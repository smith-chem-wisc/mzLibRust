//! Locating and invoking the mzLib bridge executable.
//!
//! This module is the only place in mzLibRust that knows the bridge exists. Everything above it
//! sees ordinary Rust functions and types. That boundary is deliberate: the transport (today, a
//! self-contained .NET executable invoked per call) can be replaced by an in-process binding or a
//! long-lived local server without any public API changing.
//!
//! It is the direct counterpart of pyMzLib's `_bridge.py`, and for the same reason: the wire
//! contract is language-neutral by design (decision D6), so a second binding is a transport module
//! and some typed structs, not a second implementation of mzLib.

use std::collections::HashMap;
use std::ffi::OsStr;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use serde::Deserialize;

/// Wire-format version this crate understands.
///
/// The bridge reports its own; a mismatch means the two halves were built from different sources.
pub const PROTOCOL_VERSION: u32 = 1;

/// Environment variable pointing at a bridge executable, overriding the bundled one.
///
/// The counterpart of pyMzLib's `PYMZLIB_BRIDGE`. Used during development, before a self-contained
/// binary has been staged into the package, and as the offline escape hatch.
pub const BRIDGE_ENV_VAR: &str = "MZLIB_BRIDGE";

/// The error type the bridge uses for availability failures.
///
/// Must match `Program.ServiceUnavailableType` on the C# side — the two halves agree by this
/// string, and nothing else.
pub const SERVICE_UNAVAILABLE_TYPE: &str = "ServiceUnavailable";

/// Everything mzLibRust can fail with.
///
/// Python needs a hierarchy (`PyMzLibError` → `BridgeError` → `ServiceUnavailableError`) so a
/// caller can write one `except` and be done. Rust gets that for free: one enum, matched
/// exhaustively, and the compiler will not let a new variant be forgotten at a call site.
///
/// The classification itself is made in the **bridge**, not here, so every consumer of the wire
/// format inherits it and not only this crate.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum MzLibError {
    /// A call was malformed — a missing or invalid argument. Raised before any work happens.
    ///
    /// Either this crate rejected the input, or the bridge answered `{"type": "usage"}` (exit 2).
    #[error("{0}")]
    Usage(String),

    /// An external service is unavailable — down, rate-limited, timing out, or unreachable.
    ///
    /// Deliberately a distinct variant, because the difference between "the repository is having a
    /// bad morning" and "something is broken" is the difference between retrying later and filing
    /// a bug. HTTP 408, 429 and 5xx count; 404 and 400 do not, because a wrong URL or a malformed
    /// request is our problem and excusing it as an outage would hide a real bug.
    #[error("{message}")]
    ServiceUnavailable {
        /// Always [`SERVICE_UNAVAILABLE_TYPE`]; carried so the two variants destructure alike.
        error_type: String,
        /// What the bridge said.
        message: String,
    },

    /// mzLib reported a failure.
    ///
    /// `error_type` is the .NET exception type name, e.g. `HttpRequestException` — useful for
    /// distinguishing a network failure from a bad accession without parsing prose.
    #[error("{message}")]
    Bridge {
        /// The .NET exception type name.
        error_type: String,
        /// What the bridge said.
        message: String,
    },

    /// The bridge process did not finish within the timeout.
    ///
    /// Deliberately **not** [`MzLibError::ServiceUnavailable`], and the distinction is the whole
    /// point. A subprocess timeout has several possible causes and only one of them is a slow
    /// service: the bridge may be wedged, the executable may be corrupt, antivirus may be holding
    /// it, or the caller may simply have passed a timeout that was too short. Reporting all of
    /// that as "the repository is down" is how a real bug gets skipped by every test suite and
    /// never seen again.
    #[error(
        "mzLib bridge did not finish within {seconds}s. This may mean the service is slow, but it \
         can equally mean the bridge is wedged or the timeout was too short — mzLibRust will not \
         guess which."
    )]
    Timeout {
        /// The timeout that elapsed, in seconds.
        seconds: f64,
    },

    /// The bridge executable could not be located.
    ///
    /// In a released crate this should be impossible. It normally means mzLibRust is being used
    /// from a source checkout where no bridge has been staged yet — set [`BRIDGE_ENV_VAR`].
    #[error("{0}")]
    BridgeNotFound(String),

    /// The bridge produced output this version cannot interpret.
    ///
    /// Empty output from a process that died, something that is not JSON, or a protocol version
    /// this crate does not speak.
    #[error("{0}")]
    Protocol(String),

    /// No project with that accession exists, or it has no files.
    ///
    /// PRIDE answers an unknown accession with an empty result rather than a 404, so a naive
    /// binding returns an empty list. That is a mistake worth naming: an empty list is
    /// indistinguishable from "this project genuinely has nothing matching", so a typo'd accession
    /// produces a script that reports "0 files, done" and moves on. A wrong answer that looks like
    /// a right answer is worse than an error.
    #[error("{0}")]
    ProjectNotFound(String),

    /// The bridge could not be run at all — a missing execute bit, a quarantined binary.
    ///
    /// Without this the caller would see a bare [`std::io::Error`], contradicting the promise that
    /// every failure from this crate is an [`MzLibError`].
    #[error("Could not run the mzLib bridge: {0}")]
    Io(#[from] std::io::Error),
}

/// A convenient alias: every fallible call in this crate returns one of these.
pub type Result<T> = std::result::Result<T, MzLibError>;

// ------------------------------------------------------------------ locating the payload

/// The subdirectory name a bridge build is staged under, for the current platform.
///
/// These strings must equal the .NET runtime identifiers `publish-bridge.ps1` stages under, or
/// nothing is ever found.
pub fn platform_tag() -> Result<String> {
    platform_tag_for(std::env::consts::OS, std::env::consts::ARCH)
}

/// The staging subdirectory for an explicitly named platform.
///
/// Split out from [`platform_tag`] so the mapping can be tested without a machine of each kind.
/// Both Rust's own spellings (`windows`/`macos`, `x86_64`/`aarch64`) and the ones pyMzLib's
/// `platform` module produces (`Windows`/`Darwin`, `AMD64`/`arm64`) are accepted, so the two
/// bindings can be checked against one table.
pub fn platform_tag_for(os: &str, arch: &str) -> Result<String> {
    let prefix = match os.to_ascii_lowercase().as_str() {
        "windows" => "win",
        "linux" => "linux",
        "macos" | "darwin" => "osx",
        _ => {
            return Err(MzLibError::BridgeNotFound(format!(
                "Unsupported platform: {os} {arch}"
            )))
        }
    };

    let lowered = arch.to_ascii_lowercase();
    let suffix = match lowered.as_str() {
        "amd64" | "x86_64" | "x64" => "x64",
        "arm64" | "aarch64" => "arm64",
        other => other,
    };

    Ok(format!("{prefix}-{suffix}"))
}

/// The name of the bridge executable on this platform.
fn executable_name() -> &'static str {
    if cfg!(windows) {
        "mzlib-bridge.exe"
    } else {
        "mzlib-bridge"
    }
}

/// The path of the bridge executable that will be used.
///
/// Resolution order: the [`BRIDGE_ENV_VAR`] environment variable, then the copy staged inside this
/// crate for the current platform.
///
/// # Errors
///
/// [`MzLibError::BridgeNotFound`] if neither exists.
pub fn bridge_path() -> Result<PathBuf> {
    if let Some(override_path) = std::env::var_os(BRIDGE_ENV_VAR) {
        let candidate = PathBuf::from(&override_path);
        if !candidate.is_file() {
            return Err(MzLibError::BridgeNotFound(format!(
                "{BRIDGE_ENV_VAR} points at '{}', which is not a file.",
                candidate.display()
            )));
        }
        return Ok(candidate);
    }

    // Whatever `build.rs` staged — a payload it downloaded, or one it found already in the source
    // tree. Checked before the source-tree path so a released crate, whose sources are read-only,
    // still resolves.
    if let Some(staged) = option_env!("MZLIB_BRIDGE_STAGED") {
        let candidate = PathBuf::from(staged);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }

    let candidate = staged_root().join(platform_tag()?).join(executable_name());
    if !candidate.is_file() {
        return Err(MzLibError::BridgeNotFound(missing_bridge_message(
            &candidate,
        )));
    }
    Ok(candidate)
}

/// What to tell someone who has no bridge staged.
///
/// Split out from [`bridge_path`] so it can be tested for its content without the test depending on
/// whether the machine running it happens to have a bridge staged — which made the original test
/// pass in CI and fail for any developer who had one.
fn missing_bridge_message(candidate: &Path) -> String {
    format!(
        "No mzLib bridge for this platform at '{}'.\n\
         \n\
         Three ways to fix it, cheapest first:\n\
           1. Set {BRIDGE_ENV_VAR} to a bridge executable you already have — for example the one \
         pyMzLib stages under pkg/python/src/pymzlib/_dotnet/<rid>/.\n\
           2. Run scripts/stage-bridge.ps1, or build one with pyMzLib's \
         pkg/build/publish-bridge.ps1 and stage it at '{}'.\n\
           3. Set MZLIB_BRIDGE_URL (and MZLIB_BRIDGE_SHA256) before building, and the build script \
         will download it.",
        candidate.display(),
        candidate.display()
    )
}

/// Where a staged bridge lives relative to the crate.
///
/// `_dotnet/<rid>/mzlib-bridge` mirrors the layout inside pyMzLib's wheel, so one
/// `publish-bridge.ps1` run stages a payload both bindings can consume.
fn staged_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("_dotnet")
}

// ------------------------------------------------------------------ running the process

/// What a completed bridge invocation produced.
#[derive(Debug, Clone, Default)]
pub(crate) struct Output {
    pub stdout: String,
    pub stderr: String,
    pub code: Option<i32>,
}

/// How a bridge invocation is actually executed.
///
/// This is the seam pyMzLib gets by monkeypatching `subprocess.run`. Rust has no monkeypatching,
/// so the indirection is explicit: the failure paths that matter — a process that dies silently,
/// output that is not JSON, a timeout, a binary that will not launch — are exactly the code most
/// likely to be wrong and the least convenient to provoke with a real executable.
pub(crate) trait Runner {
    /// Run the bridge once and collect everything it produced.
    fn run(
        &self,
        exe: &Path,
        args: &[String],
        stdin: Option<&str>,
        timeout: Option<Duration>,
    ) -> Result<Output>;
}

/// The real runner: spawns the executable and talks to it over pipes.
pub(crate) struct ProcessRunner;

impl Runner for ProcessRunner {
    fn run(
        &self,
        exe: &Path,
        args: &[String],
        stdin: Option<&str>,
        timeout: Option<Duration>,
    ) -> Result<Output> {
        let mut child = Command::new(exe)
            .args(args.iter().map(OsStr::new))
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|error| {
                // A missing execute bit, a quarantined binary, or a file that is not an executable
                // at all. Naming the path is the difference between a puzzling errno and an
                // actionable message.
                MzLibError::Io(std::io::Error::new(
                    error.kind(),
                    format!("'{}': {error}", exe.display()),
                ))
            })?;

        // stdin, stdout and stderr are drained on their own threads. Writing a large payload while
        // the child blocks on a full stdout pipe is the classic subprocess deadlock, and the
        // FlashLFQ verb takes a whole experiment's worth of runs on stdin. `subprocess.run` does
        // exactly this for Python; nothing about the hazard is language-specific.
        let payload = stdin.map(str::to_owned);
        let mut child_stdin = child.stdin.take();
        let writer = std::thread::spawn(move || {
            if let (Some(mut handle), Some(text)) = (child_stdin.take(), payload) {
                // A child that exits before reading everything gives a broken pipe. That is not an
                // error worth reporting: whatever it decided on the strength of what it did read
                // is in the envelope, which is the answer we came for.
                let _ = handle.write_all(text.as_bytes());
            }
            // Dropping the handle closes stdin, which is what tells the bridge the input ended.
        });

        let mut stdout_pipe = child.stdout.take();
        let stdout_reader = std::thread::spawn(move || {
            let mut buffer = String::new();
            if let Some(handle) = stdout_pipe.as_mut() {
                let _ = handle.read_to_string(&mut buffer);
            }
            buffer
        });

        let mut stderr_pipe = child.stderr.take();
        let stderr_reader = std::thread::spawn(move || {
            let mut buffer = String::new();
            if let Some(handle) = stderr_pipe.as_mut() {
                let _ = handle.read_to_string(&mut buffer);
            }
            buffer
        });

        let deadline = timeout.map(|limit| Instant::now() + limit);
        let status = loop {
            match child.try_wait()? {
                Some(status) => break status,
                None => {
                    if deadline.is_some_and(|end| Instant::now() >= end) {
                        // Kill first, then drain: the reader threads only finish once the pipes
                        // close, and they close when the process is gone.
                        let _ = child.kill();
                        let _ = child.wait();
                        let _ = writer.join();
                        let _ = stdout_reader.join();
                        let _ = stderr_reader.join();
                        return Err(MzLibError::Timeout {
                            seconds: timeout.map_or(0.0, |limit| limit.as_secs_f64()),
                        });
                    }
                    std::thread::sleep(Duration::from_millis(5));
                }
            }
        };

        let _ = writer.join();
        let stdout = stdout_reader.join().unwrap_or_default();
        let stderr = stderr_reader.join().unwrap_or_default();

        Ok(Output {
            stdout,
            stderr,
            code: status.code(),
        })
    }
}

// ------------------------------------------------------------------ the envelope

/// The single object the bridge writes to stdout, whatever the verb.
#[derive(Debug, Deserialize)]
struct Envelope {
    #[serde(default)]
    ok: bool,
    #[serde(default)]
    data: Option<serde_json::Value>,
    #[serde(default)]
    error: Option<ErrorInfo>,
}

/// A failure, described in terms a non-.NET caller can act on.
#[derive(Debug, Deserialize)]
struct ErrorInfo {
    #[serde(default)]
    r#type: Option<String>,
    #[serde(default)]
    message: Option<String>,
}

/// Reject a timeout that cannot mean anything, before spawning a process.
///
/// Python has to check for zero, negatives, `inf`, `nan`, strings and booleans, because
/// `subprocess` accepts them all and then fails somewhere unrecognisable. [`Duration`] makes every
/// one of those a compile error except zero, which still has to be caught: a zero timeout looks
/// exactly like a service that never answered.
fn validate_timeout(timeout: Option<Duration>) -> Result<()> {
    match timeout {
        Some(limit) if limit.is_zero() => Err(MzLibError::Usage(
            "timeout must be greater than zero; pass None to wait indefinitely.".to_owned(),
        )),
        _ => Ok(()),
    }
}

/// Run one bridge command and return the decoded `data` payload.
///
/// `args` is the command and its options, e.g. `["pride", "files", "--accession", "PXD000001"]`.
/// `stdin` carries a payload that would not fit on the command line — argv has a hard ceiling of
/// roughly 32 KB, and a few thousand file names exceed it. `timeout` of `None` waits indefinitely,
/// which is the right default for a large download.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the command or its arguments were malformed; [`MzLibError::Bridge`] or
/// [`MzLibError::ServiceUnavailable`] if mzLib itself failed; [`MzLibError::Protocol`] if the
/// bridge produced output this version cannot interpret.
pub fn invoke(
    args: &[String],
    stdin: Option<&str>,
    timeout: Option<Duration>,
) -> Result<serde_json::Value> {
    invoke_with(&ProcessRunner, args, stdin, timeout)
}

/// [`invoke`], against an explicit runner. The seam the transport tests use.
pub(crate) fn invoke_with(
    runner: &dyn Runner,
    args: &[String],
    stdin: Option<&str>,
    timeout: Option<Duration>,
) -> Result<serde_json::Value> {
    validate_timeout(timeout)?;
    for arg in args {
        if arg.contains('\0') {
            return Err(MzLibError::Usage(
                "Arguments may not contain a null character.".to_owned(),
            ));
        }
    }

    let exe = bridge_path()?;
    let output = runner.run(&exe, args, stdin, timeout)?;
    decode(&output)
}

/// Turn a completed invocation into the `data` payload, or the right kind of error.
fn decode(output: &Output) -> Result<serde_json::Value> {
    let stdout = output.stdout.trim();
    if stdout.is_empty() {
        // A silent non-zero exit means the process died before it could report anything — surface
        // stderr, which is the only evidence left.
        let stderr = output.stderr.trim();
        return Err(MzLibError::Protocol(format!(
            "mzLib bridge exited with code {} and no output. stderr: {}",
            output
                .code
                .map_or_else(|| "unknown".to_owned(), |code| code.to_string()),
            if stderr.is_empty() { "(empty)" } else { stderr }
        )));
    }

    let envelope: Envelope = serde_json::from_str(stdout).map_err(|_| {
        MzLibError::Protocol(format!(
            "mzLib bridge returned output that is not JSON: {}",
            truncate(stdout, 400)
        ))
    })?;

    if envelope.ok {
        return Ok(envelope.data.unwrap_or(serde_json::Value::Null));
    }

    let error = envelope.error.unwrap_or(ErrorInfo {
        r#type: None,
        message: None,
    });
    let error_type = error.r#type.unwrap_or_else(|| "Unknown".to_owned());
    let message = error
        .message
        .unwrap_or_else(|| "mzLib reported a failure with no message.".to_owned());

    Err(match error_type.as_str() {
        "usage" => MzLibError::Usage(message),
        SERVICE_UNAVAILABLE_TYPE => MzLibError::ServiceUnavailable {
            error_type,
            message,
        },
        _ => MzLibError::Bridge {
            error_type,
            message,
        },
    })
}

/// The first `limit` characters, on a character boundary so a multi-byte name cannot panic.
fn truncate(text: &str, limit: usize) -> &str {
    match text.char_indices().nth(limit) {
        Some((index, _)) => &text[..index],
        None => text,
    }
}

// ------------------------------------------------------------------ version handshake

/// What the bridge reports about itself.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct BridgeVersion {
    /// The bridge assembly's version, e.g. `"1.0.0.0"`.
    pub bridge: String,
    /// The wire-format version it speaks. Must equal [`PROTOCOL_VERSION`].
    pub protocol: u32,
    /// The .NET runtime it is carrying.
    pub runtime: String,
    /// Which mzLib the bridge was built against, as `1.0.0+<commit>`.
    ///
    /// `None` when the bridge did not report one — either because it predates this field, or
    /// because its build recorded no source commit. A crate downloads its bridge from a release
    /// rather than building it, so this is the only way to ask which mzLib is actually running.
    ///
    /// Not a compatibility check: [`protocol`](Self::protocol) is that, and it is what
    /// [`bridge_version`] verifies. This is for reporting a run, filing a bug, or pinning a
    /// result to the library that produced it.
    #[serde(default)]
    pub mzlib: Option<String>,
}

/// The bridge's own version information, with the protocol compatibility check.
///
/// This is the whole transport story end to end — locate the executable, run it, parse an
/// envelope, agree on a wire format — in one call with no network and no arguments, which is why
/// it is the first thing to make work and the first thing to check when something is wrong.
///
/// # Errors
///
/// [`MzLibError::Protocol`] if the bridge speaks a different wire format than this crate.
pub fn bridge_version() -> Result<BridgeVersion> {
    bridge_version_with(&ProcessRunner)
}

/// [`bridge_version`], against an explicit runner.
pub(crate) fn bridge_version_with(runner: &dyn Runner) -> Result<BridgeVersion> {
    let data = invoke_with(
        runner,
        &["version".to_owned()],
        None,
        Some(Duration::from_secs(60)),
    )?;

    let reported = data.get("protocol").and_then(serde_json::Value::as_u64);
    if reported != Some(u64::from(PROTOCOL_VERSION)) {
        return Err(MzLibError::Protocol(format!(
            "mzLib bridge speaks protocol {}, but this mzLibRust expects {PROTOCOL_VERSION}. The \
             Rust crate and the bridge were built from different sources.",
            reported.map_or_else(|| "nothing".to_owned(), |value| value.to_string())
        )));
    }

    serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!(
            "mzLib bridge reported an unreadable version: {error}"
        ))
    })
}

// ------------------------------------------------------------------ shared helpers

/// Render an optional float the way the bridge's invariant-culture parser expects it.
///
/// The bridge parses with [`std::f64`] semantics under `InvariantCulture`, so a comma-decimal
/// locale on the host must never reach it. Rust formats floats invariantly already; this exists so
/// the intent is stated once rather than assumed at four call sites.
pub(crate) fn format_number(value: f64, name: &str) -> Result<String> {
    if !value.is_finite() {
        return Err(MzLibError::Usage(format!(
            "{name} must be a finite number; got {value}."
        )));
    }
    Ok(format!("{value:?}"))
}

/// Read `null` as the type's default rather than failing.
///
/// The bridge emits `null` for a value mzLib had none of — a checksum a repository does not
/// publish, an organism a result file never carried. A caller wants `""` there, not a parse error
/// on a field they were not going to read.
pub(crate) fn null_to_default<'de, D, T>(deserializer: D) -> std::result::Result<T, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de> + Default,
{
    Ok(Option::<T>::deserialize(deserializer)?.unwrap_or_default())
}

/// A map whose `null` values become `0.0`.
///
/// FlashLFQ's peptide intensities are `0.0` when missing and never null; the bridge nonetheless
/// routes every double through a finite-check, so a null is theoretically reachable. Reading it as
/// zero keeps the peptide type a plain `f64` — which is the point, since the whole reason protein
/// intensities are `Option<f64>` is that only *they* can mean "could not be resolved".
pub(crate) fn deserialize_intensities<'de, D>(
    deserializer: D,
) -> std::result::Result<HashMap<String, f64>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let raw: HashMap<String, Option<f64>> = HashMap::deserialize(deserializer)?;
    Ok(raw
        .into_iter()
        .map(|(key, value)| (key, value.unwrap_or(0.0)))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    /// A runner that returns a canned result without executing anything.
    struct StubRunner {
        output: Result<Output>,
        seen: Mutex<Option<(Vec<String>, Option<String>)>>,
    }

    impl StubRunner {
        fn returning(stdout: &str, stderr: &str, code: i32) -> Self {
            Self {
                output: Ok(Output {
                    stdout: stdout.to_owned(),
                    stderr: stderr.to_owned(),
                    code: Some(code),
                }),
                seen: Mutex::new(None),
            }
        }

        fn failing(error: MzLibError) -> Self {
            Self {
                output: Err(error),
                seen: Mutex::new(None),
            }
        }
    }

    impl Runner for StubRunner {
        fn run(
            &self,
            _exe: &Path,
            args: &[String],
            stdin: Option<&str>,
            _timeout: Option<Duration>,
        ) -> Result<Output> {
            *self.seen.lock().unwrap() = Some((args.to_vec(), stdin.map(str::to_owned)));
            match &self.output {
                Ok(output) => Ok(output.clone()),
                Err(MzLibError::Timeout { seconds }) => {
                    Err(MzLibError::Timeout { seconds: *seconds })
                }
                Err(error) => Err(MzLibError::Protocol(error.to_string())),
            }
        }
    }

    /// Point `bridge_path()` at a file that exists but is never actually executed.
    ///
    /// The env var is process-global, so these tests take a lock rather than racing each other.
    static ENV_LOCK: Mutex<()> = Mutex::new(());

    struct FakeBridge {
        _guard: std::sync::MutexGuard<'static, ()>,
        _dir: tempdir::TempDir,
    }

    fn fake_bridge() -> FakeBridge {
        let guard = ENV_LOCK
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        let dir = tempdir::TempDir::new();
        let stub = dir.path().join(executable_name());
        std::fs::write(&stub, b"not a real executable").unwrap();
        std::env::set_var(BRIDGE_ENV_VAR, &stub);
        FakeBridge {
            _guard: guard,
            _dir: dir,
        }
    }

    /// A minimal scratch directory, so the crate takes no dev-dependency for four tests.
    mod tempdir {
        use std::path::{Path, PathBuf};
        use std::sync::atomic::{AtomicU32, Ordering};

        static COUNTER: AtomicU32 = AtomicU32::new(0);

        pub struct TempDir(PathBuf);

        impl TempDir {
            pub fn new() -> Self {
                let unique = COUNTER.fetch_add(1, Ordering::Relaxed);
                let path = std::env::temp_dir()
                    .join(format!("mzlibrust-test-{}-{unique}", std::process::id()));
                std::fs::create_dir_all(&path).unwrap();
                Self(path)
            }

            pub fn path(&self) -> &Path {
                &self.0
            }
        }

        impl Drop for TempDir {
            fn drop(&mut self) {
                let _ = std::fs::remove_dir_all(&self.0);
            }
        }
    }

    fn envelope(text: &str) -> Result<serde_json::Value> {
        decode(&Output {
            stdout: text.to_owned(),
            stderr: String::new(),
            code: Some(0),
        })
    }

    // ---------------------------------------------------------- locating the payload

    #[test]
    fn env_override_takes_precedence() {
        let _fake = fake_bridge();
        let expected = std::env::var_os(BRIDGE_ENV_VAR).unwrap();
        assert_eq!(bridge_path().unwrap(), PathBuf::from(expected));
    }

    #[test]
    fn env_override_pointing_at_nothing_says_so() {
        let _guard = ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        std::env::set_var(BRIDGE_ENV_VAR, "/definitely/does/not/exist/mzlib-bridge");
        let error = bridge_path().unwrap_err();
        assert!(matches!(error, MzLibError::BridgeNotFound(_)));
        assert!(error.to_string().contains("not a file"));
        std::env::remove_var(BRIDGE_ENV_VAR);
    }

    #[test]
    fn missing_payload_names_the_path_and_the_way_out() {
        // The message has to be actionable: a bare "not found" leaves the user with nowhere to go.
        //
        // Asserted against the message builder rather than by calling `bridge_path()` with the
        // environment stripped. The original version did the latter, which meant it only passed on
        // a machine with NO bridge staged — green in CI, red for any developer who had staged one.
        // A test whose result depends on the developer's machine is worse than no test.
        let message =
            missing_bridge_message(Path::new("/somewhere/_dotnet/linux-x64/mzlib-bridge"));

        assert!(message.contains("_dotnet"), "{message}");
        assert!(message.contains(BRIDGE_ENV_VAR), "{message}");
        // All three remedies, because naming only one sends the reader looking for the others.
        assert!(message.contains("MZLIB_BRIDGE_URL"), "{message}");
        assert!(message.contains("stage-bridge.ps1"), "{message}");
    }

    #[test]
    fn a_staged_bridge_is_found_without_any_environment_variable() {
        // The distribution promise: if `build.rs` staged something, the crate finds it with no
        // MZLIB_BRIDGE set. Skips rather than fails when nothing is staged, since the offline
        // suite must pass on a machine that has never seen a .NET binary.
        let _guard = ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        std::env::remove_var(BRIDGE_ENV_VAR);

        match bridge_path() {
            Ok(found) => assert!(found.is_file(), "{} should exist", found.display()),
            Err(error) => {
                let message = error.to_string();
                assert!(
                    message.contains("No mzLib bridge"),
                    "with nothing staged the failure must be the actionable one: {message}"
                );
            }
        }
    }

    #[test]
    fn platform_tags_match_dotnet_runtime_identifiers() {
        // These strings must equal the RIDs publish-bridge.ps1 stages under, or nothing is found.
        for (os, arch, expected) in [
            ("Windows", "AMD64", "win-x64"),
            ("Linux", "x86_64", "linux-x64"),
            ("Linux", "aarch64", "linux-arm64"),
            ("Darwin", "arm64", "osx-arm64"),
            ("Darwin", "x86_64", "osx-x64"),
            // Rust's own spellings must land on the same tags as Python's.
            ("windows", "x86_64", "win-x64"),
            ("macos", "aarch64", "osx-arm64"),
        ] {
            assert_eq!(platform_tag_for(os, arch).unwrap(), expected, "{os}/{arch}");
        }
    }

    #[test]
    fn unsupported_platform_is_reported_not_guessed() {
        let error = platform_tag_for("Plan9", "x86_64").unwrap_err();
        assert!(error.to_string().contains("Unsupported platform"));
    }

    // ---------------------------------------------------------- the envelope

    #[test]
    fn success_returns_only_the_data() {
        let data = envelope(r#"{"ok":true,"data":{"a":1}}"#).unwrap();
        assert_eq!(data, serde_json::json!({"a": 1}));
    }

    #[test]
    fn usage_failure_maps_to_usage_error() {
        let error = envelope(
            r#"{"ok":false,"error":{"type":"usage","message":"Missing required option --x."}}"#,
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(error.to_string().contains("Missing required option"));
    }

    #[test]
    fn runtime_failure_preserves_the_dotnet_error_type() {
        // `error_type` is what lets a caller tell a network blip from a bad accession.
        let error =
            envelope(r#"{"ok":false,"error":{"type":"HttpRequestException","message":"503"}}"#)
                .unwrap_err();
        match error {
            MzLibError::Bridge { error_type, .. } => assert_eq!(error_type, "HttpRequestException"),
            other => panic!("expected Bridge, got {other:?}"),
        }
    }

    #[test]
    fn service_unavailable_is_its_own_variant() {
        let error = envelope(
            r#"{"ok":false,"error":{"type":"ServiceUnavailable","message":"EBI is down"}}"#,
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::ServiceUnavailable { .. }));
    }

    #[test]
    fn failure_with_no_message_still_raises_something_readable() {
        let error = envelope(r#"{"ok":false}"#).unwrap_err();
        assert!(error.to_string().contains("no message"));
    }

    // ---------------------------------------------------------- when things go wrong

    #[test]
    fn silent_death_surfaces_stderr() {
        // A process that dies before writing anything leaves stderr as the only evidence.
        let error = decode(&Output {
            stdout: String::new(),
            stderr: "Segmentation fault".to_owned(),
            code: Some(139),
        })
        .unwrap_err();
        let message = error.to_string();
        assert!(message.contains("139"), "{message}");
        assert!(message.contains("Segmentation fault"), "{message}");
    }

    #[test]
    fn silent_death_with_no_stderr_does_not_produce_an_empty_message() {
        let error = decode(&Output {
            stdout: String::new(),
            stderr: String::new(),
            code: Some(1),
        })
        .unwrap_err();
        assert!(error.to_string().contains("(empty)"));
    }

    #[test]
    fn non_json_output_is_quoted_back_not_swallowed() {
        // If the bridge prints something unexpected, the user needs to see what it was.
        let error = envelope("Unhandled exception. System.Whatever").unwrap_err();
        let message = error.to_string();
        assert!(message.contains("not JSON"), "{message}");
        assert!(message.contains("System.Whatever"), "{message}");
    }

    #[test]
    fn timeout_becomes_a_typed_error_not_a_process_artifact() {
        // Nothing above `bridge` should ever have to reason about process handles.
        let _fake = fake_bridge();
        let runner = StubRunner::failing(MzLibError::Timeout { seconds: 5.0 });
        let error = invoke_with(
            &runner,
            &["anything".to_owned()],
            None,
            Some(Duration::from_secs(5)),
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Timeout { .. }));
    }

    #[test]
    fn a_timeout_is_not_reported_as_a_service_outage() {
        // The regression this guards is subtle and was live in pyMzLib: every subprocess timeout
        // used to raise ServiceUnavailable, which the canary suites turn into a skip. A wedged
        // bridge, a corrupt binary, or a caller passing too small a timeout all reported "EBI is
        // down" and the live tests passed green — precisely the failure the testing convention
        // exists to prevent.
        let _fake = fake_bridge();
        let runner = StubRunner::failing(MzLibError::Timeout { seconds: 1.0 });
        let error = invoke_with(
            &runner,
            &["anything".to_owned()],
            None,
            Some(Duration::from_secs(1)),
        )
        .unwrap_err();
        assert!(!matches!(error, MzLibError::ServiceUnavailable { .. }));
    }

    #[test]
    fn a_zero_timeout_is_rejected_before_spawning_anything() {
        // Python must also reject negatives, inf, nan, strings and booleans; `Duration` makes
        // every one of those a compile error, leaving zero as the only reachable case — and zero
        // is indistinguishable from a service that never answered.
        let _fake = fake_bridge();
        let runner = StubRunner::returning("", "", 0);
        let error = invoke_with(
            &runner,
            &["anything".to_owned()],
            None,
            Some(Duration::ZERO),
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(runner.seen.lock().unwrap().is_none(), "should not have run");
    }

    #[test]
    fn an_unlaunchable_bridge_is_an_mzlib_error() {
        // A missing execute bit or a quarantined binary must not escape as a bare io::Error.
        let _fake = fake_bridge();
        let error = invoke(&["anything".to_owned()], None, None).unwrap_err();
        // The stub file exists but is not executable code; whatever the platform makes of that,
        // it has to arrive as one of ours.
        assert!(matches!(error, MzLibError::Io(_) | MzLibError::Protocol(_)));
    }

    #[test]
    fn null_bytes_are_rejected_rather_than_reaching_the_process() {
        let _fake = fake_bridge();
        let runner = StubRunner::returning("", "", 0);
        let error = invoke_with(
            &runner,
            &["pride".to_owned(), "PX\0D".to_owned()],
            None,
            None,
        )
        .unwrap_err();
        assert!(error.to_string().contains("null character"));
        assert!(runner.seen.lock().unwrap().is_none(), "should not have run");
    }

    #[test]
    fn every_error_is_one_type() {
        // A user should be able to write `match error { … }` and be done. In Python this needs a
        // test because the hierarchy could drift; here the compiler enforces it, so the test
        // documents the guarantee rather than defending it.
        fn assert_is_error<E: std::error::Error>(_: &E) {}
        assert_is_error(&MzLibError::Usage("x".to_owned()));
        assert_is_error(&MzLibError::Protocol("x".to_owned()));
    }

    // ---------------------------------------------------------- version handshake

    #[test]
    fn matching_protocol_returns_the_version_info() {
        let _fake = fake_bridge();
        let runner = StubRunner::returning(
            r#"{"ok":true,"data":{"bridge":"1.0.0.0","protocol":1,"runtime":"8.0.27"}}"#,
            "",
            0,
        );
        let info = bridge_version_with(&runner).unwrap();
        assert_eq!(info.bridge, "1.0.0.0");
        assert_eq!(info.protocol, PROTOCOL_VERSION);
        assert_eq!(info.runtime, "8.0.27");
        // This payload carries no `mzlib`, which is exactly what a bridge built before the field
        // existed sends. It must still deserialize: an added wire field may not strand a caller on
        // an older bridge, which is the whole reason the field is optional rather than required.
        assert_eq!(info.mzlib, None);
    }

    #[test]
    fn the_mzlib_build_is_reported_when_the_bridge_sends_it() {
        let _fake = fake_bridge();
        let runner = StubRunner::returning(
            r#"{"ok":true,"data":{"bridge":"1.0.0.0","protocol":1,"runtime":"8.0.27","mzlib":"1.0.0+f6b0f0d17f32383918ef895006aaecb71cdb9a7e"}}"#,
            "",
            0,
        );
        let info = bridge_version_with(&runner).unwrap();
        assert_eq!(
            info.mzlib.as_deref(),
            Some("1.0.0+f6b0f0d17f32383918ef895006aaecb71cdb9a7e")
        );
    }

    #[test]
    fn mismatched_protocol_fails_loudly() {
        // Halves built from different sources must not silently produce wrong results.
        let _fake = fake_bridge();
        let runner = StubRunner::returning(
            r#"{"ok":true,"data":{"bridge":"1.0.0.0","protocol":2,"runtime":"8.0.27"}}"#,
            "",
            0,
        );
        let error = bridge_version_with(&runner).unwrap_err();
        assert!(error.to_string().contains("different sources"));
    }

    #[test]
    fn the_version_verb_is_what_gets_sent() {
        let _fake = fake_bridge();
        let runner = StubRunner::returning(
            r#"{"ok":true,"data":{"bridge":"1.0.0.0","protocol":1,"runtime":"8.0.27"}}"#,
            "",
            0,
        );
        bridge_version_with(&runner).unwrap();
        let seen = runner.seen.lock().unwrap().clone().unwrap();
        assert_eq!(seen.0, vec!["version".to_owned()]);
        assert_eq!(seen.1, None);
    }

    #[test]
    fn stdin_reaches_the_runner_unchanged() {
        let _fake = fake_bridge();
        let runner = StubRunner::returning(r#"{"ok":true,"data":{}}"#, "", 0);
        invoke_with(
            &runner,
            &["quant".to_owned(), "flashlfq".to_owned()],
            Some("run_1.mzML\tcondition\n"),
            None,
        )
        .unwrap();
        let seen = runner.seen.lock().unwrap().clone().unwrap();
        assert_eq!(seen.1.as_deref(), Some("run_1.mzML\tcondition\n"));
    }
}
