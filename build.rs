//! Stages the mzLib bridge executable so the crate works without a manual setup step.
//!
//! # What this does, and what it deliberately does not
//!
//! The bridge is a self-contained .NET executable of roughly 130 MB. crates.io will not carry that
//! (the limit is ~10 MB), and `include_bytes!` is obviously out, so a Rust binding has to fetch it.
//! That is the one part of the design pyMzLib's wheel did not have to solve — a wheel just carries
//! the payload.
//!
//! Resolution order, cheapest first:
//!
//! 1. **`MZLIB_BRIDGE`** — if the caller has a bridge, nothing needs staging. Checked at *runtime*
//!    too, so it always wins.
//! 2. **`_dotnet/<rid>/mzlib-bridge`** beside the crate — a source checkout where someone has run
//!    pyMzLib's `publish-bridge.ps1`. This is the current lab workflow.
//! 3. **Download** from `MZLIB_BRIDGE_URL`, verified against `MZLIB_BRIDGE_SHA256` when given.
//!
//! **This build script never fails the build.** A missing bridge is a *runtime* problem with a good
//! error message, not a compile-time one — `cargo check`, `cargo doc`, `cargo clippy` and the whole
//! offline test suite must work on a machine that has never seen a .NET binary, and they do. Failing
//! here would make the crate undocumentable and untestable for contributors, to buy nothing: the
//! program still cannot run without the payload either way.
//!
//! Staging goes to `OUT_DIR`, not to the source tree, because a crate unpacked from the registry is
//! read-only.

use std::path::{Path, PathBuf};

fn main() {
    println!("cargo:rerun-if-env-changed=MZLIB_BRIDGE");
    println!("cargo:rerun-if-env-changed=MZLIB_BRIDGE_URL");
    println!("cargo:rerun-if-env-changed=MZLIB_BRIDGE_SHA256");
    println!("cargo:rerun-if-changed=build.rs");

    // 1. The caller already has one. Nothing to stage, and the runtime check will find it.
    if std::env::var_os("MZLIB_BRIDGE").is_some() {
        return;
    }

    let Some(rid) = runtime_identifier() else {
        warn("mzLibRust does not know a .NET runtime identifier for this target; set MZLIB_BRIDGE to a bridge you built yourself.");
        return;
    };
    let exe = executable_name();

    // 2. A source checkout with a staged payload.
    let staged = PathBuf::from(env("CARGO_MANIFEST_DIR"))
        .join("_dotnet")
        .join(&rid)
        .join(exe);
    if staged.is_file() {
        emit_staged_path(&staged);
        return;
    }

    // 3. Download, if we were told where from.
    let Some(url) = std::env::var("MZLIB_BRIDGE_URL")
        .ok()
        .filter(|u| !u.is_empty())
    else {
        warn(&format!(
            "no mzLib bridge staged for {rid}. The crate will compile, and the offline tests will \
             pass, but any call that reaches mzLib will fail until you either set MZLIB_BRIDGE to a \
             bridge executable, stage one at _dotnet/{rid}/{exe}, or set MZLIB_BRIDGE_URL to \
             download from. See README.md."
        ));
        return;
    };

    let destination = PathBuf::from(env("OUT_DIR")).join(exe);
    match download(&url, &destination) {
        Ok(()) => {
            if let Some(expected) = std::env::var("MZLIB_BRIDGE_SHA256").ok().filter(|s| !s.is_empty())
            {
                match verify(&destination, &expected) {
                    Ok(()) => {}
                    Err(problem) => {
                        // A payload that fails its checksum must not be left where the runtime
                        // would pick it up and run it.
                        let _ = std::fs::remove_file(&destination);
                        warn(&format!("discarded the downloaded bridge: {problem}"));
                        return;
                    }
                }
            } else {
                warn("MZLIB_BRIDGE_URL was used without MZLIB_BRIDGE_SHA256, so the download was not verified. Pin a checksum for anything reproducible.");
            }
            make_executable(&destination);
            emit_staged_path(&destination);
        }
        Err(problem) => warn(&format!(
            "could not download the mzLib bridge from {url}: {problem}. The crate still builds; set \
             MZLIB_BRIDGE at runtime to point at one."
        )),
    }
}

/// Tell the crate where the staged bridge landed, as a compile-time fallback path.
fn emit_staged_path(path: &Path) {
    println!("cargo:rustc-env=MZLIB_BRIDGE_STAGED={}", path.display());
}

fn env(key: &str) -> String {
    std::env::var(key).unwrap_or_else(|_| panic!("cargo should set {key}"))
}

fn warn(message: &str) {
    println!("cargo:warning={message}");
}

/// The .NET runtime identifier for the target being built.
///
/// These strings have to match what `publish-bridge.ps1` stages under, or nothing is ever found.
fn runtime_identifier() -> Option<String> {
    let os = match env("CARGO_CFG_TARGET_OS").as_str() {
        "windows" => "win",
        "linux" => "linux",
        "macos" => "osx",
        _ => return None,
    };
    let arch = match env("CARGO_CFG_TARGET_ARCH").as_str() {
        "x86_64" => "x64",
        "aarch64" => "arm64",
        _ => return None,
    };
    Some(format!("{os}-{arch}"))
}

fn executable_name() -> &'static str {
    if env("CARGO_CFG_TARGET_OS") == "windows" {
        "mzlib-bridge.exe"
    } else {
        "mzlib-bridge"
    }
}

/// Fetch a file, shelling out rather than taking an HTTP stack as a build dependency.
///
/// `curl` ships with Windows 10+, macOS and essentially every Linux; PowerShell is the fallback on
/// Windows. Pulling `reqwest` and its TLS stack into the *build* graph to copy one file would cost
/// every downstream build far more than this is worth.
fn download(url: &str, destination: &Path) -> Result<(), String> {
    if let Some(parent) = destination.parent() {
        std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }

    let curl = std::process::Command::new("curl")
        .args([
            "--fail",
            "--location",
            "--silent",
            "--show-error",
            "--output",
        ])
        .arg(destination)
        .arg(url)
        .status();

    match curl {
        Ok(status) if status.success() => return Ok(()),
        Ok(status) => {
            if !cfg!(windows) {
                return Err(format!("curl exited with {status}"));
            }
        }
        Err(error) => {
            if !cfg!(windows) {
                return Err(format!("could not run curl: {error}"));
            }
        }
    }

    let script = format!(
        "$ProgressPreference='SilentlyContinue'; Invoke-WebRequest -Uri '{url}' -OutFile '{}'",
        destination.display()
    );
    let status = std::process::Command::new("powershell")
        .args(["-NoProfile", "-NonInteractive", "-Command", &script])
        .status()
        .map_err(|error| format!("could not run powershell: {error}"))?;

    if status.success() {
        Ok(())
    } else {
        Err(format!("powershell download exited with {status}"))
    }
}

/// Check a downloaded payload against an expected SHA-256, shelling out for the digest.
fn verify(path: &Path, expected: &str) -> Result<(), String> {
    let expected = expected.trim().to_ascii_lowercase();
    let actual = sha256(path)?;
    if actual == expected {
        Ok(())
    } else {
        Err(format!(
            "checksum mismatch: expected {expected}, got {actual}"
        ))
    }
}

fn sha256(path: &Path) -> Result<String, String> {
    // One of these exists on every platform we target.
    let candidates: [(&str, Vec<String>); 3] = [
        ("sha256sum", vec![path.display().to_string()]),
        (
            "shasum",
            vec!["-a".into(), "256".into(), path.display().to_string()],
        ),
        (
            "certutil",
            vec![
                "-hashfile".into(),
                path.display().to_string(),
                "SHA256".into(),
            ],
        ),
    ];

    for (program, args) in candidates {
        let Ok(output) = std::process::Command::new(program).args(&args).output() else {
            continue;
        };
        if !output.status.success() {
            continue;
        }
        let text = String::from_utf8_lossy(&output.stdout).to_ascii_lowercase();
        // sha256sum/shasum print "<digest>  <path>"; certutil puts the digest on its own line.
        if let Some(digest) = text
            .split_whitespace()
            .find(|token| token.len() == 64 && token.chars().all(|c| c.is_ascii_hexdigit()))
        {
            return Ok(digest.to_owned());
        }
    }

    Err("no usable SHA-256 tool found (tried sha256sum, shasum, certutil)".to_owned())
}

/// Downloaded files arrive without the execute bit on Unix.
fn make_executable(path: &Path) {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        if let Ok(metadata) = std::fs::metadata(path) {
            let mut permissions = metadata.permissions();
            permissions.set_mode(permissions.mode() | 0o755);
            let _ = std::fs::set_permissions(path, permissions);
        }
    }
    #[cfg(not(unix))]
    {
        let _ = path;
    }
}
