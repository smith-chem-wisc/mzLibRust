//! Fetching a bridge executable, on the caller's explicit instruction and never otherwise.
//!
//! # Why this exists at all, when `build.rs` can already download
//!
//! It can, but only when *you* tell it where from. `build.rs` resolves [`crate::BRIDGE_ENV_VAR`],
//! then a staged `_dotnet/<rid>/`, then `MZLIB_BRIDGE_URL` — all caller-supplied, with no default.
//! That is deliberate: the build script never fails the build, and `cargo check`, `cargo doc`,
//! `cargo clippy` and the whole offline suite must pass on a machine that has never seen a .NET
//! binary. Giving it a built-in URL would make the crate download roughly 130 MB the first time
//! anybody compiled it, and break vendored and air-gapped builds.
//!
//! So the fetch moves out of the build and becomes something a person does: this module.
//! **Fetching the bridge is a user action, never a build side-effect** — bridge decision
//! `BRIDGE-FETCH-IS-A-USER-ACTION` (2026-08-17), which settles link `[C]` of the release chain.
//!
//! # This is the R design, not a new one
//!
//! mzLibR reached the same shape first and for stricter reasons — CRAN forbids downloading at
//! install time and forbids writing outside `tempdir()` without consent — so `R/install-bridge.R`
//! is "an exported function the user calls, which asks first". Python is silent rather than
//! opposed: a wheel simply carries the payload, so the question never arises there. This module is
//! a port of the R one, keeping its refusals, its error classes and its checks; only the idiom
//! changes.
//!
//! The payload is the raw bridge tarball pyMzLib publishes on every release,
//! `mzlib-bridge-<rid>.tar.gz`, alongside a `SHA256SUMS` manifest — the same artefact all three
//! bindings take, per `BRIDGE-DISTRIBUTION`. How a binding *obtains* it is shared surface; where it
//! caches it and how it asks first is this crate's own business.
//!
//! ```no_run
//! # fn main() -> Result<(), mzlib::MzLibError> {
//! let bridge = mzlib::install::install_bridge()?;
//! println!("installed {}", bridge.display());
//! # Ok(())
//! # }
//! ```

use std::io::IsTerminal;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use crate::bridge::{executable_name, platform_tag, MzLibError, Result};

// The pyMzLib release the payload comes from, and the SHA-256 of each platform's tarball.
//
// Recorded here rather than fetched, so the checksum being verified does not come down the same
// connection as the thing it is verifying.
//
// The block below is REWRITTEN MECHANICALLY by a bridge-watch workflow, which reads the SHA256SUMS
// asset pyMzLib publishes on every release and regenerates these pins from it — the same treatment
// mzLibR's `R/install-bridge.R` already gets. Editing it by hand works and is fine for a one-off;
// the next bump overwrites the edit. Keep the markers: the workflow locates the region by them, and
// fails loudly rather than guessing if either is missing.
//
// Only the digest is recorded per platform. The asset NAME is derived from the runtime identifier
// (see `asset_name`), so a table that has drifted from its own key is a compile-or-test failure
// rather than a wrong-platform download that passes its checksum — mzLibR learned that one the
// expensive way.
//
// BEGIN generated bridge pins
/// The pyMzLib release whose published bridge this crate installs by default.
pub const MZLIB_BRIDGE_VERSION: &str = "0.1.0.dev5";

/// The platforms pyMzLib publishes a bridge for, and the SHA-256 of each tarball.
pub const BRIDGE_ASSETS: &[(&str, &str)] = &[
    (
        "win-x64",
        "6ff5ae776c3daa2b228ad78555e026335725396e3b72cfb965534a9476b6ad0e",
    ),
    (
        "osx-arm64",
        "397e1520df30f6cb06fcc27c114ea72fddf7044929826186c4ae9fec0519f4eb",
    ),
    (
        "osx-x64",
        "90fc1d23574cdfdccebaadb93a649a0b74a9eda05f298ac03924850d7513a58e",
    ),
    (
        "linux-x64",
        "fb6a15eaa1e9ea8cc719a245b2169ec88a9a33e2ad6b9302685ac56ee933741b",
    ),
];
// END generated bridge pins

/// How much the download weighs, for the consent prompt.
///
/// Measured, not guessed: a real `linux-x64` install of `v0.1.0.dev4` pulled 129.5 MB. The
/// platforms differ — the macOS x64 payload is larger, because mzLib copies native content files
/// regardless of runtime identifier (`bridge/UPSTREAM.md` U7) — so this is deliberately "about".
const APPROXIMATE_SIZE: &str = "about 130 MB";

/// The tarball that carries the bridge for a runtime identifier.
///
/// Derived rather than recorded. See the note on [`BRIDGE_ASSETS`].
pub fn asset_name(rid: &str) -> String {
    format!("mzlib-bridge-{rid}.tar.gz")
}

/// Where a published bridge tarball is downloaded from.
pub fn release_url(rid: &str, version: &str) -> String {
    format!(
        "https://github.com/smith-chem-wisc/pyMzLib/releases/download/v{version}/{}",
        asset_name(rid)
    )
}

/// The recorded SHA-256 for a runtime identifier, if pyMzLib publishes one for it.
pub fn published_sha256(rid: &str) -> Option<&'static str> {
    BRIDGE_ASSETS
        .iter()
        .find(|(known, _)| *known == rid)
        .map(|(_, sha)| *sha)
}

/// The per-user directory an installed bridge is cached in.
///
/// Follows each platform's convention, which is what `tools::R_user_dir()` does for mzLibR and what
/// a Rust user expects of a cache: `%LOCALAPPDATA%` on Windows, `~/Library/Caches` on macOS,
/// `$XDG_CACHE_HOME` or `~/.cache` elsewhere. Hand-rolled rather than taking the `dirs` crate,
/// because the whole dependency would be this one function.
pub fn cache_dir() -> Result<PathBuf> {
    let base = if cfg!(windows) {
        std::env::var_os("LOCALAPPDATA").map(PathBuf::from)
    } else if cfg!(target_os = "macos") {
        std::env::var_os("HOME").map(|home| PathBuf::from(home).join("Library").join("Caches"))
    } else {
        std::env::var_os("XDG_CACHE_HOME")
            .map(PathBuf::from)
            .or_else(|| std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".cache")))
    };

    base.map(|base| base.join("mzlib").join("bridge"))
        .ok_or_else(|| {
            MzLibError::BridgeNotFound(
                "Could not work out a cache directory for this platform (no HOME or LOCALAPPDATA). \
                 Pass a destination explicitly, or set MZLIB_BRIDGE to a bridge you already have."
                    .to_owned(),
            )
        })
}

/// How [`install_bridge_with`] should behave.
///
/// The R sibling expresses these as eight arguments with defaults; Rust has no keyword arguments,
/// so they become a struct with [`Default`]. Every field means exactly what the R parameter of the
/// same name means.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct InstallOptions {
    /// The pyMzLib release to take the payload from.
    pub version: String,
    /// Directory to unpack into. `None` uses [`cache_dir`].
    pub destination: Option<PathBuf>,
    /// `Some(true)` confirms the download without asking. `None` asks on a terminal and refuses
    /// otherwise — the tri-state R spells `consent = NA`.
    pub consent: Option<bool>,
    /// Whether to replace a bridge already installed at the destination.
    pub overwrite: bool,
    /// An explicit URL, overriding `version`. Requires `sha256`.
    pub url: Option<String>,
    /// The expected SHA-256 of the download, lowercase hex. Required with `url`.
    pub sha256: Option<String>,
    /// Suppress progress output.
    pub quiet: bool,
    /// How long to allow for the download.
    pub timeout: Duration,
}

impl Default for InstallOptions {
    fn default() -> Self {
        Self {
            version: MZLIB_BRIDGE_VERSION.to_owned(),
            destination: None,
            consent: None,
            overwrite: false,
            url: None,
            sha256: None,
            quiet: false,
            // R raises its 60-second default for the same reason: a 130 MB payload will not meet it.
            timeout: Duration::from_secs(1800),
        }
    }
}

/// Download the bridge for this platform into the per-user cache, asking first.
///
/// Equivalent to `install_bridge_with(InstallOptions::default())`.
///
/// # Errors
///
/// See [`install_bridge_with`].
pub fn install_bridge() -> Result<PathBuf> {
    install_bridge_with(InstallOptions::default())
}

/// Download, verify and unpack a bridge executable.
///
/// **It is never called for you.** Nothing in this crate installs a bridge as a side effect of
/// building, loading or calling anything — see the module documentation for why.
///
/// If you already have a bridge you do not need this at all: set [`crate::BRIDGE_ENV_VAR`] to point
/// at it. That is also how you relink a modified mzLib, which LGPL section 4 requires this crate to
/// permit; see `NOTICE`.
///
/// # Errors
///
/// - [`MzLibError::Usage`] — consent was refused or could not be asked for, or `url` was given
///   without `sha256`. mzLibRust will not install an executable it cannot verify.
/// - [`MzLibError::BridgeNotFound`] — no bridge is published for this platform, or the unpacked
///   payload does not contain one.
/// - [`MzLibError::ServiceUnavailable`] — the download did not complete.
/// - [`MzLibError::Protocol`] — the download did not match its checksum, or would not unpack.
///   Nothing is installed in either case.
pub fn install_bridge_with(options: InstallOptions) -> Result<PathBuf> {
    let rid = platform_tag()?;

    let (url, expected_sha256) = match (&options.url, &options.sha256) {
        (Some(url), Some(sha)) => (url.clone(), sha.clone()),
        (Some(_), None) => {
            return Err(MzLibError::Usage(
                "When you pass a url, pass sha256 as well. mzLibRust will not install an \
                 executable it cannot verify."
                    .to_owned(),
            ))
        }
        (None, explicit) => {
            // pyMzLib publishes four platforms; `platform_tag` maps more than that — linux-arm64
            // resolves happily and has nothing to fetch. Saying which platforms do have one, and
            // how to build the missing one, is the difference between a dead end and a next step.
            let Some(published) = published_sha256(&rid) else {
                let known: Vec<&str> = BRIDGE_ASSETS.iter().map(|(rid, _)| *rid).collect();
                return Err(MzLibError::BridgeNotFound(format!(
                    "No published bridge for {rid}.\n\
                     pyMzLib publishes one for {}.\n\
                     Build one with pyMzLib's pkg/build/publish-bridge.ps1 and set {} to it, or \
                     pass url and sha256 explicitly.",
                    known.join(", "),
                    crate::BRIDGE_ENV_VAR
                )));
            };
            (
                release_url(&rid, &options.version),
                explicit.clone().unwrap_or_else(|| published.to_owned()),
            )
        }
    };

    let destination = match &options.destination {
        Some(explicit) => explicit.clone(),
        None => cache_dir()?,
    };
    let target_dir = destination.join(&rid);
    let target = target_dir.join(executable_name());

    if target.is_file() && !options.overwrite {
        if !options.quiet {
            eprintln!(
                "A bridge is already installed at {}. Set overwrite to replace it.",
                target.display()
            );
        }
        return Ok(target);
    }

    if options.consent != Some(true) {
        ask_consent(&target_dir)?;
    }

    let archive = target_dir.join(format!("{}.download", asset_name(&rid)));
    std::fs::create_dir_all(&target_dir)?;
    // Whatever happens next, a part-downloaded or rejected archive must not be left behind.
    let _cleanup = ScopedFile(archive.clone());

    if !options.quiet {
        eprintln!("Downloading {url}");
    }
    download(&url, &archive, options.timeout).map_err(|problem| {
        MzLibError::ServiceUnavailable {
            error_type: crate::SERVICE_UNAVAILABLE_TYPE.to_owned(),
            message: format!("Could not download the bridge from {url}: {problem}"),
        }
    })?;

    let observed = sha256(&archive).map_err(|problem| {
        MzLibError::Usage(format!(
            "{problem}, so the download cannot be verified. mzLibRust will not install an \
             unverified executable. Install one of those tools, or download the bridge yourself \
             and set {}.",
            crate::BRIDGE_ENV_VAR
        ))
    })?;
    if observed != expected_sha256.trim().to_ascii_lowercase() {
        return Err(MzLibError::Protocol(format!(
            "The download does not match its expected checksum.\n  \
             expected {}\n  observed {observed}\n\
             Nothing was installed. This means the file was corrupted in transit, or is not the \
             file mzLibRust expected.",
            expected_sha256.trim().to_ascii_lowercase()
        )));
    }

    // The tarball's root IS the payload tree — the executable sits at the top alongside Resources/,
    // the Bruker and timsTOF native libraries, and libmmd. So it unpacks straight into place, with
    // no intermediate directory and no copy step.
    extract(&archive, &target_dir).map_err(|problem| {
        MzLibError::Protocol(format!("Could not unpack the downloaded bridge: {problem}"))
    })?;

    // Checked rather than assumed: the checksum proves the bytes are the ones published, not that
    // they contain what this platform needs. A tarball built for the wrong runtime identifier would
    // pass verification and then fail at the first call with a missing-file error.
    if !target.is_file() {
        return Err(MzLibError::BridgeNotFound(format!(
            "Unpacked the payload but found no {} at {}. The archive for {rid} does not carry the \
             executable at its root; this is a packaging problem upstream.",
            executable_name(),
            target.display()
        )));
    }

    make_runnable(&target);

    if !options.quiet {
        eprintln!("Installed {}", target.display());
    }
    Ok(target)
}

/// Delete a scratch file however the surrounding function leaves.
struct ScopedFile(PathBuf);

impl Drop for ScopedFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// Ask, when there is somebody to ask.
///
/// A download of this size should never appear on somebody's disk because they mistyped. R gets
/// `interactive()`; the equivalent here is whether stdin is a terminal, which is also what makes
/// the refusal testable — it is never a terminal under `cargo test`.
fn ask_consent(destination: &Path) -> Result<()> {
    let where_to = destination.display();

    if !std::io::stdin().is_terminal() {
        return Err(MzLibError::Usage(format!(
            "install_bridge() downloads {APPROXIMATE_SIZE} and writes it to\n  {where_to}\n\
             This session is not interactive, so it cannot ask. Set consent to Some(true) to \
             proceed, or set {} to a bridge you already have.",
            crate::BRIDGE_ENV_VAR
        )));
    }

    eprintln!("install_bridge() will download {APPROXIMATE_SIZE} and write it to\n  {where_to}");
    eprint!("Proceed? [y/N] ");
    let mut answer = String::new();
    std::io::stdin().read_line(&mut answer)?;
    if !matches!(answer.trim().to_ascii_lowercase().as_str(), "y" | "yes") {
        return Err(MzLibError::Usage(
            "Cancelled; nothing was downloaded or written.".to_owned(),
        ));
    }
    Ok(())
}

/// Fetch a file, shelling out rather than taking an HTTP stack as a dependency.
///
/// The same reasoning as `build.rs`: `curl` ships with Windows 10+, macOS and essentially every
/// Linux, with PowerShell as the Windows fallback. Pulling `reqwest` and its TLS stack in to copy
/// one file would cost every downstream build far more than this is worth.
fn download(url: &str, destination: &Path, timeout: Duration) -> std::result::Result<(), String> {
    let seconds = timeout.as_secs().max(1).to_string();

    let curl = Command::new("curl")
        .args(["--fail", "--location", "--silent", "--show-error"])
        .args(["--max-time", &seconds])
        .arg("--output")
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
    let status = Command::new("powershell")
        .args(["-NoProfile", "-NonInteractive", "-Command", &script])
        .status()
        .map_err(|error| format!("could not run powershell: {error}"))?;

    if status.success() {
        Ok(())
    } else {
        Err(format!("powershell download exited with {status}"))
    }
}

/// The SHA-256 of a file, as lowercase hex.
///
/// Verification is never skipped quietly: a 130 MB executable that is about to be run is not
/// something to accept on the strength of HTTPS alone. If no hasher exists this returns an error
/// and the install refuses, exactly as mzLibR does.
fn sha256(path: &Path) -> std::result::Result<String, String> {
    let display = path.display().to_string();
    let candidates: [(&str, Vec<String>); 3] = [
        ("sha256sum", vec![display.clone()]),
        ("shasum", vec!["-a".into(), "256".into(), display.clone()]),
        (
            "certutil",
            vec!["-hashfile".into(), display, "SHA256".into()],
        ),
    ];

    for (program, args) in candidates {
        let Ok(output) = Command::new(program).args(&args).output() else {
            continue;
        };
        if !output.status.success() {
            continue;
        }
        let text = String::from_utf8_lossy(&output.stdout).to_ascii_lowercase();
        // sha256sum/shasum print "<digest>  <path>"; certutil puts the digest on its own line,
        // sometimes spaced, which is why the whole-text fallback below exists.
        if let Some(digest) = text
            .split_whitespace()
            .find(|token| token.len() == 64 && token.chars().all(|c| c.is_ascii_hexdigit()))
        {
            return Ok(digest.to_owned());
        }
        let joined: String = text.chars().filter(|c| c.is_ascii_hexdigit()).collect();
        if joined.len() >= 64 {
            return Ok(joined[..64].to_owned());
        }
    }

    Err("no usable SHA-256 tool found (tried sha256sum, shasum, certutil)".to_owned())
}

/// Unpack the payload tree.
///
/// `tar` ships with Windows 10+, macOS and every Linux, and it is what makes this one call rather
/// than a zip's four: tar records the execute bit and restores it, where an unzip drops it.
fn extract(archive: &Path, into: &Path) -> std::result::Result<(), String> {
    let status = Command::new("tar")
        .arg("-xzf")
        .arg(archive)
        .arg("-C")
        .arg(into)
        .status()
        .map_err(|error| format!("could not run tar: {error}"))?;

    if status.success() {
        Ok(())
    } else {
        Err(format!("tar exited with {status}"))
    }
}

/// Make a freshly unpacked bridge actually runnable on this platform.
///
/// The chmod is belt and braces: tar records the mode and the published payload lists
/// `mzlib-bridge` as `-rwxr-xr-x`. It is kept because it costs nothing and covers someone pointing
/// `url` at an archive built by something less careful.
///
/// The macOS half is load-bearing and has nothing to do with archive format. A file downloaded by
/// an application inherits a quarantine attribute, so Gatekeeper refuses to run it and reports that
/// the developer cannot be verified — which sounds like a problem with mzLib rather than a side
/// effect of downloading.
fn make_runnable(path: &Path) {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        if let Ok(metadata) = std::fs::metadata(path) {
            let mut permissions = metadata.permissions();
            permissions.set_mode(permissions.mode() | 0o755);
            let _ = std::fs::set_permissions(path, permissions);
        }
    }

    if cfg!(target_os = "macos") {
        // Best effort: if the attribute is absent xattr fails, and that is the good case.
        let _ = Command::new("xattr")
            .args(["-d", "com.apple.quarantine"])
            .arg(path)
            .output();
    }

    #[cfg(not(unix))]
    {
        let _ = path;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Nothing here downloads anything. The parts worth testing are the refusals — mzLibRust must
    // never install an executable it has not verified, and must never write 130 MB anywhere without
    // being asked — and every one of those is reachable before a byte moves.

    #[test]
    fn every_published_platform_has_a_tarball_and_a_checksum() {
        let mut rids: Vec<&str> = BRIDGE_ASSETS.iter().map(|(rid, _)| *rid).collect();
        rids.sort_unstable();
        assert_eq!(rids, ["linux-x64", "osx-arm64", "osx-x64", "win-x64"]);

        for (rid, sha) in BRIDGE_ASSETS {
            // A checksum is 64 lowercase hex characters. A truncated or upper-cased one would
            // compare unequal against a correct download and produce a very confusing failure.
            assert_eq!(sha.len(), 64, "{rid}");
            assert!(
                sha.chars()
                    .all(|c| c.is_ascii_hexdigit() && !c.is_uppercase()),
                "{rid}: {sha}"
            );
        }
    }

    #[test]
    fn the_asset_name_is_derived_from_the_runtime_identifier() {
        // mzLibR recorded the asset name beside the digest and could therefore drift from its own
        // key — which fetches another platform's binary and passes its checksum. Deriving it makes
        // that unrepresentable; this pins the derivation against the names pyMzLib publishes.
        assert_eq!(asset_name("win-x64"), "mzlib-bridge-win-x64.tar.gz");
        assert_eq!(asset_name("linux-x64"), "mzlib-bridge-linux-x64.tar.gz");
    }

    #[test]
    fn the_payload_is_the_raw_tarball_not_a_python_wheel() {
        // mzLibR unzipped a release WHEEL until 2026-08-09, because no neutral artefact existed.
        // pyMzLib #31 publishes mzlib-bridge-<rid>.tar.gz for exactly this, so no binding has to
        // reach through another language's packaging format.
        let url = release_url("win-x64", MZLIB_BRIDGE_VERSION);
        assert!(url.starts_with("https://github.com/smith-chem-wisc/pyMzLib/releases/download/"));
        assert!(url.ends_with("/mzlib-bridge-win-x64.tar.gz"));
        assert!(!url.contains(".whl"));
    }

    #[test]
    fn the_pinned_version_is_in_the_url_exactly_once() {
        // The version exists in one place; the URL is built from it rather than restating it.
        let url = release_url("linux-x64", MZLIB_BRIDGE_VERSION);
        assert!(
            url.contains(&format!("/download/v{MZLIB_BRIDGE_VERSION}/")),
            "{url}"
        );
    }

    #[test]
    fn linux_arm64_has_no_published_bridge() {
        // platform_tag_for maps it happily, so there is genuinely nothing to fetch and the caller
        // must be told which platforms do have one.
        assert!(published_sha256("linux-arm64").is_none());
        assert!(published_sha256("linux-x64").is_some());
    }

    #[test]
    fn a_url_without_a_checksum_is_refused() {
        let error = install_bridge_with(InstallOptions {
            url: Some("https://example.invalid/bridge.tar.gz".to_owned()),
            consent: Some(true),
            ..Default::default()
        })
        .unwrap_err();

        assert!(matches!(error, MzLibError::Usage(_)), "{error:?}");
        assert!(error.to_string().contains("pass sha256"), "{error}");
    }

    #[test]
    fn a_non_interactive_session_without_consent_refuses_before_downloading() {
        // stdin is never a terminal under `cargo test`, which is what makes this testable.
        let destination = std::env::temp_dir().join("mzlibrust-consent-refusal");
        let error = install_bridge_with(InstallOptions {
            destination: Some(destination.clone()),
            ..Default::default()
        })
        .unwrap_err();

        assert!(matches!(error, MzLibError::Usage(_)), "{error:?}");
        let message = error.to_string();
        assert!(message.contains("not interactive"), "{message}");
        assert!(message.contains("consent"), "{message}");
        // It refused before creating anything.
        assert!(!destination.exists());
    }

    #[test]
    fn the_consent_message_names_the_size_and_the_destination() {
        // Both are what a person needs in order to decide. A bare "proceed?" is not consent to
        // anything in particular.
        let destination = std::env::temp_dir().join("mzlibrust-consent-message");
        let error = ask_consent(&destination).unwrap_err();
        let message = error.to_string();
        assert!(message.contains(APPROXIMATE_SIZE), "{message}");
        assert!(
            message.contains(&destination.display().to_string()),
            "{message}"
        );
    }

    #[test]
    fn an_already_installed_bridge_is_not_replaced_without_being_asked() {
        let destination = std::env::temp_dir().join("mzlibrust-already-installed");
        let rid = platform_tag().expect("this platform has a runtime identifier");
        let target_dir = destination.join(&rid);
        std::fs::create_dir_all(&target_dir).expect("scratch directory");
        let existing = target_dir.join(executable_name());
        std::fs::write(&existing, b"not a real executable").expect("scratch file");

        let returned = install_bridge_with(InstallOptions {
            destination: Some(destination.clone()),
            quiet: true,
            ..Default::default()
        })
        .expect("an installed bridge is returned, not re-fetched");

        assert_eq!(returned, existing);
        // Untouched: it did not download, and it did not overwrite.
        assert_eq!(
            std::fs::read(&existing).expect("still there"),
            b"not a real executable"
        );

        let _ = std::fs::remove_dir_all(&destination);
    }

    #[test]
    fn sha256_is_computable_on_this_machine() {
        // The installer refuses rather than skipping verification, so on a machine with no hasher
        // it would never install at all. Confirm this one can, and against a known answer.
        let scratch = std::env::temp_dir().join("mzlibrust-hash-abc");
        std::fs::write(&scratch, b"abc").expect("scratch file");

        let digest = sha256(&scratch).expect("some hasher exists");
        assert_eq!(
            digest,
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );

        let _ = std::fs::remove_file(&scratch);
    }

    #[test]
    fn the_cache_directory_is_under_the_platforms_own_cache_root() {
        // Not asserting an exact path — the point is that it is derived from the platform's cache
        // location and namespaced, rather than dropped in the current directory.
        let Ok(cache) = cache_dir() else {
            return; // A machine with no HOME; the error path is the documented behaviour.
        };
        assert!(
            cache.ends_with(PathBuf::from("mzlib").join("bridge")),
            "{}",
            cache.display()
        );
    }
}
