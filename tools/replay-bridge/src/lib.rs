//! A stand-in mzLib bridge that answers from recorded fixtures, so mzLibRust's doc examples run.
//!
//! Every runnable example in the crate's documentation starts with a hidden line,
//!
//! ```text
//! # mzlib_replay::activate();
//! ```
//!
//! which points `MZLIB_BRIDGE` at an executable built from `src/stub.rs`. That executable answers a
//! call from a fixture in `tests/fixtures/` recorded from the real bridge, and only when the
//! recording fits the call: the same file, the same echoed options, the same `limit`/`offset`
//! window, and a bulk recording only for a bulk call. An example whose arguments fit no recording
//! fails, and the error names the recordings it tried and why each did not fit.
//!
//! mzLibRust cannot tell the stub from the real bridge: it is found, spawned and read exactly the
//! same way. So an example is executed Rust, not decoration, and needs no .NET, no network and no
//! data file. The rules are pyMzLib's `pkg/python/tests/replay_bridge.py`, so a recording answers
//! the same call in both bindings.
//!
//! This crate is a path-only dev-dependency of `mzlib` and is never published.

/// The stub executable, built by this crate's build script.
pub const BRIDGE: &str = env!("MZLIB_REPLAY_BRIDGE");

/// Point this process at the replay bridge, and move it into a scratch directory.
///
/// The scratch directory is for examples that write a file (`out`, `output_directory`): they write
/// there rather than into the source tree. Doctests each run in their own process, so neither the
/// environment variable nor the directory leaks into anything else.
pub fn activate() {
    std::env::set_var("MZLIB_BRIDGE", BRIDGE);
    let scratch = std::env::temp_dir().join(format!("mzlib-doctest-{}", std::process::id()));
    if std::fs::create_dir_all(&scratch).is_ok() {
        let _ = std::env::set_current_dir(&scratch);
    }
}
