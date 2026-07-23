//! Shared test support.
//!
//! The important thing here is [`external_service`] — the Rust counterpart of pyMzLib's
//! `conftest.py` helper, and of mzLib's own `ExternalServiceTestHelper`. Tests that touch a live
//! service must be able to tell two failures apart:
//!
//! * the service is unavailable (down, rate-limited, 5xx, timed out) — **not** our bug, so the
//!   test should be **skipped** with a message saying so; versus
//! * the service answered but the contract is broken (wrong URL, response no longer parses, an
//!   expected value missing) — a real regression that must **fail**.
//!
//! Without that distinction a red build is ambiguous, and an ambiguous red build gets ignored,
//! which is how a genuine contract break goes unnoticed for a month.
//!
//! The classification itself lives in the bridge, not here: it reports `ServiceUnavailable` for
//! availability failures, so every consumer of the wire format benefits, not just this suite.
//!
//! **On "skip":** `cargo test` has no skip verdict — a test either passes or panics. So a skip
//! here is an early return with a loud `SKIPPED` line on stderr, visible with `--nocapture` and in
//! CI logs. That is weaker than pytest's `skip`, and deliberately noted rather than papered over:
//! the alternative, failing on an outage, is the thing this whole convention exists to prevent.

#![allow(dead_code)]

use mzlib::MzLibError;

/// Announce a skip in a form that stands out in a CI log.
pub fn skip(reason: &str) {
    eprintln!("SKIPPED: {reason}");
}

/// Ensure a bridge executable is available, or skip.
///
/// Returns `None` when there is nothing to run against, so the caller returns early. A missing
/// bridge is a *staging* problem, not a code failure — the same reasoning as an outage.
pub fn require_bridge() -> Option<()> {
    match mzlib::bridge_path() {
        Ok(_) => Some(()),
        Err(error) => {
            skip(&format!(
                "no mzLib bridge is staged ({error}). Set {} to a built bridge, e.g. the one \
                 pyMzLib stages under pkg/python/src/pymzlib/_dotnet/<rid>/.",
                mzlib::BRIDGE_ENV_VAR
            ));
            None
        }
    }
}

/// Unwrap a live call, skipping rather than failing when the service is unavailable.
///
/// Any other error is a genuine regression and panics, which is the whole point: a 404 from a
/// wrong URL, or a response that no longer parses, must not be excused as an outage.
pub fn external_service<T>(service: &str, result: mzlib::Result<T>) -> Option<T> {
    match result {
        Ok(value) => Some(value),
        Err(MzLibError::ServiceUnavailable { message, .. }) => {
            skip(&format!(
                "{service} unavailable ({message}). This is a third-party availability problem, \
                 not a code failure."
            ));
            None
        }
        Err(error) => panic!("{service} answered, but the contract is broken: {error}"),
    }
}
