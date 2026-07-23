//! The transport story, end to end, against a real bridge executable.
//!
//! This is the M0 proof: locate the executable, run it, parse an envelope, agree on a wire format.
//! It needs no network — only a staged bridge — so unlike the PRIDE and UniProt canaries it cannot
//! be excused by an outage. If this fails, the distribution story is broken and nothing above it
//! can be trusted.
//!
//! Run with `cargo test --features live`, having set `MZLIB_BRIDGE` to a staged bridge (or staged
//! one under `_dotnet/<rid>/`).

#![cfg(feature = "live")]

mod support;

#[test]
fn the_bridge_and_the_crate_agree_on_the_protocol() {
    let Some(_bridge) = support::require_bridge() else {
        return;
    };

    let info = mzlib::bridge_version().expect("the staged bridge should report its version");

    assert_eq!(
        info.protocol,
        mzlib::PROTOCOL_VERSION,
        "the crate and the bridge were built from different sources"
    );
    assert!(!info.bridge.is_empty(), "bridge version should be reported");
    assert!(!info.runtime.is_empty(), "runtime should be reported");
}

#[test]
fn an_unknown_verb_is_a_usage_error_naming_the_known_ones() {
    let Some(_bridge) = support::require_bridge() else {
        return;
    };

    let error = mzlib::bridge::invoke(&["not-a-verb".to_owned()], None, None)
        .expect_err("an unknown verb must not succeed");

    assert!(
        matches!(error, mzlib::MzLibError::Usage(_)),
        "expected a usage error, got {error:?}"
    );
    let message = error.to_string();
    assert!(message.contains("pride files"), "{message}");
    assert!(message.contains("quant flashlfq"), "{message}");
}
