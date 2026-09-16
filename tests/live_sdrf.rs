//! Live checks of the sdrf module against the real bridge and mzLib's own SDRF fixture.
//!
//! The offline suite in `src/sdrf.rs` pins the projection against payloads recorded by pyMzLib.
//! These pin that the recording still matches what the bridge emits — a fixture that has quietly
//! diverged makes every offline test prove something about a shape nothing produces any more.
//!
//! They need the bridge and `MZLIB_TEST_FILES`, and **skip** rather than fail when either is
//! missing. Run with `cargo test --features live`.

#![cfg(feature = "live")]

mod support;

use std::path::PathBuf;

use mzlib::sdrf::{self, SdrfDocument};
use support::require_bridge;

const RELATIVE: &str = "FileReadingTests/ExternalFileTypes/PXD000070.sdrf.tsv";

fn mzlib_sdrf() -> Option<PathBuf> {
    let path = PathBuf::from(std::env::var("MZLIB_TEST_FILES").ok()?).join(RELATIVE);
    if path.exists() {
        Some(path)
    } else {
        eprintln!(
            "skipping: set MZLIB_TEST_FILES to an mzLib Test directory containing {RELATIVE}"
        );
        None
    }
}

#[test]
fn the_recorded_read_fixture_still_matches_the_live_bridge() {
    let Some(()) = require_bridge() else { return };
    let Some(path) = mzlib_sdrf() else { return };

    let live = sdrf::read(&path).expect("mzLib's own SDRF fixture should read");
    let recorded: SdrfDocument =
        serde_json::from_str(include_str!("fixtures/sdrf_read_PXD000070.json")).unwrap();

    assert_eq!(live.columns, recorded.columns);
    assert_eq!(live.rows, recorded.rows);
    assert_eq!(live.row_count, recorded.row_count);
    assert_eq!(live.caveats, recorded.caveats);
}

#[test]
fn pooling_one_document_with_itself_keeps_both_copies_apart() {
    // Also the cheapest live proof that the stdin rendering survives a round trip, since a pooled
    // table of one repeated document is where a dropped line is most visible.
    let Some(()) = require_bridge() else { return };
    let Some(path) = mzlib_sdrf() else { return };

    let once = sdrf::pool_labelled(&[(&path, "first")]).expect("a one-document pool");
    let twice = sdrf::pool(&[&path, &path]).expect("a two-document pool");

    assert_eq!(once.document_count, 1);
    assert_eq!(twice.document_count, 2);
    assert_eq!(twice.document.row_count, 2 * once.document.row_count);
    assert!(
        twice
            .document
            .caveats
            .iter()
            .any(|c| c.contains("not reproducible")),
        "an unlabelled pool must say its provenance depends on where the files sit"
    );
}
