//! Live checks of the library module against the real bridge and the two `.msl` libraries the
//! recordings were made from (`tests/fixtures/library/`, shared byte for byte with pyMzLib).
//!
//! The offline suite in `src/library.rs` pins the projection against those recordings. These pin
//! that the recordings still match what the bridge emits. They need the bridge from pyMzLib 0.6.0,
//! the first with `library read`, and **skip** rather than fail without one.
//!
//! Run with `cargo test --features live`.

#![cfg(feature = "live")]

mod support;

use std::path::PathBuf;

use mzlib::library::{self, LibraryBatch, LibraryBulkOptions, LibraryReadOptions, LibraryRecords};
use mzlib::readers::BulkOptions;
use mzlib::{MzLibError, OnError};
use support::require_verb;

fn msl(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("library")
        .join(name)
}

fn recorded<T: serde::de::DeserializeOwned>(json: &str) -> T {
    serde_json::from_str(json).unwrap()
}

fn scratch(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("mzlib-live-library-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The facts a recording pins, without its machine-specific absolute path.
fn assert_same(live: &LibraryRecords, want: &LibraryRecords) {
    assert_eq!(live.columns, want.columns);
    assert_eq!(
        (live.record_count, live.returned_count, live.row_count),
        (want.record_count, want.returned_count, want.row_count)
    );
    assert_eq!(
        (live.target_count, live.decoy_count, live.proteoform_count),
        (want.target_count, want.decoy_count, want.proteoform_count)
    );
    assert_eq!(
        (live.is_index_only, live.index_only_threshold),
        (want.is_index_only, want.index_only_threshold)
    );
    assert_eq!(live.absent_fields, want.absent_fields);
    assert_eq!(live.caveats, want.caveats);
}

#[test]
fn the_recordings_are_what_the_bridge_says_now() {
    let Some(()) = require_verb("library read") else {
        return;
    };
    let cases: [(&str, LibraryRecords, LibraryReadOptions); 4] = [
        (
            "myPrositLib.msl",
            recorded(include_str!("fixtures/library_read_prosit.json")),
            LibraryReadOptions::default(),
        ),
        (
            "synthetic_library.msl",
            recorded(include_str!("fixtures/library_read_synthetic.json")),
            LibraryReadOptions::default(),
        ),
        (
            "synthetic_library.msl",
            recorded(include_str!("fixtures/library_read_fragments.json")),
            LibraryReadOptions {
                fragments: true,
                limit: Some(1),
                ..Default::default()
            },
        ),
        (
            "synthetic_library.msl",
            recorded(include_str!("fixtures/library_read_index_only.json")),
            LibraryReadOptions {
                index_only_threshold: Some(0),
                limit: Some(2),
                ..Default::default()
            },
        ),
    ];
    for (file, want, options) in cases {
        let live = library::read_with(msl(file), &options).unwrap();
        assert_same(&live, &want);
    }
}

#[test]
fn the_batch_recording_is_what_the_bridge_says_now() {
    let Some(()) = require_verb("library read") else {
        return;
    };
    let live = library::read_many(
        &[msl("synthetic_library.msl"), msl("myPrositLib.msl")],
        &LibraryBulkOptions::default(),
    )
    .unwrap();
    let want: LibraryBatch = recorded(include_str!("fixtures/library_read_many.json"));
    assert_eq!(
        (live.file_count, live.read_count, live.record_count),
        (want.file_count, want.read_count, want.record_count)
    );
    for column in live.columns.names() {
        if column != "source_path" {
            assert_eq!(
                live.columns.raw(column),
                want.columns.raw(column),
                "{column}"
            );
        }
    }
    let absent = |b: &LibraryBatch| -> Vec<Vec<String>> {
        b.files.iter().map(|f| f.absent_fields.clone()).collect()
    };
    assert_eq!(absent(&live), absent(&want));
}

#[test]
fn a_batch_is_identical_at_one_and_four_threads() {
    let Some(()) = require_verb("library read") else {
        return;
    };
    let read = |threads: i32| {
        library::read_many(
            &[msl("synthetic_library.msl"), msl("myPrositLib.msl")],
            &LibraryBulkOptions {
                bulk: BulkOptions {
                    threads,
                    ..Default::default()
                },
                fragments: true,
                ..Default::default()
            },
        )
        .unwrap()
    };
    let (one, four) = (read(1), read(4));
    assert_eq!(one.columns, four.columns);
    assert_eq!(one.files, four.files);
    assert_eq!(one.returned_count, 11);
    assert!(one.row_count > one.returned_count);
}

#[test]
fn index_only_and_full_load_give_the_same_fragment_rows() {
    let Some(()) = require_verb("library read") else {
        return;
    };
    let read = |threshold: Option<u64>| {
        library::read_with(
            msl("synthetic_library.msl"),
            &LibraryReadOptions {
                fragments: true,
                index_only_threshold: threshold,
                ..Default::default()
            },
        )
        .unwrap()
    };
    let (full, index_only) = (read(None), read(Some(0)));
    assert_eq!(
        (full.is_index_only, index_only.is_index_only),
        (false, true)
    );
    assert_eq!(index_only.columns, full.columns);
}

#[test]
fn out_writes_the_table() {
    let Some(()) = require_verb("library read") else {
        return;
    };
    let dir = scratch("out");
    let target = dir.join("nested").join("ions.tsv");
    let written = library::read_with(
        msl("synthetic_library.msl"),
        &LibraryReadOptions {
            fragments: true,
            out: Some(target.to_string_lossy().into_owned()),
            ..Default::default()
        },
    )
    .unwrap();
    let text = std::fs::read_to_string(&target).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert!(written.columns.is_empty());
    let output = written.output.as_ref().expect("out was given");
    assert_eq!(lines.len() - 1, output.row_count as usize);
    assert_eq!(output.row_count, 14);
    assert!(lines[0].split('\t').any(|c| c == "fragment_mz"));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_text_library_is_a_usage_error() {
    let Some(()) = require_verb("library read") else {
        return;
    };
    let dir = scratch("msp");
    let text = dir.join("library.msp");
    std::fs::write(&text, "Name: PEPTIDE/2\n").unwrap();
    let error = library::read(&text).unwrap_err();
    assert!(
        matches!(&error, MzLibError::Usage(m) if m.contains(".msl")),
        "{error:?}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_damaged_library_is_mzlibs_correctness_error() {
    let Some(()) = require_verb("library read") else {
        return;
    };
    let dir = scratch("damaged");
    let damaged = dir.join("damaged.msl");
    let mut bytes = std::fs::read(msl("synthetic_library.msl")).unwrap();
    bytes[200] ^= 0xFF; // inside the data the CRC-32 covers
    std::fs::write(&damaged, bytes).unwrap();
    match library::read(&damaged).unwrap_err() {
        MzLibError::Bridge { error_type, .. } => assert!(
            error_type == "InvalidDataException" || error_type == "FormatException",
            "{error_type}"
        ),
        other => panic!("expected mzLib's correctness error, got {other:?}"),
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_skipped_library_is_reported_and_the_rest_read() {
    let Some(()) = require_verb("library read") else {
        return;
    };
    let dir = scratch("skip");
    let batch = library::read_many(
        &[msl("synthetic_library.msl"), dir.join("missing.msl")],
        &LibraryBulkOptions {
            bulk: BulkOptions {
                on_error: OnError::Skip,
                ..Default::default()
            },
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!((batch.read_count, batch.failed_count), (1, 1));
    let failed = batch.failed_files();
    assert_eq!(failed.len(), 1);
    assert_eq!(failed[0].error.as_ref().unwrap().kind, "usage");
    assert_eq!(
        (failed[0].record_count, failed[0].target_count),
        (None, None)
    );
    let _ = std::fs::remove_dir_all(&dir);
}
