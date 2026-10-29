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

// ---- library calibrate-rt (pyMzLib 0.7.0) ---------------------------------------------------------
//
// test_library_calibrate_live.py, test for test: each recording's call re-run on the files in
// tests/fixtures/calibration/, and mzLib's refusals. Skips on a bridge without the verb.

use library::{calibrate_rt, RtAnchor, RtCalibration, RtCalibrationOptions, RtQuery};

/// A calibration input as one map per row.
fn calibration_rows(name: &str) -> Vec<std::collections::HashMap<String, String>> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("calibration")
        .join(name);
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_owned)
        .collect();
    lines
        .map(|line| {
            header
                .iter()
                .cloned()
                .zip(line.split('\t').map(str::to_owned))
                .collect()
        })
        .collect()
}

fn number(row: &std::collections::HashMap<String, String>, column: &str) -> f64 {
    row[column].parse().unwrap()
}

/// The guide's call: q-value 0.01 or below are anchors; the rest are queried both ways.
fn hela_call() -> mzlib::Result<RtCalibration> {
    let rows = calibration_rows("diann2_hela_rt.tsv");
    let (confident, others): (Vec<_>, Vec<_>) =
        rows.iter().partition(|r| number(r, "q_value") <= 0.01);
    let anchors: Vec<RtAnchor> = confident
        .iter()
        .map(|r| {
            RtAnchor::new(
                &r["run"],
                number(r, "retention_time"),
                number(r, "indexed_retention_time"),
            )
        })
        .collect();
    let queries: Vec<RtQuery> = others
        .iter()
        .map(|r| RtQuery::irt_to_minutes(&r["run"], number(r, "indexed_retention_time")))
        .chain(
            others
                .iter()
                .map(|r| RtQuery::minutes_to_irt(&r["run"], number(r, "retention_time"))),
        )
        .collect();
    calibrate_rt(&anchors, &queries, &RtCalibrationOptions::default())
}

/// The guide's and the function example's: one precursor held out of each run's anchors and
/// queried.
fn pbmc_inputs() -> (Vec<RtAnchor>, Vec<RtQuery>) {
    let rows = calibration_rows("diann_pbmc_rt.tsv");
    let held_out = "AAAAAAAAGAFAGR2";
    let anchors = rows
        .iter()
        .filter(|r| r["precursor_id"] != held_out)
        .map(|r| {
            RtAnchor::new(
                &r["run"],
                number(r, "retention_time"),
                number(r, "indexed_retention_time"),
            )
        })
        .collect();
    let queries = rows
        .iter()
        .filter(|r| r["precursor_id"] == held_out)
        .map(|r| RtQuery::irt_to_minutes(&r["run"], 11.886096))
        .collect();
    (anchors, queries)
}

/// Equal to the recording, except that floats may differ in their last bits across platforms
/// (1e-9 relative).
fn assert_same_calibration(live: &RtCalibration, fixture: &str) {
    let recorded: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("fixtures")
                .join(fixture),
        )
        .unwrap(),
    )
    .unwrap();
    let close = |a: &serde_json::Value, b: &serde_json::Value, at: &str| match (
        a.as_f64(),
        b.as_f64(),
        a.is_boolean(),
    ) {
        (Some(x), Some(y), false) => {
            assert!((x - y).abs() <= 1e-9 * y.abs().max(1e-300), "{at}: {x} {y}");
        }
        _ => assert_eq!(a, b, "{at}"),
    };
    let names: Vec<&str> = recorded["column_names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    assert_eq!(live.columns.names(), names, "{fixture}");
    for name in names {
        let ours = live.columns.raw(name).unwrap();
        let theirs = recorded["columns"][name].as_array().unwrap();
        assert_eq!(ours.len(), theirs.len(), "{fixture} {name}");
        for (i, (a, b)) in ours.iter().zip(theirs).enumerate() {
            close(a, b, &format!("{fixture} {name}[{i}]"));
        }
    }
    let scalars = serde_json::json!({
        "kind": live.kind,
        "minimum_anchors": live.minimum_anchors,
        "bandwidth": live.bandwidth,
        "knots": live.knots,
        "robustness_iterations": live.robustness_iterations,
        "on_error": live.on_error.as_str(),
        "group_count": live.group_count,
        "fitted_count": live.fitted_count,
        "failed_count": live.failed_count,
        "anchor_count": live.anchor_count,
        "query_count": live.query_count,
        "row_count": live.row_count,
    });
    for (key, value) in scalars.as_object().unwrap() {
        close(value, &recorded[key], &format!("{fixture} {key}"));
    }
    let groups = recorded["groups"].as_array().unwrap();
    assert_eq!(live.groups.len(), groups.len(), "{fixture}");
    for (ours, theirs) in live.groups.iter().zip(groups) {
        let at = format!("{fixture} {}", ours.group);
        assert_eq!(ours.group, theirs["group"].as_str().unwrap(), "{at}");
        assert_eq!(ours.kind.as_deref(), theirs["kind"].as_str(), "{at}");
        assert_eq!(
            (ours.anchor_count, ours.query_count),
            (
                theirs["anchor_count"].as_u64().unwrap(),
                theirs["query_count"].as_u64().unwrap()
            ),
            "{at}"
        );
        for (key, value) in [
            ("residual_sd", ours.residual_sd),
            ("first_anchor_minutes", ours.first_anchor_minutes),
            ("last_anchor_minutes", ours.last_anchor_minutes),
        ] {
            close(
                &serde_json::json!(value),
                &theirs[key],
                &format!("{at} {key}"),
            );
        }
        let error = ours
            .error
            .as_ref()
            .map(|e| serde_json::json!({"kind": e.kind, "type": e.r#type, "message": e.message}));
        assert_eq!(serde_json::json!(error), theirs["error"], "{at}");
    }
    assert_eq!(
        serde_json::to_value(&live.caveats).unwrap(),
        recorded["caveats"],
        "{fixture}"
    );
}

#[test]
fn the_calibration_recordings_are_what_the_bridge_returns_today() {
    let Some(()) = require_verb("library calibrate-rt") else {
        return;
    };
    assert_same_calibration(&hela_call().unwrap(), "library_calibrate_rt_hela.json");
    let (anchors, queries) = pbmc_inputs();
    let skip = calibrate_rt(
        &anchors,
        &queries,
        &RtCalibrationOptions {
            on_error: OnError::Skip,
            ..Default::default()
        },
    )
    .unwrap();
    assert_same_calibration(&skip, "library_calibrate_rt_pbmc_skip.json");
    let linear = calibrate_rt(
        &anchors,
        &queries,
        &RtCalibrationOptions {
            kind: Some("linear".into()),
            minimum_anchors: Some(10),
            ..Default::default()
        },
    )
    .unwrap();
    assert_same_calibration(&linear, "library_calibrate_rt_pbmc_linear.json");
}

#[test]
fn calibration_refusals_cross_as_usage_errors() {
    let Some(()) = require_verb("library calibrate-rt") else {
        return;
    };
    let usage = |result: mzlib::Result<RtCalibration>, needle: &str| match result {
        Err(MzLibError::Usage(message)) => assert!(message.contains(needle), "{needle}: {message}"),
        other => panic!("expected a usage error containing {needle:?}, got {other:?}"),
    };
    let (anchors, queries) = pbmc_inputs();
    let refused = calibrate_rt(&anchors, &queries, &RtCalibrationOptions::default());
    match &refused {
        Err(MzLibError::Usage(message)) => assert!(
            message.contains("Group '20240316_PBMCs")
                && message.contains("needs at least 20 anchors; got 13"),
            "{message}"
        ),
        other => panic!("expected the first group's refusal, got {other:?}"),
    }
    usage(
        calibrate_rt(
            &anchors,
            &queries,
            &RtCalibrationOptions {
                knots: Some(2),
                on_error: OnError::Skip,
                ..Default::default()
            },
        ),
        "Knots must be at least 3",
    );
    usage(
        calibrate_rt(
            &anchors,
            &[],
            &RtCalibrationOptions {
                bandwidth: Some(1.5),
                ..Default::default()
            },
        ),
        "Bandwidth must be in",
    );
    usage(
        calibrate_rt(
            &[
                RtAnchor::new("r", 1.0, 1.0),
                RtAnchor::new("r", 2.0, f64::NAN),
                RtAnchor::new("r", 3.0, 3.0),
            ],
            &[],
            &RtCalibrationOptions {
                minimum_anchors: Some(3),
                ..Default::default()
            },
        ),
        "must be finite",
    );
}
