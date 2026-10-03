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

// ---- mzLib 1.0.592: validate, lint, assess, samples, parse_ages ------------------------------
//
// These run against the SDRF documents the recordings were made from (tests/fixtures/, shared
// with pyMzLib), so a recording that has quietly diverged from the bridge fails here. They need
// the bridge from pyMzLib 0.2.0, and skip on an older one.

fn cohort_fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

fn recorded<T: serde::de::DeserializeOwned>(text: &str) -> T {
    serde_json::from_str(text).unwrap()
}

#[test]
fn validate_still_matches_its_recordings() {
    let Some(()) = support::require_verb("sdrf validate") else {
        return;
    };
    for (name, text) in [
        (
            "sdrf_cohort.sdrf.tsv",
            include_str!("fixtures/sdrf_validate_cohort.json"),
        ),
        (
            "sdrf_skeleton.sdrf.tsv",
            include_str!("fixtures/sdrf_validate_skeleton.json"),
        ),
    ] {
        let live = sdrf::validate(cohort_fixture(name)).expect("validate should answer");
        let want: sdrf::SdrfValidation = recorded(text);
        assert_eq!(live.is_valid, want.is_valid, "{name}");
        assert_eq!(live.error_count, want.error_count, "{name}");
        assert_eq!(live.warning_count, want.warning_count, "{name}");
        assert_eq!(live.row_count, want.row_count, "{name}");
        assert_eq!(live.messages().unwrap(), want.messages().unwrap(), "{name}");
    }
}

#[test]
fn a_missing_document_is_skipped_not_fatal() {
    let Some(()) = support::require_verb("sdrf validate") else {
        return;
    };
    let batch = sdrf::validate_many(
        &[
            cohort_fixture("sdrf_skeleton.sdrf.tsv"),
            cohort_fixture("sdrf_cohort.sdrf.tsv"),
            cohort_fixture("missing.sdrf.tsv"),
        ],
        &sdrf::BulkOptions {
            on_error: sdrf::OnError::Skip,
            threads: 2,
            ..Default::default()
        },
    )
    .expect("a skipped document is not an error");
    assert_eq!(
        (batch.file_count, batch.read_count, batch.valid_count),
        (3, 2, 1)
    );
    assert_eq!(batch.files[2].error.as_ref().unwrap().kind, "usage");
}

#[test]
fn lint_still_matches_its_recording() {
    let Some(()) = support::require_verb("sdrf lint") else {
        return;
    };
    let live = sdrf::lint_labelled(&[
        (cohort_fixture("sdrf_cohort.sdrf.tsv"), "cohort"),
        (cohort_fixture("sdrf_cohort_partner.sdrf.tsv"), "partner"),
    ])
    .expect("lint should answer");
    let want: sdrf::SdrfDrift = recorded(include_str!("fixtures/sdrf_lint_cohort.json"));
    assert_eq!(live.findings().unwrap(), want.findings().unwrap());
}

#[test]
fn assess_many_still_matches_its_recording() {
    let Some(()) = support::require_verb("sdrf assess") else {
        return;
    };
    let want: sdrf::SdrfAssessmentBatch = recorded(include_str!("fixtures/sdrf_assess_bulk.json"));
    // PXD000070 is mzLib's own fixture; skip it when mzLib's test files are not at hand.
    let Some(pxd) = std::env::var_os("MZLIB_TEST_FILES")
        .map(|root| PathBuf::from(root).join(RELATIVE))
        .filter(|p| p.exists())
    else {
        eprintln!(
            "skipping: set MZLIB_TEST_FILES to an mzLib Test directory containing {RELATIVE}"
        );
        return;
    };
    let live = sdrf::assess_many(
        &[
            cohort_fixture("sdrf_cohort.sdrf.tsv"),
            cohort_fixture("sdrf_skeleton.sdrf.tsv"),
            pxd,
        ],
        &sdrf::BulkOptions::default(),
    )
    .expect("assess_many should answer");
    assert_eq!(live.verdict_counts, want.verdict_counts);
    let verdicts = |b: &sdrf::SdrfAssessmentBatch| -> Vec<Option<String>> {
        b.files.iter().map(|f| f.verdict.clone()).collect()
    };
    assert_eq!(verdicts(&live), verdicts(&want));
    assert_eq!(live.record_count, want.record_count);
}

#[test]
fn samples_and_ages_still_match_their_recordings() {
    let Some(()) = support::require_verb("sdrf samples") else {
        return;
    };
    let live = sdrf::samples(cohort_fixture("sdrf_cohort.sdrf.tsv")).expect("samples");
    let want: sdrf::SdrfSamples = recorded(include_str!("fixtures/sdrf_samples_cohort.json"));
    assert_eq!(live.sample_count, want.sample_count);
    assert_eq!(live.conflicts().unwrap(), want.conflicts().unwrap());
    assert_eq!(live.ages().unwrap(), want.ages().unwrap());

    let Some(()) = support::require_verb("sdrf parse-age") else {
        return;
    };
    let cells = [
        "58Y",
        "30Y6M",
        "40Y-85Y",
        "40Y-40Y",
        ">=90Y",
        "<1Y",
        "6-8 weeks",
        "63",
        "not available",
        "",
        "about forty",
    ];
    let live = sdrf::parse_ages(&cells).expect("parse_ages");
    let want: sdrf::ParsedAges = recorded(include_str!("fixtures/sdrf_parse_age.json"));
    assert_eq!(live.ages().unwrap(), want.ages().unwrap());
}

// ---- sdrf design ----------------------------------------------------------------------------
//
// The SDRFs are mzLib's own design fixtures, copied byte for byte from pyMzLib into
// tests/fixtures, so these need only the bridge.

fn crate_fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

fn both_factors() -> sdrf::DesignOptions {
    sdrf::DesignOptions {
        condition_columns: vec![
            "factor value[genotype]".to_owned(),
            "factor value[treatment]".to_owned(),
        ],
        ..Default::default()
    }
}

#[test]
fn each_design_recording_still_matches_the_live_bridge() {
    let Some(()) = support::require_verb("sdrf design") else {
        return;
    };
    for (sdrf_file, recording) in [
        ("PXD067622.sdrf.tsv", "sdrf_design_PXD067622.json"),
        ("PXD049018.sdrf.tsv", "sdrf_design_PXD049018.json"),
        ("PXD067622_studywide.sdrf.tsv", "sdrf_design_studywide.json"),
    ] {
        let live = sdrf::design_with(crate_fixture(sdrf_file), &both_factors())
            .unwrap_or_else(|error| panic!("{sdrf_file}: {error}"));
        let text = std::fs::read_to_string(crate_fixture(recording)).unwrap();
        let recorded: sdrf::SdrfDesign = serde_json::from_str(&text).unwrap();

        assert_eq!(live.is_valid, recorded.is_valid, "{sdrf_file}");
        assert_eq!(live.file_count, recorded.file_count, "{sdrf_file}");
        assert_eq!(live.columns, recorded.columns, "{sdrf_file}");
        assert_eq!(live.refusals, recorded.refusals, "{sdrf_file}");
        assert_eq!(live.notes, recorded.notes, "{sdrf_file}");
        assert_eq!(live.caveats, recorded.caveats, "{sdrf_file}");
    }
}

#[test]
fn a_valid_design_writes_a_one_based_experimental_design_and_a_refused_one_writes_nothing() {
    let Some(()) = support::require_verb("sdrf design") else {
        return;
    };
    let scratch = std::env::temp_dir().join(format!("mzlib-live-design-{}", std::process::id()));
    std::fs::create_dir_all(&scratch).unwrap();
    let out = scratch.join("ExperimentalDesign.tsv");

    let d = sdrf::design_with(
        crate_fixture("PXD067622.sdrf.tsv"),
        &sdrf::DesignOptions {
            out: Some(out.clone()),
            ..both_factors()
        },
    )
    .expect("a valid design writes");
    let written = d.written.as_ref().expect("written is reported");
    assert_eq!(written.file_count, 24);
    let text = std::fs::read_to_string(&out).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some("FileName\tCondition\tBiorep\tFraction\tTechrep")
    );
    let first: Vec<&str> = lines.next().unwrap().split('\t').collect();
    // mzLib adds the one when it writes, and nowhere else.
    assert_eq!(
        first[2].parse::<u32>().unwrap(),
        d.files().unwrap()[0].biological_replicate + 1
    );

    let refused_out = scratch.join("refused.tsv");
    let refused = sdrf::design_with(
        crate_fixture("PXD049018.sdrf.tsv"),
        &sdrf::DesignOptions {
            out: Some(refused_out.clone()),
            ..both_factors()
        },
    )
    .expect("a refusal is an answer");
    assert!(!refused.is_valid);
    assert_eq!(refused.written, None);
    assert!(!refused_out.exists());

    let error = sdrf::design_with(
        crate_fixture("PXD067622.sdrf.tsv"),
        &sdrf::DesignOptions {
            out: Some(scratch.join("design.txt")),
            ..both_factors()
        },
    )
    .unwrap_err();
    assert!(
        matches!(error, mzlib::MzLibError::Usage(ref m) if m.contains(".tsv")),
        "{error}"
    );
    let _ = std::fs::remove_dir_all(&scratch);
}

#[test]
fn searched_files_restrict_the_design_and_replace_the_paths() {
    let Some(()) = support::require_verb("sdrf design") else {
        return;
    };
    let full = sdrf::design_with(crate_fixture("PXD067622.sdrf.tsv"), &both_factors()).unwrap();
    let files = full.files().unwrap();
    let searched: Vec<PathBuf> = files
        .iter()
        .map(|f| PathBuf::from("/data/runs").join(&f.full_path))
        .collect();

    let d = sdrf::design_with(
        crate_fixture("PXD067622.sdrf.tsv"),
        &sdrf::DesignOptions {
            searched_files: Some(searched.clone()),
            ..both_factors()
        },
    )
    .expect("every row's file is searched");
    assert!(d.searched_files_given);
    assert_eq!(d.file_count, 24);
    assert!(d.files().unwrap()[0].full_path.starts_with("/data/runs"));

    let undeclared = sdrf::design(crate_fixture("PXD067622.sdrf.tsv")).unwrap();
    assert!(
        !undeclared.is_valid,
        "two factor columns and none declared is a refusal"
    );
}
