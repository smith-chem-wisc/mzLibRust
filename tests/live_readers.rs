//! Live checks of the readers module against the real bridge and real mzLib fixtures.
//!
//! The offline suite in `src/readers.rs` covers the projection — how a wire payload becomes a
//! [`Table`], how a `null` becomes `None`, what an absent column reports. What it cannot cover is
//! that the bridge still emits the fields this crate reads, which is exactly the drift a binding
//! is most likely to suffer.
//!
//! These need both the bridge binary and the mzLib test files, and **skip** rather than fail when
//! either is missing, so the suite stays runnable on a machine that has neither.
//!
//! Run with `cargo test --features live`.

#![cfg(feature = "live")]

mod support;

use std::path::{Path, PathBuf};

use support::{require_bridge, require_verb};

/// mzLib's own test tree, wherever the developer put it.
///
/// Located from `MZLIB_TEST_FILES` rather than guessed, because the mzLib source is not part of
/// this crate — it is a sibling checkout whose path is a local choice.
fn mzlib_test_files() -> Option<PathBuf> {
    let root = PathBuf::from(std::env::var("MZLIB_TEST_FILES").ok()?);
    root.is_dir().then_some(root)
}

fn fixture(relative: &str) -> Option<PathBuf> {
    let path = mzlib_test_files()?.join(relative);
    (path.exists()).then_some(path)
}

macro_rules! fixture_or_skip {
    ($relative:expr) => {
        match fixture($relative) {
            Some(path) => path,
            None => {
                eprintln!(
                    "skipping: set MZLIB_TEST_FILES to an mzLib Test directory containing {}",
                    $relative
                );
                return;
            }
        }
    };
}

#[test]
fn every_format_the_bridge_lists_is_one_this_crate_can_describe() {
    let Some(()) = require_bridge() else { return };

    let formats = mzlib::readers::formats().expect("the formats verb should answer");

    assert_eq!(
        formats.len(),
        36,
        "mzLib 1.0.592 recognises 36 file types; a change here means the crate's documented \
         count is stale (the bridge's readers formats spec records the count per pin)"
    );

    // 17 of 36 belong to no cross-format family, which is the fact that makes read_records
    // necessary rather than a convenience.
    let viewless = formats.iter().filter(|f| f.views.is_empty()).count();
    assert_eq!(viewless, 17, "17 of 36 have no view at all");

    // Four, not three: mzLib 1.0.585 added DiaNnReport (mzLib #1120), which is how DIA data
    // reaches read_results and FlashLFQ at all.
    let quantifiable = formats.iter().filter(|f| f.is_quantifiable()).count();
    assert_eq!(quantifiable, 4, "exactly 4 offer the quantifiable view");
}

#[test]
fn a_format_with_no_view_at_all_still_reads() {
    // The whole point of read_records: TopPIC belongs to no family and was unreachable before it.
    let Some(()) = require_bridge() else { return };
    let path =
        fixture_or_skip!("FileReadingTests/ExternalFileTypes/ToppicPrsm_TopPICv1.6.2_prsm.tsv");

    let table = mzlib::readers::read_records(&path).expect("TopPIC should read");

    assert_eq!(table.file_type, "ToppicPrsm");
    assert_eq!(table.record_type, "ToppicPrsm");
    assert!(table.views.is_empty(), "TopPIC has no cross-format view");
    assert!(table.columns.has("e_value"));
    // A pluralising 's' belongs to the acronym before it: FixedPTMs, not fixed_pt_ms.
    assert!(table.columns.has("fixed_ptms"));
    assert!(table.columns.has("mi_score"));
}

#[test]
fn a_field_that_cannot_cross_the_wire_is_named_rather_than_dropped() {
    let Some(()) = require_bridge() else { return };
    let path =
        fixture_or_skip!("FileReadingTests/ExternalFileTypes/ToppicPrsm_TopPICv1.6.2_prsm.tsv");

    let table = mzlib::readers::read_records(&path).expect("TopPIC should read");

    let excluded = table
        .excluded_fields
        .iter()
        .find(|field| field.field == "alternative_identifications")
        .expect("the composite list field must be reported, not silently dropped");
    assert!(
        !excluded.reason.is_empty(),
        "an exclusion without a reason tells a caller nothing"
    );
}

#[test]
fn every_column_is_as_long_as_the_returned_record_count() {
    // A column shorter than the others is the silent-misalignment failure the whole projection
    // exists to refuse.
    let Some(()) = require_bridge() else { return };
    let path = fixture_or_skip!("FileReadingTests/ExternalFileTypes/crux.txt");

    let table = mzlib::readers::read_records(&path).expect("Crux should read");

    let rows = usize::try_from(table.returned_count).expect("a sane row count");
    for name in table.columns.names() {
        let column = table.columns.raw(name).expect("a named column exists");
        assert_eq!(column.len(), rows, "column '{name}' is the wrong length");
    }
}

#[test]
fn the_ms1_feature_retention_time_unit_is_unknown_and_refuses_to_convert() {
    // TopFD wrote seconds through v1.6.2 and minutes from v1.7.0 without changing the file type.
    let Some(()) = require_bridge() else { return };
    let path =
        fixture_or_skip!("FileReadingTests/ExternalFileTypes/Ms1Feature_TopFDv1.6.2_ms1.feature");

    let features = mzlib::readers::read_features(&path).expect("the feature view should read");

    assert_eq!(features.retention_time_unit, "unknown");
    let error = features
        .retention_time_start_in_minutes()
        .expect_err("an unknown unit must refuse rather than guess");
    assert!(error.to_string().contains("no basis to say"), "{error}");
}

#[test]
fn the_ms1_feature_unit_really_did_change_between_topfd_versions() {
    // The evidence for the caveat above, pinned so it cannot rot into folklore.
    let Some(()) = require_bridge() else { return };
    let older =
        fixture_or_skip!("FileReadingTests/ExternalFileTypes/Ms1Feature_TopFDv1.6.2_ms1.feature");
    let newer =
        fixture_or_skip!("FileReadingTests/ExternalFileTypes/Ms1Feature_TopFDv1.7.0_ms1.feature");

    let first_start = |path: &Path| -> f64 {
        mzlib::readers::read_features(path)
            .expect("the feature view should read")
            .columns
            .floats("retention_time_start")
            .expect("a numeric column")[0]
            .expect("a value")
    };

    assert!(
        first_start(&older) > 600.0,
        "v1.6.2 writes seconds — a value beyond any plausible gradient length in minutes"
    );
    assert!(first_start(&newer) < 600.0, "v1.7.0 writes minutes");
}

#[test]
fn casanovo_is_decoy_is_none_because_de_novo_sequencing_has_no_decoys() {
    let Some(()) = require_bridge() else { return };
    let path = fixture_or_skip!("FileReadingTests/ExternalFileTypes/Casanovo_5.0.0.mztab");

    let matches = mzlib::readers::read_matches(&path).expect("Casanovo should read");

    let decoys = matches
        .columns
        .booleans("is_decoy")
        .expect("a boolean column");
    assert!(
        decoys.iter().all(Option::is_none),
        "false would be a fabricated value someone could filter on"
    );
}

#[test]
fn spectra_read_headers_by_default_and_peaks_on_request() {
    let Some(()) = require_bridge() else { return };
    let path = fixture_or_skip!("DataFiles/sliced_ethcd.mzML");

    let headers = mzlib::readers::read_spectra(&path).expect("mzML should read");
    assert!(!headers.peaks_included);
    assert!(
        !headers.columns.has("mz"),
        "peaks must be absent by default: a mid-size mzML would otherwise serialise hundreds of \
         megabytes for the ordinary 'what is in this file' call"
    );
    assert_eq!(headers.retention_time_unit, "minutes");

    let with_peaks = mzlib::readers::read_spectra_with(
        &path,
        &mzlib::readers::SpectraOptions {
            read: mzlib::readers::ReadOptions {
                limit: Some(2),
                ..Default::default()
            },
            peaks: true,
            ..Default::default()
        },
    )
    .expect("mzML should read with peaks");

    assert!(with_peaks.peaks_included);
    let peaks = with_peaks
        .columns
        .float_arrays("mz")
        .expect("one array per scan");
    let counts = with_peaks
        .columns
        .integers("peak_count")
        .expect("a whole-number column");
    for (peaks, count) in peaks.iter().zip(counts) {
        assert_eq!(
            i64::try_from(peaks.as_ref().expect("a scan's peaks").len()).unwrap(),
            count.expect("a peak count"),
            "the peak array and the reported peak count must agree"
        );
    }
}

#[test]
fn an_ms_order_filter_reports_the_files_real_total_alongside_it() {
    let Some(()) = require_bridge() else { return };
    let path = fixture_or_skip!("DataFiles/sliced_ethcd.mzML");

    let all = mzlib::readers::read_spectra(&path).expect("mzML should read");
    let ms2 = mzlib::readers::read_spectra_with(
        &path,
        &mzlib::readers::SpectraOptions {
            ms_order: Some(2),
            ..Default::default()
        },
    )
    .expect("mzML should read filtered");

    assert_eq!(
        ms2.scan_count, all.scan_count,
        "scan_count reports the file's real total, so a filter that matched nothing can never \
         look like an empty file"
    );
    assert!(ms2.record_count <= all.record_count);
    for order in ms2.columns.integers("ms_order").expect("a column") {
        assert_eq!(order, Some(2));
    }
}

#[test]
fn asking_for_a_view_a_file_does_not_have_names_the_alternative() {
    let Some(()) = require_bridge() else { return };
    let path = fixture_or_skip!("FileReadingTests/SearchResults/ExcelEditedPeptide.psmtsv");

    let error =
        mzlib::readers::read_features(&path).expect_err("a psmtsv has no ms1_features view");

    let message = error.to_string();
    assert!(message.contains("quantifiable"), "{message}");
    assert!(message.contains("read-records"), "{message}");
}

#[test]
fn a_missing_apex_intensity_is_absent_rather_than_a_fabricated_zero() {
    // A within-type schema variant: Apex_intensity is optional and the FLASHDeconv/OpenMS
    // _ms1.feature layout omits it, so mzLib substitutes zero for every feature. A whole column of
    // zeros is indistinguishable from real measurements of nothing.
    //
    // This test used to pin the bridge passing mzLib's zero through with a FABRICATED caveat, and
    // said it would flip when the wire learned to name the gap. It has: since the mzLib 1.0.592
    // bridge (pyMzLib 0.2.0, BULK.md section 4), the column is named in absent_fields and every
    // value is None.
    let Some(()) = require_verb("readers read-features") else {
        return;
    };
    let path = fixture_or_skip!(
        "FileReadingTests/ExternalFileTypes/Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature"
    );

    let features = mzlib::readers::read_features(&path).expect("the feature view should read");

    assert!(
        features.absent_fields.iter().any(|f| f == "intensity"),
        "a column this file has no source for must be named: {:?}",
        features.absent_fields
    );
    let intensities = features
        .columns
        .floats("intensity")
        .expect("a numeric column");
    assert!(
        intensities.iter().all(Option::is_none),
        "an absent column is None in every row, never mzLib's substituted zero"
    );
}

#[test]
fn a_topfd_feature_file_still_reports_real_intensities() {
    // The counterpart, and the fixture that proves the caveat above is conditional: TopFD writes
    // Apex_intensity, so its intensities are real and nothing is claimed about them.
    let Some(()) = require_bridge() else { return };
    let path =
        fixture_or_skip!("FileReadingTests/ExternalFileTypes/Ms1Feature_TopFDv1.6.2_ms1.feature");

    let features = mzlib::readers::read_features(&path).expect("the feature view should read");

    assert!(features
        .columns
        .floats("intensity")
        .expect("a numeric column")
        .iter()
        .all(|value| value.is_some_and(|intensity| intensity > 0.0)));
    assert!(!features
        .caveats
        .iter()
        .any(|caveat| caveat.contains("FABRICATED")));
}

// ---- the mzLib 1.0.592 batch -----------------------------------------------------------------
//
// These need the bridge from pyMzLib 0.2.0 (mzLib 1.0.592) and SKIP on an older one. Each compares
// the live answer with the recording pyMzLib made from the same file, so a recording that has
// quietly diverged from what the bridge emits is caught here rather than trusted offline.

fn recording(name: &str) -> serde_json::Value {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name);
    serde_json::from_str(&std::fs::read_to_string(path).expect("the recording exists"))
        .expect("the recording is JSON")
}

#[test]
fn a_spectra_file_reports_the_run_it_came_from() {
    let Some(()) = require_verb("readers read-spectra") else {
        return;
    };
    let path = fixture_or_skip!("DataFiles/sliced_ethcd.mzML");

    let scans = mzlib::readers::read_spectra(&path).expect("mzML should read");
    let recorded = recording("readers_spectra_mzml.json");
    let source = scans.source.expect("an mzML records its source");

    assert_eq!(
        source.instrument_model.as_deref(),
        recorded["source"]["instrument_model"].as_str()
    );
    assert_eq!(
        source.instrument_serial_number.as_deref(),
        recorded["source"]["instrument_serial_number"].as_str()
    );
    assert!(matches!(
        source.acquired_at(),
        Some(mzlib::readers::AcquisitionTime::Utc(_))
    ));
}

#[test]
fn many_spectra_files_read_in_one_call_match_the_recording() {
    let Some(()) = require_verb("readers read-spectra") else {
        return;
    };
    let ethcd = fixture_or_skip!("DataFiles/sliced_ethcd.mzML");
    let mgf = fixture_or_skip!("DataFiles/withZeros.mgf");
    let missing = ethcd.with_file_name("no-such-run.mzML");

    let batch = mzlib::readers::read_spectra_many(
        &[&ethcd, &missing, &mgf],
        &mzlib::readers::SpectraBulkOptions {
            bulk: mzlib::readers::BulkOptions {
                threads: 2,
                on_error: mzlib::OnError::Skip,
                ..Default::default()
            },
            ..Default::default()
        },
    )
    .expect("a skip batch reads what it can");
    let recorded = recording("readers_many_spectra.json");

    assert_eq!(batch.file_count, 3);
    assert_eq!(batch.read_count, recorded["read_count"].as_u64().unwrap());
    assert_eq!(
        batch.record_count,
        recorded["record_count"].as_u64().unwrap()
    );
    assert_eq!(batch.columns.names()[..2], ["source_index", "source_path"]);
    assert_eq!(batch.failed_files().len(), 1);
    assert_eq!(batch.files[1].error.as_ref().unwrap().kind, "usage");
}

#[test]
fn a_batch_is_identical_at_any_thread_count() {
    // BULK.md section 2: byte-identical output at any --threads is a tested property of the
    // bridge. Asserted here too, through the binding, because a binding that reordered rows would
    // break it without the bridge knowing.
    let Some(()) = require_verb("readers read-records") else {
        return;
    };
    let plain = fixture_or_skip!("DataFiles/PXD078927_msgf_1_1_0.mzid");
    let gz = fixture_or_skip!("DataFiles/PXD078927_msgf_1_1_0.mzid.gz");

    let read = |threads| {
        mzlib::readers::read_records_many(
            &[&plain, &gz],
            &mzlib::readers::BulkOptions {
                threads,
                ..Default::default()
            },
        )
        .expect("the same record type twice reads")
    };
    let one = read(1);
    let four = read(4);
    assert_eq!(one.columns, four.columns);
    assert_eq!(one.files, four.files);
    assert_eq!(one.record_count, 24);
}

#[test]
fn identify_many_answers_each_path_in_order() {
    let Some(()) = require_verb("readers identify") else {
        return;
    };
    let mzid = fixture_or_skip!("DataFiles/PXD078927_msgf_1_1_0.mzid");
    let groups = fixture_or_skip!(
        "FileReadingTests/ExternalFileTypes/MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"
    );
    let missing = mzid.with_file_name("no-such-run.mzML");

    let batch = mzlib::readers::identify_many(
        &[&mzid, &missing, &groups],
        &mzlib::readers::BulkOptions {
            on_error: mzlib::OnError::Skip,
            ..Default::default()
        },
    )
    .expect("a skip batch identifies what it can");

    let types: Vec<&str> = batch.files.iter().map(|f| f.file_type.as_str()).collect();
    assert_eq!(
        types,
        ["MzIdentML", "", "MetaMorpheusQuantifiedProteinGroups"]
    );
    assert_eq!(batch.failed_count, 1);
}

#[test]
fn mzidentml_scores_come_back_long_as_recorded() {
    let Some(()) = require_verb("readers read-matches") else {
        return;
    };
    let path = fixture_or_skip!("DataFiles/PXD078927_msgf_1_1_0.mzid");

    let scored = mzlib::readers::read_matches_with(
        &path,
        &mzlib::readers::MatchOptions {
            read: mzlib::readers::ReadOptions {
                limit: Some(1),
                ..Default::default()
            },
            scores: true,
        },
    )
    .expect("an mzIdentML reads with scores");
    let recorded = recording("readers_matches_mzid_scores.json");

    assert_eq!(scored.row_count, recorded["row_count"].as_u64().unwrap());
    assert_eq!(
        scored.columns.strings("score_name").unwrap(),
        serde_json::from_value::<Vec<Option<String>>>(recorded["columns"]["score_name"].clone())
            .unwrap()
    );
    assert_eq!(scored.skipped_count, Some(0));
}

#[test]
fn the_quantification_tables_match_their_recordings() {
    let Some(()) = require_verb("readers read-protein-groups") else {
        return;
    };
    let groups_path = fixture_or_skip!(
        "FileReadingTests/ExternalFileTypes/MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"
    );
    let peptides_path = fixture_or_skip!(
        "FileReadingTests/ExternalFileTypes/MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"
    );
    let one = mzlib::readers::ReadOptions {
        limit: Some(1),
        ..Default::default()
    };

    let groups = mzlib::readers::read_protein_groups_with(&groups_path, &one)
        .expect("a protein-group table reads");
    let recorded = recording("readers_protein_groups.json");
    assert_eq!(groups.row_count, recorded["row_count"].as_u64().unwrap());
    assert_eq!(
        groups.columns.floats("intensity").unwrap(),
        serde_json::from_value::<Vec<Option<f64>>>(recorded["columns"]["intensity"].clone())
            .unwrap()
    );

    let peptides = mzlib::readers::read_quantified_peptides_with(&peptides_path, &one)
        .expect("a peptide table reads");
    assert_eq!(peptides.absent_fields, ["peak_order", "retention_time"]);

    let sites = mzlib::readers::read_occupancy(&groups_path).expect("occupancy reads");
    let recorded = recording("readers_occupancy.json");
    assert_eq!(sites.row_count, recorded["row_count"].as_u64().unwrap());
    assert_eq!(sites.truncated_cell_count, 0);
}
