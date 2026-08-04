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

use support::require_bridge;

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
        29,
        "mzLib recognises 29 file types; a change here means the crate's documented count is stale"
    );

    // 13 of 29 belong to no cross-format family, which is the fact that makes read_records
    // necessary rather than a convenience.
    let viewless = formats.iter().filter(|f| f.views.is_empty()).count();
    assert_eq!(viewless, 13, "13 of 29 have no view at all");

    let quantifiable = formats.iter().filter(|f| f.is_quantifiable()).count();
    assert_eq!(quantifiable, 3, "exactly 3 offer the quantifiable view");
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
fn a_fabricated_zero_intensity_crosses_as_none() {
    // A within-type schema variant: Apex_intensity is optional and the FLASHDeconv/OpenMS
    // _ms1.feature layout omits it, so mzLib substitutes zero for every feature. A whole column of
    // zeros is indistinguishable from real measurements of nothing, which is exactly what Option
    // exists to prevent — and the reason a binding must not paper over it with 0.0.
    let Some(()) = require_bridge() else { return };
    let path = fixture_or_skip!(
        "FileReadingTests/ExternalFileTypes/Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature"
    );

    let features = mzlib::readers::read_features(&path).expect("the feature view should read");

    let intensities = features
        .columns
        .floats("intensity")
        .expect("a numeric column");
    assert!(
        intensities.iter().all(Option::is_none),
        "a fabricated zero must not be handed back as a measurement"
    );
    assert!(
        features
            .caveats
            .iter()
            .any(|caveat| caveat.contains("intensity is NULL")),
        "and the reason must be stated, not left as an unexplained empty column"
    );
}

#[test]
fn a_topfd_feature_file_still_reports_real_intensities() {
    // The counterpart: nulling every intensity would be an equally serious over-correction.
    let Some(()) = require_bridge() else { return };
    let path =
        fixture_or_skip!("FileReadingTests/ExternalFileTypes/Ms1Feature_TopFDv1.6.2_ms1.feature");

    let features = mzlib::readers::read_features(&path).expect("the feature view should read");

    assert!(features
        .columns
        .floats("intensity")
        .expect("a numeric column")
        .iter()
        .all(Option::is_some));
}
