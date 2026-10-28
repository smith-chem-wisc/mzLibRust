//! The diff verbs through the real bridge, over the committed sample sheets and peptide tables.
//!
//! The offline suite replays recorded payloads, so it cannot see the bridge and this crate drift
//! apart. These run the real executable on the files in `tests/fixtures/diff/` and check that what
//! it returns today is what was recorded (so the docs' printed numbers are still true), and that
//! the files `diff design` writes are what `stats fit` reads. Local files only, no network.
//! pyMzLib's `tests/test_diff_live.py`, test for test.
//!
//! Run with `cargo test --features live`; skips on a bridge without the diff verbs.

#![cfg(feature = "live")]

mod support;

use std::path::{Path, PathBuf};

use mzlib::diff::{
    design_with, observations, observations_with, CovariateSpec, Design, DesignOptions,
    DesignStratum, FactorSpec, ObservationRun, Observations, ObservationsOptions, QuantBasis,
};
use mzlib::stats::fit;
use support::require_verb;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
}

fn diff_file(name: &str) -> PathBuf {
    fixtures().join("diff").join(name)
}

fn stats_file(name: &str) -> PathBuf {
    fixtures().join("stats").join(name)
}

fn recorded<T: serde::de::DeserializeOwned>(name: &str) -> T {
    serde_json::from_str(&std::fs::read_to_string(fixtures().join(name)).unwrap()).unwrap()
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("mzlib-live-diff-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The last part of a path written on either platform: a recording made on Windows is compared on
/// Linux, so a platform path splitter is never applied to a recorded path.
fn file_name(path: &str) -> String {
    path.rsplit(['/', '\\']).next().unwrap_or(path).to_owned()
}

fn without_paths(mut stratum: DesignStratum) -> DesignStratum {
    stratum.design_file = stratum.design_file.as_deref().map(file_name);
    stratum.contrasts_file = stratum.contrasts_file.as_deref().map(file_name);
    stratum
}

fn assert_same_design(live: &Design, recorded: &Design, name: &str) {
    let strata = |d: &Design| -> Vec<DesignStratum> {
        d.strata.iter().cloned().map(without_paths).collect()
    };
    assert_eq!(strata(live), strata(recorded), "{name}");
    assert_eq!(live.columns, recorded.columns, "{name}");
    assert_eq!(live.caveats, recorded.caveats, "{name}");
    assert_eq!(
        live.estimability_counts, recorded.estimability_counts,
        "{name}"
    );
    assert_eq!(live.unused_columns, recorded.unused_columns, "{name}");
}

fn mm_runs() -> Vec<ObservationRun> {
    let header =
        std::fs::read_to_string(diff_file("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"))
            .unwrap();
    header
        .lines()
        .next()
        .unwrap()
        .split('\t')
        .filter_map(|c| c.strip_prefix("Intensity_"))
        .map(|f| f.trim_end_matches("-calib"))
        .map(|f| ObservationRun::new(f, f))
        .collect()
}

fn spike_in() -> Vec<ObservationRun> {
    vec![
        ObservationRun::new("09-04-18_EcoliSpikeInSingleShot1x", "spike_1x"),
        ObservationRun::new("09-04-18_EcoliSpikeInSingleShot2x", "spike_2x"),
    ]
}

#[test]
fn the_design_recordings_are_still_what_the_bridge_says() {
    let Some(()) = require_verb("diff design") else {
        return;
    };
    let dir = scratch("design");
    let cases = [
        (
            "diff_design_malat.json",
            "malat_dilution_samples.tsv",
            DesignOptions {
                factors: vec![FactorSpec::new("malat_load", Some("500ng"))],
                out: Some(dir.join("malat_design")),
                ..Default::default()
            },
        ),
        (
            "diff_design_factor_batch_covariate.json",
            "design_factor_batch_covariate_samples.tsv",
            DesignOptions {
                factors: vec![FactorSpec::new("condition", Some("ctrl"))],
                covariates: vec![CovariateSpec::new("age", "years")
                    .scale(10.0)
                    .scaled_unit("decade")
                    .centre(50.0)],
                batch: Some("batch".into()),
                ..Default::default()
            },
        ),
        (
            "diff_design_confounded.json",
            "design_confounded_samples.tsv",
            DesignOptions {
                factors: vec![FactorSpec::new("condition", Some("A"))],
                batch: Some("batch".into()),
                ..Default::default()
            },
        ),
        (
            "diff_design_strata.json",
            "stratified_samples.tsv",
            DesignOptions {
                factors: vec![FactorSpec::new("condition", Some("A"))],
                strata: vec!["organism part".into()],
                ..Default::default()
            },
        ),
    ];
    for (name, samples, options) in cases {
        let live = design_with(diff_file(samples), &options).unwrap();
        assert_same_design(&live, &recorded(name), name);
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_observation_recordings_are_still_what_the_bridge_says() {
    let Some(()) = require_verb("diff observations") else {
        return;
    };
    let dir = scratch("observations");
    let groups = Some(diff_file(
        "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv",
    ));
    let mm = "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv";
    let spike = "EcoliSpikeIn_QuantifiedPeptides.tsv";
    let cases = [
        (
            "diff_observations_mm_mbr_kept.json",
            mm,
            mm_runs(),
            QuantBasis::MbrKept,
            ObservationsOptions {
                protein_groups: groups.clone(),
                ..Default::default()
            },
        ),
        (
            "diff_observations_mm_msms_only.json",
            mm,
            mm_runs(),
            QuantBasis::MsmsOnly,
            ObservationsOptions {
                protein_groups: groups.clone(),
                ..Default::default()
            },
        ),
        (
            "diff_observations_spike_in.json",
            spike,
            spike_in(),
            QuantBasis::MsmsOnly,
            ObservationsOptions::default(),
        ),
        (
            "diff_observations_spike_in_normalized.json",
            spike,
            spike_in(),
            QuantBasis::MsmsOnly,
            ObservationsOptions {
                min_peptides: Some(50),
                normalize: true,
                responses_out: Some(dir.join("spike_in_responses.tsv")),
                ..Default::default()
            },
        ),
    ];
    for (name, peptides, runs, basis, options) in cases {
        let live = observations_with(diff_file(peptides), &runs, basis, &options).unwrap();
        let expected: Observations = recorded(name);
        assert_eq!(live.columns, expected.columns, "{name}");
        assert_eq!(live.normalization, expected.normalization, "{name}");
        assert_eq!(live.state_counts, expected.state_counts, "{name}");
        assert_eq!(live.caveats, expected.caveats, "{name}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn design_then_fit_is_the_committed_fit() {
    // The design file diff design writes is the same model as the committed 0/1 design.
    let Some(()) = require_verb("diff design") else {
        return;
    };
    let dir = scratch("design-fit");
    let design = design_with(
        diff_file("malat_dilution_samples.tsv"),
        &DesignOptions {
            factors: vec![FactorSpec::new("malat_load", Some("500ng"))],
            out: Some(dir.clone()),
            ..Default::default()
        },
    )
    .unwrap();
    let ours = fit(
        stats_file("malat_dilution_log2_vs_fluc.tsv"),
        design.strata[0].design_file.as_deref().unwrap(),
        &["malat_load=250ng", "malat_load=125ng"],
    )
    .unwrap();
    let committed = fit(
        stats_file("malat_dilution_log2_vs_fluc.tsv"),
        stats_file("malat_dilution_design.tsv"),
        &["malat_250ng", "malat_125ng"],
    )
    .unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    assert_eq!(ours.status_counts, committed.status_counts);
    let (a, b) = (
        ours.columns.floats("estimate").unwrap(),
        committed.columns.floats("estimate").unwrap(),
    );
    for (a, b) in a.iter().zip(&b) {
        match (a, b) {
            (Some(a), Some(b)) => assert!((a - b).abs() <= 1e-9 * b.abs() + 1e-12, "{a} {b}"),
            _ => assert_eq!(a, b),
        }
    }
}

#[test]
fn the_responses_table_is_what_stats_fit_reads() {
    let Some(()) = require_verb("diff observations") else {
        return;
    };
    let dir = scratch("responses");
    let responses = dir.join("responses.tsv");
    let table = observations_with(
        diff_file("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"),
        &mm_runs(),
        QuantBasis::MbrKept,
        &ObservationsOptions {
            responses_out: Some(responses.clone()),
            ..Default::default()
        },
    )
    .unwrap();
    let design = dir.join("design.tsv");
    let rows: String = table.samples.iter().map(|s| format!("{s}\t1\n")).collect();
    std::fs::write(&design, format!("sample\tintercept\n{rows}")).unwrap();
    let fitted = fit(&responses, &design, &["intercept"]).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    assert_eq!(fitted.feature_count, table.peptide_count);
    assert_eq!(fitted.sample_names, table.samples);
    assert_eq!(
        table.responses_written.map(|w| file_name(&w.path)),
        Some("responses.tsv".to_owned())
    );
}

#[test]
fn a_refusal_crosses_as_a_usage_error() {
    let Some(()) = require_verb("diff design") else {
        return;
    };
    let dir = scratch("refusal");
    let sheet = dir.join("sheet.tsv");
    std::fs::write(&sheet, "sample\tcondition\ns1\tA\ns2\t\n").unwrap();
    let usage = |error: mzlib::MzLibError, needle: &str| {
        assert!(
            matches!(&error, mzlib::MzLibError::Usage(m) if m.contains(needle)),
            "{needle}: {error}"
        );
    };
    usage(
        design_with(
            &sheet,
            &DesignOptions {
                factors: vec![FactorSpec::new("condition", Some("A"))],
                ..Default::default()
            },
        )
        .unwrap_err(),
        "has no value for factor 'condition'",
    );
    usage(
        observations(
            diff_file("EcoliSpikeIn_QuantifiedPeptides.tsv"),
            &spike_in()[..1],
            QuantBasis::MsmsOnly,
        )
        .unwrap_err(),
        "not in the experimental design",
    );
    let _ = std::fs::remove_dir_all(Path::new(&dir));
}
