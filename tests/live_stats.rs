//! The stats verbs through the real bridge, over the committed reference tables.
//!
//! The offline suite replays recorded payloads, so it cannot see the bridge and this crate drifting
//! apart. These run the real executable on the files in `tests/fixtures/stats/` and check two
//! things: that what it returns today is what was recorded (so the docs' printed numbers are still
//! true), and that it still agrees with limma and metafor, through Rust, to 1e-8. This crate does
//! none of the arithmetic; these tests only compare mzLib's answers with R's. Local files only, no
//! network. pyMzLib's `tests/test_stats_live.py`, test for test.
//!
//! Run with `cargo test --features live`; skips on a bridge without the stats verbs.

#![cfg(feature = "live")]

mod support;

use std::path::PathBuf;

use mzlib::stats::{adjust, fit_with, meta, FitOptions, Study};
use serde_json::Value;
use support::require_verb;

fn stats_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("stats")
        .join(name)
}

/// A reference TSV as (header, rows).
fn tsv(name: &str) -> Vec<std::collections::HashMap<String, String>> {
    let text = std::fs::read_to_string(stats_file(name)).unwrap();
    let mut lines = text.lines().filter(|l| !l.is_empty());
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

fn worst_relative(ours: &[Option<f64>], theirs: &[f64]) -> f64 {
    assert_eq!(ours.len(), theirs.len());
    ours.iter()
        .zip(theirs)
        .map(|(a, b)| (a.expect("a fitted value") - b).abs() / b.abs().max(1e-300))
        .fold(0.0, f64::max)
}

fn recording(name: &str) -> Value {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name);
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// Equal to the recording, except that floats may differ in their last bits: the recordings were
/// made on one platform and replayed on others. 1e-9 relative is far below the 1e-8 the reference
/// comparisons promise, and far above floating-point noise.
fn assert_same_columns(live: &mzlib::readers::Table, recorded: &Value) {
    let recorded_names: Vec<&str> = recorded["column_names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    assert_eq!(live.names(), recorded_names);
    for name in recorded_names {
        let ours = live.raw(name).unwrap();
        let theirs = recorded["columns"][name].as_array().unwrap();
        assert_eq!(ours.len(), theirs.len(), "{name}");
        for (i, (a, b)) in ours.iter().zip(theirs).enumerate() {
            match (a.as_f64(), b.as_f64()) {
                (Some(a), Some(b)) => {
                    assert!(
                        (a - b).abs() <= 1e-9 * b.abs().max(1e-300),
                        "{name}[{i}]: {a} {b}"
                    );
                }
                _ => assert_eq!(a, b, "{name}[{i}]"),
            }
        }
    }
}

fn limma_fit() -> mzlib::stats::ModeratedFit {
    fit_with(
        stats_file("limma_reference_responses.tsv"),
        stats_file("limma_reference_design.tsv"),
        &["age_decades"],
        &FitOptions::default(),
    )
    .expect("the limma reference fits")
}

fn metafor_studies() -> Vec<Study> {
    tsv("metafor_dl_inputs.tsv")
        .iter()
        .map(|r| {
            Study::new(
                r["case"].clone(),
                r["yi"].parse().unwrap(),
                r["sei"].parse().unwrap(),
            )
        })
        .collect()
}

#[test]
fn fit_matches_limma() {
    let Some(()) = require_verb("stats fit") else {
        return;
    };
    let fit = limma_fit();
    let limma = tsv("limma_ebayes_notrend.tsv");
    assert_eq!(fit.feature_count as usize, limma.len());
    for (ours, theirs) in [
        ("t", "t_age"),
        ("p_value", "p_age"),
        ("bh_adjusted", "bh_age"),
        ("posterior_variance", "s2_post"),
        ("df_total", "df_total"),
        ("prior_variance", "s2_prior"),
    ] {
        let reference: Vec<f64> = limma.iter().map(|r| r[theirs].parse().unwrap()).collect();
        let worst = worst_relative(&fit.columns.floats(ours).unwrap(), &reference);
        eprintln!("stats fit {ours} vs limma {theirs}: worst relative difference {worst:e}");
        assert!(worst < 1e-8, "{ours}: {worst}");
    }
}

#[test]
fn meta_matches_metafor() {
    let Some(()) = require_verb("stats meta") else {
        return;
    };
    let pooled = meta(&metafor_studies()).expect("metafor's reference cases pool");
    let results = tsv("metafor_dl_results.tsv");
    let features = pooled.columns.strings("feature").unwrap();
    let mut worst = 0.0_f64;
    for (ours, theirs) in [
        ("estimate", "estimate"),
        ("standard_error", "se"),
        ("p_value", "pval"),
        ("confidence_low", "ci_lb"),
        ("confidence_high", "ci_ub"),
    ] {
        let column = pooled.columns.floats(ours).unwrap();
        for (i, feature) in features.iter().enumerate() {
            let expected = results
                .iter()
                .find(|r| Some(&r["case"]) == feature.as_ref())
                .unwrap();
            let difference = worst_relative(&[column[i]], &[expected[theirs].parse().unwrap()]);
            assert!(difference < 1e-8, "{feature:?} {ours}: {difference}");
            worst = worst.max(difference);
        }
    }
    eprintln!("stats meta vs metafor: worst relative difference {worst:e}");
}

#[test]
fn the_fit_recordings_are_what_the_bridge_returns_today() {
    let Some(()) = require_verb("stats fit") else {
        return;
    };
    let malat = fit_with(
        stats_file("malat_dilution_log2_vs_fluc.tsv"),
        stats_file("malat_dilution_design.tsv"),
        &["malat_250ng", "malat_125ng"],
        &FitOptions::default(),
    )
    .expect("the MALAT dilution fits");
    let recorded = recording("stats_fit_malat.json");
    assert_same_columns(&malat.columns, &recorded);
    assert_eq!(
        serde_json::to_value(&malat.caveats).unwrap(),
        recorded["caveats"]
    );

    let limma = limma_fit();
    let recorded = recording("stats_fit_limma.json");
    assert_same_columns(&limma.columns, &recorded);
    assert_eq!(
        serde_json::to_value(&limma.caveats).unwrap(),
        recorded["caveats"]
    );
}

#[test]
fn the_adjust_and_meta_recordings_are_what_the_bridge_returns_today() {
    let Some(()) = require_verb("stats adjust") else {
        return;
    };
    let adjusted = adjust(&[
        Some(0.0002),
        Some(0.004),
        Some(0.019),
        None,
        Some(0.031),
        Some(0.2),
        Some(f64::NAN),
        Some(0.74),
    ])
    .expect("the p-values adjust");
    let recorded = recording("stats_adjust.json");
    assert_same_columns(&adjusted.columns, &recorded);
    assert_eq!(
        serde_json::to_value(&adjusted.caveats).unwrap(),
        recorded["caveats"]
    );

    let pooled = meta(&metafor_studies()).expect("the studies pool");
    assert_same_columns(&pooled.columns, &recording("stats_meta_metafor.json"));
}

#[test]
fn a_trend_fit_runs_and_reports_its_basis() {
    let Some(()) = require_verb("stats fit") else {
        return;
    };
    let fit = fit_with(
        stats_file("limma_reference_responses.tsv"),
        stats_file("limma_reference_design.tsv"),
        &["age_decades"],
        &FitOptions {
            trend: true,
            ..Default::default()
        },
    )
    .expect("a trended fit runs");
    assert!(fit.prior.trended);
    assert_eq!(fit.prior.scale, None);
    assert!((1..=4).contains(&fit.prior.spline_basis_count));
}

#[test]
fn bad_input_is_a_usage_error_from_the_bridge() {
    let Some(()) = require_verb("stats fit") else {
        return;
    };
    let dir = std::env::temp_dir().join(format!("mzlib-live-stats-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let design = dir.join("design.tsv");
    std::fs::write(&design, "sample\tintercept\na\t1\nb\t1\n").unwrap();
    let responses = dir.join("responses.tsv");
    std::fs::write(&responses, "id\ta\tc\nf1\t1\t2\n").unwrap();

    let error = fit_with(&responses, &design, &["intercept"], &FitOptions::default()).unwrap_err();
    assert!(
        matches!(&error, mzlib::MzLibError::Usage(m) if m.contains("Not in the design: c")),
        "{error}"
    );
    let error = adjust(&[Some(0.5), Some(2.0)]).unwrap_err();
    assert!(
        matches!(&error, mzlib::MzLibError::Usage(m) if m.contains("outside [0, 1]")),
        "{error}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
