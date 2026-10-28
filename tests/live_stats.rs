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
fn a_small_trend_fit_sizes_its_spline_as_limma_does() {
    // Without spline_basis nothing is sent, and mzLib follows limma's fitFDist rule (mzLib #1398):
    // 1 + (n >= 3) + (n >= 6) + (n >= 30) basis functions for n fitted features. Ten features
    // therefore get 3, where the bridge used to ask for 4 on every trend fit.
    let Some(()) = require_verb("stats fit") else {
        return;
    };
    let text = std::fs::read_to_string(stats_file("limma_reference_responses.tsv")).unwrap();
    let mut lines = text.lines().filter(|l| !l.is_empty());
    let header = lines.next().unwrap();
    let complete: Vec<&str> = lines
        .filter(|line| {
            !line
                .split('\t')
                .skip(1)
                .any(|cell| cell == "NaN" || cell.is_empty())
        })
        .take(10)
        .collect();
    assert_eq!(
        complete.len(),
        10,
        "the reference table has ten complete features"
    );

    let dir = std::env::temp_dir().join(format!("mzlib-live-spline-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let responses = dir.join("ten_features.tsv");
    std::fs::write(&responses, format!("{header}\n{}\n", complete.join("\n"))).unwrap();

    let fit = fit_with(
        &responses,
        stats_file("limma_reference_design.tsv"),
        &["age_decades"],
        &FitOptions {
            trend: true,
            ..Default::default()
        },
    )
    .expect("a small trended fit runs");
    let _ = std::fs::remove_dir_all(&dir);

    assert_eq!(fit.feature_count, 10);
    assert_eq!(
        fit.prior.spline_basis_count, 3,
        "limma's rule gives 3 basis functions for 6 to 29 features"
    );
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

// ---- pyMzLib 0.6.0: contrasts, the prior estimator, and the five new stats verbs ----------------
//
// pyMzLib's test_stats_live.py (contrasts, estimator), test_stats_combine_correlate_live.py and
// test_stats_models_live.py, test for test. Each skips on a bridge without the verb.

use mzlib::stats::{
    combine_p, correlate, correlate_with, entrapment_fdp, entrapment_fdp_with, fit_logistic,
    fit_logistic_with, fit_mixed, fit_mixed_with, CorrelateOptions, EntrapmentFdpOptions,
    FitLogisticOptions, FitMixedOptions, PValueTest,
};

/// Closed-form and iteratively reweighted fits agree across platforms to 1e-9. The mixed model's
/// variance ratio and the marginal-likelihood prior are minima of a flat likelihood, found by a
/// search whose last bits follow the platform's math library: 1e-6 there, as pyMzLib holds them.
const CROSS_PLATFORM: f64 = 1e-9;
const CROSS_PLATFORM_OPTIMISED: f64 = 1e-6;

fn assert_columns_within(live: &mzlib::readers::Table, recorded: &Value, relative: f64) {
    let names: Vec<&str> = recorded["column_names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    assert_eq!(live.names(), names);
    let close = |a: &Value, b: &Value, at: &str| match (a.as_f64(), b.as_f64()) {
        (Some(a), Some(b)) => assert!(
            (a - b).abs() <= relative * b.abs().max(1e-12),
            "{at}: {a} {b}"
        ),
        _ => assert_eq!(a, b, "{at}"),
    };
    for name in names {
        let ours = live.raw(name).unwrap();
        let theirs = recorded["columns"][name].as_array().unwrap();
        assert_eq!(ours.len(), theirs.len(), "{name}");
        for (i, (a, b)) in ours.iter().zip(theirs).enumerate() {
            let at = format!("{name}[{i}]");
            match (a, b) {
                (Value::Array(a), Value::Array(b)) => {
                    assert_eq!(a.len(), b.len(), "{at}");
                    for (x, y) in a.iter().zip(b) {
                        close(x, y, &at);
                    }
                }
                _ => close(a, b, &at),
            }
        }
    }
}

fn cohorts() -> Vec<PValueTest> {
    [
        ("up_everywhere", 0.003, 4.0, 0.8),
        ("up_everywhere", 0.02, 5.0, 0.5),
        ("up_everywhere", 0.3, 3.0, 0.3),
        ("flips_sign", 0.001, 4.0, 1.2),
        ("flips_sign", 0.004, 5.0, -0.9),
        ("flips_sign", 0.4, 3.0, 0.2),
        ("no_effect", 0.6, 4.0, -0.1),
        ("no_effect", 0.45, 5.0, 0.1),
        ("no_effect", 0.85, 3.0, -0.05),
    ]
    .iter()
    .map(|&(f, p, w, e)| PValueTest::new(f, p).weight(w).effect(e))
    .collect()
}

fn combination_tests(weighted: bool) -> Vec<PValueTest> {
    tsv("combination_inputs.tsv")
        .iter()
        .map(|r| {
            let test = PValueTest::new(r["case"].as_str(), r["p"].parse::<f64>().unwrap());
            if weighted {
                test.weight(r["weight"].parse().unwrap())
            } else {
                test
            }
        })
        .collect()
}

/// `diff design` with no reference level, written to a scratch directory: every pair of loads.
fn pairs_design(dir: &std::path::Path) -> mzlib::diff::DesignStratum {
    let design = mzlib::diff::design_with(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("diff")
            .join("malat_dilution_samples.tsv"),
        &mzlib::diff::DesignOptions {
            factors: vec![mzlib::diff::FactorSpec::new("malat_load", None)],
            out: Some(dir.join("malat_pairs")),
            ..Default::default()
        },
    )
    .expect("the MALAT1 sheet designs");
    design.strata[0].clone()
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("mzlib-live-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn malat_fit(options: &FitOptions) -> mzlib::stats::ModeratedFit {
    fit_with(
        stats_file("malat_dilution_log2_vs_fluc.tsv"),
        stats_file("malat_dilution_design.tsv"),
        &["malat_250ng", "malat_125ng"],
        options,
    )
    .expect("the MALAT dilution fits")
}

fn usage_containing<T: std::fmt::Debug>(result: mzlib::Result<T>, needle: &str) {
    match result {
        Err(mzlib::MzLibError::Usage(message)) => {
            assert!(message.contains(needle), "{needle}: {message}");
        }
        other => panic!("expected a usage error containing {needle:?}, got {other:?}"),
    }
}

#[test]
fn the_marginal_recording_is_what_the_bridge_returns_today() {
    // A bridge from before pyMzLib 0.6.0 ignores --estimator; it has no stats combine-p either.
    let Some(()) = require_verb("stats combine-p") else {
        return;
    };
    let live = malat_fit(&FitOptions {
        estimator: Some("marginal_likelihood".into()),
        ..Default::default()
    });
    let recorded = recording("stats_fit_malat_marginal.json");
    assert_columns_within(&live.columns, &recorded, CROSS_PLATFORM_OPTIMISED);
    assert_eq!(
        serde_json::to_value(&live.caveats).unwrap(),
        recorded["caveats"]
    );
    assert_eq!(
        serde_json::to_value(&live.status_counts).unwrap(),
        recorded["status_counts"]
    );
    assert_eq!(live.prior.estimator.as_deref(), Some("marginal_likelihood"));
}

#[test]
fn the_contrasts_recording_is_what_the_chain_returns_today() {
    let Some(()) = require_verb("diff design") else {
        return;
    };
    let dir = scratch("contrasts");
    let stratum = pairs_design(&dir);
    let none: [&str; 0] = [];
    let live = fit_with(
        stats_file("malat_dilution_log2_vs_fluc.tsv"),
        stratum.design_file.as_deref().unwrap(),
        &none,
        &FitOptions {
            contrasts: stratum.contrasts_file.as_deref().map(Into::into),
            ..Default::default()
        },
    )
    .expect("diff design's files fit");
    let _ = std::fs::remove_dir_all(&dir);
    let recorded = recording("stats_fit_malat_contrasts.json");
    let contrasts = live.contrasts.as_ref().unwrap();
    let theirs = &recorded["contrasts"];
    assert_eq!(
        serde_json::to_value(&contrasts.tested_contrasts).unwrap(),
        theirs["tested_contrasts"]
    );
    assert_eq!(
        Some(contrasts.confidence_level),
        theirs["confidence_level"].as_f64()
    );
    for (id, weights) in theirs["weights"].as_object().unwrap() {
        let recorded: Vec<f64> = weights
            .as_array()
            .unwrap()
            .iter()
            .map(|w| w.as_f64().unwrap())
            .collect();
        assert_eq!(contrasts.weights[id], recorded, "{id}");
    }
    assert_columns_within(&contrasts.columns, theirs, CROSS_PLATFORM);
    assert_eq!(
        serde_json::to_value(&live.caveats).unwrap(),
        recorded["caveats"]
    );
}

#[test]
fn a_contrast_of_one_coefficient_is_that_coefficients_test() {
    let Some(()) = require_verb("diff design") else {
        return;
    };
    let dir = scratch("one-contrast");
    let stratum = pairs_design(&dir);
    let fit = fit_with(
        stats_file("malat_dilution_log2_vs_fluc.tsv"),
        stratum.design_file.as_deref().unwrap(),
        &["malat_load=250ng"],
        &FitOptions {
            contrasts: stratum.contrasts_file.as_deref().map(Into::into),
            ..Default::default()
        },
    )
    .unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let contrasts = fit.contrasts.as_ref().unwrap();
    assert_eq!(contrasts.weights["c1"], [0.0, 1.0, 0.0]);
    let coefficient = fit.rows("malat_load=250ng").unwrap();
    let contrast = contrasts.rows("c1").unwrap();
    for (a, b) in coefficient.iter().zip(&contrast) {
        assert_eq!(
            (a.estimate, a.standard_error, a.t, a.p_value, a.bh_adjusted),
            (b.estimate, b.standard_error, b.t, b.p_value, b.bh_adjusted)
        );
    }
}

#[test]
fn a_lower_confidence_gives_narrower_intervals_and_mzlib_refuses_one_outside_0_1() {
    let Some(()) = require_verb("diff design") else {
        return;
    };
    let dir = scratch("confidence");
    let stratum = pairs_design(&dir);
    let design = stratum.design_file.clone().unwrap();
    let contrasts = stratum.contrasts_file.as_deref().map(PathBuf::from);
    let none: [&str; 0] = [];
    let run = |confidence: Option<f64>| {
        fit_with(
            stats_file("malat_dilution_log2_vs_fluc.tsv"),
            &design,
            &none,
            &FitOptions {
                contrasts: contrasts.clone(),
                confidence,
                ..Default::default()
            },
        )
    };
    let wide = run(None).unwrap().contrasts.unwrap();
    let narrow = run(Some(0.8)).unwrap().contrasts.unwrap();
    assert_eq!(narrow.confidence_level, 0.8);
    for (w, n) in wide
        .rows("c3")
        .unwrap()
        .iter()
        .zip(narrow.rows("c3").unwrap())
    {
        if w.status == "fitted" {
            assert!(w.confidence_low < n.confidence_low);
            assert!(n.confidence_low < n.confidence_high);
            assert!(n.confidence_high < w.confidence_high);
        }
    }
    usage_containing(run(Some(1.5)), "strictly between 0 and 1");
    usage_containing(
        fit_with(
            stats_file("malat_dilution_log2_vs_fluc.tsv"),
            &design,
            &["malat_load=250ng"],
            &FitOptions {
                confidence: Some(0.8),
                ..Default::default()
            },
        ),
        "give --contrasts as well",
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_default_estimator_is_the_one_mzlib_reports() {
    let Some(()) = require_verb("stats combine-p") else {
        return;
    };
    let default = malat_fit(&FitOptions::default());
    let legacy = malat_fit(&FitOptions {
        estimator: Some("moments_legacy".into()),
        ..Default::default()
    });
    assert_eq!(default.prior.estimator.as_deref(), Some("moments_legacy"));
    assert_eq!(default.columns.raw("t"), legacy.columns.raw("t"));
}

#[test]
fn a_fit_of_fewer_than_three_features_fits_no_prior() {
    // mzLib#1398, as limma's squeezeVar: the prior df is 0, and no feature's own variance is
    // reported as the prior. Bridges before pyMzLib 0.6.0 put the first feature's variance in scale.
    let Some(()) = require_verb("stats combine-p") else {
        return;
    };
    let text = std::fs::read_to_string(stats_file("limma_reference_responses.tsv")).unwrap();
    let mut lines = text.lines().filter(|l| !l.is_empty());
    let header = lines.next().unwrap();
    let two: Vec<&str> = lines.take(2).collect();
    let dir = scratch("two-features");
    let responses = dir.join("two_features.tsv");
    std::fs::write(&responses, format!("{header}\n{}\n", two.join("\n"))).unwrap();
    let fit = fit_with(
        &responses,
        stats_file("limma_reference_design.tsv"),
        &["age_decades"],
        &FitOptions::default(),
    )
    .expect("two features fit");
    let _ = std::fs::remove_dir_all(&dir);
    assert_eq!((fit.prior.df, fit.prior.scale), (Some(0.0), None));
    let fitted: Vec<_> = fit
        .rows("age_decades")
        .unwrap()
        .into_iter()
        .filter(|r| r.status == "fitted")
        .collect();
    assert!(!fitted.is_empty());
    assert!(fitted.iter().all(|r| r.prior_variance.is_some()));
}

#[test]
fn combine_p_matches_r() {
    let Some(()) = require_verb("stats combine-p") else {
        return;
    };
    let results = tsv("combination_results.tsv");
    let fisher = combine_p(&combination_tests(false), "fisher").unwrap();
    let plain = combine_p(&combination_tests(false), "stouffer").unwrap();
    let weighted = combine_p(&combination_tests(true), "stouffer").unwrap();
    for (result, column) in [
        (&fisher, "fisher_p"),
        (&plain, "stouffer_p"),
        (&weighted, "weighted_p"),
    ] {
        let reference: Vec<f64> = results.iter().map(|r| r[column].parse().unwrap()).collect();
        let worst = worst_relative(&result.columns.floats("p_value").unwrap(), &reference);
        eprintln!("stats combine-p {column}: worst relative difference from R {worst:e}");
        assert!(worst < 1e-8, "{column}: {worst}");
    }
}

#[test]
fn the_combine_and_correlate_recordings_are_what_the_bridge_returns_today() {
    let Some(()) = require_verb("stats combine-p") else {
        return;
    };
    let unsigned: Vec<PValueTest> = cohorts()
        .into_iter()
        .map(|mut t| {
            t.effect = None;
            t
        })
        .collect();
    let combined = [
        (
            "stats_combine_p_fisher.json",
            combine_p(&combination_tests(false), "fisher"),
        ),
        (
            "stats_combine_p_weighted.json",
            combine_p(&combination_tests(true), "stouffer"),
        ),
        (
            "stats_combine_p_signed.json",
            combine_p(&cohorts(), "stouffer"),
        ),
        (
            "stats_combine_p_unsigned.json",
            combine_p(&unsigned, "stouffer"),
        ),
    ];
    for (fixture, live) in combined {
        let live = live.unwrap();
        let recorded = recording(fixture);
        assert_columns_within(&live.columns, &recorded, CROSS_PLATFORM);
        assert_eq!(
            serde_json::to_value(&live.caveats).unwrap(),
            recorded["caveats"],
            "{fixture}"
        );
    }
    let responses = stats_file("limma_reference_responses.tsv");
    let age = stats_file("limma_reference_age.tsv");
    let correlations = [
        (
            "stats_correlate_malat.json",
            correlate(
                stats_file("malat_dilution_log2_vs_fluc.tsv"),
                stats_file("malat_dilution_load_ng.tsv"),
            ),
        ),
        (
            "stats_correlate_limma_age.json",
            correlate(&responses, &age),
        ),
        (
            "stats_correlate_limma_age_given_sex.json",
            correlate_with(
                &responses,
                &age,
                &CorrelateOptions {
                    given: Some(stats_file("limma_reference_sex.tsv")),
                    ..Default::default()
                },
            ),
        ),
    ];
    for (fixture, live) in correlations {
        let live = live.unwrap();
        let recorded = recording(fixture);
        assert_columns_within(&live.columns, &recorded, CROSS_PLATFORM);
        assert_eq!(
            serde_json::to_value(&live.caveats).unwrap(),
            recorded["caveats"],
            "{fixture}"
        );
        // Its keys come from an mzLib enum; a repin that grows it must not slip past the recording.
        assert_eq!(
            serde_json::to_value(&live.method_counts).unwrap(),
            recorded["method_counts"],
            "{fixture}"
        );
    }
}

fn relative_to(ours: f64, theirs: &str) -> f64 {
    let theirs: f64 = theirs.parse().unwrap();
    (ours - theirs).abs() / theirs.abs().max(1e-300)
}

#[test]
fn fit_logistic_matches_glm() {
    let Some(()) = require_verb("stats fit-logistic") else {
        return;
    };
    let fit = fit_logistic(
        stats_file("glm_reference_responses.tsv"),
        stats_file("glm_reference_design.tsv"),
        &["Intercept", "x", "grp"],
    )
    .unwrap();
    let glm = tsv("glm_reference_fit.tsv");
    for (j, coefficient) in ["Intercept", "x", "grp"].iter().enumerate() {
        for (ours, theirs) in fit.rows(coefficient).unwrap().iter().zip(&glm) {
            if theirs["status"] == "separated" {
                assert_eq!(ours.status, "separated");
                continue;
            }
            assert_eq!(ours.status, "fitted");
            assert!(relative_to(ours.estimate.unwrap(), &theirs[&format!("b{j}")]) < 1e-8);
            assert!(relative_to(ours.standard_error.unwrap(), &theirs[&format!("se{j}")]) < 1e-6);
        }
    }
    for (ours, theirs) in fit.rows("x").unwrap().iter().zip(&glm) {
        if theirs["status"] == "fitted" {
            assert!(relative_to(ours.p_value.unwrap(), &theirs["p1"]) < 1e-6);
            assert!(relative_to(ours.deviance.unwrap(), &theirs["deviance"]) < 1e-8);
        }
    }
}

#[test]
fn fit_mixed_matches_nlme() {
    let Some(()) = require_verb("stats fit-mixed") else {
        return;
    };
    let reference = tsv("lmm_reference_fit.tsv");
    for (method, estimator) in [("REML", None), ("ML", Some("ml"))] {
        let fit = fit_mixed_with(
            stats_file("lmm_reference_responses.tsv"),
            stats_file("lmm_reference_design.tsv"),
            stats_file("lmm_reference_groups.tsv"),
            &["Intercept", "age", "wp"],
            &FitMixedOptions {
                estimator: estimator.map(str::to_owned),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(fit.estimator, method.to_lowercase());
        let nlme: Vec<_> = reference.iter().filter(|r| r["method"] == method).collect();
        for (j, coefficient) in ["Intercept", "age", "wp"].iter().enumerate() {
            for (ours, theirs) in fit.rows(coefficient).unwrap().iter().zip(&nlme) {
                assert!(ours.status == "fitted" && theirs["status"] == "fitted");
                assert!(relative_to(ours.estimate.unwrap(), &theirs[&format!("b{j}")]) < 1e-5);
                assert!(
                    relative_to(ours.standard_error.unwrap(), &theirs[&format!("se{j}")]) < 1e-5
                );
                assert_eq!(ours.df, theirs[&format!("df{j}")].parse::<f64>().ok());
            }
        }
        for (ours, theirs) in fit.rows("age").unwrap().iter().zip(&nlme) {
            assert!(relative_to(ours.log_likelihood.unwrap(), &theirs["loglik"]) < 1e-5);
            let sigma2: f64 = theirs["sigma2"].parse().unwrap();
            let tau2: f64 = theirs["tau2"].parse().unwrap();
            assert!((ours.group_variance.unwrap() - tau2).abs() < 1e-5 * sigma2);
        }
    }
}

#[test]
fn entrapment_fdp_matches_the_reference_script() {
    let Some(()) = require_verb("stats entrapment-fdp") else {
        return;
    };
    let reference = tsv("entrapment_reference_r1.tsv");
    let thresholds: Vec<f64> = reference
        .iter()
        .map(|r| r["q_value_threshold"].parse().unwrap())
        .collect();
    let fdp = entrapment_fdp_with(
        stats_file("entrapment_reference_targets.tsv"),
        stats_file("entrapment_reference_entrapments.tsv"),
        1.0,
        &EntrapmentFdpOptions {
            thresholds: Some(thresholds),
            ..Default::default()
        },
    )
    .unwrap();
    let rows = fdp.rows().unwrap();
    assert_eq!(rows.len(), reference.len());
    for (ours, theirs) in rows.iter().zip(&reference) {
        assert_eq!(ours.target_count.to_string(), theirs["target_count"]);
        assert_eq!(
            ours.entrapment_count.to_string(),
            theirs["entrapment_count"]
        );
        for (value, column) in [
            (Some(ours.lower_bound), "lower_bound"),
            (Some(ours.combined), "combined"),
            (ours.paired, "paired"),
        ] {
            let expected: f64 = theirs[column].parse().unwrap();
            assert!((value.unwrap() - expected).abs() < 1e-5, "{column}");
        }
    }
}

#[test]
fn the_model_recordings_are_what_the_bridge_returns_today() {
    let Some(()) = require_verb("stats fit-mixed") else {
        return;
    };
    let (glm_responses, glm_design) = (
        stats_file("glm_reference_responses.tsv"),
        stats_file("glm_reference_design.tsv"),
    );
    let logistic = [
        (
            "stats_fit_logistic_glm.json",
            fit_logistic(&glm_responses, &glm_design, &["Intercept", "x", "grp"]),
        ),
        (
            "stats_fit_logistic_malat.json",
            fit_logistic_with(
                stats_file("malat_dilution_detected.tsv"),
                stats_file("malat_dilution_design.tsv"),
                &["malat_250ng", "malat_125ng"],
                &FitLogisticOptions::default(),
            ),
        ),
    ];
    for (fixture, live) in logistic {
        let live = live.unwrap();
        let recorded = recording(fixture);
        assert_columns_within(&live.columns, &recorded, CROSS_PLATFORM);
        assert_eq!(
            serde_json::to_value(&live.caveats).unwrap(),
            recorded["caveats"]
        );
        // Its keys are the bridge's reading of mzLib's code (mzLib#1487); a change must not slip.
        assert_eq!(
            serde_json::to_value(&live.status_counts).unwrap(),
            recorded["status_counts"]
        );
    }
    let (lmm_responses, lmm_design, lmm_groups) = (
        stats_file("lmm_reference_responses.tsv"),
        stats_file("lmm_reference_design.tsv"),
        stats_file("lmm_reference_groups.tsv"),
    );
    let coefficients = ["Intercept", "age", "wp"];
    let mixed = [
        (
            "stats_fit_mixed_lmm_reml.json",
            fit_mixed(&lmm_responses, &lmm_design, &lmm_groups, &coefficients),
        ),
        (
            "stats_fit_mixed_lmm_ml.json",
            fit_mixed_with(
                &lmm_responses,
                &lmm_design,
                &lmm_groups,
                &coefficients,
                &FitMixedOptions {
                    estimator: Some("ml".into()),
                    ..Default::default()
                },
            ),
        ),
        (
            "stats_fit_mixed_malat.json",
            fit_mixed(
                stats_file("malat_dilution_log2_vs_fluc.tsv"),
                stats_file("malat_dilution_design.tsv"),
                stats_file("malat_dilution_day.tsv"),
                &["malat_250ng", "malat_125ng"],
            ),
        ),
    ];
    for (fixture, live) in mixed {
        let live = live.unwrap();
        let recorded = recording(fixture);
        assert_columns_within(&live.columns, &recorded, CROSS_PLATFORM_OPTIMISED);
        assert_eq!(
            serde_json::to_value(&live.caveats).unwrap(),
            recorded["caveats"]
        );
        assert_eq!(
            serde_json::to_value(&live.status_counts).unwrap(),
            recorded["status_counts"]
        );
    }
    let targets = stats_file("entrapment_reference_targets.tsv");
    let entrapments = stats_file("entrapment_reference_entrapments.tsv");
    let entrapment = [
        (
            "stats_entrapment_fdp_reference.json",
            entrapment_fdp(&targets, &entrapments, 1.0),
        ),
        (
            "stats_entrapment_fdp_reference_r2.json",
            entrapment_fdp_with(
                &targets,
                &entrapments,
                2.0,
                &EntrapmentFdpOptions {
                    thresholds: Some(vec![0.01, 0.05]),
                    ..Default::default()
                },
            ),
        ),
    ];
    for (fixture, live) in entrapment {
        let live = live.unwrap();
        let recorded = recording(fixture);
        assert_columns_within(&live.columns, &recorded, CROSS_PLATFORM);
        assert_eq!(
            serde_json::to_value(&live.caveats).unwrap(),
            recorded["caveats"]
        );
    }
}

#[test]
fn the_new_verbs_cross_their_refusals_as_usage_errors() {
    let Some(()) = require_verb("stats fit-mixed") else {
        return;
    };
    // A blank weight among weighted tests: what Rust sends for a test with no weight.
    usage_containing(
        combine_p(
            &[
                PValueTest::new("A", 0.1).weight(2.0),
                PValueTest::new("A", 0.2),
            ],
            "stouffer",
        ),
        "every line or on none",
    );
    usage_containing(
        combine_p(&[PValueTest::new("A", 0.1).weight(2.0)], "fisher"),
        "takes no weights",
    );
    usage_containing(
        combine_p(&[PValueTest::new("A", 1.5)], "fisher"),
        "outside [0, 1]",
    );

    let dir = scratch("refusals");
    let trait_table = dir.join("trait.tsv");
    std::fs::write(&trait_table, "sample\tage\tsex\n").unwrap();
    usage_containing(
        correlate(stats_file("limma_reference_responses.tsv"), &trait_table),
        "has 2 value columns",
    );
    let proportions = dir.join("proportions.tsv");
    std::fs::write(&proportions, "feature\tV1\tV2\nA\t0.5\t1\n").unwrap();
    let design = dir.join("design.tsv");
    std::fs::write(&design, "sample\tintercept\nV1\t1\nV2\t1\n").unwrap();
    usage_containing(
        fit_logistic(&proportions, &design, &["intercept"]),
        "feature 'A', sample 'V1'",
    );
    usage_containing(
        fit_logistic_with(
            stats_file("glm_reference_responses.tsv"),
            stats_file("glm_reference_design.tsv"),
            &["x"],
            &FitLogisticOptions {
                max_iterations: Some(0),
                ..Default::default()
            },
        ),
        "refuses --max-iterations 0",
    );
    let groups = dir.join("groups.tsv");
    std::fs::write(&groups, "sample\tday\tsite\n").unwrap();
    usage_containing(
        fit_mixed(
            stats_file("lmm_reference_responses.tsv"),
            stats_file("lmm_reference_design.tsv"),
            &groups,
            &["age"],
        ),
        "needs exactly two",
    );
    usage_containing(
        entrapment_fdp(
            stats_file("entrapment_reference_targets.tsv"),
            stats_file("entrapment_reference_entrapments.tsv"),
            0.0,
        ),
        "ratio must be a positive finite number",
    );
    let extra = dir.join("targets.tsv");
    std::fs::write(&extra, "peptide\tscore\tq_value\nPEPTIDE\t10\t0.01\n").unwrap();
    usage_containing(
        entrapment_fdp(
            &extra,
            stats_file("entrapment_reference_entrapments.tsv"),
            1.0,
        ),
        "Not a column this file takes: peptide",
    );
    let _ = std::fs::remove_dir_all(&dir);
}
