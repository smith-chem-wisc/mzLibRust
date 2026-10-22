//! The ptmqtl verbs through the real bridge, over mzLib's hand-computed PtmQtl cases.
//!
//! The offline suite replays recorded payloads, so it cannot see the bridge and this crate drift
//! apart. These run the real executable on the files in `tests/fixtures/ptmqtl/` and check that
//! what it returns today is what was recorded, so the docs' printed numbers are still true; that
//! the file `occupancy` writes is the one `trait_effects` reads; and that the filters and refusals
//! reach mzLib. Local files only, no network. pyMzLib's `tests/test_ptmqtl_live.py`, test for test.
//!
//! Run with `cargo test --features live`; skips on a bridge without the ptmqtl verbs.

#![cfg(feature = "live")]

mod support;

use std::path::PathBuf;

use mzlib::ptmqtl::{
    occupancy, occupancy_with, pairs, pairs_with, trait_effects, trait_effects_with,
    OccupancyOptions, PairType, PairsOptions, TraitEffectsOptions,
};
use serde_json::Value;
use support::require_verb;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
}

fn ptmqtl_file(name: &str) -> PathBuf {
    fixtures().join("ptmqtl").join(name)
}

fn recording(name: &str) -> Value {
    serde_json::from_str(&std::fs::read_to_string(fixtures().join(name)).unwrap()).unwrap()
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("mzlib-live-ptmqtl-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Equal, except that floats may differ in their last bits across platforms (1e-9 relative).
fn assert_same_columns(live: &mzlib::readers::Table, recorded: &Value) {
    let names: Vec<&str> = recorded["column_names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    assert_eq!(live.names(), names);
    let close = |a: &Value, b: &Value, at: &str| match (a.as_f64(), b.as_f64()) {
        (Some(a), Some(b)) => assert!((a - b).abs() <= 1e-9 * b.abs().max(1e-12), "{at}: {a} {b}"),
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
                    a.iter().zip(b).for_each(|(x, y)| close(x, y, &at));
                }
                _ => close(a, b, &at),
            }
        }
    }
}

fn check(
    fixture: &str,
    columns: &mzlib::readers::Table,
    caveats: &[String],
    counts: Option<(&str, Value)>,
) {
    let recorded = recording(fixture);
    assert_same_columns(columns, &recorded);
    assert_eq!(
        serde_json::to_value(caveats).unwrap(),
        recorded["caveats"],
        "{fixture}"
    );
    // Keys from mzLib enums; a repin that grows one must not slip past the recording.
    if let Some((key, live)) = counts {
        assert_eq!(live, recorded[key], "{fixture} {key}");
    }
}

#[test]
fn the_recordings_are_what_the_bridge_returns_today() {
    let Some(()) = require_verb("ptmqtl occupancy") else {
        return;
    };
    let counts = |v: &std::collections::BTreeMap<String, u64>| serde_json::to_value(v).unwrap();

    let occ = occupancy(ptmqtl_file("cases.tsv")).unwrap();
    check(
        "ptmqtl_occupancy_cases.json",
        &occ.columns,
        &occ.caveats,
        Some(("state_counts", counts(&occ.state_counts))),
    );
    let occ = occupancy_with(
        ptmqtl_file("fractions.tsv"),
        &OccupancyOptions {
            samples: Some(ptmqtl_file("fractions_samples.tsv")),
            ..Default::default()
        },
    )
    .unwrap();
    check(
        "ptmqtl_occupancy_fractions.json",
        &occ.columns,
        &occ.caveats,
        Some(("state_counts", counts(&occ.state_counts))),
    );

    let found = pairs(ptmqtl_file("cases.tsv"), PairType::Physical).unwrap();
    check(
        "ptmqtl_pairs_physical.json",
        &found.columns,
        &found.caveats,
        None,
    );
    let found = pairs(ptmqtl_file("batch.tsv"), PairType::Covarying).unwrap();
    check(
        "ptmqtl_pairs_batch.json",
        &found.columns,
        &found.caveats,
        None,
    );
    let found = pairs_with(
        ptmqtl_file("batch.tsv"),
        PairType::Covarying,
        &PairsOptions {
            remove_run_level: true,
            ..Default::default()
        },
    )
    .unwrap();
    check(
        "ptmqtl_pairs_batch_run_level.json",
        &found.columns,
        &found.caveats,
        None,
    );

    let fit = trait_effects(
        ptmqtl_file("trait_occupancy.tsv"),
        ptmqtl_file("trait_runs.tsv"),
    )
    .unwrap();
    check(
        "ptmqtl_trait_effects.json",
        &fit.columns,
        &fit.caveats,
        Some(("status_counts", counts(&fit.status_counts))),
    );
    let fit = trait_effects_with(
        ptmqtl_file("trait_occupancy.tsv"),
        ptmqtl_file("trait_runs.tsv"),
        &TraitEffectsOptions {
            min_replicates: Some(6),
            ..Default::default()
        },
    )
    .unwrap();
    check(
        "ptmqtl_trait_effects_strict.json",
        &fit.columns,
        &fit.caveats,
        Some(("status_counts", counts(&fit.status_counts))),
    );
}

#[test]
fn occupancy_out_writes_the_table_trait_effects_reads() {
    let Some(()) = require_verb("ptmqtl occupancy") else {
        return;
    };
    let dir = scratch("out");
    let out = dir.join("nested").join("occupancy.tsv");
    let occ = occupancy_with(
        ptmqtl_file("trait.tsv"),
        &OccupancyOptions {
            out: Some(out.clone()),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(occ.columns.is_empty());
    assert_eq!(
        occ.written.as_ref().map(|w| w.row_count),
        Some(occ.row_count)
    );
    // The committed file may be checked out with CRLF; the bridge writes LF.
    let lf = |text: String| text.replace("\r\n", "\n");
    assert_eq!(
        lf(std::fs::read_to_string(&out).unwrap()),
        lf(std::fs::read_to_string(ptmqtl_file("trait_occupancy.tsv")).unwrap())
    );
    let fit = trait_effects(&out, ptmqtl_file("trait_runs.tsv")).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    assert_same_columns(&fit.columns, &recording("ptmqtl_trait_effects.json"));
}

#[test]
fn filters_reach_mzlib() {
    let Some(()) = require_verb("ptmqtl occupancy") else {
        return;
    };
    let every = occupancy(ptmqtl_file("cases.tsv")).unwrap();
    let biological = occupancy_with(
        ptmqtl_file("cases.tsv"),
        &OccupancyOptions {
            categories: Some(vec!["Common Biological".into()]),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(biological
        .rows()
        .unwrap()
        .iter()
        .all(|r| r.category == "Common Biological"));
    assert_eq!(
        biological.categories,
        Some(vec!["Common Biological".to_owned()])
    );
    assert!(biological.row_count < every.row_count);
    let phospho = occupancy_with(
        ptmqtl_file("cases.tsv"),
        &OccupancyOptions {
            modifications: Some(vec!["Common Biological:Phosphorylation on S".into()]),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(phospho
        .rows()
        .unwrap()
        .iter()
        .all(|r| r.modification == "Common Biological:Phosphorylation on S"));
}

#[test]
fn refusals_cross_as_usage_errors() {
    let Some(()) = require_verb("ptmqtl pairs") else {
        return;
    };
    let usage = |error: mzlib::MzLibError, needles: &[&str]| {
        assert!(
            matches!(&error, mzlib::MzLibError::Usage(m) if needles.iter().any(|n| m.contains(n))),
            "{needles:?}: {error}"
        );
    };
    usage(
        pairs_with(
            ptmqtl_file("batch.tsv"),
            PairType::Physical,
            &PairsOptions {
                remove_run_level: true,
                ..Default::default()
            },
        )
        .unwrap_err(),
        &["shapes co-varying pairs only"],
    );
    usage(
        pairs_with(
            ptmqtl_file("batch.tsv"),
            PairType::Covarying,
            &PairsOptions {
                run_level_group: Some("category".into()),
                ..Default::default()
            },
        )
        .unwrap_err(),
        &["give --remove-run-level as well"],
    );
    let dir = scratch("refusals");
    usage(
        occupancy_with(
            ptmqtl_file("cases.tsv"),
            &OccupancyOptions {
                out: Some(dir.join("occupancy.csv")),
                ..Default::default()
            },
        )
        .unwrap_err(),
        &["must end in", "must name a .tsv"],
    );
    let repeated = dir.join("psms.tsv");
    let text = std::fs::read_to_string(ptmqtl_file("cases.tsv")).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    std::fs::write(
        &repeated,
        format!("{}\n{}\n{}\n", lines[0], lines[1], lines[1]),
    )
    .unwrap();
    usage(occupancy(&repeated).unwrap_err(), &["not one per PSM"]);
    usage(
        trait_effects_with(
            ptmqtl_file("trait_occupancy.tsv"),
            ptmqtl_file("trait_runs.tsv"),
            &TraitEffectsOptions {
                logit_epsilon: Some(0.5),
                ..Default::default()
            },
        )
        .unwrap_err(),
        &["above 0 and below 0.5"],
    );
    let _ = std::fs::remove_dir_all(&dir);
}
