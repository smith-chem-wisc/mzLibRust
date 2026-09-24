//! Live check of median polish against the real bridge.
//!
//! `quant median-polish` opens only the peptide table it is given, so this needs a bridge and
//! nothing else: no network, no mzML. The offline suite pins the projection against the recorded
//! payload; this pins that the bridge still reads a FlashLFQ-shaped table and labels samples the
//! way the docs say. Skips when no bridge is staged. Run with `cargo test --features live`.

#![cfg(feature = "live")]

mod support;

use mzlib::flashlfq::{median_polish, median_polish_with, DesignEntry, MedianPolishOptions};
use support::require_bridge;

/// A QuantifiedPeptides.tsv in FlashLFQ's layout: two runs, three proteins, one of them with a
/// single peptide seen in one run only.
const TABLE: &str = "Sequence\tBase Sequence\tProtein Groups\tGene Names\tOrganism\t\
Intensity_run_3\tIntensity_run_4\tDetection Type_run_3\tDetection Type_run_4\n\
PEPTIDEK\tPEPTIDEK\tP1\tGENE1\tHomo sapiens\t1000\t2000\tMSMS\tMSMS\n\
ACDEFGHIK\tACDEFGHIK\tP1\tGENE1\tHomo sapiens\t3000\t6000\tMSMS\tMBR\n\
LMNPQR\tLMNPQR\tP2\tGENE2\tHomo sapiens\t500\t0\tMSMS\tNotDetected\n";

fn table() -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("mzlibrust-live-quant-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("QuantifiedPeptides.tsv");
    std::fs::write(&path, TABLE).unwrap();
    path
}

#[test]
fn with_no_design_every_run_is_its_own_sample() {
    let Some(()) = require_bridge() else { return };
    let result = median_polish(table()).expect("a FlashLFQ-shaped table should roll up");

    let labels: Vec<&str> = result.samples.iter().map(|s| s.label.as_str()).collect();
    assert_eq!(labels, ["run_3", "run_4"]);
    assert_eq!(result.peptide_count, 3);
    assert_eq!(result.protein_count as usize, result.proteins.len());
    let p1 = result
        .proteins
        .iter()
        .find(|p| p.protein_group == "P1")
        .unwrap();
    assert!(p1.intensity("run_3").is_some_and(|v| v > 0.0));
}

#[test]
fn a_design_relabels_the_samples_by_condition_and_replicate() {
    let Some(()) = require_bridge() else { return };
    let result = median_polish_with(
        table(),
        &MedianPolishOptions {
            design: vec![
                DesignEntry::new("run_3")
                    .condition("control")
                    .biological_replicate(0),
                DesignEntry::new("run_4")
                    .condition("treated")
                    .biological_replicate(0),
            ],
            ..MedianPolishOptions::default()
        },
    )
    .expect("the design names every run in the table");

    let labels: Vec<&str> = result.samples.iter().map(|s| s.label.as_str()).collect();
    assert_eq!(labels, ["control_1", "treated_1"]);
    for protein in &result.proteins {
        assert!(protein
            .intensities
            .keys()
            .all(|k| labels.contains(&k.as_str())));
    }
}

#[test]
fn a_design_that_misses_a_run_is_a_usage_error() {
    let Some(()) = require_bridge() else { return };
    let error = median_polish_with(
        table(),
        &MedianPolishOptions {
            design: vec![DesignEntry::new("run_3").condition("control")],
            ..MedianPolishOptions::default()
        },
    )
    .unwrap_err();
    assert!(matches!(error, mzlib::MzLibError::Usage(_)), "{error}");
}
