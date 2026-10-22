//! Differential abundance, set up by mzLib: the experiment's design, and every peptide in every
//! sample.
//!
//! | You want to know | Call | mzLib computes it with |
//! |---|---|---|
//! | What is my design matrix, and which comparisons can these samples answer? | [`design_with`] | `AnalysisDesign.ByStratum`, `.DefaultContrasts`, `.Estimability` |
//! | What did each peptide measure in each sample, and if nothing, why? | [`observations_with`] | `FlashLfqObservations.FromPeptideTable`, `ObservationTable` |
//! | How far apart do the samples sit on the peptides they share? | [`observations_with`] | `SharedPeptideNormalization.Measure` / `.Apply` |
//!
//! Before any differential test there are two tables to get right, and both are easy to get
//! quietly wrong. The **design**: which coefficient means what, which level is the baseline,
//! whether two things you want to separate can be separated by these samples at all. And the
//! **observations**: what each peptide measured in each biological sample, with fractions summed
//! and technical replicates averaged correctly, and a missing value kept as missing rather than
//! written as 0. This module hands both to mzLib (mzLib #1450, #1455). [`design_with`] turns the
//! sample sheet you already keep into mzLib's design matrix, which equals R's `model.matrix`, and
//! says which comparisons cannot be estimated and why. [`observations_with`] turns a FlashLFQ or
//! MetaMorpheus peptide table into every peptide in every sample, as a log2 intensity or `None`
//! with the reason. Both write files that [`crate::stats::fit_with`] reads, so the analysis needs no
//! glue code of yours. This crate does no arithmetic of its own.
//!
//! ## From a sample sheet to a design
//!
//! The [`stats`](crate::stats) guide analyses a MALAT1 dilution series — 27 runs at 500, 250 or
//! 125 ng — with a hand-made design of 0/1 columns. Here is the same experiment as a scientist
//! writes it down, one row per run and its load (`run`, `malat_load`). Say which column is a factor
//! and which level is the baseline; everything else is mzLib's:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::diff::{design_with, DesignOptions, FactorSpec};
//!
//! let design = design_with(
//!     "malat_dilution_samples.tsv",
//!     &DesignOptions {
//!         factors: vec![FactorSpec::new("malat_load", Some("500ng"))],
//!         out: Some("malat_design".into()),
//!         ..Default::default()
//!     },
//! )?;
//! let stratum = design.stratum("all")?;
//! assert_eq!(stratum.column_names, ["intercept", "malat_load=125ng", "malat_load=250ng"]);
//! assert_eq!((stratum.rank, stratum.aliased_terms.len(), stratum.warnings.len()), (3, 0, 0));
//!
//! let contrasts = design.contrasts()?;
//! assert_eq!(contrasts[0].label, "malat_load=125ng vs malat_load=500ng");
//! assert_eq!(contrasts[0].estimability, "estimable");
//! assert_eq!(contrasts[0].fitted_weights, Some(vec![0.0, 1.0, 0.0]));
//! assert_eq!(contrasts[1].contrast_id, "c2");
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! What came back:
//!
//! - **The columns** are named `factor=level`. The reference comes first and has no column; the
//!   other levels follow in ordinal (character) order, so `125ng` before `250ng`. That is R's
//!   `model.matrix` with the same level order.
//! - **The default contrasts** are each level against the reference. Each is one design column
//!   here, so `fitted_weights` has a single 1.
//! - **[`DesignOptions::out`]** wrote `design_1.tsv`, the design table
//!   [`crate::stats::fit_with`] reads, and `contrasts_1.tsv`, each estimable contrast's weights on
//!   the same columns, which [`crate::stats::FitOptions::contrasts`] reads.
//!
//! Fit it. Each contrast against the reference is its design column, tested by name, and the
//! contrasts file adds a confidence interval to each answer. Name the coefficients and the file, and
//! one call returns both tables:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::diff::{design_with, DesignOptions, FactorSpec};
//! use mzlib::stats::{fit_with, FitOptions};
//!
//! # let design = design_with(
//! #     "malat_dilution_samples.tsv",
//! #     &DesignOptions {
//! #         factors: vec![FactorSpec::new("malat_load", Some("500ng"))],
//! #         out: Some("malat_design".into()),
//! #         ..Default::default()
//! #     },
//! # )?;
//! let stratum = design.stratum("all")?;
//! let design_file = stratum.design_file.as_deref().unwrap_or_default();
//! let coefficients = ["malat_load=250ng", "malat_load=125ng"];
//! let fit = fit_with("malat_dilution_log2_vs_fluc.tsv", design_file, &coefficients, &FitOptions::default())?;
//! assert_eq!(fit.status_counts["fitted"], 212);
//!
//! let median_change = |coefficient: &str| -> mzlib::Result<f64> {
//!     let mut changes: Vec<f64> = fit
//!         .rows(coefficient)?
//!         .into_iter()
//!         .filter(|r| r.status == "fitted" && r.feature.starts_with("MALAT1:"))
//!         .filter_map(|r| r.estimate)
//!         .collect();
//!     changes.sort_by(f64::total_cmp);
//!     Ok(((changes[57] + changes[58]) / 2.0 * 100.0).round() / 100.0)
//! };
//! assert_eq!((median_change("malat_load=250ng")?, median_change("malat_load=125ng")?), (-0.89, -1.01));
//!
//! // The contrasts file too: c2 is "250 ng vs 500 ng", the malat_load=250ng column.
//! let both = fit_with(
//!     "malat_dilution_log2_vs_fluc.tsv",
//!     design_file,
//!     &coefficients,
//!     &FitOptions { contrasts: stratum.contrasts_file.as_deref().map(Into::into), ..Default::default() },
//! )?;
//! let contrasts = both.contrasts.as_ref().expect("contrasts were asked for");
//! let same = both
//!     .rows("malat_load=250ng")?
//!     .iter()
//!     .zip(contrasts.rows("c2")?)
//!     .all(|(a, b)| a.t == b.t && a.bh_adjusted == b.bh_adjusted);
//! assert!(same);
//! let row = contrasts.rows("c2")?.into_iter()
//!     .find(|r| r.status == "fitted" && r.feature.starts_with("MALAT1:")).unwrap();
//! let two = |x: Option<f64>| (x.unwrap() * 100.0).round() / 100.0;
//! assert_eq!((row.feature.as_str(), row.observed), ("MALAT1:AAAAAAAAACUCG", Some(4)));
//! assert_eq!((two(row.estimate), two(row.confidence_low), two(row.confidence_high)), (-0.75, -3.68, 2.17));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! The same numbers as the hand-made design, because it is the same model: the bridge's tests hold
//! the two fits equal, feature by feature. A contrast of one column is that column's test, number
//! for number, so `c2` repeats the `malat_load=250ng` rows and adds the interval. This oligo was
//! measured in 4 runs, and its 95% interval runs from -3.7 to +2.2 log2 units: it cannot tell
//! halving from no change. The contrasts file earns its keep for comparisons that are no single
//! column, such as 250 ng against 125 ng: declare no reference level and [`design_with`] writes
//! every pair, as the [`stats`](crate::stats#contrasts-any-comparison-with-a-confidence-interval)
//! guide shows.
//!
//! ## Covariates, units and batch
//!
//! A covariate is a number in a unit. Say the unit, and how the slope should read: here age in
//! years, as a slope per decade, centred at 50 years. Batch is a fixed block. The sheet is mzLib's
//! own test design, and R's `model.matrix` output for it ships beside it:
//!
//! ```
//! # mzlib_replay::activate();
//! # let r_matrix = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/diff/design_factor_batch_covariate_matrix.tsv");
//! use mzlib::diff::{design_with, CovariateSpec, DesignOptions, FactorSpec};
//!
//! let design = design_with(
//!     "design_factor_batch_covariate_samples.tsv",
//!     &DesignOptions {
//!         factors: vec![FactorSpec::new("condition", Some("ctrl"))],
//!         covariates: vec![CovariateSpec::new("age", "years").scale(10.0).scaled_unit("decade").centre(50.0)],
//!         batch: Some("batch".into()),
//!         ..Default::default()
//!     },
//! )?;
//! let stratum = &design.strata[0];
//! assert_eq!(stratum.column_names, ["intercept", "condition=treated", "batch=b2", "age (per decade)"]);
//! let slope = &design.contrasts()?[1];
//! assert_eq!(slope.label, "age (per decade)");
//! assert_eq!(slope.covariate_unit.as_deref(), Some("decade"));
//! assert_eq!(slope.covariate_scale.as_deref(), Some("years/10, centred at 50 years"));
//!
//! // R's model.matrix(~ condition + batch + I((age - 50)/10)), value for value.
//! let r = std::fs::read_to_string(r_matrix).unwrap();
//! let r_rows: Vec<Vec<f64>> = r
//!     .lines()
//!     .skip(1)
//!     .map(|line| line.split('\t').map(|cell| cell.parse().unwrap()).collect())
//!     .collect();
//! assert_eq!(r_rows.len(), stratum.matrix.len());
//! let same = stratum
//!     .matrix
//!     .iter()
//!     .zip(&r_rows)
//!     .all(|(ours, theirs)| ours.iter().zip(theirs).all(|(a, b)| (a - b).abs() < 1e-12));
//! assert!(same);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Leave a field out and mzLib chooses.** `CovariateSpec::new("age", "years")` sends no scale,
//! scaled unit or centre, so the slope is per year and the centre is the stratum's mean;
//! [`DesignStratum::covariate_centres`] reports the centre used. Batch enters the model only when
//! you name its column, and drops out with a warning in a stratum that has only one batch.
//!
//! ## When the samples cannot answer
//!
//! Here every sample of condition A ran in batch 1 and every sample of B in batch 2:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::diff::{design_with, DesignOptions, FactorSpec};
//!
//! let design = design_with(
//!     "design_confounded_samples.tsv",
//!     &DesignOptions {
//!         factors: vec![FactorSpec::new("condition", Some("A"))],
//!         batch: Some("batch".into()),
//!         ..Default::default()
//!     },
//! )?;
//! let stratum = &design.strata[0];
//! assert_eq!(stratum.column_names, ["intercept", "condition=B", "batch=2"]);
//! assert_eq!((stratum.rank, stratum.aliased_terms.clone()), (2, vec!["condition".to_owned(), "batch".to_owned()]));
//! let contrast = &design.contrasts()?[0];
//! assert_eq!(contrast.label, "condition=B vs condition=A");
//! assert_eq!((contrast.estimability.as_str(), contrast.fitted_weights.clone()), ("not_estimable:confounded", None));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **No analysis of these samples can separate the condition from the batch**, and mzLib says so
//! instead of returning a number. The verdicts ([`ESTIMABILITIES`]) are `estimable`;
//! `not_estimable:confounded`, a comparison along a direction two or more terms share, as here; and
//! `not_estimable:rank_deficient`, one that needs a column no sample informs, such as an interaction
//! cell with no samples. The design file still holds a design that can be fitted:
//! [`DesignStratum::fitted_column_names`] leaves out the later of two aliased columns, because a fit
//! refuses a design without full rank. A contrast that is not estimable has no fitted weights, so it
//! never reaches the contrasts file.
//!
//! ## Strata: one design per tissue
//!
//! When samples are different kinds of thing — liver and muscle — a comparison is made within each
//! kind. Name the stratum column, and each stratum gets its own design:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::diff::{design_with, DesignOptions, FactorSpec, NotRunContrast};
//!
//! let design = design_with(
//!     "stratified_samples.tsv",
//!     &DesignOptions {
//!         factors: vec![FactorSpec::new("condition", Some("A"))],
//!         strata: vec!["organism part".into()],
//!         ..Default::default()
//!     },
//! )?;
//! let strata: Vec<(&str, usize)> = design.strata.iter().map(|s| (s.stratum.as_str(), s.sample_ids.len())).collect();
//! assert_eq!(strata, [("organism part=kidney", 2), ("organism part=liver", 3), ("organism part=muscle", 3)]);
//! assert_eq!(
//!     design.stratum("organism part=kidney")?.not_run,
//!     [NotRunContrast { label: "condition (one level: A)".into(), reason: "not_in_stratum".into() }]
//! );
//! let ids: Vec<(String, String)> = design.contrasts()?.into_iter().map(|c| (c.stratum, c.contrast_id)).collect();
//! assert_eq!(ids[0], ("organism part=liver".to_owned(), "c1".to_owned()));
//! assert_eq!(ids[1], ("organism part=muscle".to_owned(), "c1".to_owned()));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Kidney has only condition A, so there is nothing to compare there, and `not_run` says so. **A
//! contrast id is numbered within its stratum:** `c1` is liver's first comparison and muscle's
//! first, so match strata by `label`. With [`DesignOptions::out`], stratum *i* writes
//! `design_i.tsv`, numbered in the order of [`Design::strata`], because a stratum key can hold
//! characters a file name cannot.
//!
//! ## From a peptide table to observations
//!
//! [`observations_with`] reads the peptide table FlashLFQ (`QuantifiedPeptides.tsv`) or MetaMorpheus
//! (`AllQuantifiedPeptides.tsv`) wrote, and needs the runs: which biological sample each file
//! belongs to and, if it applies, its fraction and technical replicate. **Grouping comes only from
//! what you give here**; mzLib never reads it out of a file name. This MetaMorpheus table is mzLib's
//! test data: three peptides in 18 runs, each run its own sample. MetaMorpheus quantifies calibrated
//! files, so its columns end in `-calib`; give the file name without it, and mzLib matches either:
//!
//! ```
//! # mzlib_replay::activate();
//! # let table = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/diff/MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv");
//! use mzlib::diff::{observations_with, ObservationRun, ObservationsOptions, QuantBasis};
//!
//! let header = std::fs::read_to_string(table).unwrap();
//! let runs: Vec<ObservationRun> = header
//!     .lines()
//!     .next()
//!     .unwrap()
//!     .split('\t')
//!     .filter_map(|column| column.strip_prefix("Intensity_"))
//!     .map(|file| file.trim_end_matches("-calib"))
//!     .map(|file| ObservationRun::new(file, file))
//!     .collect();
//! let options = ObservationsOptions {
//!     protein_groups: Some("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv".into()),
//!     ..Default::default()
//! };
//! let peptides = "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv";
//! let table = observations_with(peptides, &runs, QuantBasis::MbrKept, &options)?;
//! assert_eq!((table.peptide_count, table.sample_count, table.row_count), (3, 18, 54));
//! assert_eq!((table.state_counts["quantified"], table.state_counts["mbr_transferred"]), (8, 1));
//!
//! let rows = table.rows()?;
//! let transferred = rows.iter().find(|r| r.state == "mbr_transferred").unwrap();
//! assert_eq!((transferred.base_sequence.as_str(), transferred.sample.as_str()), ("CACASHVAK", "QE-002118_GM6_a"));
//! assert_eq!(transferred.log2_intensity.map(|v| (v * 100.0).round() / 100.0), Some(17.14));
//! assert_eq!(transferred.unique_protein_group.as_deref(), Some("Q16643"));
//!
//! // Under msms_only the same cell has no value, and still says why.
//! let only = observations_with(peptides, &runs, QuantBasis::MsmsOnly, &options)?;
//! let same = only.rows()?.into_iter()
//!     .find(|r| r.sample == "QE-002118_GM6_a" && r.base_sequence == "CACASHVAK").unwrap();
//! assert_eq!((same.log2_intensity, same.state.as_str()), (None, "mbr_transferred"));
//!
//! // The protein-group table is how a contaminant or a decoy is recognised.
//! let flagged: Vec<(&str, Option<&str>)> = table
//!     .protein_groups
//!     .iter()
//!     .filter(|g| g.is_contaminant || g.is_decoy)
//!     .map(|g| (g.name.as_str(), g.genes.as_deref()))
//!     .collect();
//! assert_eq!(flagged, [("P02769", Some("ALB")), ("DECOY_P62750", Some("RPL23A"))]);
//! assert!(table.warnings[0].contains("not in the protein-group table"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **[`QuantBasis`] decides whether match-between-runs values count**, and mzLib has no default: it
//! is your scientific choice. What each `state` means ([`OBSERVATION_STATES`]):
//!
//! - `quantified`: a value from MS/MS-identified peaks only.
//! - `mbr_transferred`: the value includes a match-between-runs peak. Under
//!   [`QuantBasis::MsmsOnly`] it is left out, and `log2_intensity` is `None`.
//! - `ambiguous_peak`: a peak was found but withheld as ambiguous. No value.
//! - `identified_not_quantified`: identified by MS/MS, but no peak could be quantified. No value.
//! - `not_detected`: no value.
//!
//! A sample whose runs disagree takes the most informative state among them. With fractions, mzLib
//! sums a technical replicate's fractions on the linear scale, then averages the sample's
//! replicates in log2; give each run's [`ObservationRun::fraction`] and
//! [`ObservationRun::technical_replicate`].
//!
//! ## Normalization, and the assumption it makes
//!
//! [`observations_with`] always measures how far each sample sits from the others, and subtracts it
//! only when you ask ([`ObservationsOptions::normalize`]). mzLib's rule: over the peptides with a
//! value in every sample (decoys, contaminants and entrapments excluded), each sample's shift is
//! the median of its log2 value minus the peptide's mean. Two runs of FlashLFQ's E. coli spike-in
//! vignette, which their names say were spiked at 1x and 2x:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::diff::{observations, observations_with, ObservationRun, ObservationsOptions, QuantBasis};
//!
//! let runs = [
//!     ObservationRun::new("09-04-18_EcoliSpikeInSingleShot1x", "spike_1x"),
//!     ObservationRun::new("09-04-18_EcoliSpikeInSingleShot2x", "spike_2x"),
//! ];
//! let peptides = "EcoliSpikeIn_QuantifiedPeptides.tsv";
//! let spike = observations(peptides, &runs, QuantBasis::MsmsOnly)?;
//! assert_eq!((spike.normalization.setting.as_str(), spike.normalization.reference_set_size), ("none", 86));
//! assert_eq!(spike.normalization.minimum_reference_set_size, 100); // mzLib's floor
//!
//! // Lower the floor, subtract the shifts, and write the table stats::fit_with reads.
//! let shifted = observations_with(
//!     peptides,
//!     &runs,
//!     QuantBasis::MsmsOnly,
//!     &ObservationsOptions {
//!         min_peptides: Some(50),
//!         normalize: true,
//!         responses_out: Some("spike_in_responses.tsv".into()),
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!(shifted.normalization.setting, "shared_peptide_median");
//! let shift = |sample: &str| shifted.normalization.shifts[sample].map(|s| (s * 100.0).round() / 100.0);
//! assert_eq!((shift("spike_1x"), shift("spike_2x")), (Some(-0.47), Some(0.47)));
//! let written = shifted.responses_written.as_ref().unwrap();
//! assert_eq!((written.row_count, written.sample_count), (243, 2));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Below mzLib's floor of 100 shared peptides, nothing is shifted**, because a median of a few
//! peptides can move a whole sample by itself. And look at what lowering it did here: **the
//! normalization removed the spike-in.** This slice of the vignette holds almost only E. coli
//! peptides, so nearly every peptide was doubled; the shifts of -0.47 and +0.47 take all of it
//! away. Median normalization assumes most peptides do not change between samples. When most of
//! them do — a table of the spiked organism alone, a loading series, a pull-down — that assumption
//! is false, and normalizing deletes the result. Read the shifts before you subtract them.
//!
//! [`ObservationsOptions::responses_out`] wrote the table's own peptide-by-sample shape, which
//! [`crate::stats::fit_with`] reads as its responses: rows are peptides by full sequence, columns
//! the samples, and an empty cell means no value. For more than one stratum, narrow the call to one
//! stratum's samples with [`ObservationsOptions::samples`] (a design file from [`design_with`] will
//! do): the shifts are measured over every sample of the table, so measure each stratum on its own.
//!
//! ## What these functions will not do
//!
//! - **No protein roll-up and no row choice.** mzLib does not yet turn the observation table into a
//!   protein-by-sample table, choose which rows a test uses (decoys, contaminants and q-value cuts
//!   included), or assemble fits into its differential result table. This crate does not do it
//!   either. Every peptide is a row; filter with [`Observation::protein_groups`],
//!   [`Observation::unique_protein_group`] and [`Observations::protein_groups`].
//! - **Only the default contrasts.** [`design_with`] lists each level against its reference (every
//!   pair when there is none) and each covariate's slope. mzLib's custom contrasts have no wire
//!   form yet; write your own contrasts file for [`crate::stats::FitOptions::contrasts`].
//! - **The peptide table does not record FlashLFQ's design**, so mzLib cannot check that your runs
//!   group the files as FlashLFQ did. Give the runs the design FlashLFQ quantified with.
//! - **A table written in a comma-decimal locale does not read correctly.** FlashLFQ writes in the
//!   machine's culture, and mzLib reads with a decimal point.
//!
//! ## Cite
//!
//! The methods behind these numbers, from the verbs' specs:
//!
#![doc = include_str!("../docs/reference/cite.diff.md")]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};
use crate::readers::{Table, WrittenTable};

/// Every value of [`design_with`]'s `estimability` column, as mzLib's `AnalysisDesign` writes them.
/// Only `"estimable"` rows have `fitted_weights`.
pub const ESTIMABILITIES: [&str; 3] = [
    "estimable",
    "not_estimable:confounded",
    "not_estimable:rank_deficient",
];

/// Every value of [`observations_with`]' `state` column, mzLib's `ObservationState` from most to
/// least informative. Only `"quantified"` and `"mbr_transferred"` rows can carry a value.
pub const OBSERVATION_STATES: [&str; 5] = [
    "quantified",
    "mbr_transferred",
    "ambiguous_peak",
    "identified_not_quantified",
    "not_detected",
];

/// The wire names of [`QuantBasis`], mzLib's `QuantBasis` machine names, in its order.
pub const QUANT_BASES: [&str; 2] = ["msms_only", "mbr_kept"];

/// The pyMzLib release whose bridge first has the diff verbs (mzLib 1.0.595): their specs'
/// `since.pymzlib`.
const MZLIB_1_0_595_BRIDGE: &str = "0.6.0";

const DESIGN_TIMEOUT: Duration = Duration::from_secs(60);
const OBSERVATIONS_TIMEOUT: Duration = Duration::from_secs(600);

// ---- inputs ------------------------------------------------------------------------------------

/// Whether match-between-runs values enter [`observations_with`]: mzLib's `QuantBasis`. There is no
/// default; you choose.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QuantBasis {
    /// `"msms_only"`: no match-between-runs value enters; such a cell has no value and says why.
    MsmsOnly,
    /// `"mbr_kept"`: match-between-runs values are kept.
    MbrKept,
}

impl QuantBasis {
    /// The wire name, one of [`QUANT_BASES`].
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::MsmsOnly => "msms_only",
            Self::MbrKept => "mbr_kept",
        }
    }
}

/// A condition factor of the model: mzLib's `FactorSpec`.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct FactorSpec {
    /// The sample sheet column holding each sample's level.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub name: String,
    /// The level every other level is compared against, or `None` to compare every pair of levels
    /// (mzLib then warns, because the comparisons are many).
    #[serde(default)]
    pub reference_level: Option<String>,
}

impl FactorSpec {
    /// The factor in column `name`, compared against `reference_level` (or every pair, for `None`).
    #[must_use]
    pub fn new(name: impl Into<String>, reference_level: Option<&str>) -> Self {
        Self {
            name: name.into(),
            reference_level: reference_level.map(str::to_owned),
        }
    }
}

/// A numeric covariate of the model, entered as `(value - centre) / scale`: mzLib's `CovariateSpec`.
///
/// Leave a field `None` and nothing is sent for it, so mzLib's own default applies.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct CovariateSpec {
    /// The sample sheet column holding each sample's value, in [`Self::unit`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub name: String,
    /// The stored unit, e.g. `"years"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub unit: String,
    /// Stored units per scaled unit, e.g. `10` years per decade. `None`: mzLib's 1.
    #[serde(default)]
    pub scale: Option<f64>,
    /// The unit a slope is per, e.g. `"decade"`. `None`: [`Self::scale`] stored units (`"years"`
    /// at scale 1, `"10 years"` at scale 10).
    #[serde(default)]
    pub scaled_unit: Option<String>,
    /// The value subtracted first, in [`Self::unit`]. `None`: the stratum's mean.
    #[serde(default)]
    pub centre: Option<f64>,
}

impl CovariateSpec {
    /// The covariate in column `name`, stored in `unit`, with mzLib's scale and centre.
    #[must_use]
    pub fn new(name: impl Into<String>, unit: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            unit: unit.into(),
            scale: None,
            scaled_unit: None,
            centre: None,
        }
    }

    /// This covariate with stored units per scaled unit, e.g. `10.0` years per decade.
    #[must_use]
    pub fn scale(mut self, scale: f64) -> Self {
        self.scale = Some(scale);
        self
    }

    /// This covariate with the unit a slope is per, e.g. `"decade"`.
    #[must_use]
    pub fn scaled_unit(mut self, unit: impl Into<String>) -> Self {
        self.scaled_unit = Some(unit.into());
        self
    }

    /// This covariate centred at `centre`, in its stored unit.
    #[must_use]
    pub fn centre(mut self, centre: f64) -> Self {
        self.centre = Some(centre);
        self
    }
}

/// One spectra file and its place in the design: mzLib's `ObservationRun`.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct ObservationRun {
    /// The file as the peptide table's `Intensity_` columns name it, without its extension. A
    /// column with `-calib` after it (MetaMorpheus's calibrated files) matches.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_name: String,
    /// The biological sample the run belongs to.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_id: String,
    /// The fraction; fractions of one technical replicate are summed. `None` (as an argument): not
    /// sent, so mzLib's 0. In a result, the value mzLib holds.
    #[serde(default)]
    pub fraction: Option<i64>,
    /// The technical replicate; a sample's replicates are averaged in log2. `None` (as an
    /// argument): not sent, so mzLib's 0. In a result, the value mzLib holds.
    #[serde(default)]
    pub technical_replicate: Option<i64>,
}

impl ObservationRun {
    /// The run `file_name` (without extension) of sample `sample_id`, with mzLib's fraction and
    /// technical replicate.
    #[must_use]
    pub fn new(file_name: impl Into<String>, sample_id: impl Into<String>) -> Self {
        Self {
            file_name: file_name.into(),
            sample_id: sample_id.into(),
            fraction: None,
            technical_replicate: None,
        }
    }

    /// This run as fraction `fraction` of its technical replicate.
    #[must_use]
    pub fn fraction(mut self, fraction: i64) -> Self {
        self.fraction = Some(fraction);
        self
    }

    /// This run as technical replicate `replicate` of its sample.
    #[must_use]
    pub fn technical_replicate(mut self, replicate: i64) -> Self {
        self.technical_replicate = Some(replicate);
        self
    }
}

/// Options for [`design_with`]: the model, and where to write. Every field is optional, and one
/// that is left empty sends nothing.
#[derive(Debug, Clone)]
pub struct DesignOptions {
    /// Condition factors, in column order, each a [`FactorSpec`]. With the covariates,
    /// interactions, strata and batch below, the model: one stdin line per term.
    pub factors: Vec<FactorSpec>,
    /// Numeric covariates, in column order, each a [`CovariateSpec`] with its stored unit.
    pub covariates: Vec<CovariateSpec>,
    /// Two-factor interactions, as `(first, second)` pairs of factor names.
    pub interactions: Vec<(String, String)>,
    /// Columns whose values split the samples into strata, each designed on its own.
    pub strata: Vec<String>,
    /// The column holding each sample's batch, fitted as a fixed block. `None`: no batch.
    pub batch: Option<String>,
    /// A directory to write `design_<i>.tsv` (for [`crate::stats::fit_with`]'s `design`) and
    /// `contrasts_<i>.tsv` (for [`crate::stats::FitOptions::contrasts`]) into, per stratum.
    /// `None`: nothing is written.
    pub out: Option<PathBuf>,
    /// How long to allow; `None` waits indefinitely. One minute by default.
    pub timeout: Option<Duration>,
}

impl Default for DesignOptions {
    fn default() -> Self {
        Self {
            factors: Vec::new(),
            covariates: Vec::new(),
            interactions: Vec::new(),
            strata: Vec::new(),
            batch: None,
            out: None,
            timeout: Some(DESIGN_TIMEOUT),
        }
    }
}

/// Options for [`observations_with`]. The default sends nothing beyond the table, the runs and
/// the basis, so every rule is mzLib's.
#[derive(Debug, Clone)]
pub struct ObservationsOptions {
    /// A MetaMorpheus `AllQuantifiedProteinGroups.tsv`, for each group's q-value and decoy,
    /// contaminant and entrapment flags. `None`: none described.
    pub protein_groups: Option<PathBuf>,
    /// A TSV whose first column (after a header row) lists the samples to keep, e.g. one stratum's
    /// design file from [`design_with`]. `None`: every sample.
    pub samples: Option<PathBuf>,
    /// Subtract each sample's shift. The shifts are measured and reported either way.
    pub normalize: bool,
    /// The floor on the normalization's reference set, in peptides. `None`: not sent, so mzLib's
    /// own floor (`DefaultMinimumPeptides`, 100 at the pin).
    pub min_peptides: Option<u32>,
    /// A `.tsv` path to write the long table to instead of returning its rows.
    pub out: Option<PathBuf>,
    /// A `.tsv` path to also write the peptide-by-sample log2 table that
    /// [`crate::stats::fit_with`] reads as `responses`.
    pub responses_out: Option<PathBuf>,
    /// How long to allow; `None` waits indefinitely. Ten minutes by default.
    pub timeout: Option<Duration>,
}

impl Default for ObservationsOptions {
    fn default() -> Self {
        Self {
            protein_groups: None,
            samples: None,
            normalize: false,
            min_peptides: None,
            out: None,
            responses_out: None,
            timeout: Some(OBSERVATIONS_TIMEOUT),
        }
    }
}

// ---- design results ----------------------------------------------------------------------------

/// A two-factor interaction of the model, as [`Design::interactions`] echoes it.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct Interaction {
    /// The first factor's name.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub first: String,
    /// The second factor's name.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub second: String,
}

/// A contrast the model asks for that one stratum cannot run: mzLib's `NotRunContrast`.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct NotRunContrast {
    /// Which contrast, e.g. `"condition (one level: A)"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub label: String,
    /// Why: `"not_in_stratum"` (one side has no sample there).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub reason: String,
}

/// Where [`design_with`] wrote its files ([`DesignOptions::out`]).
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct WrittenDesign {
    /// The absolute directory written into.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub directory: String,
    /// Design and contrast files written, in files.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
}

/// One stratum's design, as mzLib's `AnalysisDesign` built it.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct DesignStratum {
    /// The stratum key: `"all"`, or the stratum columns sorted by name as `name=value` joined by
    /// `;` (`%`, `;` and `=` in a value escaped).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub stratum: String,
    /// Samples in this stratum, in samples.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// The matrix rows, ordinal by sample id.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_ids: Vec<String>,
    /// mzLib's columns: `intercept`; `factor=level` for each non-reference level; `batch=level`;
    /// `name (per unit)` for a covariate; `a=level:b=level` for an interaction.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub column_names: Vec<String>,
    /// The columns a fit takes: all of them at full rank, else each kept unless it is redundant
    /// with those before it. The design file's coefficient columns.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub fitted_column_names: Vec<String>,
    /// The numerical rank, in columns, by mzLib's `LinearModel` full-rank test.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub rank: u64,
    /// Terms that share a direction with another and cannot be told apart from it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub aliased_terms: Vec<String>,
    /// What mzLib left out or changed here, in plain words.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub warnings: Vec<String>,
    /// The contrasts this stratum cannot run.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub not_run: Vec<NotRunContrast>,
    /// The centre subtracted from each covariate here, in the covariate's stored unit: the one you
    /// gave, else the stratum's mean.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub covariate_centres: BTreeMap<String, f64>,
    /// The design matrix, one row per [`Self::sample_ids`] entry and one value per
    /// [`Self::column_names`] entry, in its coding (0/1 indicators; covariate columns in scaled
    /// units). Equal to R's `model.matrix` for the same formula and level order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub matrix: Vec<Vec<f64>>,
    /// [`Self::matrix`] with only the fitted columns, in the same coding (0/1 indicators; covariate
    /// columns in scaled units): what [`crate::stats::fit_with`] fits.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub fitted_matrix: Vec<Vec<f64>>,
    /// The design table written for this stratum (`design_<i>.tsv`), or `None` without
    /// [`DesignOptions::out`]. Give it to [`crate::stats::fit_with`] as `design`.
    #[serde(default)]
    pub design_file: Option<String>,
    /// The contrast weights written for this stratum (`contrasts_<i>.tsv`), one row per estimable
    /// contrast on the design file's columns; `None` without [`DesignOptions::out`], or when nothing
    /// here is estimable. Give it to [`crate::stats::FitOptions::contrasts`], with
    /// [`Self::design_file`] as `design`.
    #[serde(default)]
    pub contrasts_file: Option<String>,
}

/// What [`design_with`] returns: each stratum's design, and every default contrast in one table.
///
/// [`Self::columns`] holds one row per (stratum, default contrast), strata in order;
/// [`Self::contrasts`] gives them as typed [`DefaultContrast`]s.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct Design {
    /// Absolute path of the sample sheet.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub samples_file: String,
    /// Samples read, across every stratum, in samples.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// The factors given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub factors: Vec<FactorSpec>,
    /// The covariates given; a `None` field was not given, so mzLib's default applied.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub covariates: Vec<CovariateSpec>,
    /// The interactions given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub interactions: Vec<Interaction>,
    /// The stratum columns given; empty for one stratum.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub stratum_columns: Vec<String>,
    /// The batch column, or `None` when batch is not in the model.
    #[serde(default)]
    pub batch_column: Option<String>,
    /// Sheet columns the model does not name, ignored.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub unused_columns: Vec<String>,
    /// Designs built, in strata.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub stratum_count: u64,
    /// One [`DesignStratum`] per stratum, in mzLib's order (ordinal by key).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub strata: Vec<DesignStratum>,
    /// Default contrasts per verdict across strata, in contrasts, keyed by [`ESTIMABILITIES`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub estimability_counts: BTreeMap<String, u64>,
    /// Where [`DesignOptions::out`] wrote, or `None` without it: nothing was written.
    #[serde(default)]
    pub written: Option<WrittenDesign>,
    /// Default contrasts across strata, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per (stratum, default contrast).
    #[serde(flatten)]
    pub columns: Table,
    /// How the matrix is coded and what the contrast ids mean, plus counts of contrasts that
    /// cannot be estimated or run.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One default contrast of one stratum: a row of [`Design::columns`].
#[derive(Debug, Clone, PartialEq)]
pub struct DefaultContrast {
    /// The stratum key ([`DesignStratum::stratum`]).
    pub stratum: String,
    /// `c1`, `c2`, ... numbered **within the stratum**: the id of its row in the contrasts file.
    pub contrast_id: String,
    /// `factor=level vs factor=reference`, or a slope's column.
    pub label: String,
    /// `factor=level`; `None` for a slope.
    pub numerator: Option<String>,
    /// `factor=reference`; `None` for a slope.
    pub denominator: Option<String>,
    /// The covariate of a slope; `None` when it is not a slope.
    pub covariate: Option<String>,
    /// The unit a slope is per, e.g. `decade`; `None` when it is not a slope.
    pub covariate_unit: Option<String>,
    /// How the covariate was scaled, e.g. `years/10, centred at 50 years`; `None` when it is not a
    /// slope.
    pub covariate_scale: Option<String>,
    /// One of [`ESTIMABILITIES`].
    pub estimability: String,
    /// One contrast weight per design column (unitless), per [`DesignStratum::column_names`].
    pub weights: Vec<f64>,
    /// One contrast weight per fitted column (unitless), per
    /// [`DesignStratum::fitted_column_names`]: the contrasts file's row. `None` when not
    /// estimable, so no fit gives it an estimate.
    pub fitted_weights: Option<Vec<f64>>,
}

fn weights(cell: Option<Vec<Option<f64>>>) -> Option<Vec<f64>> {
    cell.map(|values| values.into_iter().map(|v| v.unwrap_or(f64::NAN)).collect())
}

impl Design {
    /// One stratum's design, by its key, e.g. `"all"` or `"organism part=liver"`.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if there is no stratum with that key.
    pub fn stratum(&self, key: &str) -> Result<&DesignStratum> {
        self.strata
            .iter()
            .find(|s| s.stratum == key)
            .ok_or_else(|| {
                MzLibError::Usage(format!(
                    "No stratum '{key}'; the strata are: {}.",
                    self.strata
                        .iter()
                        .map(|s| s.stratum.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                ))
            })
    }

    /// Every default contrast, typed, strata in order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn contrasts(&self) -> Result<Vec<DefaultContrast>> {
        let t = &self.columns;
        let stratum = t.strings("stratum")?;
        let contrast_id = t.strings("contrast_id")?;
        let label = t.strings("label")?;
        let numerator = t.strings("numerator")?;
        let denominator = t.strings("denominator")?;
        let covariate = t.strings("covariate")?;
        let covariate_unit = t.strings("covariate_unit")?;
        let covariate_scale = t.strings("covariate_scale")?;
        let estimability = t.strings("estimability")?;
        let all_weights = t.float_arrays("weights")?;
        let fitted_weights = t.float_arrays("fitted_weights")?;
        Ok((0..stratum.len())
            .map(|i| DefaultContrast {
                stratum: stratum[i].clone().unwrap_or_default(),
                contrast_id: contrast_id[i].clone().unwrap_or_default(),
                label: label[i].clone().unwrap_or_default(),
                numerator: numerator[i].clone(),
                denominator: denominator[i].clone(),
                covariate: covariate[i].clone(),
                covariate_unit: covariate_unit[i].clone(),
                covariate_scale: covariate_scale[i].clone(),
                estimability: estimability[i].clone().unwrap_or_default(),
                weights: weights(all_weights[i].clone()).unwrap_or_default(),
                fitted_weights: weights(fitted_weights[i].clone()),
            })
            .collect())
    }
}

// ---- observation results -----------------------------------------------------------------------

/// A MetaMorpheus protein group, as its protein-group table describes it: mzLib's
/// `ProteinGroupInfo`.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct ProteinGroupInfo {
    /// The group's name as MetaMorpheus writes it: member accessions joined by `|`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub name: String,
    /// The group's `Protein QValue`, a fraction (0 to 1), or `None` if not finite.
    #[serde(default)]
    pub q_value: Option<f64>,
    /// The group is a decoy.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub is_decoy: bool,
    /// The group is a contaminant.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub is_contaminant: bool,
    /// The group comes from an entrapment database.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub is_entrapment: bool,
    /// The `Gene` cell as written, or `None` when the table has none.
    #[serde(default)]
    pub genes: Option<String>,
}

/// mzLib's shared-peptide median normalization, as measured: its `NormalizationResult`.
#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
pub struct NormalizationResult {
    /// `"shared_peptide_median"` (measured on the peptides with a value in every sample),
    /// `"shared_peptide_median_half"` (in at least half; used when too few are in all) or `"none"`
    /// (too few either way: no sample shifted).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub setting: String,
    /// Peptides in the reference set used, in peptides, or in the last one tried when
    /// [`Self::setting`] is `"none"`. Decoy, contaminant and entrapment peptides never enter.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub reference_set_size: u64,
    /// The floor used, in peptides: [`ObservationsOptions::min_peptides`], else mzLib's 100.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub minimum_reference_set_size: u64,
    /// Each sample's shift, in log2 units, subtracted with [`ObservationsOptions::normalize`];
    /// `None` for a sample with no shift (not shifted).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub shifts: BTreeMap<String, Option<f64>>,
    /// Which samples have no shift, and why.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub warnings: Vec<String>,
}

/// Where [`observations_with`] wrote the peptide-by-sample responses table
/// ([`ObservationsOptions::responses_out`]).
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct WrittenResponses {
    /// The absolute path written; give it to [`crate::stats::fit_with`] as `responses`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Always `"tsv"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub format: String,
    /// Peptide rows written, in peptides.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Sample columns written, in samples.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
}

/// What [`observations_with`] returns: one row per (peptide, sample), peptides by full sequence.
///
/// [`Self::rows`] gives the rows as typed [`Observation`]s; their columns and units are in the
/// reference below.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct Observations {
    /// Absolute path of the peptide table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub peptides_file: String,
    /// Absolute path of the protein-group table, or `None` without
    /// [`ObservationsOptions::protein_groups`].
    #[serde(default)]
    pub protein_groups_file: Option<String>,
    /// Absolute path of the [`ObservationsOptions::samples`] list, or `None` when every sample was
    /// kept.
    #[serde(default)]
    pub samples_file: Option<String>,
    /// `"msms_only"` or `"mbr_kept"`, as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub basis: String,
    /// Whether the shifts were subtracted.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub normalize: bool,
    /// [`ObservationsOptions::min_peptides`] as given, in peptides, or `None`: mzLib's
    /// `DefaultMinimumPeptides` applied ([`NormalizationResult::minimum_reference_set_size`] says
    /// which).
    #[serde(default)]
    pub min_peptides: Option<u64>,
    /// Runs in the table, in runs.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub run_count: u64,
    /// The runs as mzLib holds them, ordinal by file name, with fraction and technical replicate
    /// filled in.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub runs: Vec<ObservationRun>,
    /// Biological samples, in samples.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// The samples, ordinal by id, with or without values.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub samples: Vec<String>,
    /// Samples in which no peptide has a value under this basis.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub samples_without_values: Vec<String>,
    /// Peptides, in peptides: distinct full sequences.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub peptide_count: u64,
    /// Groups the protein-group table describes, in protein groups.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub protein_group_count: u64,
    /// Those groups.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub protein_groups: Vec<ProteinGroupInfo>,
    /// Rows per state, in cells, keyed by [`OBSERVATION_STATES`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub state_counts: BTreeMap<String, u64>,
    /// The normalization, measured whether or not it was applied.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub normalization: NormalizationResult,
    /// What mzLib noticed that did not stop the build, e.g. groups the peptides name that the
    /// protein-group table does not describe.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub warnings: Vec<String>,
    /// Where [`ObservationsOptions::out`] wrote the table, or `None` when the rows are in
    /// [`Self::columns`].
    #[serde(default)]
    pub written: Option<WrittenTable>,
    /// Where [`ObservationsOptions::responses_out`] wrote, or `None` without it.
    #[serde(default)]
    pub responses_written: Option<WrittenResponses>,
    /// [`Self::peptide_count`] x [`Self::sample_count`], in rows, returned or written.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per (peptide, sample); empty with [`ObservationsOptions::out`].
    #[serde(flatten)]
    pub columns: Table,
    /// How values are combined, what the basis does, what mzLib cannot check about the table, and
    /// what the normalization assumes.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of [`Observations`]: one peptide in one biological sample.
#[derive(Debug, Clone, PartialEq)]
pub struct Observation {
    /// The engine's modified sequence: the key.
    pub full_sequence: String,
    /// The unmodified sequence.
    pub base_sequence: String,
    /// Every group the peptide table names (`UNDEFINED` = in no surviving group).
    pub protein_groups: Vec<String>,
    /// The one group the peptide is unique to; `None` when shared (two or more groups), only
    /// `UNDEFINED`, or none.
    pub unique_protein_group: Option<String>,
    /// The biological sample.
    pub sample: String,
    /// The value, in log2 of FlashLFQ intensity units: fractions summed, replicates averaged in
    /// log2, less the shift with [`ObservationsOptions::normalize`]. `None` = **no value** (see
    /// [`Self::state`]), never 0.
    pub log2_intensity: Option<f64>,
    /// One of [`OBSERVATION_STATES`].
    pub state: String,
}

impl Observations {
    /// Every row, typed, in the table's order: peptides by full sequence, samples within each.
    ///
    /// Empty when the table was written to [`ObservationsOptions::out`] instead.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self) -> Result<Vec<Observation>> {
        let t = &self.columns;
        if t.is_empty() {
            return Ok(Vec::new());
        }
        let full_sequence = t.strings("full_sequence")?;
        let base_sequence = t.strings("base_sequence")?;
        let protein_groups = crate::proteins::string_lists(t, "protein_groups")?;
        let unique_protein_group = t.strings("unique_protein_group")?;
        let sample = t.strings("sample")?;
        let log2_intensity = t.floats("log2_intensity")?;
        let state = t.strings("state")?;
        Ok((0..full_sequence.len())
            .map(|i| Observation {
                full_sequence: full_sequence[i].clone().unwrap_or_default(),
                base_sequence: base_sequence[i].clone().unwrap_or_default(),
                protein_groups: protein_groups[i].clone().unwrap_or_default(),
                unique_protein_group: unique_protein_group[i].clone(),
                sample: sample[i].clone().unwrap_or_default(),
                log2_intensity: log2_intensity[i],
                state: state[i].clone().unwrap_or_default(),
            })
            .collect())
    }
}

// ---- requests ----------------------------------------------------------------------------------

fn path_text(path: &Path, what: &str) -> Result<String> {
    let text = path
        .to_str()
        .ok_or_else(|| MzLibError::Usage(format!("The {what} path is not valid UTF-8.")))?
        .trim();
    if text.is_empty() {
        return Err(MzLibError::Usage(format!("A {what} path is required.")));
    }
    Ok(text.to_owned())
}

/// A name field of a stdin line: non-empty, no tab or line break, sent trimmed.
fn name(value: &str, what: &str) -> Result<String> {
    if value.trim().is_empty() {
        return Err(MzLibError::Usage(format!(
            "{what} must be non-empty; got {value:?}."
        )));
    }
    if value.contains(['\t', '\r', '\n']) {
        return Err(MzLibError::Usage(format!(
            "{what} contains a tab or a line break, which separate the fields: {value:?}."
        )));
    }
    Ok(value.trim().to_owned())
}

fn optional_name(value: Option<&str>, what: &str) -> Result<String> {
    value.map_or_else(|| Ok(String::new()), |v| name(v, what))
}

fn number(value: Option<f64>, what: &str) -> Result<String> {
    value.map_or_else(|| Ok(String::new()), |v| bridge::format_number(v, what))
}

/// A stdin line, with trailing fields left out when they were not given.
fn line(fields: &[String]) -> String {
    let kept = fields
        .iter()
        .rposition(|f| !f.is_empty())
        .map_or(0, |last| last + 1);
    fields[..kept].join("\t")
}

fn model_stdin(options: &DesignOptions) -> Result<String> {
    let mut lines = Vec::new();
    for factor in &options.factors {
        lines.push(line(&[
            "factor".to_owned(),
            name(&factor.name, "A factor name")?,
            optional_name(factor.reference_level.as_deref(), "A reference level")?,
        ]));
    }
    for covariate in &options.covariates {
        let of = |field: &str| format!("The {field} of covariate {:?}", covariate.name);
        lines.push(line(&[
            "covariate".to_owned(),
            name(&covariate.name, "A covariate name")?,
            name(&covariate.unit, &of("unit"))?,
            number(covariate.scale, &of("scale"))?,
            optional_name(covariate.scaled_unit.as_deref(), &of("scaled unit"))?,
            number(covariate.centre, &of("centre"))?,
        ]));
    }
    for (first, second) in &options.interactions {
        lines.push(line(&[
            "interaction".to_owned(),
            name(first, "An interaction's factor")?,
            name(second, "An interaction's factor")?,
        ]));
    }
    for column in &options.strata {
        lines.push(line(&[
            "stratum".to_owned(),
            name(column, "A stratum column")?,
        ]));
    }
    if let Some(batch) = &options.batch {
        lines.push(line(&[
            "batch".to_owned(),
            name(batch, "The batch column")?,
        ]));
    }
    Ok(lines.iter().map(|l| format!("{l}\n")).collect())
}

fn design_request(samples: &Path, options: &DesignOptions) -> Result<(Vec<String>, String)> {
    let mut args = vec![
        "diff".to_owned(),
        "design".to_owned(),
        "--samples".to_owned(),
        path_text(samples, "samples")?,
    ];
    let stdin = model_stdin(options)?;
    if let Some(out) = &options.out {
        args.push("--out".to_owned());
        args.push(path_text(out, "out")?);
    }
    Ok((args, stdin))
}

fn runs_stdin(runs: &[ObservationRun]) -> Result<String> {
    if runs.is_empty() {
        return Err(MzLibError::Usage(
            "observations needs the runs: one ObservationRun per file.".to_owned(),
        ));
    }
    let mut stdin = String::new();
    for run in runs {
        stdin.push_str(&line(&[
            name(&run.file_name, "A run's file name")?,
            name(&run.sample_id, "A run's sample")?,
            run.fraction.map_or_else(String::new, |f| f.to_string()),
            run.technical_replicate
                .map_or_else(String::new, |r| r.to_string()),
        ]));
        stdin.push('\n');
    }
    Ok(stdin)
}

fn observations_request(
    peptides: &Path,
    runs: &[ObservationRun],
    basis: QuantBasis,
    options: &ObservationsOptions,
) -> Result<(Vec<String>, String)> {
    let mut args = vec![
        "diff".to_owned(),
        "observations".to_owned(),
        "--peptides".to_owned(),
        path_text(peptides, "peptides")?,
        "--basis".to_owned(),
        basis.as_str().to_owned(),
    ];
    let stdin = runs_stdin(runs)?;
    for (option, value, what) in [
        (
            "--protein-groups",
            &options.protein_groups,
            "protein_groups",
        ),
        ("--samples", &options.samples, "samples"),
    ] {
        if let Some(path) = value {
            args.push(option.to_owned());
            args.push(path_text(path, what)?);
        }
    }
    if options.normalize {
        args.push("--normalize".to_owned());
    }
    if let Some(minimum) = options.min_peptides {
        args.push("--min-peptides".to_owned());
        args.push(minimum.to_string());
    }
    for (option, value, what) in [
        ("--out", &options.out, "out"),
        ("--responses-out", &options.responses_out, "responses_out"),
    ] {
        if let Some(path) = value {
            args.push(option.to_owned());
            args.push(path_text(path, what)?);
        }
    }
    Ok((args, stdin))
}

fn call<T: serde::de::DeserializeOwned>(
    args: &[String],
    stdin: &str,
    timeout: Option<Duration>,
) -> Result<T> {
    let data = bridge::invoke(args, Some(stdin), timeout)?;
    serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!("diff payload could not be interpreted: {error}"))
    })
}

// ---- the verbs ---------------------------------------------------------------------------------

/// Build each stratum's design matrix and default contrasts from a sample sheet, with mzLib.
///
/// mzLib's `AnalysisDesign.ByStratum` splits the samples by stratum, and in each stratum builds one
/// additive model: an intercept, each factor treatment-coded against its reference level, batch,
/// each covariate as `(value - centre) / scale`, and any interactions. The matrix equals R's
/// `model.matrix` for the same formula. It then lists `DefaultContrasts` — each level against its
/// reference, every pair when there is no reference, each covariate's slope — and judges whether
/// each can be estimated at all from these samples.
///
/// `samples` is the sample sheet: a TSV with a header row and one row per biological sample, first
/// column the sample id; the other columns are available by header. A column the model does not
/// name is ignored ([`Design::unused_columns`]). The model is [`DesignOptions`]'s factors,
/// covariates, interactions, strata and batch; they travel on stdin, one term per line.
///
/// Needs the bridge from pyMzLib 0.6.0 or later.
#[doc = include_str!("../docs/reference/diff.design.md")]
///
/// # Examples
///
/// The MALAT1 dilution series the [`stats`](crate::stats) guide analyses, described by each run's
/// load rather than by a hand-made 0/1 design:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::diff::{design_with, DesignOptions, FactorSpec};
///
/// let design = design_with(
///     "malat_dilution_samples.tsv",
///     &DesignOptions {
///         factors: vec![FactorSpec::new("malat_load", Some("500ng"))],
///         out: Some("malat_design".into()),
///         ..Default::default()
///     },
/// )?;
/// let stratum = design.stratum("all")?;
/// assert_eq!((stratum.sample_count, stratum.rank), (27, 3));
/// assert_eq!(stratum.fitted_column_names, ["intercept", "malat_load=125ng", "malat_load=250ng"]);
/// assert_eq!(design.estimability_counts["estimable"], 2);
/// assert_eq!(design.written.as_ref().map(|w| w.file_count), Some(2));
/// assert!(stratum.contrasts_file.as_deref().is_some_and(|f| f.ends_with("contrasts_1.tsv")));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// Each contrast against the reference is one design column, so [`crate::stats::fit_with`] tests it
/// by that column's name in the design file this call wrote, or as a row of the contrasts file,
/// which adds a confidence interval.
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when a path is blank or not UTF-8, a name in
/// the model is blank or holds a tab or line break, a covariate's scale or centre is not finite, or
/// the bridge predates the verb.
#[doc = include_str!("../docs/reference/diff.design.see-also.md")]
pub fn design_with(samples: impl AsRef<Path>, options: &DesignOptions) -> Result<Design> {
    let (args, stdin) = design_request(samples.as_ref(), options)?;
    bridge::require_verb("diff design", MZLIB_1_0_595_BRIDGE)?;
    call(&args, &stdin, options.timeout)
}

/// [`observations_with`] with every rule mzLib's: no protein groups, every sample, nothing
/// subtracted or written.
///
/// # Errors
///
/// As [`observations_with`].
pub fn observations(
    peptides: impl AsRef<Path>,
    runs: &[ObservationRun],
    basis: QuantBasis,
) -> Result<Observations> {
    observations_with(peptides, runs, basis, &ObservationsOptions::default())
}

/// Every peptide in every biological sample, as log2 intensity or no value with the reason, by
/// mzLib.
///
/// mzLib's `FlashLfqObservations.FromPeptideTable` reads a stored FlashLFQ or MetaMorpheus peptide
/// table, maps each run's detection type to a state, sums each technical replicate's fractions on
/// the linear scale and averages a sample's replicates in log2. Then
/// `SharedPeptideNormalization.Measure` measures how far each sample sits from the others on the
/// peptides they share; [`ObservationsOptions::normalize`] subtracts it.
///
/// `peptides` is a FlashLFQ `QuantifiedPeptides.tsv` or a MetaMorpheus `AllQuantifiedPeptides.tsv`.
/// `runs` lists every run of the table, each an [`ObservationRun`]: the run's name as the
/// `Intensity_` columns write it (without extension), its sample, and optionally its fraction and
/// technical replicate; they travel on stdin. `basis` decides whether match-between-runs values
/// enter. It is required: mzLib has no default.
///
/// Needs the bridge from pyMzLib 0.6.0 or later.
#[doc = include_str!("../docs/reference/diff.observations.md")]
///
/// # Examples
///
/// Two runs of FlashLFQ's E. coli spike-in vignette, which their names say were spiked at 1x and at
/// 2x:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::diff::{observations, ObservationRun, QuantBasis};
///
/// let table = observations(
///     "EcoliSpikeIn_QuantifiedPeptides.tsv",
///     &[
///         ObservationRun::new("09-04-18_EcoliSpikeInSingleShot1x", "spike_1x"),
///         ObservationRun::new("09-04-18_EcoliSpikeInSingleShot2x", "spike_2x"),
///     ],
///     QuantBasis::MsmsOnly,
/// )?;
/// assert_eq!((table.peptide_count, table.samples.clone()), (243, vec!["spike_1x".to_owned(), "spike_2x".to_owned()]));
/// let counts: Vec<u64> = mzlib::diff::OBSERVATION_STATES.iter().map(|s| table.state_counts[*s]).collect();
/// assert_eq!(counts, [319, 43, 9, 6, 109]);
/// assert_eq!((table.normalization.setting.as_str(), table.normalization.reference_set_size), ("none", 86));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// Only 86 peptides were measured in both runs, under mzLib's floor of 100, so mzLib measured no
/// shift rather than trust too few.
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when `runs` is empty, a run's file or sample
/// is blank or holds a tab or line break, a path is blank or not UTF-8, or the bridge predates the
/// verb.
#[doc = include_str!("../docs/reference/diff.observations.see-also.md")]
pub fn observations_with(
    peptides: impl AsRef<Path>,
    runs: &[ObservationRun],
    basis: QuantBasis,
    options: &ObservationsOptions,
) -> Result<Observations> {
    let (args, stdin) = observations_request(peptides.as_ref(), runs, basis, options)?;
    bridge::require_verb("diff observations", MZLIB_1_0_595_BRIDGE)?;
    call(&args, &stdin, options.timeout)
}

#[cfg(test)]
mod tests {
    //! Offline. The recordings are pyMzLib's, made from the real bridge.

    use super::*;

    const MALAT: &str = include_str!("../tests/fixtures/diff_design_malat.json");
    const CONFOUNDED: &str = include_str!("../tests/fixtures/diff_design_confounded.json");
    const STRATA: &str = include_str!("../tests/fixtures/diff_design_strata.json");
    const COVARIATE: &str =
        include_str!("../tests/fixtures/diff_design_factor_batch_covariate.json");
    const MM: &str = include_str!("../tests/fixtures/diff_observations_mm_mbr_kept.json");
    const MM_MSMS: &str = include_str!("../tests/fixtures/diff_observations_mm_msms_only.json");
    const SPIKE: &str =
        include_str!("../tests/fixtures/diff_observations_spike_in_normalized.json");

    #[test]
    fn design_writes_one_stdin_line_per_term_and_leaves_out_what_was_not_given() {
        let options = DesignOptions {
            factors: vec![
                FactorSpec::new("condition", Some("ctrl")),
                FactorSpec::new(" genotype ", None),
            ],
            covariates: vec![
                CovariateSpec::new("age", "years")
                    .scale(10.0)
                    .scaled_unit("decade")
                    .centre(50.0),
                CovariateSpec::new("bmi", "kg/m2"),
                CovariateSpec::new("dose", "mg").centre(1.5),
            ],
            interactions: vec![("condition".into(), "genotype".into())],
            strata: vec!["organism part".into()],
            batch: Some("batch".into()),
            ..Default::default()
        };
        let (args, stdin) = design_request(Path::new("s.tsv"), &options).unwrap();
        assert_eq!(args, ["diff", "design", "--samples", "s.tsv"]);
        assert_eq!(
            stdin,
            "factor\tcondition\tctrl\n\
             factor\tgenotype\n\
             covariate\tage\tyears\t10.0\tdecade\t50.0\n\
             covariate\tbmi\tkg/m2\n\
             covariate\tdose\tmg\t\t\t1.5\n\
             interaction\tcondition\tgenotype\n\
             stratum\torganism part\n\
             batch\tbatch\n"
        );
    }

    #[test]
    fn design_sends_out_only_when_given() {
        let (args, _) = design_request(Path::new("s.tsv"), &DesignOptions::default()).unwrap();
        assert!(!args.iter().any(|a| a == "--out"));
        let (args, _) = design_request(
            Path::new("s.tsv"),
            &DesignOptions {
                out: Some("dir".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(&args[4..], ["--out", "dir"]);
    }

    #[test]
    fn design_refuses_a_malformed_model_before_anything_is_spawned() {
        for options in [
            DesignOptions {
                factors: vec![FactorSpec::new("", None)],
                ..Default::default()
            },
            DesignOptions {
                factors: vec![FactorSpec::new("a\tb", None)],
                ..Default::default()
            },
            DesignOptions {
                covariates: vec![CovariateSpec::new("age", "")],
                ..Default::default()
            },
            DesignOptions {
                covariates: vec![CovariateSpec::new("age", "years").scale(f64::INFINITY)],
                ..Default::default()
            },
            DesignOptions {
                strata: vec!["a\nb".into()],
                ..Default::default()
            },
            DesignOptions {
                batch: Some(" ".into()),
                ..Default::default()
            },
        ] {
            let error = design_request(Path::new("s.tsv"), &options).unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{options:?}");
        }
        assert!(design_request(Path::new(" "), &DesignOptions::default()).is_err());
    }

    #[test]
    fn design_parses_the_malat_recording() {
        let design: Design = serde_json::from_str(MALAT).unwrap();
        assert_eq!(design.sample_count, 27);
        assert_eq!(
            design.factors,
            [FactorSpec::new("malat_load", Some("500ng"))]
        );
        assert!(design.batch_column.is_none() && design.covariates.is_empty());
        let stratum = design.stratum("all").unwrap();
        assert_eq!(stratum.matrix.len(), 27);
        assert_eq!(stratum.fitted_matrix[0].len(), 3);
        assert!(stratum.design_file.is_some());
        let contrasts = design.contrasts().unwrap();
        assert_eq!(contrasts.len(), 2);
        assert_eq!(contrasts[1].weights, [0.0, 0.0, 1.0]);
        assert_eq!(contrasts[0].numerator.as_deref(), Some("malat_load=125ng"));
        assert!(contrasts[0].covariate.is_none());
        assert!(matches!(design.stratum("liver"), Err(MzLibError::Usage(m)) if m.contains("all")));
    }

    #[test]
    fn design_says_why_a_contrast_cannot_be_estimated_or_run() {
        let confounded: Design = serde_json::from_str(CONFOUNDED).unwrap();
        let contrast = &confounded.contrasts().unwrap()[0];
        assert_eq!(contrast.estimability, ESTIMABILITIES[1]);
        assert_eq!(contrast.fitted_weights, None);
        assert_eq!(confounded.strata[0].aliased_terms, ["condition", "batch"]);

        let strata: Design = serde_json::from_str(STRATA).unwrap();
        assert_eq!(strata.stratum_columns, ["organism part"]);
        let kidney = strata.stratum("organism part=kidney").unwrap();
        assert_eq!(kidney.not_run[0].reason, "not_in_stratum");

        let covariate: Design = serde_json::from_str(COVARIATE).unwrap();
        assert_eq!(
            covariate.covariates,
            [CovariateSpec::new("age", "years")
                .scale(10.0)
                .scaled_unit("decade")
                .centre(50.0)]
        );
        assert_eq!(covariate.strata[0].covariate_centres["age"], 50.0);
        assert_eq!(covariate.batch_column.as_deref(), Some("batch"));
    }

    #[test]
    fn observations_write_one_line_per_run_and_send_only_what_was_given() {
        let (args, stdin) = observations_request(
            Path::new("p.tsv"),
            &[
                ObservationRun::new("f1", "s1"),
                ObservationRun::new("f2", "s1").fraction(2),
                ObservationRun::new("f3", "s2").technical_replicate(1),
            ],
            QuantBasis::MsmsOnly,
            &ObservationsOptions::default(),
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "diff",
                "observations",
                "--peptides",
                "p.tsv",
                "--basis",
                "msms_only"
            ]
        );
        assert_eq!(stdin, "f1\ts1\nf2\ts1\t2\nf3\ts2\t\t1\n");
    }

    #[test]
    fn observations_send_every_option_given() {
        let (args, _) = observations_request(
            Path::new("p.tsv"),
            &[ObservationRun::new("f1", "s1")],
            QuantBasis::MbrKept,
            &ObservationsOptions {
                protein_groups: Some("g.tsv".into()),
                samples: Some("design_1.tsv".into()),
                normalize: true,
                min_peptides: Some(50),
                out: Some("o.tsv".into()),
                responses_out: Some("r.tsv".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            &args[4..],
            [
                "--basis",
                "mbr_kept",
                "--protein-groups",
                "g.tsv",
                "--samples",
                "design_1.tsv",
                "--normalize",
                "--min-peptides",
                "50",
                "--out",
                "o.tsv",
                "--responses-out",
                "r.tsv"
            ]
        );
    }

    #[test]
    fn observations_refuse_bad_runs_before_anything_is_spawned() {
        for runs in [
            vec![],
            vec![ObservationRun::new("", "s")],
            vec![ObservationRun::new("f", "s\t1")],
        ] {
            let error = observations_request(
                Path::new("p.tsv"),
                &runs,
                QuantBasis::MsmsOnly,
                &ObservationsOptions::default(),
            )
            .unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{runs:?}");
        }
    }

    #[test]
    fn observations_parse_the_metamorpheus_recording() {
        let table: Observations = serde_json::from_str(MM).unwrap();
        assert_eq!(table.basis, "mbr_kept");
        assert_eq!((table.run_count, table.runs.len()), (18, 18));
        assert_eq!(table.runs[0].fraction, Some(0));
        assert_eq!(table.protein_group_count, 6);
        assert!(table.min_peptides.is_none() && table.written.is_none());
        let rows = table.rows().unwrap();
        assert_eq!(rows.len(), 54);
        assert!(rows
            .iter()
            .all(|r| OBSERVATION_STATES.contains(&r.state.as_str())));
        assert_eq!(rows[0].protein_groups, ["Q16643"]);
        // Missing is never 0: a row without a value says why.
        assert!(rows
            .iter()
            .filter(|r| r.log2_intensity.is_none())
            .all(|r| r.state != "quantified"));
    }

    #[test]
    fn the_two_bases_differ_only_where_a_transfer_is() {
        let kept: Observations = serde_json::from_str(MM).unwrap();
        let only: Observations = serde_json::from_str(MM_MSMS).unwrap();
        let (kept, only) = (kept.rows().unwrap(), only.rows().unwrap());
        for (a, b) in kept.iter().zip(&only) {
            assert_eq!(a.state, b.state);
            if a.log2_intensity != b.log2_intensity {
                assert_eq!(
                    (a.state.as_str(), b.log2_intensity),
                    ("mbr_transferred", None)
                );
            }
        }
    }

    #[test]
    fn the_spike_in_normalization_reports_its_shifts_and_the_responses_file() {
        let shifted: Observations = serde_json::from_str(SPIKE).unwrap();
        assert!(shifted.normalize);
        assert_eq!(shifted.min_peptides, Some(50));
        assert_eq!(shifted.normalization.minimum_reference_set_size, 50);
        let shifts: Vec<f64> = shifted
            .normalization
            .shifts
            .values()
            .map(|s| s.unwrap())
            .collect();
        assert!((shifts[0] + shifts[1]).abs() < 1e-9);
        let written = shifted.responses_written.unwrap();
        assert_eq!((written.format.as_str(), written.row_count), ("tsv", 243));
    }
}
