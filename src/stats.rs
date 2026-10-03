//! Differential abundance with no R: limma's moderated t-test, Benjamini-Hochberg, and
//! random-effects meta-analysis, computed by mzLib.
//!
//! | You want to know | Call | mzLib computes it with |
//! |---|---|---|
//! | Which features change with a coefficient, and how surely? | [`fit_with`] | `LinearModel.Fit` + `EmpiricalBayes.Moderate` |
//! | Which of these p-values survive a false discovery rate? | [`adjust`] | `MultipleTesting.BenjaminiHochberg` |
//! | What is one feature's effect, pooled across studies? | [`meta_with`] | `RandomEffectsMeta.Pool` (DerSimonian-Laird) |
//!
//! Every number on this page is mzLib's (`StatisticalModels`, mzLib #1341 and #1357). This crate
//! sends your tables to the bridge and types what comes back; it does no arithmetic of its own.
//!
//! ## It is limma, checked
//!
//! [`fit_with`] is limma's `lmFit` followed by `eBayes(legacy = TRUE)` (Smyth 2004), and mzLib
//! holds it to limma to 1e-8 relative on limma's own reference output. The crate's test data
//! carries that reference — limma's input, and limma's own `eBayes` output for it — so the
//! comparison is something you can run, not something you have to take on trust:
//!
//! ```
//! # mzlib_replay::activate();
//! # let limma_output = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/stats/limma_ebayes_notrend.tsv");
//! use mzlib::stats::{fit_with, FitOptions};
//!
//! let fit = fit_with(
//!     "limma_reference_responses.tsv",
//!     "limma_reference_design.tsv",
//!     &["age_decades"],
//!     &FitOptions::default(),
//! )?;
//! assert_eq!(fit.feature_count, 400);
//! assert!(fit.residual_df_differ);
//!
//! // limma's own eBayes(legacy = TRUE) output, column t_age.
//! let limma = std::fs::read_to_string(limma_output).unwrap();
//! let limma_t: Vec<f64> = limma
//!     .lines()
//!     .skip(1)
//!     .map(|line| line.split('\t').nth(3).unwrap().parse().unwrap())
//!     .collect();
//! let ours = fit.columns.floats("t")?;
//! let worst = ours
//!     .iter()
//!     .zip(&limma_t)
//!     .map(|(a, b)| (a.unwrap() - b).abs() / b.abs())
//!     .fold(0.0, f64::max);
//! assert!(worst < 1e-8);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Every one of the 400 moderated t-statistics agrees with limma to better than one part in 10^8.
//!
//! **What it is not.** It is limma's *legacy* estimator. When features have different residual
//! degrees of freedom — which omitting missing values causes, so in most label-free data — current
//! limma defaults to a different prior estimator, and [`ModeratedFit::residual_df_differ`] says
//! when your input is such a case. `robust = TRUE`, contrasts, the B-statistic and observation
//! weights are not implemented; write a contrast as its own design column instead.
//!
//! ## A dilution series, end to end
//!
//! mzLib's own RNA test data has MALAT1 loaded at 500, 250 and 125 ng against a constant 500 ng
//! FLuc spike. The responses are log2 intensities normalised to FLuc in each run; the design has
//! an intercept and 0/1 indicators for the two lower loads, so each coefficient is a **log2 fold
//! change against 500 ng**. Halving MALAT1 should read as about -1, quartering it as about -2:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::stats::{fit_with, FitOptions};
//!
//! let fit = fit_with(
//!     "malat_dilution_log2_vs_fluc.tsv",
//!     "malat_dilution_design.tsv",
//!     &["malat_250ng", "malat_125ng"],
//!     &FitOptions::default(),
//! )?;
//! assert_eq!((fit.feature_count, fit.sample_count), (492, 27));
//! assert_eq!(fit.status_counts["fitted"], 212);
//! assert_eq!(fit.status_counts["too_few_observations"], 214);
//! assert_eq!(fit.status_counts["rank_deficient"], 66);
//! assert_eq!(fit.prior.df.map(|df| (df * 10.0).round() / 10.0), Some(19.6));
//!
//! // The MALAT1 oligos' log2 fold change at 250 ng, over the fitted rows only.
//! let rows = fit.rows("malat_250ng")?;
//! let mut malat: Vec<f64> = rows
//!     .iter()
//!     .filter(|r| r.status == "fitted" && r.feature.starts_with("MALAT1:"))
//!     .filter_map(|r| r.estimate)
//!     .collect();
//! malat.sort_by(f64::total_cmp);
//! let median = (malat[malat.len() / 2 - 1] + malat[malat.len() / 2]) / 2.0;
//! assert_eq!((malat.len(), (median * 100.0).round() / 100.0), (116, -0.89));
//!
//! // Which oligos change, at a 5% false discovery rate?
//! let hits = rows.iter().filter(|r| r.bh_adjusted.is_some_and(|q| q < 0.05)).count();
//! assert_eq!(hits, 3);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! What this means for the experiment: the median MALAT1 oligo reads -0.89 log2 at half the load,
//! close to the -1 the dilution predicts, but only 3 oligos pass a 5% FDR at 250 ng — with nine
//! runs per load, most single oligos are too noisy to call a two-fold change on their own. And
//! 280 of the 492 oligos were never tested at all: [`ModeratedFit::status_counts`] reports them
//! (seen in too few runs, or only at some loads) rather than dropping them, and they are not in
//! the BH family.
//!
//! ## Your data, as files
//!
//! [`fit_with`] reads two tab-separated files — a feature-by-sample table of (already
//! log-transformed) responses, and a sample-by-coefficient design — because that is the shape a
//! table library writes, and it keeps a 10,000-protein table off the command line. A blank, `NA` or
//! `NaN` cell is missing and is left out of that feature's fit, never imputed. **A 0 is an
//! observation**, as in mzLib, so blank out any 0 that means "not measured".
//!
//! ## p-values from anywhere else, and several studies
//!
//! [`adjust`] Benjamini-Hochberg adjusts any list of p-values, keeping every position; `None` is an
//! untested feature and is **not counted** in the family. [`meta_with`] pools one effect size per
//! study into a random-effects estimate per feature (DerSimonian and Laird 1986, as metafor's
//! `rma(method = "DL")`), with the two checks a reader asks for first: how many studies agree in
//! direction, and how far the estimate moves when any one study is dropped.
//!
//! ## Cite
//!
//! The methods behind these numbers, from the verbs' specs:
//!
#![doc = include_str!("../docs/reference/cite.stats.md")]

use std::collections::BTreeMap;
use std::path::Path;
use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};
use crate::readers::Table;

/// Every value of the `status` column, as mzLib's `FeatureFitStatus` is written on the wire. Only
/// `"fitted"` rows carry statistics or enter the Benjamini-Hochberg family.
pub const FIT_STATUSES: [&str; 3] = ["fitted", "too_few_observations", "rank_deficient"];

const FIT_TIMEOUT: Duration = Duration::from_secs(600);
const SHORT_TIMEOUT: Duration = Duration::from_secs(60);

// ---- result types ------------------------------------------------------------------------------

/// The fitted prior of the per-feature residual variances, s²_g ~ s0² · χ²(d0)/d0.
///
/// It depends only on the fit, never on which coefficient is tested, so one prior serves every
/// coefficient in a [`ModeratedFit`].
#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
pub struct VariancePrior {
    /// Prior degrees of freedom d0. `None` means **infinite** (see [`Self::df_infinite`]). Larger
    /// means the features' variances agree, so each is shrunk harder toward [`Self::scale`].
    #[serde(default)]
    pub df: Option<f64>,
    /// `true` when the residual variances are no more dispersed than sampling alone explains, so
    /// every feature takes the prior variance.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub df_infinite: bool,
    /// Whether s0² varies with each feature's average response ([`FitOptions::trend`]).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub trended: bool,
    /// Basis functions the trend spline used, intercept included; 1 without a trend.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub spline_basis_count: u32,
    /// s0², the prior variance, in squared response units. `None` when trended: s0² is then per
    /// feature, in the `prior_variance` column.
    #[serde(default)]
    pub scale: Option<f64>,
}

/// What [`fit_with`] returns: one row per (tested coefficient, feature), with limma's moderated t.
///
/// [`Self::columns`] holds one block of [`Self::feature_count`] rows per tested coefficient, in
/// the order you named them, each block in the responses file's feature order. Its columns, and
/// their units, are in the reference below; [`Self::rows`] gives one coefficient's block as typed
/// [`FitRow`]s.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct ModeratedFit {
    /// Absolute path of the responses table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub responses_file: String,
    /// Absolute path of the design table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub design_file: String,
    /// Feature rows read, in features.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub feature_count: u64,
    /// Sample columns read, in samples.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// The samples, in the responses header's order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_names: Vec<String>,
    /// Every design column, in order: the model always fits all of them.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub coefficient_names: Vec<String>,
    /// The coefficients you asked to test, in order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub tested_coefficients: Vec<String>,
    /// Whether the prior variance was allowed to trend with average response.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub trend: bool,
    /// The thread count used, in threads; the answer is identical at any count.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub threads: i64,
    /// Response cells that were blank, NA or NaN, in cells, left out of their features' fits.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub missing_count: u64,
    /// Response cells that are exactly 0, in cells, fitted as observations.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub zero_count: u64,
    /// Features per status, in features, keyed by [`FIT_STATUSES`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub status_counts: BTreeMap<String, u64>,
    /// `true` when the fitted features do not all have the same residual degrees of freedom. Then
    /// current limma (`eBayes(legacy = FALSE)`) would use a different prior estimator and report
    /// different moderated statistics; these are `legacy = TRUE`. A caveat says so.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub residual_df_differ: bool,
    /// The fitted variance prior, shared by every tested coefficient.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub prior: VariancePrior,
    /// Rows in the table, in rows: [`Self::feature_count`] x the number of tested coefficients.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per (tested coefficient, feature), with its column order.
    #[serde(flatten)]
    pub columns: Table,
    /// What the numbers are and are not, including the ones that apply to this input (unfitted
    /// features, zeros, [`Self::residual_df_differ`], an infinite prior df).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of a [`ModeratedFit`]: one feature tested for one coefficient.
///
/// Every statistic is `None` when [`Self::status`] is not `"fitted"`.
#[derive(Debug, Clone, PartialEq)]
pub struct FitRow {
    /// The feature id from the responses file's first column.
    pub feature: String,
    /// The coefficient this row tests.
    pub coefficient: String,
    /// One of [`FIT_STATUSES`].
    pub status: String,
    /// The least-squares coefficient, in response units per coefficient unit: a **log2 fold
    /// change** for a 0/1 group column on log2 intensities. Moderation does not change it.
    pub estimate: Option<f64>,
    /// Moderated standard error, in response units per coefficient unit.
    pub standard_error: Option<f64>,
    /// Moderated t-statistic, in standard errors.
    pub t: Option<f64>,
    /// Degrees of freedom of the moderated t: df_residual + prior df, capped at the pooled df.
    pub df_total: Option<f64>,
    /// Two-sided p-value, a fraction (0 to 1).
    pub p_value: Option<f64>,
    /// Benjamini-Hochberg adjusted p-value over this coefficient's fitted features, a fraction (0
    /// to 1). **Not a target-decoy q-value.**
    pub bh_adjusted: Option<f64>,
    /// The moderated residual variance, in squared response units.
    pub posterior_variance: Option<f64>,
    /// s0² for this feature, in squared response units.
    pub prior_variance: Option<f64>,
    /// Residual standard deviation before moderation, in response units.
    pub sigma: Option<f64>,
    /// Observed samples minus coefficients.
    pub df_residual: Option<i64>,
    /// Samples with a finite response for this feature.
    pub observed: Option<i64>,
    /// Mean of the observed responses, in response units: the covariate a trend fits on. `None`
    /// when the feature has no finite response at all.
    pub average_response: Option<f64>,
}

impl ModeratedFit {
    /// The rows testing one coefficient, typed, in feature order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `coefficient` was not tested in this fit;
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self, coefficient: &str) -> Result<Vec<FitRow>> {
        if !self.tested_coefficients.iter().any(|c| c == coefficient) {
            return Err(MzLibError::Usage(format!(
                "'{coefficient}' was not tested in this fit; it tested: {}.",
                self.tested_coefficients.join(", ")
            )));
        }
        let t = &self.columns;
        let feature = t.strings("feature")?;
        let tested = t.strings("coefficient")?;
        let status = t.strings("status")?;
        let estimate = t.floats("estimate")?;
        let standard_error = t.floats("standard_error")?;
        let t_stat = t.floats("t")?;
        let df_total = t.floats("df_total")?;
        let p_value = t.floats("p_value")?;
        let bh_adjusted = t.floats("bh_adjusted")?;
        let posterior_variance = t.floats("posterior_variance")?;
        let prior_variance = t.floats("prior_variance")?;
        let sigma = t.floats("sigma")?;
        let df_residual = t.integers("df_residual")?;
        let observed = t.integers("observed")?;
        let average_response = t.floats("average_response")?;
        Ok((0..feature.len())
            .filter(|&i| tested[i].as_deref() == Some(coefficient))
            .map(|i| FitRow {
                feature: feature[i].clone().unwrap_or_default(),
                coefficient: coefficient.to_owned(),
                status: status[i].clone().unwrap_or_default(),
                estimate: estimate[i],
                standard_error: standard_error[i],
                t: t_stat[i],
                df_total: df_total[i],
                p_value: p_value[i],
                bh_adjusted: bh_adjusted[i],
                posterior_variance: posterior_variance[i],
                prior_variance: prior_variance[i],
                sigma: sigma[i],
                df_residual: df_residual[i],
                observed: observed[i],
                average_response: average_response[i],
            })
            .collect())
    }
}

/// What [`adjust`] returns: one row per p-value given, in input order.
///
/// Its `p_value` column is each p-value as read, a fraction (0 to 1), `None` where you passed
/// `None` or NaN; `bh_adjusted` is its Benjamini-Hochberg adjusted value, a fraction (0 to 1),
/// `None` for an untested entry.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct Adjusted {
    /// Values given, in lines: the length of every column.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub line_count: u64,
    /// m, the finite p-values the adjustment is over, in p-values.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub tested_count: u64,
    /// Rows in the table, in rows: equal to [`Self::line_count`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, `p_value` and `bh_adjusted`, one row per input.
    #[serde(flatten)]
    pub columns: Table,
    /// What the adjustment is and is not.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl Adjusted {
    /// The adjusted values alone, aligned with the input: a fraction (0 to 1), `None` for an entry
    /// that was not tested.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if the column is not numbers.
    pub fn bh_adjusted(&self) -> Result<Vec<Option<f64>>> {
        self.columns.floats("bh_adjusted")
    }
}

/// What [`meta_with`] returns: one pooled row per feature, in order of first appearance.
///
/// The columns and their units are in the reference below; the estimates, intervals and `tau2`
/// are in the units of the study estimates you sent.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct MetaAnalysis {
    /// Studies read, over every feature, in studies.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub study_count: u64,
    /// Distinct features, in features: the table's row count.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub feature_count: u64,
    /// The interval's coverage, as a fraction (0 to 1).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub confidence: f64,
    /// Rows in the table, in rows: equal to [`Self::feature_count`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per feature.
    #[serde(flatten)]
    pub columns: Table,
    /// What the pooling assumes.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

// ---- options -----------------------------------------------------------------------------------

/// Options for [`fit_with`]. The default is limma's: no trend, one thread.
#[derive(Debug, Clone)]
pub struct FitOptions {
    /// Let the prior variance follow average response (limma's `trend = TRUE`), for data whose
    /// low-abundance features are noisier.
    pub trend: bool,
    /// Basis functions of the trend spline, intercept included (the bridge's default is 4). Only
    /// with `trend`; mzLib may use fewer, and [`VariancePrior::spline_basis_count`] reports what
    /// it used.
    pub spline_basis: Option<u32>,
    /// Features fitted at once, in threads, or `-1` for every core. The answer is identical at any
    /// value.
    pub threads: i32,
    /// How long to allow; `None` waits indefinitely. Ten minutes by default.
    pub timeout: Option<Duration>,
}

impl Default for FitOptions {
    fn default() -> Self {
        Self {
            trend: false,
            spline_basis: None,
            threads: 1,
            timeout: Some(FIT_TIMEOUT),
        }
    }
}

/// One study's effect size for one feature, for [`meta_with`].
#[derive(Debug, Clone, PartialEq)]
pub struct Study {
    /// The feature this study measured. Studies of one feature need not be adjacent.
    pub feature: String,
    /// The study's estimate, in whatever units all studies of the feature share.
    pub estimate: f64,
    /// Its standard error, in the same units; must be finite and positive.
    pub standard_error: f64,
}

impl Study {
    /// One study of `feature`.
    #[must_use]
    pub fn new(feature: impl Into<String>, estimate: f64, standard_error: f64) -> Self {
        Self {
            feature: feature.into(),
            estimate,
            standard_error,
        }
    }
}

/// Options for [`meta_with`].
#[derive(Debug, Clone)]
pub struct MetaOptions {
    /// Coverage of the reported interval, as a fraction (0 to 1). 0.95 by default.
    pub confidence: f64,
    /// How long to allow; `None` waits indefinitely. One minute by default.
    pub timeout: Option<Duration>,
}

impl Default for MetaOptions {
    fn default() -> Self {
        Self {
            confidence: 0.95,
            timeout: Some(SHORT_TIMEOUT),
        }
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

fn one_line(text: &str, what: &str) -> Result<()> {
    if text.contains(['\n', '\r']) {
        return Err(MzLibError::Usage(format!(
            "{what} contains a newline, which separates entries: {text:?}."
        )));
    }
    Ok(())
}

fn fit_request<S: AsRef<str>>(
    responses: &Path,
    design: &Path,
    coefficients: &[S],
    options: &FitOptions,
) -> Result<(Vec<String>, String)> {
    if coefficients.is_empty() {
        return Err(MzLibError::Usage(
            "fit needs at least one coefficient to test.".to_owned(),
        ));
    }
    for name in coefficients {
        let name = name.as_ref();
        if name.trim().is_empty() {
            return Err(MzLibError::Usage(format!(
                "Every coefficient must be a non-empty design column name; got {name:?}."
            )));
        }
        one_line(name, "A coefficient name")?;
    }
    let mut args = vec![
        "stats".to_owned(),
        "fit".to_owned(),
        "--responses".to_owned(),
        path_text(responses, "responses")?,
        "--design".to_owned(),
        path_text(design, "design")?,
        "--threads".to_owned(),
        options.threads.to_string(),
    ];
    if options.trend {
        args.push("--trend".to_owned());
    }
    if let Some(basis) = options.spline_basis {
        args.push("--spline-basis".to_owned());
        args.push(basis.to_string());
    }
    let mut stdin: String = coefficients
        .iter()
        .map(AsRef::as_ref)
        .collect::<Vec<_>>()
        .join("\n");
    stdin.push('\n');
    Ok((args, stdin))
}

fn adjust_stdin(p_values: &[Option<f64>]) -> Result<String> {
    if p_values.is_empty() {
        return Err(MzLibError::Usage(
            "adjust needs at least one p-value.".to_owned(),
        ));
    }
    // A blank line is the wire's "not tested"; a NaN is sent as NaN, which the bridge reads the
    // same way. A trailing newline ends the last line, so a final untested entry survives.
    let mut stdin = String::new();
    for value in p_values {
        if let Some(value) = value {
            stdin.push_str(&value.to_string());
        }
        stdin.push('\n');
    }
    Ok(stdin)
}

fn meta_request(studies: &[Study], options: &MetaOptions) -> Result<(Vec<String>, String)> {
    if studies.is_empty() {
        return Err(MzLibError::Usage(
            "meta needs at least one study.".to_owned(),
        ));
    }
    let mut stdin = String::new();
    for (i, study) in studies.iter().enumerate() {
        if study.feature.trim().is_empty() {
            return Err(MzLibError::Usage(format!(
                "studies[{i}] needs a non-empty feature name."
            )));
        }
        if study.feature.contains('\t') {
            return Err(MzLibError::Usage(format!(
                "studies[{i}]'s feature name contains a tab: {:?}.",
                study.feature
            )));
        }
        one_line(&study.feature, &format!("studies[{i}]'s feature name"))?;
        stdin.push_str(&format!(
            "{}\t{}\t{}\n",
            study.feature, study.estimate, study.standard_error
        ));
    }
    let args = vec![
        "stats".to_owned(),
        "meta".to_owned(),
        "--confidence".to_owned(),
        options.confidence.to_string(),
    ];
    Ok((args, stdin))
}

fn call<T: serde::de::DeserializeOwned>(
    args: &[String],
    stdin: &str,
    timeout: Option<Duration>,
) -> Result<T> {
    let data = bridge::invoke(args, Some(stdin), timeout)?;
    serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!("stats payload could not be interpreted: {error}"))
    })
}

// ---- the verbs ---------------------------------------------------------------------------------

/// [`fit_with`] with every default: no trend, one thread.
///
/// # Errors
///
/// As [`fit_with`].
pub fn fit<S: AsRef<str>>(
    responses: impl AsRef<Path>,
    design: impl AsRef<Path>,
    coefficients: &[S],
) -> Result<ModeratedFit> {
    fit_with(responses, design, coefficients, &FitOptions::default())
}

/// Fit one linear model per feature and test coefficients with limma's moderated t.
///
/// Runs mzLib's `LinearModel.Fit` once — every design column is fitted — and then
/// `EmpiricalBayes.Moderate` for each coefficient you name, which also applies Benjamini-Hochberg
/// across that coefficient's fitted features. This is limma's `lmFit(...)` then
/// `eBayes(..., legacy = TRUE)`, reproduced by mzLib to 1e-8 relative.
///
/// A feature is fitted on the samples where it was observed, and is reported — never dropped —
/// when it cannot be: `too_few_observations` when it has no more observed samples than the design
/// has coefficients, `rank_deficient` when its observed samples cannot separate the coefficients
/// (a feature seen in only one group cannot have a group effect).
///
/// `responses` is a TSV, first column the feature id, every other column one sample; `design` a
/// TSV with one row per sample (matching the responses header by name) and one numeric column per
/// coefficient. `coefficients` are the design columns to test, each with its own moderated test
/// and its own Benjamini-Hochberg family; they travel on stdin.
///
/// Needs the bridge from pyMzLib 0.3.0 or later: an older one is refused before anything is read.
#[doc = include_str!("../docs/reference/stats.fit.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::stats::{fit_with, FitOptions};
///
/// let fit = fit_with(
///     "malat_dilution_log2_vs_fluc.tsv",
///     "malat_dilution_design.tsv",
///     &["malat_250ng", "malat_125ng"],
///     &FitOptions::default(),
/// )?;
/// assert_eq!(fit.tested_coefficients, ["malat_250ng", "malat_125ng"]);
/// assert_eq!(fit.row_count, 984); // 492 features x 2 coefficients
/// let quartered = fit.rows("malat_125ng")?;
/// let calls = quartered.iter().filter(|r| r.bh_adjusted.is_some_and(|q| q < 0.05)).count();
/// assert_eq!(calls, 6);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when `coefficients` is empty or holds a
/// blank name or one with a newline, or when the bridge predates the verb.
#[doc = include_str!("../docs/reference/stats.fit.see-also.md")]
pub fn fit_with<S: AsRef<str>>(
    responses: impl AsRef<Path>,
    design: impl AsRef<Path>,
    coefficients: &[S],
    options: &FitOptions,
) -> Result<ModeratedFit> {
    let (args, stdin) = fit_request(responses.as_ref(), design.as_ref(), coefficients, options)?;
    bridge::require_verb("stats fit", bridge::MZLIB_1_0_593_BRIDGE)?;
    call(&args, &stdin, options.timeout)
}

/// Benjamini-Hochberg adjust a list of p-values, keeping every position.
///
/// For the i-th smallest of m p-values the adjusted value is min over j >= i of p_(j) m / j, capped
/// at 1 (Benjamini and Hochberg 1995), computed by mzLib's `MultipleTesting.BenjaminiHochberg`.
/// `None` (or NaN) marks a feature that was not tested: it stays `None` and is **not counted in
/// m**, which is the honest family. Do not write 1 for an untested feature; that enlarges m.
///
/// [`fit_with`] already reports `bh_adjusted`; this is for p-values from anywhere else. The
/// `p_values` travel on stdin, one per line, each a fraction (0 to 1).
///
/// Needs the bridge from pyMzLib 0.3.0 or later.
#[doc = include_str!("../docs/reference/stats.adjust.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let adjusted = mzlib::stats::adjust(&[
///     Some(0.0002), Some(0.004), Some(0.019), None, Some(0.031), Some(0.2), None, Some(0.74),
/// ])?;
/// assert_eq!(adjusted.tested_count, 6); // the two None are not in the family
/// let rounded: Vec<Option<f64>> = adjusted
///     .bh_adjusted()?
///     .iter()
///     .map(|v| v.map(|v| (v * 1e4).round() / 1e4))
///     .collect();
/// assert_eq!(
///     rounded,
///     [Some(0.0012), Some(0.012), Some(0.038), None, Some(0.0465), Some(0.24), None, Some(0.74)]
/// );
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when `p_values` is empty, or when the bridge
/// predates the verb.
#[doc = include_str!("../docs/reference/stats.adjust.see-also.md")]
pub fn adjust(p_values: &[Option<f64>]) -> Result<Adjusted> {
    let stdin = adjust_stdin(p_values)?;
    bridge::require_verb("stats adjust", bridge::MZLIB_1_0_593_BRIDGE)?;
    call(
        &["stats".to_owned(), "adjust".to_owned()],
        &stdin,
        Some(SHORT_TIMEOUT),
    )
}

/// [`meta_with`] at 95% confidence.
///
/// # Errors
///
/// As [`meta_with`].
pub fn meta(studies: &[Study]) -> Result<MetaAnalysis> {
    meta_with(studies, &MetaOptions::default())
}

/// Pool one effect size per study into a random-effects estimate per feature.
///
/// Each feature is pooled on its own by mzLib's `RandomEffectsMeta.Pool`, with the
/// DerSimonian-Laird (1986) moment estimator of the between-study variance, as metafor's
/// `rma(method = "DL")` does, to 1e-8 relative. Every row also says how many studies agree in
/// direction and how far the estimate moves when any one is dropped — the two checks a reader of a
/// pooled result asks for first. The `studies` travel on stdin, one tab-separated line each.
///
/// Needs the bridge from pyMzLib 0.3.0 or later.
#[doc = include_str!("../docs/reference/stats.meta.md")]
///
/// # Examples
///
/// metafor's own DerSimonian-Laird reference cases, from mzLib's test data:
///
/// ```
/// # mzlib_replay::activate();
/// # let inputs = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/stats/metafor_dl_inputs.tsv");
/// use mzlib::stats::{meta_with, MetaOptions, Study};
///
/// let studies: Vec<Study> = std::fs::read_to_string(inputs)
///     .unwrap()
///     .lines()
///     .skip(1)
///     .map(|line| {
///         let cells: Vec<&str> = line.split('\t').collect();
///         Study::new(cells[0], cells[1].parse().unwrap(), cells[2].parse().unwrap())
///     })
///     .collect();
/// let pooled = meta_with(&studies, &MetaOptions::default())?;
/// assert_eq!(pooled.columns.strings("feature")?[1].as_deref(), Some("heterogeneous"));
/// assert_eq!(pooled.columns.integers("studies")?, [Some(5), Some(6), Some(2), Some(8)]);
/// let tau2 = pooled.columns.floats("tau2")?;
/// assert_eq!(tau2[0], Some(0.0)); // homogeneous: no between-study variance
/// assert_eq!(pooled.columns.integers("direction_disagree")?[1], Some(1));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when `studies` is empty or a feature name is
/// blank or holds a tab or newline, or when the bridge predates the verb.
#[doc = include_str!("../docs/reference/stats.meta.see-also.md")]
pub fn meta_with(studies: &[Study], options: &MetaOptions) -> Result<MetaAnalysis> {
    let (args, stdin) = meta_request(studies, options)?;
    bridge::require_verb("stats meta", bridge::MZLIB_1_0_593_BRIDGE)?;
    call(&args, &stdin, options.timeout)
}

#[cfg(test)]
mod tests {
    //! Offline. The recordings are pyMzLib's, made from the real bridge.

    use super::*;

    const MALAT: &str = include_str!("../tests/fixtures/stats_fit_malat.json");
    const LIMMA: &str = include_str!("../tests/fixtures/stats_fit_limma.json");
    const ADJUST: &str = include_str!("../tests/fixtures/stats_adjust.json");
    const META: &str = include_str!("../tests/fixtures/stats_meta_metafor.json");

    #[test]
    fn a_fit_parses_its_envelope_prior_and_rows() {
        let fit: ModeratedFit = serde_json::from_str(MALAT).unwrap();
        assert_eq!(
            (fit.feature_count, fit.sample_count, fit.row_count),
            (492, 27, 984)
        );
        assert_eq!(fit.status_counts["rank_deficient"], 66);
        assert_eq!((fit.missing_count, fit.zero_count), (10160, 15));
        assert!(fit.residual_df_differ);
        assert!(!fit.prior.trended && !fit.prior.df_infinite);
        assert_eq!(fit.prior.spline_basis_count, 1);
        assert!(fit.prior.scale.is_some());
        assert_eq!(fit.columns.rows(), 984);

        let rows = fit.rows("malat_250ng").unwrap();
        assert_eq!(rows.len(), 492);
        assert!(rows
            .iter()
            .all(|r| FIT_STATUSES.contains(&r.status.as_str())));
        // An unfitted row has no statistics, and is not in the BH family.
        let unfitted = rows.iter().find(|r| r.status != "fitted").unwrap();
        assert_eq!(
            (unfitted.t, unfitted.bh_adjusted, unfitted.df_residual),
            (None, None, None)
        );
        assert!(unfitted.observed.is_some());
    }

    #[test]
    fn a_coefficient_not_tested_is_refused() {
        let fit: ModeratedFit = serde_json::from_str(LIMMA).unwrap();
        let error = fit.rows("sex").unwrap_err();
        assert!(matches!(error, MzLibError::Usage(ref m) if m.contains("age_decades")));
    }

    #[test]
    fn a_fit_sends_its_files_as_arguments_and_its_coefficients_on_stdin() {
        let (args, stdin) = fit_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            &["a", "b"],
            &FitOptions {
                trend: true,
                spline_basis: Some(3),
                threads: -1,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "stats",
                "fit",
                "--responses",
                "r.tsv",
                "--design",
                "d.tsv",
                "--threads",
                "-1",
                "--trend",
                "--spline-basis",
                "3"
            ]
        );
        assert_eq!(stdin, "a\nb\n");

        let (args, _) = fit_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            &["a"],
            &FitOptions::default(),
        )
        .unwrap();
        assert_eq!(args.last().map(String::as_str), Some("1"));
    }

    #[test]
    fn a_fit_with_no_usable_coefficient_is_refused_before_anything_is_spawned() {
        let none: [&str; 0] = [];
        for coefficients in [&none[..], &[" "], &["a\nb"]] {
            let error = fit_request(
                Path::new("r.tsv"),
                Path::new("d.tsv"),
                coefficients,
                &FitOptions::default(),
            )
            .unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{coefficients:?}");
        }
        let error = fit_request(
            Path::new(" "),
            Path::new("d.tsv"),
            &["a"],
            &FitOptions::default(),
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    #[test]
    fn adjust_keeps_every_position_and_a_final_untested_entry() {
        assert_eq!(
            adjust_stdin(&[Some(0.0002), None, Some(0.5), None]).unwrap(),
            "0.0002\n\n0.5\n\n"
        );
        assert!(adjust_stdin(&[]).is_err());

        let adjusted: Adjusted = serde_json::from_str(ADJUST).unwrap();
        assert_eq!((adjusted.line_count, adjusted.tested_count), (8, 6));
        let values = adjusted.bh_adjusted().unwrap();
        assert_eq!(values.len(), 8);
        assert_eq!((values[3], values[6]), (None, None));
    }

    #[test]
    fn meta_sends_one_tab_separated_line_per_study() {
        let (args, stdin) = meta_request(
            &[Study::new("a", 0.3, 0.1), Study::new("a", -0.25, 0.12)],
            &MetaOptions {
                confidence: 0.9,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(args, ["stats", "meta", "--confidence", "0.9"]);
        assert_eq!(stdin, "a\t0.3\t0.1\na\t-0.25\t0.12\n");

        for bad in [Study::new("", 1.0, 1.0), Study::new("a\tb", 1.0, 1.0)] {
            assert!(meta_request(&[bad], &MetaOptions::default()).is_err());
        }
        assert!(meta_request(&[], &MetaOptions::default()).is_err());
        let (args, _) =
            meta_request(&[Study::new("x", 1.0, 1.0)], &MetaOptions::default()).unwrap();
        assert_eq!(args[3], "0.95");
    }

    #[test]
    fn a_meta_analysis_parses_with_one_row_per_feature() {
        let pooled: MetaAnalysis = serde_json::from_str(META).unwrap();
        assert_eq!(
            (pooled.study_count, pooled.feature_count, pooled.row_count),
            (21, 4, 4)
        );
        assert_eq!(pooled.confidence, 0.95);
        assert_eq!(pooled.columns.rows(), 4);
        assert_eq!(pooled.columns.floats("tau2").unwrap()[0], Some(0.0));
        assert_eq!(pooled.caveats.len(), 3);
    }
}
