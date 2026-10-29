//! Differential abundance with no R: limma's moderated t-test with contrasts, Benjamini-Hochberg,
//! random-effects meta-analysis, p-value combination, Spearman correlation, logistic and mixed
//! models, and entrapment FDP estimates, computed by mzLib.
//!
//! | You want to know | Call | mzLib computes it with |
//! |---|---|---|
//! | Which features change with a coefficient, and how surely? | [`fit_with`] | `LinearModel.Fit` + `EmpiricalBayes.Moderate` |
//! | Which features differ by a contrast, within what interval? | [`fit_with`] with [`FitOptions::contrasts`] | `EmpiricalBayes.ModerateContrasts` |
//! | Which of these p-values survive a false discovery rate? | [`adjust`] | `MultipleTesting.BenjaminiHochberg` |
//! | What is one feature's effect, pooled across studies? | [`meta_with`] | `RandomEffectsMeta.Pool` (DerSimonian-Laird) |
//! | What do several studies' p-values say together? | [`combine_p`] | `PValueCombination.Fisher` / `.Stouffer` |
//! | Which features rise or fall with a per-sample trait? | [`correlate_with`] | `SpearmanCorrelation.Correlate` / `.PartialCorrelate` |
//! | Does a 0/1 outcome (detected or not) follow a coefficient? | [`fit_logistic_with`] | `LogisticRegression.Fit` |
//! | Which features change, when samples come in groups? | [`fit_mixed_with`] | `MixedModel.Fit` (a random intercept per group) |
//! | Does a search really control its FDR (entrapment)? | [`entrapment_fdp_with`] | `EntrapmentFdp.Sweep` |
//!
//! Every number on this page is mzLib's (`StatisticalModels`, mzLib #1341, #1357, #1367, #1368,
//! #1385, #1446 and #1452). This crate sends your tables to the bridge and types what comes back; it
//! does no arithmetic of its own.
//!
//! [`fit_with`] already adjusts its own p-values, so you need [`adjust`] only for p-values from
//! somewhere else, and [`meta_with`] only when you have one estimate per study. [`combine_p`] is for
//! when you have only p-values (and the direction of each effect), and [`correlate_with`] for a trait
//! measured on a scale rather than in groups. [`fit_logistic_with`] models a yes/no outcome, such as
//! whether a protein was detected; [`fit_mixed_with`] is [`fit_with`] for samples that share a
//! baseline within groups; and [`entrapment_fdp_with`] checks a search's identifications rather than
//! its quantities.
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
//! **Which limma.** By default it is limma's *legacy* estimator. When features have different
//! residual degrees of freedom — which omitting missing values causes, so in most label-free data —
//! limma 3.61 and later default to a different prior estimator, and
//! [`ModeratedFit::residual_df_differ`] says when your input is such a case. For that case mzLib
//! has an estimator of its own, `"marginal_likelihood"` (below). `robust = TRUE`, the B-statistic
//! and observation weights are not implemented.
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
//! ## Unequal degrees of freedom: the marginal-likelihood prior
//!
//! This experiment is the case the note on limma warns about: oligos are missing from different
//! runs, so they have different residual degrees of freedom, and the default prior is limma's
//! legacy estimator. mzLib's second estimator, `"marginal_likelihood"` (mzLib #1368), fits the
//! prior by maximum likelihood with **each oligo's own degrees of freedom**, so an oligo seen in
//! five runs counts as five runs' worth of evidence about the variance, not as an average one. Ask
//! for it by name, with [`FitOptions::estimator`]:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::stats::{fit_with, FitOptions, ModeratedFit};
//!
//! let (responses, design) = ("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv");
//! let coefficients = ["malat_250ng", "malat_125ng"];
//! let legacy = fit_with(responses, design, &coefficients, &FitOptions::default())?;
//! let marginal = fit_with(
//!     responses,
//!     design,
//!     &coefficients,
//!     &FitOptions {
//!         estimator: Some("marginal_likelihood".into()),
//!         ..Default::default()
//!     },
//! )?;
//!
//! // Which estimator mzLib used, and the prior degrees of freedom it found.
//! let prior = |fit: &ModeratedFit| {
//!     let df = fit.prior.df.map(|df| (df * 10.0).round() / 10.0);
//!     (fit.prior.estimator.clone().unwrap_or_default(), df)
//! };
//! assert_eq!(prior(&legacy), ("moments_legacy".to_owned(), Some(19.6)));
//! assert_eq!(prior(&marginal), ("marginal_likelihood".to_owned(), Some(11.0)));
//!
//! // The 5% FDR calls at half load, by transcript.
//! let calls = |fit: &ModeratedFit| -> mzlib::Result<Vec<String>> {
//!     let mut calls: Vec<String> = fit
//!         .rows("malat_250ng")?
//!         .into_iter()
//!         .filter(|r| r.bh_adjusted.is_some_and(|q| q < 0.05))
//!         .map(|r| r.feature.split(':').next().unwrap_or_default().to_owned())
//!         .collect();
//!     calls.sort();
//!     Ok(calls)
//! };
//! assert_eq!(calls(&legacy)?, ["FLuc", "MALAT1", "MALAT1"]);
//! assert_eq!(calls(&marginal)?, ["MALAT1"]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **The prior df falls from 19.6 to 11.0**: once each oligo's own degrees of freedom count, the
//! data support less shrinkage toward the common variance. The estimates are the same least-squares
//! numbers; only the standard errors move. At half load the FLuc call is gone, which is right, since
//! FLuc did not change, and so is one MALAT1 call. Before you choose it:
//!
//! - **It is mzLib's own estimator, not limma's.** Default limma handles unequal degrees of freedom
//!   with a different method (`fitFDistUnequalDF1`); mzLib's tests hold this one close to it, not
//!   equal.
//! - **mzLib keeps `"moments_legacy"` the default**, so that [`fit_with`] reproduces
//!   `eBayes(legacy = TRUE)` exactly. Leaving [`FitOptions::estimator`] `None` always means mzLib's
//!   default, and [`VariancePrior::estimator`] says which one you got.
//! - **A prior df of about 0.05 is a floor, not a finding.** mzLib searches d0 from e^-3 upward.
//! - One prior serves the whole call: every coefficient and every contrast is moderated with it.
//!
//! ## Contrasts: any comparison, with a confidence interval
//!
//! A design's coefficients answer one kind of question: each level against the baseline. The
//! questions a scientist asks are often other comparisons — 500 ng against 250 ng when the baseline
//! is 125 ng, or the mean of two groups against a third. Each is a **contrast**, a weighted sum of
//! the coefficients, and [`FitOptions::contrasts`] tests every row of a contrasts file with the same
//! moderated t, plus a confidence interval: limma's `contrasts.fit`, then `eBayes`, then
//! `topTable(confint = 0.95)`. [`crate::diff::design_with`] writes the file beside the design. With
//! no reference level, every pair of levels is a contrast:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::diff::{design_with, DesignOptions, FactorSpec};
//! use mzlib::stats::{fit_with, FitOptions};
//!
//! let design = design_with(
//!     "malat_dilution_samples.tsv",
//!     &DesignOptions {
//!         factors: vec![FactorSpec::new("malat_load", None)],
//!         out: Some("malat_pairs".into()),
//!         ..Default::default()
//!     },
//! )?;
//! let stratum = &design.strata[0];
//! assert_eq!(stratum.column_names, ["intercept", "malat_load=250ng", "malat_load=500ng"]);
//! let labels: Vec<String> = design.contrasts()?.into_iter().map(|c| c.label).collect();
//! assert_eq!(labels[2], "malat_load=500ng vs malat_load=250ng");
//!
//! // No coefficient is named, so the coefficient table is empty and every answer is a contrast.
//! let no_coefficients: [&str; 0] = [];
//! let pairs = fit_with(
//!     "malat_dilution_log2_vs_fluc.tsv",
//!     stratum.design_file.as_deref().unwrap_or_default(),
//!     &no_coefficients,
//!     &FitOptions {
//!         contrasts: stratum.contrasts_file.as_deref().map(Into::into),
//!         ..Default::default()
//!     },
//! )?;
//! let contrasts = pairs.contrasts.as_ref().expect("contrasts were asked for");
//! assert_eq!(contrasts.tested_contrasts, ["c1", "c2", "c3"]);
//! assert_eq!(contrasts.confidence_level, 0.95);
//! assert_eq!(contrasts.weights["c3"], [0.0, -1.0, 1.0]);
//! assert_eq!((pairs.row_count, contrasts.row_count), (0, 1476));
//!
//! // The 116 MALAT1 oligos, against the known changes: +1, +2 and +1 in log2.
//! for (id, known, median, holds_it) in [("c1", 1.0, 0.18, 103), ("c2", 2.0, 1.01, 105), ("c3", 1.0, 0.89, 109)] {
//!     let malat: Vec<_> = contrasts
//!         .rows(id)?
//!         .into_iter()
//!         .filter(|r| r.status == "fitted" && r.feature.starts_with("MALAT1:"))
//!         .collect();
//!     let mut estimates: Vec<f64> = malat.iter().filter_map(|r| r.estimate).collect();
//!     estimates.sort_by(f64::total_cmp);
//!     let middle = (estimates[57] + estimates[58]) / 2.0;
//!     let covers = malat
//!         .iter()
//!         .filter(|r| r.confidence_low <= Some(known) && Some(known) <= r.confidence_high)
//!         .count();
//!     assert_eq!(((middle * 100.0).round() / 100.0, covers), (median, holds_it));
//! }
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! The weights are on [`ModeratedFit::coefficient_names`]. `c1` and `c2` are single coefficients;
//! **`c3` is not**: it is the 500 ng coefficient minus the 250 ng one, a comparison no column of
//! this design holds. What the intervals say:
//!
//! - **Each oligo's 95% interval is wide.** It holds the known change for 103 to 109 of the 116
//!   oligos, and holds 0 for most of them too: one oligo, measured in a handful of runs, cannot tell
//!   a doubling from no change. The medians across oligos are far more precise than any one
//!   interval.
//! - **The contrasts locate the compression.** 500 ng against 250 ng reads +0.89 against a known
//!   +1, but 250 ng against 125 ng reads only +0.18 against +1: the shortfall is between the two
//!   lower loads.
//! - **The intervals are per feature and not adjusted.** For a decision across features use
//!   `bh_adjusted`, which each contrast has over its own fitted features.
//!
//! **Your own contrasts.** A contrasts file is a header of a contrast id column and then exactly the
//! design's coefficient columns, in any order, then one row per contrast: its id and one weight per
//! coefficient. To compare the mean of 250 ng and 125 ng against 500 ng on the design above (125 ng
//! is the baseline, so its coefficient is 0 and needs no weight of its own):
//!
//! ```text
//! contrast           intercept  malat_load=250ng  malat_load=500ng
//! low_mean_vs_500    0          0.5               -1
//! ```
//!
//! The traps: a contrast id is the same comparison in every stratum — [`crate::diff::design_with`]
//! numbers them once over the whole analysis, and a stratum that cannot run one skips its id. Read
//! the labels from [`crate::diff::Design::contrasts`]. With missing values a contrast's standard
//! error is exact per feature, where limma's `contrasts.fit` approximates it, so on such data the
//! numbers differ from limma's by design; with complete data they agree. A contrast of one
//! coefficient is that coefficient's test, number for number. [`FitOptions::confidence`] changes
//! the level; leave it `None` for mzLib's default, and [`ModeratedContrasts::confidence_level`]
//! reports the level used.
//!
//! ## Your data, as files
//!
//! [`fit_with`] reads two tab-separated files — a feature-by-sample table of (already
//! log-transformed) responses, and a sample-by-coefficient design — because that is the shape a
//! table library writes, and it keeps a 10,000-protein table off the command line. A blank, `NA` or
//! `NaN` cell is missing and is left out of that feature's fit, never imputed. **A 0 is an
//! observation**, as in mzLib, so blank out any 0 that means "not measured".
//! [`crate::diff::design_with`] writes the design from a sample sheet, and
//! [`crate::diff::observations_with`] the responses from a peptide table.
//!
//! ## p-values from anywhere else, and several studies
//!
//! [`adjust`] Benjamini-Hochberg adjusts any list of p-values, keeping every position; `None` is an
//! untested feature and is **not counted** in the family. [`meta_with`] pools one effect size per
//! study into a random-effects estimate per feature (DerSimonian and Laird 1986, as metafor's
//! `rma(method = "DL")`), with the two checks a reader asks for first: how many studies agree in
//! direction, and how far the estimate moves when any one study is dropped.
//!
//! ## Combining p-values across studies
//!
//! Sometimes all you have from each study is a p-value. [`combine_p`] turns several independent
//! p-values for one feature into one, by the rule you choose — there is no default, because mzLib
//! prefers neither. `"fisher"` (X = -2 sum ln p, against chi-square on 2k degrees of freedom) is
//! driven by the smallest p-values; `"stouffer"` (Z = sum w z / sqrt(sum w²), against the standard
//! normal) weighs every study, by its weight. mzLib's own reference cases show the difference. In
//! `mixed`, one p-value of 0.001 sits among four that disagree:
//!
//! ```
//! # mzlib_replay::activate();
//! # let inputs = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/stats/combination_inputs.tsv");
//! use mzlib::stats::{combine_p, CombinedPValues, PValueTest};
//!
//! // mzLib's reference cases, one (case, p, weight) per line.
//! let rows: Vec<(String, f64, f64)> = std::fs::read_to_string(inputs)
//!     .unwrap()
//!     .lines()
//!     .skip(1)
//!     .map(|line| {
//!         let cells: Vec<&str> = line.split('\t').collect();
//!         (cells[0].to_owned(), cells[1].parse().unwrap(), cells[2].parse().unwrap())
//!     })
//!     .collect();
//! let unweighted: Vec<PValueTest> = rows.iter().map(|(case, p, _)| PValueTest::new(case.as_str(), *p)).collect();
//! let weighted: Vec<PValueTest> = rows
//!     .iter()
//!     .map(|(case, p, weight)| PValueTest::new(case.as_str(), *p).weight(*weight))
//!     .collect();
//! let fisher = combine_p(&unweighted, "fisher")?;
//! let stouffer = combine_p(&weighted, "stouffer")?;
//!
//! let mixed = |combined: &CombinedPValues| -> mzlib::Result<f64> {
//!     let p = combined.p_values()?.into_iter().find(|(case, _)| case == "mixed").unwrap().1;
//!     Ok((p.unwrap() * 1000.0).round() / 1000.0)
//! };
//! assert_eq!((mixed(&fisher)?, mixed(&stouffer)?), (0.019, 0.713));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Fisher says the five studies together are significant, because one of them is. Weighted Stouffer
//! says they are not, because the heaviest study (p = 0.95) points the other way. Both agree with
//! R's `pchisq` and `qnorm` to 1e-8 ([`combine_p`]'s example checks it).
//!
//! **The trap: a two-sided p-value has no direction.** Both rules combine *one-sided* p-values, and
//! the p-value [`fit_with`] reports is two-sided: combine those and two studies that disagree count
//! as two that agree. Give each test its effect ([`PValueTest::effect`]; only the sign is used) and
//! mzLib first turns each p into the one-sided p-value for "effect > 0" (`PValueCombination.OneSided`:
//! p/2 for a positive effect, 1 - p/2 for a negative one), so opposite effects cancel. Three
//! proteins, each tested in cohorts of 16, 25 and 9 people (illustrative numbers, written here, not
//! measured):
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::stats::{combine_p, PValueTest};
//!
//! // (protein, two-sided p, weight = sqrt(cohort size), log2 fold change)
//! let cohorts = [
//!     ("up_everywhere", 0.003, 4.0, 0.8), ("up_everywhere", 0.02, 5.0, 0.5), ("up_everywhere", 0.3, 3.0, 0.3),
//!     ("flips_sign", 0.001, 4.0, 1.2), ("flips_sign", 0.004, 5.0, -0.9), ("flips_sign", 0.4, 3.0, 0.2),
//!     ("no_effect", 0.6, 4.0, -0.1), ("no_effect", 0.45, 5.0, 0.1), ("no_effect", 0.85, 3.0, -0.05),
//! ];
//! let ignored: Vec<PValueTest> = cohorts.iter().map(|&(f, p, w, _)| PValueTest::new(f, p).weight(w)).collect();
//! let signed: Vec<PValueTest> = cohorts.iter().map(|&(f, p, w, e)| PValueTest::new(f, p).weight(w).effect(e)).collect();
//! let ignored = combine_p(&ignored, "stouffer")?;
//! let signed = combine_p(&signed, "stouffer")?;
//! assert!(signed.signed && signed.weighted);
//!
//! let two_figures = |p: Option<f64>| format!("{:.1e}", p.unwrap());
//! for ((protein, ignoring), (_, with_sign)) in ignored.p_values()?.into_iter().zip(signed.p_values()?) {
//!     let row = (protein.as_str(), two_figures(ignoring), two_figures(with_sign));
//!     match protein.as_str() {
//!         "up_everywhere" => assert_eq!(row.1, "6.2e-4"),
//!         "flips_sign" => assert_eq!((row.1.as_str(), row.2.as_str()), ("9.5e-5", "4.3e-1")),
//!         _ => assert_eq!((row.1.as_str(), row.2.as_str()), ("6.9e-1", "4.4e-1")),
//!     }
//! }
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Ignoring direction, `flips_sign` is the strongest result of the three**: it rose sharply in one
//! cohort and fell sharply in another. Signed, opposite effects cancel and it comes out at 0.43. The
//! combined `p_value` is the **upper tail**, so signed it tests effect > 0, and a Stouffer `p_value`
//! near 1 means a consistently *negative* effect; mzLib computes no two-sided combined p-value. The
//! studies must be independent (two analyses of the same samples are one study), the result is not
//! adjusted across features (send it to [`adjust`]), and it is not an effect size — when every study
//! has an estimate and a standard error on one scale, [`meta_with`] pools those; prefer it.
//!
//! ## Correlating features with a trait
//!
//! When the samples differ by a number rather than a group — age, dose, the amount loaded — the
//! question is which features rise or fall with it. [`correlate_with`] answers with Spearman's rank
//! correlation: for each feature, mzLib ranks the feature's values and the trait over the samples
//! where both are known, and tests rho = 0, two-sided. Ranks make it indifferent to a log transform
//! and to one extreme value, and it finds any relationship that only ever rises (or only ever
//! falls). The trait table has one row per sample, the sample name, then one value column whose
//! header names the trait. The dilution series, against each run's MALAT1 load:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::stats::correlate;
//!
//! let result = correlate("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_load_ng.tsv")?;
//! assert_eq!((result.trait_name.as_str(), result.feature_count, result.sample_count), ("malat_ng", 492, 27));
//! assert_eq!(
//!     (result.method_counts["exact"], result.method_counts["asymptotic"], result.method_counts["not_estimable"]),
//!     (4, 308, 180)
//! );
//!
//! let rows = result.rows()?;
//! let median_rho = |transcript: &str| {
//!     let mut rho: Vec<f64> = rows
//!         .iter()
//!         .filter(|r| r.feature.split(':').next() == Some(transcript))
//!         .filter_map(|r| r.rho)
//!         .collect();
//!     rho.sort_by(f64::total_cmp);
//!     let n = rho.len();
//!     let median = if n % 2 == 1 { rho[n / 2] } else { (rho[n / 2 - 1] + rho[n / 2]) / 2.0 };
//!     (n, (median * 100.0).round() / 100.0)
//! };
//! assert_eq!(median_rho("MALAT1"), (194, 0.43));
//! assert_eq!(median_rho("FLuc"), (111, 0.0));
//!
//! // Small n cannot be significant: three samples have six orderings.
//! let exact: Vec<(i64, Option<f64>)> = rows.iter().filter(|r| r.method == "exact").map(|r| (r.n, r.rho)).collect();
//! assert_eq!(exact, [(3, Some(-1.0)), (3, Some(-0.5)), (3, Some(0.5)), (3, Some(0.5))]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! MALAT1's oligos rise with its load, and FLuc's, the spike every run was normalised to, do not.
//! `p_value`'s `method` says how it was computed: `"exact"` (the full permutation distribution, for
//! n of 9 or fewer with no ties, as R's `cor.test`), `"asymptotic"` (t on n - 2 degrees of freedom,
//! as `cor.test(exact = FALSE)`), or `"not_estimable"`, where `rho` and `p_value` are `None`. Even a
//! perfect rank correlation on three samples has p = 1/3.
//!
//! **Partial correlation.** If something else also varies by sample — batch, run order, sex — and
//! could produce the correlation by itself, give it as [`CorrelateOptions::given`], a table of the
//! same shape with one column per covariate. mzLib then correlates what is left of both sets of
//! ranks after each is regressed on the covariates' ranks (`SpearmanCorrelation.PartialCorrelate`,
//! as `ppcor::pcor.test(method = "spearman")`). The partial test is always asymptotic, on n - 2 - k
//! degrees of freedom for k covariates. Partial out only nuisance variables: a covariate that
//! follows the trait takes the trait's signal with it. And for 10 to 1289 samples without ties, R's
//! default `cor.test` uses an Edgeworth series that mzLib does not, so p-values there differ
//! slightly from R's default output (not from `cor.test(exact = FALSE)`).
//!
//! ## Is a feature's detection related to the condition?
//!
//! Missing values in label-free data are rarely random: a protein near the detection limit is
//! missed more often where there is less of it. Before you choose how to handle them, ask the data
//! whether each feature's chance of being detected depends on the condition. That is a logistic
//! regression per feature, on a table of 1 (measured) and 0 (not measured): [`fit_logistic_with`],
//! mzLib's `LogisticRegression.Fit`, which is R's `glm(family = binomial)` per feature, tested with a
//! Wald z. The dilution series' detection table is 1 where the responses file has a value and 0
//! where it is blank:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::stats::fit_logistic;
//!
//! let detected = fit_logistic(
//!     "malat_dilution_detected.tsv",
//!     "malat_dilution_design.tsv",
//!     &["malat_250ng", "malat_125ng"],
//! )?;
//! assert_eq!((detected.status_counts["fitted"], detected.status_counts["separated"]), (191, 301));
//!
//! // MALAT1 at a quarter of its load: how many oligos could be fitted, and what changed.
//! let rows = detected.rows("malat_125ng")?;
//! let malat: Vec<_> = rows.iter().filter(|r| r.feature.starts_with("MALAT1:")).collect();
//! let fitted: Vec<_> = malat.iter().filter(|r| r.status == "fitted").collect();
//! let mut change: Vec<f64> = fitted.iter().filter_map(|r| r.estimate).collect();
//! change.sort_by(f64::total_cmp);
//! let median = change[change.len() / 2]; // 113 values: the middle one
//! let below_5_percent = fitted.iter().filter(|r| r.p_value.is_some_and(|p| p < 0.05)).count();
//! assert_eq!((malat.len(), fitted.len(), (median * 100.0).round() / 100.0, below_5_percent), (345, 113, -0.18, 4));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! `estimate` is a change in the **log-odds** of detection per unit of the coefficient; for a 0/1
//! column it is the log of the odds ratio. Quartering MALAT1 moves its oligos' odds of detection by
//! a median factor of exp(-0.18), about 0.84, and only 4 of 113 pass an unadjusted p of 0.05: here
//! MALAT1's missing values follow its load only weakly. Every oligo has an outcome in every run,
//! because 0 is an outcome here, so none is `too_few_observations`.
//!
//! **A separated feature has no estimate, and that is not the same as no effect.** When the design
//! splits a feature's 0s from its 1s, no maximum-likelihood estimate exists. R's `glm` prints a huge
//! coefficient and a warning there; mzLib reports `separated` and `None`. With three load groups, an
//! oligo seen in a handful of runs is almost always missing from every run of one group, so here
//! separation mostly means "too sparse to say". In a cleaner experiment the same status can mark the
//! strongest result in the table — a protein detected in every treated sample and in no control — so
//! read [`LogisticRow::events`] and [`LogisticRow::observed`] before you discard one. Firth's
//! penalized likelihood (R's `logistf`) is not implemented, and no likelihood-ratio test is offered.
//!
//! ## When samples come in groups
//!
//! Runs acquired on different days, replicates from one donor, studies pooled into one analysis:
//! each group can sit at its own baseline. [`fit_with`] treats every sample as independent, so a
//! day-to-day shift either inflates its residual variance or, when the days are unbalanced, leaks
//! into the effect. A **mixed model** gives every group its own baseline, a random intercept with
//! variance tau², and estimates each coefficient *within* groups: [`fit_mixed_with`], mzLib's
//! `MixedModel.Fit`, which is `nlme::lme(y ~ ..., random = ~ 1 | group)` per feature. A groups table
//! names each sample's group; the dilution series' runs were acquired on five days:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::stats::{fit, fit_mixed};
//!
//! let (responses, design) = ("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv");
//! let coefficients = ["malat_250ng", "malat_125ng"];
//! let mixed = fit_mixed(responses, design, "malat_dilution_day.tsv", &coefficients)?;
//! assert_eq!((mixed.estimator.as_str(), mixed.grouping.as_str(), mixed.group_count), ("reml", "day", 5));
//! assert_eq!((mixed.status_counts["fitted"], mixed.status_counts["too_few_groups"]), (210, 2));
//!
//! // The same MALAT1 oligos at half load, under both models: median estimate and standard error.
//! let median = |mut values: Vec<f64>| {
//!     values.sort_by(f64::total_cmp);
//!     let n = values.len();
//!     ((values[n / 2 - 1] + values[n / 2]) / 2.0 * 100.0).round() / 100.0
//! };
//! let linear = fit(responses, design, &coefficients)?;
//! let malat_linear: Vec<_> = linear.rows("malat_250ng")?.into_iter()
//!     .filter(|r| r.status == "fitted" && r.feature.starts_with("MALAT1:")).collect();
//! let malat_mixed: Vec<_> = mixed.rows("malat_250ng")?.into_iter()
//!     .filter(|r| r.status == "fitted" && r.feature.starts_with("MALAT1:")).collect();
//! assert_eq!((malat_linear.len(), malat_mixed.len()), (116, 116));
//! assert_eq!(median(malat_linear.iter().filter_map(|r| r.estimate).collect()), -0.89);
//! assert_eq!(median(malat_mixed.iter().filter_map(|r| r.estimate).collect()), -1.05);
//! assert_eq!(median(malat_linear.iter().filter_map(|r| r.standard_error).collect()), 1.0);
//! assert_eq!(median(malat_mixed.iter().filter_map(|r| r.standard_error).collect()), 0.71);
//!
//! // How much the days matter: tau² (between days) against sigma² (within a day).
//! assert_eq!(median(malat_mixed.iter().filter_map(|r| r.group_variance).collect()), 0.32);
//! assert_eq!(median(malat_mixed.iter().filter_map(|r| r.residual_variance).collect()), 0.55);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Halving MALAT1 now reads as -1.05 against a known -1, with smaller standard errors**: part of
//! what [`fit_with`] counted as noise was the days. The design is unchanged — do not add a column per
//! day, the random intercept replaces it — and `estimator` is `"reml"` because none was given and
//! that is mzLib's default. Two oligos were seen on only one day, or on no day twice, so their day
//! variance cannot be told from their residual variance: `too_few_groups`. Before you trust a mixed
//! model:
//!
//! - **The df are nlme's containment df**, not Satterthwaite's or Kenward-Roger's, so `lmerTest`
//!   and `pbkrtest` report different p-values. A `None` df means the coefficient was constant within
//!   every group the feature was seen in, and too few groups are left to test it.
//! - **One random intercept, nothing more**: no random slopes, no nested or crossed groupings.
//! - **A feature the fixed effects fit exactly is reported fitted**, with a residual variance near
//!   1e-16 and a t in the millions. Check `residual_variance` before you believe a tiny p-value.
//! - `log_likelihood` is restricted under `"reml"`; pass [`FitMixedOptions::estimator`] `"ml"` to
//!   compare designs.
//!
//! ## Does a search really control its FDR?
//!
//! A search engine reports q-values, and a 1% threshold should mean that about 1% of the
//! discoveries are false. **Entrapment** checks that claim: peptides that cannot be in the sample
//! are added to the search database, so every entrapment the search reports is a known false
//! discovery. [`entrapment_fdp_with`] turns those counts into estimates of the false discovery
//! proportion at each q-value threshold, with mzLib's `EntrapmentFdp.Sweep` (Wen et al. 2025, the
//! method FDRBench implements). mzLib's hand-built test search has seven targets and four
//! entrapments:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::stats::entrapment_fdp;
//!
//! let fdp = entrapment_fdp(
//!     "entrapment_reference_targets.tsv",
//!     "entrapment_reference_entrapments.tsv",
//!     1.0,
//! )?;
//! assert!(!fdp.thresholds_given); // mzLib's default sweep
//! assert_eq!((fdp.threshold_count, &fdp.thresholds[..3]), (100, &[0.001, 0.002, 0.003][..]));
//!
//! let at = fdp.at(0.01)?.expect("discoveries at 1%");
//! assert_eq!((at.target_count, at.entrapment_count), (4, 3));
//! let three = |x: f64| (x * 1000.0).round() / 1000.0;
//! assert_eq!((three(at.lower_bound), three(at.combined), at.paired.map(three)), (0.429, 0.857, Some(0.857)));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! At 1% the search reports seven discoveries, three of them entrapments: far from a 1% FDR. The
//! three estimates answer different questions:
//!
//! - **`lower_bound`**, N_E / (N_T + N_E): above the FDR you claimed, the search **fails** to
//!   control it. Below, it proves nothing.
//! - **`combined`**, N_E (1 + 1/r) / (N_T + N_E), with r the entrapment-to-target database size
//!   ratio: an upper bound. At or below the claimed FDR, it is **evidence of control**.
//! - **`paired`**: a tighter upper bound for 1:1 pairing, reported only when the ratio is exactly 1.
//!   For 1:1 pairing pass 1 itself, not a ratio computed from database sizes after exclusions.
//!
//! None is capped at 1. Give one row per unit you are estimating for (precursor, peptide or protein
//! group), at its best q-value, decoys removed: mzLib does not deduplicate, and does not check the
//! partner columns against the targets file. The "sample" estimator N_E / (r N_T) is deliberately
//! absent: it is invalid in both directions.
//!
//! ## What these functions will not do
//!
//! - **No imputation.** A missing value is left out of its feature's fit.
//! - **Not limma's current default** when residual df differ: `"marginal_likelihood"` is mzLib's
//!   own estimator for that case, close to it but not equal. Not `robust = TRUE`, the B-statistic
//!   or observation weights.
//! - **`bh_adjusted` is not a target-decoy q-value.** It controls the false discovery rate among
//!   the features you tested, a different question from identification FDR.
//! - **[`meta_with`] uses a normal reference** (no Hartung-Knapp adjustment), so its intervals do
//!   not widen for two or three studies.
//! - **[`combine_p`] gives no two-sided combined p-value and no effect size**, and
//!   [`correlate_with`] does not report how many covariates mzLib kept. **[`fit_logistic_with`] has
//!   no Firth correction** and **[`fit_mixed_with`] one random intercept**. None of them adjusts
//!   across features: that is [`adjust`].
//! - **[`entrapment_fdp_with`] estimates; it does not build the entrapment database**, and has no
//!   k-fold paired estimator.
//!
//! ## Cite
//!
//! The methods behind these numbers, from the verbs' specs:
//!
#![doc = include_str!("../docs/reference/cite.stats.md")]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::Deserialize;
use serde_json::Value;

use crate::bridge::{self, MzLibError, Result};
use crate::readers::Table;

/// Every value of the `status` column, as mzLib's `FeatureFitStatus` is written on the wire. Only
/// `"fitted"` rows carry statistics or enter the Benjamini-Hochberg family.
pub const FIT_STATUSES: [&str; 3] = ["fitted", "too_few_observations", "rank_deficient"];

/// The variance-prior estimators [`FitOptions::estimator`] takes: mzLib's `VariancePriorEstimator`
/// `MomentsLegacy` (limma's `eBayes(legacy = TRUE)`) and `MarginalLikelihood` (each feature's own
/// residual df, for unequal df). Leave the option `None` for mzLib's default;
/// [`VariancePrior::estimator`] says which was used.
pub const PRIOR_ESTIMATORS: [&str; 2] = ["moments_legacy", "marginal_likelihood"];

/// The rules [`combine_p`] takes as `method`: mzLib's `PValueCombinationMethod`. There is no
/// default; you choose.
pub const COMBINATION_METHODS: [&str; 2] = ["fisher", "stouffer"];

/// Every value of [`correlate_with`]'s `method` column, as mzLib's `SpearmanPValueMethod` is written
/// on the wire. Only `"not_estimable"` rows have `None` for `rho` and `p_value`.
pub const SPEARMAN_METHODS: [&str; 3] = ["exact", "asymptotic", "not_estimable"];

/// Every value of [`fit_logistic_with`]'s `status` column: the outcomes mzLib's
/// `LogisticRegression.Fit` can give a feature, read from its code (mzLib does not list them per
/// model yet). Only `"fitted"` rows carry statistics.
pub const LOGISTIC_STATUSES: [&str; 5] = [
    "fitted",
    "too_few_observations",
    "rank_deficient",
    "separated",
    "not_converged",
];

/// Every value of [`fit_mixed_with`]'s `status` column: the outcomes mzLib's `MixedModel.Fit` can
/// give a feature, read from its code. Only `"fitted"` rows carry statistics.
pub const MIXED_STATUSES: [&str; 6] = [
    "fitted",
    "too_few_observations",
    "rank_deficient",
    "too_few_groups",
    "variance_ratio_at_limit",
    "not_converged",
];

/// The variance estimators [`FitMixedOptions::estimator`] takes: mzLib's `VarianceEstimator` `Reml`
/// and `MaximumLikelihood`. Leave the option `None` for mzLib's default.
pub const MIXED_ESTIMATORS: [&str; 2] = ["reml", "ml"];

/// The pyMzLib release whose bridge first has `stats combine-p`, `correlate`, `fit-logistic`,
/// `fit-mixed` and `entrapment-fdp`, and whose `stats fit` first reads `--contrasts`, `--confidence`
/// and `--estimator` (mzLib 1.0.595): their specs' `since.pymzlib`.
const MZLIB_1_0_595_BRIDGE: &str = "0.6.0";

const FIT_TIMEOUT: Duration = Duration::from_secs(600);
const SHORT_TIMEOUT: Duration = Duration::from_secs(60);

// ---- result types ------------------------------------------------------------------------------

/// The fitted prior of the per-feature residual variances, s²_g ~ s0² · χ²(d0)/d0.
///
/// It depends only on the fit, never on which coefficient or contrast is tested, so one prior
/// serves every coefficient and every contrast in a [`ModeratedFit`].
#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
pub struct VariancePrior {
    /// How the prior was estimated, one of [`PRIOR_ESTIMATORS`]: [`FitOptions::estimator`], or
    /// mzLib's default (`"moments_legacy"` at the pin) when you left it `None`. `None` only from a
    /// bridge before pyMzLib 0.6.0, which did not report it.
    #[serde(default)]
    pub estimator: Option<String>,
    /// Prior degrees of freedom d0, in degrees of freedom. `None` means **infinite** (see
    /// [`Self::df_infinite`]). Larger means the features' variances agree, so each is shrunk harder
    /// toward [`Self::scale`]. 0 with fewer than 3 features: mzLib then fits no prior (mzLib#1398,
    /// as limma's `squeezeVar` does), so each feature keeps its own variance. Under
    /// `"marginal_likelihood"` mzLib searches d0 from about 0.05 (e^-3) upward, so a d0 there is
    /// that floor, not an optimum.
    #[serde(default)]
    pub df: Option<f64>,
    /// `true` when the residual variances are no more dispersed than sampling alone explains, so
    /// every feature takes the prior variance.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub df_infinite: bool,
    /// Whether s0² varies with each feature's average response ([`FitOptions::trend`]).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub trended: bool,
    /// Basis functions the trend spline used, in basis functions, intercept included; 1 without a
    /// trend.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub spline_basis_count: u32,
    /// s0², the prior variance, in squared response units, when every fitted feature shares one.
    /// `None` when it differs by feature, which is then in the `prior_variance` column
    /// ([`FitRow::prior_variance`]): with a trend, or with fewer than 3 features. With so few,
    /// mzLib fits no prior (mzLib#1398, as limma's `squeezeVar`): [`Self::df`] is 0 and each
    /// feature keeps its own variance.
    ///
    /// Bridges before pyMzLib 0.6.0 reported the first feature's own variance here when there were
    /// fewer than 3 features, as though it were the prior. The pinned bridge
    /// ([`MZLIB_BRIDGE_VERSION`](crate::install::MZLIB_BRIDGE_VERSION)) reports `None`.
    #[serde(default)]
    pub scale: Option<f64>,
}

/// What [`fit_with`] returns: one row per (tested coefficient, feature), with limma's moderated t.
///
/// [`Self::columns`] holds one block of [`Self::feature_count`] rows per tested coefficient, in
/// the order you named them, each block in the responses file's feature order. Its columns, and
/// their units, are in the reference below; [`Self::rows`] gives one coefficient's block as typed
/// [`FitRow`]s. Contrasts, when you asked for them, are a second table, [`Self::contrasts`], so this
/// one is the same with or without them.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct ModeratedFit {
    /// Absolute path of the responses table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub responses_file: String,
    /// Absolute path of the design table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub design_file: String,
    /// Absolute path of the contrasts table, or `None` when [`FitOptions::contrasts`] was `None`.
    #[serde(default)]
    pub contrasts_file: Option<String>,
    /// Feature rows read, in features.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub feature_count: u64,
    /// Sample columns read, in samples.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// The samples, in the responses header's order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_names: Vec<String>,
    /// Every design column, in order: the model always fits all of them, and each contrast's
    /// weights follow this order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub coefficient_names: Vec<String>,
    /// The coefficients you asked to test, in order; empty when you passed only
    /// [`FitOptions::contrasts`].
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
    /// default limma (3.61 and later, `eBayes(legacy = FALSE)`) would use a different prior
    /// estimator and report different moderated statistics; under `"moments_legacy"` these are
    /// `legacy = TRUE`, and `"marginal_likelihood"` ([`FitOptions::estimator`]) is mzLib's
    /// estimator for this case. A caveat says so.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub residual_df_differ: bool,
    /// The fitted variance prior, shared by every tested coefficient and contrast.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub prior: VariancePrior,
    /// Rows in the table, in rows: [`Self::feature_count`] x the number of tested coefficients.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per (tested coefficient, feature), with its column order.
    #[serde(flatten)]
    pub columns: Table,
    /// The tested contrasts, a [`ModeratedContrasts`] with its own table and confidence intervals,
    /// or `None` when [`FitOptions::contrasts`] was `None`.
    #[serde(default)]
    pub contrasts: Option<ModeratedContrasts>,
    /// What the numbers are and are not, including the ones that apply to this input (the
    /// estimator, unfitted features, zeros, [`Self::residual_df_differ`], an infinite prior df,
    /// contrasts).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// The contrasts [`fit_with`] tested: one row per (contrast, feature), with a confidence interval.
///
/// [`ModeratedFit::contrasts`] holds one when you set [`FitOptions::contrasts`]. Each contrast is a
/// row of your contrasts file: a weighted sum of the design's coefficients, such as "500 ng minus
/// 250 ng" (weights -1 and 1 on two level columns). mzLib moderates it with the same prior as the
/// coefficients (`EmpiricalBayes.ModerateContrasts`), as limma's `contrasts.fit`, then `eBayes`,
/// then `topTable(confint = ...)`. [`Self::columns`] holds one block of
/// [`ModeratedFit::feature_count`] rows per contrast, in the file's order, each block in the
/// responses file's feature order; [`Self::rows`] gives one contrast's block as typed
/// [`ContrastRow`]s.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct ModeratedContrasts {
    /// The two-sided level of every interval, as a fraction (0 to 1): [`FitOptions::confidence`],
    /// or mzLib's default (0.95 at the pin) when you left it `None`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub confidence_level: f64,
    /// The contrast ids, in the file's order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub tested_contrasts: Vec<String>,
    /// Contrast id -> its weights, in contrast weight per coefficient (unitless), one per
    /// [`ModeratedFit::coefficient_names`] entry, in that order: what mzLib tested, whatever order
    /// the file's columns were in.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub weights: BTreeMap<String, Vec<f64>>,
    /// Rows in the table, in rows: [`ModeratedFit::feature_count`] x the number of contrasts.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The contrasts table, one row per (contrast, feature), with its column order.
    #[serde(flatten)]
    pub columns: Table,
}

/// One row of a [`ModeratedContrasts`]: one feature tested for one contrast.
///
/// Every statistic is `None` when [`Self::status`] is not `"fitted"`.
#[derive(Debug, Clone, PartialEq)]
pub struct ContrastRow {
    /// The feature id from the responses file's first column.
    pub feature: String,
    /// The contrast id this row tests, from the contrasts file's first column.
    pub contrast: String,
    /// The feature's fit status, one of [`FIT_STATUSES`], as in the coefficient table.
    pub status: String,
    /// The weighted sum of the coefficients, in response units per contrast unit: a **log2 fold
    /// change** for one level against another on log2 intensities. Moderation does not change it.
    pub estimate: Option<f64>,
    /// Moderated standard error, in response units per contrast unit, exact per feature.
    pub standard_error: Option<f64>,
    /// Lower bound of the two-sided [`ModeratedContrasts::confidence_level`] interval, in response
    /// units per contrast unit. Per feature, not adjusted across features.
    pub confidence_low: Option<f64>,
    /// Upper bound of the same interval, in response units per contrast unit.
    pub confidence_high: Option<f64>,
    /// Moderated t-statistic, in standard errors.
    pub t: Option<f64>,
    /// Degrees of freedom of the moderated t: df_residual + prior df, capped at the pooled df.
    pub df_total: Option<f64>,
    /// Two-sided p-value, a fraction (0 to 1).
    pub p_value: Option<f64>,
    /// Benjamini-Hochberg adjusted p-value over this contrast's fitted features, a fraction (0 to
    /// 1). **Not a target-decoy q-value.**
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
    /// Mean of the observed responses, in response units. `None` when the feature has no finite
    /// response at all.
    pub average_response: Option<f64>,
}

impl ModeratedContrasts {
    /// The rows testing one contrast, typed, in feature order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `contrast` was not tested in this fit;
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self, contrast: &str) -> Result<Vec<ContrastRow>> {
        if !self.tested_contrasts.iter().any(|c| c == contrast) {
            return Err(MzLibError::Usage(format!(
                "'{contrast}' was not tested in this fit; it tested the contrasts: {}.",
                self.tested_contrasts.join(", ")
            )));
        }
        let t = &self.columns;
        let feature = t.strings("feature")?;
        let tested = t.strings("contrast")?;
        let status = t.strings("status")?;
        let estimate = t.floats("estimate")?;
        let standard_error = t.floats("standard_error")?;
        let confidence_low = t.floats("confidence_low")?;
        let confidence_high = t.floats("confidence_high")?;
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
            .filter(|&i| tested[i].as_deref() == Some(contrast))
            .map(|i| ContrastRow {
                feature: feature[i].clone().unwrap_or_default(),
                contrast: contrast.to_owned(),
                status: status[i].clone().unwrap_or_default(),
                estimate: estimate[i],
                standard_error: standard_error[i],
                confidence_low: confidence_low[i],
                confidence_high: confidence_high[i],
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
    /// s0² for this feature, in squared response units. With fewer than 3 features it is the
    /// feature's own variance, because no prior is fitted.
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

/// What [`combine_p`] returns: one combined row per feature, in order of first appearance.
///
/// Its columns: `feature`; `studies`, the tests combined, in tests (those with a finite p-value,
/// after the signed conversion when effects were given); `omitted`, this feature's tests that were
/// not combined, in tests (no p-value, or, signed, an effect of 0 or `None`); `statistic`, Fisher's
/// X = -2 sum ln p (a chi-square statistic on 2 x studies degrees of freedom) or Stouffer's
/// Z = sum w z / sqrt(sum w²) (in standard normal deviates); and `p_value`, the upper tail of the
/// statistic, a fraction (0 to 1), so ONE direction (effect > 0 when signed). A `None` statistic
/// with a p-value is an infinite one: `p_value` 0 means plus infinity (a p-value of 0 among the
/// tests), 1 means minus infinity (stouffer, a p-value of 1). A `None` `p_value` means nothing was
/// combined, or stouffer met a p of 0 and a p of 1 together.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct CombinedPValues {
    /// `"fisher"` or `"stouffer"`, as you asked.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub method: String,
    /// Whether the tests carried weights (stouffer only).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub weighted: bool,
    /// Whether the tests carried effects, so every p-value was read as two-sided and converted to
    /// the one-sided p-value for effect > 0 before combining.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub signed: bool,
    /// Tests sent, over every feature, in lines.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub line_count: u64,
    /// Distinct features, in features: the table's row count.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub feature_count: u64,
    /// Rows in the table, in rows: equal to [`Self::feature_count`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per feature.
    #[serde(flatten)]
    pub columns: Table,
    /// What the combination assumes, plus any that apply to this input (omitted tests, infinite
    /// statistics).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl CombinedPValues {
    /// Each feature with its combined p-value, a fraction (0 to 1), in table order. `None` means
    /// nothing could be combined.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn p_values(&self) -> Result<Vec<(String, Option<f64>)>> {
        let features = self.columns.strings("feature")?;
        let p_values = self.columns.floats("p_value")?;
        Ok(features
            .into_iter()
            .map(Option::unwrap_or_default)
            .zip(p_values)
            .collect())
    }
}

/// What [`correlate_with`] returns: one Spearman row per feature, in the responses file's order.
///
/// [`Self::rows`] gives the rows as typed [`Correlation`]s; their columns and units are in the
/// reference below.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct Correlations {
    /// Absolute path of the responses table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub responses_file: String,
    /// Absolute path of the trait table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub trait_file: String,
    /// Absolute path of the covariates table, or `None` without [`CorrelateOptions::given`].
    #[serde(default)]
    pub given_file: Option<String>,
    /// The trait column's header.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub trait_name: String,
    /// The covariate columns' headers, in order; empty without [`CorrelateOptions::given`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub covariate_names: Vec<String>,
    /// `true` when covariates were given, so every row is a partial correlation.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub partial: bool,
    /// Feature rows read, in features: the table's row count.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub feature_count: u64,
    /// Sample columns read, in samples.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// The samples, in the responses header's order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_names: Vec<String>,
    /// Response cells that were blank, NA or NaN, in cells.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub missing_count: u64,
    /// Samples with no trait value, in samples, so in no feature's correlation.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub trait_missing_count: u64,
    /// Samples missing a covariate, in samples, so in no feature's correlation; `None` without
    /// [`CorrelateOptions::given`].
    #[serde(default)]
    pub given_incomplete_count: Option<u64>,
    /// Features per method, in features, keyed by [`SPEARMAN_METHODS`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub method_counts: BTreeMap<String, u64>,
    /// Rows in the table, in rows: equal to [`Self::feature_count`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per feature.
    #[serde(flatten)]
    pub columns: Table,
    /// How the p-values were computed and how they differ from R's defaults, plus a count of
    /// not-estimable features when there are any.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of [`Correlations`]: one feature against the trait.
#[derive(Debug, Clone, PartialEq)]
pub struct Correlation {
    /// The feature id from the responses file's first column.
    pub feature: String,
    /// Spearman's rho, unitless (-1 to 1), ties at their average rank; the partial rho when
    /// covariates were given. `None` when [`Self::method`] is `"not_estimable"`.
    pub rho: Option<f64>,
    /// Two-sided p-value against rho = 0, a fraction (0 to 1). **Not adjusted across features.**
    /// `None` when not estimable.
    pub p_value: Option<f64>,
    /// Complete samples used, in samples: response, trait and every covariate finite.
    pub n: i64,
    /// Whether any variable was tied among those samples (`false` too when there were too few to
    /// rank).
    pub has_ties: bool,
    /// One of [`SPEARMAN_METHODS`]: `"exact"` (n <= 9, no ties, plain), `"asymptotic"` (t on n - 2
    /// df, or n - 2 - k with k covariates kept) or `"not_estimable"`.
    pub method: String,
}

impl Correlations {
    /// Every row, typed, in the responses file's order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self) -> Result<Vec<Correlation>> {
        let t = &self.columns;
        let feature = t.strings("feature")?;
        let rho = t.floats("rho")?;
        let p_value = t.floats("p_value")?;
        let n = t.integers("n")?;
        let has_ties = t.booleans("has_ties")?;
        let method = t.strings("method")?;
        Ok((0..feature.len())
            .map(|i| Correlation {
                feature: feature[i].clone().unwrap_or_default(),
                rho: rho[i],
                p_value: p_value[i],
                n: n[i].unwrap_or_default(),
                has_ties: has_ties[i].unwrap_or_default(),
                method: method[i].clone().unwrap_or_default(),
            })
            .collect())
    }
}

/// The error for a row lookup that names a coefficient a per-coefficient fit did not report.
fn not_reported(coefficient: &str, reported: &[String]) -> MzLibError {
    MzLibError::Usage(format!(
        "'{coefficient}' was not reported by this fit; it reported: {}.",
        reported.join(", ")
    ))
}

/// What [`fit_logistic_with`] returns: one row per (reported coefficient, feature).
///
/// [`Self::columns`] holds one block of [`Self::feature_count`] rows per coefficient you named, in
/// that order, each block in the responses file's feature order; [`Self::rows`] gives one
/// coefficient's block as typed [`LogisticRow`]s.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct LogisticModelFit {
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
    /// The coefficients you asked to report, in order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub tested_coefficients: Vec<String>,
    /// The iteration limit you passed, in iterations, or `None`: mzLib's default was used (100 at
    /// the pin).
    #[serde(default)]
    pub max_iterations: Option<u64>,
    /// The convergence tolerance you passed, as a relative change in deviance (unitless), or
    /// `None`: mzLib's default was used (1e-12 at the pin).
    #[serde(default)]
    pub tolerance: Option<f64>,
    /// The thread count used, in threads; the answer is identical at any count.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub threads: i64,
    /// Response cells that were blank, NA or NaN, in cells, left out of their features' fits.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub missing_count: u64,
    /// Features per status, in features, keyed by [`LOGISTIC_STATUSES`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub status_counts: BTreeMap<String, u64>,
    /// Rows in the table, in rows: [`Self::feature_count`] x the coefficients reported.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per (coefficient, feature), with its column order.
    #[serde(flatten)]
    pub columns: Table,
    /// What the numbers are and are not, with a count of unfitted features by status.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of a [`LogisticModelFit`]: one feature's logistic regression, for one coefficient.
///
/// Every statistic is `None` when [`Self::status`] is not `"fitted"`.
#[derive(Debug, Clone, PartialEq)]
pub struct LogisticRow {
    /// The feature id from the responses file's first column.
    pub feature: String,
    /// The coefficient this row reports.
    pub coefficient: String,
    /// One of [`LOGISTIC_STATUSES`]. `"separated"`: the design splits the 0s from the 1s (or every
    /// outcome is the same), so no finite estimate exists.
    pub status: String,
    /// The maximum-likelihood coefficient, in log-odds per coefficient unit: the **log odds ratio**
    /// for a 0/1 column, so `estimate.exp()` is the odds ratio.
    pub estimate: Option<f64>,
    /// Its standard error, in log-odds per coefficient unit, from the inverse Fisher information.
    pub standard_error: Option<f64>,
    /// The Wald statistic, estimate / standard_error, in standard errors.
    pub z: Option<f64>,
    /// Two-sided Wald p-value, a fraction (0 to 1). **Not adjusted across features.**
    pub p_value: Option<f64>,
    /// The feature's residual deviance (-2 x log-likelihood, natural log), the same in every
    /// coefficient's block.
    pub deviance: Option<f64>,
    /// IRLS iterations used, in iterations; 0 when no fit was attempted.
    pub iterations: i64,
    /// Samples with a 0 or 1, in samples.
    pub observed: i64,
    /// Observed samples whose outcome is 1, in samples.
    pub events: i64,
}

impl LogisticModelFit {
    /// The rows reporting one coefficient, typed, in feature order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `coefficient` was not reported by this fit;
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self, coefficient: &str) -> Result<Vec<LogisticRow>> {
        if !self.tested_coefficients.iter().any(|c| c == coefficient) {
            return Err(not_reported(coefficient, &self.tested_coefficients));
        }
        let t = &self.columns;
        let feature = t.strings("feature")?;
        let reported = t.strings("coefficient")?;
        let status = t.strings("status")?;
        let estimate = t.floats("estimate")?;
        let standard_error = t.floats("standard_error")?;
        let z = t.floats("z")?;
        let p_value = t.floats("p_value")?;
        let deviance = t.floats("deviance")?;
        let iterations = t.integers("iterations")?;
        let observed = t.integers("observed")?;
        let events = t.integers("events")?;
        Ok((0..feature.len())
            .filter(|&i| reported[i].as_deref() == Some(coefficient))
            .map(|i| LogisticRow {
                feature: feature[i].clone().unwrap_or_default(),
                coefficient: coefficient.to_owned(),
                status: status[i].clone().unwrap_or_default(),
                estimate: estimate[i],
                standard_error: standard_error[i],
                z: z[i],
                p_value: p_value[i],
                deviance: deviance[i],
                iterations: iterations[i].unwrap_or_default(),
                observed: observed[i].unwrap_or_default(),
                events: events[i].unwrap_or_default(),
            })
            .collect())
    }
}

/// What [`fit_mixed_with`] returns: one row per (reported coefficient, feature).
///
/// [`Self::columns`] holds one block of [`Self::feature_count`] rows per coefficient you named, in
/// that order, each block in the responses file's feature order; [`Self::rows`] gives one
/// coefficient's block as typed [`MixedRow`]s.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct MixedModelFit {
    /// Absolute path of the responses table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub responses_file: String,
    /// Absolute path of the design table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub design_file: String,
    /// Absolute path of the groups table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub groups_file: String,
    /// The groups file's label column header, e.g. `"day"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub grouping: String,
    /// The distinct labels, in order of first appearance over the samples.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub group_names: Vec<String>,
    /// Distinct labels over all samples, in groups.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub group_count: u64,
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
    /// The coefficients you asked to report, in order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub tested_coefficients: Vec<String>,
    /// `"reml"` or `"ml"`, one of [`MIXED_ESTIMATORS`]: the estimator mzLib used, its default when
    /// [`FitMixedOptions::estimator`] was `None`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub estimator: String,
    /// The thread count used, in threads; the answer is identical at any count.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub threads: i64,
    /// Response cells that were blank, NA or NaN, in cells, left out of their features' fits.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub missing_count: u64,
    /// Response cells that are exactly 0, in cells, fitted as observations.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub zero_count: u64,
    /// Features per status, in features, keyed by [`MIXED_STATUSES`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub status_counts: BTreeMap<String, u64>,
    /// Rows in the table, in rows: [`Self::feature_count`] x the coefficients reported.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per (coefficient, feature), with its column order.
    #[serde(flatten)]
    pub columns: Table,
    /// What the model is and is not, plus counts for this input: unfitted features, rows with no
    /// df, features at the boundary, zero cells.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of a [`MixedModelFit`]: one feature's mixed model, for one coefficient.
///
/// Every statistic is `None` when [`Self::status`] is not `"fitted"`.
#[derive(Debug, Clone, PartialEq)]
pub struct MixedRow {
    /// The feature id from the responses file's first column.
    pub feature: String,
    /// The coefficient this row reports.
    pub coefficient: String,
    /// One of [`MIXED_STATUSES`].
    pub status: String,
    /// The fixed effect within groups, in response units per coefficient unit: a **log2 fold
    /// change** for a 0/1 column on log2 intensities.
    pub estimate: Option<f64>,
    /// Its standard error, in response units per coefficient unit.
    pub standard_error: Option<f64>,
    /// nlme's containment degrees of freedom, a whole number. `None` also when there is no positive
    /// df: a group-level coefficient with too few observed groups to test it.
    pub df: Option<f64>,
    /// The Wald statistic, estimate / standard_error, in standard errors.
    pub t: Option<f64>,
    /// Two-sided t test on [`Self::df`], a fraction (0 to 1). **Not adjusted across features.**
    /// `None` also when `df` is.
    pub p_value: Option<f64>,
    /// sigma², the variance within groups, in squared response units.
    pub residual_variance: Option<f64>,
    /// tau², the variance between the groups' baselines, in squared response units; 0 is the
    /// boundary (no group effect).
    pub group_variance: Option<f64>,
    /// The maximised log-likelihood (natural log), restricted under `"reml"`.
    pub log_likelihood: Option<f64>,
    /// Samples with a finite response, in samples.
    pub observed: i64,
    /// Distinct groups among them, in groups.
    pub groups: i64,
}

impl MixedModelFit {
    /// The rows reporting one coefficient, typed, in feature order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `coefficient` was not reported by this fit;
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self, coefficient: &str) -> Result<Vec<MixedRow>> {
        if !self.tested_coefficients.iter().any(|c| c == coefficient) {
            return Err(not_reported(coefficient, &self.tested_coefficients));
        }
        let t = &self.columns;
        let feature = t.strings("feature")?;
        let reported = t.strings("coefficient")?;
        let status = t.strings("status")?;
        let estimate = t.floats("estimate")?;
        let standard_error = t.floats("standard_error")?;
        let df = t.floats("df")?;
        let t_stat = t.floats("t")?;
        let p_value = t.floats("p_value")?;
        let residual_variance = t.floats("residual_variance")?;
        let group_variance = t.floats("group_variance")?;
        let log_likelihood = t.floats("log_likelihood")?;
        let observed = t.integers("observed")?;
        let groups = t.integers("groups")?;
        Ok((0..feature.len())
            .filter(|&i| reported[i].as_deref() == Some(coefficient))
            .map(|i| MixedRow {
                feature: feature[i].clone().unwrap_or_default(),
                coefficient: coefficient.to_owned(),
                status: status[i].clone().unwrap_or_default(),
                estimate: estimate[i],
                standard_error: standard_error[i],
                df: df[i],
                t: t_stat[i],
                p_value: p_value[i],
                residual_variance: residual_variance[i],
                group_variance: group_variance[i],
                log_likelihood: log_likelihood[i],
                observed: observed[i].unwrap_or_default(),
                groups: groups[i].unwrap_or_default(),
            })
            .collect())
    }
}

/// What [`entrapment_fdp_with`] returns: one row per q-value threshold with a discovery.
///
/// [`Self::rows`] gives the rows as typed [`FdpRow`]s, and [`Self::at`] one threshold's row.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct EntrapmentFdpEstimates {
    /// Absolute path of the targets table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub targets_file: String,
    /// Absolute path of the entrapments table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub entrapments_file: String,
    /// r, the entrapment-to-target database size ratio (unitless), as you gave it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub ratio: f64,
    /// The q-value thresholds estimated at, in order, each a q-value (fraction, 0 to 1): yours, or
    /// mzLib's default sweep.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub thresholds: Vec<f64>,
    /// `true` when you passed [`EntrapmentFdpOptions::thresholds`]; `false` when they are mzLib's
    /// `EntrapmentFdp.DefaultQValueThresholds`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub thresholds_given: bool,
    /// Rows of the targets file, in candidates.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub target_candidate_count: u64,
    /// Rows of the entrapments file, in candidates.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub entrapment_candidate_count: u64,
    /// Entrapments marked foreign, in candidates.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub foreign_count: u64,
    /// Entrapments with a partner, in candidates.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub partnered_count: u64,
    /// The number of thresholds, in thresholds.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub threshold_count: u64,
    /// Thresholds with no discoveries, in thresholds, where every estimate is 0/0, so they have no
    /// row.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub empty_threshold_count: u64,
    /// Rows in the table, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per threshold with a discovery.
    #[serde(flatten)]
    pub columns: Table,
    /// What each estimate can and cannot show, plus any that apply to this input (empty
    /// thresholds, no paired estimate, no partners).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of [`EntrapmentFdpEstimates`]: the estimates at one q-value threshold.
#[derive(Debug, Clone, PartialEq)]
pub struct FdpRow {
    /// The threshold, a q-value (fraction, 0 to 1); a candidate is discovered at or below it.
    pub q_value_threshold: f64,
    /// N_T, targets discovered, in discoveries.
    pub target_count: i64,
    /// N_E, entrapments discovered, in discoveries.
    pub entrapment_count: i64,
    /// N_E / (N_T + N_E), a fraction of discoveries: a **lower bound** on the FDP. Above your
    /// claimed FDR, the search fails to control it; below, it proves nothing.
    pub lower_bound: f64,
    /// N_E (1 + 1/r) / (N_T + N_E), a fraction of discoveries: an **upper bound**. At or below your
    /// claimed FDR, it is evidence of control. Can exceed 1.
    pub combined: f64,
    /// A tighter upper bound for 1:1 pairing, a fraction of discoveries; can exceed 1. `None` when
    /// it does not apply: the ratio is not exactly 1, or every entrapment is foreign.
    pub paired: Option<f64>,
}

impl EntrapmentFdpEstimates {
    /// Every row, typed, in threshold order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self) -> Result<Vec<FdpRow>> {
        let t = &self.columns;
        let threshold = t.floats("q_value_threshold")?;
        let target_count = t.integers("target_count")?;
        let entrapment_count = t.integers("entrapment_count")?;
        let lower_bound = t.floats("lower_bound")?;
        let combined = t.floats("combined")?;
        let paired = t.floats("paired")?;
        Ok((0..threshold.len())
            .map(|i| FdpRow {
                q_value_threshold: threshold[i].unwrap_or(f64::NAN),
                target_count: target_count[i].unwrap_or_default(),
                entrapment_count: entrapment_count[i].unwrap_or_default(),
                lower_bound: lower_bound[i].unwrap_or(f64::NAN),
                combined: combined[i].unwrap_or(f64::NAN),
                paired: paired[i],
            })
            .collect())
    }

    /// The row for one q-value threshold, or `None` when nothing was discovered there.
    ///
    /// `threshold` is one of [`Self::thresholds`], as a fraction (0.01 for 1%), compared exactly.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `threshold` is not one of [`Self::thresholds`];
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    // Exact comparison, as pyMzLib's `threshold in self.thresholds`: both sides are the same
    // decimal read by the same parser.
    #[allow(clippy::float_cmp)]
    pub fn at(&self, threshold: f64) -> Result<Option<FdpRow>> {
        if !self.thresholds.contains(&threshold) {
            return Err(MzLibError::Usage(format!(
                "{threshold} is not one of the thresholds estimated at; pass it in \
                 EntrapmentFdpOptions::thresholds."
            )));
        }
        Ok(self
            .rows()?
            .into_iter()
            .find(|r| r.q_value_threshold == threshold))
    }
}

// ---- options -----------------------------------------------------------------------------------

/// Options for [`fit_with`]. The default is limma's: no trend, one thread, no contrasts, and
/// mzLib's own prior estimator.
#[derive(Debug, Clone)]
pub struct FitOptions {
    /// A TSV of contrasts to test, as [`crate::diff::design_with`] writes it beside the design
    /// (`contrasts_<i>.tsv`): a header of the contrast id column, then exactly the design's
    /// coefficient columns (in any order); then one row per contrast, its id and one weight per
    /// coefficient, e.g. `-1` and `1` for the second level minus the first. Each contrast gets its
    /// own moderated test, Benjamini-Hochberg family and confidence interval, in
    /// [`ModeratedFit::contrasts`]. `None` (the default): no contrasts table. Needs the pyMzLib
    /// 0.6.0 bridge.
    pub contrasts: Option<PathBuf>,
    /// The two-sided level of the contrasts' intervals, as a fraction (0 to 1), strictly between 0
    /// and 1; only with [`Self::contrasts`]. `None` (the default) sends nothing, so mzLib's default
    /// applies (0.95 at the pin); [`ModeratedContrasts::confidence_level`] says which level was
    /// used. Needs the pyMzLib 0.6.0 bridge.
    pub confidence: Option<f64>,
    /// How the variance prior is estimated, one of [`PRIOR_ESTIMATORS`]: `"moments_legacy"`
    /// (limma's `eBayes(legacy = TRUE)`) or `"marginal_likelihood"` (mzLib's own estimator, which
    /// uses each feature's own residual degrees of freedom, so it handles the unequal df that
    /// missing values cause exactly). `None` (the default) sends nothing, so mzLib's default
    /// applies (`"moments_legacy"` at the pin); [`VariancePrior::estimator`] says which was used.
    /// Needs the pyMzLib 0.6.0 bridge.
    pub estimator: Option<String>,
    /// Let the prior variance follow average response (limma's `trend = TRUE`), for data whose
    /// low-abundance features are noisier.
    pub trend: bool,
    /// Basis functions of the trend spline, intercept included. Only with `trend`. Leave it `None`
    /// to let mzLib choose, as limma does, for the number of features: nothing is sent for it.
    /// mzLib may also use fewer than you ask for. Either way
    /// [`VariancePrior::spline_basis_count`] reports what it used.
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
            contrasts: None,
            confidence: None,
            estimator: None,
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

/// One study's test of one feature, for [`combine_p`].
///
/// [`Self::new`] gives a feature and its p-value; [`Self::weight`] and [`Self::effect`] add the
/// optional fields. One stdin line is sent per test: when any test has a weight, every line
/// carries the weight field, blank where a test has none; when any test has an effect, every line
/// carries the effect field too, the same way. The bridge refuses a weight on some tests and not
/// others, and mzLib omits a test whose effect is blank (it has no direction).
#[derive(Debug, Clone, PartialEq)]
pub struct PValueTest {
    /// The feature this test is of. A feature's tests need not be adjacent.
    pub feature: String,
    /// The p-value, a fraction (0 to 1); `None` (or NaN) for a test that could not be run, which
    /// mzLib omits and counts in the `omitted` column. Read as two-sided when the tests carry
    /// effects.
    pub p_value: Option<f64>,
    /// Stouffer only: a positive weight (unitless), for example the square root of the study's
    /// sample size. Give one for every test or for none.
    pub weight: Option<f64>,
    /// The test's estimate, in any unit: only its sign is used. Give one for every test or for
    /// none; a 0 or NaN has no direction, so mzLib omits that test.
    pub effect: Option<f64>,
}

impl PValueTest {
    /// A test of `feature` with this p-value (a fraction, 0 to 1); `None` for one that could not
    /// be run.
    #[must_use]
    pub fn new(feature: impl Into<String>, p_value: impl Into<Option<f64>>) -> Self {
        Self {
            feature: feature.into(),
            p_value: p_value.into(),
            weight: None,
            effect: None,
        }
    }

    /// This test with a Stouffer weight (unitless, positive).
    #[must_use]
    pub fn weight(mut self, weight: f64) -> Self {
        self.weight = Some(weight);
        self
    }

    /// This test with its effect (its estimate, in any unit; only the sign is used), so its
    /// p-value is read as two-sided and converted to the one-sided p-value for effect > 0.
    #[must_use]
    pub fn effect(mut self, effect: f64) -> Self {
        self.effect = Some(effect);
        self
    }
}

/// Options for [`correlate_with`].
#[derive(Debug, Clone)]
pub struct CorrelateOptions {
    /// A TSV of the trait table's shape with one or more covariate columns (batch, run order,
    /// loading) to partial out. `None` (the default): a plain correlation.
    pub given: Option<PathBuf>,
    /// How long to allow; `None` waits indefinitely. Ten minutes by default.
    pub timeout: Option<Duration>,
}

impl Default for CorrelateOptions {
    fn default() -> Self {
        Self {
            given: None,
            timeout: Some(FIT_TIMEOUT),
        }
    }
}

/// Options for [`fit_logistic_with`]. The default leaves every model setting to mzLib.
#[derive(Debug, Clone)]
pub struct FitLogisticOptions {
    /// The IRLS iteration limit per feature, in iterations. `None` (the default) sends nothing, so
    /// mzLib's default applies (100 at the pin); a feature still changing at the limit is
    /// `not_converged`.
    pub max_iterations: Option<u32>,
    /// Stop when the relative change in deviance (unitless) is below this, R `glm`'s rule. `None`
    /// (the default) sends nothing, so mzLib's default applies (1e-12 at the pin, tighter than
    /// `glm`'s 1e-8).
    pub tolerance: Option<f64>,
    /// Features fitted at once, in threads, or `-1` for every core. The answer is identical at any
    /// value. 1 by default, the bridge's default.
    pub threads: i32,
    /// How long to allow; `None` waits indefinitely. Ten minutes by default.
    pub timeout: Option<Duration>,
}

impl Default for FitLogisticOptions {
    fn default() -> Self {
        Self {
            max_iterations: None,
            tolerance: None,
            threads: 1,
            timeout: Some(FIT_TIMEOUT),
        }
    }
}

/// Options for [`fit_mixed_with`]. The default leaves the estimator to mzLib.
#[derive(Debug, Clone)]
pub struct FitMixedOptions {
    /// `"reml"` or `"ml"`, one of [`MIXED_ESTIMATORS`]: how the variances are estimated. `None`
    /// (the default) sends nothing, so mzLib's default applies (`"reml"` at the pin);
    /// [`MixedModelFit::estimator`] says which was used. `"ml"` gives log-likelihoods you can
    /// compare between designs.
    pub estimator: Option<String>,
    /// Features fitted at once, in threads, or `-1` for every core. The answer is identical at any
    /// value. 1 by default, the bridge's default.
    pub threads: i32,
    /// How long to allow; `None` waits indefinitely. Ten minutes by default.
    pub timeout: Option<Duration>,
}

impl Default for FitMixedOptions {
    fn default() -> Self {
        Self {
            estimator: None,
            threads: 1,
            timeout: Some(FIT_TIMEOUT),
        }
    }
}

/// Options for [`entrapment_fdp_with`].
#[derive(Debug, Clone)]
pub struct EntrapmentFdpOptions {
    /// The thresholds to estimate at, each a q-value (fraction, 0 to 1): 0.01 for 1%. `None` (the
    /// default) sends nothing, so mzLib's default sweep applies, echoed in
    /// [`EntrapmentFdpEstimates::thresholds`].
    pub thresholds: Option<Vec<f64>>,
    /// How long to allow; `None` waits indefinitely. Ten minutes by default.
    pub timeout: Option<Duration>,
}

impl Default for EntrapmentFdpOptions {
    fn default() -> Self {
        Self {
            thresholds: None,
            timeout: Some(FIT_TIMEOUT),
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

/// The stdin of a per-coefficient fit: one design column name per line, each ended by a newline,
/// so an empty list is an empty stdin.
fn coefficient_lines<S: AsRef<str>>(coefficients: &[S]) -> Result<String> {
    let mut stdin = String::new();
    for name in coefficients {
        let name = name.as_ref();
        if name.trim().is_empty() {
            return Err(MzLibError::Usage(format!(
                "Every coefficient must be a non-empty design column name; got {name:?}."
            )));
        }
        one_line(name, "A coefficient name")?;
        stdin.push_str(name);
        stdin.push('\n');
    }
    Ok(stdin)
}

/// The stdin of `fit-logistic` and `fit-mixed`, which report at least one coefficient.
fn reported_coefficients<S: AsRef<str>>(coefficients: &[S], caller: &str) -> Result<String> {
    if coefficients.is_empty() {
        return Err(MzLibError::Usage(format!(
            "{caller} needs at least one coefficient to report."
        )));
    }
    coefficient_lines(coefficients)
}

fn fit_request<S: AsRef<str>>(
    responses: &Path,
    design: &Path,
    coefficients: &[S],
    options: &FitOptions,
) -> Result<(Vec<String>, String)> {
    if coefficients.is_empty() && options.contrasts.is_none() {
        return Err(MzLibError::Usage(
            "fit needs at least one coefficient to test, or a contrasts file (FitOptions::contrasts)."
                .to_owned(),
        ));
    }
    let stdin = coefficient_lines(coefficients)?;
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
    if let Some(contrasts) = &options.contrasts {
        args.push("--contrasts".to_owned());
        args.push(path_text(contrasts, "contrasts")?);
    }
    if let Some(confidence) = options.confidence {
        args.push("--confidence".to_owned());
        args.push(bridge::format_number(confidence, "confidence")?);
    }
    if let Some(estimator) = &options.estimator {
        if !PRIOR_ESTIMATORS.contains(&estimator.as_str()) {
            return Err(MzLibError::Usage(format!(
                "estimator must be one of {} (mzLib's VariancePriorEstimator MomentsLegacy or \
                 MarginalLikelihood), or None for mzLib's default; got {estimator:?}.",
                PRIOR_ESTIMATORS.join(", ")
            )));
        }
        args.push("--estimator".to_owned());
        args.push(estimator.clone());
    }
    Ok((args, stdin))
}

/// The error for a `stats fit` option the bridge in use predates, and so ignored. pyMzLib's
/// `_older_fit_bridge`.
fn older_fit_bridge(option: &str) -> MzLibError {
    let at = bridge::bridge_path()
        .map(|path| format!("the bridge at '{}'", path.display()))
        .unwrap_or_else(|_| "the bridge in use".to_owned());
    MzLibError::Usage(format!(
        "FitOptions::{option} needs the bridge from pyMzLib {MZLIB_1_0_595_BRIDGE} or later, and \
         {at} predates it, so it ignored --{option}. Point {} at a newer bridge, or install one \
         with mzlib::install::install_bridge().",
        bridge::BRIDGE_ENV_VAR
    ))
}

/// A bridge before pyMzLib 0.6.0 lists `stats fit` but ignores `--contrasts` and `--estimator`, so
/// the answer is checked for what they add. pyMzLib's check, after its call.
fn check_fit_options_reached(data: &Value, options: &FitOptions) -> Result<()> {
    if options.contrasts.is_some() && data.get("contrasts").map_or(true, Value::is_null) {
        return Err(older_fit_bridge("contrasts"));
    }
    let estimator = data
        .get("prior")
        .and_then(|prior| prior.get("estimator"))
        .filter(|e| !e.is_null());
    if options.estimator.is_some() && estimator.is_none() {
        return Err(older_fit_bridge("estimator"));
    }
    Ok(())
}

/// One numeric field of a stdin line: `None` and NaN are blank (missing), anything else the number
/// as Rust writes it, which the bridge's invariant-culture parser reads back exactly.
fn cell(value: Option<f64>) -> String {
    match value {
        Some(number) if !number.is_nan() => format!("{number:?}"),
        _ => String::new(),
    }
}

fn combine_request(tests: &[PValueTest], method: &str) -> Result<(Vec<String>, String)> {
    if !COMBINATION_METHODS.contains(&method) {
        return Err(MzLibError::Usage(format!(
            "method must be one of {} (mzLib's PValueCombinationMethod); got {method:?}.",
            COMBINATION_METHODS.join(", ")
        )));
    }
    if tests.is_empty() {
        return Err(MzLibError::Usage(
            "combine_p needs at least one test.".to_owned(),
        ));
    }
    // A test has a weight field when any test has a weight or an effect (an effect is the fourth
    // field, so a blank weight comes before it), and an effect field when any test has one: the
    // wire's "on every line or on none".
    let signed = tests.iter().any(|t| t.effect.is_some());
    let weighted = signed || tests.iter().any(|t| t.weight.is_some());
    let mut stdin = String::new();
    for (i, test) in tests.iter().enumerate() {
        if test.feature.trim().is_empty() {
            return Err(MzLibError::Usage(format!(
                "tests[{i}] needs a non-empty feature name."
            )));
        }
        if test.feature.contains('\t') {
            return Err(MzLibError::Usage(format!(
                "tests[{i}]'s feature name contains a tab: {:?}.",
                test.feature
            )));
        }
        one_line(&test.feature, &format!("tests[{i}]'s feature name"))?;
        let mut fields = vec![test.feature.clone(), cell(test.p_value)];
        if weighted {
            fields.push(cell(test.weight));
        }
        if signed {
            fields.push(cell(test.effect));
        }
        stdin.push_str(&fields.join("\t"));
        stdin.push('\n');
    }
    let args = vec![
        "stats".to_owned(),
        "combine-p".to_owned(),
        "--method".to_owned(),
        method.to_owned(),
    ];
    Ok((args, stdin))
}

fn correlate_request(
    responses: &Path,
    trait_table: &Path,
    options: &CorrelateOptions,
) -> Result<Vec<String>> {
    let mut args = vec![
        "stats".to_owned(),
        "correlate".to_owned(),
        "--responses".to_owned(),
        path_text(responses, "responses")?,
        "--trait".to_owned(),
        path_text(trait_table, "trait")?,
    ];
    if let Some(given) = &options.given {
        args.push("--given".to_owned());
        args.push(path_text(given, "given")?);
    }
    Ok(args)
}

fn logistic_request<S: AsRef<str>>(
    responses: &Path,
    design: &Path,
    coefficients: &[S],
    options: &FitLogisticOptions,
) -> Result<(Vec<String>, String)> {
    let stdin = reported_coefficients(coefficients, "fit_logistic")?;
    let mut args = vec![
        "stats".to_owned(),
        "fit-logistic".to_owned(),
        "--responses".to_owned(),
        path_text(responses, "responses")?,
        "--design".to_owned(),
        path_text(design, "design")?,
        "--threads".to_owned(),
        options.threads.to_string(),
    ];
    if let Some(limit) = options.max_iterations {
        args.push("--max-iterations".to_owned());
        args.push(limit.to_string());
    }
    if let Some(tolerance) = options.tolerance {
        args.push("--tolerance".to_owned());
        args.push(bridge::format_number(tolerance, "tolerance")?);
    }
    Ok((args, stdin))
}

fn mixed_request<S: AsRef<str>>(
    responses: &Path,
    design: &Path,
    groups: &Path,
    coefficients: &[S],
    options: &FitMixedOptions,
) -> Result<(Vec<String>, String)> {
    let stdin = reported_coefficients(coefficients, "fit_mixed")?;
    let mut args = vec![
        "stats".to_owned(),
        "fit-mixed".to_owned(),
        "--responses".to_owned(),
        path_text(responses, "responses")?,
        "--design".to_owned(),
        path_text(design, "design")?,
        "--groups".to_owned(),
        path_text(groups, "groups")?,
        "--threads".to_owned(),
        options.threads.to_string(),
    ];
    if let Some(estimator) = &options.estimator {
        if !MIXED_ESTIMATORS.contains(&estimator.as_str()) {
            return Err(MzLibError::Usage(format!(
                "estimator must be one of {} (mzLib's VarianceEstimator Reml or \
                 MaximumLikelihood), or None for mzLib's default; got {estimator:?}.",
                MIXED_ESTIMATORS.join(", ")
            )));
        }
        args.push("--estimator".to_owned());
        args.push(estimator.clone());
    }
    Ok((args, stdin))
}

fn entrapment_request(
    targets: &Path,
    entrapments: &Path,
    ratio: f64,
    options: &EntrapmentFdpOptions,
) -> Result<Vec<String>> {
    let mut args = vec![
        "stats".to_owned(),
        "entrapment-fdp".to_owned(),
        "--targets".to_owned(),
        path_text(targets, "targets")?,
        "--entrapments".to_owned(),
        path_text(entrapments, "entrapments")?,
        "--ratio".to_owned(),
        bridge::format_number(ratio, "ratio")?,
    ];
    if let Some(thresholds) = &options.thresholds {
        if thresholds.is_empty() {
            return Err(MzLibError::Usage(
                "thresholds takes a non-empty list of q-values, or None for mzLib's default."
                    .to_owned(),
            ));
        }
        let text = thresholds
            .iter()
            .enumerate()
            .map(|(i, t)| bridge::format_number(*t, &format!("thresholds[{i}]")))
            .collect::<Result<Vec<_>>>()?
            .join(",");
        args.push("--thresholds".to_owned());
        args.push(text);
    }
    Ok(args)
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
    parse(bridge::invoke(args, Some(stdin), timeout)?)
}

fn parse<T: serde::de::DeserializeOwned>(data: Value) -> Result<T> {
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

/// Fit one linear model per feature and test coefficients, and contrasts, with limma's moderated t.
///
/// Runs mzLib's `LinearModel.Fit` once — every design column is fitted — and then
/// `EmpiricalBayes.Moderate` for each coefficient you name, which also applies Benjamini-Hochberg
/// across that coefficient's fitted features. This is limma's `lmFit(...)` then
/// `eBayes(..., legacy = TRUE)`, reproduced by mzLib to 1e-8 relative. With
/// [`FitOptions::contrasts`], `EmpiricalBayes.ModerateContrasts` also tests each row of a contrasts
/// file, with its confidence interval: limma's `contrasts.fit`, then `eBayes`, then
/// `topTable(confint = ...)`. [`FitOptions::estimator`] chooses how the variance prior is
/// estimated; `None` is mzLib's default.
///
/// A feature is fitted on the samples where it was observed, and is reported — never dropped —
/// when it cannot be: `too_few_observations` when it has no more observed samples than the design
/// has coefficients, `rank_deficient` when its observed samples cannot separate the coefficients
/// (a feature seen in only one group cannot have a group effect).
///
/// `responses` is a TSV, first column the feature id, every other column one sample; `design` a
/// TSV with one row per sample (matching the responses header by name) and one numeric column per
/// coefficient; [`crate::diff::design_with`] writes one (`design_<i>.tsv`) from a sample sheet.
/// `coefficients` are the design columns to test, each with its own moderated test and its own
/// Benjamini-Hochberg family; they travel on stdin. It may be empty when you give
/// [`FitOptions::contrasts`].
///
/// Needs the bridge from pyMzLib 0.3.0 or later: an older one is refused before anything is read.
/// [`FitOptions::contrasts`], [`FitOptions::confidence`] and [`FitOptions::estimator`] need the
/// pyMzLib 0.6.0 bridge: an older one lists the verb but ignores them, so the answer is checked for
/// what they add, and a fit that ignored them is refused rather than returned.
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
/// assert_eq!(fit.prior.estimator.as_deref(), Some("moments_legacy")); // mzLib's default
/// let quartered = fit.rows("malat_125ng")?;
/// let calls = quartered.iter().filter(|r| r.bh_adjusted.is_some_and(|q| q < 0.05)).count();
/// assert_eq!(calls, 6);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// Contrasts, and the module guide's worked examples of both options, are under
/// [Contrasts](crate::stats#contrasts-any-comparison-with-a-confidence-interval).
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when `coefficients` is empty and there is no
/// [`FitOptions::contrasts`], holds a blank name or one with a newline, when
/// [`FitOptions::estimator`] is not one of [`PRIOR_ESTIMATORS`] or [`FitOptions::confidence`] is not
/// finite, or when the bridge predates the verb; and after the call, when the bridge predates
/// pyMzLib 0.6.0 and so ignored [`FitOptions::contrasts`] or [`FitOptions::estimator`].
#[doc = include_str!("../docs/reference/stats.fit.see-also.md")]
pub fn fit_with<S: AsRef<str>>(
    responses: impl AsRef<Path>,
    design: impl AsRef<Path>,
    coefficients: &[S],
    options: &FitOptions,
) -> Result<ModeratedFit> {
    let (args, stdin) = fit_request(responses.as_ref(), design.as_ref(), coefficients, options)?;
    bridge::require_verb("stats fit", bridge::MZLIB_1_0_593_BRIDGE)?;
    let data = match bridge::invoke(&args, Some(&stdin), options.timeout) {
        // A bridge from before 0.6.0 ignores --contrasts and so, with no coefficient, says it was
        // given none: the option it did not read is the real cause.
        Err(MzLibError::Usage(message))
            if options.contrasts.is_some()
                && coefficients.is_empty()
                && message.contains("No coefficient to test") =>
        {
            return Err(older_fit_bridge("contrasts"));
        }
        other => other?,
    };
    check_fit_options_reached(&data, options)?;
    parse(data)
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

/// Combine several independent p-values per feature into one, by Fisher's or Stouffer's method.
///
/// Each feature is combined on its own by mzLib's `PValueCombination`: Fisher (1925),
/// X = -2 sum ln p against chi-square on 2k degrees of freedom, or Stouffer (1949),
/// Z = sum w z / sqrt(sum w²) against the standard normal. mzLib holds both to R's `pchisq` and
/// `qnorm` to 1e-8 relative. Use it when you have one p-value per study — per dataset, cohort or
/// laboratory — for the same feature, and want one answer from all of them.
///
/// `method` is `"fisher"` or `"stouffer"`, one of [`COMBINATION_METHODS`]. It is required: mzLib
/// prefers neither rule. Fisher is driven by the strongest single result and takes no weights;
/// Stouffer weighs every study, and takes [`PValueTest::weight`]. The `tests` travel on stdin, one
/// tab-separated line each.
///
/// **Direction.** Both rules combine *one-sided* p-values. A two-sided p-value from [`fit_with`]
/// does not say which way the effect went, so two studies that disagree look like two that agree.
/// Give every test its [`PValueTest::effect`] and mzLib first turns each p into the one-sided
/// p-value for effect > 0, so opposite effects cancel; the combined `p_value` then tests effect > 0.
/// Neither rule estimates an effect size ([`meta_with`] does) or adjusts across features
/// ([`adjust`] does).
///
/// Needs the bridge from pyMzLib 0.6.0 or later.
#[doc = include_str!("../docs/reference/stats.combine-p.md")]
///
/// # Examples
///
/// mzLib's own reference cases (constructed p-values, not proteomics data), against R's answers:
///
/// ```
/// # mzlib_replay::activate();
/// # let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/stats/");
/// # let read = |name: &str| std::fs::read_to_string(format!("{dir}{name}")).unwrap();
/// use mzlib::stats::{combine_p, PValueTest};
///
/// let inputs = read("combination_inputs.tsv"); // case, p, weight
/// let tests: Vec<PValueTest> = inputs
///     .lines()
///     .skip(1)
///     .map(|line| {
///         let cells: Vec<&str> = line.split('\t').collect();
///         PValueTest::new(cells[0], cells[1].parse::<f64>().unwrap())
///     })
///     .collect();
/// let fisher = combine_p(&tests, "fisher")?;
/// assert_eq!(fisher.columns.integers("studies")?, [Some(3), Some(4), Some(1), Some(12), Some(5), Some(3)]);
///
/// // R's pchisq, column fisher_p of mzLib's reference answers.
/// let r: Vec<f64> = read("combination_results.tsv")
///     .lines()
///     .skip(1)
///     .map(|line| line.split('\t').nth(2).unwrap().parse().unwrap())
///     .collect();
/// let worst = fisher
///     .p_values()?
///     .iter()
///     .zip(&r)
///     .map(|((_, ours), theirs)| (ours.unwrap() - theirs).abs() / theirs)
///     .fold(0.0, f64::max);
/// assert!(worst < 1e-8);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when `method` is not one of
/// [`COMBINATION_METHODS`], `tests` is empty, a feature name is blank or holds a tab or newline, or
/// the bridge predates the verb.
#[doc = include_str!("../docs/reference/stats.combine-p.see-also.md")]
pub fn combine_p(tests: &[PValueTest], method: &str) -> Result<CombinedPValues> {
    let (args, stdin) = combine_request(tests, method)?;
    bridge::require_verb("stats combine-p", MZLIB_1_0_595_BRIDGE)?;
    call(&args, &stdin, Some(SHORT_TIMEOUT))
}

/// [`correlate_with`] with a plain correlation: no covariates.
///
/// # Errors
///
/// As [`correlate_with`].
pub fn correlate(
    responses: impl AsRef<Path>,
    trait_table: impl AsRef<Path>,
) -> Result<Correlations> {
    correlate_with(responses, trait_table, &CorrelateOptions::default())
}

/// Spearman-correlate every feature with a per-sample trait, optionally given covariates.
///
/// For each feature (row) of the responses table, mzLib's `SpearmanCorrelation.Correlate` ranks the
/// feature's values and the trait over the samples where both are known, and tests rho = 0,
/// two-sided. With [`CorrelateOptions::given`], `SpearmanCorrelation.PartialCorrelate` first
/// removes from both sets of ranks what the covariates' ranks explain, as
/// `ppcor::pcor.test(method = "spearman")` does. mzLib holds both to R's `cor.test` and to ppcor 1.1
/// to 1e-8 relative.
///
/// `responses` is [`fit_with`]'s responses table (**0 is an observation**, ranked like any other
/// value). `trait_table` is a TSV with a header row and one row per sample: the sample name
/// (exactly the responses' sample columns, in any order), then ONE value column, whose header
/// names the trait — age, dose, loading. A blank or `NA` cell leaves that sample out. Spearman
/// uses ranks only: it needs no log transform, one extreme value cannot drive it, and it finds any
/// relationship that only ever rises (or only ever falls).
///
/// Needs the bridge from pyMzLib 0.6.0 or later.
#[doc = include_str!("../docs/reference/stats.correlate.md")]
///
/// # Examples
///
/// Age given sex, on limma's reference input:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::stats::{correlate, correlate_with, CorrelateOptions};
///
/// let (responses, age) = ("limma_reference_responses.tsv", "limma_reference_age.tsv");
/// let plain = correlate(responses, age)?;
/// let partial = correlate_with(
///     responses,
///     age,
///     &CorrelateOptions { given: Some("limma_reference_sex.tsv".into()), ..Default::default() },
/// )?;
/// assert!(partial.partial && !plain.partial);
/// assert_eq!(partial.covariate_names, ["sex"]);
/// // The partial test is always asymptotic, so the exact p-values are gone.
/// assert_eq!((plain.method_counts["exact"], partial.method_counts["exact"]), (22, 0));
/// let g2 = |result: &mzlib::stats::Correlations| -> mzlib::Result<f64> {
///     Ok((result.rows()?[1].rho.unwrap() * 100.0).round() / 100.0)
/// };
/// assert_eq!((g2(&plain)?, g2(&partial)?), (0.93, 0.97));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when a path is blank or not UTF-8, or the
/// bridge predates the verb.
#[doc = include_str!("../docs/reference/stats.correlate.see-also.md")]
pub fn correlate_with(
    responses: impl AsRef<Path>,
    trait_table: impl AsRef<Path>,
    options: &CorrelateOptions,
) -> Result<Correlations> {
    let args = correlate_request(responses.as_ref(), trait_table.as_ref(), options)?;
    bridge::require_verb("stats correlate", MZLIB_1_0_595_BRIDGE)?;
    parse(bridge::invoke(&args, None, options.timeout)?)
}

/// [`fit_logistic_with`] with mzLib's iteration limit and tolerance, on one thread.
///
/// # Errors
///
/// As [`fit_logistic_with`].
pub fn fit_logistic<S: AsRef<str>>(
    responses: impl AsRef<Path>,
    design: impl AsRef<Path>,
    coefficients: &[S],
) -> Result<LogisticModelFit> {
    fit_logistic_with(
        responses,
        design,
        coefficients,
        &FitLogisticOptions::default(),
    )
}

/// Fit one logistic regression per feature of a 0/1 table, such as a detection table.
///
/// mzLib's `LogisticRegression.Fit` models logit P(outcome = 1) = X beta for every feature (row)
/// separately, by maximum likelihood, and tests each coefficient with a two-sided Wald z. It is R's
/// `glm(y ~ 0 + X, family = binomial)` per feature, which mzLib reproduces to 1e-8 relative in the
/// coefficients (1e-6 in standard errors and p-values).
///
/// The usual question in proteomics: **is a feature's chance of being measured related to the
/// condition?** Write 1 where a protein was quantified and 0 where it was not; a coefficient far
/// from 0 says its missing values are not missing at random. A feature is reported — never dropped
/// — when it cannot be fitted. `separated` is the one to read carefully: the design splits its 0s
/// from its 1s (or every outcome is the same), so no finite estimate exists. R's `glm` prints a
/// huge coefficient and a warning there; mzLib reports the status and `None`.
///
/// `responses` is [`fit_with`]'s layout with every cell 0, 1, or blank / `NA` / `NaN` for missing
/// (left out of that feature's fit); any other value is refused. `design` is [`fit_with`]'s design;
/// include an intercept column yourself if you want one. `coefficients` are the design columns to
/// report, each as one block of rows; every design column is fitted whatever you name. They travel
/// on stdin.
///
/// Needs the bridge from pyMzLib 0.6.0 or later.
#[doc = include_str!("../docs/reference/stats.fit-logistic.md")]
///
/// # Examples
///
/// mzLib's own reference input (simulated binary features), which R's `glm` output ships beside:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::stats::{fit_logistic_with, FitLogisticOptions};
///
/// let fit = fit_logistic_with(
///     "glm_reference_responses.tsv",
///     "glm_reference_design.tsv",
///     &["Intercept", "x", "grp"],
///     &FitLogisticOptions::default(),
/// )?;
/// assert_eq!((fit.status_counts["fitted"], fit.status_counts["separated"]), (76, 4));
/// assert_eq!((fit.max_iterations, fit.tolerance), (None, None)); // mzLib's defaults were used
///
/// let first = &fit.rows("x")?[0];
/// let four = |x: f64| (x * 1e4).round() / 1e4;
/// assert_eq!((four(first.estimate.unwrap()), four(first.p_value.unwrap())), (0.4592, 0.2724));
/// assert_eq!((first.events, first.observed), (25, 35));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when `coefficients` is empty or holds a blank
/// name or one with a newline, [`FitLogisticOptions::tolerance`] is not finite, or the bridge
/// predates the verb.
#[doc = include_str!("../docs/reference/stats.fit-logistic.see-also.md")]
pub fn fit_logistic_with<S: AsRef<str>>(
    responses: impl AsRef<Path>,
    design: impl AsRef<Path>,
    coefficients: &[S],
    options: &FitLogisticOptions,
) -> Result<LogisticModelFit> {
    let (args, stdin) =
        logistic_request(responses.as_ref(), design.as_ref(), coefficients, options)?;
    bridge::require_verb("stats fit-logistic", MZLIB_1_0_595_BRIDGE)?;
    call(&args, &stdin, options.timeout)
}

/// [`fit_mixed_with`] with mzLib's estimator, on one thread.
///
/// # Errors
///
/// As [`fit_mixed_with`].
pub fn fit_mixed<S: AsRef<str>>(
    responses: impl AsRef<Path>,
    design: impl AsRef<Path>,
    groups: impl AsRef<Path>,
    coefficients: &[S],
) -> Result<MixedModelFit> {
    fit_mixed_with(
        responses,
        design,
        groups,
        coefficients,
        &FitMixedOptions::default(),
    )
}

/// Fit one linear mixed model per feature, with a random intercept for each group of samples.
///
/// mzLib's `MixedModel.Fit` models y = X beta + u_group + e for every feature separately: the fixed
/// effects in your design, plus one baseline per group (a dataset, a donor, an acquisition day)
/// drawn from N(0, tau²). Each coefficient is tested by a two-sided Wald t on nlme's "containment"
/// degrees of freedom. It is `nlme::lme(y ~ ..., random = ~ 1 | group)` per feature, which mzLib
/// reproduces to 1e-5 relative (nlme's own precision), with identical df.
///
/// Use it when samples come in groups that share a baseline and you want each coefficient's effect
/// *within* groups. [`fit_with`] treats every sample as independent, so a day-to-day shift there
/// either inflates the residual variance or, if the days are unbalanced, leaks into the effect.
///
/// `responses` and `design` are [`fit_with`]'s (log-transform intensities first; **0 is an
/// observation**). Do not add a design column per group: the random intercept replaces it.
/// `groups` is a TSV with a header row and one row per sample: the sample name (exactly the
/// responses' samples, in any order), then ONE column with the sample's group label, whose header
/// names the grouping (`day`, `donor`, `dataset`). `coefficients` are the design columns to report;
/// every design column is fitted whatever you name. They travel on stdin.
///
/// Needs the bridge from pyMzLib 0.6.0 or later.
#[doc = include_str!("../docs/reference/stats.fit-mixed.md")]
///
/// # Examples
///
/// mzLib's own reference input (simulated features, samples in seven groups of unequal size; `wp`
/// is constant within each group, `age` varies within them):
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::stats::{fit_mixed_with, FitMixedOptions};
///
/// let fit = fit_mixed_with(
///     "lmm_reference_responses.tsv",
///     "lmm_reference_design.tsv",
///     "lmm_reference_groups.tsv",
///     &["Intercept", "age", "wp"],
///     &FitMixedOptions::default(),
/// )?;
/// assert_eq!((fit.estimator.as_str(), fit.grouping.as_str(), fit.group_count), ("reml", "group", 7));
/// assert_eq!(fit.status_counts["fitted"], 120);
///
/// // wp is a group-level covariate, so it is tested on the groups, not the samples.
/// let mut df: Vec<f64> = fit.rows("wp")?.iter().filter_map(|r| r.df).collect();
/// df.sort_by(f64::total_cmp);
/// df.dedup();
/// assert_eq!(df, [4.0, 5.0]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// Seven groups minus the intercept minus `wp` leaves 5 (4 for the one feature unobserved in a
/// group).
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when `coefficients` is empty or holds a blank
/// name or one with a newline, [`FitMixedOptions::estimator`] is not one of [`MIXED_ESTIMATORS`],
/// or the bridge predates the verb.
#[doc = include_str!("../docs/reference/stats.fit-mixed.see-also.md")]
pub fn fit_mixed_with<S: AsRef<str>>(
    responses: impl AsRef<Path>,
    design: impl AsRef<Path>,
    groups: impl AsRef<Path>,
    coefficients: &[S],
    options: &FitMixedOptions,
) -> Result<MixedModelFit> {
    let (args, stdin) = mixed_request(
        responses.as_ref(),
        design.as_ref(),
        groups.as_ref(),
        coefficients,
        options,
    )?;
    bridge::require_verb("stats fit-mixed", MZLIB_1_0_595_BRIDGE)?;
    call(&args, &stdin, options.timeout)
}

/// [`entrapment_fdp_with`] at mzLib's default thresholds.
///
/// # Errors
///
/// As [`entrapment_fdp_with`].
pub fn entrapment_fdp(
    targets: impl AsRef<Path>,
    entrapments: impl AsRef<Path>,
    ratio: f64,
) -> Result<EntrapmentFdpEstimates> {
    entrapment_fdp_with(
        targets,
        entrapments,
        ratio,
        &EntrapmentFdpOptions::default(),
    )
}

/// Estimate a search's false discovery proportion from its entrapment hits.
///
/// An entrapment experiment adds peptides that cannot be in the sample (from another species, or
/// made from the target peptides) to the search database. Every entrapment the search reports is
/// then a known false discovery, which lets you check the search's FDR claim instead of trusting
/// it. mzLib's `EntrapmentFdp.Sweep` computes the three valid estimators of Wen et al. 2025 (the
/// method FDRBench implements) at each q-value threshold: `lower_bound`, which can only show that a
/// search **fails** to control its FDR; `combined`, an upper bound and so **evidence of control**;
/// and `paired`, a tighter upper bound when each target has exactly one entrapment made from it.
/// The "sample" estimator N_E / (r N_T) is invalid in both directions and is deliberately absent.
///
/// `targets` is a TSV whose header names exactly `score` and `q_value`, in any order: one row per
/// original-target candidate (precursor, peptide or protein group, whichever you are estimating
/// for), decoys removed, each with its best q-value. Higher scores are better. `entrapments` is a
/// TSV whose header names `score`, `q_value`, `partner_score` and `partner_q_value`, and optionally
/// `foreign`: one row per entrapment candidate. The partner columns copy the score and q-value of
/// the one target this entrapment was made from, when the search identified it; leave both blank
/// otherwise. `foreign` is 1 for an entrapment with no partner by construction (another species);
/// blank or absent is mzLib's default, not foreign. Other columns are refused, so select these
/// first.
///
/// `ratio` is r, the entrapment-to-target database size ratio (unitless). It is required: mzLib
/// has no default. `paired` is reported only when it is exactly 1, so for 1:1 pairing pass `1.0`,
/// not a ratio computed from database sizes after exclusions.
///
/// Needs the bridge from pyMzLib 0.6.0 or later.
#[doc = include_str!("../docs/reference/stats.entrapment-fdp.md")]
///
/// # Examples
///
/// mzLib's hand-built search at a ratio of 2: `paired` needs a ratio of exactly 1, so it is absent.
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::stats::{entrapment_fdp_with, EntrapmentFdpOptions};
///
/// let doubled = entrapment_fdp_with(
///     "entrapment_reference_targets.tsv",
///     "entrapment_reference_entrapments.tsv",
///     2.0,
///     &EntrapmentFdpOptions { thresholds: Some(vec![0.01, 0.05]), ..Default::default() },
/// )?;
/// assert!(doubled.thresholds_given);
/// let rows: Vec<(f64, f64, Option<f64>)> = doubled
///     .rows()?
///     .iter()
///     .map(|r| (r.q_value_threshold, (r.combined * 1000.0).round() / 1000.0, r.paired))
///     .collect();
/// assert_eq!(rows, [(0.01, 0.643, None), (0.05, 0.6, None)]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when `ratio` or a threshold is not finite,
/// [`EntrapmentFdpOptions::thresholds`] is an empty list, a path is blank or not UTF-8, or the
/// bridge predates the verb.
#[doc = include_str!("../docs/reference/stats.entrapment-fdp.see-also.md")]
pub fn entrapment_fdp_with(
    targets: impl AsRef<Path>,
    entrapments: impl AsRef<Path>,
    ratio: f64,
    options: &EntrapmentFdpOptions,
) -> Result<EntrapmentFdpEstimates> {
    let args = entrapment_request(targets.as_ref(), entrapments.as_ref(), ratio, options)?;
    bridge::require_verb("stats entrapment-fdp", MZLIB_1_0_595_BRIDGE)?;
    parse(bridge::invoke(&args, None, options.timeout)?)
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
    fn a_trend_fit_without_a_spline_basis_sends_none_so_mzlib_chooses() {
        // mzLib then follows limma's fitFDist rule for the number of features (mzLib #1398); the
        // crate makes up no basis count of its own.
        let (args, _) = fit_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            &["a"],
            &FitOptions {
                trend: true,
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
                "1",
                "--trend"
            ]
        );
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

    // ---- stats fit: contrasts, confidence, estimator ------------------------------------------

    const CONTRASTS: &str = include_str!("../tests/fixtures/stats_fit_malat_contrasts.json");
    const MARGINAL: &str = include_str!("../tests/fixtures/stats_fit_malat_marginal.json");

    fn data(recording: &str) -> Value {
        serde_json::from_str(recording).unwrap()
    }

    #[test]
    fn a_fit_sends_contrasts_confidence_and_estimator_only_when_given() {
        let (args, _) = fit_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            &["a"],
            &FitOptions::default(),
        )
        .unwrap();
        for option in ["--contrasts", "--confidence", "--estimator"] {
            assert!(!args.iter().any(|a| a == option), "{option} was sent");
        }
        let (args, stdin) = fit_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            &["a"],
            &FitOptions {
                contrasts: Some("c.tsv".into()),
                confidence: Some(0.9),
                estimator: Some("marginal_likelihood".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            &args[8..],
            [
                "--contrasts",
                "c.tsv",
                "--confidence",
                "0.9",
                "--estimator",
                "marginal_likelihood"
            ]
        );
        assert_eq!(stdin, "a\n");
    }

    #[test]
    fn contrasts_alone_send_an_empty_stdin_and_no_coefficient_at_all_is_refused() {
        let none: [&str; 0] = [];
        let (_, stdin) = fit_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            &none,
            &FitOptions {
                contrasts: Some("c.tsv".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(stdin, "");
        let error = fit_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            &none,
            &FitOptions::default(),
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(ref m) if m.contains("contrasts file")));
    }

    #[test]
    fn bad_contrast_options_are_refused_before_anything_is_spawned() {
        for options in [
            FitOptions {
                estimator: Some("limma".into()),
                ..Default::default()
            },
            FitOptions {
                contrasts: Some("c.tsv".into()),
                confidence: Some(f64::NAN),
                ..Default::default()
            },
            FitOptions {
                contrasts: Some(" ".into()),
                ..Default::default()
            },
        ] {
            let error =
                fit_request(Path::new("r.tsv"), Path::new("d.tsv"), &["a"], &options).unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{options:?}");
        }
    }

    #[test]
    fn a_fit_parses_the_recorded_contrasts() {
        let fit: ModeratedFit = serde_json::from_str(CONTRASTS).unwrap();
        assert!(fit.tested_coefficients.is_empty());
        assert_eq!((fit.row_count, fit.columns.rows()), (0, 0));
        assert!(fit
            .contrasts_file
            .as_deref()
            .is_some_and(|f| f.ends_with("contrasts_1.tsv")));
        let contrasts = fit.contrasts.as_ref().unwrap();
        assert_eq!(contrasts.tested_contrasts, ["c1", "c2", "c3"]);
        assert_eq!(contrasts.confidence_level, 0.95);
        assert_eq!(contrasts.weights["c3"], [0.0, -1.0, 1.0]);
        assert_eq!(contrasts.row_count, 1476);
        let rows = contrasts.rows("c3").unwrap();
        assert_eq!(rows.len(), 492);
        let fitted = rows.iter().find(|r| r.status == "fitted").unwrap();
        let (low, estimate, high) = (
            fitted.confidence_low.unwrap(),
            fitted.estimate.unwrap(),
            fitted.confidence_high.unwrap(),
        );
        assert!(low < estimate && estimate < high);
        let unfitted = rows.iter().find(|r| r.status != "fitted").unwrap();
        assert_eq!(
            (
                unfitted.estimate,
                unfitted.confidence_low,
                unfitted.bh_adjusted
            ),
            (None, None, None)
        );
        let error = contrasts.rows("c4").unwrap_err();
        assert!(matches!(error, MzLibError::Usage(ref m) if m.contains("c1, c2, c3")));
    }

    #[test]
    fn a_fit_without_contrasts_has_none_and_reports_its_estimator() {
        let fit: ModeratedFit = serde_json::from_str(MALAT).unwrap();
        assert!(fit.contrasts.is_none() && fit.contrasts_file.is_none());
        assert_eq!(fit.prior.estimator.as_deref(), Some("moments_legacy"));
        let marginal: ModeratedFit = serde_json::from_str(MARGINAL).unwrap();
        assert_eq!(
            marginal.prior.estimator.as_deref(),
            Some("marginal_likelihood")
        );
        assert!(marginal.prior.df.unwrap() < fit.prior.df.unwrap());
    }

    #[test]
    fn an_answer_that_ignored_contrasts_or_the_estimator_is_refused() {
        // A bridge before pyMzLib 0.6.0 answers without the keys the options add.
        let mut old = data(MALAT);
        old.as_object_mut().unwrap().remove("contrasts");
        old["prior"].as_object_mut().unwrap().remove("estimator");
        let contrasts = FitOptions {
            contrasts: Some("c.tsv".into()),
            ..Default::default()
        };
        let estimator = FitOptions {
            estimator: Some("marginal_likelihood".into()),
            ..Default::default()
        };
        for (options, option) in [(contrasts, "contrasts"), (estimator, "estimator")] {
            let error = check_fit_options_reached(&old, &options).unwrap_err();
            assert!(
                matches!(&error, MzLibError::Usage(m) if m.contains(&format!("--{option}")) && m.contains("0.6.0")),
                "{error}"
            );
        }
        assert!(check_fit_options_reached(&old, &FitOptions::default()).is_ok());
        assert!(check_fit_options_reached(
            &data(CONTRASTS),
            &FitOptions {
                contrasts: Some("c.tsv".into()),
                ..Default::default()
            }
        )
        .is_ok());
    }

    // ---- stats combine-p ----------------------------------------------------------------------

    const FISHER: &str = include_str!("../tests/fixtures/stats_combine_p_fisher.json");
    const SIGNED: &str = include_str!("../tests/fixtures/stats_combine_p_signed.json");

    #[test]
    fn combine_p_writes_one_line_per_test_with_blank_missing_fields() {
        let (args, stdin) = combine_request(
            &[
                PValueTest::new("a", 0.01).weight(2.0).effect(-0.5),
                PValueTest::new("a", None).weight(3.0).effect(0.4),
                PValueTest::new("b", f64::NAN).effect(1.0),
            ],
            "stouffer",
        )
        .unwrap();
        assert_eq!(args, ["stats", "combine-p", "--method", "stouffer"]);
        assert_eq!(stdin, "a\t0.01\t2.0\t-0.5\na\t\t3.0\t0.4\nb\t\t\t1.0\n");
    }

    #[test]
    fn combine_p_sends_only_the_fields_the_tests_have() {
        let (_, plain) = combine_request(
            &[PValueTest::new("a", 0.2), PValueTest::new("b", 0.3)],
            "fisher",
        )
        .unwrap();
        assert_eq!(plain, "a\t0.2\nb\t0.3\n");
        let (_, weighted) =
            combine_request(&[PValueTest::new("a", 0.2).weight(4.0)], "stouffer").unwrap();
        assert_eq!(weighted, "a\t0.2\t4.0\n");
    }

    #[test]
    fn combine_p_refuses_bad_input_before_anything_is_spawned() {
        assert!(combine_request(&[PValueTest::new("a", 0.1)], "tippett").is_err());
        assert!(combine_request(&[], "fisher").is_err());
        for feature in ["", " ", "a\tb", "a\nb"] {
            let error = combine_request(&[PValueTest::new(feature, 0.1)], "fisher").unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{feature:?}");
        }
    }

    #[test]
    fn combine_p_parses_the_recordings() {
        let fisher: CombinedPValues = serde_json::from_str(FISHER).unwrap();
        assert_eq!(fisher.method, "fisher");
        assert!(!fisher.weighted && !fisher.signed);
        assert_eq!((fisher.line_count, fisher.feature_count), (28, 6));
        let p = fisher.p_values().unwrap();
        assert_eq!(p.len(), 6);
        assert_eq!(p[0].0, "few");
        let signed: CombinedPValues = serde_json::from_str(SIGNED).unwrap();
        assert!(signed.weighted && signed.signed);
        let flips = signed
            .p_values()
            .unwrap()
            .into_iter()
            .find(|(f, _)| f == "flips_sign")
            .unwrap();
        assert!((flips.1.unwrap() - 0.43).abs() < 0.01);
    }

    // ---- stats correlate ----------------------------------------------------------------------

    const CORRELATE: &str = include_str!("../tests/fixtures/stats_correlate_malat.json");
    const PARTIAL: &str =
        include_str!("../tests/fixtures/stats_correlate_limma_age_given_sex.json");

    #[test]
    fn correlate_sends_the_files_and_given_only_when_set() {
        let args = correlate_request(
            Path::new("r.tsv"),
            Path::new("t.tsv"),
            &CorrelateOptions::default(),
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "stats",
                "correlate",
                "--responses",
                "r.tsv",
                "--trait",
                "t.tsv"
            ]
        );
        let args = correlate_request(
            Path::new("r.tsv"),
            Path::new("t.tsv"),
            &CorrelateOptions {
                given: Some("g.tsv".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(&args[6..], ["--given", "g.tsv"]);
        assert!(correlate_request(
            Path::new(""),
            Path::new("t.tsv"),
            &CorrelateOptions::default()
        )
        .is_err());
    }

    #[test]
    fn correlate_parses_the_plain_and_partial_recordings() {
        let plain: Correlations = serde_json::from_str(CORRELATE).unwrap();
        assert_eq!(plain.trait_name, "malat_ng");
        assert!(!plain.partial && plain.given_file.is_none());
        assert_eq!(plain.given_incomplete_count, None);
        assert_eq!(plain.method_counts["not_estimable"], 180);
        let rows = plain.rows().unwrap();
        assert_eq!(rows.len(), 492);
        assert!(SPEARMAN_METHODS.contains(&rows[0].method.as_str()));
        let missing = rows.iter().find(|r| r.method == "not_estimable").unwrap();
        assert_eq!((missing.rho, missing.p_value), (None, None));

        let partial: Correlations = serde_json::from_str(PARTIAL).unwrap();
        assert!(partial.partial && partial.given_file.is_some());
        assert_eq!(partial.covariate_names, ["sex"]);
        assert_eq!(partial.given_incomplete_count, Some(0));
    }

    // ---- stats fit-logistic, fit-mixed --------------------------------------------------------

    const GLM: &str = include_str!("../tests/fixtures/stats_fit_logistic_glm.json");
    const LMM: &str = include_str!("../tests/fixtures/stats_fit_mixed_lmm_reml.json");
    const LMM_ML: &str = include_str!("../tests/fixtures/stats_fit_mixed_lmm_ml.json");

    #[test]
    fn fit_logistic_sends_only_what_was_given() {
        let (args, stdin) = logistic_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            &["x"],
            &FitLogisticOptions::default(),
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "stats",
                "fit-logistic",
                "--responses",
                "r.tsv",
                "--design",
                "d.tsv",
                "--threads",
                "1"
            ]
        );
        assert_eq!(stdin, "x\n");
        let (args, _) = logistic_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            &["x"],
            &FitLogisticOptions {
                max_iterations: Some(50),
                tolerance: Some(1e-10),
                threads: -1,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            &args[6..],
            [
                "--threads",
                "-1",
                "--max-iterations",
                "50",
                "--tolerance",
                "1e-10"
            ]
        );
    }

    #[test]
    fn fit_logistic_parses_the_glm_recording() {
        let fit: LogisticModelFit = serde_json::from_str(GLM).unwrap();
        assert_eq!(fit.tested_coefficients, ["Intercept", "x", "grp"]);
        assert_eq!((fit.max_iterations, fit.tolerance), (None, None));
        assert_eq!(fit.status_counts["separated"], 4);
        let rows = fit.rows("x").unwrap();
        assert_eq!(rows.len() as u64, fit.feature_count);
        assert!(rows
            .iter()
            .all(|r| LOGISTIC_STATUSES.contains(&r.status.as_str())));
        let separated = rows.iter().find(|r| r.status == "separated").unwrap();
        assert_eq!((separated.estimate, separated.z), (None, None));
        assert!(separated.observed > 0);
        assert!(matches!(fit.rows("y"), Err(MzLibError::Usage(_))));
    }

    #[test]
    fn fit_mixed_sends_the_estimator_only_when_given_and_refuses_an_unknown_one() {
        let (args, _) = mixed_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            Path::new("g.tsv"),
            &["age"],
            &FitMixedOptions::default(),
        )
        .unwrap();
        assert!(!args.iter().any(|a| a == "--estimator"));
        assert_eq!(&args[6..8], ["--groups", "g.tsv"]);
        let (args, _) = mixed_request(
            Path::new("r.tsv"),
            Path::new("d.tsv"),
            Path::new("g.tsv"),
            &["age"],
            &FitMixedOptions {
                estimator: Some("ml".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(&args[args.len() - 2..], ["--estimator", "ml"]);
        let none: [&str; 0] = [];
        for (coefficients, estimator) in [(&none[..], None), (&["age"][..], Some("REML"))] {
            let error = mixed_request(
                Path::new("r.tsv"),
                Path::new("d.tsv"),
                Path::new("g.tsv"),
                coefficients,
                &FitMixedOptions {
                    estimator: estimator.map(str::to_owned),
                    ..Default::default()
                },
            )
            .unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)));
        }
    }

    #[test]
    fn fit_mixed_parses_the_reference_recordings_and_a_null_df_has_no_p_value() {
        let reml: MixedModelFit = serde_json::from_str(LMM).unwrap();
        let ml: MixedModelFit = serde_json::from_str(LMM_ML).unwrap();
        assert_eq!(
            (reml.estimator.as_str(), ml.estimator.as_str()),
            ("reml", "ml")
        );
        assert_eq!(reml.group_count, 7);
        let rows = reml.rows("wp").unwrap();
        assert!(rows
            .iter()
            .all(|r| MIXED_STATUSES.contains(&r.status.as_str())));
        for row in &rows {
            if row.df.is_none() {
                assert!(row.p_value.is_none());
            }
        }
        assert!(matches!(reml.rows("sex"), Err(MzLibError::Usage(_))));
    }

    // ---- stats entrapment-fdp -----------------------------------------------------------------

    const ENTRAPMENT: &str = include_str!("../tests/fixtures/stats_entrapment_fdp_reference.json");
    const ENTRAPMENT_R2: &str =
        include_str!("../tests/fixtures/stats_entrapment_fdp_reference_r2.json");

    #[test]
    fn entrapment_fdp_sends_the_ratio_and_thresholds_only_when_given() {
        let args = entrapment_request(
            Path::new("t.tsv"),
            Path::new("e.tsv"),
            1.0,
            &EntrapmentFdpOptions::default(),
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "stats",
                "entrapment-fdp",
                "--targets",
                "t.tsv",
                "--entrapments",
                "e.tsv",
                "--ratio",
                "1.0"
            ]
        );
        let args = entrapment_request(
            Path::new("t.tsv"),
            Path::new("e.tsv"),
            2.0,
            &EntrapmentFdpOptions {
                thresholds: Some(vec![0.01, 0.05]),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(&args[7..], ["2.0", "--thresholds", "0.01,0.05"]);
        for (ratio, thresholds) in [
            (f64::NAN, None),
            (1.0, Some(vec![])),
            (1.0, Some(vec![f64::INFINITY])),
        ] {
            assert!(entrapment_request(
                Path::new("t.tsv"),
                Path::new("e.tsv"),
                ratio,
                &EntrapmentFdpOptions {
                    thresholds,
                    ..Default::default()
                },
            )
            .is_err());
        }
    }

    #[test]
    fn entrapment_fdp_parses_the_recordings_and_paired_is_none_off_a_ratio_of_one() {
        let fdp: EntrapmentFdpEstimates = serde_json::from_str(ENTRAPMENT).unwrap();
        assert!(!fdp.thresholds_given);
        assert_eq!((fdp.threshold_count, fdp.row_count), (100, 100));
        let at = fdp.at(0.01).unwrap().unwrap();
        assert_eq!((at.target_count, at.entrapment_count), (4, 3));
        assert!(at.paired.is_some());
        assert!(matches!(fdp.at(0.0123), Err(MzLibError::Usage(_))));

        let doubled: EntrapmentFdpEstimates = serde_json::from_str(ENTRAPMENT_R2).unwrap();
        assert!(doubled.thresholds_given);
        assert_eq!(doubled.ratio, 2.0);
        assert!(doubled.rows().unwrap().iter().all(|r| r.paired.is_none()));
    }
}
