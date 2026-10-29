//! PTM site occupancy and PTM-QTL statistics, computed by mzLib's `Quantification.PtmQtl`.
//!
//! | You want to know | Call | mzLib computes it with |
//! |---|---|---|
//! | What fraction of each site carries its modification, run by run? | [`occupancy_with`] | `SiteOccupancyCalculator.Calculate` |
//! | The same for a sample measured as several fractions? | [`occupancy_with`] with [`OccupancyOptions::samples`] | `SiteOccupancyCalculator.Calculate` with a run-to-sample map |
//! | Which two modifications sit on one molecule? | [`pairs_with`] with [`PairType::Physical`] | `PtmPairEngine.Physical` |
//! | Which two sites rise and fall together across runs? | [`pairs_with`] with [`PairType::Covarying`] | `PtmPairEngine.CoVarying` (Spearman, BH) |
//! | The same, counting runs where a site was seen in only one form? | [`pairs_with`] with [`PairsOptions::include_bounds`] | `PtmPairEngine.CoVarying(includeBounds)`, `PairTestStatus` |
//! | The same, given each run's depth? | [`pairs_with`] with [`PairsOptions::run_covariate`] | `PtmPairEngine.CoVarying(runCovariate)` |
//! | Does a site's occupancy change with a per-run trait? | [`trait_effects_with`] | `SiteTraitEffects.Fit` (a mixed model per site) |
//! | The same, with the effect in log2 odds? | [`trait_effects_with`] with [`TraitEffectsOptions::scale`] | `SiteTraitOptions.Scale`, `OccupancyScale` |
//!
//! A modification site is rarely all-or-nothing. In one run 15% of a serine is phosphorylated, in
//! the next 40%. That fraction, the site's **occupancy**, is what changes with age, dose or
//! genotype, and what two sites can share. Measuring it from label-free data is easy to get wrong:
//! a modified form that was identified but not quantified is not a zero, a site seen only in its
//! modified form is not truly 100%, and two sites read off one peptide are correlated by
//! construction. mzLib's `Quantification.PtmQtl` engine (mzLib #1369, #1430, #1446) keeps those
//! cases apart (mzLib #1464 and #1468 added the peptidoform counts, the reason a pair was not
//! tested, and the log2 scale). This crate sends your table to the bridge and types what comes
//! back; it does no arithmetic of its own.
//!
//! The examples on this page read mzLib's own hand-computed test cases, so each number below is
//! one mzLib's tests assert.
//!
//! ## The input: one row per peptidoform, per run, per protein
//!
//! Every function reads a tab-separated table with a header row: mzLib's own
//! `PeptidoformObservation`, as columns.
//!
//! | column | what it holds |
//! |---|---|
//! | `run` | the run (raw file) the peptidoform was identified in |
//! | `full_sequence` | the peptidoform in MetaMorpheus notation: `PEPS[Common Biological:Phosphorylation on S]K` |
//! | `protein_accession` | one protein it maps to |
//! | `start_residue`, `end_residue` | the peptide's first and last residue in that protein, counting from 1 |
//! | `intensity` | its MS/MS-quantified intensity in that run, in instrument units; **blank** when it was identified but not quantified |
//! | `previous_residue` | optional: the protein residue before the peptide, one letter |
//!
//! Four rules decide whether the numbers mean anything:
//!
//! - **One row per (run, peptidoform, protein), never one per PSM.** Two rows for one peptidoform
//!   in one run would count its intensity twice; mzLib refuses them.
//! - **MS/MS identifications only.** A match-between-runs transfer is not an identification.
//! - **Blank, not 0, for "identified but not quantified".** FlashLFQ writes 0 for a peptide it did
//!   not quantify. Blank those cells.
//! - **A shared peptide is one row per protein.** Each protein then gets its full intensity.
//!   MetaMorpheus instead splits a shared peptide across protein groups, so its
//!   `IntensityOccupancy_` columns can differ.
//!
//! **mzLib does not build this table for you yet.** It has no reader that turns a MetaMorpheus or
//! FlashLFQ result into these rows, because that step makes choices — which PSMs count, how a
//! shared peptide is placed — and this crate will not make them on mzLib's behalf. Assemble it from
//! the PSMs ([`crate::readers::read_records`] on `AllPSMs.psmtsv`: runs, proteins and positions) and
//! the FlashLFQ peptide table ([`crate::readers::read_quantified_peptides`]: intensities), keeping
//! one row per (run, full sequence, protein) and a blank intensity where FlashLFQ has no MS/MS
//! value. Check that the run names of the two files match: FlashLFQ labels a calibrated file
//! `<run>-calib`.
//!
//! ## Occupancy, run by run
//!
//! For a site (protein, position, modification) in a run, occupancy is the summed intensity of the
//! quantified peptidoforms that carry the modification there, over the summed intensity of every
//! quantified peptidoform that covers the position — the estimator MetaMorpheus writes as
//! `IntensityOccupancy_`. The cases file holds several of mzLib's cases, each on its own protein,
//! named for what it tests:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::ptmqtl::occupancy;
//!
//! let occ = occupancy("cases.tsv")?;
//! assert_eq!((occ.observation_count, occ.site_count), (22, 9));
//! let rows = occ.rows()?;
//! let ratio = rows.iter().find(|r| r.protein_accession == "RATIO").unwrap();
//! assert_eq!(ratio.site_key, "RATIO:S13:Common Biological:Phosphorylation on S");
//! assert_eq!((ratio.position, ratio.fraction, ratio.unmodified_quantified), (13, Some(0.15), true));
//!
//! // Each cell has a state, and only a quantified cell has a fraction.
//! let states: Vec<(&str, &str, Option<f64>)> = rows
//!     .iter()
//!     .filter(|r| r.protein_accession == "STATES")
//!     .map(|r| (r.run.as_str(), r.state.as_str(), r.fraction))
//!     .collect();
//! assert_eq!(
//!     states,
//!     [
//!         ("countOnly", "count_only", None),
//!         ("floor", "floor", None),
//!         ("notDetected", "not_detected", None),
//!         ("quant", "quantified", Some(0.1)),
//!     ]
//! );
//!
//! // A 1 can be a ceiling: the site was seen only in its modified form.
//! let ceiling = rows.iter().find(|r| r.protein_accession == "CEILING").unwrap();
//! assert_eq!((ceiling.fraction, ceiling.unmodified_quantified), (Some(1.0), false));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! `RATIO` has 30 units of the phosphopeptide, 70 of the same peptide unmodified, and 100 of a
//! missed-cleavage peptide that also covers the serine: 30 / 200 = 0.15. Every quantified form that
//! covers the position counts in the denominator, whatever its own sequence. The states
//! ([`OCCUPANCY_STATES`]):
//!
//! | state | what was seen | read it as |
//! |---|---|---|
//! | `quantified` | the modified form was quantified | a measured fraction in (0, 1] |
//! | `floor` | the modified form was identified but not quantified, while forms covering the site were | low, **not zero** |
//! | `count_only` | the modified form was identified, but nothing covering the site was quantified | present, amount unknown |
//! | `not_detected` | peptidoforms covering the site were identified in the run, and none carried the modification | a bound below what the run could detect, **not zero** |
//!
//! Neither a floor nor an undetected cell has a fraction, so neither can be averaged in as 0 by
//! accident. Since mzLib #1468 a `not_detected` cell is mzLib's covered-zero: the forms covering
//! the site need only be identified, not quantified. Each row also counts identifications, not
//! intensity: [`OccupancyRow::covering_peptidoforms`] distinct peptidoforms covered the position in
//! that run, and [`OccupancyRow::modified_peptidoforms`] of them carried the modification there —
//! the k of n of a cell with no fraction. A `count_only` cell with k = n was seen only modified;
//! with k < n both forms were identified. Two modifications at one position share a denominator: on
//! `SHARED_POSITION`, 20 units carry phosphorylation, 30 HexNAc and 50 neither, so they are two
//! sites, at 0.2 and 0.3. **Positions count from 1, in the protein**: a protein N-terminal
//! modification sits at the residue it is on, residue 2 when the initiator methionine was removed.
//! So **do not join this table with [`crate::readers::read_occupancy`] on position**:
//! MetaMorpheus's occupancy file puts the N-terminus at position 0.
//!
//! **What mzLib leaves out.** mzLib always skips the `Common Variable` and `Common Fixed`
//! categories, peptide-terminal modifications and every C-terminal modification (the table does not
//! carry the protein's length). An N-terminal modification is a site only on a peptide that starts
//! at residue 1, or at 2 after the initiator Met is removed. A modification name that mzLib's
//! registry does not know is treated as protein N-terminal. To keep only some modifications, set
//! [`OccupancyOptions::modifications`] (exact names) or [`OccupancyOptions::categories`].
//!
//! ## A sample measured as several fractions
//!
//! Fractions of one sample are summed before the ratio, as FlashLFQ sums them. Give the
//! run-to-sample map as a file:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::ptmqtl::{occupancy_with, OccupancyOptions};
//!
//! let occ = occupancy_with(
//!     "fractions.tsv",
//!     &OccupancyOptions { samples: Some("fractions_samples.tsv".into()), ..Default::default() },
//! )?;
//! assert_eq!((occ.run_count, occ.sample_count, occ.run_names.clone()), (3, Some(2), vec!["A".to_owned(), "B".to_owned()]));
//! let cells: Vec<(String, String, String, Option<f64>)> = occ
//!     .rows()?
//!     .into_iter()
//!     .map(|r| (r.protein_accession, r.run, r.state, r.fraction))
//!     .collect();
//! assert_eq!(cells[1], ("SUMMED".to_owned(), "A".to_owned(), "quantified".to_owned(), Some(0.2)));
//! assert_eq!(cells[0].2, "floor"); // FLOOR_ACROSS
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! `SUMMED` reads 0.3 in fraction 1 and 0.1 in fraction 2. Sample A is (30 + 10) / (100 + 100) =
//! 0.2, not the mean of 0.3 and 0.1. **Map only fractions together.** Technical replicates are
//! separate measurements of one sample; give each its own name, and let [`trait_effects_with`] treat
//! them as replicates of one biological sample.
//!
//! ## Pairs of modifications
//!
//! A peptidoform carrying two modifications is direct evidence that both sit on one molecule:
//! [`PairType::Physical`]. `statistic` is then the median **co-occupancy** over the runs where the
//! doubly modified form was quantified, and there is no test:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::ptmqtl::{pairs, PairType};
//!
//! let found = pairs("cases.tsv", PairType::Physical)?;
//! let pair = &found.rows()?[0];
//! assert_eq!(pair.site_a, "TWO_ON_ONE:K14:Common Biological:Acetylation on K");
//! assert_eq!(pair.site_b, "TWO_ON_ONE:S13:Common Biological:Phosphorylation on S");
//! assert_eq!((pair.statistic, pair.n, pair.statistic_n, pair.p_value), (Some(0.25), 2, 1, None));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Identified together in two runs, quantified in one, where 25 of the 100 intensity units covering
//! both positions carried both marks.
//!
//! Two sites that rise and fall together across runs are co-varying ([`PairType::Covarying`]).
//! mzLib correlates their occupancies with Spearman's rho on the runs where both were quantified,
//! and adjusts the p-values with Benjamini-Hochberg within each family. **The trap: a batch
//! correlates everything.** These twelve runs come in two batches, and every site reads 0.2 higher
//! in the second. [`PairsOptions::remove_run_level`] correlates the two sites given each run's
//! **level**, the median occupancy of the other sites in that run:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::ptmqtl::{pairs, pairs_with, PairType, PairsOptions, PtmPairs};
//!
//! let plain = pairs("batch.tsv", PairType::Covarying)?;
//! // mzLib's default guard, reported because none was set. Ceilings are left out by mzLib's own
//! // rule, which the crate does not restate: exclude_ceiling is None unless you chose.
//! assert_eq!((plain.covarying_run_count, plain.min_quantified_fraction, plain.exclude_ceiling), (Some(12), Some(0.7), None));
//!
//! let given = pairs_with(
//!     "batch.tsv",
//!     PairType::Covarying,
//!     &PairsOptions { remove_run_level: true, ..Default::default() },
//! )?;
//! assert_eq!(given.run_level_group.as_deref(), Some("none"));
//! assert_eq!((given.min_runs_for_run_level, given.min_sites_for_run_level), (Some(10), Some(3)));
//!
//! let rho = |result: &PtmPairs| -> mzlib::Result<f64> {
//!     let pair = result
//!         .rows()?
//!         .into_iter()
//!         .find(|r| r.site_a.starts_with("PAIR_S") && r.site_b.starts_with("PAIR_T"))
//!         .unwrap();
//!     Ok((pair.statistic.unwrap() * 100.0).round() / 100.0)
//! };
//! assert_eq!((rho(&plain)?, rho(&given)?), (0.76, 0.03));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! The batch alone gives `PAIR_S` and `PAIR_T` a rho of 0.76; given each run's level, the pair is
//! unrelated. A site enters co-varying pairs only if it is quantified in at least
//! [`PtmPairs::min_quantified_fraction`] of the runs, and ceilings are left out by default. **Remove
//! the run level for a batch or loading shift only**: if the level follows your design — one
//! condition reads higher everywhere — removing it removes the effect you are studying. Sites of
//! different chemistry can have different batch shifts: [`PairsOptions::run_level_group`] gives each
//! category or modification its own level, and [`PairsOptions::run_level_groups`] any grouping, from
//! a two-column file (`modification`, `group`).
//!
//! What to read in the table: **`overlapping` pairs have no `q`** — when one peptide covers both
//! sites, their occupancies share intensities and are correlated by construction, so mzLib reports
//! them but leaves them out of the adjustment. `fdr_family` is the type (`P` physical, `A`
//! co-varying) and whether the two sites are on one protein (`intra`) or two (`inter`); `q` is
//! adjusted within it, and is not a target-decoy q-value. [`PairsOptions::keep_ceiling`] keeps runs
//! where a site was seen only modified: two unrelated sites that both read 1 in low-load runs then
//! correlate through detection alone, which is why mzLib leaves them out by default.
//! [`PtmPair::test_status`] says why a co-varying pair has no statistic: `not_estimable` (too few
//! complete runs, or a site constant over them), `too_few_shared_runs` (with the run level or a
//! run covariate removed, fewer than [`PtmPairs::min_runs_for_run_level`] runs) or
//! `both_sites_mostly_bounds` (next section), one of [`PAIR_TEST_STATUSES`]. An untested pair has
//! no `statistic`, `p_value` or `q`, and is not in its family's adjustment.
//!
//! ### Counting one-form runs as bounds
//!
//! A site seen in only one form in a run has no fraction, but it is not missing either. A
//! `not_detected` cell says the occupancy was below what the run could see; a ceiling says it was
//! at the top. By default mzLib leaves such runs out of the correlation. A site that is often below
//! detection then never reaches the 70% guard, and its pairs are never tested.
//!
//! [`PairsOptions::include_bounds`] lets those runs in as **tied ranks**, never as invented values
//! (mzLib #1468). A `not_detected` or `floor` cell ties below every fraction of its site. A ceiling,
//! or a `count_only` cell whose covering peptidoforms all carry the mark (k of n with k = n), ties
//! above. A `count_only` cell where both forms were identified has no rank and stays out. The
//! example is mzLib's `BoundsEnterTypeAAsTiedRanks`: two sites with a value in six of ten runs, and
//! a bound in the other four.
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::ptmqtl::{pairs, pairs_with, PairType, PairsOptions};
//!
//! let plain = pairs("bounds.tsv", PairType::Covarying)?;
//! assert_eq!((plain.row_count, plain.include_bounds), (0, Some(false)));
//!
//! let bounded = pairs_with(
//!     "bounds.tsv",
//!     PairType::Covarying,
//!     &PairsOptions { include_bounds: true, ..Default::default() },
//! )?;
//! let pair = &bounded.rows()?[0];
//! assert_eq!((pair.n, pair.test_status.as_deref(), pair.bounds_included), (9, Some("tested"), true));
//! assert_eq!(format!("{:.3}", pair.statistic.unwrap()), "0.004");
//! assert_eq!((bounded.min_untied_share, bounded.exclude_ceiling), (Some(0.5), None));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Six of ten is below the guard, so the plain call tests nothing. With bounds, `BOUND_S` has a
//! value or a bound in all ten runs and `BOUND_T` in nine, and the pair is tested on the nine they
//! share. The tenth is `BOUND_T`'s `count_only` cell with both forms identified (1 of 2): it has no
//! rank. The two sites turn out unrelated (rho 0.004).
//!
//! **Too many ties make the test unreliable.** When both sites are mostly one bound, Spearman's t
//! test calls far too many pairs: in mzLib's simulation, 2 to 5 times the nominal rate at
//! p < 0.001 with no relation at all. So a pair is tested only if one of its sites has at least
//! [`PtmPairs::min_untied_share`] (half) of the shared runs outside its largest tie. Otherwise its
//! [`PtmPair::test_status`] is `both_sites_mostly_bounds`, and it has no p-value.
//! [`PtmPairs::exclude_ceiling`] stays `None`: with bounds, mzLib treats a ceiling as a high bound,
//! and [`PairsOptions::keep_ceiling`] gives the same answer.
//!
//! ### Given each run's depth
//!
//! Shallow runs miss low-abundance forms everywhere at once. With bounds in, that makes
//! covered-zeros co-occur across unrelated sites, and any site that rises with depth correlates with
//! every other. [`PairsOptions::run_covariate`] takes a per-run value, such as the run's target
//! PSMs, and correlates the two sites **given** it: a partial Spearman correlation, as for the run
//! level. The file is tab-separated, with a header row whose second field names the covariate,
//! then one row per run (`depth_psms.tsv`, tabs shown as spaces):
//!
//! ```text
//! run   target_psms
//! r00   9000
//! r01   12000
//! ```
//!
//! The example is mzLib's `ARunCovariateIsPartialledOut`: twelve runs in two batches, and two sites
//! that also rise a little with each run's depth.
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::ptmqtl::{pairs, pairs_with, PairType, PairsOptions, PtmPairs};
//!
//! let depth_pair = |result: &PtmPairs| -> mzlib::Result<String> {
//!     let r = result
//!         .rows()?
//!         .into_iter()
//!         .find(|r| r.site_a.starts_with("DEPTH_S") && r.site_b.starts_with("DEPTH_T"))
//!         .unwrap();
//!     Ok(format!("rho {:.2}  p {:.3}  n {}", r.statistic.unwrap(), r.p_value.unwrap(), r.n))
//! };
//! assert_eq!(depth_pair(&pairs("depth.tsv", PairType::Covarying)?)?, "rho 0.77  p 0.003  n 12");
//!
//! let given = pairs_with(
//!     "depth.tsv",
//!     PairType::Covarying,
//!     &PairsOptions { run_covariate: Some("depth_psms.tsv".into()), ..Default::default() },
//! )?;
//! assert_eq!(depth_pair(&given)?, "rho 0.76  p 0.006  n 12");
//! assert_eq!(given.run_covariate_name.as_deref(), Some("target_psms"));
//! assert_eq!(given.run_covariate_runs_missing.as_deref(), Some(&[][..]));
//! assert_eq!(given.min_runs_for_run_level, Some(10));
//!
//! let both = pairs_with(
//!     "depth.tsv",
//!     PairType::Covarying,
//!     &PairsOptions {
//!         run_covariate: Some("depth_psms.tsv".into()),
//!         remove_run_level: true,
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!(depth_pair(&both)?, "rho 0.08  p 0.832  n 12");
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Depth alone barely moves rho: in these runs most of the shared rise is the batch. Given the
//! batch's run level as well, the pair is unrelated. A covariate removes only what it explains.
//!
//! - **Any increasing transform gives the same answer.** The covariate is ranked, so target PSMs
//!   and log2 target PSMs are the same input.
//! - **Leave a run out for "unknown".** It is then left out of every pair, and listed in
//!   [`PtmPairs::run_covariate_runs_missing`]. A blank or `NA` value is not "unknown": mzLib refuses
//!   it, and the call returns [`crate::MzLibError::Usage`].
//! - **The same floor of shared runs applies.** A pair partial on the covariate is tested only on
//!   at least [`PtmPairs::min_runs_for_run_level`] shared runs; below that it is
//!   `too_few_shared_runs`.
//!
//! ## Does occupancy change with a trait?
//!
//! [`trait_effects_with`] fits, for each site, logit(occupancy) against a per-run trait, with a
//! random intercept for each biological replicate, so repeated injections of one animal or donor
//! are not counted as independent evidence. mzLib's `TraitDesign` has six replicates at three
//! doses, two injections each, built so that logit(occupancy) = -2 + 0.5 x dose + a replicate
//! offset. Compute occupancy and write it to a file first; the runs table gives each run its
//! replicate, then the trait, then any covariates, each column named by its header (`run`,
//! `replicate`, `dose`, `injection`):
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::ptmqtl::{occupancy_with, trait_effects, trait_effects_with, OccupancyOptions, TraitEffectsOptions};
//!
//! let occ = occupancy_with(
//!     "trait.tsv",
//!     &OccupancyOptions { out: Some("trait_occupancy.tsv".into()), ..Default::default() },
//! )?;
//! assert_eq!((occ.written.as_ref().map(|w| w.row_count), occ.columns.is_empty()), (Some(24), true));
//!
//! let fit = trait_effects("trait_occupancy.tsv", "trait_runs.tsv")?;
//! assert_eq!((fit.trait_name.as_str(), fit.covariate_names.clone()), ("dose", vec!["injection".to_owned()]));
//! assert_eq!((fit.min_runs, fit.min_replicates), (6, 4)); // mzLib's defaults
//! let effects: Vec<(String, f64, Option<f64>, i64, i64)> = fit
//!     .rows()?
//!     .into_iter()
//!     .map(|r| (r.protein_accession, (r.effect.unwrap() * 100.0).round() / 100.0, r.df, r.runs, r.replicates))
//!     .collect();
//! assert_eq!(effects[0], ("CENSORED_SITE".to_owned(), 0.54, Some(3.0), 10, 5));
//! assert_eq!(effects[1], ("TRAIT_SITE".to_owned(), 0.5, Some(4.0), 12, 6));
//!
//! // Ask for six replicates, and CENSORED_SITE is below support.
//! let strict = trait_effects_with(
//!     "trait_occupancy.tsv",
//!     "trait_runs.tsv",
//!     &TraitEffectsOptions { min_replicates: Some(6), ..Default::default() },
//! )?;
//! let statuses: Vec<(String, String, Option<f64>)> =
//!     strict.rows()?.into_iter().map(|r| (r.protein_accession, r.status, r.effect.map(|e| (e * 10.0).round() / 10.0))).collect();
//! assert_eq!(statuses[0], ("CENSORED_SITE".to_owned(), "below_support".to_owned(), None));
//! assert_eq!(statuses[1], ("TRAIT_SITE".to_owned(), "fitted".to_owned(), Some(0.5)));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! `TRAIT_SITE` recovers the 0.5 it was built with. **`effect` is in logit units per unit of the
//! trait**, so `0.5_f64.exp()` = 1.65: each dose unit multiplies the odds of the modified form by
//! about 1.65. [`TraitEffects::scale`] names the scale, `natural_log_odds` by mzLib's default. The
//! degrees of freedom are between replicates (six replicates, less the intercept and the dose),
//! because the dose never varies within a replicate. `CENSORED_SITE` is the same site with one
//! replicate's two runs replaced by a floor and a ceiling. Both leave the fit — only quantified,
//! non-ceiling cells are responses — so it rests on five replicates and three degrees of freedom.
//!
//! **In log2 odds.** Natural-log odds are mzLib's default. [`TraitEffectsOptions::scale`] set to
//! [`OccupancyScale::Log2Odds`] reports the same fit in log2 odds (mzLib #1464), where an effect of
//! +1 doubles the odds of the modified form per trait unit:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::ptmqtl::{trait_effects, trait_effects_with, OccupancyScale, TraitEffectsOptions};
//!
//! let fit = trait_effects("trait_occupancy.tsv", "trait_runs.tsv")?;
//! let log2 = trait_effects_with(
//!     "trait_occupancy.tsv",
//!     "trait_runs.tsv",
//!     &TraitEffectsOptions { scale: Some(OccupancyScale::Log2Odds), ..Default::default() },
//! )?;
//! assert_eq!((fit.scale.as_str(), log2.scale.as_str()), ("natural_log_odds", "log2_odds"));
//! for (natural, base2) in fit.rows()?.iter().zip(log2.rows()?) {
//!     println!(
//!         "{:<14}{:.3} -> {:.3}  p {:.5} -> {:.5}",
//!         natural.protein_accession,
//!         natural.effect.unwrap(),
//!         base2.effect.unwrap(),
//!         natural.p_value.unwrap(),
//!         base2.p_value.unwrap()
//!     );
//! }
//! //   CENSORED_SITE 0.543 -> 0.783  p 0.00361 -> 0.00361
//! //   TRAIT_SITE    0.501 -> 0.723  p 0.00082 -> 0.00082
//! let trait_site = &log2.rows()?[1];
//! assert_eq!(format!("{:.3}", trait_site.effect.unwrap()), "0.723");
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! `TRAIT_SITE`'s 0.50 natural-log units are 0.72 log2 units: 2^0.72 = 1.65, the same odds ratio
//! per dose unit. mzLib does not refit. It divides the effects and standard errors by ln 2 and the
//! two variances by (ln 2)², so `df`, `p_value` and `q` are unchanged. Pick the scale your other
//! results use; the science does not change. What it does not tell you:
//!
//! - **Only the modified share.** Floors, count-only and undetected cells are left out, never
//!   entered as 0. Whether the modification is *detected* more often with the trait is a separate
//!   question this does not answer.
//! - **`q` covers one call**: Benjamini-Hochberg over the sites fitted in this call, not across
//!   datasets.
//! - **A trait that never varies, or a covariate that repeats it, cannot be fitted.** mzLib refuses
//!   the design, and the call returns [`crate::MzLibError::Usage`] with mzLib's reason.
//!
//! ## What these functions will not do
//!
//! - **Read MetaMorpheus or FlashLFQ output directly.** Build the table as above; mzLib has no
//!   adapter yet.
//! - **Combine datasets.** Pairs and trait effects are per dataset. mzLib can pool pairs across
//!   datasets (`GlobalPairEngine`), which no binding calls yet.
//! - **Place a C-terminal modification**, since the table carries no protein length.
//!
//! ## Cite
//!
//! The methods behind these numbers, from the verbs' specs:
//!
#![doc = include_str!("../docs/reference/cite.ptmqtl.md")]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};
use crate::readers::{Table, WrittenTable};

/// Every value of the occupancy `state` column, as mzLib's `OccupancyState` is written on the wire.
/// Only `"quantified"` cells have a `fraction`.
pub const OCCUPANCY_STATES: [&str; 4] = ["quantified", "floor", "count_only", "not_detected"];

/// The wire names of [`PairType`], in its order.
pub const PAIR_TYPES: [&str; 3] = ["physical", "covarying", "both"];

/// The named rules [`PairsOptions::run_level_group`] takes. For any other grouping, map each
/// modification to a group in a file and pass it as [`PairsOptions::run_level_groups`].
pub const RUN_LEVEL_GROUPS: [&str; 2] = ["category", "modification"];

/// Every value of the [`pairs_with`] `test_status` column (mzLib's `PairTestStatus`, #1468), for
/// co-varying pairs; a physical pair has `None`. Only a `"tested"` pair carries a statistic.
pub const PAIR_TEST_STATUSES: [&str; 4] = [
    "tested",
    "not_estimable",
    "too_few_shared_runs",
    "both_sites_mostly_bounds",
];

/// The wire names of [`OccupancyScale`], in its order: what [`TraitEffectsOptions::scale`] sends
/// and [`TraitEffects::scale`] reports (mzLib's `OccupancyScale`).
pub const OCCUPANCY_SCALES: [&str; 2] = ["natural_log_odds", "log2_odds"];

/// Every value of the [`trait_effects_with`] `status` column. Only `"fitted"` rows carry
/// statistics.
pub const TRAIT_STATUSES: [&str; 7] = [
    "fitted",
    "below_support",
    "too_few_observations",
    "rank_deficient",
    "too_few_groups",
    "variance_ratio_at_limit",
    "not_converged",
];

/// The pyMzLib release whose bridge first has the ptmqtl verbs (mzLib 1.0.595): their specs'
/// `since.pymzlib`.
const MZLIB_1_0_595_BRIDGE: &str = "0.6.0";

const TIMEOUT: Duration = Duration::from_secs(600);

// ---- options -----------------------------------------------------------------------------------

/// Which pairs [`pairs_with`] finds. There is no default: mzLib has none between its two methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PairType {
    /// `"physical"`: both modifications identified on one peptidoform, so on one molecule
    /// (`PtmPairEngine.Physical`). No test.
    Physical,
    /// `"covarying"`: two sites whose occupancies rise and fall together across runs
    /// (`PtmPairEngine.CoVarying`), with Spearman's rho and a Benjamini-Hochberg `q`.
    Covarying,
    /// `"both"`: physical rows first, then co-varying.
    Both,
}

impl PairType {
    /// The wire name, one of [`PAIR_TYPES`].
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Physical => "physical",
            Self::Covarying => "covarying",
            Self::Both => "both",
        }
    }
}

/// The unit [`trait_effects_with`] reports effects and variances in (mzLib's `OccupancyScale`,
/// #1464). There is no default here: leave [`TraitEffectsOptions::scale`] `None` and mzLib's
/// applies (`SiteTraitOptions.Scale`, natural-log odds at the pin).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum OccupancyScale {
    /// `"natural_log_odds"`: ln(f / (1 - f)), the logit, so `effect.exp()` is the odds ratio per
    /// trait unit.
    NaturalLogOdds,
    /// `"log2_odds"`: log2(f / (1 - f)), so an effect of +1 doubles the odds per trait unit. mzLib
    /// divides the natural-log fit's effects and standard errors by ln 2 and its variances by
    /// (ln 2)², and leaves df, p and q unchanged: a change of unit, not a refit.
    Log2Odds,
}

impl OccupancyScale {
    /// The wire name, one of [`OCCUPANCY_SCALES`].
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NaturalLogOdds => "natural_log_odds",
            Self::Log2Odds => "log2_odds",
        }
    }
}

/// Options for [`occupancy_with`]. The default sends nothing beyond the observations.
#[derive(Debug, Clone)]
pub struct OccupancyOptions {
    /// A TSV with a header row, then `run`, `sample`, to sum each sample's fractions before the
    /// ratio; the `run` column of the result then holds the sample. Map only fractions together,
    /// never technical replicates. `None`: per-run occupancy.
    pub samples: Option<PathBuf>,
    /// Count only these modifications as sites: exact names, `Category:IdWithMotif`. Not with
    /// [`Self::categories`]. `None`: no filter.
    pub modifications: Option<Vec<String>>,
    /// Count only modifications of these categories (the text before the first `:`), e.g.
    /// `"Common Biological"`. Not with [`Self::modifications`]. `None`: no filter.
    pub categories: Option<Vec<String>>,
    /// Write the table to this `.tsv` instead of returning it — the file [`trait_effects_with`]
    /// reads. Then [`SiteOccupancy::columns`] is empty and [`SiteOccupancy::written`] says where.
    pub out: Option<PathBuf>,
    /// How long to allow; `None` waits indefinitely. Ten minutes by default.
    pub timeout: Option<Duration>,
}

impl Default for OccupancyOptions {
    fn default() -> Self {
        Self {
            samples: None,
            modifications: None,
            categories: None,
            out: None,
            timeout: Some(TIMEOUT),
        }
    }
}

/// Options for [`pairs_with`]. The default sends nothing beyond the observations and the type, so
/// every setting is mzLib's.
#[derive(Debug, Clone)]
pub struct PairsOptions {
    /// [`OccupancyOptions::samples`]'s run-to-sample table; pairs are then across samples.
    pub samples: Option<PathBuf>,
    /// As [`OccupancyOptions::modifications`]: only these modifications are sites.
    pub modifications: Option<Vec<String>>,
    /// As [`OccupancyOptions::categories`]: only these categories are sites.
    pub categories: Option<Vec<String>>,
    /// Co-varying only. A site enters pairs only if quantified in at least this fraction of the
    /// runs (fraction of the runs, above 0 and at most 1), rounded up. `None` sends nothing, so
    /// mzLib's default applies (`PtmPairEngine.DefaultMinQuantifiedFraction`, 0.7 at the pin),
    /// which [`PtmPairs::min_quantified_fraction`] reports.
    pub min_quantified_fraction: Option<f64>,
    /// Co-varying only. Keep runs where a site was seen only modified. `false` sends nothing, so
    /// mzLib's default (leave them out) applies.
    pub keep_ceiling: bool,
    /// Co-varying only. Let a run where a site was seen in only one form enter as a tied bound,
    /// never an invented value (mzLib's `includeBounds`, #1468): `not_detected` and `floor` cells
    /// rank below every fraction of the site, ceilings and count-only cells seen only modified
    /// above; a count-only cell with both forms identified has no rank. A pair whose two sites are
    /// both mostly one bound is not tested (`"both_sites_mostly_bounds"`). `false` sends nothing, so
    /// mzLib's default applies (`includeBounds` false at the pin: no bounds), which
    /// [`PtmPairs::include_bounds`] reports.
    pub include_bounds: bool,
    /// Co-varying only. A TSV with a header row whose second field names the covariate (e.g.
    /// `target_psms`), then `run`, value: each run's (or with [`Self::samples`], each sample's)
    /// depth, in any unit. Pairs become the partial Spearman correlation given it (mzLib's
    /// `runCovariate`, #1468); it is ranked, so log2 gives the same answer. Leave a run out for
    /// "unknown": it is then left out of every pair. A blank, `NA` or `NaN` value is refused by
    /// mzLib. `None` sends nothing: no covariate.
    pub run_covariate: Option<PathBuf>,
    /// Co-varying only. Correlate given each run's level — the median occupancy of the other tested
    /// sites of the pair's group — to remove a batch or loading shift. **Batch only**: where the
    /// level follows the design, it removes the effect you study. `false` sends nothing, so mzLib's
    /// default (not removed) applies.
    pub remove_run_level: bool,
    /// With [`Self::remove_run_level`]: `"category"` or `"modification"`, one of
    /// [`RUN_LEVEL_GROUPS`], for one level per group of sites. `None`: one level for every site,
    /// mzLib's default. Not with [`Self::run_level_groups`].
    pub run_level_group: Option<String>,
    /// With [`Self::remove_run_level`]: a TSV with a header row, then `modification`, `group`, for
    /// any grouping, such as chemistry. Every modification among the sites needs a group.
    pub run_level_groups: Option<PathBuf>,
    /// How long to allow; `None` waits indefinitely. Ten minutes by default.
    pub timeout: Option<Duration>,
}

impl Default for PairsOptions {
    fn default() -> Self {
        Self {
            samples: None,
            modifications: None,
            categories: None,
            min_quantified_fraction: None,
            keep_ceiling: false,
            include_bounds: false,
            run_covariate: None,
            remove_run_level: false,
            run_level_group: None,
            run_level_groups: None,
            timeout: Some(TIMEOUT),
        }
    }
}

/// Options for [`trait_effects_with`]. The default sends nothing beyond the two files, so every
/// threshold is mzLib's.
#[derive(Debug, Clone)]
pub struct TraitEffectsOptions {
    /// Fewest runs with a response for a site to be fitted, in runs. `None` sends nothing, so
    /// mzLib's default applies (`SiteTraitOptions.MinRuns`, 6 at the pin), which
    /// [`TraitEffects::min_runs`] reports.
    pub min_runs: Option<u32>,
    /// Fewest distinct replicates among them, in replicates. `None` sends nothing, so mzLib's
    /// default applies (`SiteTraitOptions.MinReplicates`, 4 at the pin).
    pub min_replicates: Option<u32>,
    /// Fractions are clamped to [epsilon, 1 - epsilon] before the logit, a fraction (0 to 1) above 0
    /// and below 0.5. `None` sends nothing, so mzLib's default applies
    /// (`SiteTraitOptions.LogitEpsilon`, 5e-5 at the pin).
    pub logit_epsilon: Option<f64>,
    /// Use runs where the site was seen only modified. `false` sends nothing, so mzLib's default
    /// (leave them out) applies.
    pub keep_ceiling: bool,
    /// The unit of the effects and variances. [`OccupancyScale::NaturalLogOdds`]: `effect.exp()` is
    /// the odds ratio per trait unit. [`OccupancyScale::Log2Odds`]: an effect of +1 doubles the
    /// odds; mzLib divides the natural-log fit by ln 2 (variances by (ln 2)²), and df, p and q do
    /// not change. `None` sends nothing, so mzLib's default applies (`SiteTraitOptions.Scale`,
    /// natural-log odds at the pin), which [`TraitEffects::scale`] reports.
    pub scale: Option<OccupancyScale>,
    /// How long to allow; `None` waits indefinitely. Ten minutes by default.
    pub timeout: Option<Duration>,
}

impl Default for TraitEffectsOptions {
    fn default() -> Self {
        Self {
            min_runs: None,
            min_replicates: None,
            logit_epsilon: None,
            keep_ceiling: false,
            scale: None,
            timeout: Some(TIMEOUT),
        }
    }
}

// ---- result types ------------------------------------------------------------------------------

/// What [`occupancy_with`] returns: one row per (site, run) where the site was identified or
/// covered.
///
/// Rows are in mzLib's order: by protein, then run, then position and modification. [`Self::rows`]
/// gives them as typed [`OccupancyRow`]s; their columns and units are in the reference below.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct SiteOccupancy {
    /// Absolute path of the observations table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub observations_file: String,
    /// Absolute path of the run-to-sample table, or `None` without [`OccupancyOptions::samples`]:
    /// occupancy is per run.
    #[serde(default)]
    pub samples_file: Option<String>,
    /// Observation rows read, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub observation_count: u64,
    /// Distinct runs in the observations, in runs.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub run_count: u64,
    /// Distinct samples, in samples, or `None` without [`OccupancyOptions::samples`].
    #[serde(default)]
    pub sample_count: Option<u64>,
    /// The distinct values of the `run` column, in ordinal order: the runs, or the samples with
    /// [`OccupancyOptions::samples`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub run_names: Vec<String>,
    /// Distinct sites in the table, in sites.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub site_count: u64,
    /// Cells per state, in cells, keyed by [`OCCUPANCY_STATES`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub state_counts: BTreeMap<String, u64>,
    /// The modification names the sites were limited to, or `None` for no filter.
    #[serde(default)]
    pub modifications: Option<Vec<String>>,
    /// The categories the sites were limited to, or `None` for no filter.
    #[serde(default)]
    pub categories: Option<Vec<String>>,
    /// Columns of the observations file that are not fields of mzLib's observation, so were not
    /// read. A misspelled `previous_residue` shows up here.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub ignored_columns: Vec<String>,
    /// Cells in the table, in rows, whether returned or written.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per (site, run); empty when it was written to [`OccupancyOptions::out`]
    /// instead. Its column order is also the header `out` writes.
    #[serde(flatten)]
    pub columns: Table,
    /// Where [`OccupancyOptions::out`] wrote the table, or `None` without it.
    #[serde(default)]
    pub written: Option<WrittenTable>,
    /// What the numbers are and are not, plus any that apply to this input.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of [`SiteOccupancy`]: one site in one run (or sample).
#[derive(Debug, Clone, PartialEq)]
pub struct OccupancyRow {
    /// `accession:residue position:modification`, the key [`pairs_with`] and
    /// [`trait_effects_with`] use.
    pub site_key: String,
    /// The protein, as given.
    pub protein_accession: String,
    /// The 1-based residue in the protein; a protein N-terminal modification at the residue it sits
    /// on.
    pub position: i64,
    /// The residue, one letter.
    pub residue: String,
    /// `Category:IdWithMotif`, as in the full sequence.
    pub modification: String,
    /// The text before the first `:` (`""` if none).
    pub category: String,
    /// The run, or with [`OccupancyOptions::samples`] the sample.
    pub run: String,
    /// One of [`OCCUPANCY_STATES`].
    pub state: String,
    /// Summed intensity of the quantified forms carrying the modification here, in intensity
    /// (instrument units); 0 unless quantified.
    pub modified_intensity: f64,
    /// Summed intensity of every quantified form covering the position, in intensity (instrument
    /// units); 0 for `count_only`.
    pub covering_intensity: f64,
    /// Modified over covering intensity, a fraction (0 to 1). `None` unless [`Self::state`] is
    /// `"quantified"`: the occupancy is censored (floor) or unknown, **never 0**.
    pub fraction: Option<f64>,
    /// `false` with a fraction of 1 is a **ceiling**: the site was seen only modified.
    pub unmodified_quantified: bool,
    /// Distinct peptidoforms covering the position identified in this run (or sample), quantified
    /// or not, in peptidoforms: mzLib's `SiteRunOccupancy.CoveringPeptidoforms` (#1468), the n of a
    /// cell's k of n. It counts identifications, not intensity. `None`: not supplied (mzLib always
    /// supplies it; a bridge before pyMzLib 0.7.0 has no such column).
    pub covering_peptidoforms: Option<i64>,
    /// Of [`Self::covering_peptidoforms`], those carrying this modification at the position, in
    /// peptidoforms: the k. 0 for `not_detected`. `None`: not supplied, as
    /// [`Self::covering_peptidoforms`].
    pub modified_peptidoforms: Option<i64>,
}

impl SiteOccupancy {
    /// Every row, typed, in mzLib's order. Empty when the table was written to
    /// [`OccupancyOptions::out`] instead.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self) -> Result<Vec<OccupancyRow>> {
        let t = &self.columns;
        if t.is_empty() {
            return Ok(Vec::new());
        }
        let site_key = t.strings("site_key")?;
        let protein_accession = t.strings("protein_accession")?;
        let position = t.integers("position")?;
        let residue = t.strings("residue")?;
        let modification = t.strings("modification")?;
        let category = t.strings("category")?;
        let run = t.strings("run")?;
        let state = t.strings("state")?;
        let modified_intensity = t.floats("modified_intensity")?;
        let covering_intensity = t.floats("covering_intensity")?;
        let fraction = t.floats("fraction")?;
        let unmodified_quantified = t.booleans("unmodified_quantified")?;
        let covering_peptidoforms = optional_integers(t, "covering_peptidoforms", site_key.len())?;
        let modified_peptidoforms = optional_integers(t, "modified_peptidoforms", site_key.len())?;
        Ok((0..site_key.len())
            .map(|i| OccupancyRow {
                site_key: site_key[i].clone().unwrap_or_default(),
                protein_accession: protein_accession[i].clone().unwrap_or_default(),
                position: position[i].unwrap_or_default(),
                residue: residue[i].clone().unwrap_or_default(),
                modification: modification[i].clone().unwrap_or_default(),
                category: category[i].clone().unwrap_or_default(),
                run: run[i].clone().unwrap_or_default(),
                state: state[i].clone().unwrap_or_default(),
                modified_intensity: modified_intensity[i].unwrap_or_default(),
                covering_intensity: covering_intensity[i].unwrap_or_default(),
                fraction: fraction[i],
                unmodified_quantified: unmodified_quantified[i].unwrap_or_default(),
                covering_peptidoforms: covering_peptidoforms[i],
                modified_peptidoforms: modified_peptidoforms[i],
            })
            .collect())
    }
}

/// A column a bridge before pyMzLib 0.7.0 (mzLib 1.0.596) does not write: all `None` when absent.
fn optional_integers(t: &Table, name: &str, rows: usize) -> Result<Vec<Option<i64>>> {
    if t.has(name) {
        t.integers(name)
    } else {
        Ok(vec![None; rows])
    }
}

/// As [`optional_integers`], for a text column.
fn optional_strings(t: &Table, name: &str, rows: usize) -> Result<Vec<Option<String>>> {
    if t.has(name) {
        t.strings(name)
    } else {
        Ok(vec![None; rows])
    }
}

/// As [`optional_integers`], for a true/false column.
fn optional_booleans(t: &Table, name: &str, rows: usize) -> Result<Vec<Option<bool>>> {
    if t.has(name) {
        t.booleans(name)
    } else {
        Ok(vec![None; rows])
    }
}

/// What [`pairs_with`] returns: one row per pair of sites, physical pairs first.
///
/// [`Self::rows`] gives the rows as typed [`PtmPair`]s; their columns and units are in the
/// reference below. Every co-varying setting is `None` for [`PairType::Physical`].
#[derive(Debug, Clone, Default, Deserialize)]
pub struct PtmPairs {
    /// Absolute path of the observations table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub observations_file: String,
    /// Absolute path of the run-to-sample table, or `None` without [`PairsOptions::samples`]: pairs
    /// are across runs.
    #[serde(default)]
    pub samples_file: Option<String>,
    /// `"physical"`, `"covarying"` or `"both"`, as you asked.
    #[serde(rename = "type", default, deserialize_with = "bridge::null_to_default")]
    pub r#type: String,
    /// Observation rows read, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub observation_count: u64,
    /// Distinct runs in the observations, in runs.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub run_count: u64,
    /// Distinct samples, in samples, or `None` without [`PairsOptions::samples`].
    #[serde(default)]
    pub sample_count: Option<u64>,
    /// The runs, or samples, the pairs were computed over, in ordinal order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub run_names: Vec<String>,
    /// The runs (or samples) where any site was identified or covered, in runs: the denominator of
    /// [`Self::min_quantified_fraction`]. `None` for [`PairType::Physical`].
    #[serde(default)]
    pub covarying_run_count: Option<u64>,
    /// The guard used, as a fraction of the runs (0 to 1): yours, or mzLib's default. `None` for
    /// [`PairType::Physical`].
    #[serde(default)]
    pub min_quantified_fraction: Option<f64>,
    /// `Some(false)` with [`PairsOptions::keep_ceiling`]. `None` when you did not choose: mzLib's
    /// `CoVarying` then resolves its own default (ceilings left out) since mzLib #1468, and the
    /// crate does not restate that rule to name the value. Also `None` for [`PairType::Physical`].
    #[serde(default)]
    pub exclude_ceiling: Option<bool>,
    /// Whether runs where a site was seen in only one form entered as tied bounds: `true` with
    /// [`PairsOptions::include_bounds`], else mzLib's default (`false`). `None` for
    /// [`PairType::Physical`].
    #[serde(default)]
    pub include_bounds: Option<bool>,
    /// mzLib's `PtmPairEngine.MinUntiedShare`, a fraction of the shared runs (0 to 1): with bounds,
    /// a pair is tested only if one of its sites has at least this share of the shared runs outside
    /// its largest tie. `None` when bounds were not included (or for [`PairType::Physical`]).
    #[serde(default)]
    pub min_untied_share: Option<f64>,
    /// Absolute path of the run covariate table, or `None` without
    /// [`PairsOptions::run_covariate`].
    #[serde(default)]
    pub run_covariate_file: Option<String>,
    /// The covariate's name, the run covariate table's second header, or `None` without one.
    #[serde(default)]
    pub run_covariate_name: Option<String>,
    /// Runs (or samples) with no value in the run covariate table, in ordinal order: mzLib left
    /// them out of every pair. Empty when every run has one; `None` without a run covariate.
    #[serde(default)]
    pub run_covariate_runs_missing: Option<Vec<String>>,
    /// Whether the run level was removed; `None` for [`PairType::Physical`].
    #[serde(default)]
    pub remove_run_level: Option<bool>,
    /// `"none"` (one group), `"category"`, `"modification"` or `"file"`; `None` when the run level
    /// was not removed (or for [`PairType::Physical`]).
    #[serde(default)]
    pub run_level_group: Option<String>,
    /// Absolute path of the modification-to-group table, or `None` without
    /// [`PairsOptions::run_level_groups`].
    #[serde(default)]
    pub run_level_groups_file: Option<String>,
    /// mzLib's `MinRunsForRunLevel`, in runs: with the run level or a run covariate removed, a pair
    /// sharing fewer runs is not tested ([`PtmPair::test_status`] `"too_few_shared_runs"`). `None`
    /// when neither was removed.
    #[serde(default)]
    pub min_runs_for_run_level: Option<u64>,
    /// mzLib's `MinSitesForRunLevel`, in sites: a run whose group has fewer other sites there has no
    /// level. `None` when the run level was not removed.
    #[serde(default)]
    pub min_sites_for_run_level: Option<u64>,
    /// The modification names the sites were limited to, or `None` for no filter.
    #[serde(default)]
    pub modifications: Option<Vec<String>>,
    /// The categories the sites were limited to, or `None` for no filter.
    #[serde(default)]
    pub categories: Option<Vec<String>>,
    /// Columns of the observations file that were not read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub ignored_columns: Vec<String>,
    /// Physical pairs found, in pairs, or `None` for [`PairType::Covarying`].
    #[serde(default)]
    pub physical_count: Option<u64>,
    /// Co-varying pairs reported, in pairs, or `None` for [`PairType::Physical`].
    #[serde(default)]
    pub covarying_count: Option<u64>,
    /// Pairs in the table, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per pair.
    #[serde(flatten)]
    pub columns: Table,
    /// What each pair type means, and the run level when it was removed, plus any that apply to
    /// this input, among them one counting untested pairs by [`PtmPair::test_status`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of [`PtmPairs`]: one pair of sites.
#[derive(Debug, Clone, PartialEq)]
pub struct PtmPair {
    /// `"physical"` (both on one peptidoform, so one molecule) or `"covarying"` (across runs).
    pub result_type: String,
    /// The first site's key, first by ordinal order; a pair is written once.
    pub site_a: String,
    /// The second site's key.
    pub site_b: String,
    /// Both sites are on one protein accession.
    pub same_protein: bool,
    /// One peptide covers both sites; such a co-varying pair is correlated by construction.
    pub overlapping: bool,
    /// Physical: the median co-occupancy, a fraction (0 to 1); co-varying: Spearman's rho, a
    /// unitless correlation (-1 to 1). `None`: physical pairs never quantified together in any run,
    /// or a co-varying pair that is not estimable.
    pub statistic: Option<f64>,
    /// Two-sided Spearman p-value, a fraction (0 to 1). `None` for a physical pair (no test) or one
    /// that is not estimable.
    pub p_value: Option<f64>,
    /// Benjamini-Hochberg within [`Self::fdr_family`], a fraction (0 to 1). `None` = not tested
    /// (physical, overlapping, not estimable). Not a target-decoy q-value.
    pub q: Option<f64>,
    /// mzLib's family: `P` or `A`, then `intra` or `inter` protein, e.g. `"A:inter"`.
    pub fdr_family: String,
    /// Physical: runs identified together; co-varying: runs where both were quantified. In runs.
    pub n: i64,
    /// Runs [`Self::statistic`] was computed on, in runs.
    pub statistic_n: i64,
    /// `"exact"`, `"asymptotic"` or `"not_estimable"` (not tested, for the reason
    /// [`Self::test_status`] gives); `None` for a physical pair.
    pub spearman_method: Option<String>,
    /// Whether rho is partial given the run level.
    pub run_level_removed: bool,
    /// One of [`PAIR_TEST_STATUSES`], mzLib's `PairTestStatus` (#1468): `"tested"` (statistic,
    /// p-value and, unless overlapping, `q` are set) or why not — `"not_estimable"` (too few
    /// complete runs, a site constant over them, or a covariate that carries a site entirely),
    /// `"too_few_shared_runs"` (partial on the run level or a run covariate, on fewer than
    /// [`PtmPairs::min_runs_for_run_level`] runs) or `"both_sites_mostly_bounds"`. `None` for a
    /// physical pair, which has no test, and from a bridge before pyMzLib 0.7.0.
    pub test_status: Option<String>,
    /// Whether one-form cells entered as tied bounds. Always `false` for a physical pair.
    pub bounds_included: bool,
    /// Whether rho is partial given the run covariate. Always `false` for a physical pair.
    pub run_covariate_removed: bool,
}

impl PtmPairs {
    /// Every row, typed, physical pairs first.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self) -> Result<Vec<PtmPair>> {
        let t = &self.columns;
        if t.is_empty() {
            return Ok(Vec::new());
        }
        let result_type = t.strings("result_type")?;
        let site_a = t.strings("site_a")?;
        let site_b = t.strings("site_b")?;
        let same_protein = t.booleans("same_protein")?;
        let overlapping = t.booleans("overlapping")?;
        let statistic = t.floats("statistic")?;
        let p_value = t.floats("p_value")?;
        let q = t.floats("q")?;
        let fdr_family = t.strings("fdr_family")?;
        let n = t.integers("n")?;
        let statistic_n = t.integers("statistic_n")?;
        let spearman_method = t.strings("spearman_method")?;
        let run_level_removed = t.booleans("run_level_removed")?;
        let test_status = optional_strings(t, "test_status", site_a.len())?;
        let bounds_included = optional_booleans(t, "bounds_included", site_a.len())?;
        let run_covariate_removed = optional_booleans(t, "run_covariate_removed", site_a.len())?;
        Ok((0..site_a.len())
            .map(|i| PtmPair {
                result_type: result_type[i].clone().unwrap_or_default(),
                site_a: site_a[i].clone().unwrap_or_default(),
                site_b: site_b[i].clone().unwrap_or_default(),
                same_protein: same_protein[i].unwrap_or_default(),
                overlapping: overlapping[i].unwrap_or_default(),
                statistic: statistic[i],
                p_value: p_value[i],
                q: q[i],
                fdr_family: fdr_family[i].clone().unwrap_or_default(),
                n: n[i].unwrap_or_default(),
                statistic_n: statistic_n[i].unwrap_or_default(),
                spearman_method: spearman_method[i].clone(),
                run_level_removed: run_level_removed[i].unwrap_or_default(),
                test_status: test_status[i].clone(),
                bounds_included: bounds_included[i].unwrap_or_default(),
                run_covariate_removed: run_covariate_removed[i].unwrap_or_default(),
            })
            .collect())
    }
}

/// What [`trait_effects_with`] returns: one row per site with a usable response, by site key.
///
/// [`Self::rows`] gives the rows as typed [`TraitEffect`]s; their columns and units are in the
/// reference below.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct TraitEffects {
    /// Absolute path of the occupancy table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub occupancy_file: String,
    /// Absolute path of the runs table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub runs_file: String,
    /// The runs table's third header: the trait.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub trait_name: String,
    /// The covariate headers, in order; empty without covariates.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub covariate_names: Vec<String>,
    /// Occupancy rows read, in cells.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub cell_count: u64,
    /// Runs in the runs table, in runs.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub run_count: u64,
    /// Runs with a trait value, in runs: the runs the design is built over.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub trait_run_count: u64,
    /// Occupancy runs with no trait value, so not used.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub occupancy_runs_not_used: Vec<String>,
    /// The support threshold used, in runs: yours, or mzLib's default.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub min_runs: u64,
    /// The replicate threshold used, in replicates: yours, or mzLib's default.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub min_replicates: u64,
    /// The clamp used before the logit, a fraction (0 to 1).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub logit_epsilon: f64,
    /// Whether ceiling cells were left out.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub exclude_ceiling: bool,
    /// The scale every effect, standard error and variance is in: `"natural_log_odds"` (mzLib's
    /// default, `SiteTraitOptions.Scale`; the logit, so `effect.exp()` is the odds ratio) or
    /// `"log2_odds"` (an effect of +1 doubles the odds; mzLib #1464), one of [`OCCUPANCY_SCALES`]:
    /// [`TraitEffectsOptions::scale`], or mzLib's default when it was `None`. Empty from a bridge
    /// before pyMzLib 0.7.0.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub scale: String,
    /// Sites per status, in sites, keyed by [`TRAIT_STATUSES`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub status_counts: BTreeMap<String, u64>,
    /// Columns of the occupancy file that were not read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub ignored_columns: Vec<String>,
    /// Sites in the table, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per site.
    #[serde(flatten)]
    pub columns: Table,
    /// The model and what it leaves out, plus any that apply to this input.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of [`TraitEffects`]: one site's mixed model.
///
/// Every statistic is `None` when [`Self::status`] is not `"fitted"`.
#[derive(Debug, Clone, PartialEq)]
pub struct TraitEffect {
    /// The site key, as [`occupancy_with`] writes it.
    pub site_key: String,
    /// The protein, as [`occupancy_with`] gives it.
    pub protein_accession: String,
    /// The 1-based residue in the protein.
    pub position: i64,
    /// The residue, one letter.
    pub residue: String,
    /// `Category:IdWithMotif`.
    pub modification: String,
    /// One of [`TRAIT_STATUSES`].
    pub status: String,
    /// The change in the log odds of occupancy per trait unit, on [`TraitEffects::scale`]:
    /// natural-log odds (the logit; `effect.exp()` is the odds ratio) by default.
    pub effect: Option<f64>,
    /// Its standard error, in log odds of occupancy per trait unit (as [`TraitEffects::scale`]).
    pub standard_error: Option<f64>,
    /// Degrees of freedom by the containment rule: between-replicate for a trait constant within
    /// each replicate.
    pub df: Option<f64>,
    /// Two-sided p-value for the trait, a fraction (0 to 1).
    pub p_value: Option<f64>,
    /// Benjamini-Hochberg over this call's fitted sites, a fraction (0 to 1); `None` when not
    /// fitted, so not in the family.
    pub q: Option<f64>,
    /// Runs contributing a response, in runs.
    pub runs: i64,
    /// Distinct replicates among them, in replicates.
    pub replicates: i64,
    /// Median occupancy over those runs, a fraction (0 to 1).
    pub median_fraction: Option<f64>,
    /// One per covariate, in [`TraitEffects::covariate_names`] order, in log odds of occupancy per
    /// covariate unit (as [`TraitEffects::scale`]); each `None` when not fitted.
    pub covariate_effects: Vec<Option<f64>>,
    /// The between-replicate variance (random intercept), in squared log odds of occupancy (as
    /// [`TraitEffects::scale`]).
    pub replicate_variance: Option<f64>,
    /// The within-replicate, run-to-run variance, in squared log odds of occupancy (as
    /// [`TraitEffects::scale`]).
    pub residual_variance: Option<f64>,
}

impl TraitEffects {
    /// Every row, typed, by site key.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self) -> Result<Vec<TraitEffect>> {
        let t = &self.columns;
        if t.is_empty() {
            return Ok(Vec::new());
        }
        let site_key = t.strings("site_key")?;
        let protein_accession = t.strings("protein_accession")?;
        let position = t.integers("position")?;
        let residue = t.strings("residue")?;
        let modification = t.strings("modification")?;
        let status = t.strings("status")?;
        let effect = t.floats("effect")?;
        let standard_error = t.floats("standard_error")?;
        let df = t.floats("df")?;
        let p_value = t.floats("p_value")?;
        let q = t.floats("q")?;
        let runs = t.integers("runs")?;
        let replicates = t.integers("replicates")?;
        let median_fraction = t.floats("median_fraction")?;
        let covariate_effects = t.float_arrays("covariate_effects")?;
        let replicate_variance = t.floats("replicate_variance")?;
        let residual_variance = t.floats("residual_variance")?;
        Ok((0..site_key.len())
            .map(|i| TraitEffect {
                site_key: site_key[i].clone().unwrap_or_default(),
                protein_accession: protein_accession[i].clone().unwrap_or_default(),
                position: position[i].unwrap_or_default(),
                residue: residue[i].clone().unwrap_or_default(),
                modification: modification[i].clone().unwrap_or_default(),
                status: status[i].clone().unwrap_or_default(),
                effect: effect[i],
                standard_error: standard_error[i],
                df: df[i],
                p_value: p_value[i],
                q: q[i],
                runs: runs[i].unwrap_or_default(),
                replicates: replicates[i].unwrap_or_default(),
                median_fraction: median_fraction[i],
                covariate_effects: covariate_effects[i].clone().unwrap_or_default(),
                replicate_variance: replicate_variance[i],
                residual_variance: residual_variance[i],
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

/// A list of names as the wire's one tab-separated option.
fn names(values: &[String], what: &str) -> Result<String> {
    if values.is_empty() {
        return Err(MzLibError::Usage(format!(
            "{what} takes a non-empty list of names, or None for no filter."
        )));
    }
    for name in values {
        if name.trim().is_empty() {
            return Err(MzLibError::Usage(format!(
                "Every name in {what} must be non-empty; got {name:?}."
            )));
        }
        if name.contains(['\t', '\r', '\n']) {
            return Err(MzLibError::Usage(format!(
                "A name in {what} contains a tab or a line break: {name:?}."
            )));
        }
    }
    Ok(values.join("\t"))
}

fn filter_args(
    modifications: Option<&Vec<String>>,
    categories: Option<&Vec<String>>,
) -> Result<Vec<String>> {
    match (modifications, categories) {
        (Some(_), Some(_)) => Err(MzLibError::Usage(
            "Give modifications or categories, not both: mzLib takes one filter on the \
             modification name."
                .to_owned(),
        )),
        (Some(values), None) => Ok(vec![
            "--modifications".to_owned(),
            names(values, "modifications")?,
        ]),
        (None, Some(values)) => Ok(vec![
            "--categories".to_owned(),
            names(values, "categories")?,
        ]),
        (None, None) => Ok(Vec::new()),
    }
}

fn occupancy_request(observations: &Path, options: &OccupancyOptions) -> Result<Vec<String>> {
    let mut args = vec![
        "ptmqtl".to_owned(),
        "occupancy".to_owned(),
        "--observations".to_owned(),
        path_text(observations, "observations")?,
    ];
    if let Some(samples) = &options.samples {
        args.push("--samples".to_owned());
        args.push(path_text(samples, "samples")?);
    }
    args.extend(filter_args(
        options.modifications.as_ref(),
        options.categories.as_ref(),
    )?);
    if let Some(out) = &options.out {
        args.push("--out".to_owned());
        args.push(path_text(out, "out")?);
    }
    Ok(args)
}

fn pairs_request(
    observations: &Path,
    pair_type: PairType,
    options: &PairsOptions,
) -> Result<Vec<String>> {
    let mut args = vec![
        "ptmqtl".to_owned(),
        "pairs".to_owned(),
        "--observations".to_owned(),
        path_text(observations, "observations")?,
        "--type".to_owned(),
        pair_type.as_str().to_owned(),
    ];
    if let Some(samples) = &options.samples {
        args.push("--samples".to_owned());
        args.push(path_text(samples, "samples")?);
    }
    args.extend(filter_args(
        options.modifications.as_ref(),
        options.categories.as_ref(),
    )?);
    if let Some(fraction) = options.min_quantified_fraction {
        args.push("--min-quantified-fraction".to_owned());
        args.push(bridge::format_number(fraction, "min_quantified_fraction")?);
    }
    if options.keep_ceiling {
        args.push("--keep-ceiling".to_owned());
    }
    if options.include_bounds {
        args.push("--include-bounds".to_owned());
    }
    if let Some(covariate) = &options.run_covariate {
        args.push("--run-covariate".to_owned());
        args.push(path_text(covariate, "run_covariate")?);
    }
    if options.remove_run_level {
        args.push("--remove-run-level".to_owned());
    }
    if let Some(group) = &options.run_level_group {
        if !RUN_LEVEL_GROUPS.contains(&group.as_str()) {
            return Err(MzLibError::Usage(format!(
                "run_level_group must be one of {}; got {group:?}. For any other grouping, pass a \
                 modification-to-group file as run_level_groups.",
                RUN_LEVEL_GROUPS.join(", ")
            )));
        }
        args.push("--run-level-group".to_owned());
        args.push(group.clone());
    }
    if let Some(groups) = &options.run_level_groups {
        args.push("--run-level-groups".to_owned());
        args.push(path_text(groups, "run_level_groups")?);
    }
    Ok(args)
}

fn trait_effects_request(
    occupancy: &Path,
    runs: &Path,
    options: &TraitEffectsOptions,
) -> Result<Vec<String>> {
    let mut args = vec![
        "ptmqtl".to_owned(),
        "trait-effects".to_owned(),
        "--occupancy".to_owned(),
        path_text(occupancy, "occupancy")?,
        "--runs".to_owned(),
        path_text(runs, "runs")?,
    ];
    if let Some(minimum) = options.min_runs {
        args.push("--min-runs".to_owned());
        args.push(minimum.to_string());
    }
    if let Some(minimum) = options.min_replicates {
        args.push("--min-replicates".to_owned());
        args.push(minimum.to_string());
    }
    if let Some(epsilon) = options.logit_epsilon {
        args.push("--logit-epsilon".to_owned());
        args.push(bridge::format_number(epsilon, "logit_epsilon")?);
    }
    if options.keep_ceiling {
        args.push("--keep-ceiling".to_owned());
    }
    if let Some(scale) = options.scale {
        args.push("--scale".to_owned());
        args.push(scale.as_str().to_owned());
    }
    Ok(args)
}

fn call<T: serde::de::DeserializeOwned>(args: &[String], timeout: Option<Duration>) -> Result<T> {
    let data = bridge::invoke(args, None, timeout)?;
    serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!("ptmqtl payload could not be interpreted: {error}"))
    })
}

// ---- the verbs ---------------------------------------------------------------------------------

/// [`occupancy_with`] per run, with no filter, returning the table.
///
/// # Errors
///
/// As [`occupancy_with`].
pub fn occupancy(observations: impl AsRef<Path>) -> Result<SiteOccupancy> {
    occupancy_with(observations, &OccupancyOptions::default())
}

/// Compute every site's occupancy in every run, with mzLib's `SiteOccupancyCalculator`.
///
/// For a site (protein, position, modification) in a run, occupancy is the summed intensity of the
/// quantified peptidoforms that carry the modification there, over the summed intensity of every
/// quantified peptidoform that covers the position. Each cell keeps its state, so a censored cell
/// is never mistaken for a 0: see [`OCCUPANCY_STATES`]. mzLib always leaves out the categories
/// `Common Variable` and `Common Fixed`, peptide-terminal modifications and every C-terminal
/// modification.
///
/// `observations` is a TSV with a header row naming its columns, one row per (run, peptidoform,
/// protein) identified by MS/MS: `run`, `full_sequence` (MetaMorpheus notation),
/// `protein_accession`, `start_residue` and `end_residue` (1-based, in the protein), `intensity`
/// (instrument units; blank or `NA` for identified but not quantified), and optionally
/// `previous_residue` (one letter). A shared peptide is one row per protein. Never one row per PSM,
/// and never a match-between-runs transfer.
///
/// Needs the bridge from pyMzLib 0.6.0 or later.
#[doc = include_str!("../docs/reference/ptmqtl.occupancy.md")]
///
/// # Examples
///
/// Written to a file, for [`trait_effects_with`] to read:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::ptmqtl::{occupancy_with, OccupancyOptions};
///
/// let occ = occupancy_with(
///     "trait.tsv",
///     &OccupancyOptions { out: Some("trait_occupancy.tsv".into()), ..Default::default() },
/// )?;
/// let written = occ.written.as_ref().expect("out was given");
/// assert_eq!((written.format.as_str(), written.row_count), ("tsv", 24));
/// assert!(occ.columns.is_empty() && occ.rows()?.is_empty());
/// assert_eq!((occ.state_counts["quantified"], occ.state_counts["floor"]), (23, 1));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when both filters are given, a filter is an
/// empty list or holds a blank name or one with a tab or line break, a path is blank or not UTF-8,
/// or the bridge predates the verb.
#[doc = include_str!("../docs/reference/ptmqtl.occupancy.see-also.md")]
pub fn occupancy_with(
    observations: impl AsRef<Path>,
    options: &OccupancyOptions,
) -> Result<SiteOccupancy> {
    let args = occupancy_request(observations.as_ref(), options)?;
    bridge::require_verb("ptmqtl occupancy", MZLIB_1_0_595_BRIDGE)?;
    call(&args, options.timeout)
}

/// [`pairs_with`] with every setting mzLib's.
///
/// # Errors
///
/// As [`pairs_with`].
pub fn pairs(observations: impl AsRef<Path>, pair_type: PairType) -> Result<PtmPairs> {
    pairs_with(observations, pair_type, &PairsOptions::default())
}

/// Find PTM-PTM pairs within one dataset, with mzLib's `PtmPairEngine`.
///
/// Two kinds of evidence, which you choose between with `pair_type` (it is required: mzLib has no
/// default between its two methods):
///
/// - **physical** (`PtmPairEngine.Physical`): both modifications were identified on one
///   peptidoform, so on one molecule. `statistic` is the median co-occupancy over the runs where
///   both were quantified. There is no test.
/// - **co-varying** (`PtmPairEngine.CoVarying`): the two sites' occupancies rise and fall together
///   (or oppositely) across runs. `statistic` is Spearman's rho on the runs where both are
///   quantified, with a Benjamini-Hochberg `q` per family. Floors, undetected cells and, by
///   default, ceilings are left out: two unrelated sites both read 1 in low-load runs, which would
///   correlate them through detection alone. [`PairsOptions::include_bounds`] lets them in as tied
///   bounds instead, and [`PairsOptions::run_covariate`] correlates given each run's depth.
///
/// `observations` is [`occupancy_with`]'s table, with the same rules; occupancy is computed inside
/// by the same mzLib call.
///
/// Needs the bridge from pyMzLib 0.6.0 or later; [`PairsOptions::include_bounds`] and
/// [`PairsOptions::run_covariate`] need pyMzLib 0.7.0's.
#[doc = include_str!("../docs/reference/ptmqtl.pairs.md")]
///
/// # Examples
///
/// The physical pair on mzLib's cases file, with every co-varying setting `None`:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::ptmqtl::{pairs, PairType};
///
/// let found = pairs("cases.tsv", PairType::Physical)?;
/// assert_eq!((found.r#type.as_str(), found.physical_count, found.covarying_count), ("physical", Some(1), None));
/// assert_eq!(found.min_quantified_fraction, None); // a co-varying setting
/// let pair = &found.rows()?[0];
/// assert_eq!((pair.result_type.as_str(), pair.fdr_family.as_str()), ("physical", "P:intra"));
/// assert_eq!((pair.same_protein, pair.q, pair.spearman_method.clone()), (true, None, None));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when both filters are given, a filter is an
/// empty list or holds a blank name or one with a tab or line break,
/// [`PairsOptions::run_level_group`] is not one of [`RUN_LEVEL_GROUPS`],
/// [`PairsOptions::min_quantified_fraction`] is not finite, a path is blank or not UTF-8, or the
/// bridge predates the verb. A co-varying option with [`PairType::Physical`] (among them
/// [`PairsOptions::include_bounds`] and [`PairsOptions::run_covariate`]), a run-level grouping
/// without [`PairsOptions::remove_run_level`], and both groupings are the bridge's refusals; a run
/// covariate value that is not finite is mzLib's.
#[doc = include_str!("../docs/reference/ptmqtl.pairs.see-also.md")]
pub fn pairs_with(
    observations: impl AsRef<Path>,
    pair_type: PairType,
    options: &PairsOptions,
) -> Result<PtmPairs> {
    let args = pairs_request(observations.as_ref(), pair_type, options)?;
    bridge::require_verb("ptmqtl pairs", MZLIB_1_0_595_BRIDGE)?;
    call(&args, options.timeout)
}

/// [`trait_effects_with`] with mzLib's thresholds and ceiling rule.
///
/// # Errors
///
/// As [`trait_effects_with`].
pub fn trait_effects(occupancy: impl AsRef<Path>, runs: impl AsRef<Path>) -> Result<TraitEffects> {
    trait_effects_with(occupancy, runs, &TraitEffectsOptions::default())
}

/// Fit, per site, the effect of a per-run trait on the log odds of occupancy, with mzLib's
/// `SiteTraitEffects`.
///
/// The model, per site, is logit(occupancy) = b0 + b_trait x trait + covariates + a random
/// intercept per biological replicate + error, fitted by REML (mzLib's `MixedModel`). The random
/// intercept is the point: repeated runs of one animal or donor are not counted as independent, so
/// a trait that is constant within each replicate is tested at the between-replicate degrees of
/// freedom. Only quantified cells are responses: floors, count-only and undetected cells are left
/// out, never entered as 0, so this answers "does the modified share change?", not "is the
/// modification detected more often?". Ceilings are left out by default.
///
/// `occupancy` is the table [`occupancy_with`] writes with [`OccupancyOptions::out`] (columns found
/// by header name); a `fraction` cell is used as written. `runs` is a TSV with a header row, then
/// one row per run: the run, its biological replicate, the trait (its header names it; any unit),
/// then any covariates (each header names one). A run with a blank trait is not used.
///
/// Needs the bridge from pyMzLib 0.6.0 or later; [`TraitEffectsOptions::scale`] needs pyMzLib
/// 0.7.0's.
#[doc = include_str!("../docs/reference/ptmqtl.trait-effects.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::ptmqtl::trait_effects;
///
/// let fit = trait_effects("trait_occupancy.tsv", "trait_runs.tsv")?;
/// assert_eq!(fit.status_counts["fitted"], 2);
/// assert_eq!((fit.logit_epsilon, fit.exclude_ceiling), (5e-5, true)); // mzLib's defaults
/// assert_eq!(fit.scale, "natural_log_odds"); // the logit
/// let site = &fit.rows()?[1];
/// assert_eq!(site.protein_accession, "TRAIT_SITE");
/// assert_eq!(site.covariate_effects.len(), fit.covariate_names.len());
/// assert!(site.q.is_some_and(|q| q < 0.01));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`], before anything is spawned, when [`TraitEffectsOptions::logit_epsilon`]
/// is not finite, a path is blank or not UTF-8, or the bridge predates the verb.
#[doc = include_str!("../docs/reference/ptmqtl.trait-effects.see-also.md")]
pub fn trait_effects_with(
    occupancy: impl AsRef<Path>,
    runs: impl AsRef<Path>,
    options: &TraitEffectsOptions,
) -> Result<TraitEffects> {
    let args = trait_effects_request(occupancy.as_ref(), runs.as_ref(), options)?;
    bridge::require_verb("ptmqtl trait-effects", MZLIB_1_0_595_BRIDGE)?;
    call(&args, options.timeout)
}

#[cfg(test)]
mod tests {
    //! Offline. The recordings are pyMzLib's, made from the real bridge.

    use super::*;

    const CASES: &str = include_str!("../tests/fixtures/ptmqtl_occupancy_cases.json");
    const FRACTIONS: &str = include_str!("../tests/fixtures/ptmqtl_occupancy_fractions.json");
    const WRITTEN: &str = include_str!("../tests/fixtures/ptmqtl_occupancy_trait_out.json");
    const PHYSICAL: &str = include_str!("../tests/fixtures/ptmqtl_pairs_physical.json");
    const BATCH: &str = include_str!("../tests/fixtures/ptmqtl_pairs_batch.json");
    const RUN_LEVEL: &str = include_str!("../tests/fixtures/ptmqtl_pairs_batch_run_level.json");
    const EFFECTS: &str = include_str!("../tests/fixtures/ptmqtl_trait_effects.json");
    const STRICT: &str = include_str!("../tests/fixtures/ptmqtl_trait_effects_strict.json");
    const LOG2: &str = include_str!("../tests/fixtures/ptmqtl_trait_effects_log2.json");
    const BOUNDS_PLAIN: &str = include_str!("../tests/fixtures/ptmqtl_pairs_bounds_plain.json");
    const BOUNDS: &str = include_str!("../tests/fixtures/ptmqtl_pairs_bounds.json");
    const DEPTH: &str = include_str!("../tests/fixtures/ptmqtl_pairs_depth.json");
    const DEPTH_COVARIATE: &str =
        include_str!("../tests/fixtures/ptmqtl_pairs_depth_covariate.json");
    const DEPTH_BOTH: &str =
        include_str!("../tests/fixtures/ptmqtl_pairs_depth_covariate_run_level.json");

    #[test]
    fn occupancy_sends_only_what_was_given() {
        let args = occupancy_request(Path::new("o.tsv"), &OccupancyOptions::default()).unwrap();
        assert_eq!(args, ["ptmqtl", "occupancy", "--observations", "o.tsv"]);
        let args = occupancy_request(
            Path::new("o.tsv"),
            &OccupancyOptions {
                samples: Some("s.tsv".into()),
                categories: Some(vec!["Common Biological".into(), "UniProt".into()]),
                out: Some("occ.tsv".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            &args[4..],
            [
                "--samples",
                "s.tsv",
                "--categories",
                "Common Biological\tUniProt",
                "--out",
                "occ.tsv"
            ]
        );
    }

    #[test]
    fn occupancy_refuses_bad_filters_before_anything_is_spawned() {
        for options in [
            OccupancyOptions {
                modifications: Some(vec!["a".into()]),
                categories: Some(vec!["b".into()]),
                ..Default::default()
            },
            OccupancyOptions {
                modifications: Some(vec![]),
                ..Default::default()
            },
            OccupancyOptions {
                categories: Some(vec!["a\tb".into()]),
                ..Default::default()
            },
            OccupancyOptions {
                out: Some(" ".into()),
                ..Default::default()
            },
        ] {
            let error = occupancy_request(Path::new("o.tsv"), &options).unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{options:?}");
        }
    }

    #[test]
    fn occupancy_parses_the_cases_and_keeps_every_censored_cell_without_a_fraction() {
        let occ: SiteOccupancy = serde_json::from_str(CASES).unwrap();
        assert_eq!((occ.site_count, occ.row_count), (9, 14));
        assert_eq!(occ.sample_count, None);
        assert!(occ.written.is_none() && occ.modifications.is_none());
        let rows = occ.rows().unwrap();
        assert_eq!(rows.len(), 14);
        for row in &rows {
            assert!(OCCUPANCY_STATES.contains(&row.state.as_str()));
            assert_eq!(row.fraction.is_some(), row.state == "quantified");
        }
        let fractions: SiteOccupancy = serde_json::from_str(FRACTIONS).unwrap();
        assert_eq!(fractions.sample_count, Some(2));
        assert!(fractions.samples_file.is_some());
    }

    #[test]
    fn occupancy_carries_mzlibs_peptidoform_counts() {
        // mzLib #1468: each cell counts the peptidoforms covering the site, and those modified.
        let occ: SiteOccupancy = serde_json::from_str(CASES).unwrap();
        let rows = occ.rows().unwrap();
        let ratio = rows
            .iter()
            .find(|r| r.protein_accession == "RATIO")
            .unwrap();
        assert_eq!(
            (ratio.covering_peptidoforms, ratio.modified_peptidoforms),
            (Some(3), Some(1))
        );
        let undetected = rows.iter().find(|r| r.state == "not_detected").unwrap();
        assert_eq!(undetected.modified_peptidoforms, Some(0));
        // A count_only cell seen only modified: k = n.
        for row in rows.iter().filter(|r| r.state == "count_only") {
            assert_eq!(row.modified_peptidoforms, row.covering_peptidoforms);
        }
        assert!(occ
            .caveats
            .iter()
            .any(|c| c.starts_with("covering_peptidoforms and modified_peptidoforms")));

        // A table without the columns (a bridge before pyMzLib 0.7.0) reads them as not supplied.
        let mut wire: serde_json::Value = serde_json::from_str(CASES).unwrap();
        let columns = wire["columns"].as_object_mut().unwrap();
        columns.remove("covering_peptidoforms");
        columns.remove("modified_peptidoforms");
        let names = wire["column_names"].as_array_mut().unwrap();
        names.retain(|n| !n.as_str().unwrap().ends_with("_peptidoforms"));
        let old: SiteOccupancy = serde_json::from_value(wire).unwrap();
        assert!(old
            .rows()
            .unwrap()
            .iter()
            .all(|r| r.covering_peptidoforms.is_none() && r.modified_peptidoforms.is_none()));
    }

    #[test]
    fn occupancy_written_to_out_has_no_rows() {
        let occ: SiteOccupancy = serde_json::from_str(WRITTEN).unwrap();
        assert!(occ.columns.is_empty());
        assert!(occ.rows().unwrap().is_empty());
        assert_eq!(occ.written.unwrap().row_count, 24);
        assert_eq!(occ.columns.names().len(), 0);
    }

    #[test]
    fn pairs_sends_only_what_was_given() {
        let args =
            pairs_request(Path::new("o.tsv"), PairType::Both, &PairsOptions::default()).unwrap();
        assert_eq!(
            args,
            [
                "ptmqtl",
                "pairs",
                "--observations",
                "o.tsv",
                "--type",
                "both"
            ]
        );
        let args = pairs_request(
            Path::new("o.tsv"),
            PairType::Covarying,
            &PairsOptions {
                modifications: Some(vec!["Common Biological:Phosphorylation on S".into()]),
                min_quantified_fraction: Some(0.5),
                keep_ceiling: true,
                remove_run_level: true,
                run_level_group: Some("category".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            &args[6..],
            [
                "--modifications",
                "Common Biological:Phosphorylation on S",
                "--min-quantified-fraction",
                "0.5",
                "--keep-ceiling",
                "--remove-run-level",
                "--run-level-group",
                "category"
            ]
        );
        let args = pairs_request(
            Path::new("o.tsv"),
            PairType::Covarying,
            &PairsOptions {
                remove_run_level: true,
                run_level_groups: Some("groups.tsv".into()),
                samples: Some("s.tsv".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            &args[6..],
            [
                "--samples",
                "s.tsv",
                "--remove-run-level",
                "--run-level-groups",
                "groups.tsv"
            ]
        );
        let args = pairs_request(
            Path::new("depth.tsv"),
            PairType::Covarying,
            &PairsOptions {
                include_bounds: true,
                run_covariate: Some("psms.tsv".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "ptmqtl",
                "pairs",
                "--observations",
                "depth.tsv",
                "--type",
                "covarying",
                "--include-bounds",
                "--run-covariate",
                "psms.tsv"
            ]
        );
        // Left out, neither is sent: mzLib's defaults (no bounds, no covariate) apply.
        let plain = pairs_request(
            Path::new("o.tsv"),
            PairType::Covarying,
            &PairsOptions::default(),
        )
        .unwrap();
        assert!(!plain
            .iter()
            .any(|a| a == "--include-bounds" || a == "--run-covariate"));
    }

    #[test]
    fn pairs_refuses_bad_arguments_before_anything_is_spawned() {
        for options in [
            PairsOptions {
                run_level_group: Some("file".into()),
                ..Default::default()
            },
            PairsOptions {
                min_quantified_fraction: Some(f64::NAN),
                ..Default::default()
            },
            PairsOptions {
                modifications: Some(vec!["a".into()]),
                categories: Some(vec!["b".into()]),
                ..Default::default()
            },
            PairsOptions {
                run_covariate: Some(" ".into()),
                ..Default::default()
            },
        ] {
            let error =
                pairs_request(Path::new("o.tsv"), PairType::Covarying, &options).unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{options:?}");
        }
    }

    #[test]
    fn pairs_with_bounds_tests_what_the_guard_kept_out() {
        // mzLib's BoundsEnterTypeAAsTiedRanks: six values in ten runs is below the 70% guard.
        let plain: PtmPairs = serde_json::from_str(BOUNDS_PLAIN).unwrap();
        assert_eq!(
            (
                plain.row_count,
                plain.include_bounds,
                plain.min_untied_share
            ),
            (0, Some(false), None)
        );
        let bounded: PtmPairs = serde_json::from_str(BOUNDS).unwrap();
        assert_eq!(
            (bounded.include_bounds, bounded.min_untied_share),
            (Some(true), Some(0.5))
        );
        assert_eq!(bounded.exclude_ceiling, None);
        let rows = bounded.rows().unwrap();
        assert_eq!(rows.len(), 1);
        let pair = &rows[0];
        assert_eq!(
            (pair.n, pair.test_status.as_deref(), pair.bounds_included),
            (9, Some("tested"), true)
        );
        assert!(pair.statistic.is_some() && !pair.run_covariate_removed);
        assert!(bounded
            .caveats
            .iter()
            .any(|c| c.starts_with("With --include-bounds, a pair is tested only if")));
    }

    #[test]
    fn pairs_given_a_run_covariate() {
        // mzLib's ARunCovariateIsPartialledOut: two sites that rise with each run's depth and a batch.
        let pair = |result: &PtmPairs| -> PtmPair {
            result
                .rows()
                .unwrap()
                .into_iter()
                .find(|r| r.site_a.starts_with("DEPTH_S") && r.site_b.starts_with("DEPTH_T"))
                .unwrap()
        };
        let plain: PtmPairs = serde_json::from_str(DEPTH).unwrap();
        assert_eq!(
            (
                plain.run_covariate_file.as_deref(),
                plain.run_covariate_name.as_deref()
            ),
            (None, None)
        );
        assert!(!pair(&plain).run_covariate_removed);

        let given: PtmPairs = serde_json::from_str(DEPTH_COVARIATE).unwrap();
        assert!(given
            .run_covariate_file
            .as_deref()
            .is_some_and(|f| f.ends_with("depth_psms.tsv")));
        assert_eq!(given.run_covariate_name.as_deref(), Some("target_psms"));
        assert_eq!(given.run_covariate_runs_missing, Some(Vec::new()));
        assert_eq!(
            (given.min_runs_for_run_level, given.min_sites_for_run_level),
            (Some(10), None)
        );
        assert!(pair(&given).run_covariate_removed && !pair(&given).run_level_removed);
        assert_ne!(pair(&given).statistic, pair(&plain).statistic);

        let both: PtmPairs = serde_json::from_str(DEPTH_BOTH).unwrap();
        assert!(pair(&both).run_covariate_removed && pair(&both).run_level_removed);
        assert!(pair(&both).statistic.unwrap().abs() < 0.1);
    }

    #[test]
    fn pairs_parses_the_physical_recording_and_mzlibs_defaults() {
        let physical: PtmPairs = serde_json::from_str(PHYSICAL).unwrap();
        assert_eq!(physical.r#type, "physical");
        assert_eq!(
            (physical.min_quantified_fraction, physical.exclude_ceiling),
            (None, None)
        );
        let pair = &physical.rows().unwrap()[0];
        assert_eq!(
            (pair.statistic, pair.p_value, pair.q),
            (Some(0.25), None, None)
        );

        assert_eq!(physical.include_bounds, None);
        assert_eq!(
            (pair.test_status.as_deref(), pair.bounds_included),
            (None, false)
        );

        let plain: PtmPairs = serde_json::from_str(BATCH).unwrap();
        assert_eq!(plain.min_quantified_fraction, Some(0.7));
        // mzLib #1468: mzLib resolves its own nullable excludeCeiling; nothing is restated here.
        assert_eq!(plain.exclude_ceiling, None);
        assert_eq!(
            (plain.include_bounds, plain.min_untied_share),
            (Some(false), None)
        );
        assert_eq!(
            (
                plain.run_covariate_file.as_deref(),
                plain.run_covariate_runs_missing.as_deref()
            ),
            (None, None)
        );
        assert!(plain
            .rows()
            .unwrap()
            .iter()
            .all(|r| r.test_status.as_deref() == Some("tested") && !r.run_covariate_removed));
        assert_eq!(
            (plain.remove_run_level, plain.run_level_group.as_deref()),
            (Some(false), None)
        );
        assert_eq!(plain.covarying_count, Some(15));

        let given: PtmPairs = serde_json::from_str(RUN_LEVEL).unwrap();
        assert_eq!(given.run_level_group.as_deref(), Some("none"));
        assert!(given.rows().unwrap().iter().all(|r| r.run_level_removed));
        // Every untested pair says why; spearman_method still reads not_estimable for each.
        for pair in given.rows().unwrap() {
            let status = pair.test_status.as_deref().unwrap();
            assert!(PAIR_TEST_STATUSES.contains(&status));
            assert_eq!(
                status == "tested",
                pair.spearman_method.as_deref() != Some("not_estimable")
            );
        }
        // An overlapping pair is reported but never tested.
        for pair in given.rows().unwrap() {
            if pair.overlapping {
                assert!(pair.q.is_none());
            }
        }
    }

    #[test]
    fn trait_effects_sends_only_what_was_given() {
        let args = trait_effects_request(
            Path::new("occ.tsv"),
            Path::new("runs.tsv"),
            &TraitEffectsOptions::default(),
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "ptmqtl",
                "trait-effects",
                "--occupancy",
                "occ.tsv",
                "--runs",
                "runs.tsv"
            ]
        );
        let args = trait_effects_request(
            Path::new("occ.tsv"),
            Path::new("runs.tsv"),
            &TraitEffectsOptions {
                min_runs: Some(4),
                min_replicates: Some(3),
                logit_epsilon: Some(0.001),
                keep_ceiling: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            &args[6..],
            [
                "--min-runs",
                "4",
                "--min-replicates",
                "3",
                "--logit-epsilon",
                "0.001",
                "--keep-ceiling"
            ]
        );
        let args = trait_effects_request(
            Path::new("occ.tsv"),
            Path::new("runs.tsv"),
            &TraitEffectsOptions {
                scale: Some(OccupancyScale::Log2Odds),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(&args[6..], ["--scale", "log2_odds"]);
        assert_eq!(
            [OccupancyScale::NaturalLogOdds, OccupancyScale::Log2Odds].map(OccupancyScale::as_str),
            OCCUPANCY_SCALES
        );
        assert!(trait_effects_request(
            Path::new("occ.tsv"),
            Path::new(""),
            &TraitEffectsOptions::default()
        )
        .is_err());
    }

    #[test]
    fn trait_effects_parses_the_recording_and_below_support_has_null_statistics() {
        let fit: TraitEffects = serde_json::from_str(EFFECTS).unwrap();
        assert_eq!((fit.min_runs, fit.min_replicates), (6, 4));
        assert!(fit.exclude_ceiling);
        assert_eq!(fit.scale, "natural_log_odds"); // mzLib's default
        let rows = fit.rows().unwrap();
        assert!(rows
            .iter()
            .all(|r| TRAIT_STATUSES.contains(&r.status.as_str())));
        assert_eq!(rows[1].covariate_effects.len(), 1);

        let strict: TraitEffects = serde_json::from_str(STRICT).unwrap();
        assert_eq!(strict.min_replicates, 6);
        let below = &strict.rows().unwrap()[0];
        assert_eq!(below.status, "below_support");
        assert_eq!(
            (below.effect, below.df, below.p_value, below.q),
            (None, None, None, None)
        );
    }

    #[test]
    fn trait_effects_on_log2_odds_is_the_natural_fit_over_ln2() {
        // mzLib's TraitFitOnLog2OddsIsTheNaturalLogFitDividedByLn2: a change of unit, not a refit.
        let natural: TraitEffects = serde_json::from_str(EFFECTS).unwrap();
        let log2: TraitEffects = serde_json::from_str(LOG2).unwrap();
        assert_eq!(log2.scale, "log2_odds");
        assert!(log2
            .caveats
            .iter()
            .any(|c| c.starts_with("scale is log2_odds")));
        let close = |a: Option<f64>, b: Option<f64>| {
            let (a, b) = (a.unwrap(), b.unwrap());
            assert!((a - b).abs() <= 1e-12 * b.abs(), "{a} {b}");
        };
        let ln2 = std::f64::consts::LN_2;
        let natural_rows = natural.rows().unwrap();
        let log2_rows = log2.rows().unwrap();
        assert_eq!(natural_rows.len(), log2_rows.len());
        for (n, r) in natural_rows.iter().zip(&log2_rows) {
            assert_eq!(n.site_key, r.site_key);
            close(r.effect, n.effect.map(|e| e / ln2));
            close(r.standard_error, n.standard_error.map(|e| e / ln2));
            close(
                r.replicate_variance,
                n.replicate_variance.map(|v| v / (ln2 * ln2)),
            );
            assert_eq!((r.df, r.p_value, r.q), (n.df, n.p_value, n.q));
        }
    }
}
