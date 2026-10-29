# Changelog

All notable changes to mzLibRust are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the crate follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

**No version has been released yet.** `Cargo.toml` says 0.1.0 and no tag exists, so everything
below is unreleased, and the crate as it stood on `main` before this file began is the baseline:
the `pride`, `peptidoform`, `flashlfq`, `readers` and `sdrf` modules and the bridge installer.

A wire verb's first appearance here is the version its spec in the bridge records as
`since.mzlibrust`; that is how the three bindings' reference pages agree on when a verb arrived.

## [Unreleased]

### Added

- **`ptmqtl` takes pyMzLib 0.7.0's three new options** (mzLib #1464, #1468), each sent only
  when set, so mzLib's default applies otherwise; the echoes the repin projected now report them:
  - `PairsOptions::include_bounds` (wire `--include-bounds`): co-varying pairs count runs where a
    site was seen in only one form as tied bounds, never invented values. A pair whose two sites
    are mostly one bound is `both_sites_mostly_bounds`, untested. `false` is mzLib's default.
  - `PairsOptions::run_covariate` (wire `--run-covariate`): a per-run covariate table, such as
    each run's target PSMs; pairs become the partial Spearman correlation given it.
  - `TraitEffectsOptions::scale` (wire `--scale`), an `OccupancyScale` (`OCCUPANCY_SCALES`):
    `Log2Odds` reports the same fit in log2 odds, a change of unit, not a refit. `None` is
    mzLib's natural-log odds.

  The ptmqtl guide gains "Counting one-form runs as bounds", "Given each run's depth" and the log2
  example, on mzLib's `BoundsEnterTypeAAsTiedRanks`, `ARunCovariateIsPartialledOut` and
  `TraitDesign` cases (`bounds.tsv`, `depth.tsv` and `depth_psms.tsv` vendored byte for byte from
  pyMzLib 0.7.0), and `PENDING` in `tests/spec_docs.rs` is empty.

- **The bridge pin is pyMzLib v0.7.0** (mzLib 1.0.596, tag commit `cc67b7a2`; bridge-watch's
  repin, mzLibRust #37). `install_bridge()` fetches it, and its bridge dispatches 57 verbs. Every
  answer it changes for the same call is in `docs/upgrading.md`'s new pyMzLib 0.7.0 section:
  `sdrf::design_with` keeps the SDRF's replicate numbers (mzLib #1426), `diff::design_with`
  numbers contrast ids once over the analysis (#1476), `ptmqtl` reports more `not_detected` cells
  and leaves `exclude_ceiling` to mzLib (#1468), `peptidoform::convert_with` fails malformed
  modification text by row (#1457) and resolves some terminal modifications (#1442), a locked
  `.parquet` is `MzLibError::Bridge` (#1478), TMT10's `TMT131` reads (#1467), and a peptide with
  a permanently charged modification weighs one proton less (pyMzLib #99).
- **New output fields of verbs the crate already has**, projected from the 0.7.0 bridge:
  - `isobaric::ReporterChannel::sdrf_label` and `IsobaricKit::sdrf_labels()`: the SDRF
    `comment[label]` of each channel, by mzLib's `SdrfIsobaricLabels` (#1467), `TMT131` for
    TMT10's 131N. `None` only from an older bridge.
  - `diff::NotRunContrast::id`: the default contrast's id a stratum skipped, the same in every
    stratum (#1476).
  - `ptmqtl::OccupancyRow::covering_peptidoforms` and `modified_peptidoforms`: the k of n of a
    cell, counted in identifications (#1468); `None` from an older bridge.
  - `ptmqtl::PtmPair::test_status` (`PAIR_TEST_STATUSES`: why a co-varying pair was not tested),
    `bounds_included` and `run_covariate_removed`; and on `PtmPairs`, `include_bounds`,
    `min_untied_share`, `run_covariate_file`, `run_covariate_name` and
    `run_covariate_runs_missing` (#1468). Each echo is mzLib's default, or `None`, until the
    crate sends the options.
  - `ptmqtl::TraitEffects::scale`: the scale of every effect and variance, `natural_log_odds` by
    mzLib's default (#1464).
- **Specs re-vendored at bridge `45833ca`** (pyMzLib 0.7.0's specs, `4c9eb13`, plus the
  `stats fit` Rust example that compiles), every spec but `modifications resolve` as before. The
  four new verbs' specs come in with them, so the front page lists `genes map-residues`,
  `library calibrate-rt`, `stats fit-latent-mixed` and `stats fit-link` as `(not yet)`, and the
  three new options of verbs the crate has (`ptmqtl pairs --include-bounds` and
  `--run-covariate`, `ptmqtl trait-effects --scale`) are `PENDING` in `tests/spec_docs.rs`. Every
  recording a spec names is vendored byte for byte from pyMzLib 0.7.0, those of the not-yet
  verbs and options included, because the replay bridge reads each spec's examples; the replay
  bridge takes pyMzLib 0.7.0's new `RECORDED_WITH` and `OPTIONAL_FILES` entries.

- **`stats` gains five verbs, ported from pyMzLib 0.6.0** (mzLib 1.0.595; mzLib #1367, #1385,
  #1446), each with its guide section, an executed example and live canaries against R:
  - **`stats::combine_p`**: one p-value per feature from several independent studies, by
    Fisher's or Stouffer's method (`COMBINATION_METHODS`; no default, you choose), optionally
    weighted and optionally signed by each test's effect (`PValueTest::weight`,
    `PValueTest::effect`), so opposite effects cancel instead of reinforcing. Held to R's
    `pchisq` and `qnorm` to 1e-8. Returns `CombinedPValues`, with `p_values()`. Wire verb
    `stats combine-p`.
  - **`stats::correlate` / `correlate_with`**: Spearman's rho of every feature with a per-sample
    trait, optionally partial given covariates (`CorrelateOptions::given`), as R's `cor.test` and
    ppcor's `pcor.test`. Returns `Correlations`, with typed `rows()` (`Correlation`); the
    argument is `trait_table` because `trait` is a Rust keyword. Wire verb `stats correlate`.
  - **`stats::fit_logistic` / `fit_logistic_with`**: one logistic regression per feature of a 0/1
    table, such as a detection table, as R's `glm(family = binomial)`, with a Wald z.
    `FitLogisticOptions::max_iterations` and `tolerance` are `None` (mzLib's) unless you set
    them. A separated feature is reported, with `None` statistics, never dropped.
    `LogisticModelFit::rows` gives `LogisticRow`s; `LOGISTIC_STATUSES`. Wire verb
    `stats fit-logistic`.
  - **`stats::fit_mixed` / `fit_mixed_with`**: one linear mixed model per feature with a random
    intercept per group of samples, as `nlme::lme`, tested on nlme's containment df.
    `FitMixedOptions::estimator` is `None` (mzLib's REML) unless you set it; `MIXED_ESTIMATORS`,
    `MIXED_STATUSES`. `MixedModelFit::rows` gives `MixedRow`s. Wire verb `stats fit-mixed`.
  - **`stats::entrapment_fdp` / `entrapment_fdp_with`**: a search's false discovery proportion
    from its entrapment hits, as the lower bound, the combined upper bound and the paired upper
    bound of Wen et al. 2025, at each q-value threshold. `ratio` is required (mzLib has no
    default); `EntrapmentFdpOptions::thresholds` is `None` for mzLib's default sweep.
    `EntrapmentFdpEstimates::rows` and `at(threshold)` give `FdpRow`s. Wire verb
    `stats entrapment-fdp`.
- **`stats::fit_with` tests contrasts, with confidence intervals, and takes mzLib's prior
  estimator** (mzLib #1452, #1368). `FitOptions::contrasts` names a contrasts file — the one
  `diff::design_with` writes beside the design, or your own — and each contrast is moderated with
  the same prior as the coefficients, has its own Benjamini-Hochberg family, and carries
  `confidence_low` and `confidence_high`: limma's `contrasts.fit`, `eBayes` and
  `topTable(confint = ...)`. The answers are a second table, `ModeratedFit::contrasts` (a
  `ModeratedContrasts`, with `rows()` giving `ContrastRow`s and `weights` per contrast), so the
  coefficient table is unchanged; `coefficients` may now be empty when contrasts are given, and
  `ModeratedFit::contrasts_file` echoes the file. `FitOptions::confidence` sets the level, as a
  fraction; `None` is mzLib's 0.95. `FitOptions::estimator` takes `"moments_legacy"` or
  `"marginal_likelihood"` (`PRIOR_ESTIMATORS`); `None` is mzLib's default, and the new
  `VariancePrior::estimator` reports which was used. On the MALAT1 dilution series the
  marginal-likelihood prior lowers the prior df from 19.6 to 11.0 and drops the false FLuc call
  at half load. A bridge before pyMzLib 0.6.0 lists `stats fit` but ignores the three options;
  `fit_with` then returns `MzLibError::Usage` naming the 0.6.0 bridge instead of a fit that
  ignored them.
- **A new module, `diff`: the design and the observations, set up by mzLib** (mzLib #1450,
  #1455), with a guide that runs:
  - **`diff::design_with`**: a sample sheet and a model (`FactorSpec`, `CovariateSpec` with its
    unit, interactions, strata, batch, in `DesignOptions`) into mzLib's design per stratum — the
    matrix, equal to R's `model.matrix`, its rank and aliased terms, and the default contrasts
    with their estimability (`ESTIMABILITIES`), as `Design`, `DesignStratum` and
    `Design::contrasts()` (`DefaultContrast`). `DesignOptions::out` writes `design_<i>.tsv` and
    `contrasts_<i>.tsv`, which `stats::fit_with` reads as they are. Wire verb `diff design`.
  - **`diff::observations` / `observations_with`**: every peptide in every biological sample of a
    FlashLFQ or MetaMorpheus peptide table, as a log2 intensity or `None` with the reason
    (`OBSERVATION_STATES`), fractions summed and technical replicates averaged by mzLib's rules
    from the `ObservationRun`s you give, with the shared-peptide normalization shifts
    (`NormalizationResult`). `QuantBasis` decides whether match-between-runs values enter; there
    is no default. `ObservationsOptions::responses_out` writes the peptide-by-sample table
    `stats::fit_with` reads. Wire verb `diff observations`.
- **A new module, `ptmqtl`: PTM site occupancy and PTM-QTL statistics, computed by mzLib's
  `Quantification.PtmQtl`** (mzLib #1369, #1430, #1446), with a guide that runs on mzLib's own
  hand-computed cases:
  - **`ptmqtl::occupancy` / `occupancy_with`**: every site's intensity occupancy in every run (or
    sample, with `OccupancyOptions::samples`), each cell with its censoring state
    (`OCCUPANCY_STATES`): a floor or an undetected cell has no fraction, never 0.
    `OccupancyOptions::out` writes the table `trait_effects` reads. `SiteOccupancy::rows` gives
    `OccupancyRow`s. Wire verb `ptmqtl occupancy`.
  - **`ptmqtl::pairs` / `pairs_with`**: PTM-PTM pairs, physical (both on one peptidoform) or
    co-varying (Spearman across runs, with a Benjamini-Hochberg `q` per family), chosen with
    `PairType` (no default). `PairsOptions::remove_run_level` correlates given each run's level,
    for a batch shift. Every co-varying threshold is `None` (mzLib's) unless you set it, and the
    result reports the value used. `PtmPairs::rows` gives `PtmPair`s. Wire verb `ptmqtl pairs`.
  - **`ptmqtl::trait_effects` / `trait_effects_with`**: per site, the effect of a per-run trait on
    logit occupancy, from a mixed model with a random intercept per biological replicate.
    `TraitEffectsOptions`' thresholds are `None` (mzLib's) unless you set them.
    `TraitEffects::rows` gives `TraitEffect`s; `TRAIT_STATUSES`. Wire verb
    `ptmqtl trait-effects`.
- **The replay bridge learns pyMzLib 0.6.0's matching rules for these verbs**: `NESTED_ECHOES`
  (`stats fit`'s estimator and confidence level), `OPTIONAL_FILES` (`--contrasts`, `--given`,
  `--samples`, `--protein-groups`, `--run-level-groups`), `OPTIONAL_OUTPUTS`
  (`--responses-out`), `GIVEN_ECHOES` (`--thresholds`), a numeric echo matching a call by value
  (`--ratio 1.0` fits a recorded `1`), the new `RECORDED_WITH` entries (flags included), and the
  stdin echoes of `stats fit`, `fit-logistic`, `fit-mixed`, `combine-p`, `diff design` and
  `diff observations`. The pyMzLib 0.6.0 recordings and inputs these examples replay are
  vendored byte for byte under `tests/fixtures/` (`stats/`, `diff/`, `ptmqtl/`).
- **The bridge pin is pyMzLib v0.6.0** (still mzLib 1.0.595, `d10017ce`; bridge-watch's repin).
  Its bridge dispatches 53 verbs. **Specs re-vendored in full at bridge `60b1d25`**: every spec
  but `modifications resolve`, which has no verb (`not_vendored` in `docs/specs/SOURCE`). With
  the stats and data ports, this crate builds every one of them, so the front page lists no verb
  as `(not yet)`. Every spec's
  recorded examples are vendored byte for byte from pyMzLib 0.6.0, because the replay bridge
  reads each spec's examples, and so are the three recordings that changed since 0.5.0
  (`bridge_version.json`, `stats_fit_malat.json`, `stats_fit_limma.json`: `stats fit`'s first
  caveat now names its estimator). The replay bridge carries pyMzLib's whole `RECORDED_WITH`
  table, flags included.
- **`pride::file_checksums`, `pride::file_checksums_with`: the MD5 and exact size PRIDE recorded for
  each file deposited in a project** (`pride checksums`, mzLib #1427; spec `pride.checksums.yaml`).
  It is what a download is checked against, and `download` and `download_files` check neither
  yet. A `BTreeMap<String, PrideFileChecksum>` keyed by file name, in the wire's file-name order,
  as pyMzLib's dict is; each row is read strictly, so a row missing its MD5 or size is
  `MzLibError::Protocol`, not a blank that blames the download later. An unknown or unpublished
  accession is `MzLibError::ProjectNotFound`; an existing project with no checksums is an empty
  map; a list mzLib cannot read stays `MzLibError::Bridge` (`MzLibException`). No time limit of
  the crate's own, as for every PRIDE call. The `pride` guide gains *Checking a download against
  PRIDE's checksums*: the check itself (size, then MD5, with the `md-5` crate, now a
  dev-dependency), what the list covers (not `generated/`, files the manifest omits, non-deposits)
  and its two traps. The tests' `spec_docs` lint reads a map return as a list of its values.
- **`sdrf::isobaric_design_with`: the isobaric (TMT, TMTpro, iTRAQ, DiLeu) design an SDRF
  describes, or every reason it cannot** (`sdrf isobaric-design`, mzLib #1380; spec
  `sdrf.isobaric-design.yaml`). `SdrfIsobaricDesign` holds one row per file and annotated channel,
  `channels()` gives them as `DesignedChannel`, and **every number is 1-based**, as the SDRF and
  `TmtDesign.txt` write them, where `design_with`'s are 0-based: the split is mzLib's, and neither
  is converted. Nothing is guessed: `IsobaricDesignOptions::kit` and exactly one of `plex_column`,
  `plex_file_name_pattern` and `single_plex` must be given, and a missing one is mzLib's refusal,
  a result, not an error. `out` names a folder, where mzLib writes MetaMorpheus's `TmtDesign.txt`.
- **`sdrf::restrict_with`, `sdrf::restrict`: an SDRF cut to the files one search read** (`sdrf
  restrict`, mzLib #1451; spec `sdrf.restrict.yaml`). mzLib joins `X-calib.mzML` to the row naming
  `X.raw` by stem; a `SearchedFile::acquired_from` overrides the rule for one file.
  `RestrictedSdrf` is the restricted `SdrfDocument` with `comment[searched data file]`
  (`SEARCHED_DATA_FILE_COLUMN`) after `comment[data file]`, and mzLib's four reports:
  `searched_without_row`, `dropped_data_files`, `ambiguous`, `shared_stem`. Restrict first, then
  read a design: both design readers key on the searched name. The `sdrf` guide gains *Read an
  isobaric (TMT) design* and *Restrict an SDRF to one search's files*, every example replayed from
  pyMzLib's recordings over mzLib's own PXD008841 and PXD061609 subsets (copied byte for byte), and
  the SDRF-B5 trap: mzLib's design readers split a path on the running platform's separator only,
  so the live tests build native paths. The replay bridge holds a `sdrf restrict` call to the
  recording made from the same searched files, as pyMzLib's does.
- **`proteins::align_with`, `proteins::align`: mzLib's global affine-gap aligner, and every
  residue's counterpart** (`proteins align`, mzLib #1428; spec `proteins.align.yaml`). BLOSUM62,
  Gotoh's Needleman-Wunsch; `Alignments` holds the pairs table and the residue map, and
  `alignments()` gives one `Alignment` per pair with `map_source_to_target` and
  `map_target_to_source` (1-based; `None` for a residue facing a gap, never a neighbour).
  `AlignOptions`' `gap_open`, `gap_extend`, `end_gaps` (`EndGaps`) and `max_cells` are `None` by
  default and sent only when set, so mzLib's defaults (11, 1, free, 200,000,000 cells) apply and
  the result echoes them; `residue_map: false` sends `--no-residue-map`.
- **`proteins::parse_accessions`: accessions read back into entry, isoform or version, and
  applied variants, with mzLib's grammars** (`proteins parse-accession`, mzLib #1382; spec
  `proteins.parse-accession.yaml`). Each accession is sent exactly as given, never trimmed; what
  fits neither UniProt's nor RefSeq's protein grammar is `unrecognized`, never repaired.
  `ParsedAccessions::entry_of` groups isoforms and variants under their entry; `NAMESPACES` is
  `resolve_genes_with`'s vocabulary.
- **`proteins::build_orthology_with`: an Ensembl Compara orthology snapshot, and optionally its
  Parquet** (`genes orthology`, mzLib #1381; spec `genes.orthology.yaml`). `OrthologyOptions`
  takes the species' GTFs (or compact gene sets), **every** species' homology dump, the gene
  trees, and `pair_status` for every gene's status in one other species (`PairStatus`,
  `ORTHOLOGY_STATUSES`). `OrthologySnapshot` says what went in (`GeneSet`, `GeneTrees`,
  `HomologyDump`), what each pair holds (`SpeciesPair`), what was written (`WrittenSnapshot`,
  `SnapshotFile`) and the `snapshot_id`. mzLib's three-species test world is vendored byte for
  byte under `tests/fixtures/orthology/`, marked `-text` so its sha256s, and the snapshot id the
  live test holds, survive a Windows checkout. The `proteins` guide gains *Map a residue from one
  sequence onto another*, *Read an accession into entry, isoform and variants* and *Orthologs
  across species*, each replayed from pyMzLib's recordings; the DuckDB queries over the Parquet
  are shown as SQL, since this crate carries no Parquet reader. The replay bridge holds each of
  the three to the stdin it was recorded from, `proteins read` to its `--sequences` choice, and
  `genes orthology` to its `--pair-status` direction and `--restrict-gene-trees` choice, as
  pyMzLib's does; it reads one level into an echoed object (`pair_status.from_species`).
- **A new module, `library`: mzLib's binary spectral libraries** (`library read`; spec
  `library.read.yaml`). `library::read`, `library::read_with` and `library::read_many` read a
  `.msl` with mzLib's `MslLibrary` into one row per precursor, in stored order (`precursor_index`,
  not m/z), or with `LibraryReadOptions::fragments` one row per fragment ion; mzLib's "not
  available" values (ion mobility 0.0, a NaN q-value, NCE 0, an empty protein) are `None`, and
  `retention_time_unit` is per row, iRT or minutes. `index_only_threshold` is `None` by default and
  sent only when set, so mzLib's own 1 GiB rule applies; `limit` and `offset` count precursors.
  `LibraryRecords`, `LibraryBatch` and `LibraryFileReport` carry the counts, the load mode and the
  format version; a skipped library's unknown facts are `None`, never 0. The module guide covers
  every column, the iRT trap, targets and entrapment, fragment ions, index-only loading and batches,
  every example replayed from pyMzLib's recordings of the two libraries vendored under
  `tests/fixtures/library/` (marked binary in `.gitattributes`). `tests/live_library.rs` holds the
  recordings to the live bridge.
- **The bridge pin is pyMzLib v0.5.0** (mzLib 1.0.595, tag commit `d10017ce`, which carries 1.0.594
  as well; mzLibRust #33). What it changes for a caller is under *Changed* below, and
  `docs/upgrading.md` lists every answer that differs from the same call.
- **DIA-NN 2.x `report.parquet` reads** (mzLib #1384, file type `DiaNnReportParquet`), through the
  same `read_results`, `read_records` and `flashlfq` calls as the 1.8/1.9 `report.tsv`, with
  retention time in minutes. **`readers::formats()` lists 39 types, up from 38, and 5 offer the
  `quantifiable` view, up from 4**; 19 still offer none. It is recognised by its columns
  (`Precursor.Id`, `Stripped.Sequence`, `Run`), not its name, and any other `.parquet` is
  `MzLibError::Usage`.
- **Specs re-vendored at bridge `98ee566`** (pyMzLib #82 to #85): each fragments digestion option,
  the trend's spline size and both PRIDE page sizes have no bridge default and name the mzLib
  symbol that owns it; the ProForma caveat is retired; `pride ftp-files` has a recorded example,
  so the replay bridge reads it from the spec; and the bridge's re-check of its 93 open questions
  (37 answered and folded into the docs, 25 reworded) is carried in. The reference pages are
  re-rendered. The bridge's `modifications resolve`
  spec is left out (`not_vendored` in `docs/specs/SOURCE`): no binding has built that verb, and a
  vendored spec would list it on the front page. `scripts/sync-specs.ps1 -OnlyVendored` does
  that, as pyMzLib's `sync_specs.py --only-vendored` does.
- The seven pyMzLib recordings that changed since 0.4.0 are re-vendored byte for byte: both
  albumin digests (`peptidoform_P02768.json`, `peptidoform_P02768_small.json`, now trypsin and
  HCD), the three `peptidoform convert` recordings, `readers_formats.json` and
  `readers_results_psmtsv.json`. The replay bridge learns pyMzLib's `RECORDED_WITH`: a call that
  sends nothing cannot be told from one that sent `--max-mods 1` by the echo alone, so the trimmed
  albumin recording answers only a call that sends it.

- **`peptidoform::convert_with`, `peptidoform::convert`: rewrite full sequences in another notation
  with mzLib** (`peptidoform convert`, pyMzLib #75; spec `peptidoform.convert.yaml`). The main use
  is MetaMorpheus full sequences to Unimod accessions: `[UniProt:N-acetylserine on S]SEQK` becomes
  `[UNIMOD:1]SEQK`, and dimethyllysine becomes `UNIMOD:36`. `SequenceConversions` holds one row per
  input, in order, with mzLib's own verdict (`converted`, `converted_with_warnings` or `failed`),
  the modifications it could not write, and its warnings; `sequences()` and `not_converted()` give
  them as `ConvertedSequence` rows. `ConvertOptions::source` and `target` take any notation mzLib
  has registered, and every result lists them; `mode` is mzLib's `SequenceConversionHandlingMode`
  (`ConversionMode`). **The ProForma target does not resolve UniProt modifications** in this mzLib
  build and writes them back by name; convert to Unimod for those. The verb needs the bridge from
  pyMzLib 0.4.0, which this crate now pins; an older bridge is refused before anything is spawned,
  and the live tests skip on one. The replay bridge holds a `peptidoform convert` call to the recording made from the same
  sequences, as pyMzLib's does.
- **The bridge pin moves to pyMzLib v0.4.0** (mzLib 1.0.593, `0a808fec`): `install_bridge()`
  fetches it, verified against that release's `SHA256SUMS`. Its bridge dispatches
  `peptidoform convert`, and its `--out` accepts only `.tsv` (pyMzLib #76; `ReadOptions::out` and `BulkOptions::out` say so); `read_records` on an
  SDRF leaves the `header` and `cells` lists out (in `excluded_fields`).
- **Specs re-vendored at bridge `e761578`**: adds `peptidoform.convert.yaml`; records pyMzLib 0.3.0
  as the first release of four verbs, the isobaric citation note and mzLibR #29; and documents the
  0.4.0 bridge's `.tsv`-only `out` (PYB-1) and its too-large-answer usage error (PYB-3 part 1) on
  the reader verbs and `sdrf design`.

- **Every guide opens with a question table, runs on real data, and says what to cite** (pyMzLib
  #72). A module's `//!` documentation is its guide; `tests/docs_lint.rs` holds each one to a
  `| You want to | Call | mzLib type |` table in its first lines, a `## Cite` section rendered from
  the specs' DOIs (`docs/reference/cite.<module>.md`), examples that run or say `Not run:` and why,
  and no counts or mzLib versions in prose (README and the crate front page included). New worked
  examples replay pyMzLib's guide recordings: the K562 MBR run and its median polish, the full
  albumin digest, PXD000001's FTP tree and an unknown accession, the MSFragger identify, the RNA
  transcript groups, MS1 peaks, a filtered protein read and the SDRF assessments. The crate front
  page gains *Citing* and a generated *Which function do I want?* table; `MzLibError` lists every
  error each verb can return; `docs/upgrading.md` lists the changes that give a different answer
  from the same file; `CITATION.cff` gives the software citation.

- **Gene Ontology for protein groups** (pyMzLib #70; specs `proteins.annotate-go.yaml`,
  `proteins.update-go.yaml`). `proteins::annotate_go_with` (+ `annotate_go`, `GoAnnotateOptions`)
  annotates a stored MetaMorpheus protein-group table against a UniProt XML and a go.obo you keep,
  with mzLib's `GoGroupAnnotator`: one row per (group, GO term) that any member holds, directly or
  through an ancestor, in mzLib's `GoAnnotationTsv` columns, with every input's sha256, mzLib's
  header counters, an optional category map (`GoCategories`), and `out` / `categories_out` written
  by mzLib's own writers. `proteins::update_go` fetches the current release on purpose, keeping a
  different file as a timestamped backup (`GoUpdate`). `string_list_maps` reads the
  `evidence_by_member` column. Both need the bridge from pyMzLib 0.3.0 and are refused, before
  anything is spawned, by an older one.
- **`sdrf::design_with`, `sdrf::design`** (`sdrf design`, pyMzLib #69, mzLib #1363): the
  label-free experimental design read out of an SDRF by mzLib's `SdrfLabelFreeDesign`, as an
  `SdrfDesign` with `is_valid`, every refusal at once, and `notes` recording each renumbering.
  `DesignOptions { condition_columns, searched_files, out, timeout }`. `SdrfDesign::files` gives
  the 0-based `DesignedFile` rows; `spectra` and `run_design` hand them to `flashlfq::quantify_with`
  and `median_polish_with`, and refuse a refused design. `out` writes MetaMorpheus's 1-based
  `ExperimentalDesign.tsv` through mzLib's writer.
- **The `isobaric` module** (`isobaric kits`, pyMzLib #69): `isobaric::kits(None)` lists every
  kit mzLib can name and `kits(Some("TMT18"))` one, matched by mzLib's whole-name rule. Each
  channel's label, theoretical reporter-ion m/z and matching window as an `IsobaricKits` table,
  grouped by `by_kit` and `kit`.

- **`stats`, a new module: differential abundance with no R** (pyMzLib #71; specs
  `stats.fit.yaml`, `stats.adjust.yaml`, `stats.meta.yaml`). `fit`/`fit_with` is limma's `lmFit`
  then `eBayes(legacy = TRUE)` (`FitOptions { trend, spline_basis, threads, timeout }`, returning
  `ModeratedFit` with its `VariancePrior` and typed `FitRow`s per coefficient); `adjust` is
  Benjamini-Hochberg over `&[Option<f64>]`, `None` left out of the family; `meta`/`meta_with` pools
  `Study` effect sizes per feature with DerSimonian-Laird (`MetaOptions { confidence, timeout }`).
  Every number is mzLib's. Against limma's and metafor's own reference output, through this crate,
  the worst relative difference is below 1e-10. Needs the bridge from pyMzLib 0.3.0: an older one is
  refused before anything is read. The reference tables are pyMzLib's `tests/fixtures/stats/`, byte
  for byte.
- **The bridge pin moves to pyMzLib v0.3.0** (mzLib 1.0.593, `0a808fec`): `install_bridge()`
  fetches it, verified against that release's `SHA256SUMS`.
- **The mzLib 1.0.593 tables** (pyMzLib #68, mzLib #1388, #1365).
  - An RNA search's `AllQuantifiedTranscriptGroups.tsv` reads through `read_protein_groups` and
    `read_occupancy`, and its `AllQuantifiedOligos.tsv` through `read_quantified_peptides`. mzLib
    reads them with subclasses of the protein-group and peptide readers, so the columns keep their
    protein names (`protein_group_name` is the transcript group, `sequence` the oligonucleotide),
    and occupancy names RNA modifications. `formats()` lists two more types, neither with a view.
  - MetaMorpheus protein-group tables written without quantification (`AllProteinGroups.tsv`,
    `<file>_ProteinGroups.tsv`) now read through `read_protein_groups`, with `intensity` in
    `absent_fields`; they failed with `Tsv file type not supported`.
- **The specs at bridge `078d5f8`** in `docs/specs/`, with the recordings they name, byte for byte
  from pyMzLib: the new verbs' specs (`proteins annotate-go`/`update-go`, `sdrf design`,
  `isobaric kits`, `stats fit`/`adjust`/`meta`) arrive ahead of their functions, which the spec lint
  skips until they exist.
- **The replay bridge follows pyMzLib's rules again**: a recording with `written` answers only an
  `out` call (and the reverse); an input file echoed as `<name>_file`, or under another name
  (`--max-mods` as `max_modifications`, `--psms` as `psm_file`), must be the recording's file; a
  filtered `proteins read` recording answers only a filtered call; a single-kit `isobaric kits`
  recording answers only a call naming a kit; and `quant flashlfq` must send the recording's runs
  on stdin.

- **The mzLib 1.0.592 batch** (pyMzLib #63, #64, #65), each verb against its spec in
  `docs/specs/`. It needs the bridge from pyMzLib 0.2.0 (mzLib 1.0.592): on the pinned 0.1.1
  bridge the new verbs and the many-file form are not dispatched, and the new fields read as empty
  or `None`.
  - **Many files in one call.** Every reader has a `_many` twin — `identify_many`,
    `read_results_many`, `read_records_many`, `read_features_many`, `read_matches_many`,
    `read_spectra_many`, `read_protein_groups_many`, `read_quantified_peptides_many`,
    `read_occupancy_many` — reading a list in one bridge process into one long table
    (`ReadBatch`: `source_index` and `source_path` first, one `FileReport` per file), with
    `BulkOptions { threads, on_error, out, timeout }` and `OnError::{Fail, Skip}` (BULK.md).
    `ReadBatch::in_minutes` converts each file's rows by its own unit.
  - **`absent_fields`** on every reader: a column the view defines but this file's format has no
    source for, so every value is `None` rather than mzLib's default — beside `failed_fields` and
    `excluded_fields`, whose entries now name the verb that carries them (`ExcludedField::verb`).
  - **The run a spectra file came from**: `ScanRecords::source` (`SpectraSource`: instrument model
    and PSI-MS accession, serial number, acquisition start; mzLib #1349), with
    `SpectraSource::acquired_at` typing the time as an instant or a local clock reading.
  - **`read_matches`**: mzIdentML's `q_value`, `rank` and `pass_threshold` columns, the items mzLib
    skipped (`skipped_count`, `skipped`; #1313), `row_count`, and the engine's scores as long rows
    with `MatchOptions::scores` (#1306).
  - **`read_protein_groups`, `read_quantified_peptides`, `read_occupancy`** (#1347): MetaMorpheus
    protein groups, FlashLFQ peptides and PTM site occupancy as long tables, one row per record
    per sample.
  - **`sdrf::validate`, `validate_many`, `lint_labelled`, `lint`, `assess`, `assess_many`,
    `samples`, `samples_many`, `parse_ages`**: mzLib's structural findings (`SdrfValidator`),
    cross-document drift (`SdrfDriftLint`), Informative / Partial / Skeleton verdicts with their
    evidence (`SdrfSampleInformativeness`), per-sample values with conflicting columns withheld
    (`SdrfSampleBlock`), and ages read into years with the reason a cell was refused
    (`SdrfAge.TryParse`), with `sdrf::BulkOptions` for a corpus.
  - **A `proteins` module**: `read` / `read_with` (UniProt XML or FASTA to one row per protein,
    with GO-term and Ensembl-gene long tables on request and an exact accession filter),
    `resolve_genes` / `resolve_genes_with` (stable Ensembl gene ids against a GTF you pin, one
    outcome per protein, every input's sha256), and `classify_peptides` /
    `classify_peptides_with` (Unique, SharedWithinGene, SharedAcrossGenes or NotInDatabase, with
    I and L one residue).
  - **`BridgeVersion::verbs`** and `has_verb`: the verbs a bridge dispatches. The quantification
    readers check it first and refuse an older bridge with the release they need, instead of
    spawning a process to be told "Unknown command".
- **Parity debt closed**: `pride::search` / `search_with` + `SearchOptions` (find PRIDE projects by
  keyword, every page fetched, no accession repeated, dates as calendar dates), and
  `flashlfq::median_polish` / `median_polish_with` + `MedianPolishOptions` / `DesignEntry`
  (re-quantify proteins from a `QuantifiedPeptides.tsv` under a new design, returning
  `MedianPolishResults` with the samples that key its intensities). Both verbs were already in the
  pinned bridge and in pyMzLib.
- **Reference facts rendered from the bridge's per-verb specs.** Every function that calls a wire
  verb now documents, in rustdoc, the parameters with their units, defaults and ranges, the result
  fields with their units and what a null means, the error kinds and the `MzLibError` variant each
  becomes, the caveats, the mzLib code it wraps (linked at the pin), the same verb in pyMzLib and
  mzLibR, references and since. They are rendered from `docs/specs/` (vendored from the bridge by
  `scripts/sync-specs.ps1`) into `docs/reference/`, and included with `#[doc = include_str!(..)]`.
  See `docs/reference-facts.md`.
- **A doc lint** (`tests/spec_docs.rs`) that fails when a spec parameter or result field has no doc
  on the Rust options or result struct, or a doc that does not name the spec's unit, and when a
  rendered fragment is stale. Deliberate differences from a spec are declared, with a reason, in one
  table.
- **Examples that run.** The doc examples are doctests against a stand-in bridge
  (`tools/replay-bridge`) that answers from the recorded fixtures pyMzLib and mzLibR share, and
  only when a recording fits the call. 60 of the crate's 69 examples execute, up from 0 of 10; the
  nine left `no_run` download from EBI, or read many files of a kind no recording exists for yet,
  and say so.
- The recorded fixtures of pyMzLib's mzLib 1.0.592 batch, byte for byte, in `tests/fixtures/`.
- `[package.metadata.docs.rs]`, so docs.rs builds the documentation with every feature.
- This changelog.

### Changed

- **`PtmPairs::exclude_ceiling` is `None` when you did not choose** (the pyMzLib 0.7.0 bridge,
  mzLib #1468). mzLib made `CoVarying`'s `excludeCeiling` nullable and resolves it itself, and the
  bridge does not restate the rule to name the value; it was `Some(true)`. `Some(false)` still
  means `PairsOptions::keep_ceiling`, and ceilings are still left out by default.
- **`sdrf::design_with`'s `DesignedFile::biological_replicate` is the SDRF's own number** (the
  pyMzLib 0.7.0 bridge, mzLib #1426), gaps included, where it was ranked within each condition;
  `SdrfDesign::notes` flags the gaps. The guide's example and the offline test now hold the kept
  numbers (22, 23, 24 stay 21, 22, 23 0-based) instead of the ranked design.
- **`diff::DefaultContrast::contrast_id` is one comparison in every stratum** (mzLib #1476), and the
  `diff` and `stats` guides say so; they told you to match strata by label.
- **`sdrf::pool_with(.., PoolOptions { out, .. })` reports the file mzLib wrote** (the pyMzLib
  0.6.0 bridge). mzLib adds `.sdrf.tsv` to a name that does not already end in it, so an `out` of
  `x.tsv` writes `x.tsv.sdrf.tsv`; `WrittenSdrf::path` said `x.tsv`, a file that does not exist,
  and now says `x.tsv.sdrf.tsv`. The file written is the same. `PoolOptions::out` and
  `WrittenSdrf::path` say so, and `docs/upgrading.md` has the entry.
- **PRIDE's metadata calls add no timeout of their own** (pyMzLib #87). `pride::list_files` and
  `pride::search` defaulted, through `ListOptions::default()` and `SearchOptions::default()`, to a
  300 s timeout, and `pride::list_ftp_files` passed one, which turned a PRIDE outage that hangs into
  `MzLibError::Timeout`. mzLib already bounds each request (#1412: three retries after 5, 20 and
  60 s, each attempt up to 100 s), so they now default to `None` and a hang is reported as
  `MzLibError::ServiceUnavailable`, after about 8 minutes per request. A timeout you set is still
  sent: `Some(Duration::from_secs(300))` keeps the old limit. `ListOptions` and `SearchOptions`
  now derive `Default`. The `pride` guide gains *Retries and timeouts*: mzLib retries, how long an
  outage takes to report, why to wait minutes rather than seconds before trying again, and the two
  EBI failures that arrive as `MzLibError::Bridge`.
- **`peptidoform::fragments` takes its defaults from mzLib, so a call that names no protease or
  dissociation type gives different peptides and ions** (pyMzLib #83). Every digestion option of
  `FragmentOptions` is now an `Option` (`protease`, `dissociation`, `terminus`: `Option<String>`;
  `missed_cleavages`, `min_length`, `max_modifications`, `max_isoforms`: `Option<u32>`), `None` by
  default, and only the ones set are sent; an option left `None` takes mzLib's `DigestionParams`
  default. `fragments("P02768")` now sends nothing but the accession. Two of the old defaults
  were this crate's own and change:
  - **protease:** mzLib's default is `"trypsin"`, which also cuts before proline. The old default
    was `"trypsin|P"`, which does not. Albumin's default digest goes from 303 to 310 peptidoforms,
    and from 195 to 202 distinct peptides.
  - **dissociation:** mzLib names no default, so the bridge uses HCD (b and y ions) until it does
    (smith-chem-wisc/mzLib#1488). The old default was ETD (c and z• ions).

  The other defaults (missed cleavages, length window, modification and isoform caps, terminus)
  have the same values as before, now set by mzLib. `max_length: None` is no longer sent as `0`;
  `Some(0)` still asks for unbounded. **To keep the old results, pass
  `FragmentOptions { protease: Some("trypsin|P".into()), dissociation: Some("ETD".into()),
  ..Default::default() }`.** The docs' histone counts are re-measured at `d10017ce`; UniProt's
  annotations had moved them even under the old protease. Code that set a field must wrap it in
  `Some`: `FragmentOptions { max_modifications: Some(1), .. }`.
- **`pride::ListOptions::page_size` and `pride::SearchOptions::page_size` are `Option<u32>`, `None`
  by default, and leave the page size to mzLib** (pyMzLib #85). `None` sends no `--page-size`, so
  mzLib's `GetProjectFilesAsync` and `SearchProjectsAsync` defaults apply (100 today, as before).
  The results never depended on it. `Some(0)` and a size over `i32::MAX` are still refused before
  anything is spawned.
- **`stats::fit_with` with `trend: true` on fewer than 30 features sizes its trend as limma does**
  (mzLib 1.0.594, #1398; pyMzLib #82). This crate already sent `--spline-basis` only when
  `spline_basis` was set, but the bridge then asked mzLib for 4 basis functions every time. It now
  asks for none, and mzLib follows limma's `fitFDist` rule, so a small trend fit's prior and
  moderated statistics change to limma's. `VariancePrior::spline_basis_count` says what was used.
  A table of 30 or more features is unchanged, and so is any call that sets `spline_basis`.
  `FitOptions::spline_basis` no longer says the bridge's default is 4.
- **`peptidoform::convert_with` with `target: "ProForma"` names UniProt modifications by UNIMOD
  accession** (mzLib 1.0.594, #1403). `[UniProt:N-acetylserine on S]SEQK` becomes
  `[UNIMOD:1]-SEQK`, where it was written back by name; the caveat that said so is retired, so a
  result carries five caveats, not six. The guide's "ProForma does not resolve UniProt
  modifications yet" section now says what still differs: a modification no list mzLib loads knows
  is written back by name in ProForma. The `pro_forma` column `read_records` computes for a
  `.psmtsv` written without one gains the same accessions.
- PXD000001's file counts in the `pride` docs are dated: on 2026-10-21 PRIDE's REST manifest listed
  9 of the 14 files in the FTP tree, omitting five including the two largest; the recording the
  examples replay (2026-10-03) lists 8. The docs said 8 of 13.

- **`QuantifyOptions::use_pep_q_value` says what it does**: it filters nothing, and changes the
  q-value FlashLFQ carries, which match-between-runs reads. `bayesian_protein_quant`'s results are
  written to `output_directory` only (pyMzLib #72).
- **The hand-written `pride_ftp_PXD000001.json` is gone**: the FTP examples and tests replay
  pyMzLib's recording of the real tree, which holds 14 files.

- **Four shipped modifications now write their Unimod accession in `pro_forma`** (mzLib #1328):
  `GG (Ubiquitination Site)` as `[UNIMOD:121]`, both `Myristoylation` entries as `[UNIMOD:45]` and
  `EQIGG` as `[UNIMOD:846]`, instead of by name, in `read_records` on a MetaMorpheus `.psmtsv` or
  `.osmtsv`. Masses are unchanged.
- **A read fault on an existing `.mzid` is `MzLibError::Bridge` of type `MzLibException` naming the
  file** (mzLib #1362), where it was an `IOException`. A missing file is still a usage error.
- **PRIDE download errors name the file and the host, never the URL** (mzLib #1350), so a reviewer
  token in a query string cannot reach a log. They were, and still are, `ServiceUnavailable`.
- **Counts and mzLib versions are gone from the readers prose**: the `formats()` example shows the
  counts as executed output instead, so the next mzLib release cannot leave them stale.
- **`readers::read_matches_with` takes `MatchOptions`** (which nests `ReadOptions`, and adds
  `scores`) instead of `ReadOptions`: the spelling the bridge's spec records, and the only change
  here to an existing signature.
- `ProteinGroup::intensities` documents both keyings: run names from `quantify`, sample labels from
  `median_polish`.
- The `readers` result types (`ResultRecords`, `NativeRecords`, `FeatureRecords`, `MatchRecords`,
  `ScanRecords`) are explicit structs rather than macro-generated, so each field states its own unit:
  `record_count` counts scans for `read_spectra`, features for `read_features`, matches for
  `read_matches`. `Table` now implements `Deserialize` from the wire's `column_names` and `columns`.
  Their fields are unchanged.

### Deprecated

- **`peptidoform::Peptide::fixed_charges` and `Modification::formal_charge`** are always 0 from
  the pyMzLib 0.7.0 bridge on, and will be removed when the bridge's protocol next changes
  (pyMzLib #99). mzLib's corrected UniProt modifications carry each permanent charge in their
  mass, so `Peptide::monoisotopic_mass`, `Modification::mass` and `Fragment::neutral_mass` are one
  proton lower per charge for such a peptide, and `Peptide::mz` gives the same value as before.
  The docs say so; `Fragment::neutral_mass` no longer warns that a charged fragment's m/z cannot be
  computed, since `(neutral_mass + z * PROTON_MASS) / z` is now right for every fragment.

### Fixed

- Stale caveats: `QuantifyOptions::max_threads` warned that FlashLFQ's roll-up changed results
  between runs (mzLib#1111), and the peptidoform notes described spurious ETD `y` ions (#1109).
  Both were fixed inside the pinned mzLib, by #1155 and #1114. `max_threads: 1` is now described
  as the conservative choice until the bridge project re-measures the K562 case at `-1`, as the
  `quant flashlfq` spec says; `STATUS.md`, `docs/findings.md` and `docs/test-parity.md` record the
  fixes.
- The live test for a FLASHDeconv feature file asserted the fabricated zero intensity the bridge
  used to pass through. On the 1.0.592 bridge the column is named in `absent_fields` and every
  value is `None`, and the test says so.
- Stale counts: mzLib 1.0.592 recognises 36 file types, of which 17 offer no cross-format view and
  6 offer `spectral_match` (the crate said 31, 14 and 4). The README heading no longer counts
  modules.
- The CI job for the minimum supported Rust version resolves dependencies that support Rust 1.74
  rather than failing on the newest `thiserror`, which needs 1.77.

### Documentation

- **`docs/upgrading.md` has a pyMzLib 0.7.0 section** (mzLib 1.0.596), entry for entry with
  pyMzLib's, naming each changed answer by its Rust function and field. The `ptmqtl` guide
  describes `not_detected` as mzLib's covered-zero, the peptidoform counts and `test_status`; the
  `sdrf` guide keeps the SDRF's replicate numbers and reads `TMT131`; the `isobaric` guide shows
  the SDRF labels. The 0.7.0 recordings the repin changed are re-vendored from pyMzLib 0.7.0, byte
  for byte: `bridge_version.json`, the five `diff_design_*`, both `isobaric_kits*`, the three
  `ptmqtl_occupancy_*`, the three `ptmqtl_pairs_*`, both `ptmqtl_trait_effects*`,
  `ptmqtl/trait_occupancy.tsv` and the three `sdrf_design_*`, and `bridge_version`'s example
  names mzLib `cc67b7a2`.

- **Every result mzLib 1.0.594 and 1.0.595 changed is now documented** (as pyMzLib #87 does), in
  `docs/upgrading.md`'s pyMzLib 0.5.0 bridge section and on the reference pages: ProForma input
  refusals (#1462), SUMO and TMTpro names (#1432, #1419), the Lactylation mass (#1417), entrapment
  flags and `entrapment_fraction` (#1371), FlashLFQ normalization and gapped designs (#1422,
  #1424), `|P` proteases at sequence ends (#1320), `stats::fit` below 3 features (#1398), and
  PRIDE's retries, resumable `.partial` files and missing download checks (#1412, #1427, #1449).
  The specs are re-vendored at bridge `f85fef7` and the reference pages re-rendered. The `pride`
  guide says a failed transfer's `.partial` pair resumes and to check large files yourself, and
  `list_ftp_files` says an empty FTP listing now arrives as `MzLibError::Bridge`, not
  `ProjectNotFound`.
- **`VariancePrior::scale` says what it means below 3 features.** mzLib fits no prior there
  (#1398), so `VariancePrior::df` is 0 and each feature keeps its own variance, in
  `FitRow::prior_variance`. The pyMzLib 0.5.0 bridge still reported the first feature's own
  variance as `scale`; the pinned pyMzLib 0.6.0 bridge reports `None`, as the spec says, and
  `a_fit_of_fewer_than_three_features_fits_no_prior` holds it live.
- **The stats guide teaches contrasts, the marginal-likelihood prior, p-value combination,
  correlation, logistic and mixed models and entrapment**, each section with an executed example
  ported from pyMzLib's guide. Its "what it is not" no longer says contrasts are not implemented.
- **`stats::fit_with`'s recordings are re-vendored from pyMzLib 0.6.0**
  (`stats_fit_malat.json`, `stats_fit_limma.json`): the 0.6.0 bridge's first caveat names the
  estimator (`moments_legacy: ...`), its `residual_df_differ` caveat points to
  `marginal_likelihood`, and its `bh_adjusted` caveat says "per coefficient and per contrast", so
  `the_fit_recordings_are_what_the_bridge_returns_today` failed against the repinned bridge until
  they were. Every number is unchanged. `bridge_version.json` is re-vendored too (the 0.6.0
  bridge's verb list), and `bridge_version`'s example names its mzLib commit, `d10017ce`.

### Tested

- **stdin is never inherited** (pyMzLib #73). The Python bridge call let a call with nothing to send
  inherit the caller's stdin, so `median_polish(path)` hung in a terminal. This crate always pipes
  stdin and closes it, so it never had the bug; a test now re-runs the test binary with a stdin
  pipe that never closes and fails if a bridge call waits on it.
