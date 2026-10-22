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

- **The bridge pin is pyMzLib v0.6.0** (still mzLib 1.0.595, `d10017ce`; bridge-watch's repin).
  Its bridge dispatches 53 verbs. **Specs re-vendored in full at bridge `60b1d25`**: every spec
  but `modifications resolve`, which has no verb (`not_vendored` in `docs/specs/SOURCE`), so the
  verbs this crate has not built yet are listed on the front page as `(not yet)`. Every spec's
  recorded examples are vendored byte for byte from pyMzLib 0.6.0, because the replay bridge
  reads each spec's examples, and so are the three recordings that changed since 0.5.0
  (`bridge_version.json`, `stats_fit_malat.json`, `stats_fit_limma.json`: `stats fit`'s first
  caveat now names its estimator). The replay bridge carries pyMzLib's whole `RECORDED_WITH`
  table, flags included. `stats fit`'s contrasts and estimator are `PENDING` in
  `tests/spec_docs.rs` until they are projected.
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
  `FitRow::prior_variance`. The pinned pyMzLib 0.5.0 bridge still reports the first feature's own
  variance as `scale`; bridges from pyMzLib 0.5.1 report `None`, as the spec now says.

### Tested

- **stdin is never inherited** (pyMzLib #73). The Python bridge call let a call with nothing to send
  inherit the caller's stdin, so `median_polish(path)` hung in a terminal. This crate always pipes
  stdin and closes it, so it never had the bug; a test now re-runs the test binary with a stdin
  pipe that never closes and fails if a bridge call waits on it.
