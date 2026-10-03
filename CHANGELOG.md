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

### Tested

- **stdin is never inherited** (pyMzLib #73). The Python bridge call let a call with nothing to send
  inherit the caller's stdin, so `median_polish(path)` hung in a terminal. This crate always pipes
  stdin and closes it, so it never had the bug; a test now re-runs the test binary with a stdin
  pipe that never closes and fails if a bridge call waits on it.
