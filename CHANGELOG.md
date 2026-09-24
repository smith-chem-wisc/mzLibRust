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
