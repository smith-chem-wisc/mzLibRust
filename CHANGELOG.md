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
  only when a recording fits the call. 27 of the crate's 30 examples now execute, up from 0 of 10;
  the other three download from EBI and say so.
- The recorded fixtures of pyMzLib's mzLib 1.0.592 batch, byte for byte, in `tests/fixtures/`.
- `[package.metadata.docs.rs]`, so docs.rs builds the documentation with every feature.
- This changelog.

### Changed

- The `readers` result types (`ResultRecords`, `NativeRecords`, `FeatureRecords`, `MatchRecords`,
  `ScanRecords`) are explicit structs rather than macro-generated, so each field states its own unit:
  `record_count` counts scans for `read_spectra`, features for `read_features`, matches for
  `read_matches`. `Table` now implements `Deserialize` from the wire's `column_names` and `columns`.
  Their fields are unchanged.

### Fixed

- Stale counts: mzLib 1.0.592 recognises 36 file types, of which 17 offer no cross-format view and
  6 offer `spectral_match` (the crate said 31, 14 and 4). The README heading no longer counts
  modules.
- The CI job for the minimum supported Rust version resolves dependencies that support Rust 1.74
  rather than failing on the newest `thiserror`, which needs 1.77.
