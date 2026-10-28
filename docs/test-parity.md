# Test parity with pyMzLib

The goal stated for this crate was **"exactly the same three functions and tests in the Rust repo
that we have in pyMzLib."** This table is how that claim stays checkable rather than asserted.

pyMzLib has **123** tests (107 offline, 16 live). mzLibRust has **136** (118 offline, 18 live).
Every Python test maps to a Rust test, is eliminated by Rust's type system, or is listed below with
the reason it is not portable. Nothing was silently dropped.

> **These counts, and the table below, cover the first three modules** (PRIDE, peptidoform,
> FlashLFQ), when they were ported. The `readers`, `sdrf` and `proteins` modules and the mzLib
> 1.0.592 batch came later and are mapped by their specs instead: `tests/spec_docs.rs` fails when a
> spec's parameter or field is not documented here, and every doc example replays the recording
> pyMzLib's own example replays. At the 1.0.592 port the crate has 231 offline unit tests, the
> 7-test spec lint, 69 doctests (60 executed; 9 are `no_run` because they download from EBI or no
> many-file recording exists yet) and 52 live canaries.

## Counts

| Suite | pyMzLib | mzLibRust | Notes |
|---|---|---|---|
| transport | `test_bridge.py` 20 | `src/bridge.rs` 23 | +3: service-unavailable variant, version-verb assertion, stdin passthrough |
| PRIDE | `test_pride.py` 38 | `src/pride.rs` 38 | 1:1 |
| peptidoform | `test_peptidoform.py` 25 | `src/peptidoform.rs` 29 | +4: the default recording's HCD product set (b and y, no c or z•), an absent max-length is not sent while 0 is unbounded, modification mass/charge |
| FlashLFQ | `test_flashlfq.py` 24 | `src/flashlfq.rs` 28 | +4: null peptide intensity, None/zero typing, every-flag-reachable, invariant formatting |
| PRIDE live | `test_pride_live.py` 7 | `tests/live_pride.rs` 7 | 1:1 (protocol handshake moved to `live_bridge.rs`) |
| peptidoform live | `test_peptidoform_live.py` 9 | `tests/live_peptidoform.rs` 9 | 1:1; 2 histone tests `#[ignore]` for slowness |
| transport live | — | `tests/live_bridge.rs` 2 | new: the M0 end-to-end proof |
| peptidoform convert | `test_peptidoform_convert.py` 18 (8 functions; one parametrized over 11 refused arguments) | `src/peptidoform.rs` 8 | 1:1 by function; 4 of the 11 refused arguments cannot be expressed in Rust (below), and `\r` is refused beside `\n` |
| peptidoform convert live | — | `tests/live_peptidoform.rs` 5 | new: the recordings against the live bridge, the four modes, format names, thread-count invariance. Skip on a bridge without the verb (older than pyMzLib 0.4.0) |
| stats, pyMzLib 0.6.0 verbs | `test_stats.py` 25, `test_stats_combine_correlate.py` 10, `test_stats_models.py` 13 | `src/stats.rs` 26 | by behaviour, not by function: every request shape (options sent only when given), every refusal before spawning, every recording parse and the old-bridge guard for `contrasts` and `estimator` has a test; one Rust test often holds several pyMzLib ones. The type checks and `__all__` are eliminated (below) |
| diff | `test_diff.py` 16 | `src/diff.rs` 11 | by behaviour, as stats; `test_the_example_inputs_are_vendored` is the doctests themselves, which fail without the inputs |
| ptmqtl | `test_ptmqtl.py` 16 | `src/ptmqtl.rs` 9 | by behaviour, as stats |
| stats live, pyMzLib 0.6.0 | `test_stats_live.py` 10, `test_stats_combine_correlate_live.py` 3, `test_stats_models_live.py` 5 | `tests/live_stats.rs` 20 | every parametrized recording case, every R comparison (limma, metafor, `pchisq`/`qnorm`, `glm`, `nlme`, the entrapment reference script) and every refusal; +1: a fit of fewer than three features fits no prior (`scale` is `None`) |
| diff live | `test_diff_live.py` 5 | `tests/live_diff.rs` 5 | 1:1 |
| ptmqtl live | `test_ptmqtl_live.py` 4 | `tests/live_ptmqtl.rs` 4 | 1:1 |
| PRIDE checksums | `test_pride.py` 10 checksum tests, and `file_checksums` in the no-timeout parametrization | `src/pride.rs` 11 | 1:1 by function, plus the manifest-versus-checksum size. The verb guard's `("pride checksums", "0.6.0")` is held by `bridge.rs`'s `require_verb` tests, not per verb; a non-string accession cannot reach `&str` (above) |
| sdrf isobaric design | `test_sdrf_isobaric_design.py` 11 offline (one parametrized over 10 refused arguments) | `src/sdrf.rs` 9 | 1:1 by function; the blank path is one of the refused arguments here, and the older-bridge refusal is `bridge.rs`'s `require_verb` test. `kit=10` and a bare-string `condition_columns` cannot be expressed in Rust; an empty `condition_columns` is this crate's "not sent", as for `design_with` |
| sdrf isobaric design live | `test_sdrf_isobaric_design.py` 2 live | `tests/live_sdrf.rs` 2 | 1:1: the recording and `TmtDesign.txt` against the 0.6.0 bridge, and an unknown kit is `isobaric kits`' own usage error |
| sdrf restrict | `test_sdrf_restrict.py` 7 offline (two parametrized over 9 and 3 refused arguments) | `src/sdrf.rs` 5, and `written_sdrf_reads_the_path_mzlib_reports` for `pool` | 1:1 by function. A bare string, a number and an empty mapping cannot reach `&[SearchedFile]`; a negative window cannot reach `u64` |
| sdrf restrict live | `test_sdrf_restrict.py` 4 live | `tests/live_sdrf.rs` 4, and `pool_out_reports_the_file_mzlib_wrote` | 1:1: the recording, `out`'s `.sdrf.tsv` rename, restrict then design with native paths (SDRF-B5), and no `comment[data file]` is `ArgumentException` |
| proteins align and parse-accession | `test_proteins_align.py` 17 offline (five of them parametrized) | `src/proteins.rs` 15 | 1:1 by function. A non-int position, a non-string pair or accession, and a bare string cannot reach `usize`, `&[(S, T)]` or `&[S]`; the verb guard is `bridge.rs`'s |
| proteins align and parse-accession live | `test_proteins_live.py` 6 (one parametrized over 4 refusals) | `tests/live_proteins.rs` 6 | 1:1 |
| genes orthology | `test_proteins_orthology.py` 15 offline (two parametrized over 15 and 5 cases) | `src/proteins.rs` 12 | 1:1 by function. An int release, a dict instead of a list and a non-bool flag cannot be expressed in Rust (`&str`, `Vec<(String, PathBuf)>`, `bool`) |
| genes orthology live | `test_proteins_live.py` 5 | `tests/live_proteins.rs` 5 | 1:1: the snapshot id over the byte-exact inputs, Parquet written once, mzLib's three refusals |
| library read | `test_library.py` 18 offline (two parametrized over refused arguments) | `src/library.rs` 14 | By function, with three exceptions: a list passed to `read` cannot reach `impl AsRef<Path>` (nor a bool or negative threshold, limit or offset reach `Option<u64>`/`u64`), the `__dir__` test has no Rust counterpart (a module's `pub` items are its API), and the per-row labels and retention-time units are pinned by the module guide's doctests |
| library read live | `test_library_live.py` 9 (one parametrized over 4 recordings) | `tests/live_library.rs` 8 | 1:1; the verb-list check is `require_verb` itself |
| PRIDE checksums live | `test_pride_live.py` 3 | `tests/live_pride.rs` 3 | 1:1: the list covers the FTP root's deposits, one download's length and MD5 match its row, an unknown accession is `ProjectNotFound`. Skip on a bridge older than pyMzLib 0.6.0 |

## Tests eliminated by Rust's type system

These Python tests defend against a wrong *type* reaching a function. In Rust that is a compile
error, so the runtime test has nothing to assert. Each is listed so the coverage is visibly
accounted for, not quietly missing.

| pyMzLib test | Why it does not port |
|---|---|
| `test_unusable_timeouts_are_rejected_before_spawning_anything` | `Duration` cannot be negative, `inf`, `nan`, a string, or a bool. Only zero remains reachable, and `a_zero_timeout_is_rejected_before_spawning_anything` covers it. |
| `test_non_string_accession_gives_a_usage_error_not_an_attribute_error` | `&str` parameters cannot receive `123`, `None`, a list, or bytes. |
| `test_a_bare_string_of_extensions_is_refused` | `&[String]` is not `&str`; the confusion the test defends against cannot be expressed. |
| `test_bad_arguments_are_refused_before_the_bridge`: one string, a non-string sequence, `mode=None`, `threads=True` | `convert` takes `&[S: AsRef<str>]`, so neither a bare `&str` nor a number can reach it; `ConvertOptions::mode` is a `ConversionMode` and `threads` an `i32`. An unknown mode name is refused by `ConversionMode::from_str`. |
| `test_spectra_must_be_a_list_not_a_string` | Same: `&[SpectraFile]` is not a string. |
| `test_download_files_rejects_things_that_are_not_pride_files` | `&[PrideFile]` is typed. |
| `test_bad_page_sizes_are_usage_errors` (partial) | `Option<u32>` excludes `"100"` and `2.5`; `None` is not a bad size but mzLib's default, as in pyMzLib since #85. The two reachable cases (`Some(0)`, over `i32::MAX`) each have a test. |
| `test_negative_replicate_raises` | design fields are `u32`. `a_negative_replicate_is_impossible_by_construction` documents it. |
| `test_boolean_max_threads_raises` | `i32` is not `bool`. |
| `test_non_integer_max_length_is_refused` | `Option<u32>`. |
| `test_mz_rejects_charges_that_are_not_positive_whole_numbers` (partial) | `i32` excludes `1.5`, `"2"`, `True`, `None`; `0` and `-1` are tested. |
| `test_combine_p_has_no_default_method`, `test_entrapment_fdp_has_no_default_ratio`, `test_pairs_requires_a_type`, `test_trait_effects_needs_both_paths` | `method`, `ratio`, `pair_type` and both paths are positional arguments, so leaving one out is a compile error. `observations_with`'s `basis` is a `QuantBasis` for the same reason. |
| `test_module_lists_only_its_public_api`, `test_the_module_is_exported`, `test_module_is_exported` (stats, diff, ptmqtl) | Rust's `pub` is the public API; `lib.rs` declares the modules, and the doctests use them by path. |
| the stats, diff and ptmqtl `test_*_refuses_*` cases for a wrong type (a bare string for a list, a bool for an int, `threads=True`, a non-path) | The arguments are `&[S: AsRef<str>]`, `&[PValueTest]`, `&[ObservationRun]`, `Option<u32>`, `i32`, `f64` and `impl AsRef<Path>`: the confusion cannot be expressed. The reachable cases (an empty list, a blank or tab-holding name, a non-finite number, an unknown estimator or grouping, both filters) each have a test. |

## Deliberately not ported

| pyMzLib | Reason |
|---|---|
| `PrideFile.as_dict()` and `test_as_dict_includes_the_computed_properties` | Exists so a pandas `DataFrame` picks up the computed properties, which `vars()` and `asdict()` skip. Rust has no equivalent hazard — a caller maps the struct explicitly — so the method and its test would be ceremony. |
| `test_bridge_failure_surfaces_as_bridge_error` (in `test_pride.py`) | Covered once at the transport layer in `src/bridge.rs` rather than repeated per tranche. |

## Where the Rust suite is stronger

- **`Option<f64>` for protein intensity is enforced by the compiler**, not by documentation. The
  Python docs warn that a protein intensity may be `None`; here you cannot read one without
  handling it. `the_none_and_zero_distinction_is_carried_by_the_types` asserts the shape.
- **A `null` peptide intensity is resolved at the boundary.** pyMzLib documents peptide intensities
  as "0.0 when missing, never None" but stores the wire value verbatim — and its own
  `flashlfq_small.json` fixture contains `"run_4": null`, so `intensity("run_4")` returns `None`
  there. No Python test covers it. See [findings.md](findings.md).
- **`etd_produces_c_and_z_ions_not_b_and_y`** (live) replaced the two tests that pinned
  smith-chem-wisc/mzLib#1109 while it was open; #1114 fixed it, and the test now fails if the y
  ions ever come back. It was an offline test on the albumin recording until that recording was
  re-made at mzLib's defaults (pyMzLib #83), which fragment by HCD; the offline
  `hcd_produces_b_and_y_but_no_c_or_zdot_ions` now pins that series instead.

## Running them

```bash
cargo test                     # the offline suite: no network, no bridge needed
cargo fmt --check && cargo clippy --all-targets

# Live canaries. MZLIB_BRIDGE must point at a staged bridge; they SKIP on an outage.
MZLIB_BRIDGE=/path/to/mzlib-bridge cargo test --features live
MZLIB_BRIDGE=/path/to/mzlib-bridge cargo test --features live -- --ignored   # the slow histones
```

`cargo test` has no skip verdict, so a skipped live canary prints `SKIPPED: …` on stderr and
passes. That is weaker than pytest's `skip` and is noted rather than papered over; the alternative —
failing on an outage — is exactly what this convention exists to prevent.
