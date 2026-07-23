# Test parity with pyMzLib

The goal stated for this crate was **"exactly the same three functions and tests in the Rust repo
that we have in pyMzLib."** This table is how that claim stays checkable rather than asserted.

pyMzLib has **123** tests (107 offline, 16 live). mzLibRust has **136** (118 offline, 18 live).
Every Python test maps to a Rust test, is eliminated by Rust's type system, or is listed below with
the reason it is not portable. Nothing was silently dropped.

## Counts

| Suite | pyMzLib | mzLibRust | Notes |
|---|---|---|---|
| transport | `test_bridge.py` 20 | `src/bridge.rs` 23 | +3: service-unavailable variant, version-verb assertion, stdin passthrough |
| PRIDE | `test_pride.py` 38 | `src/pride.rs` 38 | 1:1 |
| peptidoform | `test_peptidoform.py` 25 | `src/peptidoform.rs` 29 | +4: ETD y-ion defect (mzLib#1109) ×2, unbounded max-length, modification mass/charge |
| FlashLFQ | `test_flashlfq.py` 24 | `src/flashlfq.rs` 28 | +4: null peptide intensity, None/zero typing, every-flag-reachable, invariant formatting |
| PRIDE live | `test_pride_live.py` 7 | `tests/live_pride.rs` 7 | 1:1 (protocol handshake moved to `live_bridge.rs`) |
| peptidoform live | `test_peptidoform_live.py` 9 | `tests/live_peptidoform.rs` 9 | 1:1; 2 histone tests `#[ignore]` for slowness |
| transport live | — | `tests/live_bridge.rs` 2 | new: the M0 end-to-end proof |

## Tests eliminated by Rust's type system

These Python tests defend against a wrong *type* reaching a function. In Rust that is a compile
error, so the runtime test has nothing to assert. Each is listed so the coverage is visibly
accounted for, not quietly missing.

| pyMzLib test | Why it does not port |
|---|---|
| `test_unusable_timeouts_are_rejected_before_spawning_anything` | `Duration` cannot be negative, `inf`, `nan`, a string, or a bool. Only zero remains reachable, and `a_zero_timeout_is_rejected_before_spawning_anything` covers it. |
| `test_non_string_accession_gives_a_usage_error_not_an_attribute_error` | `&str` parameters cannot receive `123`, `None`, a list, or bytes. |
| `test_a_bare_string_of_extensions_is_refused` | `&[String]` is not `&str`; the confusion the test defends against cannot be expressed. |
| `test_spectra_must_be_a_list_not_a_string` | Same: `&[SpectraFile]` is not a string. |
| `test_download_files_rejects_things_that_are_not_pride_files` | `&[PrideFile]` is typed. |
| `test_bad_page_sizes_are_usage_errors` (partial) | `u32` excludes `"100"`, `None`, `2.5`; the two reachable cases (0, over `i32::MAX`) each have a test. |
| `test_negative_replicate_raises` | design fields are `u32`. `a_negative_replicate_is_impossible_by_construction` documents it. |
| `test_boolean_max_threads_raises` | `i32` is not `bool`. |
| `test_non_integer_max_length_is_refused` | `Option<u32>`. |
| `test_mz_rejects_charges_that_are_not_positive_whole_numbers` (partial) | `i32` excludes `1.5`, `"2"`, `True`, `None`; `0` and `-1` are tested. |

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
- **`etd_produces_c_and_z_ions_and_also_y_which_it_should_not`** pins a live mzLib defect
  (smith-chem-wisc/mzLib#1109) and fails the moment it is fixed upstream.

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
