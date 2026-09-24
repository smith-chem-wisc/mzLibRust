# Reference facts: specs, fragments, the doc lint and the replay bridge

Every function in this crate that calls a wire verb documents that verb's **facts** — its
parameters with their units and valid ranges, its result fields with their units and what a null
means, its error kinds, caveats, the mzLib code it wraps, citations, the same verb's spelling in
Python and R, and the version that shipped it. Those facts are not written here. They are written
once, per verb, in the bridge repository, so that pyMzLib, mzLibRust and mzLibR cannot disagree
about them. This page is how that works, and what to do when you add or change a verb.

The rule from the bridge (`design/verbs/README.md`) in one line: **the spec owns the facts, the
binding owns the idiom.** Where this crate's docs contradict a spec, this crate is wrong. Where a
spec contradicts what the bridge emits, the spec is wrong, and it is fixed in the bridge.

## The pieces

| Piece | Where | What it does |
|---|---|---|
| Specs | `docs/specs/*.yaml`, `docs/specs/SOURCE` | One YAML file per wire verb, copied from the bridge's `design/verbs/`. Never edited here. |
| Sync script | `scripts/sync-specs.ps1` | Copies the bridge's committed specs byte for byte, deletes ones that are gone, and records the bridge commit in `SOURCE`. |
| Fragments | `docs/reference/<verb>.md`, `.bulk.md`, `.see-also.md` | Markdown rendered from each spec, included into rustdoc with `#[doc = include_str!(...)]`. Committed. |
| Renderer and lint | `tests/spec_docs.rs` | Renders the fragments, fails when a committed one is stale, and holds every documented parameter and field to its spec. |
| Replay bridge | `tools/replay-bridge/` | A stand-in bridge the doc examples run against, answering from `tests/fixtures/`. |

## What a function's page looks like

The shared reference-page order, the same in all three bindings, is: summary, wraps, parameters,
returns, errors, caveats, an **executed** example, performance, the same verb in other bindings,
references, since. A function that binds a verb is written as:

```rust,ignore
/// One-line summary, in this crate's words.
///
/// Hand-written prose: when to use it, the traps, links to the neighbours.
#[doc = include_str!("../docs/reference/readers.read-spectra.md")]   // wraps .. caveats
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let scans = mzlib::readers::read_spectra_with("sliced_ethcd.mzML", &options)?;
/// assert_eq!(scans.scan_count, 6);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-spectra.see-also.md")]  // other bindings .. since
pub fn read_spectra_with(/* .. */) { /* .. */ }
```

A `_many` function includes `<verb>.bulk.md` the same way. `every_shipped_verb_includes_its_reference_facts`
fails when a function that binds a verb is missing its includes.

Rust-specific errors that a spec cannot know about (`ProjectNotFound`, a refusal before anything is
spawned) go under a hand-written `# Errors this crate adds` heading, after the included facts.

## When the bridge's specs change

```powershell
.\scripts\sync-specs.ps1 -From E:\CodeReview\bridge\design\verbs
$env:MZLIB_RENDER_SPEC_DOCS = 1; cargo test --test spec_docs; Remove-Item Env:MZLIB_RENDER_SPEC_DOCS
cargo test --test spec_docs
```

The first command vendors what is **committed** at the bridge's HEAD; `-WorkingTree` takes
uncommitted specs too, and `SOURCE` then says so. The second re-renders every fragment. The third is
what CI runs: it fails on a stale fragment and on any documentation that no longer matches a spec.

Commit `docs/specs/` and `docs/reference/` together, with the code change that needed them.

## The doc lint

`every_param_and_field_is_documented_with_its_unit` reads this crate's source with `syn` and, for
every spec whose Rust function exists:

- every spec **parameter** must be a documented field of the options struct the spec names (searched
  through nested option structs, so `SpectraOptions::read::limit` counts), or an argument of the
  function;
- every **result field** must be a documented field of the struct the function returns (searched
  through `#[serde(flatten)]` and nested structs), and every entry field of a record list
  (`Vec<Format>`, `Vec<PrideFile>`) or a typed sub-table must be one on the entry's struct;
- a parameter or field whose spec `unit` is not null must **name its unit** in its doc — the unit,
  its singular, or a spelled-out alias (`min` → minute). The rule is pyMzLib's, so the two
  bindings agree on what counts.

A spec that says mzLibRust ships a verb (`since.mzlibrust` set) whose function does not exist fails.

### Deviations: a deliberate difference, with its reason

This crate keeps its own idiom. `QuantifyOptions` uses FlashLFQ's parameter names (`ppm_tolerance`,
not the wire's `--ppm`); `pride::list_files` returns the list itself rather than an envelope around
it. Each such difference is one entry in `DEVIATIONS` in `tests/spec_docs.rs`:

```rust,ignore
dev("quant flashlfq", "param.ppm", Some("ppm_tolerance"), "FlashLFQ's name, echoed by parameters.ppm_tolerance"),
dev("readers formats", "field.format_count", None, "formats() returns Vec<Format>; the count is its len()"),
```

`Some(name)` is the Rust spelling; `None` means the fact is not a named field or argument at all.
The reason is required. `every_deviation_names_a_real_param_or_field` fails on an entry that names
nothing in its spec, so a stale entry cannot survive a spec change. The rendered tables show the
Rust spelling, with the wire name beside it.

Deviations are reported back to the bridge on its thread, so its specs can record each binding's
spelling in `bindings.rust`.

### Pending: a gap, not a choice

`PENDING` lists spec facts this crate does not project **yet**. The lint skips them, and fails on an
entry that is no longer a gap, so the list can only shrink. A new verb whose function does not exist
yet, and whose spec does not claim mzLibRust ships it, is skipped without an entry.

## The replay bridge

Examples are doctests, and doctests run in CI. They need a bridge, and a real one needs .NET, the
network, and data files. So every runnable example starts with a hidden line:

```rust,ignore
# mzlib_replay::activate();
```

which points `MZLIB_BRIDGE` at an executable built by `tools/replay-bridge` and moves the process
into a scratch directory (so an example that writes `out` writes there). The crate cannot tell it
from the real bridge. It answers a call from a fixture in `tests/fixtures/` — recorded from the real
bridge, and shared byte for byte with pyMzLib and mzLibR — **only when the recording fits the call**:

- `--path` must name the file the recording was made from (compared by file name), and any other
  `--option value` whose name is a top-level key of the recording must equal it;
- `--limit`/`--offset` must reproduce the recording's `returned_count` from its `record_count` (or
  `row_count`), and a `--flag` with a `<flag>_included` key must match it;
- a `--paths-stdin` call fits only a bulk recording (a `files[]` whose entries carry a `path`, and no
  top-level `path`), and a one-path call only a one-document recording.

Otherwise the example fails with a usage error naming every recording it tried and why each did not
fit. These are the rules of pyMzLib's `pkg/python/tests/replay_bridge.py`, line for line.

Which recordings a verb may answer from is the `examples` of its spec, plus `REPLAY_EXTRA` in
`tools/replay-bridge/build.rs` for recordings no spec lists. The build script parses the fixtures and
compiles `src/stub.rs` — std only — with the table baked in, using the same `rustc` that builds the
crate. It is a path-only dev-dependency, so it never reaches crates.io.

An example that cannot run — one that downloads from EBI — is `no_run`, and the prose above it says
why.

## Adding a verb

1. The spec lands in the bridge first, with its fixture recorded from the live bridge.
2. Copy the fixture into `tests/fixtures/` byte for byte from pyMzLib, then sync the specs.
3. Write the function, its options struct and its result type. Document every field with its unit.
4. Include the fragments, write an example that replays the fixture, and run `cargo test`.
5. Add a `CHANGELOG.md` entry, and report the Rust spelling and version back to the bridge.
