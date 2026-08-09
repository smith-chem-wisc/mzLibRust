# mzLibRust — status

Written 2026-07-23, at the end of the session that built it.

## Where it stands

**The crate exists, works, and is public:** <https://github.com/smith-chem-wisc/mzLibRust>.

The three capabilities pyMzLib has — PRIDE, peptidoforms, FlashLFQ — over the same language-neutral
bridge. **137 tests** (119 offline, 18 live). `cargo test` is green with no network and no .NET on
the machine; `cargo clippy --all-targets` and `cargo fmt --check` are clean.

Anyone can now use it. `build.rs` resolves a bridge from `MZLIB_BRIDGE`, from `_dotnet/<rid>/`, or by
download, and never fails the build when it finds none — a missing payload is a runtime error with
three named remedies. `scripts/stage-bridge.ps1` stages one from a pyMzLib checkout and probes it.

**What is not done:** publishing to crates.io, which waits on pyMzLib's CI uploading the raw bridge
binaries as release assets (mzLib#1103 makes that cheaper by shrinking the payload). The machinery in
`build.rs` is already there.

## What the port produced beyond the crate

Seven defects, all verified against ground truth before filing, two with fixes proposed:

| | | |
|---|---|---|
| [mzLib#1109](https://github.com/smith-chem-wisc/mzLib/issues/1109) | ETD/ECD emit `y` with no `b` — no fragmentation mechanism | **[PR #1114](https://github.com/smith-chem-wisc/mzLib/pull/1114)** |
| [mzLib#1110](https://github.com/smith-chem-wisc/mzLib/issues/1110) | z• suppressed at proline, complementary c ions still emitted | open |
| [mzLib#1111](https://github.com/smith-chem-wisc/mzLib/issues/1111) | FlashLFQ roll-up nondeterministically drops MBR intensities | open |
| [mzLib#1112](https://github.com/smith-chem-wisc/mzLib/issues/1112) | Glycosylation annotations dropped silently (corrected: exclusion right, silence wrong) | open |
| [mzLib#1113](https://github.com/smith-chem-wisc/mzLib/issues/1113) | Modifications applied after digestion → peptides trypsin cannot make | open |
| [pyMzLib#7](https://github.com/smith-chem-wisc/pyMzLib/issues/7) | `Peptide.intensity()` returned `None` against its own invariant | **[PR #9](https://github.com/smith-chem-wisc/pyMzLib/pull/9)** |
| [pyMzLib#8](https://github.com/smith-chem-wisc/pyMzLib/issues/8) | `--no-modifications` also discards proteolysis products | **closed / fixed** |

**Two deserve attention beyond the issue tracker:**

- **#1111 threatens reproducibility.** With default threading, identical inputs give different
  protein-level answers roughly 1 run in 6. Any FlashLFQ figure produced multithreaded may not
  reproduce. `max_threads = 1` is the workaround and is now documented in both bindings.
- **#1112 affects MetaMorpheus, not just the bindings.** Every UniProt XML load silently drops
  glycosylation-site annotations, and reports nothing about it.

## What the bake-off actually showed

Six biologist-persona arms, three tasks × two toolchains
([design](design/bakeoff/DESIGN.md), [results](design/bakeoff/RESULTS.md)):

| | ecosystem | mzLibRust |
|---|---|---|
| external lookups | 13 | **2** |
| dead ends | 8 | **1** |
| answers they would publish | 4 of 11 | **10 of 11** |

The quant arm is the clearest: **284 / 173 / 8** against a truth of **257 / 140 / 2** — wrong in the
plausible direction, which is the dangerous one. Their diagnosis was not about Rust:

> The stuck part wasn't Rust… it was realizing an hour in that **there is no tool here**. I now have
> to be the person who decides what "10 ppm" and "±1 minute" mean, and I have no way to tell if I
> chose wrong.

The strongest single result was a trap the mzLibRust arm **hit** rather than avoided — filtering on
`.mgf` for a file named `.pride.mgf.gz`. The library refused instead of reporting success with an
empty list. *"In a Python script I'd have gotten an empty list, shrugged, and found out three weeks
later."* A trap avoided proves the docs are readable; a trap hit and caught proves the guard works.

## Three things worth carrying to the next binding

**1. The persona is the experiment.** The first bake-off used an expert-programmer persona and
produced "nobody was confidently wrong" — an artifact, since careful experts self-correct regardless
of documentation quality. It had to be discarded and re-run. Expert arms are good **defect
scanners**; biologist arms are the only valid **comparison**. Run both, in that order. Pinned in
`DESIGN.md` so it cannot drift again.

**2. Put the warning where the mistake is made, not where the concept is defined.** Confirmed three
times. The `.gz` trap was documented on `PrideFile::extension` — but the person about to make the
mistake is reading `DownloadOptions::extensions`, whose examples were both uncompressed. The
biologist walked straight in. Conversely `Peptide::intensities` carries its warning on the field
itself, and the biologist reported 257 instead of 169 because of it.

**3. A disclosure that states a false reason is worse than no disclosure.**
`ModificationCensus::explain()` was wrong twice in one day — first claiming excluded sites had no
defined composition, then implying they should have been loaded. Both times the feature built to
prevent silent wrongness *was* the wrong thing. The correction came from a domain expert, not from
tests, tooling, or two agents independently agreeing (they agreed, and were both wrong).

## Immediate next steps

1. **Review [mzLib#1114](https://github.com/smith-chem-wisc/mzLib/pull/1114).** It changes ETD search
   results, so it is a scientific call as much as a code one and wants a second opinion from someone
   who runs ETD routinely. Full suite: 5334 passed.
2. **Review [pyMzLib#9](https://github.com/smith-chem-wisc/pyMzLib/pull/9).** One real bug fix plus
   six doc corrections, 154 passed.
3. **Publish the bridge binaries** as pyMzLib release assets — the last thing between this crate and
   crates.io.
4. **`ModificationCensus` should carry the annotations, not a tally.** It cannot currently
   distinguish "no defined composition" from "in vitro" from "variant-only", which are three
   exclusions needing three different judgements. Needs a bridge change, so both bindings gain it.
5. **mzLibR** is the same pattern and about as cheap — `processx` + `jsonlite` over the same wire.
   Every documentation lesson above transfers.
