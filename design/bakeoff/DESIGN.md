# mzLibRust bake-off — design and ground truth

Methodology follows pyMzLib's `design/bakeoff-flashlfq/DESIGN.md`: give an independent agent a real
task, one toolchain, and no hints, then measure whether they got a **correct** answer — not whether
they got *an* answer.

The product claim for these bindings is not "faster" or "smaller". It is **"you are less likely to
be confidently wrong."** So trap-avoidance is the primary metric, and everything else (steps, lines
of code, docs consulted) is secondary colour.

## Persona (all arms) — the "biologist phenotype"

**This is not optional, and getting it wrong invalidates the run.** Inherited verbatim from
pyMzLib's `design/bakeoff-flashlfq/DESIGN.md`, so results are comparable across bindings:

> A proteomics biologist: ~8 years at the bench and on an Orbitrap, comfortable with label-free
> quantification and **knows what match-between-runs is in principle**. **Not an experienced coder
> or Rust user** — leans on docs and examples, gets stuck on language mechanics, **does not read
> library source**.

The first attempt at this bake-off used "a proteomics researcher who writes some Rust" instead. The
agents read library source, reverse-engineered a dependency's internal types, and wrote 375–598
lines of verification scaffolding each — and the comparison became meaningless, because a careful
expert self-corrects regardless of how good the documentation is. It measured *can an expert reach
the right answer* rather than *do the docs carry a biologist through*.

The biologist phenotype is the harder and more honest test: it is the only one where the
documentation has to do the work.

## Instrumentation (mandatory, all arms)

One JSON object per attempt, appended **as you go**, to `<task>_arm{A,B}_log.jsonl`:

```json
{"n": 1, "action": "what you tried", "outcome": "worked|deadend|external_lookup", "note": "short"}
```

`external_lookup` = the arm had to leave the tool for something the tool should have supplied — a
function name, a parameter's meaning, how to switch MBR on, what a missing value means. **This is
the cleanest single signal**, and in pyMzLib's peptidoform run it separated the arms 3 vs 2.

Also required per arm: **"would you put it in a figure?"**, dead-end count, and the decisions the
arm had to make alone.

## Arms

Six agents, three tasks × two toolchains. Arm A is kept **blind** to mzLibRust — it is never
mentioned in an Arm A prompt. Each works alone, in its own scratch directory, and is
forbidden from reading the other arm's work. The mzLibRust arms are additionally forbidden from
reading `#[cfg(test)]` blocks, anything under `tests/`, and `docs/findings.md` — a real user does
not have the author's test suite, and letting an agent read the tests would leak every answer.

| Task | mzLibRust arm | Existing-tool arm |
|---|---|---|
| Label-free quant across two runs, MBR on | `mzlib::flashlfq` | `mzdata` / `sage` / whatever crates.io offers |
| Digest and fragment a protein | `mzlib::peptidoform` | `rustyms` |
| PRIDE manifest and filtered download | `mzlib::pride` | a PRIDE client crate if one exists, else hand-rolled HTTP |

The agents are **not** told what the traps are. They are asked for numbers, their confidence, their
judgement calls, and their surprises. Scoring happens afterwards against the ground truth below.

## Ground truth

Established independently of the Rust crate, by driving the bridge executable directly from the
shell. That matters: if the crate were wrong, ground truth derived *through* the crate would be
wrong in the same direction and the bake-off would score itself.

### Task 1 — quant, two K562 runs, match-between-runs on

Inputs: `AllPSMs.psmtsv` and `20100614_Velos1_TaGe_SA_K562_{3,4}.mzML` from mzLib's own FlashLFQ
test data — the same pair the pyMzLib bake-off used, so the two are directly comparable.

594 identifications → 354 peptides, 943 protein groups, 647 chromatographic peaks.

| Question | Correct | What the obvious-but-wrong route gives | Error |
|---|---|---|---|
| Peptides quantified in **both** runs | **257** (from `peaks`) | 169 (from the peptide roll-up) | −34% |
| Peptides rescued by MBR | **140** (from `peaks`) | 52 (from `detection_types`) | **−63%** |
| Protein groups with no obtainable intensity | **2** (`None`) | 0 (if `None` read as zero) or 849 (if `None` and `0.0` conflated) | — |

**The headline trap, verified on real data.** MBR transfers per run:

| | run_3 | run_4 |
|---|---|---|
| from `peaks` (correct) | 62 | 78 |
| from the peptide roll-up | **0** | 52 |

An entire run's 62 transfers are absent from `QuantifiedPeptides.tsv`. A user who builds their
matrix from the peptide table does not get a slightly low number — they get a result in which
match-between-runs appears not to have worked at all in half the experiment. This is precisely why
`FlashLfqResults::peaks` is documented as the MBR surface and `mbr_rescued_peptide_count` exists.

**The protein trap.** 2 protein groups are `None` — FlashLFQ's median-polish could not resolve them.
847 are `0.0` in both runs, which means *measured as absent*, an entirely different statement.
Conflating them mislabels 847 proteins; coercing `None` to `0.0` fabricates a number for 2. In
mzLibRust the type system does not permit either mistake silently: `ProteinGroup::intensities` is
`HashMap<String, Option<f64>>`.

### Task 2 — digest and fragment human serum albumin (P02768)

Peptide counts vary with two defaults, and both defaults move the answer:

| protease | min_length | peptides |
|---|---|---|
| `trypsin\|P` (Keil rule — what a mass spectrometrist means) | 7 | **195** |
| `trypsin` (also cleaves before proline) | 7 | 202 |
| `trypsin\|P` | 1 | 243 |
| `trypsin` | 1 | 257 |

- **The protease trap.** mzLib's `trypsin|P` *applies* the proline rule; its plain `trypsin` does
  not. That is the **reverse** of MaxQuant/Mascot, where `Trypsin/P` means the rule is *ignored*.
  Someone reaching for the familiar-looking name gets the opposite of their intent, in either
  direction, silently. (smith-chem-wisc/mzLib#1106)
- **The min-length trap.** The default of 7 discards 48 peptides here — a fifth of the digest — with
  no indication. On a histone it is roughly a third.
- **The census.** UniProt annotates **38** modification-like features; **14** are usable; **24** are
  glycosylation sites with no defined composition and therefore no mass. Reporting "14
  modifications" without the denominator is correct-but-misleading, which is what
  `ModificationCensus::explain()` exists to prevent.
- **The ETD trap.** mzLib emits `c`, `zDot` **and `y`** for ETD — the y ions are spurious and about
  a third of the fragment list. See smith-chem-wisc/mzLib#1109. An arm that reports "c and z ions"
  is repeating the textbook rather than reading its output; an arm that notices the y ions has done
  something genuinely good.

### Task 3 — PRIDE project PXD000001

8 files, 514,278,049 bytes (0.514 GB).

**The extension trap.** The MGF file is `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz`. Its extension is
**`.gz`**, not `.mgf`. Filtering on `.mgf` matches **zero** files. The failure mode that matters is
not the empty result — it is that an empty result plus a zero exit code reads as success, so a batch
script reports "done" having downloaded nothing. mzLibRust raises `MzLibError::Usage` naming the
`.gz` gotcha rather than returning an empty vector.

**The unknown-accession trap.** PRIDE answers `PXD0000019999` with an empty result, not a 404. A
binding that passes that through returns an empty list, and a typo becomes "0 files, done".
mzLibRust raises `MzLibError::ProjectNotFound`.

## Scoring

Per arm, per question:

- **Correct** — matches ground truth, and the arm can say why.
- **Confidently wrong** — a plausible number that is wrong, reported without hedging. The worst
  outcome, and the one the bindings exist to prevent.
- **Hedged wrong** — wrong, but the arm flagged uncertainty.
- **Could not answer** — a legitimate outcome for the existing-tool arms, and itself a finding.

Secondary: steps taken, lines written, docs consulted, dead ends, and — most useful for improving
the docs — **what each arm wished the documentation had told it**.

## What happens to the results

Per the standing project rules:

1. Any **mzLib** defect surfaced goes to `smith-chem-wisc/mzLib` as an issue. A workaround in a
   binding is not a fix.
2. Any **documentation** lesson is back-ported to pyMzLib, and later to mzLibR. The traps live in
   mzLib's behaviour, not in one language, so wording that rescues a Rust user rescues a Python one
   too.

Results land in `RESULTS.md` beside this file.
