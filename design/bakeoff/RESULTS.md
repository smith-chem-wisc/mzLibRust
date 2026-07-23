# Bake-off results

> ## ⚠ Wrong persona — the comparison below is not valid as a bake-off
>
> This run used **"a proteomics researcher who writes some Rust"**. The lab's established
> phenotype, defined in `pyMzLib/design/bakeoff-flashlfq/DESIGN.md` and used for every previous
> bake-off, is a biologist who is **"not an experienced coder … leans on docs, does not read library
> source."**
>
> The agents here did the opposite: they read library source, reverse-engineered sage's internal
> types, and wrote 375–598 lines of verification scaffolding apiece. This run also omitted the
> mandatory `external_lookup` JSONL instrumentation and the "would-you-put-it-in-a-figure" metric.
>
> **What that invalidates:** the head-to-head comparison. "Nobody was confidently wrong" is an
> artifact of using careful experts — they self-corrected by construction. The effort counts too: a
> biologist does not write 806 lines, they give up or ship something wrong. This tested *can an
> expert reach the right answer* (yes, on both toolchains) rather than *do the docs carry a
> biologist through*, which is the question the methodology exists to answer.
>
> **What still stands:** every defect in "Defects filed" and every entry in the documentation table.
> Those were verified independently against ground truth and do not depend on the persona.
>
> A corrected MBR head-to-head with the true biologist phenotype is running; results in
> `mbr_armA_log.jsonl` / `mbr_armB_log.jsonl` and a revision below.


Six independent agents, three tasks × two toolchains, scored against the ground truth in
[DESIGN.md](DESIGN.md) — which was established by driving the bridge executable directly, so the
crate could not score itself.

## Scorecard

| Task | Question | Truth | mzLibRust | Existing tools |
|---|---|---|---|---|
| **Quant** | peptides in both runs | 257 | **257** ✅ | 333 — different engine |
| | peptides rescued by MBR | 140 / 135 strict | **135** ✅ (found the distinction) | 211, self-assessed as meaningless |
| | protein groups with no intensity | 2 | **2** ✅ | 1 |
| **Peptidoform** | distinct tryptic peptides | 195 / 193 | **193** ✅ | 184 — mature chain, defensible |
| | ETD ion series | c, zDot, **y** (defect) | **found the y ions** ✅ | **found the y ions** ✅ |
| | usable modifications | 14 of 38 *per mzLib* | 14/38 ✅ *per mzLib* | **36 of 38** — and it was right |
| **PRIDE** | file count | 8 *per API* | **8** ✅ | **13** — and it was right |
| | total size | 0.514 GB *per API* | **0.514 GB** ✅ | **1.44 GB** — and it was right |
| | compressed MGF | 1 file | **1** ✅, avoided the trap | **1** ✅, caught the trap by luck |

**Nobody was confidently wrong.** Both toolchains were driven by agents who distrusted their first
answer and verified it. That is the honest headline, and it is not the headline we expected.

## What mzLibRust won on

**Effort, by a wide margin.** Lines written to answer the same questions:

| Task | mzLibRust | Existing tools |
|---|---|---|
| Quant | 186 (60 load-bearing) | **806** |
| Peptidoform | 222 (~25 load-bearing) | **598** |
| PRIDE | 117 (~30 load-bearing) | **604** |

The existing-tool arms did not spend that on ceremony. They spent it on a MetaMorpheus psmtsv
parser, a UniProt XML feature scraper, a `ptmlist.txt` parser, an RT-alignment replacement, and a
UniProt→PSI-MOD→Unimod name bridge — real domain code, each piece a place to be silently wrong.

**Trap avoidance where it counted most.** The quant arm reported 257 rather than 169 for
"peptides in both runs", and said why: the module header's third paragraph told it the peptide
roll-up drops MBR transfers. Its own words — *"had I reached for the peptide table by analogy with
`QuantifiedPeptides.tsv`, which is what a FlashLFQ user would do, I would have reported 169 and
never known."* The PRIDE arm never hit the `.mgf` filter that matches nothing, and named the
doc-comment that stopped it.

That is the product claim working exactly as designed: **the documentation prevented the error
rather than explaining it afterwards.**

## What the existing tools won on — and this is the important part

**On two questions the existing-tool arm was right and mzLibRust was wrong**, because mzLib itself
is wrong:

1. **Modifications: 36 usable, not 14.** The rustyms arm joined UniProt's `ptmlist.txt` on `CF`/`MM`
   and found that 22 of albumin's 24 "glycosylation site" annotations are
   `N-linked (Glc) (glycation) lysine` — formula `C6H10O5`, mass `162.052823`, in both PSI-MOD and
   Unimod. Entirely usable. mzLib drops them because `ProteinXmlEntry.ParseFeatureEndElement` only
   handles `modified residue` and `lipid moiety-binding region`. Filed as
   [mzLib#1112](https://github.com/smith-chem-wisc/mzLib/issues/1112).

   **Worse: our census told the user the wrong reason.** `explain()` said the excluded sites *"have
   no defined chemical composition, so no mass can be assigned"* — false for 22 of 24. We built a
   feature whose entire purpose is refusing to let a caller be silently wrong, and the disclosure
   itself was the misinformation. Corrected, with a test pinning the correct reason.

2. **PRIDE holds 13 files, not 8.** The hand-rolled arm cross-checked the API against the FTP tree
   and found five files missing from the API — including `…01-20141210.mzML` (450 MB) and the
   matching `.mzXML` (472 MB), the modern open-format conversions most users want. True project size
   is **1.44 GB**, not the 0.514 GB the API reports. This is PRIDE's omission, not mzLib's — mzLib
   faithfully reports the v3 API — but "the complete file manifest" was our claim, and it was wrong
   by 65% of the bytes. Now documented on `list_files` and `total_size_bytes`.

The lesson generalises: **a binding that faithfully reports a single upstream source inherits that
source's blind spots and presents them with unearned authority.** The hand-rolled arm was *forced*
to look at two sources and therefore saw the discrepancy. We were not, and did not.

## The state of Rust proteomics

- **No PRIDE client exists.** Searching `pride` returns LGBTQ+ flag utilities; `proteomexchange`
  returns zero crates.
- **`cargo add sage-core` installs an LLM agent framework.** The name on crates.io belongs to an
  unrelated "Core library for Sage Agent". The proteomics Sage is git-only. A silent misroute on the
  single most important crate for the quant task.
- **`rustyms` is genuinely excellent** at what it covers — ProForma 2.0, literature-cited
  fragmentation models with DOIs in the doc comments, five bundled ontologies. It does peptide
  chemistry better than mzLib does in places.
- **But nothing connects the chemistry to the databases.** rustyms bundles Unimod, PSI-MOD, RESID,
  XLMOD and 186,811 GNOme entries, and **0 of 55** of albumin's UniProt annotations could be looked
  up by the names UniProt uses. `find_name("Phosphoserine")` → `None` across every ontology; Unimod
  calls it `Phospho`, PSI-MOD calls it `O-phospho-L-serine`. That join is the gap, and it is a gap
  in the ecosystem, not in any one crate.
- **No LFQ quant exists at all.** Closing the quant gap properly was estimated at 1,100–1,700 lines
  dominated by MBR FDR control and protein rollup — with the arm's own MBR numbers correlating at
  **r = −0.13** between technical replicates (versus r = 0.75 for MS2-anchored peptides), i.e. mostly
  noise without the FDR machinery that no Rust crate provides.

So the value proposition holds, and is sharper than "pyteomics vs pyMzLib" was: for quant there is
no Rust alternative at all, and for peptidoforms the alternative has better chemistry but cannot
reach the annotations.

## Defects filed

Six issues, all found through this exercise:

| Issue | Repo | What |
|---|---|---|
| [#1109](https://github.com/smith-chem-wisc/mzLib/issues/1109) | mzLib | ETD/ECD emit `y` ions with no `b` ions — no fragmentation mechanism |
| [#1110](https://github.com/smith-chem-wisc/mzLib/issues/1110) | mzLib | z• suppressed N-terminal to proline, complementary c ions still emitted |
| [#1111](https://github.com/smith-chem-wisc/mzLib/issues/1111) | mzLib | FlashLFQ peptide roll-up nondeterministically drops MBR intensities |
| [#1112](https://github.com/smith-chem-wisc/mzLib/issues/1112) | mzLib | Glycosylation-site annotations dropped despite defined formula and mass |
| [#7](https://github.com/smith-chem-wisc/pyMzLib/issues/7) | pyMzLib | `Peptide.intensity()` returns `None` against its documented invariant |
| [#8](https://github.com/smith-chem-wisc/pyMzLib/issues/8) | pyMzLib | `--no-modifications` also discards proteolysis products |

One claim was **investigated and dropped**: the extra `zDot` numbered `L` looked like an off-by-one,
and is in fact the deliberate, correct N–Cα cleavage at residue 1 (`M − NH₂`), clearly commented in
mzLib. Filing it would have been wrong. Recorded because "verify before filing" is the practice, not
the exception.

## Documentation changes, and where they must be back-ported

Every item below applies to **pyMzLib** verbatim, and later to mzLibR — the traps live in mzLib's
behaviour, not in any one language.

| Fix | Why |
|---|---|
| `total_size_bytes`: PRIDE reports *decompressed* size for some `.gz` files | The MGF reports 16,448,103 and downloads 5,984,662 — 2.75×. An agent's first hypothesis was a truncated download, which would have been a wrong bug report. |
| `list_files`: the API manifest is incomplete | 8 of 13 for PXD000001. |
| `DownloadOptions::extensions`: cross-reference the `.gz` trap | The warning lived on `extension()`; the person about to err is reading `extensions`, whose examples were both uncompressed. **Put the warning where the mistake is made, not where the concept is defined.** |
| `FragmentOptions::dissociation`: ETD returns c, **y**, zDot | The doc said "c and z• ions". An agent nearly published "16,066 c and z• ions"; 5,338 were y. |
| `Fragment::fragment_number`: zDot runs `1..=L`; proline asymmetry | Neither was discoverable from the API. |
| `FragmentOptions::modifications`: `false` is not a clean control | It also disables signal-peptide processing. |
| `ModificationCensus::explain()`: state the real exclusion reason | The old reason was false for 22 of 24. |
| `QuantifyOptions::max_threads`: pin to 1 for reproducibility | Presented as a performance knob; it changes results. |
| `mbr_rescued_peptide_count`: state the exact predicate | Prose and implementation diverged, 140 vs 135. |
| `Peptide::intensities`: the two surfaces aggregate differently | The docs tell you to switch surfaces without saying they disagree on magnitude. |
| `ProteinGroup::intensities`: `Some(0.0)` is the common case | 847 of 943 versus 2 `None`. `Option` marks the rare case, not the one you probably mean. |
| `protease`: "37 peptides out of about 200" → **7** | The number was stale by 5×, and it is the kind of figure a reader quotes. |
| Two new methods: `fragments_by_series()`, `distinct_base_sequences()` | `fragment_count()` invited a wrong answer; `peptides.len()` counts peptidoforms, not sequences (303 vs 195). |

## Methodological notes

- **Verification dominated.** Every arm wrote far more code checking its answer than producing it —
  quant 375 of 568 lines, peptidoform 377 of 599. The agents that got things right were the ones
  that distrusted a plausible number.
- **The `#[cfg(test)]` prohibition mattered.** Agents reported navigating around test modules
  deliberately. Without it they would have read the answers.
- **Ground truth had to be independent.** Derived through the bridge, not the crate. Two of the nine
  answers turned out to be wrong *in the ground truth's own source* (mzLib), which only became
  visible because a second toolchain disagreed.
- **What we did not test:** a human user, or a user who does not verify. Both agents that avoided
  the big traps did so partly because the docs warned them — but they were also unusually careful.
  A hurried user reading `fragment_count()` still gets 16,066.
