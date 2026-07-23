# Bake-off results — biologist phenotype

Six agents, three tasks × two toolchains, run under the persona and instrumentation pinned in
[DESIGN.md](DESIGN.md): ~8 years at the bench and on an Orbitrap, knows the science, **not an
experienced coder, leans on docs, does not read library source.**

Ground truth was established by driving the bridge executable directly from the shell, so the crate
could not score itself.

*(A first attempt used an expert-programmer persona and had to be discarded — see
[the appendix](#appendix-the-invalidated-expert-run).)*

---

## Headline

| | Arm A — existing Rust ecosystem | Arm B — mzLibRust |
|---|---|---|
| **external lookups** | **13** | **2** |
| **dead ends** | **8** | **1** |
| **answers the biologist would publish** | 4 of 11 | 10 of 11 |
| **finished** | 3 of 3 — by hand-writing the tool | 3 of 3 |

The lookup ratio is the cleanest signal, and it is **6.5×**.

---

## Task 1 — label-free quant with match-between-runs

| | Arm A | Arm B | Truth |
|---|---|---|---|
| peptides quantified in both runs | **284** ✗ | **257** ✓ | 257 |
| peptides rescued by MBR | **173** ✗ | **140 / 135** ✓ | 140 / 135 |
| protein groups with no intensity | **8** ✗ | **2** (+847 zeros) ✓ | 2 |
| external lookups / dead ends | 3 / 2 | **0 / 0** | |
| would publish? | counts cautiously; **intensities no** | yes | |

Arm A is wrong in the dangerous direction — 284 against 257, 173 against 140. Plausible numbers. No
reviewer would blink.

The reason was not the language:

> The stuck part wasn't Rust. I braced for the borrow checker and it barely bit me. The stuck part
> was realizing an hour in that **there is no tool here**… I now have to be the person who decides
> what "10 ppm" and "±1 minute" mean, and I have no way to tell if I chose wrong.
>
> Every one of those changes the three numbers I just reported and **I can't defend any of them from
> documentation.**

They hand-wrote ~250 lines of quantifier — XIC construction, apex finding, integration, RT transfer,
protein rollup — against no reference implementation. To their credit they invented a control
(shift every peptide mass by +11.5 Da and re-run: 173 "transfers" collapsed to 9), which is good
instinct. Their own verdict: *"I finished by writing my own quantifier, not by using a tool."*

Arm B, by contrast, reported **zero** external lookups and never got stuck:

> The place I *should* have gotten stuck was Q1. My instinct as a biologist is "peptide table →
> count non-zero cells", and that instinct is wrong here in a way that produces a plausible-looking
> number (169). **The doc comment on that field says it in plain English, in bold, and points me at
> `peaks`.**

They also set `max_threads: 1` because the docs warned them results are otherwise nondeterministic —
a warning that only exists because the earlier expert run found it. That loop closing is the single
clearest demonstration that these bake-offs pay for themselves.

## Task 2 — digest and fragment albumin

| | Arm A | Arm B | Truth |
|---|---|---|---|
| tryptic peptides | 48 / 208 / 249 / 261 — **wouldn't publish** | **195 seq / 303 peptidoforms** ✓ | 195 / 303 |
| ETD fragments | 22,882 / 5,850 / 924 — **wouldn't publish** | **10,174 real** (16,066 raw) ✓ | — |
| searchable modifications | 31 of 55 → "53 if you patch it" ✗ | 14 loaded, flagged as incomplete | **exclusion correct; 14 usable** |
| pSer89 peptide at 2+ | **549.25499** ✓ | **549.255** ✓ | agree |
| external lookups / dead ends | 7 / 3 | **1 / 0** | |
| would publish? | 2 of 4 | 4 of 4 | |

**Both arms independently landed on 549.255 for SLHTLFGDK + phospho-Ser89.** Two entirely different
chemistry stacks agreeing to 5 decimal places is a good sign for both.

**Both arms independently concluded the 22 glycation sites are searchable** (+162.0528) — and
**both were wrong, as was this report's first version.** A mass spectrometrist reviewing it pointed
out that glycation is ambiguous in precisely the way the exclusion exists to catch: the Amadori
product is labile and heterogeneous, progresses to AGEs, and dissociates in preference to the
backbone, so an exact mass plus a clean fragment ladder describes an unobservable species. The
annotation qualifiers agree — **14 of the 22 are `; in vitro`** and both `GlcNAc` sites exist only in
disease variants.

Two independent arms converging on the same answer was not corroboration; it was two arms making the
same reasonable inference from `ptmlist.txt` having a `CF`/`MM` line. Agreement between arms is
evidence, not proof, and this is the case that shows the difference.

The defect is the **silence**, not the exclusion —
[#1112](https://github.com/smith-chem-wisc/mzLib/issues/1112) has been corrected accordingly.

Arm A's most alarming finding is about the *ecosystem*, not about us. The digestion crate they found
(69 downloads, no documentation) **cleaves before proline and never says so**:

> I only caught the proline thing because I've stared at enough albumin peptide lists to know
> `LVRPEVDVM…` doesn't look right. Someone with less bench time would have published a wrong number
> and never known.

And on the modification vocabulary gap:

> **Left to the documented API I would have concluded that zero of albumin's 55 modifications have a
> mass, which is spectacularly wrong.**

`Ontology::find_name("Phosphoserine")` returns `None` in every bundled ontology; the fuzzy matcher
confidently returns **"homoserine"**. The correct synonym is sitting in the PSI-MOD record, but no
documented API reaches it.

## Task 3 — PRIDE manifest and filtered download

| | Arm A | Arm B | Truth |
|---|---|---|---|
| file count | 8 — unaware of the other 5 | 8, **flagged incomplete → 13** ✓ | 13 |
| download size | **0.503 GB** ✗ | **1.44 GB** ✓ | 1.44 GB |
| the MGF | ✓ | ✓ | |
| typo `PXD0000019999` | caught — via a guard added *after being burned* | hard error | |
| external lookups / dead ends | 3 / 3 | **1 / 1** | |

Arm B got the size right **only because of a doc fix made earlier the same day** from the expert
run's finding.

Arm A's summary of the raw API is the case for the library existing:

> The API lies quietly in three different ways and none of them throw an error: a nonexistent
> accession looks like an empty dataset, a size field means something different than you'd assume,
> and an array's order changes between records. Every one of those produces a program that runs,
> exits 0, and prints a confident wrong answer. I caught all three, but I caught the size one **by
> accident** and the PEAK-filter one **by paranoia**, not by any process.

Two traps mzLib absorbs silently:

- **`fileCategory == "PEAK"` matches 2 files, not 1** — the MGF *and* a 243 MB mzXML. Taking the
  first would have pulled **40× the intended download**.
- **`publicFileLocations` order is not stable.** The mztab lists FTP first; the MGF lists **Aspera**
  first. Indexing `[0]` yields `prd_ascp@fasp.ebi.ac.uk:…`, unfetchable by any HTTP client. mzLib's
  `TryGetHttpsDownloadUrl` searches rather than indexes, so Arm B never saw this.

### The one place Arm B genuinely hit the trap

Arm B wrote the natural filter, `extensions: vec![".mgf"]`, and got **zero files** — the file is
`…pride.mgf.gz`, extension `.gz`. The library refused rather than reporting success:

> **How I noticed: the library refused.** It didn't create an empty folder and report success… In a
> Python script I'd have gotten an empty list, shrugged, and found out three weeks later.

That is better evidence than the expert run produced. A trap avoided proves the docs are readable; a
trap **hit and caught** proves the guard works.

---

## What this says to do next

**1. Warnings must live where the mistake is made.** Confirmed three times now — `Peptide::intensities`
(caught the 169 error), `PrideFile::extension` (caught the `.gz` error), and the counter-example
where `DownloadOptions::extensions` had no warning and the biologist walked straight in. Apply this
systematically to pyMzLib rather than case by case.

**2. `total_size_bytes` is defended only by prose, and prose is thin.** Arm B got it right but said
they'd have been wrong by 3× if they'd skimmed: *"a bar chart doesn't carry a footnote."* Worth
considering a name or return type that carries the caveat to the call site — against the cost of
breaking the follow-mzLib-names convention.

**3. `ModificationCensus` should expose the excluded annotations, not just a tally.** Arm B's single
external lookup was leaving the tool to read the 24 excluded annotation strings from UniProt, which
is how they learned 22 were glycation. The census reports *that* something was excluded but not
*what*, so the user cannot check the reason. Needs a bridge change.

**4. There is no FTP cross-check.** Arm B had to `curl` the directory listing to answer "how much
disk space do I need" — the second question anyone asks.

---

## Defects filed

Seven, all found through these exercises and all verified against ground truth before filing:

| Issue | Repo | What |
|---|---|---|
| [#1109](https://github.com/smith-chem-wisc/mzLib/issues/1109) | mzLib | ETD/ECD emit `y` ions with no `b` ions — no fragmentation mechanism |
| [#1110](https://github.com/smith-chem-wisc/mzLib/issues/1110) | mzLib | z• suppressed N-terminal to proline, complementary c ions still emitted |
| [#1111](https://github.com/smith-chem-wisc/mzLib/issues/1111) | mzLib | FlashLFQ roll-up nondeterministically drops MBR intensities |
| [#1112](https://github.com/smith-chem-wisc/mzLib/issues/1112) | mzLib | Glycosylation annotations dropped despite defined formula and mass |
| [#1113](https://github.com/smith-chem-wisc/mzLib/issues/1113) | mzLib | Modifications applied after digestion → peptides ending in a modified K trypsin cannot cleave |
| [#7](https://github.com/smith-chem-wisc/pyMzLib/issues/7) | pyMzLib | `Peptide.intensity()` returns `None` against its documented invariant |
| [#8](https://github.com/smith-chem-wisc/pyMzLib/issues/8) | pyMzLib | `--no-modifications` also discards proteolysis products |

One claim was **investigated and dropped**: the extra `z•` numbered `L` is the deliberate, correct
N–Cα cleavage at residue 1, clearly commented upstream. Filing it would have been wrong.

---

## Appendix: the invalidated expert run

The first attempt used *"a proteomics researcher who writes some Rust."* Those agents read library
source, reverse-engineered a dependency's internal types, and wrote 375–598 lines of verification
scaffolding each. It also omitted the `external_lookup` instrumentation entirely.

Its conclusion — "nobody was confidently wrong" — was an artifact of the persona. A careful expert
self-corrects regardless of documentation quality, so the comparison could not have come out any
other way. It measured *can an expert reach the right answer* rather than *do the docs carry a
biologist through*.

**It was not wasted.** Every defect above except #1113 came from it, and two documentation fixes it
prompted (the `max_threads` determinism warning and the PRIDE manifest incompleteness) are directly
why Arm B answered two questions correctly in the valid run. Expert arms are good **defect
scanners**; biologist arms are the only valid **comparison**. Run both, in that order.
