# Findings

Defects and documentation gaps surfaced while building mzLibRust.

Two standing rules govern what happens to anything on this page:

1. **A bug in mzLib is filed as an issue on mzLib.** A workaround in the bridge or in a binding is
   not the fix — mzLib bugs need to be fixed in mzLib, or MetaMorpheus and every other consumer
   stays broken and the next binding rediscovers it from scratch.
2. **A documentation lesson is back-ported to every binding.** The traps live in mzLib's behaviour,
   not in any one language, so wording that rescues a Rust user rescues a Python (and later an R)
   user too.

---

## 1. ETD and ECD generate y ions but not b ions — mzLib

**Status:** filed as [smith-chem-wisc/mzLib#1109](https://github.com/smith-chem-wisc/mzLib/issues/1109).
**Found by:** a parity test asserting "ETD produces c and z• ions, not b and y", which failed against
the recorded albumin digest.

`ProductsFromDissociationType` in `Omics/Fragmentation/Peptide/DissociationTypeCollection.cs` maps:

```csharp
{ DissociationType.ECD,   new List<ProductType>{ ProductType.c, ProductType.y, ProductType.zDot } },
{ DissociationType.ETD,   new List<ProductType>{ ProductType.c, ProductType.y, ProductType.zDot } },
{ DissociationType.EThcD, new List<ProductType>{ ProductType.b, ProductType.y, ProductType.c, ProductType.zDot } },
```

ETD and ECD cleave the **N–Cα** bond, giving c and z• ions; b and y come from **amide** cleavage
under vibrational activation. If the y ions modelled residual activation then b ions would accompany
them — which is exactly what the `EThcD` row does, correctly. **y without b has no mechanism.** The
rows also read like an edit of a `{ b, y }` row in which `b` became `c, zDot` and `y` was left
behind.

**Impact:** roughly a third of every ETD/ECD theoretical fragment list is ions ETD does not make. On
albumin's 31-mer tryptic peptide: 30 c, 30 z•, **30 y**. In a search these are matchable theoretical
ions, so they can inflate matched-ion counts and scores.

**Handling here:** `etd_produces_c_and_z_ions_and_also_y_which_it_should_not` and
`the_scale_of_the_spurious_etd_y_ions_is_recorded` assert what mzLib *currently does*, so they fail
when the upstream fix lands rather than making the fix look like a regression.

**Back-port:** pyMzLib's `peptidoform` docs describe ETD as "c and z• ions" without mentioning the y
ions its users are actually receiving. Worth a sentence there too, until #1109 is fixed.

---

## 2. A `null` peptide intensity contradicts pyMzLib's documented invariant — pyMzLib

**Status:** not yet filed. **Found by:** porting `flashlfq_small.json`, whose own contents contradict
the docs it ships with.

`pkg/python/src/pymzlib/flashlfq.py` documents peptide intensities as:

> Missing is `0.0`, **never `None`** (unlike proteins).

and `Peptide.intensity()` promises the same. But `_from_wire` stores the wire value verbatim
(`intensities=dict(payload.get("intensities") or {})`), and the shipped fixture contains:

```json
"sequence": "AC[Carbamidomethyl]DEFR",
"intensities": { "run_3": 2000.0, "run_4": null }
```

so `peptide.intensity("run_4")` returns `None` — the documented invariant is not enforced. The
Python test only covers the *missing-key* case (`intensity("never_seen") == 0.0`), never the
*null-value* case its own fixture provides.

Why it matters: the None-vs-0.0 distinction is the headline safety property of the FlashLFQ tranche
— `None` is supposed to mean "FlashLFQ could not resolve this", and that is supposed to be a
*protein-only* condition. A peptide that also returns `None` erodes exactly the signal a caller is
being told to rely on, and code branching on `is None` to detect unresolvable proteins will
misclassify peptides.

**Handling here:** `Peptide::intensities` is `HashMap<String, f64>` and the null is resolved to `0.0`
at the deserialization boundary, so the documented invariant is true by construction. Only
`ProteinGroup::intensities` is `Option<f64>`. `a_null_peptide_intensity_reads_as_zero_not_as_unresolvable`
pins it.

**Back-port:** pyMzLib should either coerce `None → 0.0` in `Peptide._from_wire` (making the docs
true, and matching this crate) or change the documentation and the type hint. The first is almost
certainly right — FlashLFQ's own `GetIntensity` returns 0 for a peptide it did not measure.

---

## 3. PRIDE reports decompressed sizes for some compressed files — documentation, both bindings

**Status:** documented in this crate; **needs back-porting to pyMzLib**.
**Found by:** the PRIDE bake-off arm, which saw a 5.98 MB file land where the manifest promised
16.4 MB and correctly suspected a truncated download before checking.

For PXD000001, `file_size_bytes` versus what actually arrives:

| File | PRIDE reports | bytes written | `gzip -l` uncompressed |
|---|---|---|---|
| `…pride.mgf.gz` | 16,448,103 | 5,984,662 | **16,448,103** |
| `…pride.mztab.gz` | 497,985 | 103,845 | **497,985** |
| `…xml.gz` | 10,677,205 | 10,668,000 | 48,038,607 |

For two of the three, PRIDE's reported size is *exactly* the decompressed length. All three pass
`gzip -t`, so these are complete downloads, not truncations. PRIDE's metadata is simply inconsistent,
and neither mzLib nor a binding can correct it.

Why it matters: `total_size_bytes()` is documented as the way to see what a download will cost
before starting it, and for the file a user most likely wants here it overstates the transfer by
**2.75×**. The near-miss is the real point — the bake-off agent's first hypothesis was a truncated
download, which would have been a serious and wrong bug report.

**Handling here:** [`total_size_bytes`](../src/pride.rs) now carries the caveat with these numbers,
and is described as an upper bound rather than a cost estimate.

**Back-port:** pyMzLib's `pride.total_size_bytes()` and `PrideFile.size_mb` have the identical
docstring promise and the identical problem.

## 4. The `.gz` trap is prevented by documentation, and it works — evidence

Not a defect; recorded because it is the clearest evidence the "disclose the traps" convention pays
for itself, and because it says where the disclosure has to *live*.

The bake-off agent was asked to download "only the compressed MGF". The natural filter,
`extensions: [".mgf"]`, matches **zero** files — PXD000001's peak list is
`PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz`, extension `.gz`. The agent did not hit it, and said why:
`PrideFile::extension`'s doc-comment volunteers the `x.mgf.gz` example unprompted.

It then ran the wrong filter deliberately, and got:

> `No file in PXD000001 matched extensions [".mgf"]. Use list_files() to see what the project actually contains — note that compressed files such as 'x.mgf.gz' have the extension '.gz'.`

An error rather than an empty success, **with the fix inside the error**, and no stray empty
destination directory left behind.

The lesson for documentation placement: the warning was on `extension()`, but *the person about to
make the mistake is reading `DownloadOptions::extensions`* — whose examples were `.raw` and `.mzML`,
both uncompressed, warning nobody. The trap is now cross-referenced from the field that causes it.
That generalises: **put the warning where the mistake is made, not where the concept is defined.**
Worth applying to pyMzLib's `download(extensions=…)` docstring too.

## 5. Glycosylation-site annotations are dropped silently — and the exclusion is *right*

**Status:** [mzLib#1112](https://github.com/smith-chem-wisc/mzLib/issues/1112), with a correction
posted to the issue.
**Found by:** a bake-off arm that resolved 22 of albumin's 24 excluded annotations to a defined mass
and concluded they were usable. **That conclusion was wrong, and so was our first report of it.**

`ProteinXmlEntry.ParseFeatureEndElement` turns a UniProt feature into a modification only for
`modified residue` and `lipid moiety-binding region`. `glycosylation site` is discarded before any
lookup, so it never reaches `unknownModifications` and nothing reports it.

The first version of this finding argued mzLib should load them, because UniProt's `ptmlist.txt`
gives `N-linked (Glc) (glycation) lysine` a formula (`C6H10O5`) and a mass (`162.052823`).
**A mass spectrometrist corrected that, and they were right:**

- **Glycation is ambiguous in exactly the way the exclusion exists to catch.** The Amadori product is
  labile and heterogeneous, progresses to advanced glycation end products, and dissociates in
  preference to the peptide backbone. Emitting one exact mass and a clean c/z ladder would describe a
  species you cannot observe. The presence of `CF`/`MM` in `ptmlist.txt` is **not** sufficient
  evidence that a modification is usable. The same argument applies to `lipid moiety-binding region`,
  which mzLib *does* load.
- **Most of these are not modifications of the real protein anyway.** By qualifier: **14** of the 22
  glycation annotations are `; in vitro`, 1 is `; alternate`, and both `N-linked (GlcNAc...)` sites
  exist only in the Redhill and Casebrook variants. mzLib reads none of those qualifiers, for **any**
  feature type — including the `modified residue` ones it does load.

**So the defect is the silence, not the exclusion.** A caller gets 14 with no indication that 24 more
were annotated, that they were dropped on feature type rather than chemistry, or that most were in
vitro. The revised ask upstream is to *report* the excluded annotations — type, description and
position — not to load them.

**Lesson for this crate:** our `explain()` originally said the excluded sites "have no defined
chemical composition", which was false; we then corrected it to imply they should have been loaded,
which was also false. Both times the trap-disclosure was itself the misinformation. It now states
what actually happened (dropped on feature type), why that is usually correct, and that the census
can only give a count — so the reader must go and read the annotations. See finding 6.

## 6. `ModificationCensus` reports a tally, not the annotations

**Status:** open; needs a bridge change (`CensusAnnotatedFeatures` in `pkg/bridge/Peptidoform.cs`).
**Found by:** the biologist bake-off arm, whose *single* external lookup in the whole task was
leaving the crate to read the excluded annotation strings from UniProt.

`by_type` gives a feature-type name, a count, and a loaded flag. It cannot distinguish "no defined
composition" from "in vitro" from "variant-only" — three exclusions needing three different
decisions. The user is told *that* something was excluded but has no way to check *what*, so they
cannot evaluate the reason and must leave the tool.

This is the same failure mode as findings 4 and 5 in miniature: disclosure that stops one level short
of actionable. The fix is to carry the description string and position through the census, which
requires the bridge to emit them — so both bindings gain it at once.

## Already known, carried over from pyMzLib

These were found during the pyMzLib build and are relevant to anyone reading this crate's source:

- **[mzLib#1103](https://github.com/smith-chem-wisc/mzLib/issues/1103)** — TorchSharp/libtorch is
  ~238 MB of the bridge payload, dragged in transitively. Making it optional shrinks *both*
  bindings, and matters more here because the Rust distribution story is download-at-build.
- **[mzLib#1106](https://github.com/smith-chem-wisc/mzLib/issues/1106)** — `trypsin` and `trypsin|P`
  mean the reverse of the MaxQuant/Mascot convention. Documented on
  [`FragmentOptions::protease`](../src/peptidoform.rs).
- **[mzLib#1108](https://github.com/smith-chem-wisc/mzLib/issues/1108)** — `Protein.Digest` returns
  duplicate peptidoforms where a chain boundary meets the initiator-Met cleavage site. Worked around
  *in the bridge*, so both bindings inherit the fix; the issue tracks the real repair.
