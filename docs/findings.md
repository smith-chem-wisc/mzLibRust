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
