# mzLibRust — mzLib for Rust

[mzLib](https://github.com/smith-chem-wisc/mzLib) is a mass-spectrometry and proteomics library
written in C#. **mzLibRust** makes its functionality callable from Rust, with no .NET installation.

It is the sibling of [pyMzLib](https://github.com/smith-chem-wisc/pyMzLib) and speaks the same
**language-neutral bridge**: a self-contained executable exchanging a versioned JSON envelope over
stdin/stdout, which assumes nothing about the language calling it. Everything genuinely hard already
lives there — the mzLib interop, the composition of mzLib's own methods, and the
availability-versus-correctness error classification. This crate is the thin, idiomatic Rust surface
over it, which is why it exists at all: a second binding costs a transport module and some typed
structs, not a second implementation of mzLib.

## Four capabilities, the same four pyMzLib has

```rust
// PRIDE Archive — what is in a project, and pull it down.
let files = mzlib::pride::list_files("PXD000001")?;
let small: Vec<_> = files.iter()
    .filter(|f| f.size_mb() < 5.0 && f.downloadable())
    .cloned()
    .collect();
mzlib::pride::download_files(&small, "downloads", &Default::default())?;

// Peptidoforms — digest an annotated protein and fragment its peptides.
let digest = mzlib::peptidoform::fragments("P02768")?;
println!("{}", digest.modification_census.explain());
//   14 of 38 annotated modifications were applied, across 14 residue positions.
//   Excluded by type: 24 × glycosylation site — these have no defined chemical composition…

// FlashLFQ — label-free quantification across runs.
use mzlib::flashlfq::{quantify_with, QuantifyOptions, SpectraFile};
let result = quantify_with(
    "AllPSMs.psmtsv",
    &[SpectraFile::from("run_3.mzML"), SpectraFile::from("run_4.mzML")],
    &QuantifyOptions { match_between_runs: true, ..Default::default() },
)?;
println!("{} peptides rescued by MBR", result.mbr_rescued_peptide_count());

// Readers — identify any of mzLib's 29 result-file types, and read ALL of them.
let info = mzlib::readers::identify("psm.tsv")?;
println!("{} {:?}", info.file_type, info.views);   // MsFraggerPsm ["quantifiable"]

let table = mzlib::readers::read_records("toppic_prsm.tsv")?;   // works on all 29
let e_values = table.columns.floats("e_value")?;                // Vec<Option<f64>>
```

### Reading: one universal function, four typed views

What differs between the 29 formats is not *whether* you can read them but *what the columns mean*.

| function | reads | columns |
|---|---|---|
| `read_records` | **all 29** | **that format's own fields**, under mzLib's names |
| `read_results` | 3 | uniform `quantifiable` view |
| `read_features` | 2 | uniform `ms1_features` view |
| `read_matches` | 4 | uniform `spectral_match` view |
| `read_spectra` | 7 | scan headers; peaks opt-in |

**13 of the 29 belong to no cross-format family at all** — TopPIC, Crux, MSFragger's peptide and
protein tables, the FlashDeconv formats. mzLib parses them into a format-specific shape and there
is no uniform view to project them onto, so `read_records` is what reaches them; it is a necessity,
not a convenience.

Because the column set depends on the format, a read returns a `Table` rather than a struct with
named fields — with typed accessors that project a wire `null` onto `Option`, so a missing cell can
never silently become a zero and can never shorten a column:

```rust
let t = mzlib::readers::read_records("crux.txt")?;
for (sequence, score) in t.columns.strings("base_sequence")?
    .iter()
    .zip(t.columns.floats("x_corr_score")?)
{
    if let (Some(sequence), Some(score)) = (sequence, score) { println!("{sequence}\t{score}"); }
}
```

## Two conventions worth knowing

**Names follow mzLib.** A field here means exactly what it means in the mzLib source, the
MetaMorpheus output columns, and the papers — `match_between_runs`, `ppm_tolerance`,
`protein_groups`, `detection_type`. Nothing is renamed to look more Rust-like, because renaming
forces every reader to hold a translation table in their head, and makes "what is `ppm_tolerance`?"
answerable straight from FlashLFQ's own docs.

**The types disclose the traps.** Where mzLib applies an invisible rule, this crate surfaces it
rather than swallowing it. The clearest case is quantification:

```rust
struct Peptide      { intensities: HashMap<String, f64>,         /* 0.0 = not measured here */ }
struct ProteinGroup { intensities: HashMap<String, Option<f64>>, /* None = could not be resolved */ }
```

FlashLFQ's median-polish protein quant emits NaN when a protein's peptide matrix is degenerate. In
Python that distinction lives in the documentation and has to be remembered. Here the compiler makes
you handle it — and a `0.0` peptide can never be mistaken for an unresolvable one.

The same doctrine runs through the rest: [`Digest::truncated`] tells you the silent isoform cap
bound, [`ModificationCensus::explain`] tells you what UniProt annotated versus what could be used,
and `FlashLfqResults::peaks` is documented as the surface to read for match-between-runs because the
peptide roll-up drops most transfers.

## Getting the bridge

**You do not need one to build, test, document or lint the crate.** `cargo test` runs 118 offline
tests with no network and no .NET on the machine. A bridge is required only for calls that actually
reach mzLib, and a missing one is a runtime error with instructions, never a build failure — so
contributors are never blocked by a 130 MB payload they may not want.

`build.rs` resolves it in this order:

1. **`MZLIB_BRIDGE`** — a path to a bridge you already have. Always wins, checked at runtime too.
2. **`_dotnet/<runtime-identifier>/mzlib-bridge[.exe]`** beside the crate.
3. **Downloaded** from `MZLIB_BRIDGE_URL`, verified against `MZLIB_BRIDGE_SHA256` when you set it.

The quickest route, if you have a pyMzLib checkout — it already stages a bridge for its wheel:

```powershell
.\scripts\stage-bridge.ps1 -PyMzLibRoot ..\pyMzLib          # copy the one pyMzLib staged
.\scripts\stage-bridge.ps1 -PyMzLibRoot ..\pyMzLib -Build   # or build a fresh one (needs .NET)
```

The script probes the staged binary by asking it for its version, because a payload that cannot
report that will certainly fail when the crate calls it.

Or simply:

```bash
export MZLIB_BRIDGE=/path/to/mzlib-bridge
cargo test --features live
```

`MZLIB_BRIDGE` is also a **licence affordance**, not only a convenience: it is how you point this
crate at a bridge built from a modified mzLib, exercising your LGPL §4 right to relink without
rebuilding anything here. See `NOTICE`.

Download-at-build becomes the default path once pyMzLib's CI publishes the raw bridge binaries as
release assets; the machinery for it is already in `build.rs`.

## Testing

```bash
cargo test          # 118 offline tests: no network, no bridge needed
MZLIB_BRIDGE=… cargo test --features live    # 18 live canaries; they SKIP on an outage
```

The offline suite is the default because it must pass anywhere, in milliseconds. Live canaries are
the ones that would catch mzLib, PRIDE or UniProt changing under us, and they **skip rather than
fail** when a service is down — an ambiguous red build gets ignored, which is how a genuine contract
break goes unnoticed for a month.

See [docs/test-parity.md](docs/test-parity.md) for the test-by-test mapping against pyMzLib, and
[docs/findings.md](docs/findings.md) for defects this port surfaced upstream.

## Licence

Same as mzLib. See `LICENSE`.
