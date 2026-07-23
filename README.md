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

## Three capabilities, the same three pyMzLib has

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

The crate needs a staged bridge executable. Resolution order:

1. `MZLIB_BRIDGE` — a path to a bridge binary. The dev and offline escape hatch.
2. `_dotnet/<runtime-identifier>/mzlib-bridge[.exe]` beside the crate.

Build one with pyMzLib's `pkg/build/publish-bridge.ps1`, which cross-publishes any .NET runtime
identifier. Download-at-build via `build.rs` lands once pyMzLib's CI publishes the raw bridge
binaries as release assets.

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
