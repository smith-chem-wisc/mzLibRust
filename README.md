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

## What it does

<!-- This heading used to read "the same three pyMzLib has", then "Four capabilities" over five
     modules. A parity claim is a claim about someone else's repository, and a count is a claim that
     goes stale the day a module lands, so neither is made here: state what this crate does and let
     the reader compare. -->

```rust
// PRIDE Archive — find projects by keyword, see what is in one, and pull it down.
let hits = mzlib::pride::search("plasmodium falciparum schizont")?;
println!("{} {:?}", hits[0].accession, hits[0].matched_fields());   // PXD070842 ["references", "title"]
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
//   Excluded by type: 24 × glycosylation site — mzLib loads only 'modified residue' and …

// FlashLFQ — label-free quantification across runs.
use mzlib::flashlfq::{quantify_with, QuantifyOptions, SpectraFile};
let result = quantify_with(
    "AllPSMs.psmtsv",
    &[SpectraFile::from("run_3.mzML"), SpectraFile::from("run_4.mzML")],
    &QuantifyOptions { match_between_runs: true, ..Default::default() },
)?;
println!("{} peptides rescued by MBR", result.mbr_rescued_peptide_count());

// Median polish — re-roll proteins under a new design, without re-reading any spectra.
use mzlib::flashlfq::{median_polish_with, DesignEntry, MedianPolishOptions};
let proteins = median_polish_with("QuantifiedPeptides.tsv", &MedianPolishOptions {
    design: vec![DesignEntry::new("run_3").condition("control"),
                 DesignEntry::new("run_4").condition("treated")],
    ..Default::default()
})?;

// Readers — spectra as well as search output; identify any of mzLib's 36 types and read ALL of them.
// mzML, Thermo .raw, Bruker .d, timsTOF .d, MGF and msalign all read through read_spectra.
let scans = mzlib::readers::read_spectra("run.mzML")?;
println!("{} scans", scans.scan_count);

let info = mzlib::readers::identify("psm.tsv")?;
println!("{} {:?}", info.file_type, info.views);   // MsFraggerPsm ["quantifiable"]

let table = mzlib::readers::read_records("toppic_prsm.tsv")?;   // works on all 36
let e_values = table.columns.floats("e_value")?;                // Vec<Option<f64>>

// Hundreds of runs, one bridge process, one long table; each run's instrument in its report.
let runs = mzlib::readers::read_spectra_many(&paths, &Default::default())?;
println!("{:?}", runs.files[0].source.as_ref().map(|s| &s.instrument_serial_number));

// MetaMorpheus protein groups and FlashLFQ peptides as long tables: one row per sample.
let groups = mzlib::readers::read_protein_groups("AllQuantifiedProteinGroups.tsv")?;

// SDRF experimental design — what was searched, pooled across experiments with provenance.
let design = mzlib::sdrf::pool_labelled(&[("PXD000070.sdrf.tsv", "malaria"),
                                          ("PXD026824.sdrf.tsv", "colon")])?;
println!("{:?}", design.document.value("characteristics[organism part]"));

// Is it well-formed, do the files agree, does it describe its samples at all? mzLib's answers.
let findings = mzlib::sdrf::validate("PXD000070.sdrf.tsv")?;
let verdicts = mzlib::sdrf::assess_many(&corpus, &Default::default())?;   // Informative / Partial / Skeleton
let ages = mzlib::sdrf::parse_ages(&["58Y", "40Y-85Y", ">=90Y", "63"])?;   // years, or why not

// Protein databases — what an accession is, which gene it is, whether a peptide is unique.
let db = mzlib::proteins::read(&["human.xml"])?;
println!("{:?}", db.taxonomy()?.get("P04406"));                            // Some(Some("9606"))
let calls = mzlib::proteins::classify_peptides(&["YLYEIAR"], &["human.xml", "bovine.fasta"])?;
```

### Reading: one universal function, four typed views

What differs between mzLib 1.0.592's 36 formats is not *whether* you can read them but *what the
columns mean*. The counts are the pin's; `mzlib::readers::formats()` gives the live ones.

| function | reads | columns |
|---|---|---|
| `read_records` | **all 36** | **that format's own fields**, under mzLib's names |
| `read_results` | 4 | uniform `quantifiable` view |
| `read_features` | 2 | uniform `ms1_features` view |
| `read_matches` | 6 | uniform `spectral_match` view |
| `read_spectra` | 7 | scan headers; peaks opt-in |

**17 of the 36 belong to no cross-format family at all** — TopPIC, Crux, MSFragger's peptide and
protein tables, the FlashDeconv formats, SDRF, and the MetaMorpheus and FlashLFQ quantification
tables. mzLib parses them into a format-specific shape and there is no uniform view to project them
onto, so `read_records` is what reaches them; it is a necessity, not a convenience.

**Many files are one call, not a loop.** Every reader has a `_many` twin that takes a list and
returns one long table whose first two columns name the file each row came from, plus one report
per file (its instrument, its absent fields, why it failed under `OnError::Skip`). The list is read
by one bridge process, `threads` files at a time; the answer is identical at any thread count.

**Three quantification tables have functions of their own** (mzLib 1.0.592): MetaMorpheus protein
groups (`read_protein_groups`), FlashLFQ peptides (`read_quantified_peptides`) and PTM site
occupancy (`read_occupancy`), each as a long table, one row per record per sample. `read_records`
reads those files too, but their per-sample values are dictionaries it names and cannot project.

**SDRF is the exception: read it with the `sdrf` module, not `read_records`.** `read_records` joins
each SDRF row into one semicolon-separated string, and SDRF's own `NT=…;AC=…` grammar puts semicolons
inside cells, so the string cannot be split back apart. `mzlib::sdrf::read` and `mzlib::sdrf::pool`
return every cell intact, in a row-major shape that keeps SDRF's repeated column names. The same
module asks mzLib's three questions about a document — `validate` (is it well-formed?),
`lint_labelled` (do several files write the same thing the same way?) and `assess` (does it
describe its samples at all?) — which are blind in different places, which is why there are three.
`samples` lifts each sample's characteristics, and `parse_ages` reads `characteristics[age]` into
years, refusing any cell that would need a guess.

**Protein databases have their own module**, `proteins`: `read` gives one row per protein —
organism, NCBI taxon, genes, mass — with GO terms and Ensembl gene links on request;
`resolve_genes_with` resolves proteins to stable Ensembl gene ids against a gene set you pin;
`classify_peptides` sorts peptides into Unique, SharedWithinGene, SharedAcrossGenes or
NotInDatabase, treating I and L as the same residue. A FASTA's silence about GO and Ensembl is
reported in `absent_fields`, never as an empty answer.

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

**You do not need one to build, test, document or lint the crate.** `cargo test` runs the whole
offline suite with no network and no .NET on the machine. A bridge is required only for calls that
actually reach mzLib, and a missing one is a runtime error with instructions, never a build failure
— so contributors are never blocked by a 130 MB payload they may not want.

The simplest way to get one:

```rust
let bridge = mzlib::install::install_bridge()?;
```

It downloads the bridge pyMzLib published for your platform, verifies it against a checksum recorded
in this crate, unpacks it and caches it per-user. **It asks first, and nothing calls it for you** —
fetching the bridge is a user action, never a side effect of building. A built-in default download
would make a first `cargo check` pull 130 MB and break vendored and air-gapped builds, so `build.rs`
has no default URL and never will. mzLibR's `mzlibr_install_bridge()` is the same function for the
same reasons.

Once it is installed the crate finds it with nothing set. Full resolution order:

1. **`MZLIB_BRIDGE`** — a path to a bridge you already have. Always wins, checked at runtime too.
2. **`_dotnet/<runtime-identifier>/mzlib-bridge[.exe]`** beside the crate, staged by `build.rs` or by
   `scripts/stage-bridge.ps1`.
3. **The per-user cache**, where `install_bridge()` puts one.

The version it fetches is **pinned**, not "whichever is newest" — a crate that silently followed the
newest bridge could not be reproducible, and would change what it runs without anybody merging
anything. `.github/workflows/bridge-watch.yml` checks weekly and opens a pull request when pyMzLib
publishes a newer one, regenerating the pinned digests from that release's `SHA256SUMS`.

`build.rs` will also download at build time from **`MZLIB_BRIDGE_URL`** (verified against
`MZLIB_BRIDGE_SHA256`) if you set it. That takes a URL for a **bare executable** — it unpacks
nothing, so it is not how to consume the published `mzlib-bridge-<rid>.tar.gz`, whose payload is a
whole tree. Point it at an archive and the build script says so and stages nothing.

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

This section used to end by promising that "download-at-build becomes the default path once
pyMzLib's CI publishes the raw bridge binaries as release assets". Those assets exist as of
`v0.1.0.dev4`, and the promise was not kept — deliberately. Downloading by default is the one thing
this crate will not do, so the fetch became `install_bridge()` instead.

## Documentation

The API reference is rustdoc (`cargo doc --open`). Two things about it are deliberate:

- **Every example runs.** The examples are doctests, executed in CI against a stand-in bridge that
  answers each call from a fixture recorded from the real one — the same recordings pyMzLib's and
  mzLibR's examples replay — and only when the recording fits the call. The few that cannot run
  download from EBI, and say so.
- **The facts come from one spec per wire verb.** Parameters with their units, result fields with
  their units and what a null means, error kinds, caveats, the mzLib code wrapped, and the same
  verb's spelling in Python and R are rendered from language-neutral specs all three bindings share,
  and a test fails when a documented field drifts from its spec. See
  [docs/reference-facts.md](docs/reference-facts.md).

Changes are recorded in [CHANGELOG.md](CHANGELOG.md).

## Testing

```bash
cargo test          # the offline suite, the spec lint and every doc example: no network, no bridge
MZLIB_BRIDGE=… cargo test --features live    # the live canaries; they SKIP on an outage
```

The offline suite is the default because it must pass anywhere, in milliseconds. Live canaries are
the ones that would catch mzLib, PRIDE or UniProt changing under us, and they **skip rather than
fail** when a service is down — an ambiguous red build gets ignored, which is how a genuine contract
break goes unnoticed for a month.

See [docs/test-parity.md](docs/test-parity.md) for the test-by-test mapping against pyMzLib, and
[docs/findings.md](docs/findings.md) for defects this port surfaced upstream.

## Licence

Same as mzLib. See `LICENSE`.
