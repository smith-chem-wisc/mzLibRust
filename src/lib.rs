//! # mzLibRust — mzLib for Rust
//!
//! [mzLib](https://github.com/smith-chem-wisc/mzLib) is a mass-spectrometry and proteomics library
//! written in C#. mzLibRust makes its functionality callable from Rust with **no .NET installation**
//! — what it needs is a self-contained bridge executable, which
//! [`install::install_bridge`] fetches for your platform on request. The crate cannot carry it:
//! crates.io allows about 10 MB and the payload is roughly 130 MB. Nothing downloads it for you,
//! and the whole offline test suite passes without it.
//!
//! It is the sibling of [pyMzLib](https://github.com/smith-chem-wisc/pyMzLib) and speaks the same
//! **language-neutral bridge** — a self-contained executable exchanging a versioned JSON envelope
//! over stdin/stdout, which assumes nothing about the language calling it. Everything genuinely
//! hard already lives there: the mzLib interop, the composition of mzLib's own methods, and the
//! availability-versus-correctness error classification. This crate is the thin, idiomatic Rust
//! surface over it.
//!
//! ```
//! # mzlib_replay::activate();
//! let files = mzlib::pride::list_files("PXD000001")?;
//! println!(
//!     "{} files, {:.2} GB",
//!     files.len(),
//!     mzlib::pride::total_size_bytes(&files) as f64 / 1e9
//! );
//! # assert_eq!(files.len(), 8);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Reading is the widest surface, and it covers instrument data as well as search output.
//! [`readers::read_spectra`] reads **mzML**, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF and
//! msalign — scan headers always, peaks opt-in:
//!
//! ```
//! # mzlib_replay::activate();
//! # use mzlib::readers::{ReadOptions, SpectraOptions};
//! # let first_three = SpectraOptions { read: ReadOptions { limit: Some(3), ..Default::default() }, ..Default::default() };
//! let scans = mzlib::readers::read_spectra_with("sliced_ethcd.mzML", &first_three)?;
//! println!("{} scans", scans.scan_count);
//! # assert_eq!(scans.scan_count, 6);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! mzLib 1.0.592 recognises **36 file types** and this crate reads all of them;
//! [`readers::formats`] lists them from the mzLib the bridge carries.
//! [`readers::read_records`] reads any format into that format's own fields;
//! [`readers::read_results`], [`readers::read_features`], [`readers::read_matches`] and
//! [`readers::read_spectra`] project the four cross-format views; and
//! [`readers::read_protein_groups`], [`readers::read_quantified_peptides`] and
//! [`readers::read_occupancy`] read the MetaMorpheus and FlashLFQ quantification tables as long
//! tables. Every reader has a `_many` twin that reads a list of files in **one** bridge process
//! into one table — [`readers::read_spectra_many`] for every run of an experiment. See the
//! [`readers`] module.
//! SDRF experimental-design files have their own module, [`sdrf`], because `read_records`
//! cannot carry them without loss. Protein databases — what an accession is, which Ensembl gene
//! it resolves to, and whether a peptide is unique — are the [`proteins`] module.
//!
//! ## Two conventions worth knowing up front
//!
//! **Names follow mzLib.** A field here means exactly what it means in the mzLib source, the
//! MetaMorpheus output columns, and the papers — `match_between_runs`, `ppm_tolerance`,
//! `protein_groups`, `detection_type`. Nothing is renamed to look more Rust-like, because renaming
//! forces every reader to hold a translation table in their head.
//!
//! **The types disclose the traps.** Where mzLib applies an invisible rule, this crate surfaces it
//! rather than swallowing it. The clearest case is quantification: a *peptide* intensity is `f64`
//! and `0.0` means "not measured here", while a *protein* intensity is [`Option<f64>`] and `None`
//! means FlashLFQ could not resolve a number at all. In Python that distinction has to live in the
//! documentation and be remembered; here the compiler makes you handle it.
//!
//! ## Every example on these pages runs
//!
//! The examples are doctests, executed in CI against a stand-in bridge that answers each call from
//! a fixture recorded from the real one, and only when the recording fits the call. They are the
//! same recordings pyMzLib's and mzLibR's examples replay. The few that are not run say why: they
//! download from EBI.
//!
//! ## The reference facts come from the bridge's specs
//!
//! Each function that calls a wire verb carries that verb's facts — parameters with their units,
//! result fields with their units and what a null means, error kinds, caveats, the mzLib code it
//! wraps, and the same verb's spelling in Python and R — rendered from one language-neutral spec
//! per verb that all three bindings share. See `docs/reference-facts.md` in the repository.

#![forbid(unsafe_code)]

pub mod bridge;
pub mod flashlfq;
pub mod install;
pub mod peptidoform;
pub mod pride;
pub mod proteins;
pub mod readers;
pub mod sdrf;

pub use bridge::{
    bridge_path, bridge_version, BridgeVersion, MzLibError, OnError, Result, BRIDGE_ENV_VAR,
    PROTOCOL_VERSION, SERVICE_UNAVAILABLE_TYPE,
};
pub use install::{install_bridge, InstallOptions};

/// This crate's version, as declared in `Cargo.toml`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
