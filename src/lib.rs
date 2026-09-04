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
//! ```no_run
//! # fn main() -> Result<(), mzlib::MzLibError> {
//! let files = mzlib::pride::list_files("PXD000001")?;
//! println!(
//!     "{} files, {:.2} GB",
//!     files.len(),
//!     mzlib::pride::total_size_bytes(&files) as f64 / 1e9
//! );
//! # Ok(())
//! # }
//! ```
//!
//! Reading is the widest surface, and it covers instrument data as well as search output.
//! [`readers::read_spectra`] reads **mzML**, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF and
//! msalign — scan headers always, peaks opt-in:
//!
//! ```no_run
//! # fn main() -> Result<(), mzlib::MzLibError> {
//! let scans = mzlib::readers::read_spectra("run.mzML")?;
//! println!("{} scans", scans.scan_count);
//! # Ok(())
//! # }
//! ```
//!
//! mzLib recognises **31 file types** in all and this crate reads all of them.
//! [`readers::read_records`] reads any format into that format's own fields;
//! [`readers::read_results`], [`readers::read_features`], [`readers::read_matches`] and
//! [`readers::read_spectra`] project the four cross-format views. See the [`readers`] module.
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

#![forbid(unsafe_code)]

pub mod bridge;
pub mod flashlfq;
pub mod install;
pub mod peptidoform;
pub mod pride;
pub mod readers;

pub use bridge::{
    bridge_path, bridge_version, BridgeVersion, MzLibError, Result, BRIDGE_ENV_VAR,
    PROTOCOL_VERSION, SERVICE_UNAVAILABLE_TYPE,
};
pub use install::{install_bridge, InstallOptions};

/// This crate's version, as declared in `Cargo.toml`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
