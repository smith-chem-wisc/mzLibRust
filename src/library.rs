//! Spectral libraries: read an mzLib `.msl` library, every precursor, and its fragment ions on
//! request.
//!
//! A spectral library says, for each **precursor** (a peptide or proteoform at one charge), where a
//! search should expect to find it: its m/z, its retention time, its ion mobility, and the fragment
//! ions its MS2 spectrum should hold, with their relative intensities. MetaMorpheus and mzLib write
//! libraries in mzLib's binary `.msl` format, predicted (Prosit, through Koina) or built from a
//! search. This module reads them with mzLib's own reader, so a library reads in Rust exactly as
//! MetaMorpheus reads it.
//!
//! | You want to | Call | mzLib type |
//! |---|---|---|
//! | See every precursor in a library | [`read`], [`read_with`] | `MslFileTypeHandler.Open`, `MslLibrary.GetEntry` |
//! | Know whether a retention time is iRT or minutes | the `retention_time_unit` column | `MslLibraryEntry.RtIsCalibrated` |
//! | Count targets, decoys and entrapment precursors | the `decoy_contaminant_target` column | `DecoyContaminantTargetLabel.For` |
//! | See a precursor's fragment ions | [`LibraryReadOptions::fragments`] | `MslFragmentIon` |
//! | Read a library too large for memory | [`LibraryReadOptions::index_only_threshold`] | `MslLibrary.LoadIndexOnly` |
//! | Read several libraries into one table | [`read_many`] | the same reader, in one bridge process |
//!
//! Every example on this page runs in CI against output recorded from the real bridge. Two small
//! libraries are used, the ones pyMzLib's examples read, in `tests/fixtures/library/`.
//! `myPrositLib.msl` is mzLib's own test library `myPrositLib.msp` (five Prosit-predicted
//! precursors) converted to `.msl` by mzLib. `synthetic_library.msl` holds six precursors written to
//! cover every field: a target, its decoy, an entrapment target whose retention time is in minutes,
//! an entrapment decoy, a predicted precursor with nothing measured, and a proteoform.
//!
//! ## Every precursor
//!
//! ```
//! # mzlib_replay::activate();
//! let library = mzlib::library::read("synthetic_library.msl")?;
//! assert_eq!(
//!     (library.record_count, library.target_count, library.decoy_count, library.proteoform_count),
//!     (6, 4, 2, 1)
//! );
//! let c = &library.columns;
//! let (index, sequence, charge) = (c.integers("precursor_index")?, c.strings("full_sequence")?, c.integers("charge_state")?);
//! let (mz, label) = (c.floats("precursor_mz")?, c.strings("decoy_contaminant_target")?);
//! for i in 0..library.columns.rows() {
//!     println!("{:?} {:?} {:?} {:.4} {:?}", index[i], sequence[i], charge[i], mz[i].unwrap(), label[i]);
//! }
//! //   0 PEPTIDEK 2 465.7400 T
//! //   1 KEDITPEP 2 465.7401 D
//! //   2 PEPTLDEK 2 458.7462 ET
//! //   3 KEDLTPEP 2 458.7463 ED
//! //   4 ELVISLIVESK 3 406.2420 T
//! //   5 MAHVRSGKHQIPAEALQ 2 952.4930 T
//! assert_eq!(sequence[5].as_deref(), Some("MAHVRSGKHQIPAEALQ"));
//! assert_eq!(label.iter().map(|l| l.as_deref().unwrap()).collect::<Vec<_>>(), ["T", "D", "ET", "ED", "T", "T"]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! One row per precursor, **in the library's stored order**. `precursor_index` is the position in
//! the file and the key mzLib's `MslLibrary.GetEntry` takes. It is not m/z order: mzLib's own
//! `GetAllEntries` enumerates by ascending m/z, so sort on `precursor_mz` if that is the order you
//! want. [`LibraryRecords::columns`] holds one value per row for each column; the typed accessors
//! of [`Table`] project them.
//!
//! ### What each column is
//!
//! The names are mzLib's `MslLibraryEntry` properties. Where mzLib stores a value it documents as
//! "not available", the column holds `None`, and the table says which:
//!
//! | column | what it is | `None` means |
//! |---|---|---|
//! | `precursor_index` | 0-based position in the library | never `None` |
//! | `full_sequence`, `base_sequence` | mzLib bracket notation, e.g. `PEPTM[Common Variable:Oxidation on M]IDE`; the bare residues | never `None` |
//! | `precursor_mz` | precursor m/z, (M + zH)/z | never `None` |
//! | `charge_state` | precursor charge | never `None` |
//! | `retention_time`, `retention_time_unit` | the time, and whether it is `"iRT"` or `"minutes"` | never `None` |
//! | `ion_mobility` | 1/K0, in V·s/cm² | mzLib stored 0.0, "not available" |
//! | `is_decoy`, `is_entrapment`, `decoy_contaminant_target` | target/decoy/entrapment flags, and mzLib's label `T`, `D`, `ET` or `ED` | never `None` |
//! | `protein_accession`, `protein_name`, `gene_name` | the source protein | mzLib stored an empty string: no protein |
//! | `is_proteotypic` | the library writer's proteotypic flag | never `None` |
//! | `q_value` | the library's q-value for the precursor | mzLib stored NaN: no q-value, as in a predicted library |
//! | `elution_group_id` | shared by every precursor of one base sequence | never `None` |
//! | `molecule_type` | `Peptide`, `Proteoform`, `Oligonucleotide` or `Glycopeptide` | never `None` |
//! | `source_type` | `Predicted`, `Empirical`, or `EmpiricalRefined` (predicted intensities, measured time or mobility) | never `None` |
//! | `dissociation_type` | e.g. `HCD`; `Unknown` when the library does not say | never `None` |
//! | `nce` | normalized collision energy, a whole number | mzLib stored 0: unknown |
//!
//! The predicted precursor (row 4) shows the "not available" values, and the empirical target
//! (row 0) shows them filled in:
//!
//! ```
//! # mzlib_replay::activate();
//! let library = mzlib::library::read("synthetic_library.msl")?;
//! let c = &library.columns;
//! assert_eq!(c.floats("ion_mobility")?[0], Some(0.949999988079071));
//! assert_eq!(c.floats("ion_mobility")?[4], None);
//! assert_eq!((c.floats("q_value")?[0], c.floats("q_value")?[4]), (Some(0.001), None));
//! assert_eq!((c.integers("nce")?[0], c.integers("nce")?[4]), (Some(28), None));
//! assert_eq!(c.strings("protein_accession")?[0].as_deref(), Some("P02768"));
//! assert_eq!(c.strings("protein_accession")?[4], None);
//! let source = c.strings("source_type")?;
//! assert_eq!((source[0].as_deref(), source[4].as_deref()), (Some("Empirical"), Some("Predicted")));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! A whole column can have no source at all. `myPrositLib.msl` has no protein table and no ion
//! mobility, so those columns are named in [`LibraryRecords::absent_fields`] and are `None` in
//! every row:
//!
//! ```
//! # mzlib_replay::activate();
//! let prosit = mzlib::library::read("myPrositLib.msl")?;
//! assert_eq!(prosit.absent_fields, ["ion_mobility", "protein_accession", "protein_name", "gene_name"]);
//! assert_eq!(prosit.columns.integers("nce")?, [None; 5]);
//! assert_eq!(prosit.columns.strings("dissociation_type")?[0].as_deref(), Some("Unknown"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! The `.msp` it was converted from records `Collision_energy=30.0`, yet `nce` is `None`. That is
//! not a reading error. mzLib converts a text library through `LibrarySpectrum`, which has no field
//! for NCE, dissociation type, q-value, ion mobility or protein, so
//! `MslLibraryEntry.FromLibrarySpectrum` writes mzLib's "not available" value for each. A converted
//! library has none of them.
//!
//! ## Retention time: iRT or minutes, row by row
//!
//! A library stores each precursor's retention time either as an **indexed retention time** (iRT,
//! a unitless scale anchored on reference peptides) or as **minutes** in one particular run. Which
//! one is a flag on each precursor (`MslLibraryEntry.RtIsCalibrated`), and one library can hold
//! both. So the unit is a column, not a property of the file:
//!
//! ```
//! # mzlib_replay::activate();
//! let library = mzlib::library::read("synthetic_library.msl")?;
//! let time = library.columns.floats("retention_time")?;
//! let unit = library.columns.strings("retention_time_unit")?;
//! assert_eq!(
//!     unit[..4].iter().map(|u| u.as_deref().unwrap()).collect::<Vec<_>>(),
//!     ["iRT", "iRT", "minutes", "iRT"]
//! );
//! assert_eq!(format!("{:.2}", time[2].unwrap()), "21.40");   // PEPTLDEK, in run minutes
//! // Keep the rows of one unit before you use the numbers.
//! let irt: Vec<f64> = time
//!     .iter()
//!     .zip(&unit)
//!     .filter(|(_, u)| u.as_deref() == Some("iRT"))
//!     .filter_map(|(t, _)| *t)
//!     .collect();
//! assert_eq!((irt.len(), library.columns.rows()), (5, 6));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Never plot or compare `retention_time` without `retention_time_unit`.** This crate converts
//! nothing: iRT becomes minutes only through a calibration against a run, which is not in the
//! library. iRT can be negative; Prosit predicts `-15.77` for one of the `myPrositLib.msl`
//! precursors. A `retention_time` of exactly 0 deserves suspicion: mzLib's conversion from a text
//! library writes 0 when the spectrum had no retention time, and documents no sentinel for it, so
//! 0 is passed through as a value.
//!
//! ## Targets, decoys and entrapment
//!
//! **The library is unfiltered.** Decoys and entrapment precursors are rows like any other.
//! `decoy_contaminant_target` is mzLib's own label, written by `DecoyContaminantTargetLabel.For`:
//! `T` target, `D` decoy, `ET` entrapment target, `ED` entrapment decoy. The `.msl` format has no
//! contaminant flag, so `C` never appears. Entrapment is a target for the search: it competes with
//! the real targets, so [`LibraryRecords::target_count`] includes the entrapment targets and
//! [`LibraryRecords::decoy_count`] the entrapment decoys. To count only the real targets, keep
//! `decoy_contaminant_target == "T"`:
//!
//! ```
//! # mzlib_replay::activate();
//! let library = mzlib::library::read("synthetic_library.msl")?;
//! let decoy = library.columns.booleans("is_decoy")?;
//! let sequence = library.columns.strings("full_sequence")?;
//! let targets: Vec<&str> = (0..decoy.len())
//!     .filter(|&i| decoy[i] == Some(false))
//!     .filter_map(|i| sequence[i].as_deref())
//!     .collect();
//! assert_eq!(targets, ["PEPTIDEK", "PEPTLDEK", "ELVISLIVESK", "MAHVRSGKHQIPAEALQ"]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! ## Fragment ions
//!
//! With [`LibraryReadOptions::fragments`] the table is **one row per fragment ion**, the
//! precursor's columns repeated on each row, then the fragment's own. `limit` and `offset` still
//! count precursors:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::library::{read_with, LibraryReadOptions};
//!
//! let ions = read_with(
//!     "synthetic_library.msl",
//!     &LibraryReadOptions { fragments: true, limit: Some(1), ..Default::default() },
//! )?;
//! assert_eq!((ions.returned_count, ions.row_count), (1, 6));
//! let c = &ions.columns;
//! let (annotation, mz) = (c.strings("annotation")?, c.floats("fragment_mz")?);
//! let (intensity, position, loss) =
//!     (c.floats("fragment_intensity")?, c.integers("residue_position")?, c.floats("neutral_loss")?);
//! for i in 0..annotation.len() {
//!     println!("{:11} {:9.4} {:5.2} {:?} {:?}", annotation[i].as_deref().unwrap(), mz[i].unwrap(),
//!              intensity[i].unwrap(), position[i], loss[i]);
//! }
//! //   y1           147.1128  0.40 Some(8) Some(0.0)
//! //   b2           227.1026  0.30 Some(2) Some(0.0)
//! //   y2           276.1554  1.00 Some(7) Some(0.0)
//! //   bIy[2-5]     326.1710  0.10 Some(2) Some(0.0)
//! //   y3-18.0106   373.2082  0.20 Some(6) Some(-18.010565)
//! //   y4-43.0058   461.2242  0.15 None Some(-43.005814)
//! assert_eq!(
//!     annotation.iter().map(|a| a.as_deref().unwrap()).collect::<Vec<_>>(),
//!     ["y1", "b2", "y2", "bIy[2-5]", "y3-18.0106", "y4-43.0058"]
//! );
//! assert_eq!((position[5], loss[5]), (None, Some(-43.005814)));   // a custom loss: no position
//! assert_eq!(loss[0], Some(0.0));                                  // no loss is 0, not None
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! What the numbers are:
//!
//! - **`fragment_intensity` is relative within its precursor, 0 to 1.** mzLib's writer scales each
//!   precursor's most intense fragment to 1, so intensities are comparable within a precursor and
//!   not across precursors.
//! - **`neutral_loss` is in Da, negative for a loss, and 0 means no loss.** It is 0, not `None`,
//!   because 0 Da is the value of "no loss", so `fragment_mz` and `neutral_loss` stay usable
//!   together.
//! - **`annotation` is mzLib's library rendering** (`MslFragmentIon.Annotation`): the loss signed
//!   to four decimals and, for a charge above 1, the charge appended (`y4^2+`). It differs from
//!   `Product.Annotation` elsewhere in mzLib, so do not match the two strings against each other.
//! - **An internal ion** (`bIy[2-5]`) has `is_internal` true, its C-terminal type in
//!   `secondary_product_type` and its end residue in `secondary_fragment_number`. Both are `None`
//!   on terminal ions, where mzLib stores 0 and null.
//! - **`residue_position` is `None` for a fragment with a custom neutral loss**, like the
//!   `-43.0058` above. The `.msl` format reuses that field to index its table of custom losses, and
//!   mzLib reads it back as 0. The loss itself survives. Named losses (water, ammonia, phosphate)
//!   keep their position.
//! - **`fragment_count`** on each row is how many fragment ions the precursor has. A precursor
//!   with none is still one row, with `fragment_count` 0 and every fragment cell `None`.
//!
//! Fragments are stored, and returned, in m/z order: mzLib's writer sorts them, so the table is a
//! spectrum to plot as is. A library has tens of fragment ions per precursor. For a whole library,
//! write the table to disk with [`LibraryReadOptions::out`] (a `.tsv`) and read that, rather than
//! returning it as one document.
//!
//! ## Large libraries: index-only
//!
//! mzLib opens a library one of two ways. **Loaded whole**, every precursor and every fragment is
//! in memory. **Index-only**, the precursors are in memory and each precursor's fragment block is
//! read from disk when it is asked for. mzLib chooses by file size (`MslFileTypeHandler.Open`): a
//! library at or above `MslFileTypeHandler.DefaultIndexOnlyThresholdBytes`, 1 GiB in the mzLib this
//! crate's bridge is built from, opens index-only. Nothing is sent unless you set
//! [`LibraryReadOptions::index_only_threshold`], a size in bytes; `Some(0)` opens every library
//! index-only:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::library::{read_with, LibraryReadOptions};
//!
//! let small = read_with(
//!     "synthetic_library.msl",
//!     &LibraryReadOptions { index_only_threshold: Some(0), limit: Some(2), ..Default::default() },
//! )?;
//! assert_eq!((small.is_index_only, small.index_only_threshold), (true, 0));
//! assert_eq!(
//!     small.columns.strings("full_sequence")?,
//!     [Some("PEPTIDEK".to_owned()), Some("KEDITPEP".to_owned())]
//! );
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **The rows are the same in both modes**; only memory and speed differ. Two things to know:
//!
//! - A **compressed** library ([`LibraryRecords::is_compressed`]) always loads whole: its fragment
//!   offsets point into the decompressed data, not into the file, so mzLib decompresses it all and
//!   `is_index_only` stays `false` whatever the threshold.
//! - Index-only **with fragments** reads every fragment block from disk twice: once to count the
//!   rows, once to write them. mzLib keeps each precursor's stored fragment count to itself, so
//!   there is no cheaper way to size the table. Without fragments, index-only reads no fragment
//!   block at all, which is the fast way to list the precursors of a very large library.
//!
//! ## Many libraries at once
//!
//! [`read_many`] reads a list of libraries into one long table, in one bridge process. Its first
//! two columns say where each row came from (`source_index`, `source_path`), and
//! [`LibraryBatch::files`] holds one report per library:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::library::{read_many, LibraryBulkOptions};
//!
//! let batch = read_many(&["synthetic_library.msl", "myPrositLib.msl"], &LibraryBulkOptions::default())?;
//! assert_eq!(batch.record_count, 11);
//! assert_eq!(batch.files.iter().map(|f| f.record_count).collect::<Vec<_>>(), [Some(6), Some(5)]);
//! assert_eq!(
//!     batch.columns.integers("source_index")?,
//!     [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1].map(Some)
//! );
//! assert!(batch.files[0].absent_fields.is_empty());
//! assert_eq!(batch.files[1].absent_fields, ["ion_mobility", "protein_accession", "protein_name", "gene_name"]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! [`BulkOptions::threads`](crate::readers::BulkOptions::threads) reads that many libraries at
//! once. The answer is identical at any value; the default is 1 because a library loaded whole is
//! held in memory while its rows are written.
//! [`OnError::Skip`](crate::OnError::Skip) records an unreadable library in its report's `error`
//! and reads the rest.
//!
//! ## Traps
//!
//! - **Stored order, not m/z order.** Rows follow `precursor_index`. Sort on `precursor_mz`
//!   yourself.
//! - **The extra digits are the float.** `precursor_mz`, `retention_time` and `ion_mobility` are
//!   stored as 32-bit floats and widened by mzLib, so 465.74 reads back as `465.739990234375`.
//!   Round them for display; they carry no more precision than about seven significant digits.
//! - **`nce` is a whole number.** The format stores NCE × 10 and mzLib reads it back with integer
//!   division, so an NCE of 27.5 reads as 27.
//! - **Converted libraries are thin.** A library converted from `.msp` has no NCE, dissociation
//!   type, q-value, ion mobility or protein (see above).
//! - **A broken file is an error, not a short table.** mzLib checks the whole file (magic number,
//!   version, footer, precursor count and a CRC-32 checksum) before it returns a single precursor.
//!   A damaged library is [`MzLibError::Bridge`] with `error_type` `"InvalidDataException"`, and a
//!   file that is not a library is one with `"FormatException"`. A file whose name does not end in
//!   `.msl` is [`MzLibError::Usage`] before anything is opened.
//!
//! ## What this does not do
//!
//! - **Read a text library** (`.msp`, pDeep, ms2pip). mzLib can turn one into a library in memory
//!   (`MslConverter.FromMspFile`), but no wire verb projects it yet.
//! - **Write a library.** mzLib writes them (`MslLibrary.Save`); no wire verb projects that yet.
//! - **Query by m/z or retention-time window.** mzLib's `QueryMzWindow` and `QueryWindow` are what
//!   a search engine calls. To find a precursor here, read the table and filter it.
//!
//! Ported from pyMzLib's `pymzlib.library`, which decided the verb, the wire fields and the
//! caveats. mzLib #1383 added the entrapment flag and made `rt_is_calibrated` work, #1421 added
//! `GetEntry` without fragments, and #1420 made libraries over 2 GB open.
//!
//! Wire verb: `library read`.
//!
//! ## Cite
//!
#![doc = include_str!("../docs/reference/cite.library.md")]

use std::path::Path;
use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, OnError, Result};
use crate::readers::{BulkOptions, ExcludedField, ReadError, Table, WrittenTable};

/// The first pyMzLib whose bridge has `library read`: its spec's `since.pymzlib`.
const LIBRARY_SINCE: &str = "0.6.0";

// ---------------------------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------------------------

/// How [`read_with`] reads one library. Every default is mzLib's, and nothing is sent for it.
#[derive(Debug, Clone, Default)]
pub struct LibraryReadOptions {
    /// One row per **fragment ion** instead of one per precursor, with the precursor's columns
    /// repeated on each row. Off by default: a library has tens of fragments per precursor, so
    /// pair it with [`Self::limit`] or [`Self::out`] for a large one.
    pub fragments: bool,
    /// A library this many **bytes** or larger opens index-only (precursors in memory, fragment
    /// blocks read from disk as needed); a smaller one loads whole. `None` (the default) sends
    /// nothing, so mzLib's own default applies: `MslFileTypeHandler.DefaultIndexOnlyThresholdBytes`,
    /// 1073741824 bytes (1 GiB) at the pin. `Some(0)` opens every library index-only. The rows are
    /// the same either way.
    pub index_only_threshold: Option<u64>,
    /// Return at most this many **precursors**; with [`Self::fragments`] each gives one row per
    /// fragment ion. `None` returns all of them.
    pub limit: Option<u64>,
    /// Skip this many **precursors**, in `precursor_index` order. A window, not a cursor: the
    /// library is opened again on every call.
    pub offset: u64,
    /// Write the table to this path as **tab-separated** text and return only a summary
    /// ([`LibraryRecords::output`]). Must differ from the input; parent directories are created.
    /// The path must end in `.tsv` (any case); any other is [`MzLibError::Usage`] from the bridge
    /// before the library is opened.
    pub out: Option<String>,
    /// Time to allow. `None` (the default) waits indefinitely.
    pub timeout: Option<Duration>,
}

/// How [`read_many`] reads its libraries. Every default is mzLib's.
#[derive(Debug, Clone, Default)]
pub struct LibraryBulkOptions {
    /// Libraries read at once, what one unreadable library does, and where the long table goes.
    /// The answer is identical at any thread count; the default is 1 because a library loaded
    /// whole is held in memory while its rows are written.
    pub bulk: BulkOptions,
    /// One row per fragment ion, in every library. Off by default: with it, write the batch to
    /// [`BulkOptions::out`].
    pub fragments: bool,
    /// The size in **bytes** at and above which a library opens index-only, for every library.
    /// `None` (the default) leaves it to mzLib (1 GiB at the pin).
    pub index_only_threshold: Option<u64>,
}

// ---------------------------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------------------------

/// An mzLib spectral library read into one table: what [`read_with`] returns.
///
/// Without fragments, [`Self::columns`] is one row per precursor; with
/// [`LibraryReadOptions::fragments`], one row per fragment ion. The columns and what a `None` in
/// each means are tabled in the [module guide](crate::library#what-each-column-is).
#[derive(Debug, Clone, Deserialize)]
pub struct LibraryRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Always `None`: `.msl` is not one of mzLib's `SupportedFileType` members (mzLib reads it with
    /// `MslLibrary`, outside `FileReader`), so there is no type name to give.
    #[serde(default)]
    pub file_type: Option<String>,
    /// `"MslLibrary"`, the mzLib class that read the file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub reader: String,
    /// The file's `.msl` format version (`MslFileHeader.FormatVersion`), 1 to 3 at the pin.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub format_version: u32,
    /// Whether the fragment section is zstd-compressed. A compressed library always loads whole.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub is_compressed: bool,
    /// Whether mzLib opened it index-only: at or above [`Self::index_only_threshold`], and not
    /// compressed.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub is_index_only: bool,
    /// Precursors that are not decoys, entrapment targets included.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub target_count: u64,
    /// Decoy precursors, entrapment decoys included.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub decoy_count: u64,
    /// Precursors whose `molecule_type` is `"Proteoform"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub proteoform_count: u64,
    /// The smallest precursor m/z. `None`: the library has no precursors.
    #[serde(default)]
    pub min_precursor_mz: Option<f64>,
    /// The largest precursor m/z. `None`: the library has no precursors.
    #[serde(default)]
    pub max_precursor_mz: Option<f64>,
    /// Precursors in the whole library, before the offset and limit.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Always `None`: data rows are not counted for a binary library. mzLib checks the whole file
    /// (magic, version, footer, precursor count, CRC-32) before returning any precursor, so every
    /// one is read or the call fails.
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// Always `None`: the unit varies by precursor, so it is the `retention_time_unit` column.
    #[serde(default)]
    pub retention_time_unit: Option<String>,
    /// What this table cannot be trusted to mean, each citing the mzLib source.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Columns this library has no source for, from its file flags, so every value is `None`:
    /// `ion_mobility` when no precursor has one, `protein_accession` and `protein_name` when there
    /// is no protein table, `gene_name` when there are no gene names.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"field: ExceptionType"` for each column whose read threw on some returned rows; those
    /// cells are `None`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Always empty: every field crosses, the fragments with [`LibraryReadOptions::fragments`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// Whether the rows are fragment ions.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub fragments_included: bool,
    /// The size, in bytes, at and above which mzLib opened the library index-only: the one you
    /// gave, or mzLib's default (1073741824 bytes, 1 GiB, at the pin).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub index_only_threshold: u64,
    /// Precursors returned in [`Self::columns`]: precursors, not rows, the unit the limit and
    /// offset count in. Zero when the table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// Rows in [`Self::columns`]: [`Self::returned_count`] without fragments, one per fragment ion
    /// (or one for a precursor with none) with them. Zero when the table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Precursors skipped: the offset applied, in precursors.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// Whether precursors were left behind, by either the limit or the offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// The table; its column order is [`Table::names`]. Empty when the table went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
}

/// The per-library facts of one input to [`read_many`]: one entry of [`LibraryBatch::files`].
///
/// The facts a single [`read_with`] reports at its top level, under the same names. For an input
/// that failed under [`OnError::Skip`], [`Self::error`] says why, and every fact that could not be
/// established is `None` or empty: **not zero**, because nothing was counted.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct LibraryFileReport {
    /// The absolute path of the input.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Always `None`: `.msl` is not an mzLib `SupportedFileType`.
    #[serde(default)]
    pub file_type: Option<String>,
    /// `"MslLibrary"`; `None` for a failed input.
    #[serde(default)]
    pub reader: Option<String>,
    /// The `.msl` format version. `None`: a failed input.
    #[serde(default)]
    pub format_version: Option<u32>,
    /// Whether the fragment section is compressed. `None`: a failed input.
    #[serde(default)]
    pub is_compressed: Option<bool>,
    /// Whether mzLib opened it index-only. `None`: a failed input.
    #[serde(default)]
    pub is_index_only: Option<bool>,
    /// Precursors that are not decoys, entrapment targets included. `None`: a failed input.
    #[serde(default)]
    pub target_count: Option<u64>,
    /// Decoy precursors, entrapment decoys included. `None`: a failed input.
    #[serde(default)]
    pub decoy_count: Option<u64>,
    /// Proteoform precursors. `None`: a failed input.
    #[serde(default)]
    pub proteoform_count: Option<u64>,
    /// The smallest precursor m/z. `None`: no precursors, or a failed input.
    #[serde(default)]
    pub min_precursor_mz: Option<f64>,
    /// The largest precursor m/z. `None`: no precursors, or a failed input.
    #[serde(default)]
    pub max_precursor_mz: Option<f64>,
    /// Precursors in the library. `None`: a failed input.
    #[serde(default)]
    pub record_count: Option<u64>,
    /// Always `None` (see [`LibraryRecords::rows_not_read`]).
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// Always `None`: the unit is the `retention_time_unit` column.
    #[serde(default)]
    pub retention_time_unit: Option<String>,
    /// What this library's rows cannot be trusted to mean.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// This library's columns, without the `source_index`/`source_path` pair the batch adds.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub column_names: Vec<String>,
    /// Columns this library has no source for.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"field: ExceptionType"` for each column whose read threw on some rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Always empty.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// Why this input could not be read, under [`OnError::Skip`]; `None` when it was.
    #[serde(default)]
    pub error: Option<ReadError>,
}

impl LibraryFileReport {
    /// Whether this input was read ([`Self::error`] is `None`).
    #[must_use]
    pub fn ok(&self) -> bool {
        self.error.is_none()
    }
}

/// Many libraries read into ONE long table: what [`read_many`] returns.
///
/// The first two columns say where each row came from: `source_index`, the library's 0-based
/// position in the list you passed, and `source_path`, its absolute path. The rest are
/// [`LibraryRecords`]' columns. Rows are grouped by library in input order, each library's rows in
/// `precursor_index` order, **whatever the thread count was**.
#[derive(Debug, Clone, Deserialize)]
pub struct LibraryBatch {
    /// Libraries in the list you passed, in files.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Libraries read, in files.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Libraries that could not be read, in files: always 0 unless [`OnError::Skip`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// [`OnError::Fail`] or [`OnError::Skip`], as requested. The thread count is deliberately not
    /// echoed: the result is identical at every thread count.
    #[serde(default)]
    pub on_error: OnError,
    /// Precursors, summed over the libraries read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Precursors carried back in [`Self::columns`], summed over the libraries. Zero when the
    /// table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// Rows in [`Self::columns`], one per fragment ion with fragments. Zero when the table was
    /// written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Whether the rows are fragment ions.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub fragments_included: bool,
    /// The size in bytes at and above which a library opened index-only: the one given, or
    /// mzLib's default.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub index_only_threshold: u64,
    /// One long table: `source_index`, `source_path`, then the library columns. Empty when the
    /// table was written to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, one library at a time, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
    /// One [`LibraryFileReport`] per input, in input order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<LibraryFileReport>,
}

impl LibraryBatch {
    /// The inputs that could not be read, under [`OnError::Skip`]. Empty otherwise.
    #[must_use]
    pub fn failed_files(&self) -> Vec<&LibraryFileReport> {
        self.files.iter().filter(|file| !file.ok()).collect()
    }
}

// ---------------------------------------------------------------------------------------------
// Argument assembly
// ---------------------------------------------------------------------------------------------

/// `--fragments` and `--index-only-threshold`, only when given.
fn library_options(fragments: bool, index_only_threshold: Option<u64>) -> Vec<String> {
    let mut args = Vec::new();
    if fragments {
        args.push("--fragments".to_owned());
    }
    if let Some(threshold) = index_only_threshold {
        args.push("--index-only-threshold".to_owned());
        args.push(threshold.to_string());
    }
    args
}

fn out_args(out: Option<&str>) -> Result<Vec<String>> {
    match out {
        None => Ok(Vec::new()),
        Some(out) if out.trim().is_empty() => Err(MzLibError::Usage(
            "out must be a non-empty path ending in .tsv, or None to return the rows.".to_owned(),
        )),
        Some(out) => Ok(vec!["--out".to_owned(), out.trim().to_owned()]),
    }
}

/// The arguments of a one-library `library read` call, validated before anything is spawned.
fn read_request(path: &Path, options: &LibraryReadOptions) -> Result<Vec<String>> {
    let text = path.to_str().ok_or_else(|| {
        MzLibError::Usage(
            "The library path is not valid UTF-8, which the bridge requires.".to_owned(),
        )
    })?;
    if text.trim().is_empty() {
        return Err(MzLibError::Usage(
            "A library path is required, e.g. 'myLibrary.msl'.".to_owned(),
        ));
    }
    let mut args = vec![
        "library".to_owned(),
        "read".to_owned(),
        "--path".to_owned(),
        text.trim().to_owned(),
    ];
    if let Some(limit) = options.limit {
        args.push("--limit".to_owned());
        args.push(limit.to_string());
    }
    if options.offset > 0 {
        args.push("--offset".to_owned());
        args.push(options.offset.to_string());
    }
    args.extend(out_args(options.out.as_deref())?);
    args.extend(library_options(
        options.fragments,
        options.index_only_threshold,
    ));
    Ok(args)
}

/// The arguments and stdin of a `library read --paths-stdin` call. The whole list goes to ONE
/// bridge process; nothing here loops over files (the bridge's `design/PARALLELISM.md`).
fn read_many_request<P: AsRef<Path>>(
    paths: &[P],
    options: &LibraryBulkOptions,
) -> Result<(Vec<String>, String)> {
    let lines = bridge::path_lines(paths, "library")?;
    let mut args = vec![
        "library".to_owned(),
        "read".to_owned(),
        "--paths-stdin".to_owned(),
        "--threads".to_owned(),
        bridge::threads_arg(options.bulk.threads)?,
        "--on-error".to_owned(),
        options.bulk.on_error.as_str().to_owned(),
    ];
    args.extend(out_args(options.bulk.out.as_deref())?);
    args.extend(library_options(
        options.fragments,
        options.index_only_threshold,
    ));
    Ok((args, lines.join("\n") + "\n"))
}

fn call<T: serde::de::DeserializeOwned>(
    args: &[String],
    stdin: Option<&str>,
    timeout: Option<Duration>,
) -> Result<T> {
    let data = bridge::invoke(args, stdin, timeout)?;
    serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!("library payload could not be interpreted: {error}"))
    })
}

// ---------------------------------------------------------------------------------------------
// The public surface
// ---------------------------------------------------------------------------------------------

/// Read an mzLib spectral library (`.msl`), one row per precursor, with every default.
///
/// See [`read_with`] for the reference.
///
/// # Errors
///
/// As [`read_with`].
pub fn read(path: impl AsRef<Path>) -> Result<LibraryRecords> {
    read_with(path, &LibraryReadOptions::default())
}

/// Read an mzLib spectral library (`.msl`) as one row per precursor, or with fragments one row per
/// fragment ion.
///
/// Wraps mzLib's `MslFileTypeHandler.Open` and `MslLibrary.GetEntry`. Rows come in the library's
/// stored order (`precursor_index`), not by m/z; sort on `precursor_mz` for that. mzLib's "not
/// available" values (an ion mobility of 0.0, a NaN q-value, an NCE of 0, an empty protein
/// accession) arrive as `None`. Any extension other than `.msl` is [`MzLibError::Usage`] before the
/// file is opened: text libraries (`.msp`) are not read here. The
/// [module guide](crate::library) walks through every column.
#[doc = include_str!("../docs/reference/library.read.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for a blank path or a blank `out`, and when
/// the bridge in use is older than pyMzLib 0.6.0's, which has no `library read`.
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::library::{read_with, LibraryReadOptions};
///
/// let library = read_with("synthetic_library.msl", &LibraryReadOptions::default())?;
/// assert_eq!(
///     library.columns.strings("decoy_contaminant_target")?,
///     ["T", "D", "ET", "ED", "T", "T"].map(|s| Some(s.to_owned()))
/// );
/// let unit = library.columns.strings("retention_time_unit")?;
/// assert_eq!(unit[2].as_deref(), Some("minutes"));
/// // The predicted precursor: no q-value and no ion mobility are "not available", so None.
/// assert_eq!(library.columns.floats("q_value")?[4], None);
/// assert_eq!(library.columns.floats("ion_mobility")?[4], None);
///
/// // One precursor's fragment ions.
/// let ions = read_with(
///     "synthetic_library.msl",
///     &LibraryReadOptions { fragments: true, limit: Some(1), ..Default::default() },
/// )?;
/// assert_eq!((ions.returned_count, ions.row_count), (1, 6));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/library.read.see-also.md")]
pub fn read_with(path: impl AsRef<Path>, options: &LibraryReadOptions) -> Result<LibraryRecords> {
    // Every argument is checked before any process is started, the version probe included.
    let args = read_request(path.as_ref(), options)?;
    bridge::require_verb("library read", LIBRARY_SINCE)?;
    call(&args, None, options.timeout)
}

/// Read many spectral libraries into one long table, in one bridge process.
///
/// [`read_with`] for a list: a predicted library and the empirical one refined from it, say. Each
/// library's facts (counts, load mode, format version) are in its [`LibraryFileReport`]. A
/// repeated path is refused rather than read twice. With [`OnError::Fail`] (the default) the first
/// library that cannot be read stops the batch, its message starting `Input <i> (<path>)`; with
/// [`OnError::Skip`] its report's `error` says why and the rest are read.
#[doc = include_str!("../docs/reference/library.read.bulk.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for an empty list, a blank path or one holding
/// a tab or line break, a thread count of 0 or below -1, or a blank `out`; and when the bridge in
/// use is older than pyMzLib 0.6.0's.
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::library::{read_many, LibraryBulkOptions};
///
/// let batch = read_many(&["synthetic_library.msl", "myPrositLib.msl"], &LibraryBulkOptions::default())?;
/// assert_eq!((batch.file_count, batch.read_count, batch.record_count), (2, 2, 11));
/// assert_eq!(batch.files[1].record_count, Some(5));
/// assert!(batch.failed_files().is_empty());
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_many<P: AsRef<Path>>(
    paths: &[P],
    options: &LibraryBulkOptions,
) -> Result<LibraryBatch> {
    let (args, stdin) = read_many_request(paths, options)?;
    bridge::require_verb("library read", LIBRARY_SINCE)?;
    call(&args, Some(&stdin), options.bulk.timeout)
}

#[cfg(test)]
mod tests {
    //! Offline. The recordings are pyMzLib's, byte for byte, made from the two libraries in
    //! `tests/fixtures/library/`.

    use super::*;

    fn recorded<T: serde::de::DeserializeOwned>(name: &str) -> T {
        let text = std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("fixtures")
                .join(name),
        )
        .unwrap();
        serde_json::from_str(&text).unwrap()
    }

    #[test]
    fn a_plain_read_sends_only_the_path() {
        let args = read_request(
            Path::new("synthetic_library.msl"),
            &LibraryReadOptions::default(),
        )
        .unwrap();
        assert_eq!(args, ["library", "read", "--path", "synthetic_library.msl"]);
        assert_eq!(LibraryReadOptions::default().timeout, None);
    }

    #[test]
    fn every_option_is_sent_under_its_wire_name() {
        let options = LibraryReadOptions {
            fragments: true,
            index_only_threshold: Some(0),
            limit: Some(2),
            offset: 1,
            out: Some("ions.tsv".to_owned()),
            timeout: None,
        };
        let args = read_request(Path::new("synthetic_library.msl"), &options).unwrap();
        assert_eq!(
            args,
            [
                "library",
                "read",
                "--path",
                "synthetic_library.msl",
                "--limit",
                "2",
                "--offset",
                "1",
                "--out",
                "ions.tsv",
                "--fragments",
                "--index-only-threshold",
                "0"
            ]
        );
    }

    #[test]
    fn the_threshold_is_left_to_mzlib_when_not_given() {
        let args = read_request(
            Path::new("synthetic_library.msl"),
            &LibraryReadOptions {
                fragments: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(!args.iter().any(|a| a == "--index-only-threshold"));
        let library: LibraryRecords = recorded("library_read_synthetic.json");
        assert_eq!(library.index_only_threshold, 1_073_741_824, "mzLib's 1 GiB");
    }

    #[test]
    fn bad_arguments_are_refused_before_the_bridge_is_spawned() {
        assert!(read_request(Path::new(""), &LibraryReadOptions::default()).is_err());
        let blank_out = LibraryReadOptions {
            out: Some(" ".to_owned()),
            ..Default::default()
        };
        let error = read_request(Path::new("a.msl"), &blank_out).unwrap_err();
        assert!(error.to_string().contains(".tsv"));
        let empty: [&str; 0] = [];
        assert!(read_many_request(&empty, &LibraryBulkOptions::default()).is_err());
        let zero_threads = LibraryBulkOptions {
            bulk: BulkOptions {
                threads: 0,
                ..Default::default()
            },
            ..Default::default()
        };
        assert!(read_many_request(&["a.msl"], &zero_threads).is_err());
    }

    #[test]
    fn read_many_sends_the_list_on_stdin() {
        let options = LibraryBulkOptions {
            bulk: BulkOptions {
                threads: 4,
                on_error: OnError::Skip,
                ..Default::default()
            },
            fragments: true,
            ..Default::default()
        };
        let (args, stdin) =
            read_many_request(&["synthetic_library.msl", "myPrositLib.msl"], &options).unwrap();
        assert_eq!(
            args,
            [
                "library",
                "read",
                "--paths-stdin",
                "--threads",
                "4",
                "--on-error",
                "skip",
                "--fragments"
            ]
        );
        assert_eq!(stdin, "synthetic_library.msl\nmyPrositLib.msl\n");
    }

    #[test]
    fn the_envelope_is_projected_field_for_field() {
        let library: LibraryRecords = recorded("library_read_synthetic.json");
        assert_eq!(library.reader, "MslLibrary");
        assert_eq!(library.file_type, None);
        assert!(!library.is_compressed && !library.is_index_only);
        assert_eq!(
            (
                library.record_count,
                library.target_count,
                library.decoy_count,
                library.proteoform_count
            ),
            (6, 4, 2, 1)
        );
        assert_eq!((library.returned_count, library.row_count), (6, 6));
        assert_eq!(
            (library.rows_not_read, library.retention_time_unit),
            (None, None)
        );
        assert!(library.min_precursor_mz.unwrap() < library.max_precursor_mz.unwrap());
        assert!(!library.fragments_included && !library.truncated);
        assert!(library.excluded_fields.is_empty() && library.failed_fields.is_empty());
        assert!(!library.caveats.is_empty());
        assert_eq!(library.output, None);
    }

    #[test]
    fn rows_are_the_stored_order_not_mz_order() {
        let library: LibraryRecords = recorded("library_read_synthetic.json");
        let index = library.columns.integers("precursor_index").unwrap();
        assert_eq!(index, (0..6).map(Some).collect::<Vec<_>>());
        let mz = library.columns.floats("precursor_mz").unwrap();
        assert!(mz.windows(2).any(|w| w[0] > w[1]), "not sorted by m/z");
    }

    #[test]
    fn mzlibs_not_available_sentinels_arrive_as_none() {
        let library: LibraryRecords = recorded("library_read_synthetic.json");
        let c = &library.columns;
        assert_eq!(c.floats("ion_mobility").unwrap()[4], None);
        assert_eq!(c.floats("q_value").unwrap()[4], None);
        assert_eq!(c.integers("nce").unwrap()[4], None);
        assert_eq!(c.strings("protein_accession").unwrap()[4], None);
        assert_eq!(c.integers("nce").unwrap()[0], Some(28));
    }

    #[test]
    fn a_library_with_no_proteins_names_those_fields_absent() {
        let prosit: LibraryRecords = recorded("library_read_prosit.json");
        assert_eq!(
            prosit.absent_fields,
            [
                "ion_mobility",
                "protein_accession",
                "protein_name",
                "gene_name"
            ]
        );
        assert_eq!(
            prosit.columns.strings("protein_accession").unwrap(),
            [None, None, None, None, None]
        );
    }

    #[test]
    fn fragment_rows_repeat_the_precursor_and_keep_mzlibs_values() {
        let ions: LibraryRecords = recorded("library_read_fragments.json");
        assert!(ions.fragments_included);
        assert_eq!((ions.returned_count, ions.row_count), (1, 6));
        let c = &ions.columns;
        assert_eq!(
            c.strings("full_sequence").unwrap(),
            vec![Some("PEPTIDEK".to_owned()); 6]
        );
        assert_eq!(c.integers("fragment_count").unwrap(), vec![Some(6); 6]);
        assert_eq!(c.integers("residue_position").unwrap()[5], None);
        assert_eq!(c.floats("neutral_loss").unwrap()[0], Some(0.0));
        assert_eq!(
            c.strings("secondary_product_type").unwrap()[3].as_deref(),
            Some("y")
        );
        assert_eq!(c.booleans("is_internal").unwrap()[3], Some(true));
        assert_eq!(c.strings("secondary_product_type").unwrap()[0], None);
    }

    #[test]
    fn an_index_only_read_says_so_and_gives_the_same_rows() {
        let index_only: LibraryRecords = recorded("library_read_index_only.json");
        let whole: LibraryRecords = recorded("library_read_synthetic.json");
        assert!(index_only.is_index_only);
        assert_eq!(index_only.index_only_threshold, 0);
        assert_eq!(
            index_only.columns.strings("full_sequence").unwrap(),
            whole.columns.strings("full_sequence").unwrap()[..2]
        );
    }

    #[test]
    fn a_batch_projects_one_report_per_library() {
        let batch: LibraryBatch = recorded("library_read_many.json");
        assert_eq!(
            (batch.file_count, batch.read_count, batch.failed_count),
            (2, 2, 0)
        );
        assert_eq!(batch.on_error, OnError::Fail);
        assert_eq!(batch.record_count, 11);
        assert_eq!(
            batch
                .files
                .iter()
                .map(|f| f.record_count)
                .collect::<Vec<_>>(),
            [Some(6), Some(5)]
        );
        assert!(batch.files.iter().all(LibraryFileReport::ok));
        assert_eq!(batch.columns.names()[..2], ["source_index", "source_path"]);
        assert_eq!(batch.files[0].reader.as_deref(), Some("MslLibrary"));
    }

    #[test]
    fn a_skipped_library_keeps_its_unknown_facts_as_none() {
        let mut data: serde_json::Value = recorded("library_read_many.json");
        let entry = &mut data["files"][1];
        for key in [
            "reader",
            "format_version",
            "is_compressed",
            "is_index_only",
            "target_count",
            "decoy_count",
            "proteoform_count",
            "min_precursor_mz",
            "max_precursor_mz",
            "record_count",
        ] {
            entry[key] = serde_json::Value::Null;
        }
        entry["error"] =
            serde_json::json!({"kind": "usage", "type": "usage", "message": "File not found"});
        let batch: LibraryBatch = serde_json::from_value(data).unwrap();
        let failed = batch.failed_files();
        assert_eq!(failed.len(), 1);
        assert_eq!(
            failed[0].record_count, None,
            "not zero: nothing was counted"
        );
        assert_eq!(failed[0].target_count, None);
        assert_eq!(failed[0].error.as_ref().unwrap().kind, "usage");
    }

    #[test]
    fn a_table_written_to_out_returns_no_columns() {
        let mut data: serde_json::Value = recorded("library_read_fragments.json");
        data["columns"] = serde_json::Value::Null;
        data["returned_count"] = 0.into();
        data["row_count"] = 0.into();
        data["output"] =
            serde_json::json!({"path": "C:/x/ions.tsv", "format": "tsv", "row_count": 6});
        let written: LibraryRecords = serde_json::from_value(data).unwrap();
        assert!(written.columns.is_empty());
        assert_eq!(written.output.unwrap().row_count, 6);
    }
}
