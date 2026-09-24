//! Reading mass-spectrometry data files and proteomics search results: what a file *is*, and
//! every one of them read — one file at a time, or hundreds in one call.
//!
//! **Spectra files are read here, not just search output.** [`read_spectra`] reads **mzML**,
//! Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF and msalign — scan headers always, peaks opt-in —
//! and reports what the file records about the run it came from:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::readers::{read_spectra_with, ReadOptions, SpectraOptions};
//!
//! let scans = read_spectra_with(
//!     "sliced_ethcd.mzML",
//!     &SpectraOptions { read: ReadOptions { limit: Some(3), ..Default::default() }, ..Default::default() },
//! )?;
//! assert_eq!(scans.scan_count, 6);                  // the whole file
//! assert_eq!(scans.columns.rows(), 3);              // what came back
//! assert!(scans.truncated);                         // and it says so
//! let minutes = scans.columns.floats("retention_time")?;
//! assert_eq!(minutes[0], Some(38.92571663975));
//!
//! let source = scans.source.as_ref().expect("an mzML records its source");
//! assert_eq!(source.instrument_model.as_deref(), Some("Orbitrap Fusion"));
//! assert_eq!(source.instrument_serial_number.as_deref(), Some("FSN10189"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! mzLib 1.0.592 recognises **36 file types**: those instrument and deconvolution formats, plus
//! the output of a dozen search tools — MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT,
//! Crux, Casanovo, FlashDeconv, Dinosaur, DIA-NN, FlashLFQ, Pytheas, and any mzIdentML writer — and
//! maintains a parser for each. **Every one is readable here.** Ask [`formats`] rather than
//! trusting that number: it is enumerated from the mzLib the bridge carries, so it cannot drift.
//!
//! ```
//! # mzlib_replay::activate();
//! let info = mzlib::readers::identify("PXD078927_msgf_1_1_0.mzid")?;
//! assert_eq!(info.file_type, "MzIdentML");
//! assert_eq!(info.views, ["spectral_match"]);        // so read_matches reads it
//!
//! let table = mzlib::readers::read_records("ToppicPrsm_TopPICv1.6.2_prsm.tsv")?;
//! assert_eq!(table.record_type, "ToppicPrsm");
//! assert_eq!(table.columns.names().len(), 36);        // TopPIC's own fields
//! let e_values = table.columns.floats("e_value")?;    // Vec<Option<f64>>
//! # assert_eq!(e_values.len(), 4);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! ## Choosing a function
//!
//! What differs between formats is not *whether* you can read them but *what the columns mean*.
//! The counts are mzLib 1.0.592's; [`formats`] gives the live ones.
//!
//! | function | reads | columns |
//! |---|---|---|
//! | [`read_records`] | **all 36** | **that format's own fields**, under mzLib's names |
//! | [`read_results`] | 4 | uniform `quantifiable` view: sequence, RT, charge, mass, protein groups |
//! | [`read_features`] | 2 | uniform `ms1_features` view: m/z, charge, RT range, intensity |
//! | [`read_matches`] | 6 | uniform `spectral_match` view: scan, sequences, accession, decoy flag, q-value, rank |
//! | [`read_spectra`] | 7 | scan headers; peaks opt-in; the run's instrument and start time |
//! | [`read_protein_groups`] | 1 | MetaMorpheus protein groups, **long**: one row per group per sample group |
//! | [`read_quantified_peptides`] | 1 | FlashLFQ peptides, **long**: one row per peptide per sample |
//! | [`read_occupancy`] | 1 | MetaMorpheus PTM site occupancy: one row per group, sample group, basis and site |
//!
//! The rule of thumb: **a typed view when you need numbers that mean the same thing across files,
//! and [`read_records`] when you need everything one file has.** A `.psmtsv` through
//! [`read_results`] gives 10 comparable columns; the same file through [`read_records`] gives 73,
//! including the q-values and scores the uniform view does not carry.
//!
//! An empty [`FileInfo::views`] is a real and common answer — **17 of the 36** have it, meaning
//! mzLib parses the file into a shape that shares nothing with any other format. Those are exactly
//! the files [`read_records`] exists for. Three of them — the MetaMorpheus and FlashLFQ
//! quantification tables — also have functions of their own, because their per-sample values are
//! dictionaries that [`read_records`] cannot project (mzLib #1347): it names them in
//! [`NativeRecords::excluded_fields`] and points at the function that carries them.
//!
//! ## Everything comes back as a [`Table`]
//!
//! One array per column rather than one struct per row, because the column set is not knowable at
//! compile time: it depends on the format. [`Table`] gives typed accessors that project a wire
//! `null` onto [`Option`], so a missing cell can never silently become a zero and can never
//! shorten a column.
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::readers::{read_records_with, ReadOptions};
//!
//! let t = read_records_with("crux.txt", &ReadOptions { limit: Some(3), ..Default::default() })?;
//! for (sequence, score) in t
//!     .columns
//!     .strings("base_sequence")?
//!     .iter()
//!     .zip(t.columns.floats("x_corr_score")?)
//! {
//!     if let (Some(sequence), Some(score)) = (sequence, score) {
//!         println!("{sequence}\t{score}");
//!     }
//! }
//! assert_eq!((t.record_count, t.returned_count), (14, 3));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! ## Four ways a field can have no value, and each is named
//!
//! A `None` in a table can mean four different things, and conflating them is how a missing column
//! turns into a measured zero. Every result keeps them apart:
//!
//! - **`absent_fields`** — the function defines the field, but *this file's format has no column
//!   for it*: MBR Score in a current FlashLFQ peaks table (mzLib #1345), `q_value` in an
//!   MsPathFinderT targets file, apex intensity in a FLASHDeconv feature file. Every value in the
//!   column is `None`, whatever default mzLib filled in.
//! - **`failed_fields`** — the field exists, but *reading it threw* for some rows (a UniProt-shaped
//!   accession parsed from a plain FASTA header). `None` in those rows, named `"field: Exception"`.
//! - **`excluded_fields`** — the field has *no column shape* (a dictionary, a nested object), so it
//!   does not cross in this function at all; each [`ExcludedField::verb`] names the command that
//!   carries it.
//! - a `None` in one row, with the field in none of those lists — the value is *genuinely missing
//!   for that row*: the precursor of an MS1 scan, a blank intensity cell.
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::readers::{read_results_with, ReadOptions};
//!
//! let psms = read_results_with(
//!     "FraggerPsm_FragPipev21.1_psm.tsv",
//!     &ReadOptions { limit: Some(2), ..Default::default() },
//! )?;
//! // MSFragger writes no target/decoy label: the column is absent, so every value is unknown.
//! assert_eq!(psms.absent_fields, ["is_decoy"]);
//! assert_eq!(psms.columns.booleans("is_decoy")?, [None, None]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! ## Many files in one call
//!
//! Every reader has a `_many` twin — [`read_spectra_many`], [`read_records_many`],
//! [`identify_many`] and the rest — that takes a list of paths and returns ONE long table
//! ([`ReadBatch`]) whose first two columns say which file each row came from (`source_index`,
//! `source_path`), plus one [`FileReport`] per file. The list is read by **one** bridge process,
//! [`BulkOptions::threads`] files at a time, so two hundred runs cost one .NET start-up rather than
//! two hundred. The answer is identical at any thread count; the default is 1 because each file in
//! flight is held whole in memory. There is deliberately no Rust-side thread pool or loop over files
//! anywhere in this module: parallelism is the bridge's to express, once, where it can be counted.
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::readers::{read_spectra_many, BulkOptions, SpectraBulkOptions};
//! use mzlib::OnError;
//!
//! let batch = read_spectra_many(
//!     &["sliced_ethcd.mzML", "no-such-run.mzML", "withZeros.mgf"],
//!     &SpectraBulkOptions {
//!         bulk: BulkOptions { threads: 2, on_error: OnError::Skip, ..Default::default() },
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!((batch.file_count, batch.read_count, batch.failed_count), (3, 2, 1));
//! assert_eq!(batch.columns.names()[..2], ["source_index", "source_path"]);
//! let missing = &batch.failed_files()[0];
//! assert_eq!(missing.error.as_ref().unwrap().kind, "usage");       // file not found, skipped
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! ## Units are not normalised across formats
//!
//! mzLib's result-file readers pass through whatever the writing tool wrote. MetaMorpheus and
//! MSFragger retention times are minutes; TopFD's `_ms1.feature` changed from seconds to minutes
//! at v1.7.0 without changing the file type. Every typed view therefore reports a
//! `retention_time_unit`, and the `*_in_minutes` methods convert — or refuse, when mzLib gives no
//! basis to say. [`ReadBatch::in_minutes`] converts each file's rows by that file's own unit.
//! Spectra readers are the exception: they convert to minutes at the boundary.
//!
//! **Nothing here is FDR-filtered.** Every result format records confidence somewhere, and none of
//! these functions filters on it. [`read_matches`] carries mzIdentML's q-value, rank and threshold
//! as columns; filter before you report.

use std::collections::BTreeMap;
use std::path::Path;
use std::time::Duration;

use serde::Deserialize;
use serde_json::Value;

use crate::bridge::{self, MzLibError, Result};

pub use crate::bridge::OnError;

/// The view name for the cross-format record shape [`crate::flashlfq::quantify`] consumes.
pub const QUANTIFIABLE: &str = "quantifiable";

/// The view name for deconvolved MS1 features — [`read_features`].
pub const MS1_FEATURES: &str = "ms1_features";

/// The view name for records that are identifications — [`read_matches`].
pub const SPECTRAL_MATCH: &str = "spectral_match";

/// The view name for files that are spectra rather than results — [`read_spectra`].
pub const SPECTRA: &str = "spectra";

// ---------------------------------------------------------------------------------------------
// The columnar table
// ---------------------------------------------------------------------------------------------

/// A columnar table: one array per field, with typed accessors.
///
/// The column set depends on the format — a TopPIC file has 36 columns and a Crux file 23 — so
/// this cannot be a struct with named fields. What it can be is a map whose accessors do the
/// projection properly: every cell is an [`Option`], a wire `null` becomes [`None`], and a column
/// whose values are not the type you asked for is an error rather than a silent default.
///
/// It carries the wire's `column_names` as [`Table::names`], so the order mzLib declares is kept.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Table {
    names: Vec<String>,
    columns: BTreeMap<String, Vec<Value>>,
}

impl Table {
    /// The column names, **in the order mzLib declares them** — base-class fields first, then
    /// each format's own. The wire's `column_names`.
    #[must_use]
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// The number of rows, or zero when the table went to disk instead.
    #[must_use]
    pub fn rows(&self) -> usize {
        self.names
            .first()
            .and_then(|name| self.columns.get(name))
            .map_or(0, Vec::len)
    }

    /// Whether the table carries any rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.rows() == 0
    }

    /// Whether a column is present.
    #[must_use]
    pub fn has(&self, name: &str) -> bool {
        self.columns.contains_key(name)
    }

    /// A column's raw JSON values, for a shape the typed accessors do not cover.
    #[must_use]
    pub fn raw(&self, name: &str) -> Option<&[Value]> {
        self.columns.get(name).map(Vec::as_slice)
    }

    /// A column as floating-point numbers, with `null` as [`None`].
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if the column is absent, [`MzLibError::Protocol`] if a value is
    /// present but is not a number.
    pub fn floats(&self, name: &str) -> Result<Vec<Option<f64>>> {
        self.project(name, "a number", |value| match value {
            Value::Null => Some(None),
            Value::Number(number) => number.as_f64().map(Some),
            _ => None,
        })
    }

    /// A column as whole numbers, with `null` as [`None`].
    ///
    /// # Errors
    ///
    /// As [`Table::floats`], and also when a value is a number with a fractional part — silently
    /// truncating one would turn a protocol change into wrong data.
    pub fn integers(&self, name: &str) -> Result<Vec<Option<i64>>> {
        self.project(name, "a whole number", |value| match value {
            Value::Null => Some(None),
            Value::Number(number) => number.as_i64().map(Some),
            _ => None,
        })
    }

    /// A column as strings, with `null` as [`None`].
    ///
    /// # Errors
    ///
    /// As [`Table::floats`], for a value that is not a string.
    pub fn strings(&self, name: &str) -> Result<Vec<Option<String>>> {
        self.project(name, "a string", |value| match value {
            Value::Null => Some(None),
            Value::String(text) => Some(Some(text.clone())),
            _ => None,
        })
    }

    /// A column as booleans, with `null` as [`None`].
    ///
    /// `None` genuinely means *unknown* in this library rather than *false* — MSFragger's
    /// `is_decoy` is the case that matters, because its `psm.tsv` has no target/decoy label at all
    /// and a `false` there would be a fabricated value someone could filter on.
    ///
    /// # Errors
    ///
    /// As [`Table::floats`], for a value that is not a boolean.
    pub fn booleans(&self, name: &str) -> Result<Vec<Option<bool>>> {
        self.project(name, "a boolean", |value| match value {
            Value::Null => Some(None),
            Value::Bool(flag) => Some(Some(*flag)),
            _ => None,
        })
    }

    /// A column whose every cell is itself an array of numbers.
    ///
    /// The shape [`read_spectra`] returns for `mz` and `intensity` under
    /// [`SpectraOptions::peaks`]: one array per scan, not one number per scan.
    ///
    /// # Errors
    ///
    /// As [`Table::floats`], for a value that is not an array of numbers.
    pub fn float_arrays(&self, name: &str) -> Result<Vec<Option<Vec<Option<f64>>>>> {
        self.project(name, "an array of numbers", |value| match value {
            Value::Null => Some(None),
            Value::Array(items) => items
                .iter()
                .map(|item| match item {
                    Value::Null => Some(None),
                    Value::Number(number) => number.as_f64().map(Some),
                    _ => None,
                })
                .collect::<Option<Vec<_>>>()
                .map(Some),
            _ => None,
        })
    }

    /// The shared body of every typed accessor: absent is a usage error, wrong-typed is a protocol
    /// error, and the two are never conflated.
    fn project<T>(
        &self,
        name: &str,
        expected: &str,
        convert: impl Fn(&Value) -> Option<T>,
    ) -> Result<Vec<T>> {
        let column = self.columns.get(name).ok_or_else(|| {
            // Naming what IS there, because the column set is per-format and a caller who guessed
            // a name from another format has no other way to find out.
            MzLibError::Usage(format!(
                "No column '{name}' in this table. Its columns are: {}.",
                self.names.join(", ")
            ))
        })?;

        column
            .iter()
            .enumerate()
            .map(|(row, value)| {
                convert(value).ok_or_else(|| {
                    MzLibError::Protocol(format!(
                        "Column '{name}' row {row} is not {expected}: {value}"
                    ))
                })
            })
            .collect()
    }

    /// Builds a table from the wire's `column_names` and `columns`.
    ///
    /// Order comes from `column_names` rather than from the map, because a JSON object has no
    /// order a parser must preserve — reading it back from the map would give alphabetical
    /// columns, which is not the order mzLib declares its fields in.
    fn from_wire(names: Vec<String>, columns: Option<BTreeMap<String, Vec<Value>>>) -> Self {
        let columns = columns.unwrap_or_default();
        // Only names that actually have a column, so `names()` cannot promise a column `raw()`
        // then fails to return.
        let names = if names.is_empty() {
            columns.keys().cloned().collect()
        } else {
            names
                .into_iter()
                .filter(|name| columns.contains_key(name))
                .collect()
        };
        Self { names, columns }
    }
}

/// A table deserializes from the two wire keys that describe it, `column_names` and `columns`, so
/// a result type can hold one with `#[serde(flatten)]` and read the rest of the payload itself.
impl<'de> Deserialize<'de> for Table {
    fn deserialize<D: serde::Deserializer<'de>>(
        deserializer: D,
    ) -> std::result::Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct Wire {
            #[serde(default, deserialize_with = "bridge::null_to_default")]
            column_names: Vec<String>,
            #[serde(default)]
            columns: Option<BTreeMap<String, Vec<Value>>>,
        }
        let wire = Wire::deserialize(deserializer)?;
        Ok(Self::from_wire(wire.column_names, wire.columns))
    }
}

// ---------------------------------------------------------------------------------------------
// Wire types
// ---------------------------------------------------------------------------------------------

/// One file type mzLib can recognise: an entry of [`formats`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct Format {
    /// mzLib's `SupportedFileType` member name, e.g. `"MsFraggerPsm"`, `"psmtsv"`,
    /// `"MzIdentMLGz"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The extension or filename suffix mzLib dispatches on, e.g. `".mzid.gz"` or
    /// `"QuantifiedProteinGroups.tsv"`.
    ///
    /// **Not unique across file types.** `BrukerD` and `BrukerTimsTof` are both `.d`, told apart
    /// by which analysis file the directory holds, and several formats share `.tsv`; `DiaNnReport`
    /// and `PytheasResult` dispatch on content. `None` only for a broken mzLib build with no
    /// extension mapping for the member — the listing reports it rather than failing.
    #[serde(default)]
    pub extension: Option<String>,
    /// The name of the mzLib class that parses it, for cross-referencing the mzLib source. `None`
    /// only for a broken mzLib build with no reader mapping for the member.
    #[serde(default)]
    pub reader: Option<String>,
    /// The cross-format views this format offers, from `quantifiable`, `ms1_features`,
    /// `spectra` and `spectral_match`. **Empty is a real answer**: 17 of 36 at mzLib 1.0.592
    /// offer none and are readable only through [`read_records`].
    #[serde(default)]
    pub views: Vec<String>,
}

impl Format {
    /// Whether this format offers the cross-format record view, and so feeds FlashLFQ.
    #[must_use]
    pub fn is_quantifiable(&self) -> bool {
        self.views.iter().any(|view| view == QUANTIFIABLE)
    }
}

/// Why one input of a `_many` call could not be read, under [`OnError::Skip`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct ReadError {
    /// `"usage"` — the input was wrong: not found, not a type mzLib recognises, or not one this
    /// function reads; `"correctness"` — mzLib recognised the file and failed to parse it. The same
    /// split every [`MzLibError`] makes (`Usage` against `Bridge`).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub kind: String,
    /// What the failure would have been on its own: `"usage"`, or the .NET exception type, e.g.
    /// `"MzLibException"`, `"HeaderValidationException"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub r#type: String,
    /// mzLib's own message.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub message: String,
}

/// What a particular file is, and what can be done with it: what [`identify`] returns, and one
/// entry of [`IdentifyBatch::files`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct FileInfo {
    /// The absolute path that was identified.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// mzLib's `SupportedFileType` name. Empty for an input of [`identify_many`] that could not be
    /// identified under [`OnError::Skip`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The extension mzLib dispatched on. `None` when mzLib maps no extension to this type, or the
    /// input could not be identified.
    #[serde(default)]
    pub extension: Option<String>,
    /// The mzLib reader class that would parse it. `None` for an input that could not be
    /// identified.
    #[serde(default)]
    pub reader: Option<String>,
    /// The uniform views this file supports. **Empty is a real answer**, and means mzLib can read
    /// the file but offers no cross-format projection of it — use [`read_records`].
    #[serde(default)]
    pub views: Vec<String>,
    /// Always `None` from [`identify`], which returns an error instead. From [`identify_many`] under
    /// [`OnError::Skip`], why this path could not be identified.
    #[serde(default)]
    pub error: Option<ReadError>,
}

impl FileInfo {
    /// Whether this file offers the cross-format record view.
    ///
    /// `true` is the precondition for [`crate::flashlfq::quantify`] — but it is **not
    /// permission**. It reports what mzLib's *interface* offers, not that the numbers are
    /// comparable. See [`read_results`] and the caveats it returns.
    #[must_use]
    pub fn is_quantifiable(&self) -> bool {
        self.views.iter().any(|view| view == QUANTIFIABLE)
    }

    /// Whether this file offers a given view, e.g. [`MS1_FEATURES`].
    #[must_use]
    pub fn has_view(&self, view: &str) -> bool {
        self.views.iter().any(|present| present == view)
    }
}

/// Where a read wrote its table, when asked to write one instead of returning it.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct WrittenTable {
    /// The absolute path written.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Always `"tsv"`. **Tab-separated, not comma-separated**, because these fields contain
    /// commas — MSFragger's mapped proteins are a comma-separated list inside one field.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub format: String,
    /// Rows written, excluding the header: the window that was selected, not the whole file. For a
    /// `_many` read, every file's rows, written one file at a time.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
}

/// A field of a record type that could not become a column, why, and where it is carried instead.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct ExcludedField {
    /// The field's wire name, as it would have been.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub field: String,
    /// Its .NET type, e.g. `"List<AlternativeToppicId>"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub r#type: String,
    /// Why it could not cross the wire: a nested object, a list of composites, or a dictionary.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub reason: String,
    /// The wire verb that **does** carry this field, or `None` if none does: `sample_groups` points
    /// at `readers read-protein-groups` ([`read_protein_groups`]), `samples` at
    /// `readers read-quantified-peptides` ([`read_quantified_peptides`]), mzIdentML's `scores` at
    /// `readers read-matches` ([`read_matches_with`] with [`MatchOptions::scores`]).
    #[serde(default)]
    pub verb: Option<String>,
}

/// An mzIdentML identification item mzLib did not turn into a record, and why (mzLib #1313).
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct SkippedMatch {
    /// The item's `id` in the document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub spectrum_identification_item_id: String,
    /// The spectrum's nativeID, as written.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub spectrum_id: String,
    /// Why it was skipped: a crosslink identification, a modification that does not resolve to a
    /// Unimod entry, a substitution.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub reason: String,
}

/// What a spectra file records about the run it came from: mzLib's `SourceFile` (#1349).
///
/// The facts a batch or instrument confound is built from — which instrument, which unit of that
/// model, and when. Each is `None` when the file does not record it; MGF and msalign record none.
#[derive(Debug, Clone, PartialEq, Eq, Default, Deserialize)]
pub struct SpectraSource {
    /// The instrument model's **name**, e.g. `"Orbitrap Fusion Lumos"`. `None`: the file does not
    /// record it.
    #[serde(default)]
    pub instrument_model: Option<String>,
    /// The model's PSI-MS accession, e.g. `"MS:1002416"`. mzML carries one; Thermo `.raw` records
    /// only the name, so it is `None` there rather than looked up. **Match on the accession, never
    /// the name.**
    #[serde(default)]
    pub instrument_model_accession: Option<String>,
    /// The instrument's serial number, verbatim, e.g. `"FSN10189"`. With the model it tells two
    /// instruments of one model apart. Free text: a converter's placeholder such as
    /// `"Serial Number N/A"` is reported as written. `None`: the file does not record it.
    #[serde(default)]
    pub instrument_serial_number: Option<String>,
    /// When acquisition started, as ISO-8601 text, e.g. `"2021-03-16T17:09:07Z"`. It ends in `Z`
    /// only when the source fixed the instant (an mzML `startTimeStamp` with an offset); otherwise
    /// it is the acquisition computer's wall-clock time with no offset, because the instant is
    /// unknown. `None`: the file does not record it. [`Self::acquired_at`] parses it.
    #[serde(default)]
    pub acquisition_start_time: Option<String>,
    /// `true` when [`Self::acquisition_start_time`] is a UTC instant; `false` when it is local
    /// instrument-PC time — **always, for Thermo `.raw`** — or absent. A `.raw` and ProteoWizard's
    /// mzML of the same run can differ by the site's UTC offset, because ProteoWizard assumes the
    /// converting machine's time zone.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub acquisition_start_time_is_utc: bool,
}

/// When acquisition started, typed by what the file actually knows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AcquisitionTime {
    /// A fixed instant: the file recorded an offset.
    Utc(chrono::DateTime<chrono::Utc>),
    /// A reading of the instrument PC's clock, with **no known offset** — which a timezone-aware
    /// type would have to invent.
    Local(chrono::NaiveDateTime),
}

impl SpectraSource {
    /// [`Self::acquisition_start_time`] parsed: an instant when the file fixed one, a local clock
    /// reading when it did not, `None` when it recorded nothing or the text does not parse.
    #[must_use]
    pub fn acquired_at(&self) -> Option<AcquisitionTime> {
        let text = self.acquisition_start_time.as_deref()?;
        if let Ok(instant) = chrono::DateTime::parse_from_rfc3339(text) {
            return Some(AcquisitionTime::Utc(instant.with_timezone(&chrono::Utc)));
        }
        ["%Y-%m-%dT%H:%M:%S%.f", "%Y-%m-%dT%H:%M:%S"]
            .iter()
            .find_map(|format| chrono::NaiveDateTime::parse_from_str(text, format).ok())
            .map(AcquisitionTime::Local)
    }
}

// ---------------------------------------------------------------------------------------------
// One file's result
// ---------------------------------------------------------------------------------------------

/// The uniform `quantifiable` view of a result file — what [`read_results`] returns.
#[derive(Debug, Clone, Deserialize)]
pub struct ResultRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// The mzLib `SupportedFileType` that was dispatched.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The mzLib reader class that parsed the file, e.g. `"PsmFromTsvFile"`.
    #[serde(default)]
    pub reader: Option<String>,
    /// Records in the **whole file**, before the window — whatever [`ReadOptions::limit`] and
    /// [`ReadOptions::offset`] selected.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Records carried back in [`Self::columns`], the unit the window counts in. Zero when the
    /// table was written to disk instead.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// The offset applied, in records.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether records were left behind**, by either the limit or the offset. A short answer and
    /// a complete one must never look alike, so check this rather than comparing counts yourself.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// The unit [`Self::columns`]' `retention_time` carries for this format: `"minutes"`,
    /// `"seconds"`, or `"unknown"`. mzLib does not normalise it, so it differs per format.
    /// Convert with [`ResultRecords::retention_time_in_minutes`] rather than by hand.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub retention_time_unit: String,
    /// Data rows that did not become records. mzLib drops a malformed row **silently**, so a
    /// non-zero value here means the file is partly unreadable and the table is incomplete.
    /// `None` when the count is not meaningful for this format: only the psmtsv family and
    /// MSFragger are one line per record.
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// **What the uniform view cannot be trusted to mean for this file**, each citing the mzLib
    /// source it came from. Worth reading before comparing anything across formats.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Columns this view defines but **this file's format has no column for**, so every value is
    /// `None`: `["is_decoy"]` for MSFragger and DIA-NN, which write no target/decoy label.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"field: ExceptionType"` for each column whose read threw on some returned rows; those
    /// cells are `None`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Fields with no column shape. Always empty for this view.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// The table, one row per record: `file_name`, `base_sequence`, `full_sequence`,
    /// `retention_time`, `charge_state`, `monoisotopic_mass`, `is_decoy`, `protein_accessions`,
    /// `gene_name`, `organism`. Empty when the records went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
}

/// Every field of one format, whatever format it is — what [`read_records`] returns.
///
/// The columns here are **not uniform**: they are this format's own mzLib record fields, under
/// mzLib's own names in `snake_case`, which makes them cross-referenceable against the mzLib
/// source. A column called `e_value` is `ToppicPrsm.EValue`, and [`Self::record_type`] names the
/// class to look in.
#[derive(Debug, Clone, Deserialize)]
pub struct NativeRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// The `SupportedFileType` mzLib dispatched — not a guess from the extension.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The mzLib reader class that parsed the file.
    #[serde(default)]
    pub reader: Option<String>,
    /// The mzLib record class the columns are properties of, e.g. `"ToppicPrsm"`,
    /// `"MzIdentMLRecord"`, `"ProteinGroupFromTsv"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_type: String,
    /// The cross-format views this file *also* offers, if any. Often empty.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub views: Vec<String>,
    /// Records mzLib parsed from the **whole file**, before the window. Records, not lines: what
    /// mzLib dropped is reported beside it, in [`Self::rows_not_read`] and [`Self::skipped`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Records carried back in [`Self::columns`]. Zero when the table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// The offset applied, in records.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether records were left behind**, by either the limit or the offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// **Fields that could not become columns**, each with the reason and the verb that carries it.
    /// A nested object or a dictionary has no faithful column shape, and inventing one would mean
    /// publishing a schema mzLib does not have — so they are named rather than dropped, because a
    /// column that simply vanished is indistinguishable from a field the format does not have.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// Fields that **raised** while being read, as `"field: ExceptionType"`. Several mzLib
    /// properties are computed and assume a UniProt-style FASTA header — Crux's and
    /// MsPathFinderT's `accession` are both `protein_id` split on `|` — so on other databases
    /// they throw. Those cells arrive as `None` rather than failing the whole read, but a
    /// failure must not look like missing data.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Columns **this file has no source for**: the record type reads them from an optional column
    /// the file lacks, so every value is mzLib's default and crosses as `None` here —
    /// `["mbr_score"]` on a current FlashLFQ peaks table (mzLib #1345), `["q_value", "pep_q_value"]`
    /// on an MsPathFinderT targets file, where mzLib would otherwise report 0. Empty means none
    /// found, or no basis to judge (the psmtsv family parses its own header).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// Data rows that did not become records, counted where one line is one record (the psmtsv
    /// family and MSFragger); mzLib drops a malformed line silently. `None`: the count is not
    /// meaningful for this format.
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// Always `None` for this function: the columns are the format's own, and several formats carry
    /// more than one time column in different units. Use a typed view when units matter.
    #[serde(default)]
    pub retention_time_unit: Option<String>,
    /// Always empty for this function; the caveats in its reference apply.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// mzIdentML items mzLib did not represent as records — crosslinks, unresolvable
    /// modifications, substitutions (mzLib #1313) — so `record_count + skipped_count` is the items
    /// in the file. `None` for every other format, which keeps no skip list.
    #[serde(default)]
    pub skipped_count: Option<u64>,
    /// The same items, each with its reason. `None` with [`Self::skipped_count`].
    #[serde(default)]
    pub skipped: Option<Vec<SkippedMatch>>,
    /// The table: the record type's projectable public properties in `snake_case`, base class
    /// first, in declaration order. Empty when the records went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline. Its `row_count` is the
    /// window written, not the whole file.
    #[serde(default)]
    pub output: Option<WrittenTable>,
}

/// Deconvolved MS1 features, in the uniform `ms1_features` view — what [`read_features`] returns.
#[derive(Debug, Clone, Deserialize)]
pub struct FeatureRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// The mzLib `SupportedFileType` that was dispatched.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The mzLib reader class that parsed the file, e.g. `"Ms1FeatureFile"`.
    #[serde(default)]
    pub reader: Option<String>,
    /// Features in the **whole file**, before the window. **For `_ms1.feature` this exceeds the
    /// file's line count**: mzLib expands each deconvolved feature into one feature per charge.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Features carried back in [`Self::columns`], the unit the window counts in. Zero when the
    /// table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// The offset applied, in features.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether features were left behind**, by either the limit or the offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// `"minutes"` for Dinosaur, and `"unknown"` for `_ms1.feature`.
    ///
    /// **It is genuinely `"unknown"` for `_ms1.feature`.** TopFD wrote seconds through v1.6.2
    /// and minutes from v1.7.0 without changing the file type, and mzLib normalises neither.
    /// That is not a gap in this crate; it is the honest state of the format.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub retention_time_unit: String,
    /// Always `None`: rows that did not become records are never counted here, because one line of
    /// an `_ms1.feature` is not one feature (rows are expanded per charge).
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// What this view cannot be trusted to mean for this file, each citing the mzLib source.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Columns this file has no source for, so every value is `None`: `number_of_isotopes` for
    /// every `_ms1.feature` (mzLib never sets it there), and `intensity` for a FLASHDeconv file,
    /// whose schema has no apex intensity — mzLib would otherwise report zero for every feature.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"field: ExceptionType"` for each column whose read threw on some returned rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Fields with no column shape. Always empty for this view.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// The table, one row per single-charge feature: `mz`, `charge`, `retention_time_start`,
    /// `retention_time_end`, `intensity`, `number_of_isotopes`. Empty when it went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
}

/// Identifications, in the uniform `spectral_match` view — what [`read_matches`] returns.
///
/// **Nothing here is FDR-filtered.** `q_value` is filled only by mzIdentML and by an MsPathFinderT
/// file with a `QValue` column; [`Self::absent_fields`] names it for every other file. mzIdentML
/// also lists every candidate, lower ranks and failed thresholds included: filter on `rank == 1`
/// and `pass_threshold` before counting.
#[derive(Debug, Clone, Deserialize)]
pub struct MatchRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// The mzLib `SupportedFileType` that was dispatched.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The mzLib reader class that parsed the file, e.g. `"MzIdentMLResultFile"`.
    #[serde(default)]
    pub reader: Option<String>,
    /// Matches in the **whole file**, before the window.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Matches carried back in [`Self::columns`], the unit the window counts in. Zero when the
    /// table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// Rows in [`Self::columns`]: one per match, or with [`MatchOptions::scores`] one per match and
    /// score, so it can exceed [`Self::returned_count`]. Zero when the table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The offset applied, in matches.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether matches were left behind**, by either the limit or the offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// Whether [`MatchOptions::scores`] made the table long by score (`match_index`, `score_name`,
    /// `score_value`).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub scores_included: bool,
    /// mzIdentML only: identification items mzLib did not turn into rows — crosslinks,
    /// unresolvable modifications, substitutions (mzLib #1313) — so `record_count +
    /// skipped_count` is the items in the file. `None` for every other format, which keeps no
    /// such list.
    #[serde(default)]
    pub skipped_count: Option<u64>,
    /// The same items, each with its reason. `None` with [`Self::skipped_count`].
    #[serde(default)]
    pub skipped: Option<Vec<SkippedMatch>>,
    /// Always `None`: rows that did not become records are not counted for this view.
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// Always `None`: this view has no time column.
    #[serde(default)]
    pub retention_time_unit: Option<String>,
    /// What this view cannot be trusted to mean for this file — that MsPathFinderT infers decoys
    /// from an `XXX` name prefix, that Casanovo's scan numbers are mzTab indices, that mzIdentML
    /// lists every candidate rank.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Columns this file's format has no source for, so every value is `None`: `is_decoy` for
    /// every format but MsPathFinderT; `q_value` for Casanovo and an MsPathFinderT targets file
    /// (mzLib would otherwise report 0, a perfect q-value); `rank` and `pass_threshold` for
    /// everything but mzIdentML; with scores, `score_name` and `score_value` for everything but
    /// mzIdentML.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"field: ExceptionType"` for each column whose read threw on some rows — MsPathFinderT's
    /// `accession` on a FASTA without UniProt-style headers.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Fields with no column shape. Always empty for this view.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// The table, one row per match (or per match and score): `file_name_without_extension`,
    /// `one_based_scan_number`, `base_sequence`, `full_sequence`, `accession`, `is_decoy`,
    /// `modifications`, `modification_count`, `q_value`, `rank`, `pass_threshold`. Empty when it
    /// went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
}

/// Scan headers, and optionally peaks — what [`read_spectra`] returns.
///
/// Retention times here **are** minutes for every format: mzLib's spectra readers convert at the
/// boundary, unlike its result-file readers.
#[derive(Debug, Clone, Deserialize)]
pub struct ScanRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// The mzLib `SupportedFileType` name.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The mzLib `MsDataFile` class that read it, e.g. `"Mzml"`, `"ThermoRawFileReader"`.
    #[serde(default)]
    pub reader: Option<String>,
    /// Scans in the **whole file**, before any MS-level filter. Reported alongside
    /// [`Self::record_count`] so a filter that matched nothing can never look like an empty file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub scan_count: u64,
    /// The MS level filtered to, or `None` when no MS-level filter was applied.
    #[serde(default)]
    pub ms_order: Option<i64>,
    /// Scans that passed the MS-level filter, before the window.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Scans carried back in [`Self::columns`]. Zero when the table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// The offset applied, in scans, counted after the MS-level filter.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether scans were left behind**, by either the limit or the offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// Whether `mz` and `intensity` are present — read them with [`Table::float_arrays`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub peaks_included: bool,
    /// Always `"minutes"` for this view: every mzLib `MsDataFile` reader converts at the boundary.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub retention_time_unit: String,
    /// What the file records about the run — instrument model and serial number, and when
    /// acquisition started (mzLib #1349). `None` only when the reader built no description at all;
    /// a fact the file does not record is `None` inside it.
    #[serde(default)]
    pub source: Option<SpectraSource>,
    /// Always `None`: rows that did not become records are never counted for this view.
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// Format-specific traps for this file — that msalign holds deconvolved neutral masses rather
    /// than m/z, that MGF scan numbers may be file order, that Bruker needs Windows x64.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Always empty for this view: a scan field a file does not record is `None` per scan, as for
    /// the precursor of an MS1 scan.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"field: ExceptionType"` for each column whose read threw on some scans; those cells are
    /// `None`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Always empty for this view.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// The table, one row per scan. With [`SpectraOptions::peaks`], `mz` and `intensity` are each
    /// one array per scan. Empty when it went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
}

/// A MetaMorpheus protein-group table in **long** form: one row per group per sample group —
/// what [`read_protein_groups`] returns.
///
/// Columns: `protein_group_name`, `gene`, `organism`, `decoy_contaminant_target`, `is_decoy`,
/// `is_contaminant`, `is_entrapment`, `q_value`, then `sample_label`, `spectral_count` and
/// `intensity`.
///
/// **Unfiltered.** Decoys, contaminants and groups above 1% FDR are all rows, as MetaMorpheus
/// writes them. Filter on `q_value` and `decoy_contaminant_target` before counting.
#[derive(Debug, Clone, Deserialize)]
pub struct ProteinGroupRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `"MetaMorpheusQuantifiedProteinGroups"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// `"ProteinGroupFromTsvFile"`.
    #[serde(default)]
    pub reader: Option<String>,
    /// The file's sample-group labels, in header order, verbatim — e.g.
    /// `"QE-002106_GM1_a-calib"`. Condition and replicate cannot be recovered from them; map them
    /// to your design yourself.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_labels: Vec<String>,
    /// Protein groups in the whole file — groups, not rows — before the window.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Groups carried back in [`Self::columns`] — groups, not rows; the unit the window counts in.
    /// Zero when the table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// Rows in [`Self::columns`]: returned groups times sample groups. Zero when the table was
    /// written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The offset applied, in groups.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether groups were left behind**, by either the limit or the offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// Always `None`: rows that did not become records are not counted for this table — one line
    /// is one group, and every group is read or the file fails.
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// Always `None`: this table has no time column.
    #[serde(default)]
    pub retention_time_unit: Option<String>,
    /// What this table cannot be trusted to mean, each citing the mzLib source.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Columns this file has no column for, so every value is `None`: `gene` or `organism` when
    /// MetaMorpheus did not write them; `sample_label`, `spectral_count` and `intensity` when the
    /// file has no per-sample columns, in which case each group is one row.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"field: ExceptionType"` for each column whose read threw on some rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// `count_occupancy` and `intensity_occupancy`, whose [`ExcludedField::verb`] is
    /// `readers read-occupancy`: [`read_occupancy`] has them, one row per site.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// The table. `intensity` is `None` where the cell was blank — not quantified in that sample
    /// group, which is not zero. Empty when it went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
}

/// A FlashLFQ peptide table in **long** form: one row per peptide per sample — what
/// [`read_quantified_peptides`] returns.
///
/// Columns: `sequence` (the full, modified sequence), `base_sequence`, `peak_order`,
/// `protein_groups`, `gene_names`, `organism`, then `sample_label`, `intensity`, `detection_type`
/// and `retention_time`.
///
/// **An intensity of 0 is not a measurement.** FlashLFQ writes a literal 0 for a peptide it did
/// not quantify in a sample; `detection_type` (`"MSMS"`, `"MBR"`, `"NotDetected"`, …) is what tells
/// the two apart.
#[derive(Debug, Clone, Deserialize)]
pub struct QuantifiedPeptideRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `"FlashLFQQuantifiedPeptide"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// `"QuantifiedPeptideFile"`.
    #[serde(default)]
    pub reader: Option<String>,
    /// The file's sample labels, in header order, verbatim.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_labels: Vec<String>,
    /// Peptides in the whole file — peptides, not rows — before the window.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Peptides carried back in [`Self::columns`] — peptides, not rows. Zero when the table was
    /// written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// Rows in [`Self::columns`]: returned peptides times samples. Zero when the table was written
    /// to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The offset applied, in peptides.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether peptides were left behind**, by either the limit or the offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// Always `None`: rows that did not become records are not counted for this table.
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// `"minutes"`: the per-sample header is `RetentionTime (min)_`. The column is filled only by
    /// IsoTracker output; [`Self::absent_fields`] names it otherwise.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub retention_time_unit: String,
    /// What this table cannot be trusted to mean, each citing the mzLib source.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Columns this file has no column for, so every value is `None` — `peak_order` and
    /// `retention_time` for anything but IsoTracker output.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"field: ExceptionType"` for each column whose read threw on some rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Always empty for this table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// The table. Empty when it went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
}

/// PTM site occupancy from a MetaMorpheus protein-group table: one row per modified site — what
/// [`read_occupancy`] returns.
///
/// Each row is one site in one occupancy cell: one group, one sample group, one `basis`. Columns:
/// `protein_group_name`, `sample_label`, `basis` (`"count"` or `"intensity"`), `entity_index`,
/// `position`, `is_n_terminus`, `modification`, `fraction`, `numerator`, `denominator` and
/// `cell_is_truncated`.
///
/// **Which number to trust depends on the basis.** A `"count"` cell prints its fraction to two
/// decimals, so use `numerator / denominator` (PSMs modified over PSMs covering the site); an
/// `"intensity"` cell prints its numerator and denominator to four significant digits, so use
/// `fraction`. `numerator` and `denominator` are read with [`Table::floats`] for both bases.
#[derive(Debug, Clone, Deserialize)]
pub struct OccupancyRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `"MetaMorpheusQuantifiedProteinGroups"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// `"ProteinGroupFromTsvFile"`.
    #[serde(default)]
    pub reader: Option<String>,
    /// The file's sample-group labels, in header order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_labels: Vec<String>,
    /// Occupancy cells the writer cut short or replaced with `"Output too long for Excel"`. Their
    /// complete sites are rows with `cell_is_truncated` set; a cut cell with no complete site gives
    /// no rows and is counted only here.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated_cell_count: u64,
    /// Protein groups in the whole file — groups, not rows — before the window.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Groups carried back — groups, not rows. Zero when the table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// Site rows in [`Self::columns`]: one per group, sample group, basis and modified site. Zero
    /// when the table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The offset applied, in groups.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// Whether groups were left behind by the limit or the offset — not the same thing as a cell
    /// the writer cut short; see [`Self::truncated_cell_count`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// Always `None`: rows that did not become records are not counted for this table.
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// Always `None`: this table has no time column.
    #[serde(default)]
    pub retention_time_unit: Option<String>,
    /// What this table cannot be trusted to mean.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Every site column, when the file has no occupancy columns at all (then there are no rows).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"count_occupancy: FormatException"` (or `intensity_occupancy`) when a cell was not an
    /// occupancy cell at all: mzLib refuses to guess, and that cell gives no rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Always empty for this table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// The table. Empty when it went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
}

// ---------------------------------------------------------------------------------------------
// Many files' result
// ---------------------------------------------------------------------------------------------

/// The per-file facts of one input to a `_many` read: one entry of [`ReadBatch::files`].
///
/// One per input, **in input order**; a file's rows in the batch table are those whose
/// `source_index` is its position. The fields are the ones a single-file read reports at its top
/// level, under the same names (BULK.md §3), so what you learn about one reads the same for many.
/// Those that belong to particular readers are `None` elsewhere. For an input that failed under
/// [`OnError::Skip`], [`Self::error`] says why and the facts that could not be established are
/// `None` or empty.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct FileReport {
    /// The absolute path of the input.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// mzLib's `SupportedFileType` name. `None`: the type could not be determined (a failed
    /// input).
    #[serde(default)]
    pub file_type: Option<String>,
    /// The mzLib class that parsed the file. `None`: a failed input.
    #[serde(default)]
    pub reader: Option<String>,
    /// Records mzLib read from the whole file, in the reader's own unit — records, scans passing
    /// the MS-level filter, features, matches, groups or peptides. `None`: a failed input — not
    /// zero, because nothing was counted.
    #[serde(default)]
    pub record_count: Option<u64>,
    /// Data rows that did not become records, counted where one line is one record (the psmtsv
    /// family and MSFragger). `None`: not established for this format, or a failed input.
    #[serde(default)]
    pub rows_not_read: Option<i64>,
    /// The unit of the reader's time columns for this file: `"minutes"`, `"seconds"` or
    /// `"unknown"`. `None`: the reader has no time column, or it is [`read_records_many`], whose
    /// columns are the format's own. [`ReadBatch::in_minutes`] applies it.
    #[serde(default)]
    pub retention_time_unit: Option<String>,
    /// What the reader cannot be trusted to mean for this file, each citing the mzLib source.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// This file's columns, without the `source_index`/`source_path` pair the batch adds.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub column_names: Vec<String>,
    /// Columns the reader defines but **this file has no source for**, so every value in them is
    /// `None`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// `"field: ExceptionType"` for each column whose read threw on some rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
    /// Fields that have no column shape and so do not cross, each with the verb that carries it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// Why this input could not be read, under [`OnError::Skip`]; `None` when it was.
    #[serde(default)]
    pub error: Option<ReadError>,
    /// [`read_records_many`] only: the mzLib record class the columns belong to.
    #[serde(default)]
    pub record_type: Option<String>,
    /// [`read_records_many`] only: the uniform views this file also offers.
    #[serde(default)]
    pub views: Option<Vec<String>>,
    /// [`read_spectra_many`] only: scans in the whole file, before the MS-level filter.
    #[serde(default)]
    pub scan_count: Option<u64>,
    /// [`read_spectra_many`] only: the instrument and acquisition facts.
    #[serde(default)]
    pub source: Option<SpectraSource>,
    /// mzIdentML items mzLib did not represent ([`read_matches_many`], [`read_records_many`]).
    /// `None`: the format keeps no such list.
    #[serde(default)]
    pub skipped_count: Option<u64>,
    /// The same items, each with its reason.
    #[serde(default)]
    pub skipped: Option<Vec<SkippedMatch>>,
    /// The quantification readers only: the file's sample labels, in header order.
    #[serde(default)]
    pub sample_labels: Option<Vec<String>>,
    /// [`read_occupancy_many`] only: occupancy cells the writer cut short.
    #[serde(default)]
    pub truncated_cell_count: Option<u64>,
}

impl FileReport {
    /// Whether this input was read ([`Self::error`] is `None`).
    #[must_use]
    pub fn ok(&self) -> bool {
        self.error.is_none()
    }
}

/// Many files read into ONE long table — what every table-reading `_many` function returns.
///
/// The table's first two columns say where each row came from: `source_index`, the file's 0-based
/// position in the list you passed, and `source_path`, its absolute path. The rest are the
/// single-file function's columns. Rows are grouped by file in input order, each file's rows in
/// the file's own order — **whatever [`BulkOptions::threads`] was**, so the same list always gives
/// the same table. The per-file facts are in [`Self::files`].
#[derive(Debug, Clone, Deserialize)]
pub struct ReadBatch {
    /// The wire verb that produced it, e.g. `"readers read-spectra"`. Not on the wire: set by the
    /// function that made the call, so [`Self::in_minutes`] can name it.
    #[serde(skip)]
    pub verb: String,
    /// Paths given, in files.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Files read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Files that could not be read — always 0 unless [`OnError::Skip`], since otherwise the first
    /// failure is an error. Each one's [`FileReport::error`] says why.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// [`OnError::Fail`] or [`OnError::Skip`], as requested. The thread count is deliberately not
    /// echoed: the result is identical at every thread count.
    #[serde(default)]
    pub on_error: OnError,
    /// Records read, summed over the files read, in the reader's own unit — records, scans,
    /// features, matches, groups or peptides. Not the row count: a long reader gives several rows
    /// per record.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Records carried back in [`Self::columns`], summed over the files, in the reader's own unit
    /// — records, scans, features, matches, groups or peptides. Zero when the table was written to
    /// disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// Rows in [`Self::columns`] — more than [`Self::returned_count`] for a long reader (protein
    /// groups, peptides, occupancy, matches with scores). Zero when the table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// [`read_spectra_many`]: the MS level every file was filtered to, or `None` when no MS-level
    /// filter was applied.
    #[serde(default)]
    pub ms_order: Option<i64>,
    /// [`read_spectra_many`]: whether `mz` and `intensity` are present.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub peaks_included: bool,
    /// [`read_matches_many`]: whether the table is long by score.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub scores_included: bool,
    /// One long table: `source_index`, `source_path`, then the single-file columns. Empty when the
    /// table was written to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline. Written one file at a time,
    /// so memory held at most [`BulkOptions::threads`] files.
    #[serde(default)]
    pub output: Option<WrittenTable>,
    /// One [`FileReport`] per input, in input order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<FileReport>,
}

impl ReadBatch {
    /// The inputs that could not be read, under [`OnError::Skip`]. Empty otherwise.
    #[must_use]
    pub fn failed_files(&self) -> Vec<&FileReport> {
        self.files.iter().filter(|file| !file.ok()).collect()
    }

    /// A time column in minutes, converting each file's rows by that file's own
    /// [`FileReport::retention_time_unit`], so a batch that mixes formats comes back on one axis.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] when the table went to disk, has no such column, or a file's unit is
    /// `"unknown"` or not reported — **raised rather than guessed**.
    pub fn in_minutes(&self, column: &str) -> Result<Vec<Option<f64>>> {
        let values = self.columns.floats(column)?;
        let sources = self.columns.integers("source_index")?;
        values
            .into_iter()
            .zip(sources)
            .map(|(value, source)| {
                let file = source
                    .and_then(|index| usize::try_from(index).ok())
                    .and_then(|index| self.files.get(index));
                let unit = file.and_then(|f| f.retention_time_unit.as_deref());
                match unit {
                    Some("minutes") => Ok(value),
                    Some("seconds") => Ok(value.map(|seconds| seconds / 60.0)),
                    _ => Err(MzLibError::Usage(format!(
                        "Cannot convert {column} from {} for '{}': mzLib gives no basis to say \
                         what unit it is in. TopFD changed from seconds to minutes at v1.7.0 \
                         without changing the file type, so check the values against your \
                         gradient length before comparing them.",
                        self.verb,
                        file.map_or("", |f| f.path.as_str())
                    ))),
                }
            })
            .collect()
    }
}

/// What [`identify_many`] returns: one [`FileInfo`] per path, in the order given.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct IdentifyBatch {
    /// Paths given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Paths identified.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Paths that could not be identified — always 0 unless [`OnError::Skip`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// [`OnError::Fail`] or [`OnError::Skip`], as requested.
    #[serde(default)]
    pub on_error: OnError,
    /// One [`FileInfo`] per path, in input order; a failed one has [`FileInfo::error`] set.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<FileInfo>,
}

impl IdentifyBatch {
    /// The paths that could not be identified, under [`OnError::Skip`].
    #[must_use]
    pub fn failed_files(&self) -> Vec<&FileInfo> {
        self.files
            .iter()
            .filter(|info| info.error.is_some())
            .collect()
    }
}

impl ResultRecords {
    /// `retention_time` converted to minutes, whatever unit the format wrote.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] when the unit is `"unknown"` — **raised rather than guessed**, because
    /// a silently unconverted time axis is the specific mistake this module exists to prevent.
    pub fn retention_time_in_minutes(&self) -> Result<Vec<Option<f64>>> {
        convert_minutes(
            &self.columns,
            "retention_time",
            &self.retention_time_unit,
            &self.file_type,
        )
    }
}

impl FeatureRecords {
    /// `retention_time_start` converted to minutes.
    ///
    /// # Errors
    ///
    /// As [`ResultRecords::retention_time_in_minutes`], and for `_ms1.feature` it **will** raise:
    /// the unit is genuinely unknown there.
    pub fn retention_time_start_in_minutes(&self) -> Result<Vec<Option<f64>>> {
        convert_minutes(
            &self.columns,
            "retention_time_start",
            &self.retention_time_unit,
            &self.file_type,
        )
    }

    /// `retention_time_end` converted to minutes.
    ///
    /// # Errors
    ///
    /// As [`FeatureRecords::retention_time_start_in_minutes`].
    pub fn retention_time_end_in_minutes(&self) -> Result<Vec<Option<f64>>> {
        convert_minutes(
            &self.columns,
            "retention_time_end",
            &self.retention_time_unit,
            &self.file_type,
        )
    }
}

fn convert_minutes(
    table: &Table,
    column: &str,
    unit: &str,
    file_type: &str,
) -> Result<Vec<Option<f64>>> {
    let values = table.floats(column)?;
    match unit {
        "minutes" => Ok(values),
        "seconds" => Ok(values
            .into_iter()
            .map(|value| value.map(|seconds| seconds / 60.0))
            .collect()),
        _ => Err(MzLibError::Usage(format!(
            "Cannot convert retention time for '{file_type}': mzLib gives no basis to say what \
             unit it is in. TopFD changed from seconds to minutes at v1.7.0 without changing the \
             file type, so check the values against your gradient length before comparing them."
        ))),
    }
}

// ---------------------------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------------------------

/// How much of one file to read, and where to put it.
///
/// `limit` and `offset` count in the reading function's own unit, which is what its
/// `returned_count` reports: **records** for [`read_records`] and [`read_results`], **scans** for
/// [`read_spectra`], **features** for [`read_features`], **matches** for [`read_matches`],
/// **groups** for [`read_protein_groups`] and [`read_occupancy`], and **peptides** for
/// [`read_quantified_peptides`]. A long reader gives several rows per group or peptide; the window
/// never splits one.
#[derive(Debug, Clone, Default)]
pub struct ReadOptions {
    /// Return at most this many records (scans for [`read_spectra`], features for
    /// [`read_features`], matches for [`read_matches`], groups for [`read_protein_groups`] and
    /// [`read_occupancy`], peptides for [`read_quantified_peptides`]). `None` returns all of them.
    ///
    /// **There is no default limit**, deliberately: a result file can carry a million rows, and a
    /// library whose default answer is "here's some of it" eventually puts a truncated table in a
    /// paper. `truncated` reports whether anything was left behind.
    pub limit: Option<u64>,
    /// Skip this many records (scans for [`read_spectra`], counted after the MS-level filter;
    /// features for [`read_features`]; matches for [`read_matches`]; groups for
    /// [`read_protein_groups`] and [`read_occupancy`]; peptides for [`read_quantified_peptides`]).
    ///
    /// **A window, not a cursor.** mzLib materialises the whole file on every call — its readers
    /// look lazy and are not — so paging re-reads and re-parses the file once per page, which
    /// makes a loop over pages quadratic. For a large file use [`Self::out`] in one call.
    pub offset: u64,
    /// Write the selected window here as a **tab-separated** table (header = the column names)
    /// and return only a summary. Must differ from the input; parent directories are created. The
    /// intended path for large files, not an escape hatch.
    pub out: Option<String>,
    /// Time to allow. `None` waits indefinitely, which a large file legitimately needs.
    pub timeout: Option<Duration>,
}

/// [`ReadOptions`], plus the two choices only a spectra read has.
#[derive(Debug, Clone, Default)]
pub struct SpectraOptions {
    /// The window and destination, as for every other read. Its `limit` and `offset` count scans.
    pub read: ReadOptions,
    /// Keep only scans at this MS level — `1` for survey scans, `2` for fragment scans. `None`
    /// keeps every level.
    ///
    /// Applied **before** the offset and limit, so `ms_order: Some(2), limit: Some(10)` means the
    /// first ten MS2 scans rather than the MS2 scans among the first ten.
    pub ms_order: Option<u32>,
    /// Include the `mz` and `intensity` arrays.
    ///
    /// **Off by default, and worth leaving so unless you need them.** A scan header is tens of
    /// bytes; its peak list is thousands, and a mid-size mzML holds tens of thousands of scans.
    /// `peak_count` still reports how many peaks each scan has. Written to [`ReadOptions::out`],
    /// each cell is a `;`-joined list.
    pub peaks: bool,
}

/// [`ReadOptions`], plus the one choice only an identification read has.
#[derive(Debug, Clone, Default)]
pub struct MatchOptions {
    /// The window and destination, as for every other read. Its `limit` and `offset` count
    /// matches, even when [`Self::scores`] makes several rows per match.
    pub read: ReadOptions,
    /// Make the table long by score: one row per match and engine score, adding `match_index`
    /// (the match's 0-based position among the file's matches), `score_name` (e.g.
    /// `"MS-GF:SpecEValue"`) and `score_value`. Only mzIdentML records scores (mzLib #1306); other
    /// formats keep one row per match and name the two columns in
    /// [`MatchRecords::absent_fields`].
    pub scores: bool,
}

/// How a `_many` read runs: how many files at once, what one unreadable file does, and where the
/// long table goes.
#[derive(Debug, Clone)]
pub struct BulkOptions {
    /// Files read at once: 1 (the default), more, or `-1` for one per core.
    ///
    /// The result is **identical at any value** — files are always returned in input order — so
    /// this trades memory for speed and never changes an answer. The default is 1 because every
    /// reader holds a whole file in memory, so `threads: 8` can mean eight whole files at once.
    pub threads: i32,
    /// [`OnError::Fail`] (the default) stops at the first file that cannot be read, naming it;
    /// [`OnError::Skip`] records the failure in that file's [`FileReport::error`] and reads the rest.
    pub on_error: OnError,
    /// Write the long table here as **tab-separated** text and return only a summary
    /// ([`ReadBatch::output`]). Files are written one at a time, in order, so memory holds at most
    /// `threads` files however long the list is — the way to read hundreds. A batch that stops on
    /// an error removes its partial table. Must differ from every input.
    pub out: Option<String>,
    /// Time to allow for the whole batch. `None` (the default) waits indefinitely.
    pub timeout: Option<Duration>,
}

impl Default for BulkOptions {
    fn default() -> Self {
        Self {
            threads: 1,
            on_error: OnError::Fail,
            out: None,
            timeout: None,
        }
    }
}

/// [`BulkOptions`], plus the spectra choices, applied to every file.
#[derive(Debug, Clone, Default)]
pub struct SpectraBulkOptions {
    /// Threads, what one unreadable file does, and where the table goes.
    pub bulk: BulkOptions,
    /// Keep only scans at this MS level, in every file. `None` keeps every level.
    pub ms_order: Option<u32>,
    /// Include the `mz` and `intensity` arrays. Off by default, and worth keeping off for a list:
    /// with peaks, write the batch to [`BulkOptions::out`].
    pub peaks: bool,
}

/// [`BulkOptions`], plus the scores choice, applied to every file.
#[derive(Debug, Clone, Default)]
pub struct MatchBulkOptions {
    /// Threads, what one unreadable file does, and where the table goes.
    pub bulk: BulkOptions,
    /// As [`MatchOptions::scores`]: one row per match and score, for the files that record
    /// scores (mzIdentML).
    pub scores: bool,
}

// ---------------------------------------------------------------------------------------------
// Argument assembly
// ---------------------------------------------------------------------------------------------

fn window_args(verb: &str, path: &Path, options: &ReadOptions) -> Result<Vec<String>> {
    let path = path.to_str().ok_or_else(|| {
        MzLibError::Usage("The file path is not valid UTF-8, which the bridge requires.".to_owned())
    })?;
    if path.trim().is_empty() {
        return Err(MzLibError::Usage(
            "A file path is required, e.g. 'AllPSMs.psmtsv'.".to_owned(),
        ));
    }

    let mut args = vec![
        "readers".to_owned(),
        verb.to_owned(),
        "--path".to_owned(),
        path.trim().to_owned(),
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
    Ok(args)
}

fn out_args(out: Option<&str>) -> Result<Vec<String>> {
    match out {
        None => Ok(Vec::new()),
        Some(out) if out.trim().is_empty() => Err(MzLibError::Usage(
            "out must be a non-empty path, or None to return the records.".to_owned(),
        )),
        Some(out) => Ok(vec!["--out".to_owned(), out.trim().to_owned()]),
    }
}

/// The arguments and stdin of a `_many` read: the list goes on stdin, one path per line.
///
/// The one implementation behind every `_many` function. **No loop over files happens here, and
/// none may be added.** The whole list is handed to ONE bridge process, which reads `threads`
/// files at once and returns one table. A Rust-side thread pool would pay .NET start-up once per
/// file and leave nobody owning the total thread count — see the bridge's `design/PARALLELISM.md`.
fn batch_args<P: AsRef<Path>>(
    verb: &str,
    paths: &[P],
    options: &BulkOptions,
) -> Result<(Vec<String>, String)> {
    let lines = bridge::path_lines(paths, "path")?;
    let mut args = vec![
        "readers".to_owned(),
        verb.to_owned(),
        "--paths-stdin".to_owned(),
        "--threads".to_owned(),
        bridge::threads_arg(options.threads)?,
        "--on-error".to_owned(),
        options.on_error.as_str().to_owned(),
    ];
    args.extend(out_args(options.out.as_deref())?);
    Ok((args, lines.join("\n") + "\n"))
}

fn ms_order_args(ms_order: Option<u32>, peaks: bool) -> Result<Vec<String>> {
    let mut args = Vec::new();
    if let Some(ms_order) = ms_order {
        // Rejected here rather than at the bridge so the caller gets the error without a process
        // launch, and gets it in Rust's vocabulary.
        if ms_order == 0 {
            return Err(MzLibError::Usage(
                "ms_order must be 1 or greater; pass None to keep every scan.".to_owned(),
            ));
        }
        args.push("--ms-order".to_owned());
        args.push(ms_order.to_string());
    }
    if peaks {
        args.push("--peaks".to_owned());
    }
    Ok(args)
}

fn read<T: serde::de::DeserializeOwned>(args: &[String], timeout: Option<Duration>) -> Result<T> {
    let data = bridge::invoke(args, None, timeout)?;
    serde_json::from_value(data).map_err(protocol)
}

fn read_many<P: AsRef<Path>>(
    verb: &str,
    paths: &[P],
    options: &BulkOptions,
    extra: Vec<String>,
) -> Result<ReadBatch> {
    let (mut args, stdin) = batch_args(verb, paths, options)?;
    args.extend(extra);
    let data = bridge::invoke(&args, Some(&stdin), options.timeout)?;
    let mut batch: ReadBatch = serde_json::from_value(data).map_err(protocol)?;
    batch.verb = format!("readers {verb}");
    Ok(batch)
}

// ---------------------------------------------------------------------------------------------
// What a file is
// ---------------------------------------------------------------------------------------------

/// Every file type mzLib can recognise.
///
/// Enumerated from mzLib itself rather than from a list maintained here, so it reflects the mzLib
/// the bridge carries and cannot go stale.
#[doc = include_str!("../docs/reference/readers.formats.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let formats = mzlib::readers::formats()?;
/// let quantifiable: Vec<&str> = formats
///     .iter()
///     .filter(|f| f.is_quantifiable())
///     .map(|f| f.file_type.as_str())
///     .collect();
/// assert_eq!(quantifiable, ["psmtsv", "osmtsv", "MsFraggerPsm", "DiaNnReport"]);
/// assert_eq!(formats.iter().filter(|f| f.views.is_empty()).count(), 17);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.formats.see-also.md")]
pub fn formats() -> Result<Vec<Format>> {
    #[derive(Deserialize)]
    struct Payload {
        #[serde(default)]
        formats: Vec<Format>,
    }

    let payload: Payload = read(
        &["readers".to_owned(), "formats".to_owned()],
        Some(Duration::from_secs(60)),
    )?;
    Ok(payload.formats)
}

/// Identify a file — its mzLib type, the extension it dispatched on, its reader and the views it
/// offers — without parsing its contents.
///
/// Cheap by design: mzLib resolves the type and stops, so identifying a million-row file costs no
/// more than identifying an empty one. mzLib has no "unknown" result, so a file is dispatchable
/// or it is an error. For a list, [`identify_many`] answers in one call.
#[doc = include_str!("../docs/reference/readers.identify.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let info = mzlib::readers::identify("PXD078927_msgf_1_1_0.mzid")?;
/// assert_eq!(info.file_type, "MzIdentML");
/// assert_eq!(info.reader.as_deref(), Some("MzIdentMLResultFile"));
/// assert!(info.has_view(mzlib::readers::SPECTRAL_MATCH));
/// assert!(!info.is_quantifiable());
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.identify.see-also.md")]
pub fn identify(path: impl AsRef<Path>) -> Result<FileInfo> {
    let args = window_args("identify", path.as_ref(), &ReadOptions::default())?;
    read(&args, Some(Duration::from_secs(60)))
}

/// Identify many files in one call, without parsing their contents: [`identify`] for a list.
///
/// One bridge process, the answers in the order given — the cheap way to sort a directory of
/// unknown files before reading any of them. There is no table and no `out`; setting
/// [`BulkOptions::out`] is refused before anything is spawned.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the list is empty or holds a blank path, a path is repeated, `threads`
/// is 0 or below -1, `out` is set; or, under [`OnError::Fail`], a file is missing or unrecognised —
/// the message then starts `Input <i> (<path>)`.
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{identify_many, BulkOptions};
/// use mzlib::OnError;
///
/// let batch = identify_many(
///     &[
///         "PXD078927_msgf_1_1_0.mzid",
///         "no-such-run.mzML",
///         "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv",
///     ],
///     &BulkOptions { on_error: OnError::Skip, ..Default::default() },
/// )?;
/// let types: Vec<&str> = batch.files.iter().map(|f| f.file_type.as_str()).collect();
/// assert_eq!(types, ["MzIdentML", "", "MetaMorpheusQuantifiedProteinGroups"]);
/// assert_eq!(batch.failed_files()[0].error.as_ref().unwrap().kind, "usage");
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn identify_many<P: AsRef<Path>>(paths: &[P], options: &BulkOptions) -> Result<IdentifyBatch> {
    if options.out.is_some() {
        return Err(MzLibError::Usage(
            "identify_many returns no table, so it takes no `out`.".to_owned(),
        ));
    }
    let (args, stdin) = batch_args("identify", paths, options)?;
    let data = bridge::invoke(&args, Some(&stdin), options.timeout)?;
    serde_json::from_value(data).map_err(protocol)
}

// ---------------------------------------------------------------------------------------------
// The typed views
// ---------------------------------------------------------------------------------------------

/// Read a result file through the uniform `quantifiable` view, with every default.
///
/// See [`read_results_with`] for the reference: parameters, fields, errors and caveats.
///
/// # Errors
///
/// As [`read_results_with`].
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_results_with, ReadOptions};
///
/// let psms = read_results_with(
///     "FraggerPsm_FragPipev21.1_psm.tsv",
///     &ReadOptions { limit: Some(2), ..Default::default() },
/// )?;
/// assert_eq!(psms.retention_time_unit, "minutes");
/// // MSFragger writes no target/decoy label, so is_decoy is unknown - not false.
/// assert_eq!(psms.columns.booleans("is_decoy")?, [None, None]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_results(path: impl AsRef<Path>) -> Result<ResultRecords> {
    read_results_with(path, &ReadOptions::default())
}

/// Read a result file through the uniform `quantifiable` view: the same columns — sequence,
/// retention time, charge, theoretical mass, decoy flag, protein groups — for every format that
/// offers it.
///
/// Four file types offer the view: MetaMorpheus `.psmtsv`/`.osmtsv`, MSFragger `psm.tsv` and
/// DIA-NN `report.tsv`. Use [`read_records`] for any other format, and [`read_results_many`] for a
/// list.
#[doc = include_str!("../docs/reference/readers.read-results.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_results_with, ReadOptions};
///
/// let psms = read_results_with(
///     "FraggerPsm_FragPipev21.1_psm.tsv",
///     &ReadOptions { limit: Some(2), ..Default::default() },
/// )?;
/// assert_eq!((psms.record_count, psms.returned_count, psms.truncated), (5, 2, true));
/// assert_eq!(psms.rows_not_read, Some(0));        // every data line became a record
/// let minutes = psms.retention_time_in_minutes()?;
/// assert_eq!(minutes[0], Some(0.03233));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-results.see-also.md")]
pub fn read_results_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<ResultRecords> {
    let args = window_args("read-results", path.as_ref(), options)?;
    read(&args, options.timeout)
}

/// Read many result files into one long table of the uniform `quantifiable` view.
///
/// [`read_results_with`] for a list, read by one bridge process. Each file's
/// `retention_time_unit` is in its [`FileReport`], and [`ReadBatch::in_minutes`] converts every
/// file's rows by its own.
#[doc = include_str!("../docs/reference/readers.read-results.bulk.md")]
///
/// # Examples
///
/// Not run: no recording of a many-file read of this view exists yet.
///
/// ```no_run
/// use mzlib::readers::{read_results_many, BulkOptions};
///
/// let batch = read_results_many(&["search1/AllPSMs.psmtsv", "search2/psm.tsv"], &BulkOptions::default())?;
/// let minutes = batch.in_minutes("retention_time")?;   // one axis across both formats
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_results_many<P: AsRef<Path>>(paths: &[P], options: &BulkOptions) -> Result<ReadBatch> {
    read_many("read-results", paths, options, Vec::new())
}

/// Read **any** file mzLib recognises into that format's own fields, with every default.
///
/// See [`read_records_with`] for the reference.
///
/// # Errors
///
/// As [`read_records_with`].
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let table = mzlib::readers::read_records("ToppicPrsm_TopPICv1.6.2_prsm.tsv")?;
/// assert_eq!((table.record_count, table.truncated), (4, false));
/// // TopPIC's alternative identifications are a list of objects: named, not dropped.
/// assert_eq!(table.excluded_fields[0].field, "alternative_identifications");
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_records(path: impl AsRef<Path>) -> Result<NativeRecords> {
    read_records_with(path, &ReadOptions::default())
}

/// Read **any** file mzLib recognises into a table of that format's own record fields, naming
/// every field that could not become a column.
///
/// The exhaustive verb: if [`identify`] succeeds on a path, this reads it — including the 17 file
/// types that belong to no cross-format view at all (TopPIC, Crux, MSFragger's peptide and
/// protein tables, the FlashDeconv formats, the MetaMorpheus and FlashLFQ quantification tables),
/// which no other function in this module can touch. The columns are **not uniform**; see
/// [`NativeRecords`].
///
/// **For SDRF, use [`crate::sdrf::read`] instead.** This verb joins each SDRF row's cells into one
/// semicolon-separated string, and SDRF's `NT=…;AC=…` grammar puts semicolons inside cells, so the
/// joined string cannot be split back apart.
#[doc = include_str!("../docs/reference/readers.read-records.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_records_with, ReadOptions};
///
/// let crux = read_records_with("crux.txt", &ReadOptions { limit: Some(3), ..Default::default() })?;
/// assert_eq!(crux.record_type, "CruxResult");
/// assert_eq!(crux.columns.floats("x_corr_score")?[0], Some(6.4364114));
///
/// // A current FlashLFQ peaks table has no MBR Score column (mzLib #1345): absent, not zero.
/// let peaks = read_records_with(
///     "FlashLFQ_MzLib1.0.591_QuantifiedPeaks.tsv",
///     &ReadOptions { limit: Some(2), ..Default::default() },
/// )?;
/// assert_eq!(peaks.absent_fields, ["mbr_score"]);
///
/// // A MetaMorpheus protein-group table keeps its per-sample values in a dictionary, which
/// // cannot be a column: named, with the verb that does carry it.
/// let groups = read_records_with(
///     "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv",
///     &ReadOptions { limit: Some(2), ..Default::default() },
/// )?;
/// let sample_groups = &groups.excluded_fields[0];
/// assert_eq!(sample_groups.field, "sample_groups");
/// assert_eq!(sample_groups.verb.as_deref(), Some("readers read-protein-groups"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-records.see-also.md")]
pub fn read_records_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<NativeRecords> {
    let args = window_args("read-records", path.as_ref(), options)?;
    read(&args, options.timeout)
}

/// Read many files of **one** record type into one long table of their own fields.
///
/// [`read_records_with`] for a list — a directory of mzIdentML submissions, say. The long table has
/// one column set, so every input must have the same record type: a mixed list is refused before
/// any file is parsed, naming the groups, whatever [`BulkOptions::on_error`] says.
#[doc = include_str!("../docs/reference/readers.read-records.bulk.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_records_many, BulkOptions};
///
/// // The same search, plain and gzipped: mzLib reads .mzid.gz without unpacking it (#1313).
/// let batch = read_records_many(
///     &["PXD078927_msgf_1_1_0.mzid", "PXD078927_msgf_1_1_0.mzid.gz"],
///     &BulkOptions { threads: 2, ..Default::default() },
/// )?;
/// assert_eq!((batch.record_count, batch.row_count), (24, 24));
/// let types: Vec<_> = batch.files.iter().map(|f| f.file_type.as_deref()).collect();
/// assert_eq!(types, [Some("MzIdentML"), Some("MzIdentMLGz")]);
/// assert_eq!(batch.files[1].skipped_count, Some(0));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_records_many<P: AsRef<Path>>(paths: &[P], options: &BulkOptions) -> Result<ReadBatch> {
    read_many("read-records", paths, options, Vec::new())
}

/// Read deconvolved MS1 features through the uniform `ms1_features` view, with every default.
///
/// See [`read_features_with`] for the reference.
///
/// # Errors
///
/// As [`read_features_with`].
pub fn read_features(path: impl AsRef<Path>) -> Result<FeatureRecords> {
    read_features_with(path, &ReadOptions::default())
}

/// Read deconvolved MS1 features through the uniform `ms1_features` view: m/z, charge,
/// retention-time range, apex intensity and isotope count.
///
/// Two file types offer it: TopFD/FLASHDeconv `_ms1.feature` and Dinosaur `.feature.tsv`.
///
/// **One row is not one line of the file for `_ms1.feature`**: mzLib expands each deconvolved
/// feature into one single-charge feature per charge in its recorded range. Dinosaur is
/// one-for-one.
#[doc = include_str!("../docs/reference/readers.read-features.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_features_with, ReadOptions};
///
/// let features = read_features_with(
///     "Ms1Feature_TopFDv1.6.2_ms1.feature",
///     &ReadOptions { limit: Some(5), ..Default::default() },
/// )?;
/// assert_eq!(features.record_count, 25);
/// assert_eq!(features.absent_fields, ["number_of_isotopes"]);   // mzLib never sets it here
/// // TopFD changed its time unit at v1.7.0 without changing the format, so this refuses.
/// assert_eq!(features.retention_time_unit, "unknown");
/// assert!(features.retention_time_start_in_minutes().is_err());
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-features.see-also.md")]
pub fn read_features_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<FeatureRecords> {
    let args = window_args("read-features", path.as_ref(), options)?;
    read(&args, options.timeout)
}

/// Read many feature files into one long table of the `ms1_features` view.
///
/// [`read_features_with`] for a list. Each file's `retention_time_unit` is in its [`FileReport`];
/// for `_ms1.feature` it is `"unknown"`, so [`ReadBatch::in_minutes`] refuses those rows.
#[doc = include_str!("../docs/reference/readers.read-features.bulk.md")]
///
/// # Examples
///
/// Not run: no recording of a many-file read of this view exists yet.
///
/// ```no_run
/// use mzlib::readers::{read_features_many, BulkOptions};
///
/// let batch = read_features_many(&["a.feature.tsv", "b.feature.tsv"], &BulkOptions::default())?;
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_features_many<P: AsRef<Path>>(paths: &[P], options: &BulkOptions) -> Result<ReadBatch> {
    read_many("read-features", paths, options, Vec::new())
}

/// Read identifications through the uniform `spectral_match` view, with every default.
///
/// See [`read_matches_with`] for the reference.
///
/// # Errors
///
/// As [`read_matches_with`].
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let matches = mzlib::readers::read_matches("MsPathFinderT_WithMods_IcTda.tsv")?;
/// assert_eq!(matches.record_count, 5);
/// // No QValue column in this file: q_value is absent, never mzLib's default of 0.
/// assert!(matches.absent_fields.iter().any(|f| f == "q_value"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_matches(path: impl AsRef<Path>) -> Result<MatchRecords> {
    read_matches_with(path, &MatchOptions::default())
}

/// Read identifications through the uniform `spectral_match` view — scan, sequences, accession,
/// decoy flag, modifications and the q-value, rank and threshold a format records — optionally
/// with each engine's scores as long rows.
///
/// Six file types offer it: MsPathFinderT's targets, decoys and combined results, Casanovo's
/// `.mztab`, and mzIdentML `.mzid` / `.mzid.gz`.
#[doc = include_str!("../docs/reference/readers.read-matches.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_matches_with, MatchOptions, ReadOptions};
///
/// let matches = read_matches_with(
///     "PXD078927_msgf_1_1_0.mzid",
///     &MatchOptions { read: ReadOptions { limit: Some(3), ..Default::default() }, ..Default::default() },
/// )?;
/// assert_eq!((matches.record_count, matches.skipped_count), (12, Some(0)));
/// assert_eq!(matches.columns.strings("base_sequence")?[0].as_deref(), Some("HSNLNDATYQRT"));
/// // Every candidate is a row: filter on rank 1 and pass_threshold before counting.
/// assert_eq!(matches.columns.integers("rank")?, [Some(1), Some(2), Some(3)]);
///
/// // The engine's own scores, one row per match and score.
/// let scored = read_matches_with(
///     "PXD078927_msgf_1_1_0.mzid",
///     &MatchOptions { read: ReadOptions { limit: Some(1), ..Default::default() }, scores: true },
/// )?;
/// assert_eq!((scored.returned_count, scored.row_count), (1, 7));
/// assert_eq!(scored.columns.strings("score_name")?[2].as_deref(), Some("MS-GF:SpecEValue"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-matches.see-also.md")]
pub fn read_matches_with(path: impl AsRef<Path>, options: &MatchOptions) -> Result<MatchRecords> {
    let mut args = window_args("read-matches", path.as_ref(), &options.read)?;
    if options.scores {
        args.push("--scores".to_owned());
    }
    read(&args, options.read.timeout)
}

/// Read many identification files into one long table of the `spectral_match` view.
///
/// [`read_matches_with`] for a list — a directory of mzIdentML submissions, say. Each file's
/// skipped items and absent fields are in its [`FileReport`].
#[doc = include_str!("../docs/reference/readers.read-matches.bulk.md")]
///
/// # Examples
///
/// Not run: no recording of a many-file read of this view exists yet.
///
/// ```no_run
/// use mzlib::readers::{read_matches_many, MatchBulkOptions};
///
/// let batch = read_matches_many(&["a.mzid", "b.mzid.gz"], &MatchBulkOptions { scores: true, ..Default::default() })?;
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_matches_many<P: AsRef<Path>>(
    paths: &[P],
    options: &MatchBulkOptions,
) -> Result<ReadBatch> {
    let extra = if options.scores {
        vec!["--scores".to_owned()]
    } else {
        Vec::new()
    };
    read_many("read-matches", paths, &options.bulk, extra)
}

/// Read the scans of a spectra file with every default: every scan's header, no peaks.
///
/// See [`read_spectra_with`] for the reference.
///
/// # Errors
///
/// As [`read_spectra_with`].
pub fn read_spectra(path: impl AsRef<Path>) -> Result<ScanRecords> {
    read_spectra_with(path, &SpectraOptions::default())
}

/// Read the scans of a spectra file: every scan's header always, and its peak arrays only on
/// request.
///
/// Seven file types offer the `spectra` view. **Two of them need Windows**: Bruker `.d` and
/// timsTOF `.d` are read through vendor native libraries and are Windows-x64 only. Thermo `.raw` is
/// managed and reads everywhere. Each read also reports what the file records about the run
/// ([`ScanRecords::source`], mzLib #1349).
#[doc = include_str!("../docs/reference/readers.read-spectra.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_spectra_with, AcquisitionTime, ReadOptions, SpectraOptions};
///
/// // The first two MS2 scans: the MS-level filter applies before the window.
/// let ms2 = read_spectra_with(
///     "sliced_ethcd.mzML",
///     &SpectraOptions {
///         ms_order: Some(2),
///         read: ReadOptions { limit: Some(2), ..Default::default() },
///         ..Default::default()
///     },
/// )?;
/// assert_eq!((ms2.scan_count, ms2.record_count, ms2.returned_count), (6, 5, 2));
/// assert_eq!(ms2.columns.floats("selected_ion_mz")?, [Some(548.453918457031), Some(796.765197753906)]);
///
/// // An mzML with an offset fixes the instant; a Thermo .raw would give a local clock reading.
/// let source = ms2.source.unwrap();
/// assert!(matches!(source.acquired_at(), Some(AcquisitionTime::Utc(_))));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-spectra.see-also.md")]
pub fn read_spectra_with(path: impl AsRef<Path>, options: &SpectraOptions) -> Result<ScanRecords> {
    let mut args = window_args("read-spectra", path.as_ref(), &options.read)?;
    args.extend(ms_order_args(options.ms_order, options.peaks)?);
    read(&args, options.read.timeout)
}

/// Read the scans of many spectra files into one long table.
///
/// [`read_spectra_with`] for a list — every run of an experiment, with each file's instrument,
/// serial number and acquisition start in its [`FileReport::source`]. Retention times are minutes
/// for every format, so the table is on one time axis.
#[doc = include_str!("../docs/reference/readers.read-spectra.bulk.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_spectra_many, BulkOptions, SpectraBulkOptions};
/// use mzlib::OnError;
///
/// let runs = read_spectra_many(
///     &["sliced_ethcd.mzML", "no-such-run.mzML", "withZeros.mgf"],
///     &SpectraBulkOptions {
///         bulk: BulkOptions { threads: 2, on_error: OnError::Skip, ..Default::default() },
///         ..Default::default()
///     },
/// )?;
/// assert_eq!((runs.record_count, runs.row_count), (8, 8));
/// let serials: Vec<_> = runs
///     .files
///     .iter()
///     .map(|f| f.source.as_ref().and_then(|s| s.instrument_serial_number.as_deref()))
///     .collect();
/// assert_eq!(serials, [Some("FSN10189"), None, None]);      // MGF records no instrument
/// assert!(!runs.files[1].ok());
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_spectra_many<P: AsRef<Path>>(
    paths: &[P],
    options: &SpectraBulkOptions,
) -> Result<ReadBatch> {
    let extra = ms_order_args(options.ms_order, options.peaks)?;
    read_many("read-spectra", paths, &options.bulk, extra)
}

// ---------------------------------------------------------------------------------------------
// The quantification tables (mzLib 1.0.592, #1347)
// ---------------------------------------------------------------------------------------------

/// Read a MetaMorpheus protein-group table with every default.
///
/// See [`read_protein_groups_with`] for the reference.
///
/// # Errors
///
/// As [`read_protein_groups_with`].
pub fn read_protein_groups(path: impl AsRef<Path>) -> Result<ProteinGroupRecords> {
    read_protein_groups_with(path, &ReadOptions::default())
}

/// Read a MetaMorpheus `AllQuantifiedProteinGroups.tsv` as one row per protein group per sample
/// group, with each sample group's spectral count and intensity as columns.
///
/// mzLib 1.0.592 reads the per-sample columns into a dictionary (#1347) that [`read_records`]
/// cannot project; this is that dictionary as a **long** table beside the fields you filter on.
/// The group's other fields — coverage, masses, member counts — are in [`read_records`], joined on
/// `protein_group_name`; PTM site occupancy is [`read_occupancy`].
///
/// Needs the bridge from pyMzLib 0.2.0 or later: an older one is refused before anything is read.
#[doc = include_str!("../docs/reference/readers.read-protein-groups.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_protein_groups_with, ReadOptions};
///
/// let groups = read_protein_groups_with(
///     "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv",
///     &ReadOptions { limit: Some(1), ..Default::default() },
/// )?;
/// // One group, eighteen sample groups: eighteen rows.
/// assert_eq!((groups.record_count, groups.returned_count, groups.row_count), (6, 1, 18));
/// assert_eq!(groups.sample_labels[0], "QE-002106_GM1_a-calib");
/// assert_eq!(groups.columns.floats("intensity")?[0], Some(1838317.1674502683));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-protein-groups.see-also.md")]
pub fn read_protein_groups_with(
    path: impl AsRef<Path>,
    options: &ReadOptions,
) -> Result<ProteinGroupRecords> {
    bridge::require_verb("readers read-protein-groups", bridge::MZLIB_1_0_592_BRIDGE)?;
    let args = window_args("read-protein-groups", path.as_ref(), options)?;
    read(&args, options.timeout)
}

/// Read many MetaMorpheus protein-group tables into one long table.
///
/// [`read_protein_groups_with`] for a list — one search per condition, say. Each file's sample
/// labels are in its [`FileReport::sample_labels`].
#[doc = include_str!("../docs/reference/readers.read-protein-groups.bulk.md")]
///
/// # Examples
///
/// Not run: no recording of a many-file read of this table exists yet.
///
/// ```no_run
/// use mzlib::readers::{read_protein_groups_many, BulkOptions};
///
/// let batch = read_protein_groups_many(
///     &["search1/AllQuantifiedProteinGroups.tsv", "search2/AllQuantifiedProteinGroups.tsv"],
///     &BulkOptions::default(),
/// )?;
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_protein_groups_many<P: AsRef<Path>>(
    paths: &[P],
    options: &BulkOptions,
) -> Result<ReadBatch> {
    bridge::require_verb("readers read-protein-groups", bridge::MZLIB_1_0_592_BRIDGE)?;
    read_many("read-protein-groups", paths, options, Vec::new())
}

/// Read a FlashLFQ peptide table with every default.
///
/// See [`read_quantified_peptides_with`] for the reference.
///
/// # Errors
///
/// As [`read_quantified_peptides_with`].
pub fn read_quantified_peptides(path: impl AsRef<Path>) -> Result<QuantifiedPeptideRecords> {
    read_quantified_peptides_with(path, &ReadOptions::default())
}

/// Read a FlashLFQ `QuantifiedPeptides.tsv` (or MetaMorpheus `AllQuantifiedPeptides.tsv`) as one
/// row per peptide per sample, with the sample's intensity and detection type as columns.
///
/// **An intensity of 0 is not a measured zero**: FlashLFQ writes 0 for a peptide it did not
/// quantify in a sample. Filter on `detection_type` before a mean or a log.
///
/// Needs the bridge from pyMzLib 0.2.0 or later: an older one is refused before anything is read.
#[doc = include_str!("../docs/reference/readers.read-quantified-peptides.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_quantified_peptides_with, ReadOptions};
///
/// let peptides = read_quantified_peptides_with(
///     "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv",
///     &ReadOptions { limit: Some(1), ..Default::default() },
/// )?;
/// assert_eq!((peptides.returned_count, peptides.row_count), (1, 18));
/// assert_eq!(peptides.absent_fields, ["peak_order", "retention_time"]);   // not IsoTracker output
/// // A 0 here was not measured: detection_type says so.
/// assert_eq!(peptides.columns.floats("intensity")?[0], Some(0.0));
/// assert_eq!(peptides.columns.strings("detection_type")?[0].as_deref(), Some("NotDetected"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-quantified-peptides.see-also.md")]
pub fn read_quantified_peptides_with(
    path: impl AsRef<Path>,
    options: &ReadOptions,
) -> Result<QuantifiedPeptideRecords> {
    bridge::require_verb(
        "readers read-quantified-peptides",
        bridge::MZLIB_1_0_592_BRIDGE,
    )?;
    let args = window_args("read-quantified-peptides", path.as_ref(), options)?;
    read(&args, options.timeout)
}

/// Read many FlashLFQ peptide tables into one long table.
///
/// [`read_quantified_peptides_with`] for a list.
#[doc = include_str!("../docs/reference/readers.read-quantified-peptides.bulk.md")]
///
/// # Examples
///
/// Not run: no recording of a many-file read of this table exists yet.
///
/// ```no_run
/// use mzlib::readers::{read_quantified_peptides_many, BulkOptions};
///
/// let batch = read_quantified_peptides_many(
///     &["a/QuantifiedPeptides.tsv", "b/QuantifiedPeptides.tsv"],
///     &BulkOptions::default(),
/// )?;
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_quantified_peptides_many<P: AsRef<Path>>(
    paths: &[P],
    options: &BulkOptions,
) -> Result<ReadBatch> {
    bridge::require_verb(
        "readers read-quantified-peptides",
        bridge::MZLIB_1_0_592_BRIDGE,
    )?;
    read_many("read-quantified-peptides", paths, options, Vec::new())
}

/// Read the PTM site occupancy of a MetaMorpheus protein-group table with every default.
///
/// See [`read_occupancy_with`] for the reference.
///
/// # Errors
///
/// As [`read_occupancy_with`].
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let sites = mzlib::readers::read_occupancy("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")?;
/// assert_eq!((sites.record_count, sites.row_count, sites.truncated_cell_count), (6, 95, 0));
/// // A count cell: trust numerator / denominator, which the printed fraction rounds.
/// assert_eq!(sites.columns.strings("basis")?[0].as_deref(), Some("count"));
/// assert_eq!(sites.columns.floats("numerator")?[0], Some(1.0));
/// assert_eq!(sites.columns.floats("denominator")?[0], Some(2.0));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_occupancy(path: impl AsRef<Path>) -> Result<OccupancyRecords> {
    read_occupancy_with(path, &ReadOptions::default())
}

/// Read the PTM site occupancy of a MetaMorpheus `AllQuantifiedProteinGroups.tsv` as one row per
/// group, sample group, basis and modified site.
///
/// MetaMorpheus writes two occupancy cells per group per sample group — one from PSM counts, one
/// from intensities — each a list of modified sites encoded as text. mzLib 1.0.592 parses them
/// (`ModificationOccupancyCell`, #1347); this returns every site of every cell as a row, with
/// `basis` saying which cell it came from. Check [`OccupancyRecords::truncated_cell_count`] and
/// [`OccupancyRecords::failed_fields`]: a cut or malformed cell shortens the table.
///
/// Needs the bridge from pyMzLib 0.2.0 or later: an older one is refused before anything is read.
#[doc = include_str!("../docs/reference/readers.read-occupancy.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_occupancy_with, ReadOptions};
///
/// let sites = read_occupancy_with(
///     "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv",
///     &ReadOptions::default(),
/// )?;
/// assert_eq!(sites.columns.strings("modification")?[0].as_deref(), Some("Deamidation on N"));
/// assert_eq!(sites.columns.integers("position")?[0], Some(329));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-occupancy.see-also.md")]
pub fn read_occupancy_with(
    path: impl AsRef<Path>,
    options: &ReadOptions,
) -> Result<OccupancyRecords> {
    bridge::require_verb("readers read-occupancy", bridge::MZLIB_1_0_592_BRIDGE)?;
    let args = window_args("read-occupancy", path.as_ref(), options)?;
    read(&args, options.timeout)
}

/// Read the PTM site occupancy of many MetaMorpheus protein-group tables into one table.
///
/// [`read_occupancy_with`] for a list. Each file's `truncated_cell_count` is in its
/// [`FileReport`].
#[doc = include_str!("../docs/reference/readers.read-occupancy.bulk.md")]
///
/// # Examples
///
/// Not run: no recording of a many-file read of this table exists yet.
///
/// ```no_run
/// use mzlib::readers::{read_occupancy_many, BulkOptions};
///
/// let sites = read_occupancy_many(
///     &["s1/AllQuantifiedProteinGroups.tsv", "s2/AllQuantifiedProteinGroups.tsv"],
///     &BulkOptions::default(),
/// )?;
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn read_occupancy_many<P: AsRef<Path>>(
    paths: &[P],
    options: &BulkOptions,
) -> Result<ReadBatch> {
    bridge::require_verb("readers read-occupancy", bridge::MZLIB_1_0_592_BRIDGE)?;
    read_many("read-occupancy", paths, options, Vec::new())
}

fn protocol(error: serde_json::Error) -> MzLibError {
    MzLibError::Protocol(format!("readers payload could not be interpreted: {error}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(names: &[&str], columns: &[(&str, Vec<Value>)]) -> Table {
        Table::from_wire(
            names.iter().map(|n| (*n).to_owned()).collect(),
            Some(
                columns
                    .iter()
                    .map(|(name, values)| ((*name).to_owned(), values.clone()))
                    .collect(),
            ),
        )
    }

    #[test]
    fn a_wire_null_becomes_none_rather_than_shortening_a_column() {
        // The hazard that runs through every binding: a shortened column silently misaligns every
        // other column in the table.
        let t = table(
            &["rt"],
            &[("rt", vec![Value::from(1.5), Value::Null, Value::from(3.5)])],
        );

        assert_eq!(t.floats("rt").unwrap(), vec![Some(1.5), None, Some(3.5)]);
        assert_eq!(t.rows(), 3);
    }

    #[test]
    fn a_null_boolean_stays_none_because_none_means_unknown() {
        // Casanovo's is_decoy: de novo sequencing has no target/decoy label, so false would be a
        // fabricated value someone could filter on.
        let t = table(
            &["is_decoy"],
            &[("is_decoy", vec![Value::Null, Value::Null])],
        );
        assert_eq!(t.booleans("is_decoy").unwrap(), vec![None, None]);
    }

    #[test]
    fn an_absent_column_names_the_columns_that_are_there() {
        // The column set is per-format, so a caller who guessed a name from another format has no
        // other way to find out what this one has.
        let t = table(&["a", "b"], &[("a", vec![]), ("b", vec![])]);
        let error = t.floats("c").unwrap_err();

        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(error.to_string().contains("a, b"), "{error}");
    }

    #[test]
    fn a_wrongly_typed_value_is_a_protocol_error_not_a_default() {
        let t = table(&["n"], &[("n", vec![Value::from("not a number")])]);
        assert!(matches!(
            t.floats("n").unwrap_err(),
            MzLibError::Protocol(_)
        ));
    }

    #[test]
    fn column_order_follows_column_names_not_the_map() {
        // A JSON object has no order a parser must preserve, so reading it back from the map would
        // give alphabetical columns rather than mzLib's declaration order.
        let t = table(
            &["z", "a"],
            &[("a", vec![Value::Null]), ("z", vec![Value::Null])],
        );
        assert_eq!(t.names(), ["z", "a"]);
    }

    #[test]
    fn a_name_without_a_column_is_dropped_from_names() {
        // Otherwise names() promises a column that raw() then fails to return.
        let t = table(&["present", "missing"], &[("present", vec![Value::Null])]);
        assert_eq!(t.names(), ["present"]);
        assert!(t.raw("missing").is_none());
    }

    #[test]
    fn per_scan_peak_arrays_stay_arrays() {
        let t = table(
            &["mz"],
            &[(
                "mz",
                vec![Value::from(vec![100.0, 200.5]), Value::from(vec![300.0])],
            )],
        );

        assert_eq!(
            t.float_arrays("mz").unwrap(),
            vec![
                Some(vec![Some(100.0), Some(200.5)]),
                Some(vec![Some(300.0)])
            ]
        );
    }

    fn result_records(file_type: &str, unit: &str, rt: f64) -> ResultRecords {
        serde_json::from_value(serde_json::json!({
            "file_type": file_type,
            "retention_time_unit": unit,
            "column_names": ["retention_time"],
            "columns": {"retention_time": [rt]},
        }))
        .unwrap()
    }

    #[test]
    fn seconds_convert_and_unknown_refuses() {
        let seconds = result_records("MsFraggerPsm", "seconds", 120.0);
        assert_eq!(
            seconds.retention_time_in_minutes().unwrap(),
            vec![Some(2.0)]
        );

        let unknown = result_records("Ms1Feature", "unknown", 2372.27);
        // Raised rather than guessed: mzLib's own deconvolution code guesses here, and this does
        // not.
        let error = unknown.retention_time_in_minutes().unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(error.to_string().contains("no basis to say"), "{error}");
    }

    #[test]
    fn a_result_holds_its_table_and_the_rest_of_the_payload() {
        // The Table is flattened into each result type, so one deserialize reads both the table
        // and the envelope fields around it, from the recording pyMzLib shares.
        let scans: ScanRecords =
            serde_json::from_str(include_str!("../tests/fixtures/readers_spectra_mzml.json"))
                .unwrap();
        assert_eq!(scans.scan_count, 6);
        assert_eq!(scans.columns.rows(), 3);
        assert_eq!(scans.columns.names()[0], "one_based_scan_number");
        assert_eq!(scans.retention_time_unit, "minutes");
        assert!(scans.truncated);
    }

    #[test]
    fn a_view_constant_matches_what_the_bridge_emits() {
        let format = Format {
            file_type: "psmtsv".to_owned(),
            extension: Some(".psmtsv".to_owned()),
            reader: Some("PsmFromTsvFile".to_owned()),
            views: vec![QUANTIFIABLE.to_owned()],
        };
        assert!(format.is_quantifiable());

        let viewless = Format {
            views: vec![],
            ..format
        };
        assert!(!viewless.is_quantifiable());
    }

    #[test]
    fn a_blank_path_is_refused_before_anything_is_spawned() {
        let error = window_args("read-records", Path::new("   "), &ReadOptions::default())
            .expect_err("a blank path must not reach the bridge");
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    #[test]
    fn a_zero_offset_is_not_sent() {
        // A default that is sent explicitly is a default the bridge can later disagree with.
        let args =
            window_args("read-records", Path::new("a.tsv"), &ReadOptions::default()).unwrap();
        assert!(!args.iter().any(|arg| arg == "--offset"));
    }

    #[test]
    fn the_window_is_assembled_in_the_documented_order() {
        let args = window_args(
            "read-records",
            Path::new("a.tsv"),
            &ReadOptions {
                limit: Some(5),
                offset: 2,
                out: Some("out.tsv".to_owned()),
                timeout: None,
            },
        )
        .unwrap();

        assert_eq!(
            args,
            vec![
                "readers",
                "read-records",
                "--path",
                "a.tsv",
                "--limit",
                "5",
                "--offset",
                "2",
                "--out",
                "out.tsv"
            ]
        );
    }

    #[test]
    fn an_empty_out_is_refused() {
        let error = window_args(
            "read-records",
            Path::new("a.tsv"),
            &ReadOptions {
                out: Some("  ".to_owned()),
                ..ReadOptions::default()
            },
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    // ---- the 1.0.592 batch: recordings pyMzLib shares ------------------------------------------

    fn fixture<T: serde::de::DeserializeOwned>(text: &str) -> T {
        serde_json::from_str(text).expect("the recording should deserialize")
    }

    #[test]
    fn a_spectra_read_carries_the_runs_instrument_and_start_time() {
        let scans: ScanRecords =
            fixture(include_str!("../tests/fixtures/readers_spectra_mzml.json"));
        let source = scans.source.expect("an mzML records its source");
        assert_eq!(
            source.instrument_model_accession.as_deref(),
            Some("MS:1002416")
        );
        assert!(source.acquisition_start_time_is_utc);
        match source.acquired_at() {
            Some(AcquisitionTime::Utc(instant)) => {
                assert_eq!(instant.to_rfc3339(), "2021-03-16T17:09:07+00:00")
            }
            other => panic!("an mzML with Z is an instant: {other:?}"),
        }
        assert!(scans.absent_fields.is_empty());
    }

    #[test]
    fn a_local_clock_reading_stays_local() {
        // Thermo .raw records the acquisition PC's wall clock with no offset; a timezone-aware
        // type would have to invent one.
        let source = SpectraSource {
            acquisition_start_time: Some("2021-03-16T12:09:07".to_owned()),
            ..SpectraSource::default()
        };
        assert!(matches!(
            source.acquired_at(),
            Some(AcquisitionTime::Local(_))
        ));
        assert_eq!(SpectraSource::default().acquired_at(), None);
    }

    #[test]
    fn an_mzidentml_read_reports_what_it_skipped_and_what_is_absent() {
        let matches: MatchRecords =
            fixture(include_str!("../tests/fixtures/readers_matches_mzid.json"));
        assert_eq!(matches.skipped_count, Some(0));
        assert_eq!(matches.skipped, Some(vec![]));
        assert_eq!(matches.absent_fields, ["is_decoy"]);
        assert!(!matches.scores_included);
        assert_eq!(matches.row_count, matches.returned_count);

        // Every format but mzIdentML keeps no skip list, and says so with None, not zero.
        let casanovo: MatchRecords = fixture(include_str!(
            "../tests/fixtures/readers_matches_casanovo.json"
        ));
        assert_eq!(casanovo.skipped_count, None);
        assert!(casanovo.absent_fields.iter().any(|f| f == "q_value"));
    }

    #[test]
    fn scores_make_the_match_table_long() {
        let scored: MatchRecords = fixture(include_str!(
            "../tests/fixtures/readers_matches_mzid_scores.json"
        ));
        assert!(scored.scores_included);
        assert_eq!((scored.returned_count, scored.row_count), (1, 7));
        assert_eq!(
            scored.columns.integers("match_index").unwrap(),
            vec![Some(0); 7]
        );
    }

    #[test]
    fn an_excluded_field_names_the_verb_that_carries_it() {
        let groups: NativeRecords = fixture(include_str!(
            "../tests/fixtures/readers_records_mm_protein_groups.json"
        ));
        assert_eq!(groups.excluded_fields[0].field, "sample_groups");
        assert_eq!(
            groups.excluded_fields[0].verb.as_deref(),
            Some("readers read-protein-groups")
        );
        let toppic: NativeRecords = fixture(include_str!(
            "../tests/fixtures/readers_records_toppic.json"
        ));
        assert_eq!(toppic.excluded_fields[0].verb, None);
    }

    #[test]
    fn a_missing_optional_column_is_absent_not_zero() {
        let peaks: NativeRecords = fixture(include_str!(
            "../tests/fixtures/readers_records_flashlfq_peaks.json"
        ));
        assert_eq!(peaks.absent_fields, ["mbr_score"]);
        assert_eq!(
            peaks.columns.floats("mbr_score").unwrap(),
            vec![None; peaks.columns.rows()]
        );
    }

    #[test]
    fn the_quantification_tables_are_long() {
        let groups: ProteinGroupRecords = fixture(include_str!(
            "../tests/fixtures/readers_protein_groups.json"
        ));
        assert_eq!(groups.sample_labels.len(), 18);
        assert_eq!(groups.row_count, 18);
        assert_eq!(groups.columns.rows(), 18);
        assert!(groups
            .excluded_fields
            .iter()
            .all(|f| f.verb.as_deref() == Some("readers read-occupancy")));

        let peptides: QuantifiedPeptideRecords = fixture(include_str!(
            "../tests/fixtures/readers_quantified_peptides.json"
        ));
        assert_eq!(peptides.retention_time_unit, "minutes");
        assert_eq!(peptides.absent_fields, ["peak_order", "retention_time"]);

        let sites: OccupancyRecords =
            fixture(include_str!("../tests/fixtures/readers_occupancy.json"));
        assert_eq!((sites.row_count, sites.truncated_cell_count), (95, 0));
        // A count-basis numerator crosses as a whole number and must still read as a float.
        assert_eq!(sites.columns.floats("numerator").unwrap()[0], Some(1.0));
    }

    #[test]
    fn a_batch_keeps_each_files_facts_and_its_failures_apart() {
        let batch: ReadBatch = fixture(include_str!("../tests/fixtures/readers_many_spectra.json"));
        assert_eq!(
            (batch.file_count, batch.read_count, batch.failed_count),
            (3, 2, 1)
        );
        assert_eq!(batch.on_error, OnError::Skip);
        assert_eq!(batch.columns.names()[..2], ["source_index", "source_path"]);
        assert_eq!(batch.files.len(), 3);
        let failed = batch.failed_files();
        assert_eq!(failed.len(), 1);
        // A failed input's facts are None, not zero: nothing was counted.
        assert_eq!(failed[0].record_count, None);
        assert_eq!(failed[0].error.as_ref().unwrap().r#type, "usage");
        assert_eq!(batch.files[2].file_type.as_deref(), Some("Mgf"));
        assert_eq!(
            batch.files[2].source.as_ref().unwrap().instrument_model,
            None
        );
    }

    #[test]
    fn a_batch_converts_each_files_rows_by_its_own_unit() {
        let mut batch: ReadBatch =
            fixture(include_str!("../tests/fixtures/readers_many_spectra.json"));
        batch.verb = "readers read-spectra".to_owned();
        let minutes = batch.in_minutes("retention_time").unwrap();
        assert_eq!(minutes, batch.columns.floats("retention_time").unwrap());

        // A file whose unit is unknown refuses rather than guessing, naming the file.
        batch.files[0].retention_time_unit = Some("unknown".to_owned());
        let error = batch.in_minutes("retention_time").unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(error.to_string().contains("sliced_ethcd.mzML"), "{error}");
    }

    #[test]
    fn an_identify_batch_answers_every_path_in_order() {
        let batch: IdentifyBatch =
            fixture(include_str!("../tests/fixtures/readers_many_identify.json"));
        assert_eq!(batch.files.len(), 3);
        assert_eq!(batch.files[1].file_type, "");
        assert_eq!(batch.failed_files().len(), 1);
        assert_eq!(
            batch.files[2].extension.as_deref(),
            Some("QuantifiedProteinGroups.tsv")
        );
    }

    #[test]
    fn a_batch_sends_its_list_on_stdin_and_its_options_as_arguments() {
        let (args, stdin) = batch_args(
            "read-spectra",
            &["a.mzML", "b.mzML"],
            &BulkOptions {
                threads: -1,
                on_error: OnError::Skip,
                out: Some("all.tsv".to_owned()),
                timeout: None,
            },
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "readers",
                "read-spectra",
                "--paths-stdin",
                "--threads",
                "-1",
                "--on-error",
                "skip",
                "--out",
                "all.tsv"
            ]
        );
        assert_eq!(stdin, "a.mzML\nb.mzML\n");
        // No window: offset and limit are refused with --paths-stdin, so they are never sent.
        assert!(!args.iter().any(|a| a == "--limit" || a == "--offset"));
    }

    #[test]
    fn the_default_batch_reads_one_file_at_a_time_and_fails_on_the_first_error() {
        let options = BulkOptions::default();
        assert_eq!((options.threads, options.on_error), (1, OnError::Fail));
        let (args, _) = batch_args("read-records", &["a.tsv"], &options).unwrap();
        assert_eq!(args[3..], ["--threads", "1", "--on-error", "fail"]);
    }

    #[test]
    fn a_bad_batch_is_refused_before_anything_is_spawned() {
        let empty: [&str; 0] = [];
        assert!(batch_args("read-records", &empty, &BulkOptions::default()).is_err());
        let zero = BulkOptions {
            threads: 0,
            ..BulkOptions::default()
        };
        assert!(batch_args("read-records", &["a.tsv"], &zero).is_err());
        assert!(batch_args("read-records", &["a\nb.tsv"], &BulkOptions::default()).is_err());
        let with_out = BulkOptions {
            out: Some("x.tsv".to_owned()),
            ..BulkOptions::default()
        };
        let error = identify_many(&["a.tsv"], &with_out).unwrap_err();
        assert!(error.to_string().contains("no table"), "{error}");
    }

    #[test]
    fn scores_and_the_ms_level_reach_the_bridge() {
        assert_eq!(
            ms_order_args(Some(2), true).unwrap(),
            ["--ms-order", "2", "--peaks"]
        );
        assert!(ms_order_args(None, false).unwrap().is_empty());
        assert!(matches!(
            ms_order_args(Some(0), false),
            Err(MzLibError::Usage(_))
        ));
    }
}
