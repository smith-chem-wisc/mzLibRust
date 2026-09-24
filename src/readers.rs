//! Reading mass-spectrometry data files and proteomics search results: what a file *is*, and
//! every one of them read.
//!
//! **Spectra files are read here, not just search output.** [`read_spectra`] reads **mzML**,
//! Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF and msalign — scan headers always, peaks opt-in:
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
//! | [`read_matches`] | 6 | uniform `spectral_match` view: scan, sequences, accession, decoy flag |
//! | [`read_spectra`] | 7 | scan headers; peaks opt-in |
//!
//! The rule of thumb: **a typed view when you need numbers that mean the same thing across files,
//! and [`read_records`] when you need everything one file has.** A `.psmtsv` through
//! [`read_results`] gives 10 comparable columns; the same file through [`read_records`] gives 73,
//! including the q-values and scores the uniform view does not carry.
//!
//! An empty [`FileInfo::views`] is a real and common answer — **17 of the 36** have it, meaning
//! mzLib parses the file into a shape that shares nothing with any other format. Those are exactly
//! the files [`read_records`] exists for.
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
//! ## Units are not normalised across formats
//!
//! mzLib's result-file readers pass through whatever the writing tool wrote. MetaMorpheus and
//! MSFragger retention times are minutes; TopFD's `_ms1.feature` changed from seconds to minutes
//! at v1.7.0 without changing the file type. Every typed view therefore reports a
//! `retention_time_unit`, and the `*_in_minutes` methods convert — or refuse, when mzLib gives no
//! basis to say. Spectra readers are the exception: they convert to minutes at the boundary.
//!
//! **Nothing here is FDR-filtered.** Every result format records confidence somewhere, and none of
//! the uniform views filters on it. Filter before you report.

use std::collections::BTreeMap;
use std::path::Path;
use std::time::Duration;

use serde::Deserialize;
use serde_json::Value;

use crate::bridge::{self, MzLibError, Result};

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

/// What a particular file is, and what can be done with it: what [`identify`] returns.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct FileInfo {
    /// The absolute path that was identified.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// mzLib's `SupportedFileType` name.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The extension mzLib dispatched on. `None` when mzLib maps no extension to this type.
    #[serde(default)]
    pub extension: Option<String>,
    /// The mzLib reader class that would parse it.
    #[serde(default)]
    pub reader: Option<String>,
    /// The uniform views this file supports. **Empty is a real answer**, and means mzLib can read
    /// the file but offers no cross-format projection of it — use [`read_records`].
    #[serde(default)]
    pub views: Vec<String>,
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
    /// Rows written, excluding the header: the window that was selected, not the whole file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
}

/// A field of a record type that could not become a column, and why.
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
}

/// The uniform `quantifiable` view of a result file — what [`read_results`] returns.
#[derive(Debug, Clone, Deserialize)]
pub struct ResultRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// The mzLib `SupportedFileType` that was dispatched.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
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
    /// Records mzLib parsed from the **whole file**, before the window. Records, not lines: see
    /// the caveats on [`read_records_with`].
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
    /// **Fields that could not become columns**, each with the reason. A nested object or a
    /// dictionary has no faithful column shape, and inventing one would mean publishing a
    /// schema mzLib does not have — so they are named rather than dropped, because a column
    /// that simply vanished is indistinguishable from a field the format does not have.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub excluded_fields: Vec<ExcludedField>,
    /// Fields that **raised** while being read, as `"field: ExceptionType"`. Several mzLib
    /// properties are computed and assume a UniProt-style FASTA header — Crux's and
    /// MsPathFinderT's `accession` are both `protein_id` split on `|` — so on other databases
    /// they throw. Those cells arrive as `null` rather than failing the whole read, but a
    /// failure must not look like missing data.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_fields: Vec<String>,
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
    /// What this view cannot be trusted to mean for this file, each citing the mzLib source.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
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
/// **Nothing here is FDR-filtered.** mzLib's `ISpectralMatch` carries identity fields; every
/// format offering this view records an E-value or q-value that [`read_records`] will give you.
#[derive(Debug, Clone, Deserialize)]
pub struct MatchRecords {
    /// The absolute path that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// The mzLib `SupportedFileType` that was dispatched.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// Matches in the **whole file**, before the window.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Matches carried back in [`Self::columns`], the unit the window counts in. Zero when the
    /// table was written to disk.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// The offset applied, in matches.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether matches were left behind**, by either the limit or the offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// What this view cannot be trusted to mean for this file — that MsPathFinderT infers decoys
    /// from an `XXX` name prefix, that Casanovo's scan numbers are mzTab indices, that mzIdentML
    /// lists every candidate rank.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// The table, one row per match. Empty when it went to disk.
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
    /// Format-specific traps for this file — that msalign holds deconvolved neutral masses rather
    /// than m/z, that MGF scan numbers may be file order, that Bruker needs Windows x64.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// The table, one row per scan. With [`SpectraOptions::peaks`], `mz` and `intensity` are each
    /// one array per scan. Empty when it went to disk.
    #[serde(flatten)]
    pub columns: Table,
    /// Where the table was written, or `None` if it came back inline.
    #[serde(default)]
    pub output: Option<WrittenTable>,
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

/// How much of a file to read, and where to put it.
///
/// `limit` and `offset` count in the reading function's own unit, which is what its
/// `returned_count` reports: **records** for [`read_records`] and [`read_results`], **scans** for
/// [`read_spectra`], **features** for [`read_features`], **matches** for [`read_matches`].
#[derive(Debug, Clone, Default)]
pub struct ReadOptions {
    /// Return at most this many records (scans for [`read_spectra`], features for
    /// [`read_features`], matches for [`read_matches`]). `None` returns all of them.
    ///
    /// **There is no default limit**, deliberately: a result file can carry a million rows, and a
    /// library whose default answer is "here's some of it" eventually puts a truncated table in a
    /// paper. `truncated` reports whether anything was left behind.
    pub limit: Option<u64>,
    /// Skip this many records (scans for [`read_spectra`], counted after the MS-level filter;
    /// features for [`read_features`]; matches for [`read_matches`]).
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
    /// `peak_count` still reports how many peaks each scan has.
    pub peaks: bool,
}

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
    if let Some(out) = &options.out {
        if out.trim().is_empty() {
            return Err(MzLibError::Usage(
                "out must be a non-empty path, or None to return the records.".to_owned(),
            ));
        }
        args.push("--out".to_owned());
        args.push(out.trim().to_owned());
    }

    Ok(args)
}

fn read<T: serde::de::DeserializeOwned>(args: &[String], timeout: Option<Duration>) -> Result<T> {
    let data = bridge::invoke(args, None, timeout)?;
    serde_json::from_value(data).map_err(protocol)
}

// ---------------------------------------------------------------------------------------------
// The public surface
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
/// or it is an error.
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
/// DIA-NN `report.tsv`. Use [`read_records`] for any other format.
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
/// let minutes = psms.retention_time_in_minutes()?;
/// assert_eq!(minutes[0], Some(0.03233));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-results.see-also.md")]
pub fn read_results_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<ResultRecords> {
    let args = window_args("read-results", path.as_ref(), options)?;
    read(&args, options.timeout)
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
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-records.see-also.md")]
pub fn read_records_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<NativeRecords> {
    let args = window_args("read-records", path.as_ref(), options)?;
    read(&args, options.timeout)
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

/// Read identifications through the uniform `spectral_match` view, with every default.
///
/// See [`read_matches_with`] for the reference.
///
/// # Errors
///
/// As [`read_matches_with`].
pub fn read_matches(path: impl AsRef<Path>) -> Result<MatchRecords> {
    read_matches_with(path, &ReadOptions::default())
}

/// Read identifications through the uniform `spectral_match` view: scan, sequences, accession,
/// decoy flag and modifications.
///
/// Six file types offer it: MsPathFinderT's targets, decoys and combined results, Casanovo's
/// `.mztab`, and mzIdentML `.mzid` / `.mzid.gz`.
#[doc = include_str!("../docs/reference/readers.read-matches.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_matches_with, ReadOptions};
///
/// let matches = read_matches_with(
///     "PXD078927_msgf_1_1_0.mzid",
///     &ReadOptions { limit: Some(3), ..Default::default() },
/// )?;
/// assert_eq!(matches.record_count, 12);
/// assert_eq!(matches.columns.strings("base_sequence")?[0].as_deref(), Some("HSNLNDATYQRT"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-matches.see-also.md")]
pub fn read_matches_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<MatchRecords> {
    let args = window_args("read-matches", path.as_ref(), options)?;
    read(&args, options.timeout)
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
/// managed and reads everywhere.
#[doc = include_str!("../docs/reference/readers.read-spectra.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::readers::{read_spectra_with, ReadOptions, SpectraOptions};
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
/// assert_eq!(ms2.columns.integers("ms_order")?, [Some(2), Some(2)]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/readers.read-spectra.see-also.md")]
pub fn read_spectra_with(path: impl AsRef<Path>, options: &SpectraOptions) -> Result<ScanRecords> {
    let mut args = window_args("read-spectra", path.as_ref(), &options.read)?;

    if let Some(ms_order) = options.ms_order {
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
    if options.peaks {
        args.push("--peaks".to_owned());
    }

    read(&args, options.read.timeout)
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
}
