//! Reading proteomics result files: what a file *is*, and every one of them read.
//!
//! mzLib recognises **31 file types** written by a dozen search and deconvolution tools —
//! MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT, Crux, Casanovo, FlashDeconv, Dinosaur,
//! DIA-NN, FlashLFQ — and maintains a parser for each. **All 31 are readable here.**
//!
//! ```no_run
//! # fn main() -> Result<(), mzlib::MzLibError> {
//! let info = mzlib::readers::identify("psm.tsv")?;
//! println!("{} {:?}", info.file_type, info.views);   // MsFraggerPsm ["quantifiable"]
//!
//! let table = mzlib::readers::read_records("toppic_prsm.tsv")?;
//! let e_values = table.columns.floats("e_value")?;   // Vec<Option<f64>>
//! # Ok(())
//! # }
//! ```
//!
//! ## Choosing a function
//!
//! What differs between formats is not *whether* you can read them but *what the columns mean*.
//!
//! | function | reads | columns |
//! |---|---|---|
//! | [`read_records`] | **all 31** | **that format's own fields**, under mzLib's names |
//! | [`read_results`] | 4 | uniform `quantifiable` view |
//! | [`read_features`] | 2 | uniform `ms1_features` view |
//! | [`read_matches`] | 4 | uniform `spectral_match` view |
//! | [`read_spectra`] | 7 | scan headers; peaks opt-in |
//!
//! The rule of thumb: **a typed view when you need numbers that mean the same thing across files,
//! and [`read_records`] when you need everything one file has.** A `.psmtsv` through
//! [`read_results`] gives 10 comparable columns; the same file through [`read_records`] gives 73,
//! including the q-values and scores the uniform view does not carry.
//!
//! An empty [`FileInfo::views`] is a real and common answer — **14 of the 31** have it, meaning
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
//! ```no_run
//! # fn main() -> Result<(), mzlib::MzLibError> {
//! let t = mzlib::readers::read_records("crux.txt")?;
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
//! # Ok(())
//! # }
//! ```

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
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Table {
    names: Vec<String>,
    columns: BTreeMap<String, Vec<Value>>,
}

impl Table {
    /// The column names, **in the order mzLib declares them** — base-class fields first, then
    /// each format's own.
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
    /// `None` genuinely means *unknown* in this library rather than *false* — Casanovo's
    /// `is_decoy` is the case that matters, because de novo sequencing has no target/decoy label
    /// at all and a `false` there would be a fabricated value someone could filter on.
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

// ---------------------------------------------------------------------------------------------
// Wire types
// ---------------------------------------------------------------------------------------------

/// One file type mzLib can recognise.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct Format {
    /// mzLib's `SupportedFileType` name, e.g. `"MsFraggerPsm"`, `"psmtsv"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The extension or filename suffix mzLib dispatches on, e.g. `"psm.tsv"`, `"_ms1.feature"`.
    ///
    /// **Not unique across file types.** `BrukerD` and `BrukerTimsTof` are both `.d`, told apart
    /// by which analysis file the directory holds, and several formats share `.tsv`.
    #[serde(default)]
    pub extension: Option<String>,
    /// The name of the mzLib class that parses it, for cross-referencing the mzLib source.
    #[serde(default)]
    pub reader: Option<String>,
    /// The uniform views this format supports. Often empty — 14 of 31 have none.
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

/// What a particular file is, and what can be done with it.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct FileInfo {
    /// The absolute path that was identified.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// mzLib's `SupportedFileType` name.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_type: String,
    /// The extension mzLib dispatched on.
    #[serde(default)]
    pub extension: Option<String>,
    /// The mzLib class that would parse it.
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
    /// comparable; `MsFraggerPsm` is quantifiable by interface and its retention times used to be
    /// in seconds. See [`read_results`] and the caveats it returns.
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
    /// Rows written, excluding the header.
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
    /// Why it could not cross the wire.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub reason: String,
}

/// The fields every read verb reports.
#[derive(Debug, Clone, Default, Deserialize)]
struct Common {
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    path: String,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    file_type: String,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    record_count: u64,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    returned_count: u64,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    offset: u64,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    truncated: bool,
    #[serde(default)]
    column_names: Vec<String>,
    #[serde(default)]
    columns: Option<BTreeMap<String, Vec<Value>>>,
    #[serde(default)]
    output: Option<WrittenTable>,
}

/// Generates a public result type over [`Common`] plus its own fields.
///
/// Written as a macro rather than five hand-rolled structs so the shared fields — and the
/// invariant that `record_count` counts the whole file while `returned_count` counts what came
/// back — cannot drift between the five verbs.
macro_rules! record_type {
    (
        $(#[$meta:meta])*
        $name:ident { $( $(#[$field_meta:meta])* $field:ident : $ty:ty ),* $(,)? }
    ) => {
        $(#[$meta])*
        #[derive(Debug, Clone)]
        pub struct $name {
            /// The absolute path that was read.
            pub path: String,
            /// mzLib's `SupportedFileType` name.
            pub file_type: String,
            /// Records in the **whole file**, regardless of any limit or offset.
            pub record_count: u64,
            /// Records actually carried back in [`Self::columns`]. Zero when the table was
            /// written to disk instead.
            pub returned_count: u64,
            /// The offset that was applied.
            pub offset: u64,
            /// **Whether records were left behind**, by either the limit or the offset. A short
            /// answer and a complete one must never look alike, so check this rather than
            /// comparing counts yourself.
            pub truncated: bool,
            /// The table. Empty when the records went to disk.
            pub columns: Table,
            /// Where the table was written, or `None` if it came back inline.
            pub output: Option<WrittenTable>,
            $( $(#[$field_meta])* pub $field : $ty ),*
        }

        impl $name {
            fn build(common: Common $(, $field: $ty)*) -> Self {
                Self {
                    path: common.path,
                    file_type: common.file_type,
                    record_count: common.record_count,
                    returned_count: common.returned_count,
                    offset: common.offset,
                    truncated: common.truncated,
                    columns: Table::from_wire(common.column_names, common.columns),
                    output: common.output,
                    $( $field ),*
                }
            }
        }
    };
}

record_type! {
    /// The uniform record view of a result file — what [`read_results`] returns.
    ResultRecords {
        /// The unit [`Self::columns`]' `retention_time` carries for this format: `"minutes"`,
        /// `"seconds"`, or `"unknown"`. mzLib does not normalise it, so it differs per format.
        /// Convert with [`ResultRecords::retention_time_in_minutes`] rather than by hand.
        retention_time_unit: String,
        /// Data rows that did not become records. mzLib drops a malformed row **silently**, so a
        /// non-zero value here means the file is partly unreadable and the table is incomplete.
        /// `None` when the count could not be established meaningfully.
        rows_not_read: Option<i64>,
        /// **What the uniform view cannot be trusted to mean for this format**, each citing the
        /// mzLib source it came from. Worth reading before comparing anything across formats.
        caveats: Vec<String>,
    }
}

record_type! {
    /// Every field of one format, whatever format it is — what [`read_records`] returns.
    ///
    /// The columns here are **not uniform**: they are this format's own mzLib record fields, under
    /// mzLib's own names in `snake_case`, which makes them cross-referenceable against the mzLib
    /// source. A column called `e_value` is `ToppicPrsm.EValue`, and [`Self::record_type`] names
    /// the class to look in.
    NativeRecords {
        /// The mzLib class that parsed the file.
        reader: Option<String>,
        /// The mzLib record class the columns came from, e.g. `"ToppicPrsm"`.
        record_type: String,
        /// The uniform views this file *also* supports, if any. Often empty.
        views: Vec<String>,
        /// **Fields that could not become columns**, each with the reason. A nested object or a
        /// dictionary has no faithful column shape, and inventing one would mean publishing a
        /// schema mzLib does not have — so they are named rather than dropped, because a column
        /// that simply vanished is indistinguishable from a field the format does not have.
        excluded_fields: Vec<ExcludedField>,
        /// Fields that **raised** while being read, with the exception type. Several mzLib
        /// properties are computed and assume a UniProt-style FASTA header — Crux's and
        /// MsPathFinderT's `accession` are both `protein_id` split on `|` — so on other databases
        /// they throw. Those cells arrive as `null` rather than failing the whole read, but a
        /// failure must not look like missing data.
        failed_fields: Vec<String>,
    }
}

record_type! {
    /// Deconvolved MS1 features — what [`read_features`] returns.
    FeatureRecords {
        /// `"minutes"` or `"unknown"`.
        ///
        /// **It is genuinely `"unknown"` for `_ms1.feature`.** TopFD wrote seconds through v1.6.2
        /// and minutes from v1.7.0 without changing the file type, and mzLib normalises neither.
        /// That is not a gap in this crate; it is the honest state of the format.
        retention_time_unit: String,
        /// What this view cannot be trusted to mean for this format.
        caveats: Vec<String>,
    }
}

record_type! {
    /// Identifications — what [`read_matches`] returns.
    ///
    /// **Nothing here is FDR-filtered, and there is no confidence column to filter on**: mzLib's
    /// `ISpectralMatch` carries identity fields only. Every format offering this view records an
    /// E-value or q-value that [`read_records`] will give you.
    MatchRecords {
        /// What this view cannot be trusted to mean for this format — that MsPathFinderT infers
        /// decoys from an `XXX` name prefix, and that Casanovo's `is_decoy` is `null` because de
        /// novo sequencing has no target/decoy label at all.
        caveats: Vec<String>,
    }
}

record_type! {
    /// Scan headers, and optionally peaks — what [`read_spectra`] returns.
    ScanRecords {
        /// The mzLib class that parsed it, e.g. `"Mzml"`, `"ThermoRawFileReader"`.
        reader: Option<String>,
        /// Scans in the **whole file**, before any MS-level filter. Reported alongside
        /// `record_count` so a filter that matched nothing can never look like an empty file.
        scan_count: u64,
        /// The MS level filtered to, or `None` if unfiltered.
        ms_order: Option<i64>,
        /// Whether `mz` and `intensity` are present — read them with [`Table::float_arrays`].
        peaks_included: bool,
        /// Always `"minutes"` for this view: mzLib's spectra readers convert at the boundary,
        /// unlike its result-file readers.
        retention_time_unit: String,
        /// What this view cannot be trusted to mean for this format.
        caveats: Vec<String>,
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

/// How much of a file to read, and where to put it.
#[derive(Debug, Clone, Default)]
pub struct ReadOptions {
    /// Maximum records to return. `None` returns all of them.
    ///
    /// **There is no default limit**, deliberately: a result file can carry a million rows, and a
    /// library whose default answer is "here's some of it" eventually puts a truncated table in a
    /// paper. `truncated` reports whether anything was left behind.
    pub limit: Option<u64>,
    /// Records to skip.
    ///
    /// **A window, not a cursor.** mzLib materialises the whole file on every call — its readers
    /// look lazy and are not — so paging re-reads and re-parses the file once per page, which
    /// makes a loop over pages quadratic. For a large file use [`Self::out`] in one call.
    pub offset: u64,
    /// Write the records here as a **tab-separated** table and return only a summary. The intended
    /// path for large files, not an escape hatch.
    pub out: Option<String>,
    /// Seconds to allow. `None` waits indefinitely, which a large file legitimately needs.
    pub timeout: Option<Duration>,
}

/// [`ReadOptions`], plus the two choices only a spectra read has.
#[derive(Debug, Clone, Default)]
pub struct SpectraOptions {
    /// The window and destination, as for every other read.
    pub read: ReadOptions,
    /// Keep only scans at this MS level — `1` for survey scans, `2` for fragment scans.
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

// ---------------------------------------------------------------------------------------------
// The public surface
// ---------------------------------------------------------------------------------------------

/// Every file type mzLib can recognise.
///
/// Enumerated from mzLib itself rather than from a list maintained here, so it reflects the
/// installed version and cannot go stale.
///
/// # Errors
///
/// [`MzLibError::Bridge`] if mzLib itself failed; [`MzLibError::Protocol`] if the payload cannot
/// be interpreted.
pub fn formats() -> Result<Vec<Format>> {
    #[derive(Deserialize)]
    struct Payload {
        #[serde(default)]
        formats: Vec<Format>,
    }

    let data = bridge::invoke(
        &["readers".to_owned(), "formats".to_owned()],
        None,
        Some(Duration::from_secs(60)),
    )?;
    let payload: Payload = serde_json::from_value(data).map_err(protocol)?;
    Ok(payload.formats)
}

/// Identify a result file without parsing its contents.
///
/// Cheap by design: mzLib resolves the type and stops, so identifying a million-row file costs no
/// more than identifying an empty one. It is not, however, *pure* — mzLib disambiguates a bare
/// `.tsv` by reading its first line, a `.mztab` by its first five, and a Bruker `.d` by which
/// analysis file the directory holds.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the path is blank, does not exist, or is not a file type mzLib
/// recognises — mzLib has no "unknown" result, so a file is dispatchable or it is an error.
pub fn identify(path: impl AsRef<Path>) -> Result<FileInfo> {
    let args = window_args("identify", path.as_ref(), &ReadOptions::default())?;
    let data = bridge::invoke(&args, None, Some(Duration::from_secs(60)))?;
    serde_json::from_value(data).map_err(protocol)
}

/// Read a result file into the uniform `quantifiable` record view.
///
/// Only the four file types offering that view can be read this way. Use [`read_records`] for any
/// other format.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the path is blank, missing, unrecognised, or has no `quantifiable`
/// view — the message names the views it does have.
pub fn read_results(path: impl AsRef<Path>) -> Result<ResultRecords> {
    read_results_with(path, &ReadOptions::default())
}

/// [`read_results`], with an explicit window.
///
/// # Errors
///
/// As [`read_results`].
pub fn read_results_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<ResultRecords> {
    #[derive(Deserialize)]
    struct Extra {
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        retention_time_unit: String,
        #[serde(default)]
        rows_not_read: Option<i64>,
        #[serde(default)]
        caveats: Vec<String>,
    }

    let args = window_args("read-results", path.as_ref(), options)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    let common: Common = serde_json::from_value(data.clone()).map_err(protocol)?;
    let extra: Extra = serde_json::from_value(data).map_err(protocol)?;
    Ok(ResultRecords::build(
        common,
        extra.retention_time_unit,
        extra.rows_not_read,
        extra.caveats,
    ))
}

/// Read **any** file mzLib recognises, into that format's own fields.
///
/// The exhaustive verb: if [`identify`] succeeds on a path, this reads it. All 31 file types,
/// including the 13 that belong to no cross-format view at all — TopPIC, Crux, MSFragger's peptide
/// and protein tables, the FlashDeconv formats — which no other function here can touch.
///
/// The columns are **not uniform**; see [`NativeRecords`].
///
/// # Errors
///
/// [`MzLibError::Usage`] if the path is blank, missing, or not a file type mzLib recognises.
pub fn read_records(path: impl AsRef<Path>) -> Result<NativeRecords> {
    read_records_with(path, &ReadOptions::default())
}

/// [`read_records`], with an explicit window.
///
/// # Errors
///
/// As [`read_records`].
pub fn read_records_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<NativeRecords> {
    #[derive(Deserialize)]
    struct Extra {
        #[serde(default)]
        reader: Option<String>,
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        record_type: String,
        #[serde(default)]
        views: Vec<String>,
        #[serde(default)]
        excluded_fields: Vec<ExcludedField>,
        #[serde(default)]
        failed_fields: Vec<String>,
    }

    let args = window_args("read-records", path.as_ref(), options)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    let common: Common = serde_json::from_value(data.clone()).map_err(protocol)?;
    let extra: Extra = serde_json::from_value(data).map_err(protocol)?;
    Ok(NativeRecords::build(
        common,
        extra.reader,
        extra.record_type,
        extra.views,
        extra.excluded_fields,
        extra.failed_fields,
    ))
}

/// Read deconvolved MS1 features, in the cross-format `ms1_features` view.
///
/// Two file types offer it: TopFD/FLASHDeconv `_ms1.feature` and Dinosaur `.feature.tsv`.
///
/// **One row is not one line of the file for `_ms1.feature`**: mzLib expands each deconvolved
/// feature into one single-charge feature per charge in its recorded range. Dinosaur is
/// one-for-one. Both facts are in [`FeatureRecords::caveats`].
///
/// # Errors
///
/// [`MzLibError::Usage`] if the file has no `ms1_features` view — the message names the views it
/// does have, and points at [`read_records`].
pub fn read_features(path: impl AsRef<Path>) -> Result<FeatureRecords> {
    read_features_with(path, &ReadOptions::default())
}

/// [`read_features`], with an explicit window.
///
/// # Errors
///
/// As [`read_features`].
pub fn read_features_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<FeatureRecords> {
    #[derive(Deserialize)]
    struct Extra {
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        retention_time_unit: String,
        #[serde(default)]
        caveats: Vec<String>,
    }

    let args = window_args("read-features", path.as_ref(), options)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    let common: Common = serde_json::from_value(data.clone()).map_err(protocol)?;
    let extra: Extra = serde_json::from_value(data).map_err(protocol)?;
    Ok(FeatureRecords::build(
        common,
        extra.retention_time_unit,
        extra.caveats,
    ))
}

/// Read identifications, in the cross-format `spectral_match` view.
///
/// Four file types offer it: MsPathFinderT's targets, decoys and combined results, and Casanovo's
/// `.mztab`.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the file has no `spectral_match` view.
pub fn read_matches(path: impl AsRef<Path>) -> Result<MatchRecords> {
    read_matches_with(path, &ReadOptions::default())
}

/// [`read_matches`], with an explicit window.
///
/// # Errors
///
/// As [`read_matches`].
pub fn read_matches_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<MatchRecords> {
    #[derive(Deserialize)]
    struct Extra {
        #[serde(default)]
        caveats: Vec<String>,
    }

    let args = window_args("read-matches", path.as_ref(), options)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    let common: Common = serde_json::from_value(data.clone()).map_err(protocol)?;
    let extra: Extra = serde_json::from_value(data).map_err(protocol)?;
    Ok(MatchRecords::build(common, extra.caveats))
}

/// Read the scans of a spectra file: headers always, peaks on request.
///
/// Seven file types offer the `spectra` view. **Two of them need Windows**: Bruker `.d` and
/// timsTOF `.d` are read through vendor native libraries and are Windows-x64 only.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the file has no `spectra` view, or `ms_order` is zero.
pub fn read_spectra(path: impl AsRef<Path>) -> Result<ScanRecords> {
    read_spectra_with(path, &SpectraOptions::default())
}

/// [`read_spectra`], with an explicit window, MS-level filter and peak choice.
///
/// # Errors
///
/// As [`read_spectra`].
pub fn read_spectra_with(path: impl AsRef<Path>, options: &SpectraOptions) -> Result<ScanRecords> {
    #[derive(Deserialize)]
    struct Extra {
        #[serde(default)]
        reader: Option<String>,
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        scan_count: u64,
        #[serde(default)]
        ms_order: Option<i64>,
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        peaks_included: bool,
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        retention_time_unit: String,
        #[serde(default)]
        caveats: Vec<String>,
    }

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

    let data = bridge::invoke(&args, None, options.read.timeout)?;
    let common: Common = serde_json::from_value(data.clone()).map_err(protocol)?;
    let extra: Extra = serde_json::from_value(data).map_err(protocol)?;
    Ok(ScanRecords::build(
        common,
        extra.reader,
        extra.scan_count,
        extra.ms_order,
        extra.peaks_included,
        extra.retention_time_unit,
        extra.caveats,
    ))
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

    #[test]
    fn seconds_convert_and_unknown_refuses() {
        let seconds = ResultRecords::build(
            Common {
                file_type: "MsFraggerPsm".to_owned(),
                column_names: vec!["retention_time".to_owned()],
                columns: Some(
                    [("retention_time".to_owned(), vec![Value::from(120.0)])]
                        .into_iter()
                        .collect(),
                ),
                ..Common::default()
            },
            "seconds".to_owned(),
            None,
            vec![],
        );
        assert_eq!(
            seconds.retention_time_in_minutes().unwrap(),
            vec![Some(2.0)]
        );

        let unknown = ResultRecords::build(
            Common {
                file_type: "Ms1Feature".to_owned(),
                column_names: vec!["retention_time".to_owned()],
                columns: Some(
                    [("retention_time".to_owned(), vec![Value::from(2372.27)])]
                        .into_iter()
                        .collect(),
                ),
                ..Common::default()
            },
            "unknown".to_owned(),
            None,
            vec![],
        );
        // Raised rather than guessed: mzLib's own deconvolution code guesses here, and this does
        // not.
        let error = unknown.retention_time_in_minutes().unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(error.to_string().contains("no basis to say"), "{error}");
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
