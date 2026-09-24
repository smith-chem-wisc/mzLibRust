//! SDRF-Proteomics experimental-design files: read one, pool several, and ask mzLib whether
//! they are well-formed, consistent and informative.
//!
//! Every other reader in this crate answers *what did the search find*. SDRF answers *what was
//! searched* — which sample, which organism part, which replicate, which instrument settings — and
//! that is the half you need to group results across experiments.
//!
//! ```
//! # mzlib_replay::activate();
//! let doc = mzlib::sdrf::read("PXD000070.sdrf.tsv")?;
//! assert_eq!((doc.row_count, doc.columns.len()), (6, 31));
//! let organisms = doc.value("characteristics[organism]");              // Vec<Option<&str>>
//! assert_eq!(organisms[0], Some("plasmodium falciparum"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Pool several experiments into one analysis table, naming each one yourself:
//!
//! ```
//! # mzlib_replay::activate();
//! use std::path::PathBuf;
//! use mzlib::sdrf::{pool_with, PoolInput, PoolOptions, ReadOptions};
//!
//! let documents = PoolInput::Labelled(vec![
//!     (PathBuf::from("PXD000070.sdrf.tsv"), "malaria".to_owned()),
//!     (PathBuf::from("PXD026824.sdrf.tsv"), "colon".to_owned()),
//! ]);
//! let options = PoolOptions {
//!     read: ReadOptions { limit: Some(4), ..Default::default() },
//!     ..Default::default()
//! };
//! let pooled = pool_with(&documents, &options)?;
//! assert_eq!((pooled.document_count, pooled.document.row_count), (2, 24));
//! assert_eq!(pooled.labels, ["malaria", "colon"]);
//! assert!(pooled.document.truncated);                 // 4 of the 24 rows came back
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Use this, not [`crate::readers::read_records`], for SDRF.** `read_records` recognises SDRF,
//! but it joins each row's cells into one semicolon-separated string, and SDRF's own key=value
//! grammar (`NT=Oxidation;AC=UNIMOD:35`) puts semicolons *inside* cells. The joined string cannot
//! be split back apart.
//!
//! **This module is row-major, and every other reader here is columnar.** That is not a style
//! choice. [`crate::readers::Table`] is keyed by column name because those names are a schema.
//! SDRF's are *data*, and they **repeat**: 649 files in the curated corpus carry
//! `comment[modification parameters]` more than once, up to eight times in one file, and one file
//! repeats an *empty* name 23 times. A map keyed by name would silently keep one occurrence and
//! drop the rest. So [`SdrfDocument::columns`] is a `Vec` that may contain duplicates,
//! [`SdrfDocument::rows`] is a `Vec` of cell `Vec`s, and position is what links them. Use
//! [`SdrfDocument::value`] for the first cell under a name and [`SdrfDocument::all`] for every one.
//!
//! **Three things worth knowing before you index anything:**
//!
//! * *Rows are ragged.* A row may be shorter than `columns`. Real files are like this — PXD059974
//!   in mzLib's own fixtures has a 46-column header with 17 of its 22 rows carrying 42 cells — and
//!   mzLib preserves it rather than padding, so the file round-trips byte for byte.
//!   [`SdrfDocument::value`] returns `None` for a position a row does not reach.
//! * *Cells are raw strings, never interpreted.* The key=value grammar arrives exactly as written.
//!   It is not decoded, because it cannot be told apart from a cell that merely contains `=` and
//!   `;` — `comment[file uri]` routinely carries pre-signed download URLs whose query strings
//!   contain `Signature=` and `Expires=`.
//! * *A reserved word is a real value.* `"not available"` and `"not applicable"` mean *the
//!   experiment stated an absence*, which is not the same as a column the document does not have.
//!   `None` means the latter. Do not collapse the two.
//!
//! ## Three questions about a document, and one about ages
//!
//! Each is mzLib's own answer (mzLib 1.0.592), projected once in the bridge for all three
//! bindings, never reimplemented here:
//!
//! | question | function | mzLib |
//! |---|---|---|
//! | *Is this file well-formed?* | [`validate`], [`validate_many`] | `SdrfValidator` |
//! | *Do these files write the same thing the same way?* | [`lint_labelled`], [`lint`] | `SdrfDriftLint` |
//! | *Does this file describe its samples at all?* | [`assess`], [`assess_many`] | `SdrfSampleInformativeness` |
//! | *What does it say about each sample?* | [`samples`], [`samples_many`] | `SdrfSampleBlock` |
//! | *How old, in years?* | [`parse_ages`] | `SdrfAge.TryParse` |
//!
//! **The first three are blind in different places, which is why there are three.** A file of
//! `"not available"` validates cleanly and lints clean, and only [`assess`] sees that it says
//! nothing. Twenty individually valid files can still disagree, and only [`lint_labelled`] looks
//! across them. Each result's `caveats` names the blind spot it leaves.
//!
//! ```
//! # mzlib_replay::activate();
//! let cohort = "sdrf_cohort.sdrf.tsv";
//! assert!(mzlib::sdrf::validate(cohort)?.is_valid);                 // well-formed...
//! assert_eq!(mzlib::sdrf::assess(cohort)?.verdict, "Informative");  // ...and says something
//! let drift = mzlib::sdrf::lint_labelled(&[
//!     (cohort, "cohort"),
//!     ("sdrf_cohort_partner.sdrf.tsv", "partner"),
//! ])?;
//! assert_eq!(drift.finding_count, 4);                                // ...but not like its partner
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! [`validate_many`], [`assess_many`] and [`samples_many`] read a whole corpus in **one** bridge
//! call, [`BulkOptions::threads`] documents at a time, and return one long table whose first two
//! columns say which document each row came from. The answer is identical at any thread count;
//! there is deliberately no Rust-side loop, because parallelism is the bridge's to express once.
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::sdrf::{assess_many, BulkOptions};
//!
//! let corpus = ["sdrf_cohort.sdrf.tsv", "sdrf_skeleton.sdrf.tsv", "PXD000070.sdrf.tsv"];
//! let batch = assess_many(&corpus, &BulkOptions::default())?;
//! assert_eq!(batch.paths_with(&["Skeleton"])?, ["sdrf_skeleton.sdrf.tsv"]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Ages are read into **years** — a month is 1/12, a week 7/365.25, a day 1/365.25 — and a cell
//! that would need a guess is refused with its reason, never assumed. A bare `63` is refused,
//! because 63 years and 63 days are both plausible in one study.
//!
//! Ported from pyMzLib's `pymzlib.sdrf`, which decided the verbs, the wire fields and the caveats.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};
use crate::readers::Table;

pub use crate::bridge::OnError;

/// The column [`pool`] adds to record which document each row came from.
pub const SOURCE_DOCUMENT_COLUMN: &str = "comment[source document]";

/// The default timeout, matching pyMzLib's `sdrf` module. SDRF files are small.
const DEFAULT_TIMEOUT: Duration = Duration::from_secs(60);

// ---------------------------------------------------------------------------------------------
// The document
// ---------------------------------------------------------------------------------------------

/// One SDRF-Proteomics document: an ordered header, and rows of raw cells.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct SdrfDocument {
    /// The path that was read. Empty for a pooled table, which is not a file that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// The column names, **verbatim and in document order**. Names may repeat, and the order is
    /// part of the document. Names are never case-normalised: the corpus contains
    /// `"comment[MS min charge]"` and `"Material Type"`, and rewriting them would silently alter
    /// someone's file.
    #[serde(
        rename = "column_names",
        default,
        deserialize_with = "bridge::null_to_default"
    )]
    pub columns: Vec<String>,
    /// One `Vec` of cells per row. **Ragged**: a row may be shorter than [`Self::columns`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub rows: Vec<Vec<String>>,
    /// Data rows in the **whole document** (the header is not a row), regardless of any limit or
    /// offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Data rows actually carried back in [`Self::rows`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// The offset that was applied, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether rows were left behind**, by either the limit or the offset. A short answer and a
    /// complete one must never look alike, so check this rather than comparing counts yourself.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// What this document's data cannot tell you about itself — raggedness, repeated names,
    /// reserved words. Worth printing the first time you read an unfamiliar file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfDocument {
    /// The position of the first column with this name, or `None` if the document has none.
    ///
    /// Comparison is exact and case-**sensitive**, matching the SDRF specification and mzLib.
    pub fn index_of(&self, column: &str) -> Option<usize> {
        self.columns.iter().position(|name| name == column)
    }

    /// Every position carrying this name, in document order. Empty when the column is absent.
    pub fn indexes_of(&self, column: &str) -> Vec<usize> {
        self.columns
            .iter()
            .enumerate()
            .filter(|(_, name)| *name == column)
            .map(|(i, _)| i)
            .collect()
    }

    /// The first cell under `column`, one entry per returned row.
    ///
    /// `None` means **the document does not have this column, or the row is too short to reach
    /// it**. It does not mean "empty": the SDRF reserved words `"not available"` and
    /// `"not applicable"` are real values that an experiment chose to write, and they come back as
    /// themselves.
    pub fn value(&self, column: &str) -> Vec<Option<&str>> {
        match self.index_of(column) {
            None => vec![None; self.rows.len()],
            Some(i) => self
                .rows
                .iter()
                .map(|row| row.get(i).map(String::as_str))
                .collect(),
        }
    }

    /// Every cell under `column`, one `Vec` per returned row.
    ///
    /// The accessor for a multi-cardinality column such as `comment[modification parameters]`,
    /// which legitimately repeats — up to eight times in one corpus file. Positions a row is too
    /// short to reach are skipped rather than reported, so each inner `Vec` holds only cells that
    /// exist.
    pub fn all(&self, column: &str) -> Vec<Vec<&str>> {
        let indexes = self.indexes_of(column);
        self.rows
            .iter()
            .map(|row| {
                indexes
                    .iter()
                    .filter_map(|&i| row.get(i).map(String::as_str))
                    .collect()
            })
            .collect()
    }

    /// The rows as name-to-cell maps, for the common case of a document with no repeated names.
    ///
    /// **Lossy when a name repeats** — a later position overwrites an earlier one — which is
    /// exactly why it is not the primary shape. [`Self::has_repeated_columns`] says whether that
    /// applies to this document; for a repeating document use [`Self::all`] instead. A position a
    /// row is too short to reach is absent from that row's map.
    pub fn records(&self) -> Vec<HashMap<&str, &str>> {
        self.rows
            .iter()
            .map(|row| {
                self.columns
                    .iter()
                    .zip(row)
                    .map(|(name, cell)| (name.as_str(), cell.as_str()))
                    .collect()
            })
            .collect()
    }

    /// Whether any column name appears more than once — see [`Self::records`].
    pub fn has_repeated_columns(&self) -> bool {
        let mut seen = std::collections::HashSet::with_capacity(self.columns.len());
        !self.columns.iter().all(|name| seen.insert(name))
    }

    /// How many returned rows carry fewer cells than there are columns.
    pub fn ragged_row_count(&self) -> usize {
        self.rows
            .iter()
            .filter(|row| row.len() < self.columns.len())
            .count()
    }
}

/// Where [`pool_with`] wrote the merged document, when asked to write one.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct WrittenSdrf {
    /// The path written.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Rows written — the **whole** merged document, not the windowed slice. The limit and offset
    /// shape what comes back over the wire; they never shorten the file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
}

/// Several SDRF documents merged into one table.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct PooledSdrf {
    /// The merged table. Its `path` is empty: a pooled table is not a file that was read.
    #[serde(flatten)]
    pub document: SdrfDocument,
    /// How many documents were pooled.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub document_count: u64,
    /// The paths pooled, in the order given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub paths: Vec<String>,
    /// The provenance label used for each, in the same order — either what you supplied or
    /// mzLib's `containing-folder/file-stem` default.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub labels: Vec<String>,
    /// Where the merged document was written, or `None` if [`PoolOptions::out`] was not set.
    #[serde(default)]
    pub written: Option<WrittenSdrf>,
}

impl PooledSdrf {
    /// The provenance label of each returned row — which document it came from.
    pub fn source_documents(&self) -> Vec<Option<&str>> {
        self.document.value(SOURCE_DOCUMENT_COLUMN)
    }
}

// ---------------------------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------------------------

/// How much of one document to return.
#[derive(Debug, Clone)]
pub struct ReadOptions {
    /// Return at most this many rows. `None` returns all of them. `Some(0)` returns the header
    /// alone.
    pub limit: Option<u64>,
    /// Skip this many rows.
    pub offset: u64,
    /// Time to allow. `None` waits indefinitely. Defaults to 60 seconds.
    pub timeout: Option<Duration>,
}

impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            limit: None,
            offset: 0,
            timeout: Some(DEFAULT_TIMEOUT),
        }
    }
}

/// How much of a pooled table to return, and where to write the whole of it.
#[derive(Debug, Clone, Default)]
pub struct PoolOptions {
    /// The window and the timeout, as for [`read_with`].
    pub read: ReadOptions,
    /// Write the merged document here as SDRF. The **whole** document is written regardless of
    /// the limit and offset.
    pub out: Option<String>,
}

/// The documents to pool: bare paths, or each path with the provenance label to stamp on its rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PoolInput {
    /// Paths only. mzLib falls back to `containing-folder/file-stem` for provenance, which depends
    /// on where the files happen to sit, so the result carries a caveat saying it is not
    /// reproducible elsewhere. Prefer [`PoolInput::Labelled`].
    Paths(Vec<PathBuf>),
    /// Each path with the label you chose for it. Every label must be non-blank.
    Labelled(Vec<(PathBuf, String)>),
}

// ---------------------------------------------------------------------------------------------
// Argument assembly
// ---------------------------------------------------------------------------------------------

fn path_text<'a>(path: &'a Path, what: &str) -> Result<&'a str> {
    let text = path.to_str().ok_or_else(|| {
        MzLibError::Usage(format!(
            "The {what} is not valid UTF-8, which the bridge requires."
        ))
    })?;
    Ok(text.trim())
}

fn window_args(args: &mut Vec<String>, options: &ReadOptions) {
    if let Some(limit) = options.limit {
        args.push("--limit".to_owned());
        args.push(limit.to_string());
    }
    if options.offset > 0 {
        args.push("--offset".to_owned());
        args.push(options.offset.to_string());
    }
}

fn read_args(path: &Path, options: &ReadOptions) -> Result<Vec<String>> {
    let path = path_text(path, "file path")?;
    if path.is_empty() {
        return Err(MzLibError::Usage(
            "A file path is required, e.g. 'PXD000070.sdrf.tsv'.".to_owned(),
        ));
    }
    let mut args = vec![
        "sdrf".to_owned(),
        "read".to_owned(),
        "--path".to_owned(),
        path.to_owned(),
    ];
    window_args(&mut args, options);
    Ok(args)
}

/// The arguments and the stdin for a pool, validated before anything is spawned.
///
/// Documents cross on stdin, one per line, as `path[\tlabel]` — the bridge keeps named options in
/// a dictionary, so a repeated `--path` would silently keep only the last.
fn pool_request(documents: &PoolInput, options: &PoolOptions) -> Result<(Vec<String>, String)> {
    let lines = document_lines(documents)?;

    let mut args = vec!["sdrf".to_owned(), "pool".to_owned()];
    window_args(&mut args, &options.read);
    if let Some(out) = &options.out {
        if out.trim().is_empty() {
            return Err(MzLibError::Usage(
                "out must be a non-empty path, or None to skip writing.".to_owned(),
            ));
        }
        args.push("--out".to_owned());
        args.push(out.trim().to_owned());
    }

    Ok((args, lines.join("\n")))
}

/// The `path[\tlabel]` stdin lines `sdrf pool` and `sdrf lint` share, validated.
fn document_lines(documents: &PoolInput) -> Result<Vec<String>> {
    let pairs: Vec<(&Path, Option<&str>)> = match documents {
        PoolInput::Paths(paths) => paths.iter().map(|p| (p.as_path(), None)).collect(),
        PoolInput::Labelled(pairs) => pairs
            .iter()
            .map(|(p, label)| (p.as_path(), Some(label.as_str())))
            .collect(),
    };
    if pairs.is_empty() {
        return Err(MzLibError::Usage(
            "At least one SDRF document is required.".to_owned(),
        ));
    }

    let mut lines = Vec::with_capacity(pairs.len());
    for (path, label) in pairs {
        let path = path_text(path, "document path")?;
        if path.is_empty() {
            return Err(MzLibError::Usage(
                "A document path may not be blank.".to_owned(),
            ));
        }
        let label = label.map(str::trim);
        if label.is_some_and(str::is_empty) {
            return Err(MzLibError::Usage(format!(
                "The label for '{path}' is blank. Give every document a label, or use \
                 PoolInput::Paths to accept mzLib's path-derived default for all of them."
            )));
        }
        let separator = |text: &str| text.contains(['\t', '\n', '\r']);
        if separator(path) || label.is_some_and(separator) {
            return Err(MzLibError::Usage(format!(
                "A path or label contains a tab or newline, which the bridge uses to separate \
                 them: {path:?} / {:?}.",
                label.unwrap_or("")
            )));
        }
        lines.push(match label {
            Some(label) => format!("{path}\t{label}"),
            None => path.to_owned(),
        });
    }
    Ok(lines)
}

// ---------------------------------------------------------------------------------------------
// The public surface
// ---------------------------------------------------------------------------------------------

/// Read one SDRF-Proteomics file, every row.
///
/// See [`read_with`] for the reference.
///
/// # Errors
///
/// As [`read_with`].
pub fn read(path: impl AsRef<Path>) -> Result<SdrfDocument> {
    read_with(path, &ReadOptions::default())
}

/// Read one SDRF-Proteomics file with every cell intact, in a row-major shape that keeps the
/// document's repeated column names.
///
/// A blank path is refused before anything is spawned.
#[doc = include_str!("../docs/reference/sdrf.read.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// // PXD059974: a 46-column header over rows of 42 cells. mzLib keeps the raggedness.
/// let doc = mzlib::sdrf::read("PXD059974.sdrf.tsv")?;
/// assert_eq!(doc.ragged_row_count(), 17);
/// let last = doc.columns.last().unwrap();
/// assert!(doc.value(last).iter().any(Option::is_none)); // a short row reaches no cell there
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.read.see-also.md")]
pub fn read_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<SdrfDocument> {
    let args = read_args(path.as_ref(), options)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    serde_json::from_value(data).map_err(protocol)
}

/// Merge several SDRF documents into one analysis table, with provenance from their paths.
///
/// Prefer [`pool_labelled`]: see [`PoolInput::Paths`] for why.
///
/// # Errors
///
/// As [`pool_with`].
pub fn pool<P: AsRef<Path>>(paths: &[P]) -> Result<PooledSdrf> {
    pool_with(
        &PoolInput::Paths(paths.iter().map(|p| p.as_ref().to_path_buf()).collect()),
        &PoolOptions::default(),
    )
}

/// Merge several SDRF documents into one analysis table, naming each one yourself.
///
/// # Errors
///
/// As [`pool_with`].
pub fn pool_labelled<P: AsRef<Path>, L: AsRef<str>>(documents: &[(P, L)]) -> Result<PooledSdrf> {
    pool_with(
        &PoolInput::Labelled(
            documents
                .iter()
                .map(|(p, l)| (p.as_ref().to_path_buf(), l.as_ref().to_owned()))
                .collect(),
        ),
        &PoolOptions::default(),
    )
}

/// Merge several SDRF documents into one analysis table, with provenance.
///
/// Columns are the union of every document's, ordered by SDRF's own block structure, and a name
/// that repeats is carried at the highest multiplicity any single document used, so nothing is
/// dropped. A cell a document did not have is filled with the reserved word `"not available"`,
/// and a [`SOURCE_DOCUMENT_COLUMN`] records which document each row came from.
///
/// **The result is an analysis table, not something to deposit.** `source name` + `assay name` +
/// `comment[label]` is unique within one document, but two experiments may both have a
/// `"Sample 1"`, so a pooled table will usually violate SDRF's uniqueness rule. Use the
/// source-document column as part of any key.
///
/// No documents, a blank path or label, or a path or label containing a tab or newline, is refused
/// before anything is spawned.
#[doc = include_str!("../docs/reference/sdrf.pool.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{pool_with, PoolInput, PoolOptions, ReadOptions, SOURCE_DOCUMENT_COLUMN};
///
/// let pooled = pool_with(
///     &PoolInput::Labelled(vec![
///         ("PXD000070.sdrf.tsv".into(), "malaria".to_owned()),
///         ("PXD026824.sdrf.tsv".into(), "colon".to_owned()),
///     ]),
///     &PoolOptions { read: ReadOptions { limit: Some(4), ..Default::default() }, out: None },
/// )?;
/// assert!(pooled.document.columns.iter().any(|c| c == SOURCE_DOCUMENT_COLUMN));
/// assert_eq!(pooled.source_documents()[0], Some("malaria"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.pool.see-also.md")]
pub fn pool_with(documents: &PoolInput, options: &PoolOptions) -> Result<PooledSdrf> {
    let (args, stdin) = pool_request(documents, options)?;
    let data = bridge::invoke(&args, Some(&stdin), options.read.timeout)?;
    serde_json::from_value(data).map_err(protocol)
}

// ---------------------------------------------------------------------------------------------
// Validate, lint, assess, samples and ages (mzLib 1.0.592)
// ---------------------------------------------------------------------------------------------

/// The three verdicts [`assess`] can return, best first, as mzLib spells them
/// (`SdrfSampleVerdict`).
pub const VERDICTS: [&str; 3] = ["Informative", "Partial", "Skeleton"];

/// How much of an age a cell pins down, as mzLib spells it (`SdrfAgePrecision`). See
/// [`parse_ages`].
pub const AGE_PRECISIONS: [&str; 4] = ["Exact", "Range", "LowerBound", "UpperBound"];

/// Why an age cell was refused. mzLib's `SdrfAge.TryParse` answers only yes or no; the bridge
/// names the reason from the cell: `empty`; `reserved_word` (`"not available"` and friends, any
/// case); `no_unit` (a bare number such as `63`, where years and days cannot be told apart);
/// `unreadable`.
pub const AGE_REFUSALS: [&str; 4] = ["empty", "reserved_word", "no_unit", "unreadable"];

/// How a many-document call reads its corpus: how many documents at once, and what one
/// unreadable document does to the rest (BULK.md §1).
///
/// The whole list goes to **one** bridge process, which reads `threads` documents at once and
/// returns one table in input order. There is deliberately no Rust-side loop or thread pool:
/// parallelism is the bridge's to express, once, where it can be counted.
#[derive(Debug, Clone)]
pub struct BulkOptions {
    /// Documents read at once: 1 or more, or `-1` for one per core. The default is 1. The answer is
    /// byte-identical at any value (tested at 1, 4 and -1), so this is a resource choice only.
    pub threads: i32,
    /// [`OnError::Fail`] (the default) returns the first unreadable document's error, in input
    /// order. [`OnError::Skip`] records the failure in that document's `files` entry and reads
    /// the rest.
    pub on_error: OnError,
    /// Time to allow for the whole batch. `None`, the default, waits indefinitely.
    pub timeout: Option<Duration>,
}

impl Default for BulkOptions {
    fn default() -> Self {
        Self {
            threads: 1,
            on_error: OnError::Fail,
            timeout: None,
        }
    }
}

/// Why one document of a `*_many` call produced no result, under [`OnError::Skip`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct FileError {
    /// `"usage"` (the file is missing) or `"correctness"` (mzLib could not read it; the message
    /// then starts with the .NET exception type, e.g. `"MzLibException: ..."`).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub kind: String,
    /// What went wrong, naming the file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub message: String,
}

// ---- validate -------------------------------------------------------------------------------

/// One finding from mzLib's `SdrfValidator`, located as precisely as its rule allows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidationMessage {
    /// `"Error"` (the document cannot be reliably consumed: a ragged row, a missing required
    /// column, two indistinguishable rows) or `"Warning"` (it deviates from the specification but
    /// can still be joined: casing, an empty cell, a missing recommended column).
    pub severity: String,
    /// The stable rule id, e.g. `"RequiredColumn"`, `"RowWidth"`, `"RowKeyUniqueness"`,
    /// `"ReservedWordCase"`. Stable, so you can count or suppress by it.
    pub rule: String,
    /// Human-readable, including the offending value.
    pub message: String,
    /// 0-based index into the document's data rows; `None` for a finding about the whole
    /// document (a missing column, an empty header, no rows).
    pub row_index: Option<i64>,
    /// The 1-based line in the file, `row_index + 2` because the header is line 1; `None` exactly
    /// when `row_index` is.
    pub line_number: Option<i64>,
    /// The column involved; `None` when the finding is not about one column.
    pub column_name: Option<String>,
}

fn messages(table: &Table) -> Result<Vec<ValidationMessage>> {
    let severity = table.strings("severity")?;
    let rule = table.strings("rule")?;
    let message = table.strings("message")?;
    let row_index = table.integers("row_index")?;
    let line_number = table.integers("line_number")?;
    let column_name = table.strings("column_name")?;
    Ok((0..rule.len())
        .map(|i| ValidationMessage {
            severity: severity[i].clone().unwrap_or_default(),
            rule: rule[i].clone().unwrap_or_default(),
            message: message[i].clone().unwrap_or_default(),
            row_index: row_index[i],
            line_number: line_number[i],
            column_name: column_name[i].clone(),
        })
        .collect())
}

/// The structural findings for one SDRF document — what [`validate`] returns.
///
/// The findings are a table, one row per finding: `severity`, `rule`, `message`, `row_index`,
/// `line_number` and `column_name`, in mzLib's order. [`SdrfValidation::messages`] gives the same
/// rows as [`ValidationMessage`] values.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfValidation {
    /// The absolute path validated.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `true` when there is no `Error` finding. **Warnings never make a document invalid**: mzLib
    /// calibrated every severity against the 1,236-file curated corpus, and a rule that fired on
    /// most curated files was judged wrong rather than the files.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub is_valid: bool,
    /// Findings of severity `Error`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub error_count: u64,
    /// Findings of severity `Warning`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub warning_count: u64,
    /// All findings: the rows of [`Self::columns`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub message_count: u64,
    /// Data rows in the document (the header is not a row).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The findings table, one row per finding, in mzLib's order.
    #[serde(flatten)]
    pub columns: Table,
    /// What a clean result does *not* mean. Read these once.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfValidation {
    /// Every finding as a [`ValidationMessage`], in mzLib's order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn messages(&self) -> Result<Vec<ValidationMessage>> {
        messages(&self.columns)
    }

    /// The findings that make the document invalid.
    ///
    /// # Errors
    ///
    /// As [`Self::messages`].
    pub fn errors(&self) -> Result<Vec<ValidationMessage>> {
        Ok(self
            .messages()?
            .into_iter()
            .filter(|m| m.severity == "Error")
            .collect())
    }

    /// The findings worth fixing that do not make the document invalid.
    ///
    /// # Errors
    ///
    /// As [`Self::messages`].
    pub fn warnings(&self) -> Result<Vec<ValidationMessage>> {
        Ok(self
            .messages()?
            .into_iter()
            .filter(|m| m.severity == "Warning")
            .collect())
    }
}

/// One document's summary in a [`validate_many`] result. Every count is `None` when the document
/// could not be read, and [`Self::error`] then says why.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct ValidatedFile {
    /// The path as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// No `Error` findings; `None` if unread.
    #[serde(default)]
    pub is_valid: Option<bool>,
    /// `Error` findings; `None` if unread.
    #[serde(default)]
    pub error_count: Option<u64>,
    /// `Warning` findings; `None` if unread.
    #[serde(default)]
    pub warning_count: Option<u64>,
    /// All findings; `None` if unread.
    #[serde(default)]
    pub message_count: Option<u64>,
    /// Data rows in the document; `None` if unread.
    #[serde(default)]
    pub row_count: Option<u64>,
    /// Why it was not read, or `None` when it was.
    #[serde(default)]
    pub error: Option<FileError>,
}

/// The findings for many SDRF documents as one long table — what [`validate_many`] returns.
///
/// The table has the columns of [`SdrfValidation`], preceded by `source_index` (the document's
/// 0-based position in the list you passed) and `source_path`. Rows are in input order, then
/// mzLib's order within a document, whatever [`BulkOptions::threads`] was.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfValidationBatch {
    /// Documents given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Documents read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Documents that could not be read; non-zero only under [`OnError::Skip`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// Rows in [`Self::columns`] over every document read — the count every bulk verb reports
    /// (BULK.md §2). Equal to [`Self::message_count`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Documents read with no `Error`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub valid_count: u64,
    /// Findings over every document read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub message_count: u64,
    /// One entry per input, in input order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<ValidatedFile>,
    /// The findings table: `source_index`, `source_path`, then the single-document columns.
    #[serde(flatten)]
    pub columns: Table,
    /// As for one document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfValidationBatch {
    /// Every finding as a [`ValidationMessage`]; read [`Self::columns`]' `source_index` to know
    /// which document each came from.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn messages(&self) -> Result<Vec<ValidationMessage>> {
        messages(&self.columns)
    }

    /// The inputs that could not be read, under [`OnError::Skip`].
    #[must_use]
    pub fn failed_files(&self) -> Vec<&ValidatedFile> {
        self.files.iter().filter(|f| f.error.is_some()).collect()
    }
}

// ---- lint -----------------------------------------------------------------------------------

/// One spelling of one inconsistently written concept: a row of an [`SdrfDrift`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DriftVariant {
    /// Which finding the row belongs to, 0-based, most impactful first.
    pub finding_index: i64,
    /// `AccessionNameConflict`, `NameAccessionConflict`, `MixedTermAndFreeText`,
    /// `ColumnNameVariant` or `ValueCaseVariant` (mzLib's `SdrfDriftKind`).
    pub kind: String,
    /// What was written inconsistently: an accession, a name, a normalised value, or a normalised
    /// column name, depending on `kind`.
    pub concept: String,
    /// The column whose values drift; `None` for `ColumnNameVariant`, a finding about a name.
    pub column_name: Option<String>,
    /// 0 for the majority spelling, then by frequency. Descriptive, not advice.
    pub variant_rank: i64,
    /// This spelling exactly as written.
    pub value: String,
    /// How widely this spelling was used, in documents.
    pub occurrences: i64,
    /// The labels of the documents that used it.
    pub documents: Vec<String>,
}

/// Concepts a set of SDRF documents wrote inconsistently, from mzLib's `SdrfDriftLint` — what
/// [`lint_labelled`] and [`lint`] return.
///
/// The table is long, **one row per finding × spelling**: a finding with three spellings is three
/// rows sharing a `finding_index`. Its columns are `finding_index`, `kind`, `concept`,
/// `column_name`, `variant_rank`, `value`, `occurrences` and `documents` — the last a **list** per
/// row, since labels are your text and no delimiter could be undone.
/// [`SdrfDrift::findings`] groups the rows.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfDrift {
    /// Documents linted.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub document_count: u64,
    /// The paths, in the order given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub paths: Vec<String>,
    /// The label each document is named by in the `documents` column: yours, or mzLib's
    /// `containing-folder/file-stem` default.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub labels: Vec<String>,
    /// Distinct findings (not rows).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub finding_count: u64,
    /// The findings table, one row per finding × spelling, findings in mzLib's impact order and
    /// spellings majority first.
    #[serde(flatten)]
    pub columns: Table,
    /// Including why the majority spelling is not advice.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfDrift {
    /// Every row as a [`DriftVariant`], in the table's order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn variants(&self) -> Result<Vec<DriftVariant>> {
        let t = &self.columns;
        let finding_index = t.integers("finding_index")?;
        let kind = t.strings("kind")?;
        let concept = t.strings("concept")?;
        let column_name = t.strings("column_name")?;
        let variant_rank = t.integers("variant_rank")?;
        let value = t.strings("value")?;
        let occurrences = t.integers("occurrences")?;
        let documents = t.raw("documents").unwrap_or(&[]);
        (0..finding_index.len())
            .map(|i| {
                let labels = match documents.get(i) {
                    Some(serde_json::Value::Array(items)) => items
                        .iter()
                        .map(|item| item.as_str().map(str::to_owned))
                        .collect::<Option<Vec<_>>>(),
                    Some(serde_json::Value::Null) | None => Some(Vec::new()),
                    Some(_) => None,
                }
                .ok_or_else(|| {
                    MzLibError::Protocol(format!(
                        "Column 'documents' row {i} is not a list of labels."
                    ))
                })?;
                Ok(DriftVariant {
                    finding_index: finding_index[i].unwrap_or_default(),
                    kind: kind[i].clone().unwrap_or_default(),
                    concept: concept[i].clone().unwrap_or_default(),
                    column_name: column_name[i].clone(),
                    variant_rank: variant_rank[i].unwrap_or_default(),
                    value: value[i].clone().unwrap_or_default(),
                    occurrences: occurrences[i].unwrap_or_default(),
                    documents: labels,
                })
            })
            .collect()
    }

    /// The rows grouped by finding: one `Vec` of spellings per finding, majority first.
    ///
    /// # Errors
    ///
    /// As [`Self::variants`].
    pub fn findings(&self) -> Result<Vec<Vec<DriftVariant>>> {
        let mut grouped: Vec<Vec<DriftVariant>> = Vec::new();
        for variant in self.variants()? {
            match grouped.last_mut() {
                Some(group) if group[0].finding_index == variant.finding_index => {
                    group.push(variant);
                }
                _ => grouped.push(vec![variant]),
            }
        }
        Ok(grouped)
    }
}

// ---- assess ---------------------------------------------------------------------------------

/// Whether one SDRF document's sample half describes an experimental design — what [`assess`]
/// returns.
///
/// mzLib's `SdrfSampleInformativeness` asks three questions, and the verdict is how many pass: all
/// three is `Informative`, none is `Skeleton`, anything between is `Partial`. The table is the
/// evidence, one row per column a check read: `role` (`factor_value`, `sample_characteristic` or
/// `biological_replicate`), `column_name`, `rows`, `filled` (rows with a real answer), `absent`
/// (empty or a reserved word), `distinct_values` (compared ignoring case and surrounding space)
/// and `fill_rate` (`filled / rows`, a fraction from 0 to 1, which may cross as a whole number:
/// read it with [`Table::floats`]).
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfAssessment {
    /// The absolute path assessed.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `"Informative"`, `"Partial"` or `"Skeleton"` — see [`VERDICTS`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub verdict: String,
    /// Some `factor value[...]` column holds two or more distinct real answers.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub factor_value_varies: bool,
    /// Some `characteristics[...]` column other than organism and biological replicate holds a
    /// real answer. Organism is left out because a search can fill it without a human.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_is_described: bool,
    /// `characteristics[biological replicate]` holds two or more distinct real answers; `false`
    /// also when the column is missing (the table then has no `biological_replicate` row).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub biological_replicate_varies: bool,
    /// Data rows in the document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The evidence table: factor values, then sample characteristics, then biological
    /// replicate, each in mzLib's order.
    #[serde(flatten)]
    pub columns: Table,
    /// Including when `Partial` is perfectly legitimate.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// How many documents of an [`assess_many`] batch got each verdict.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Deserialize)]
pub struct VerdictCounts {
    /// Documents judged `Informative`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub informative: u64,
    /// Documents judged `Partial`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub partial: u64,
    /// Documents judged `Skeleton`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub skeleton: u64,
}

/// One document's verdict in an [`assess_many`] result; `None` fields mean it was not read.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct AssessedFile {
    /// The path as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `"Informative"`, `"Partial"` or `"Skeleton"`; `None` if unread.
    #[serde(default)]
    pub verdict: Option<String>,
    /// As [`SdrfAssessment::factor_value_varies`]; `None` if unread.
    #[serde(default)]
    pub factor_value_varies: Option<bool>,
    /// As [`SdrfAssessment::sample_is_described`]; `None` if unread.
    #[serde(default)]
    pub sample_is_described: Option<bool>,
    /// As [`SdrfAssessment::biological_replicate_varies`]; `None` if unread.
    #[serde(default)]
    pub biological_replicate_varies: Option<bool>,
    /// Data rows; `None` if unread.
    #[serde(default)]
    pub row_count: Option<u64>,
    /// Why it was not read, or `None`.
    #[serde(default)]
    pub error: Option<FileError>,
}

/// Verdicts for many SDRF documents, with the evidence as one long table — what [`assess_many`]
/// returns: the gate for a corpus.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfAssessmentBatch {
    /// Documents given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Documents assessed.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Documents that could not be read; non-zero only under [`OnError::Skip`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// Rows in [`Self::columns`] over every document read — the count every bulk verb reports
    /// (BULK.md §2).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Documents with each verdict, over the documents read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub verdict_counts: VerdictCounts,
    /// One entry per input, in input order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<AssessedFile>,
    /// The evidence table: `source_index`, `source_path`, then the columns of [`SdrfAssessment`].
    #[serde(flatten)]
    pub columns: Table,
    /// As for one document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfAssessmentBatch {
    /// The paths whose verdict is one of `verdicts`, in input order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] for a verdict that is not one of [`VERDICTS`]: a typo would otherwise
    /// silently select nothing.
    pub fn paths_with(&self, verdicts: &[&str]) -> Result<Vec<&str>> {
        if let Some(unknown) = verdicts.iter().find(|v| !VERDICTS.contains(v)) {
            return Err(MzLibError::Usage(format!(
                "Unknown verdict '{unknown}'; expected some of {VERDICTS:?}."
            )));
        }
        Ok(self
            .files
            .iter()
            .filter(|f| f.verdict.as_deref().is_some_and(|v| verdicts.contains(&v)))
            .map(|f| f.path.as_str())
            .collect())
    }
}

// ---- samples --------------------------------------------------------------------------------

/// One sample's `characteristics[age]`, read into years: a row of [`SdrfSamples::ages`].
#[derive(Debug, Clone, PartialEq)]
pub struct SampleAge {
    /// The sample, as its first row spelled it.
    pub source_name: String,
    /// The cell, verbatim.
    pub value: Option<String>,
    /// The age, a range's midpoint, or a bound, in years; `None` when refused.
    pub years: Option<f64>,
    /// The youngest the cell allows, in years; `None` when refused.
    pub min_years: Option<f64>,
    /// The oldest the cell allows, in years; `None` when refused **or unbounded** (a
    /// `LowerBound` such as `>=90Y`).
    pub max_years: Option<f64>,
    /// One of [`AGE_PRECISIONS`]; `None` when refused.
    pub precision: Option<String>,
    /// `true` for the specification's grammar, `false` for unambiguous words; `None` when refused.
    pub follows_specification: Option<bool>,
    /// `None` when the cell was read, else one of [`AGE_REFUSALS`].
    pub refusal: Option<String>,
}

/// Every sample in one SDRF document and what its rows agree on, as a long table — what
/// [`samples`] returns.
///
/// A sample is a `source name`, matched ignoring case and surrounding space, exactly as mzLib's
/// `SdrfSampleBlock.BySourceName` keys it: one sample measured in three fractions is one sample.
/// The table has one row per sample × column × position, samples in document order, with
/// `source_name`, `sample_row_count`, `column_name` (the header's own spelling), `column_kind`
/// (`source_name`, `characteristic` or `factor_value`), `position` (`None` on a conflicting row),
/// `status` (`agreed`, or `conflicting` when the sample's rows disagree and mzLib withheld the
/// column), `value` (verbatim; `None` when withheld), and the six `age_*` columns, filled only on
/// `characteristics[age]` rows with the meanings [`ParsedAges`] documents.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfSamples {
    /// The absolute path read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Distinct samples (source names).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// Data rows in the document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Sample × column pairs mzLib withheld because the sample's rows disagree.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub conflict_count: u64,
    /// Rows mzLib could not place: no `source name` column, or a blank source name. Empty for a
    /// well-formed document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub problems: Vec<String>,
    /// The samples table.
    #[serde(flatten)]
    pub columns: Table,
    /// Including how to key samples across documents.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

fn sample_ages(table: &Table) -> Result<Vec<SampleAge>> {
    let source_name = table.strings("source_name")?;
    let value = table.strings("value")?;
    let years = table.floats("age_years")?;
    let min_years = table.floats("age_min_years")?;
    let max_years = table.floats("age_max_years")?;
    let precision = table.strings("age_precision")?;
    let follows = table.booleans("age_follows_specification")?;
    let refusal = table.strings("age_refusal")?;
    Ok((0..source_name.len())
        .filter(|&i| years[i].is_some() || refusal[i].is_some())
        .map(|i| SampleAge {
            source_name: source_name[i].clone().unwrap_or_default(),
            value: value[i].clone(),
            years: years[i],
            min_years: min_years[i],
            max_years: max_years[i],
            precision: precision[i].clone(),
            follows_specification: follows[i],
            refusal: refusal[i].clone(),
        })
        .collect())
}

impl SdrfSamples {
    /// Only the parsed-age rows: one per sample whose `characteristics[age]` was agreed.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn ages(&self) -> Result<Vec<SampleAge>> {
        sample_ages(&self.columns)
    }

    /// `(source_name, column_name)` for every column mzLib withheld as conflicting.
    ///
    /// # Errors
    ///
    /// As [`Self::ages`].
    pub fn conflicts(&self) -> Result<Vec<(String, String)>> {
        let status = self.columns.strings("status")?;
        let source_name = self.columns.strings("source_name")?;
        let column_name = self.columns.strings("column_name")?;
        Ok((0..status.len())
            .filter(|&i| status[i].as_deref() == Some("conflicting"))
            .map(|i| {
                (
                    source_name[i].clone().unwrap_or_default(),
                    column_name[i].clone().unwrap_or_default(),
                )
            })
            .collect())
    }
}

/// One document's summary in a [`samples_many`] result; `None` fields mean it was not read.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct SampledFile {
    /// The path as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Distinct samples; `None` if unread.
    #[serde(default)]
    pub sample_count: Option<u64>,
    /// Data rows; `None` if unread.
    #[serde(default)]
    pub row_count: Option<u64>,
    /// Columns withheld as conflicting; `None` if unread.
    #[serde(default)]
    pub conflict_count: Option<u64>,
    /// Rows mzLib could not place; `None` if unread.
    #[serde(default)]
    pub problems: Option<Vec<String>>,
    /// Why it was not read, or `None`.
    #[serde(default)]
    pub error: Option<FileError>,
}

/// The samples of many SDRF documents, as one long table — what [`samples_many`] returns.
///
/// The table has the columns of [`SdrfSamples`], preceded by `source_index` and `source_path`.
/// **A source name is only unique within one document** — two studies may both have a
/// `"Sample 1"` — so key a sample on `(source_index, source_name)`.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfSamplesBatch {
    /// Documents given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Documents read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Documents that could not be read; non-zero only under [`OnError::Skip`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// Rows in [`Self::columns`] over every document read — the count every bulk verb reports
    /// (BULK.md §2).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Samples over every document read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// One entry per input, in input order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<SampledFile>,
    /// The samples table: `source_index`, `source_path`, then the columns of [`SdrfSamples`].
    #[serde(flatten)]
    pub columns: Table,
    /// As for one document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

// ---- parse_ages -----------------------------------------------------------------------------

/// One age cell read into years: a row of [`ParsedAges::ages`].
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedAge {
    /// The cell exactly as given.
    pub cell: String,
    /// The single figure to place the sample on an age axis, in years: the age, a range's
    /// **midpoint**, or a bound. `None` when refused.
    pub years: Option<f64>,
    /// The youngest the cell allows, in years; `0` for an upper bound. `None` when refused.
    pub min_years: Option<f64>,
    /// The oldest the cell allows, in years. `None` when refused **or** when there is no upper
    /// bound (`precision == "LowerBound"`, e.g. `>=90Y`): mzLib's +infinity, which JSON cannot
    /// carry. [`Self::precision`] tells the two apart.
    pub max_years: Option<f64>,
    /// One of [`AGE_PRECISIONS`]; `None` when refused.
    pub precision: Option<String>,
    /// `true` for the specification's `nYnMnD` / `nW` grammar, `false` for unambiguous words
    /// (`"3 year"`, `"6-8 weeks"`); `None` when refused.
    pub follows_specification: Option<bool>,
    /// `None` when the cell was read, else one of [`AGE_REFUSALS`].
    pub refusal: Option<String>,
}

/// Age cells read into years by mzLib's `SdrfAge.TryParse`, one row per cell given, in input
/// order — what [`parse_ages`] returns.
///
/// The table's columns are `cell`, `years`, `min_years`, `max_years`, `precision`,
/// `follows_specification` and `refusal`, with the meanings of [`ParsedAge`]. Every age is in
/// **years**: a month is 1/12 year, a week 7/365.25, a day 1/365.25, an hour 1/(365.25 × 24).
#[derive(Debug, Clone, Deserialize)]
pub struct ParsedAges {
    /// Cells given: the rows of [`Self::columns`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub cell_count: u64,
    /// Cells read as an age.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub parsed_count: u64,
    /// The ages table, one row per input cell.
    #[serde(flatten)]
    pub columns: Table,
    /// The unit conventions and refusal rules.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl ParsedAges {
    /// Every row as a [`ParsedAge`], row `i` for cell `i`.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn ages(&self) -> Result<Vec<ParsedAge>> {
        let t = &self.columns;
        let cell = t.strings("cell")?;
        let years = t.floats("years")?;
        let min_years = t.floats("min_years")?;
        let max_years = t.floats("max_years")?;
        let precision = t.strings("precision")?;
        let follows = t.booleans("follows_specification")?;
        let refusal = t.strings("refusal")?;
        Ok((0..cell.len())
            .map(|i| ParsedAge {
                cell: cell[i].clone().unwrap_or_default(),
                years: years[i],
                min_years: min_years[i],
                max_years: max_years[i],
                precision: precision[i].clone(),
                follows_specification: follows[i],
                refusal: refusal[i].clone(),
            })
            .collect())
    }
}

// ---- argument assembly for the analysis verbs -----------------------------------------------

fn one_path_args(verb: &str, path: &Path) -> Result<Vec<String>> {
    let path = path_text(path, "file path")?;
    if path.is_empty() {
        return Err(MzLibError::Usage(
            "A file path is required, e.g. 'PXD000070.sdrf.tsv'.".to_owned(),
        ));
    }
    Ok(vec![
        "sdrf".to_owned(),
        verb.to_owned(),
        "--path".to_owned(),
        path.to_owned(),
    ])
}

/// The arguments and stdin of a `*_many` call: the list on stdin, one path per line, and the
/// bulk options as flags. One process reads the whole corpus; nothing here loops over files.
fn bulk_request<P: AsRef<Path>>(
    verb: &str,
    paths: &[P],
    options: &BulkOptions,
) -> Result<(Vec<String>, String)> {
    let lines = bridge::path_lines(paths, "SDRF document")?;
    let args = vec![
        "sdrf".to_owned(),
        verb.to_owned(),
        "--paths-stdin".to_owned(),
        "--threads".to_owned(),
        bridge::threads_arg(options.threads)?,
        "--on-error".to_owned(),
        options.on_error.as_str().to_owned(),
    ];
    Ok((args, lines.join("\n")))
}

/// The stdin of `parse_ages`: one cell per line, blank cells kept, and a trailing newline so a
/// final blank cell survives the trip.
fn age_lines<S: AsRef<str>>(cells: &[S]) -> Result<String> {
    if cells.is_empty() {
        return Err(MzLibError::Usage(
            "parse_ages needs at least one cell.".to_owned(),
        ));
    }
    let mut text = String::new();
    for (index, cell) in cells.iter().enumerate() {
        let cell = cell.as_ref();
        if cell.contains(['\n', '\r']) {
            return Err(MzLibError::Usage(format!(
                "Cell {index} contains a line break, which separates cells: {cell:?}."
            )));
        }
        text.push_str(cell);
        text.push('\n');
    }
    Ok(text)
}

fn call<T: serde::de::DeserializeOwned>(
    args: &[String],
    stdin: Option<&str>,
    timeout: Option<Duration>,
) -> Result<T> {
    let data = bridge::invoke(args, stdin, timeout)?;
    serde_json::from_value(data).map_err(protocol)
}

// ---- the verbs ------------------------------------------------------------------------------

/// Check one SDRF file against the SDRF-Proteomics specification's structural rules, one row per
/// finding.
///
/// Calls mzLib's `SdrfValidator.Validate`. **Structure only** — required and recommended columns,
/// column order, casing, malformed names, ragged rows, integer replicate and fraction columns,
/// reserved-word casing, and the one hard row rule: `source name` + `assay name` +
/// `comment[label]` must be unique. Controlled-vocabulary accessions are **not** resolved. A
/// document of reserved words validates cleanly: use [`assess`] for whether it says anything, and
/// [`lint_labelled`] for whether several agree.
#[doc = include_str!("../docs/reference/sdrf.validate.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let result = mzlib::sdrf::validate("sdrf_cohort.sdrf.tsv")?;
/// assert_eq!((result.is_valid, result.error_count, result.warning_count), (true, 0, 2));
/// for m in result.warnings()? {
///     println!("line {:?}: {} in {:?}", m.line_number, m.rule, m.column_name);
/// }
/// assert_eq!(result.warnings()?[0].line_number, Some(10));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// A file generated from a list of data files fails on structure too:
///
/// ```
/// # mzlib_replay::activate();
/// let skeleton = mzlib::sdrf::validate("sdrf_skeleton.sdrf.tsv")?;
/// let first = &skeleton.errors()?[0];
/// assert!(!skeleton.is_valid);
/// assert_eq!((first.rule.as_str(), first.line_number), ("RequiredColumn", None)); // whole-document
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.validate.see-also.md")]
pub fn validate(path: impl AsRef<Path>) -> Result<SdrfValidation> {
    let args = one_path_args("validate", path.as_ref())?;
    call(&args, None, Some(DEFAULT_TIMEOUT))
}

/// Validate many SDRF files in one bridge call.
///
/// One process for the whole corpus instead of one per file — roughly 120 ms saved per document —
/// with the parallelism inside the bridge. The result does not depend on
/// [`BulkOptions::threads`]: rows are always in input order.
#[doc = include_str!("../docs/reference/sdrf.validate.bulk.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{validate_many, BulkOptions, OnError};
///
/// let batch = validate_many(
///     &["sdrf_skeleton.sdrf.tsv", "sdrf_cohort.sdrf.tsv", "missing.sdrf.tsv"],
///     &BulkOptions { on_error: OnError::Skip, ..Default::default() },
/// )?;
/// assert_eq!((batch.file_count, batch.read_count, batch.valid_count), (3, 2, 1));
/// let valid: Vec<Option<bool>> = batch.files.iter().map(|f| f.is_valid).collect();
/// assert_eq!(valid, [Some(false), Some(true), None]);        // None: not read
/// assert_eq!(batch.failed_files()[0].error.as_ref().unwrap().kind, "usage");
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors
///
/// As [`validate`], per the bulk rules above: an empty list, a blank path or a bad `threads` is
/// refused before anything is spawned.
pub fn validate_many<P: AsRef<Path>>(
    paths: &[P],
    options: &BulkOptions,
) -> Result<SdrfValidationBatch> {
    let (args, stdin) = bulk_request("validate", paths, options)?;
    call(&args, Some(&stdin), options.timeout)
}

/// Find the concepts a set of SDRF documents wrote inconsistently, naming each document by
/// mzLib's `containing-folder/file-stem` default.
///
/// Prefer [`lint_labelled`]: the default label depends on where the files sit. See it for the
/// reference.
///
/// # Errors
///
/// As [`lint_labelled`].
pub fn lint<P: AsRef<Path>>(paths: &[P]) -> Result<SdrfDrift> {
    lint_documents(&PoolInput::Paths(
        paths.iter().map(|p| p.as_ref().to_path_buf()).collect(),
    ))
}

/// Find the concepts a set of SDRF documents wrote inconsistently, one row per finding and
/// spelling, naming each document yourself.
///
/// Validity is a property of one file; comparability is a property of the relationship between
/// files. Every document can pass [`validate`] and the pooled table still be unusable because one
/// file writes `"Homo sapiens"` and another `"homo sapiens"`. This is mzLib's `SdrfDriftLint`
/// check for that. **Label the documents**: the labels are how each row's `documents` names them.
#[doc = include_str!("../docs/reference/sdrf.lint.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let drift = mzlib::sdrf::lint_labelled(&[
///     ("sdrf_cohort.sdrf.tsv", "cohort"),
///     ("sdrf_cohort_partner.sdrf.tsv", "partner"),
/// ])?;
/// assert_eq!(drift.finding_count, 4);
/// for spellings in drift.findings()? {
///     let values: Vec<&str> = spellings.iter().map(|v| v.value.as_str()).collect();
///     println!("{}: {values:?}", spellings[0].kind);
/// }
/// let last = drift.findings()?.pop().unwrap();
/// assert_eq!(last[0].kind, "ValueCaseVariant");
/// assert_eq!((last[0].value.as_str(), last[1].value.as_str()), ("Homo sapiens", "homo sapiens"));
/// assert_eq!(last[1].documents, ["partner"]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.lint.see-also.md")]
pub fn lint_labelled<P: AsRef<Path>, L: AsRef<str>>(documents: &[(P, L)]) -> Result<SdrfDrift> {
    lint_documents(&PoolInput::Labelled(
        documents
            .iter()
            .map(|(p, l)| (p.as_ref().to_path_buf(), l.as_ref().to_owned()))
            .collect(),
    ))
}

fn lint_documents(documents: &PoolInput) -> Result<SdrfDrift> {
    let lines = document_lines(documents)?;
    let args = ["sdrf".to_owned(), "lint".to_owned()];
    call(&args, Some(&lines.join("\n")), Some(DEFAULT_TIMEOUT))
}

/// Decide whether one SDRF file's sample half describes an experimental design — `Informative`,
/// `Partial` or `Skeleton` — with the per-column counts behind the verdict.
///
/// A file generated from a list of data files — reserved words in every sample column, one
/// replicate number, no factor — passes [`validate`], because reserved words are the
/// specification's correct way to say nothing. For grouping results by biology it is the same as
/// having no SDRF at all. mzLib's `SdrfSampleInformativeness.Assess` is the gate for that.
#[doc = include_str!("../docs/reference/sdrf.assess.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let a = mzlib::sdrf::assess("sdrf_cohort.sdrf.tsv")?;
/// assert_eq!(a.verdict, "Informative");
/// assert!(a.factor_value_varies && a.sample_is_described && a.biological_replicate_varies);
/// assert_eq!(a.columns.strings("column_name")?[0].as_deref(), Some("factor value[disease]"));
/// assert_eq!(a.columns.integers("distinct_values")?[0], Some(2));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.assess.see-also.md")]
pub fn assess(path: impl AsRef<Path>) -> Result<SdrfAssessment> {
    let args = one_path_args("assess", path.as_ref())?;
    call(&args, None, Some(DEFAULT_TIMEOUT))
}

/// Assess many SDRF files in one bridge call — the gate for a corpus.
#[doc = include_str!("../docs/reference/sdrf.assess.bulk.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{assess_many, BulkOptions};
///
/// let batch = assess_many(
///     &["sdrf_cohort.sdrf.tsv", "sdrf_skeleton.sdrf.tsv", "PXD000070.sdrf.tsv"],
///     &BulkOptions::default(),
/// )?;
/// let counts = batch.verdict_counts;
/// assert_eq!((counts.informative, counts.partial, counts.skeleton), (1, 1, 1));
/// let usable = batch.paths_with(&["Informative", "Partial"])?;
/// assert_eq!(usable.len(), 2);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors
///
/// As [`assess`], with the bulk rules of [`validate_many`].
pub fn assess_many<P: AsRef<Path>>(
    paths: &[P],
    options: &BulkOptions,
) -> Result<SdrfAssessmentBatch> {
    let (args, stdin) = bulk_request("assess", paths, options)?;
    call(&args, Some(&stdin), options.timeout)
}

/// One row per sample, column and position: what each source name's rows agree on, the columns
/// they disagree about named and withheld, and `characteristics[age]` read into years.
///
/// Calls mzLib's `SdrfSampleBlock.BySourceName`: the sample half of the document (`source name`,
/// every `characteristics[...]`, every `factor value[...]`) per sample, merged over the sample's
/// rows. Where those rows **disagree** — one fraction says `normal` and another `COVID-19` —
/// mzLib names the column and withholds it rather than let the first row win.
#[doc = include_str!("../docs/reference/sdrf.samples.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let s = mzlib::sdrf::samples("sdrf_cohort.sdrf.tsv")?;
/// assert_eq!(s.sample_count, 6);
/// assert_eq!(s.conflicts()?, [("S6".to_owned(), "characteristics[disease]".to_owned())]);
/// for age in s.ages()? {
///     println!("{} {:?} {:?} {:?}", age.source_name, age.value, age.years, age.refusal);
/// }
/// let ages = s.ages()?;
/// assert_eq!((ages[1].years, ages[1].precision.as_deref()), (Some(62.5), Some("Range")));
/// assert_eq!(ages[3].refusal.as_deref(), Some("no_unit"));            // "63": years or days?
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.samples.see-also.md")]
pub fn samples(path: impl AsRef<Path>) -> Result<SdrfSamples> {
    let args = one_path_args("samples", path.as_ref())?;
    call(&args, None, Some(DEFAULT_TIMEOUT))
}

/// Lift the samples of many SDRF files in one bridge call.
#[doc = include_str!("../docs/reference/sdrf.samples.bulk.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{samples_many, BulkOptions};
///
/// let batch = samples_many(
///     &["sdrf_cohort.sdrf.tsv", "sdrf_cohort_partner.sdrf.tsv"],
///     &BulkOptions::default(),
/// )?;
/// let per_file: Vec<Option<u64>> = batch.files.iter().map(|f| f.sample_count).collect();
/// assert_eq!((batch.sample_count, per_file), (8, vec![Some(6), Some(2)]));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors
///
/// As [`samples`], with the bulk rules of [`validate_many`].
pub fn samples_many<P: AsRef<Path>>(
    paths: &[P],
    options: &BulkOptions,
) -> Result<SdrfSamplesBatch> {
    let (args, stdin) = bulk_request("samples", paths, options)?;
    call(&args, Some(&stdin), options.timeout)
}

/// Read `characteristics[age]` cells into years with an honest precision, refusing any cell that
/// would need a guess.
///
/// Calls mzLib's `SdrfAge.TryParse` on each cell. It reads the specification's grammar (`58Y`,
/// `30Y6M`, `16W`), ranges (`40Y-85Y`, `6-8 weeks`), bounds (`>=90Y`, `<1Y`) and unambiguous
/// words (`3 year`, `4 hour`). It **refuses** a bare number — `63` is 11% of real age cells, and 63
/// years and 63 days are both plausible in one study — as well as reserved words and free text.
/// One bridge call for any number of cells, so pass them all at once rather than looping. An empty
/// string is sent as an empty cell and comes back refused as `"empty"`, so the result stays
/// aligned with the input.
#[doc = include_str!("../docs/reference/sdrf.parse-age.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let ages = mzlib::sdrf::parse_ages(&[
///     "58Y", "30Y6M", "40Y-85Y", "40Y-40Y", ">=90Y", "<1Y", "6-8 weeks",
///     "63", "not available", "", "about forty",
/// ])?;
/// assert_eq!((ages.parsed_count, ages.cell_count), (7, 11));
/// let rows = ages.ages()?;
/// assert_eq!((rows[2].years, rows[2].min_years, rows[2].max_years), (Some(62.5), Some(40.0), Some(85.0)));
/// assert_eq!((rows[4].max_years, rows[4].precision.as_deref()), (None, Some("LowerBound")));
/// let refusals: Vec<Option<&str>> = rows[7..].iter().map(|r| r.refusal.as_deref()).collect();
/// assert_eq!(refusals, [Some("no_unit"), Some("reserved_word"), Some("empty"), Some("unreadable")]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.parse-age.see-also.md")]
pub fn parse_ages<S: AsRef<str>>(cells: &[S]) -> Result<ParsedAges> {
    let stdin = age_lines(cells)?;
    let args = ["sdrf".to_owned(), "parse-age".to_owned()];
    call(&args, Some(&stdin), Some(DEFAULT_TIMEOUT))
}

fn protocol(error: serde_json::Error) -> MzLibError {
    MzLibError::Protocol(format!("sdrf payload could not be interpreted: {error}"))
}

#[cfg(test)]
mod tests {
    //! Offline. The fixtures are recorded from the real bridge by pyMzLib and shared verbatim, so
    //! they carry the shape mzLib actually produces — PXD059974's raggedness and PXD000070's eight
    //! repetitions of `comment[modification parameters]` are the point, not an accident.

    use super::*;

    const READ_FIXTURE: &str = include_str!("../tests/fixtures/sdrf_read_PXD000070.json");
    const RAGGED_FIXTURE: &str = include_str!("../tests/fixtures/sdrf_read_ragged.json");
    const POOL_FIXTURE: &str = include_str!("../tests/fixtures/sdrf_pool_two.json");

    fn document() -> SdrfDocument {
        serde_json::from_str(READ_FIXTURE).unwrap()
    }

    fn ragged() -> SdrfDocument {
        serde_json::from_str(RAGGED_FIXTURE).unwrap()
    }

    fn pooled() -> PooledSdrf {
        serde_json::from_str(POOL_FIXTURE).unwrap()
    }

    fn labelled(pairs: &[(&str, &str)]) -> PoolInput {
        PoolInput::Labelled(
            pairs
                .iter()
                .map(|(p, l)| (PathBuf::from(p), (*l).to_owned()))
                .collect(),
        )
    }

    // ---- the shape ----------------------------------------------------------------------------

    #[test]
    fn columns_is_a_vec_because_sdrf_names_repeat() {
        // The whole reason this module is row-major. A map keyed by name would keep one of eight.
        let doc = document();
        assert!(doc.has_repeated_columns());
        assert_eq!(doc.indexes_of("comment[modification parameters]").len(), 8);
    }

    #[test]
    fn all_returns_every_occurrence_where_value_returns_the_first() {
        let doc = document();
        let every = &doc.all("comment[modification parameters]")[0];
        let first = doc.value("comment[modification parameters]")[0].unwrap();

        assert_eq!(every.len(), 8);
        assert_eq!(every[0], first);
        assert!(first.starts_with("NT=Carbamidomethyl"));
    }

    #[test]
    fn a_cv_cell_crosses_intact_rather_than_being_split_on_its_semicolons() {
        // The defect this module exists to close: read_records joins an SdrfRow's cells with ";".
        let doc = document();
        assert_eq!(
            doc.value("comment[modification parameters]")[0],
            Some("NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed")
        );
    }

    #[test]
    fn an_absent_column_is_none_and_a_reserved_word_is_itself() {
        // None means "no such column"; "not applicable" is a real answer a curator wrote.
        let doc = document();
        assert_eq!(
            doc.value("characteristics[nonesuch]"),
            vec![None; doc.rows.len()]
        );
        assert_eq!(doc.index_of("characteristics[nonesuch]"), None);
        assert_eq!(
            doc.value("characteristics[disease]")[0],
            Some("not applicable")
        );
    }

    #[test]
    fn a_short_row_reports_none_rather_than_shifting_later_columns() {
        let doc = ragged();
        assert_eq!(doc.ragged_row_count(), 17);
        let short = doc
            .rows
            .iter()
            .position(|row| row.len() < doc.columns.len())
            .unwrap();
        let full = doc
            .rows
            .iter()
            .position(|row| row.len() == doc.columns.len())
            .unwrap();

        let last = doc.columns.last().unwrap();
        assert!(doc.value(last)[full].is_some());
        assert!(doc.value(last)[short].is_none());
        // A position both rows reach still lines up, which is what "no shifting" means.
        assert!(doc.value("characteristics[organism]")[short].is_some());
    }

    #[test]
    fn caveats_name_the_raggedness_rather_than_leaving_it_to_be_discovered() {
        assert!(ragged().caveats.iter().any(|c| c.contains("SHORT")));
    }

    #[test]
    fn truncated_says_rows_were_left_behind() {
        // Recorded with limit=4 against a 24-row merge.
        let pooled = pooled();
        assert!(pooled.document.truncated);
        assert_eq!(pooled.document.returned_count, 4);
        assert_eq!(pooled.document.row_count, 24);
    }

    #[test]
    fn a_complete_read_is_not_marked_truncated() {
        let doc = ragged();
        assert!(!doc.truncated);
        assert_eq!(doc.returned_count, doc.row_count);
    }

    #[test]
    fn records_is_offered_but_lossy_when_names_repeat() {
        let doc = document();
        assert!(doc.has_repeated_columns());
        assert!(doc.records()[0].len() < doc.columns.len());
    }

    #[test]
    fn a_null_list_on_the_wire_is_empty_rather_than_an_error() {
        let doc: SdrfDocument = serde_json::from_str(
            r#"{"path": "x", "column_names": null, "rows": null, "caveats": null}"#,
        )
        .unwrap();
        assert!(doc.columns.is_empty() && doc.rows.is_empty() && doc.caveats.is_empty());
    }

    // ---- pool ---------------------------------------------------------------------------------

    #[test]
    fn pool_carries_provenance_and_the_labels_it_was_given() {
        let pooled = pooled();
        assert_eq!(pooled.document_count, 2);
        assert_eq!(pooled.labels, ["malaria", "colon"]);
        assert!(pooled.written.is_none());
        assert!(pooled
            .document
            .columns
            .iter()
            .any(|c| c == SOURCE_DOCUMENT_COLUMN));
        assert!(pooled
            .source_documents()
            .iter()
            .all(|s| matches!(s, Some("malaria" | "colon"))));
    }

    #[test]
    fn pool_renders_one_tab_separated_stdin_line_per_document() {
        let (args, stdin) = pool_request(
            &labelled(&[("a.sdrf.tsv", "malaria"), ("b.sdrf.tsv", "colon")]),
            &PoolOptions::default(),
        )
        .unwrap();
        assert_eq!(stdin, "a.sdrf.tsv\tmalaria\nb.sdrf.tsv\tcolon");
        assert_eq!(args, ["sdrf", "pool"]);
    }

    #[test]
    fn pool_without_labels_sends_bare_paths() {
        let (_, stdin) = pool_request(
            &PoolInput::Paths(vec!["a.sdrf.tsv".into(), "b.sdrf.tsv".into()]),
            &PoolOptions::default(),
        )
        .unwrap();
        assert_eq!(stdin, "a.sdrf.tsv\nb.sdrf.tsv");
    }

    #[test]
    fn pool_passes_the_fallback_warning_through() {
        // The recorded fixture was pooled WITH labels, so it carries no warning; the bridge adds
        // one when labels are absent, and it must reach the caller.
        assert!(!pooled()
            .document
            .caveats
            .iter()
            .any(|c| c.contains("reproducible")));

        let mut payload: serde_json::Value = serde_json::from_str(POOL_FIXTURE).unwrap();
        payload["caveats"]
            .as_array_mut()
            .unwrap()
            .push("... not reproducible on another machine ...".into());
        let unlabelled: PooledSdrf = serde_json::from_value(payload).unwrap();
        assert!(unlabelled
            .document
            .caveats
            .iter()
            .any(|c| c.contains("reproducible")));
    }

    #[test]
    fn pool_writes_the_whole_document_and_reports_it() {
        let (args, _) = pool_request(
            &PoolInput::Paths(vec!["a.sdrf.tsv".into()]),
            &PoolOptions {
                read: ReadOptions {
                    limit: Some(4),
                    ..ReadOptions::default()
                },
                out: Some("merged.sdrf.tsv".to_owned()),
            },
        )
        .unwrap();
        assert_eq!(
            args,
            ["sdrf", "pool", "--limit", "4", "--out", "merged.sdrf.tsv"]
        );

        let mut payload: serde_json::Value = serde_json::from_str(POOL_FIXTURE).unwrap();
        payload["written"] = serde_json::json!({"path": "merged.sdrf.tsv", "row_count": 24});
        let written: PooledSdrf = serde_json::from_value(payload).unwrap();
        assert_eq!(
            written.written,
            Some(WrittenSdrf {
                path: "merged.sdrf.tsv".to_owned(),
                row_count: 24
            })
        );
    }

    // ---- argument validation, before anything is spawned --------------------------------------

    #[test]
    fn an_empty_selection_is_refused() {
        let error = pool_request(&PoolInput::Paths(vec![]), &PoolOptions::default()).unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(error.to_string().contains("At least one"), "{error}");
    }

    #[test]
    fn a_blank_label_is_refused_rather_than_silently_defaulted() {
        let error = pool_request(
            &labelled(&[("a.sdrf.tsv", "malaria"), ("b.sdrf.tsv", "   ")]),
            &PoolOptions::default(),
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(
            error.to_string().contains("label for 'b.sdrf.tsv'"),
            "{error}"
        );
    }

    #[test]
    fn a_tab_in_a_label_is_refused_because_it_is_the_field_separator() {
        let error = pool_request(
            &labelled(&[("a.sdrf.tsv", "mal\taria")]),
            &PoolOptions::default(),
        )
        .unwrap_err();
        assert!(error.to_string().contains("tab or newline"), "{error}");
    }

    #[test]
    fn a_blank_path_is_refused() {
        for path in ["", "   "] {
            let error = read_args(Path::new(path), &ReadOptions::default()).unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)));
            assert!(
                error.to_string().contains("file path is required"),
                "{error}"
            );
        }
        let error =
            pool_request(&PoolInput::Paths(vec![" ".into()]), &PoolOptions::default()).unwrap_err();
        assert!(error.to_string().contains("may not be blank"), "{error}");
    }

    #[test]
    fn an_empty_out_is_refused() {
        let error = pool_request(
            &PoolInput::Paths(vec!["a.sdrf.tsv".into()]),
            &PoolOptions {
                out: Some("  ".to_owned()),
                ..PoolOptions::default()
            },
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    #[test]
    fn the_window_options_reach_the_bridge() {
        let args = read_args(
            Path::new("a.sdrf.tsv"),
            &ReadOptions {
                limit: Some(2),
                offset: 1,
                ..ReadOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "sdrf",
                "read",
                "--path",
                "a.sdrf.tsv",
                "--limit",
                "2",
                "--offset",
                "1"
            ]
        );
    }

    #[test]
    fn a_zero_limit_is_sent_and_a_zero_offset_is_not() {
        // The bridge accepts --limit 0 (a header-only answer); a default offset sent explicitly is
        // a default the bridge could later disagree with.
        let args = read_args(
            Path::new("a.sdrf.tsv"),
            &ReadOptions {
                limit: Some(0),
                ..ReadOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            args,
            ["sdrf", "read", "--path", "a.sdrf.tsv", "--limit", "0"]
        );
    }

    #[test]
    fn the_default_timeout_matches_pymzlib() {
        assert_eq!(
            ReadOptions::default().timeout,
            Some(Duration::from_secs(60))
        );
    }

    // ---- validate, lint, assess, samples, parse_ages ------------------------------------------

    fn fixture<T: serde::de::DeserializeOwned>(text: &str) -> T {
        serde_json::from_str(text).unwrap()
    }

    #[test]
    fn a_validation_reads_into_typed_findings() {
        let cohort: SdrfValidation =
            fixture(include_str!("../tests/fixtures/sdrf_validate_cohort.json"));
        assert!(cohort.is_valid);
        assert_eq!((cohort.error_count, cohort.warning_count), (0, 2));
        assert_eq!(cohort.columns.rows() as u64, cohort.message_count);
        let warnings = cohort.warnings().unwrap();
        assert_eq!(warnings.len(), 2);
        assert!(warnings.iter().all(|m| m.rule == "ReservedWordCase"));
        // line_number is row_index + 2: the header is line 1.
        for m in &warnings {
            assert_eq!(m.line_number, m.row_index.map(|r| r + 2));
        }
        assert!(cohort.errors().unwrap().is_empty());
    }

    #[test]
    fn a_document_level_finding_has_no_row_and_no_line() {
        let skeleton: SdrfValidation = fixture(include_str!(
            "../tests/fixtures/sdrf_validate_skeleton.json"
        ));
        assert!(!skeleton.is_valid);
        assert_eq!(
            skeleton.errors().unwrap().len() as u64,
            skeleton.error_count
        );
        let first = &skeleton.errors().unwrap()[0];
        assert_eq!(first.rule, "RequiredColumn");
        assert_eq!((first.row_index, first.line_number), (None, None));
    }

    #[test]
    fn a_skipped_document_is_named_and_counts_nothing() {
        let batch: SdrfValidationBatch =
            fixture(include_str!("../tests/fixtures/sdrf_validate_bulk.json"));
        assert_eq!(
            (batch.file_count, batch.read_count, batch.failed_count),
            (3, 2, 1)
        );
        assert_eq!(batch.record_count, batch.message_count);
        let missing = &batch.files[2];
        assert_eq!(missing.is_valid, None); // not read, which is not "invalid"
        assert_eq!(missing.error.as_ref().unwrap().kind, "usage");
        assert_eq!(batch.failed_files().len(), 1);
        assert_eq!(batch.columns.names()[..2], ["source_index", "source_path"]);
        assert_eq!(batch.messages().unwrap().len() as u64, batch.message_count);
    }

    #[test]
    fn drift_groups_by_finding_and_keeps_the_label_lists() {
        let drift: SdrfDrift = fixture(include_str!("../tests/fixtures/sdrf_lint_cohort.json"));
        assert_eq!(drift.labels, ["cohort", "partner"]);
        let findings = drift.findings().unwrap();
        assert_eq!(findings.len() as u64, drift.finding_count);
        assert!(findings.iter().all(|f| f[0].variant_rank == 0));
        let names = &findings[2];
        assert_eq!(names[0].kind, "ColumnNameVariant");
        assert_eq!(names[0].column_name, None); // a finding about a name, not a column
        assert_eq!(findings[0][0].documents, ["partner"]);
    }

    #[test]
    fn an_assessment_reads_its_verdict_and_integer_fill_rates_as_floats() {
        let a: SdrfAssessment = fixture(include_str!("../tests/fixtures/sdrf_assess_cohort.json"));
        assert_eq!(a.verdict, "Informative");
        assert!(VERDICTS.contains(&a.verdict.as_str()));
        // fill_rate crosses as 1, not 1.0: a float accessor must still read it.
        let rates = a.columns.floats("fill_rate").unwrap();
        assert!(rates
            .iter()
            .all(|r| r.is_some_and(|r| (0.0..=1.0).contains(&r))));
    }

    #[test]
    fn paths_with_refuses_a_verdict_that_does_not_exist() {
        let batch: SdrfAssessmentBatch =
            fixture(include_str!("../tests/fixtures/sdrf_assess_bulk.json"));
        assert_eq!(
            batch.paths_with(&["Skeleton"]).unwrap(),
            ["sdrf_skeleton.sdrf.tsv"]
        );
        let error = batch.paths_with(&["Informativ"]).unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert_eq!(
            batch.verdict_counts,
            VerdictCounts {
                informative: 1,
                partial: 1,
                skeleton: 1
            }
        );
    }

    #[test]
    fn samples_withhold_a_conflict_and_parse_ages() {
        let s: SdrfSamples = fixture(include_str!("../tests/fixtures/sdrf_samples_cohort.json"));
        assert_eq!(s.conflict_count as usize, s.conflicts().unwrap().len());
        let ages = s.ages().unwrap();
        assert_eq!(ages.len(), 6);
        assert_eq!(ages[0].years, Some(58.0));
        // >=90Y has no upper bound: None, told apart from a refusal by its precision.
        assert_eq!(
            (ages[2].max_years, ages[2].precision.as_deref()),
            (None, Some("LowerBound"))
        );
        assert_eq!(ages[4].refusal.as_deref(), Some("reserved_word"));
        let batch: SdrfSamplesBatch =
            fixture(include_str!("../tests/fixtures/sdrf_samples_bulk.json"));
        assert_eq!(batch.sample_count, 8);
        assert_eq!(batch.columns.rows() as u64, batch.record_count);
    }

    #[test]
    fn parsed_ages_stay_aligned_with_the_cells() {
        let ages: ParsedAges = fixture(include_str!("../tests/fixtures/sdrf_parse_age.json"));
        let rows = ages.ages().unwrap();
        assert_eq!(rows.len() as u64, ages.cell_count);
        assert_eq!(rows[9].cell, "");
        assert_eq!(rows[9].refusal.as_deref(), Some("empty"));
        for row in &rows {
            assert_eq!(
                row.refusal.is_none(),
                row.precision.is_some(),
                "{}",
                row.cell
            );
            if let Some(refusal) = &row.refusal {
                assert!(AGE_REFUSALS.contains(&refusal.as_str()));
            }
        }
    }

    #[test]
    fn a_bulk_call_sends_its_options_and_one_path_per_line() {
        let (args, stdin) = bulk_request(
            "validate",
            &["a.sdrf.tsv", "b.sdrf.tsv"],
            &BulkOptions {
                threads: -1,
                on_error: OnError::Skip,
                timeout: None,
            },
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "sdrf",
                "validate",
                "--paths-stdin",
                "--threads",
                "-1",
                "--on-error",
                "skip"
            ]
        );
        assert_eq!(stdin, "a.sdrf.tsv\nb.sdrf.tsv");
        assert_eq!(BulkOptions::default().threads, 1);
        assert_eq!(BulkOptions::default().on_error, OnError::Fail);
    }

    #[test]
    fn bulk_calls_refuse_bad_input_before_anything_is_spawned() {
        let options = BulkOptions::default();
        assert!(bulk_request::<&str>("assess", &[], &options).is_err());
        assert!(bulk_request("assess", &["a", "  "], &options).is_err());
        assert!(bulk_request("assess", &["a\nb"], &options).is_err());
        let zero = BulkOptions {
            threads: 0,
            ..BulkOptions::default()
        };
        assert!(matches!(
            bulk_request("assess", &["a"], &zero),
            Err(MzLibError::Usage(_))
        ));
        assert!(one_path_args("samples", Path::new(" ")).is_err());
        assert_eq!(
            one_path_args("samples", Path::new("x.sdrf.tsv")).unwrap(),
            ["sdrf", "samples", "--path", "x.sdrf.tsv"]
        );
    }

    #[test]
    fn age_cells_keep_blanks_and_end_with_a_newline() {
        assert_eq!(age_lines(&["58Y", "", "63"]).unwrap(), "58Y\n\n63\n");
        // A final blank cell must survive: without the trailing newline it would vanish.
        assert_eq!(age_lines(&["58Y", ""]).unwrap(), "58Y\n\n");
        assert!(age_lines::<&str>(&[]).is_err());
        assert!(age_lines(&["58Y\n63"]).is_err());
    }

    #[test]
    fn lint_shares_pools_document_lines() {
        let lines =
            document_lines(&labelled(&[("a.sdrf.tsv", "one"), ("b.sdrf.tsv", "two")])).unwrap();
        assert_eq!(lines, ["a.sdrf.tsv\tone", "b.sdrf.tsv\ttwo"]);
        assert!(document_lines(&labelled(&[("a.sdrf.tsv", " ")])).is_err());
    }
}
