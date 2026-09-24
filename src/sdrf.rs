//! SDRF-Proteomics experimental-design files: read one, or pool several into one table.
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
//! **What this module does not do yet is validate.** mzLib models SDRF's structural rules in
//! `SdrfValidator` and its vocabulary-drift rules in `SdrfDriftLint`. Both have been public since
//! mzLib #1207, which the pinned bridge includes, but the bridge does not expose them yet. When it
//! does, they will be projected here from mzLib rather than reimplemented. Until then, this module
//! reads, pools and reports honestly, and makes no claim about whether a document is *correct*.
//!
//! Ported from pyMzLib's `pymzlib.sdrf`, which decided the verbs, the wire fields and the caveats.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};

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
}
