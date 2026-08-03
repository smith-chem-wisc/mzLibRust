//! PRIDE Archive access, backed by mzLib's `PrideArchiveClient`.
//!
//! The [PRIDE Archive](https://www.ebi.ac.uk/pride/archive/) is EBI's public proteomics data
//! repository. This module lets a Rust user list what is in a project and pull files down, using
//! the same paging, URL-resolution, and safe-download logic that mzLib uses in C#.
//!
//! ```no_run
//! # fn main() -> Result<(), mzlib::MzLibError> {
//! let files = mzlib::pride::list_files("PXD000001")?;
//! println!("{} files", files.len());
//!
//! // Filter however you like — the full expressiveness of Rust — then fetch exactly that.
//! let small: Vec<_> = files
//!     .iter()
//!     .filter(|f| f.size_mb() < 5.0 && f.downloadable())
//!     .cloned()
//!     .collect();
//! mzlib::pride::download_files(&small, "downloads", &Default::default())?;
//! # Ok(())
//! # }
//! ```

use std::path::{Path, PathBuf};
use std::time::Duration;

use chrono::{DateTime, FixedOffset};
use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};

/// One published location of a file, as a controlled-vocabulary term.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct CvParam {
    /// The CV accession, e.g. `"PRIDE:0000469"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub accession: String,
    /// The term's name, e.g. `"FTP Protocol"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub name: String,
    /// The term's value — for a location term, the URL.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub value: String,
}

/// One file belonging to a PRIDE Archive project.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct PrideFile {
    /// The file's name, e.g. `"run1.raw"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_name: String,
    /// Size in bytes as reported by PRIDE.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_size_bytes: u64,
    /// The repository's checksum, or `""` if it provides none.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub checksum: String,
    /// The file category, e.g. `"RAW"`, `"PEAK"`, `"SEARCH"`, `"OTHER"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub category: String,
    /// The category's controlled-vocabulary accession.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub category_accession: String,
    /// A direct HTTPS download URL, or `None` when the file is only reachable by a protocol that
    /// cannot be fetched over HTTPS (Aspera-only files).
    #[serde(default)]
    pub https_url: Option<String>,
    /// Every published location, for callers that want the raw controlled-vocabulary terms.
    #[serde(default)]
    pub locations: Vec<CvParam>,
    /// When the file was submitted, if the repository said.
    #[serde(default, deserialize_with = "lenient_timestamp")]
    pub submission_date: Option<DateTime<FixedOffset>>,
    /// When the file was published, if the repository said.
    #[serde(default, deserialize_with = "lenient_timestamp")]
    pub publication_date: Option<DateTime<FixedOffset>>,
    /// When the file was last updated, if the repository said.
    #[serde(default, deserialize_with = "lenient_timestamp")]
    pub updated_date: Option<DateTime<FixedOffset>>,
    /// Which project this file came from.
    ///
    /// Not on the wire — stamped on by [`list_files`], so [`download_files`] can tell which project
    /// to fetch from without the caller having to carry the accession alongside the list.
    #[serde(skip)]
    pub project_accession: String,
}

impl PrideFile {
    /// The file size in megabytes, for the common case of eyeballing a manifest.
    #[must_use]
    pub fn size_mb(&self) -> f64 {
        self.file_size_bytes as f64 / 1_000_000.0
    }

    /// The file's lowercase extension including the dot, e.g. `".raw"`. Empty if none.
    ///
    /// Note that a compressed file such as `x.mgf.gz` has the extension `".gz"` — which is what
    /// trips people up when a filter matches nothing.
    #[must_use]
    pub fn extension(&self) -> String {
        Path::new(&self.file_name)
            .extension()
            .map(|ext| format!(".{}", ext.to_string_lossy().to_lowercase()))
            .unwrap_or_default()
    }

    /// Whether this file can be fetched by [`download`] — i.e. it has an HTTPS location.
    #[must_use]
    pub fn downloadable(&self) -> bool {
        self.https_url.is_some()
    }
}

/// One file found by walking a PRIDE project's FTP directory tree — the COMPLETE listing.
///
/// PRIDE's REST manifest ([`list_files`]) is knowingly incomplete: for PXD000001 it returns 8
/// files while the FTP tree holds 13, omitting the two largest. When completeness or a true project
/// size matters, [`list_ftp_files`] walks the directory (subdirectories included) and returns
/// everything (mzLib #1121). The trade-off is the size: [`PrideFtpFile::approximate_size_bytes`]
/// is PRIDE's rounded index value, not the exact transfer size.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct PrideFtpFile {
    /// Path relative to the project's FTP root, e.g. `"run1.raw"` or, for a file in a subdirectory,
    /// `"generated/summary.mztab"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub relative_path: String,
    /// The bare file name — the last segment of [`PrideFtpFile::relative_path`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_name: String,
    /// The HTTPS URL the file can be downloaded from.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub url: String,
    /// PRIDE's rounded index size in bytes — good for a project-size estimate, but not exact. For
    /// the precise transfer size of one file, issue an HTTP HEAD against [`PrideFtpFile::url`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub approximate_size_bytes: u64,
}

impl PrideFtpFile {
    /// The approximate size in megabytes, for eyeballing a project's footprint.
    #[must_use]
    pub fn approximate_size_mb(&self) -> f64 {
        self.approximate_size_bytes as f64 / 1_000_000.0
    }

    /// The file's lowercase extension including the dot, e.g. `".raw"`. Empty if none.
    #[must_use]
    pub fn extension(&self) -> String {
        Path::new(&self.file_name)
            .extension()
            .map(|ext| format!(".{}", ext.to_string_lossy().to_lowercase()))
            .unwrap_or_default()
    }
}

/// How a manifest is fetched. Defaults match pyMzLib's: 100 files per API call, 300 s.
#[derive(Debug, Clone)]
pub struct ListOptions {
    /// How many files to request per underlying API call. Only affects how the manifest is
    /// fetched, never what you get back.
    pub page_size: u32,
    /// Seconds to allow for the whole fetch.
    pub timeout: Option<Duration>,
}

impl Default for ListOptions {
    fn default() -> Self {
        Self {
            page_size: 100,
            timeout: Some(Duration::from_secs(300)),
        }
    }
}

/// How a download is filtered and how long it may take.
#[derive(Debug, Clone, Default)]
pub struct DownloadOptions {
    /// Keep only files of this category, e.g. `"RAW"`. `None` keeps all.
    pub category: Option<String>,
    /// Keep only files with these extensions, e.g. `[".raw", ".mzML"]`. Empty keeps all.
    /// Combined with `category` as AND.
    ///
    /// **A compressed file's extension is `.gz`, not what it is compressed from.** This is the
    /// most common way to get an empty result: PXD000001's peak list is
    /// `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz`, so `[".mgf"]` matches **nothing** — while
    /// `[".gz"]` over-matches to three unrelated files. To select one compressed type, combine
    /// [`Self::category`] with the extension (`category: "PEAK"`, `extensions: [".gz"]`), or skip
    /// these filters entirely and pass the files you want to [`download_files`].
    ///
    /// A filter that matched nothing is an error here, not an empty success — see
    /// [`PrideFile::extension`].
    pub extensions: Vec<String>,
    /// When `false`, a file already present at the destination is left alone and not re-fetched —
    /// a cheap resume for a large project.
    pub overwrite: Option<bool>,
    /// Seconds to allow. `None` waits as long as it takes, which is usually what you want for
    /// multi-gigabyte projects.
    pub timeout: Option<Duration>,
}

impl DownloadOptions {
    /// Whether existing files should be replaced. Defaults to `true`, matching pyMzLib.
    fn overwrite_or_default(&self) -> bool {
        self.overwrite.unwrap_or(true)
    }
}

// ------------------------------------------------------------------ validation

/// Validate and canonicalise an accession, failing loudly rather than returning nothing.
///
/// Accessions are upper-cased, because PRIDE's API is case-sensitive on the accession while the
/// category matching is case-insensitive — two rules pointing opposite ways is a trap, and this is
/// the one that can be fixed without surprising anybody.
///
/// The grammar is a short letter prefix and a run of digits, hand-rolled rather than pulled from
/// `regex`: `^[A-Z]{2,4}[0-9]{4,}$`.
fn normalise_accession(accession: &str) -> Result<String> {
    let candidate = accession.trim().to_ascii_uppercase();
    if candidate.is_empty() {
        return Err(MzLibError::Usage(
            "A PRIDE project accession is required, e.g. 'PXD000001'.".to_owned(),
        ));
    }

    let letters = candidate
        .chars()
        .take_while(char::is_ascii_uppercase)
        .count();
    let digits = candidate[letters..].len();
    let all_digits = candidate[letters..].chars().all(|c| c.is_ascii_digit());

    if !(2..=4).contains(&letters) || digits < 4 || !all_digits {
        return Err(MzLibError::Usage(format!(
            "'{accession}' is not a valid repository accession. Expected a short letter prefix \
             followed by digits, e.g. 'PXD000001'."
        )));
    }

    Ok(candidate)
}

/// Reject a blank destination instead of quietly writing into the current directory.
///
/// An empty path is the current directory, so `dest = config.outdir.unwrap_or_default()` would
/// spray a multi-gigabyte project across the working directory.
fn normalise_destination(destination: &Path) -> Result<&Path> {
    if destination.as_os_str().is_empty() || destination.to_string_lossy().trim().is_empty() {
        return Err(MzLibError::Usage(
            "A destination directory is required; got an empty path.".to_owned(),
        ));
    }
    Ok(destination)
}

/// Accept a list of extensions, refusing one that normalises to nothing.
///
/// A caller who asked for extensions and whose list normalises to nothing would have had `--ext`
/// omitted entirely, which the bridge reads as "no filter" and downloads the whole project. Asking
/// for a filter and getting everything is never right.
fn normalise_extensions(extensions: &[String]) -> Result<Vec<String>> {
    for value in extensions {
        if value.contains(',') {
            return Err(MzLibError::Usage(format!(
                "An extension may not contain a comma; got '{value}'. Pass separate list items."
            )));
        }
    }

    let kept: Vec<String> = extensions
        .iter()
        .map(|value| value.trim().to_owned())
        .filter(|value| !value.is_empty())
        .collect();

    if !extensions.is_empty() && kept.is_empty() {
        return Err(MzLibError::Usage(format!(
            "extensions was given but names no extensions; got {extensions:?}. Omit it to \
             download every file type."
        )));
    }
    Ok(kept)
}

/// Refuse a value that would be parsed as another option by the bridge.
///
/// The bridge's parser treats `--a --b` as two flags, so a value beginning with `-` silently
/// discards the option it belonged to — and can smuggle in a flag the caller never intended.
fn reject_flag_like(name: &str, value: &str) -> Result<String> {
    if value.starts_with('-') {
        return Err(MzLibError::Usage(format!(
            "{name} may not begin with '-'; got '{value}'. That would be read as another option."
        )));
    }
    Ok(value.to_owned())
}

/// Convert an ISO-8601 timestamp from the bridge, treating anything unreadable as absent.
///
/// A timestamp that will not parse is not worth failing a whole manifest over — the caller asked
/// for a file list, not a date. It becomes `None`, exactly as pyMzLib's `_parse_timestamp` does.
fn lenient_timestamp<'de, D>(
    deserializer: D,
) -> std::result::Result<Option<DateTime<FixedOffset>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let raw: Option<String> = Option::deserialize(deserializer)?;
    Ok(raw.as_deref().and_then(parse_timestamp))
}

/// Parse one ISO-8601 timestamp, with or without an offset. `None` if it will not parse.
fn parse_timestamp(value: &str) -> Option<DateTime<FixedOffset>> {
    if value.trim().is_empty() {
        return None;
    }
    if let Ok(parsed) = DateTime::parse_from_rfc3339(value) {
        return Some(parsed);
    }
    // System.Text.Json writes a DateTime with no offset when the source had none.
    for format in ["%Y-%m-%dT%H:%M:%S%.f", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"] {
        if let Ok(naive) = chrono::NaiveDateTime::parse_from_str(value, format) {
            return Some(naive.and_utc().fixed_offset());
        }
        if let Ok(date) = chrono::NaiveDate::parse_from_str(value, format) {
            return Some(date.and_hms_opt(0, 0, 0)?.and_utc().fixed_offset());
        }
    }
    None
}

// ------------------------------------------------------------------ argument assembly

/// The argv for `pride files`.
fn build_list_args(accession: &str, page_size: u32) -> Result<Vec<String>> {
    let canonical = normalise_accession(accession)?;
    if page_size == 0 {
        return Err(MzLibError::Usage(format!(
            "page_size must be positive; got {page_size}."
        )));
    }
    if page_size > i32::MAX as u32 {
        return Err(MzLibError::Usage(format!(
            "page_size is larger than the API allows; got {page_size}."
        )));
    }

    Ok(vec![
        "pride".to_owned(),
        "files".to_owned(),
        "--accession".to_owned(),
        canonical,
        "--page-size".to_owned(),
        page_size.to_string(),
    ])
}

/// The argv for `pride download`.
fn build_download_args(
    accession: &str,
    destination: &Path,
    options: &DownloadOptions,
) -> Result<Vec<String>> {
    let canonical = normalise_accession(accession)?;
    let target = normalise_destination(destination)?;
    let wanted = normalise_extensions(&options.extensions)?;

    let mut args = vec![
        "pride".to_owned(),
        "download".to_owned(),
        "--accession".to_owned(),
        canonical,
        "--dest".to_owned(),
        target.to_string_lossy().into_owned(),
    ];

    if let Some(category) = &options.category {
        if category.trim().is_empty() {
            return Err(MzLibError::Usage(
                "category is empty. Omit it to download every category, rather than passing a \
                 blank value — a filter that selects nothing must not silently select everything."
                    .to_owned(),
            ));
        }
        args.push("--category".to_owned());
        args.push(reject_flag_like("category", category.trim())?);
    }

    if !wanted.is_empty() {
        args.push("--ext".to_owned());
        args.push(reject_flag_like("extensions", &wanted.join(","))?);
    }

    if !options.overwrite_or_default() {
        args.push("--no-overwrite".to_owned());
    }

    Ok(args)
}

/// The argv and stdin payload for `pride download --names-from-stdin`.
fn build_download_files_args(
    files: &[PrideFile],
    destination: &Path,
    overwrite: bool,
) -> Result<(Vec<String>, String)> {
    let target = normalise_destination(destination)?;

    if files.is_empty() {
        return Err(MzLibError::Usage(
            "No files selected. An empty selection is almost always a filter that did not match \
             what you expected, so mzLibRust refuses it rather than reporting success."
                .to_owned(),
        ));
    }

    let unreachable: Vec<&str> = files
        .iter()
        .filter(|file| !file.downloadable())
        .map(|file| file.file_name.as_str())
        .collect();
    if !unreachable.is_empty() {
        return Err(MzLibError::Usage(format!(
            "{} of {} selected files have no HTTPS location and cannot be downloaded (e.g. '{}'). \
             Filter on `downloadable()` first.",
            unreachable.len(),
            files.len(),
            unreachable[0]
        )));
    }

    let mut accessions: Vec<&str> = files
        .iter()
        .map(|file| file.project_accession.as_str())
        .filter(|accession| !accession.is_empty())
        .collect();
    accessions.sort_unstable();
    accessions.dedup();

    if accessions.len() > 1 {
        return Err(MzLibError::Usage(format!(
            "All files must come from one project; got {accessions:?}."
        )));
    }
    let Some(accession) = accessions.first() else {
        return Err(MzLibError::Usage(
            "These PrideFile values carry no project accession, so mzLibRust cannot tell which \
             project to fetch from. Obtain them from list_files()."
                .to_owned(),
        ));
    };

    // The selection travels on stdin rather than argv: a few thousand names would blow the ~32 KB
    // command-line ceiling. The framing is newline-delimited, which is *almost* general — a POSIX
    // file name may legally contain a newline, so such a name would split into two and silently
    // select the wrong files. PRIDE has never published one, but "never seen it" is not a
    // contract, so it is refused explicitly rather than mis-parsed quietly.
    if let Some(offender) = files
        .iter()
        .find(|file| file.file_name.contains('\n') || file.file_name.contains('\r'))
    {
        return Err(MzLibError::Usage(format!(
            "Cannot select '{}': the file name contains a line break, which the selection format \
             cannot represent. Please open an issue — this is a limitation worth fixing properly \
             if a real repository ever publishes such a name.",
            offender.file_name.escape_debug()
        )));
    }

    let mut args = vec![
        "pride".to_owned(),
        "download".to_owned(),
        "--accession".to_owned(),
        (*accession).to_owned(),
        "--dest".to_owned(),
        target.to_string_lossy().into_owned(),
        "--names-from-stdin".to_owned(),
    ];
    if !overwrite {
        args.push("--no-overwrite".to_owned());
    }

    let stdin = files
        .iter()
        .map(|file| file.file_name.as_str())
        .collect::<Vec<_>>()
        .join("\n");

    Ok((args, stdin))
}

// ------------------------------------------------------------------ parsing

/// Turn the `pride files` payload into typed files, refusing an empty manifest.
fn parse_manifest(data: &serde_json::Value, accession: &str) -> Result<Vec<PrideFile>> {
    let raw = data
        .get("files")
        .cloned()
        .unwrap_or(serde_json::Value::Null);
    let mut files: Vec<PrideFile> = if raw.is_null() {
        Vec::new()
    } else {
        serde_json::from_value(raw).map_err(|error| {
            MzLibError::Protocol(format!("PRIDE manifest could not be read: {error}"))
        })?
    };

    for file in &mut files {
        file.project_accession = accession.to_owned();
    }

    if files.is_empty() {
        return Err(MzLibError::ProjectNotFound(format!(
            "PRIDE returned no files for '{accession}'. Either the accession does not exist \
             (check for a typo) or the project is private. PRIDE does not distinguish the two, so \
             neither can mzLibRust."
        )));
    }

    Ok(files)
}

/// Parse the `files` array of a `pride ftp-files` payload. Simpler than [`parse_manifest`]: there is
/// no per-file project accession to stamp, and the empty/not-found decision is made by the caller.
fn parse_ftp_files(data: &serde_json::Value) -> Result<Vec<PrideFtpFile>> {
    let raw = data
        .get("files")
        .cloned()
        .unwrap_or(serde_json::Value::Null);
    if raw.is_null() {
        return Ok(Vec::new());
    }
    serde_json::from_value(raw).map_err(|error| {
        MzLibError::Protocol(format!("PRIDE FTP listing could not be read: {error}"))
    })
}

/// Read the written paths out of a `pride download` payload.
fn parse_paths(data: &serde_json::Value) -> Vec<PathBuf> {
    data.get("paths")
        .and_then(serde_json::Value::as_array)
        .map(|paths| {
            paths
                .iter()
                .filter_map(serde_json::Value::as_str)
                .map(PathBuf::from)
                .collect()
        })
        .unwrap_or_default()
}

/// Refuse to report success for a filter that matched nothing.
///
/// A filter that matched nothing is nearly always a filter that does not mean what its author
/// thought, and reporting success with an empty list lets a batch script carry on as though the
/// work had been done.
fn check_filter_matched(
    written: &[PathBuf],
    accession: &str,
    options: &DownloadOptions,
) -> Result<()> {
    let filtered = options.category.is_some() || !options.extensions.is_empty();
    if !written.is_empty() || !filtered {
        return Ok(());
    }

    let mut described = Vec::new();
    if let Some(category) = &options.category {
        described.push(format!("category '{category}'"));
    }
    if !options.extensions.is_empty() {
        described.push(format!("extensions {:?}", options.extensions));
    }

    Err(MzLibError::Usage(format!(
        "No file in {accession} matched {}. Use list_files() to see what the project actually \
         contains — note that compressed files such as 'x.mgf.gz' have the extension '.gz'.",
        described.join(" and ")
    )))
}

// ------------------------------------------------------------------ the public surface

/// The file manifest of a PRIDE Archive project, with pyMzLib's default paging.
///
/// Paging is handled for you: however many pages the project spans, you get one list.
///
/// **This is what PRIDE's REST API publishes, which is not always everything in the project.**
/// For PXD000001 the API returns **8** files while the FTP tree holds **13** — and the five it
/// omits include the two largest, `…60min_01-20141210.mzML` (450 MB) and the matching `.mzXML`
/// (472 MB), which are exactly the modern open-format conversions most people want. The omission is
/// PRIDE's, not mzLib's: mzLib faithfully reports the v3 API, and the API simply does not list
/// them.
///
/// So a manifest that looks short may be short. If completeness matters — mirroring a project,
/// budgeting a download, proving you analysed everything — cross-check the FTP directory at
/// `https://ftp.pride.ebi.ac.uk/pride/data/archive/<year>/<month>/<accession>/`.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the accession is malformed; [`MzLibError::ProjectNotFound`] if PRIDE
/// returns nothing (a typo and a private project are indistinguishable to PRIDE, so also to us);
///
/// Note that the accession check is **grammatical only** — a letter prefix and a run of digits.
/// `"PXD0000019999"` is well-formed and costs a live round trip before it fails, so no offline
/// check will catch a transposed digit. That is deliberate: PXD accessions are not fixed-width
/// forever, and rejecting a valid future accession would be worse than one wasted request.
///
/// [`MzLibError::ServiceUnavailable`] if EBI was unreachable.
pub fn list_files(accession: &str) -> Result<Vec<PrideFile>> {
    list_files_with(accession, &ListOptions::default())
}

/// [`list_files`], with the paging and timeout stated explicitly.
///
/// # Errors
///
/// As [`list_files`], plus [`MzLibError::Usage`] if `page_size` is zero or larger than the API
/// allows.
pub fn list_files_with(accession: &str, options: &ListOptions) -> Result<Vec<PrideFile>> {
    let args = build_list_args(accession, options.page_size)?;
    let canonical = normalise_accession(accession)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    parse_manifest(&data, &canonical)
}

/// Return the COMPLETE file list of a PRIDE project, read from its FTP directory tree.
///
/// The authoritative counterpart to [`list_files`]: where that returns PRIDE's REST manifest —
/// incomplete for some projects, omitting for PXD000001 the two largest of 13 files — this walks
/// the FTP directory (subdirectories included) and returns everything the project holds (mzLib
/// #1121). Sizes are PRIDE's rounded index sizes, so [`approximate_total_size_bytes`] is an
/// estimate over the *whole* project, unlike [`total_size_bytes`].
///
/// This is a **listing** surface only. [`download`] and [`download_files`] operate on the REST
/// manifest, so a file that appears *only* here — the whole point of this function — is fetched
/// directly from its [`PrideFtpFile::url`] with an ordinary HTTPS client.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the accession is malformed; [`MzLibError::ProjectNotFound`] if no
/// project has that accession (or it lacks the publication date that locates its FTP directory), or
/// the directory listed no files; [`MzLibError::ServiceUnavailable`] if EBI was unreachable.
pub fn list_ftp_files(accession: &str) -> Result<Vec<PrideFtpFile>> {
    list_ftp_files_with(accession, Some(Duration::from_secs(300)))
}

/// [`list_ftp_files`], with the timeout stated explicitly.
///
/// # Errors
///
/// As [`list_ftp_files`].
pub fn list_ftp_files_with(
    accession: &str,
    timeout: Option<Duration>,
) -> Result<Vec<PrideFtpFile>> {
    let canonical = normalise_accession(accession)?;
    let args = [
        "pride".to_owned(),
        "ftp-files".to_owned(),
        "--accession".to_owned(),
        canonical.clone(),
    ];

    let data = match bridge::invoke(&args, None, timeout) {
        Ok(data) => data,
        // mzLib resolves the project (via the REST API) before walking the tree, so an unknown
        // accession comes back as an MzLibException, not an empty list. Map it to the same
        // ProjectNotFound that list_files raises; a mid-walk transport failure keeps its Bridge error.
        Err(MzLibError::Bridge { error_type, .. }) if error_type == "MzLibException" => {
            return Err(MzLibError::ProjectNotFound(format!(
                "PRIDE has no project '{canonical}' (or it lacks the publication date needed to \
                 locate its FTP directory). Check for a typo - a private project looks the same."
            )));
        }
        Err(other) => return Err(other),
    };

    let files = parse_ftp_files(&data)?;
    if files.is_empty() {
        return Err(MzLibError::ProjectNotFound(format!(
            "The FTP directory for '{canonical}' listed no files. Either the project is genuinely \
             empty or PRIDE's directory-index format has changed; list_files() may still return \
             its REST manifest."
        )));
    }
    Ok(files)
}

/// Download a project's files, optionally filtered, and return where they landed.
///
/// Files are streamed to a temporary name and moved into place only once complete, so an
/// interrupted download never leaves a truncated file behind.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the accession or destination is blank, or if a filter was asked for
/// and matched nothing; [`MzLibError::Bridge`] if a request failed.
pub fn download(
    accession: &str,
    destination: impl AsRef<Path>,
    options: &DownloadOptions,
) -> Result<Vec<PathBuf>> {
    let destination = destination.as_ref();
    let args = build_download_args(accession, destination, options)?;
    let canonical = normalise_accession(accession)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    let written = parse_paths(&data);
    check_filter_matched(&written, &canonical, options)?;
    Ok(written)
}

/// Download exactly the files you selected, and nothing else.
///
/// This is the counterpart to [`list_files`], and usually the one you want. [`download`]'s
/// `category` and `extensions` filters can only express what they were built to express; "under
/// 5 MB", "the three newest", or "everything except the MGF" cannot be said in that vocabulary at
/// all. They can all be said with an iterator.
///
/// The returned paths say **where each file is**, not what was transferred just now: with
/// [`DownloadOptions::overwrite`] set to `false`, a file already present is left alone and its path
/// is still returned. Do not read the length of this vector as a byte-count of work done.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the selection is empty, spans several projects, or includes a file
/// with no HTTPS location.
pub fn download_files(
    files: &[PrideFile],
    destination: impl AsRef<Path>,
    options: &DownloadOptions,
) -> Result<Vec<PathBuf>> {
    let (args, stdin) =
        build_download_files_args(files, destination.as_ref(), options.overwrite_or_default())?;
    let data = bridge::invoke(&args, Some(&stdin), options.timeout)?;
    Ok(parse_paths(&data))
}

/// Sum the sizes of some files.
///
/// **This is the size PRIDE reports, which is not always the number of bytes you will transfer.**
/// For compressed files PRIDE frequently reports the *decompressed* size: in PXD000001 the reported
/// size of `PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz` is 16,448,103 bytes, which is exactly what
/// `gzip -l` reports as its uncompressed length — the actual download is 5,984,662 bytes, 2.75×
/// smaller. The `.mztab.gz` in the same project behaves the same way; the `.xml.gz` does not.
/// PRIDE's own metadata is inconsistent here, so neither this crate nor mzLib can correct it.
///
/// Use this to answer "how much data is in this project", and treat it as an upper bound on
/// transfer time. If you need the real figure for a compressed file, only the download knows.
///
/// **It is also a sum over an incomplete manifest.** For PXD000001 this returns 0.51 GB; the
/// project on disk is **1.44 GB**, because PRIDE's API omits five files including the two largest
/// (see [`list_files`]). The two errors happen to run in opposite directions and do **not** cancel:
/// compressed sizes are over-reported, whole files are missing entirely.
#[must_use]
pub fn total_size_bytes(files: &[PrideFile]) -> u64 {
    files.iter().map(|file| file.file_size_bytes).sum()
}

/// Sum the approximate sizes of some FTP files — the honest project-size estimate.
///
/// The counterpart to [`total_size_bytes`] with the trade-offs reversed: it sums over the COMPLETE
/// FTP listing ([`list_ftp_files`]), so no files are missing, but each size is PRIDE's rounded
/// directory-index value, so the total is an estimate, not an exact byte count. For PXD000001 it
/// lands near the true 1.44 GB, where [`total_size_bytes`] reports 0.51 GB over the incomplete REST
/// manifest. For the exact bytes of one file, HTTP HEAD its [`PrideFtpFile::url`].
#[must_use]
pub fn approximate_total_size_bytes(files: &[PrideFtpFile]) -> u64 {
    files.iter().map(|file| file.approximate_size_bytes).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The recorded PXD000001 manifest, the same fixture pyMzLib's offline suite uses.
    const FIXTURE: &str = include_str!("../tests/fixtures/pride_PXD000001_files.json");

    fn recorded_manifest() -> serde_json::Value {
        serde_json::from_str(FIXTURE).expect("fixture should be valid JSON")
    }

    fn recorded_files() -> Vec<PrideFile> {
        parse_manifest(&recorded_manifest(), "PXD000001").expect("fixture should parse")
    }

    fn download_args(accession: &str, dest: &str, options: &DownloadOptions) -> Vec<String> {
        build_download_args(accession, Path::new(dest), options).expect("should assemble")
    }

    // ---------------------------------------------------------- parsing

    #[test]
    fn list_files_parses_every_file() {
        let manifest = recorded_manifest();
        let expected = manifest["files"].as_array().unwrap().len();
        assert_eq!(recorded_files().len(), expected);
    }

    #[test]
    fn file_fields_are_typed() {
        let files = recorded_files();
        let first = &files[0];
        assert!(!first.file_name.is_empty());
        assert!(first.file_size_bytes > 0);
        assert!(!first.category.is_empty());
    }

    #[test]
    fn derived_properties() {
        let files = recorded_files();
        let first = &files[0];
        assert!(
            (first.size_mb() - first.file_size_bytes as f64 / 1_000_000.0).abs() < f64::EPSILON
        );
        assert!(first.extension().starts_with('.'));
        assert_eq!(first.extension(), first.extension().to_lowercase());
        assert_eq!(first.downloadable(), first.https_url.is_some());
    }

    #[test]
    fn total_size_matches_sum() {
        let files = recorded_files();
        let expected: u64 = files.iter().map(|f| f.file_size_bytes).sum();
        assert_eq!(total_size_bytes(&files), expected);
    }

    // ---------------------------------------------------------- ftp-files (mzLib #1121)

    const FTP_FIXTURE: &str = include_str!("../tests/fixtures/pride_ftp_PXD000001.json");

    fn recorded_ftp() -> Vec<PrideFtpFile> {
        let data: serde_json::Value =
            serde_json::from_str(FTP_FIXTURE).expect("ftp fixture should be valid JSON");
        parse_ftp_files(&data).expect("ftp fixture should parse")
    }

    #[test]
    fn ftp_files_parse_every_file_including_the_rest_hidden_one() {
        let files = recorded_ftp();
        assert_eq!(files.len(), 4);
        // The whole point of the verb: a subdirectory file the REST manifest hides is present.
        assert!(files
            .iter()
            .any(|f| f.relative_path == "generated/summary.mztab"));
    }

    #[test]
    fn ftp_nested_file_keeps_its_path_but_a_bare_leaf_name() {
        let nested = recorded_ftp()
            .into_iter()
            .find(|f| f.relative_path.contains('/'))
            .expect("a nested file");
        assert_eq!(nested.relative_path, "generated/summary.mztab");
        assert_eq!(nested.file_name, "summary.mztab");
        assert_eq!(nested.extension(), ".mztab");
        assert!(nested.url.starts_with("https://"));
    }

    #[test]
    fn approximate_total_size_sums_the_complete_listing() {
        let files = recorded_ftp();
        let expected: u64 = files.iter().map(|f| f.approximate_size_bytes).sum();
        assert_eq!(approximate_total_size_bytes(&files), expected);

        let run = files
            .iter()
            .find(|f| f.file_name == "run1.raw")
            .expect("run1.raw");
        assert!(
            (run.approximate_size_mb() - run.approximate_size_bytes as f64 / 1_000_000.0).abs()
                < f64::EPSILON
        );
    }

    #[test]
    fn aspera_only_file_is_not_downloadable() {
        // A file reachable only by Aspera has no HTTPS URL, and `download` cannot fetch it. Saying
        // so up front beats a confusing failure halfway through a multi-gigabyte project.
        let aspera = PrideFile {
            file_name: "only-aspera.raw".to_owned(),
            https_url: None,
            ..recorded_files()[0].clone()
        };
        assert!(!aspera.downloadable());
    }

    #[test]
    fn files_know_which_project_they_came_from() {
        // download_files needs this to tell which project to fetch from; without it the caller
        // would have to carry the accession alongside the list.
        assert!(recorded_files()
            .iter()
            .all(|file| file.project_accession == "PXD000001"));
    }

    #[test]
    fn malformed_timestamp_becomes_none_rather_than_crashing() {
        // The caller asked for a file list, not a date. One unreadable timestamp must not fail the
        // whole manifest.
        assert_eq!(parse_timestamp("not-a-date"), None);
        assert_eq!(parse_timestamp(""), None);
        assert!(parse_timestamp("2012-02-07T00:00:00").is_some());
        assert!(parse_timestamp("2012-02-07T00:00:00Z").is_some());
        assert!(parse_timestamp("2012-02-07").is_some());
    }

    #[test]
    fn a_valid_but_unknown_accession_raises_rather_than_returning_empty() {
        // PRIDE answers an unknown accession with an empty result, not a 404. Returning an empty
        // list would let a typo produce "0 files, done" and a green exit.
        let empty = serde_json::json!({"accession": "PXD999999", "files": []});
        let error = parse_manifest(&empty, "PXD999999").unwrap_err();
        assert!(matches!(error, MzLibError::ProjectNotFound(_)));
        assert!(error.to_string().contains("typo"));
    }

    // ---------------------------------------------------------- accessions

    #[test]
    fn blank_accession_rejected_before_any_work() {
        for accession in ["", "   "] {
            let error = normalise_accession(accession).unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{accession:?}");
        }
    }

    #[test]
    fn malformed_accessions_are_rejected_not_silently_empty() {
        for accession in ["banana", "PXD", "12345", "PXD00", "PXD000001x", "-PXD1"] {
            assert!(
                normalise_accession(accession).is_err(),
                "{accession} should be rejected"
            );
        }
    }

    #[test]
    fn accession_case_and_whitespace_are_normalised() {
        // PRIDE's API is case-sensitive on the accession while category matching is not; two rules
        // pointing opposite ways is a trap, and this is the one fixable without surprise.
        for accession in ["pxd000001", "  PXD000001  ", "Pxd000001"] {
            assert_eq!(normalise_accession(accession).unwrap(), "PXD000001");
        }
    }

    #[test]
    fn nonpositive_page_size_rejected() {
        let error = build_list_args("PXD000001", 0).unwrap_err();
        assert!(error.to_string().contains("must be positive"));
    }

    #[test]
    fn a_page_size_larger_than_the_api_allows_is_refused() {
        let error = build_list_args("PXD000001", 2_147_483_648).unwrap_err();
        assert!(error.to_string().contains("larger than the API allows"));
    }

    #[test]
    fn list_args_carry_the_canonical_accession_and_page_size() {
        let args = build_list_args("  pxd000001 ", 100).unwrap();
        assert_eq!(
            args,
            vec![
                "pride",
                "files",
                "--accession",
                "PXD000001",
                "--page-size",
                "100"
            ]
        );
    }

    #[test]
    fn list_defaults_match_pymzlib() {
        let defaults = ListOptions::default();
        assert_eq!(defaults.page_size, 100);
        assert_eq!(defaults.timeout, Some(Duration::from_secs(300)));
    }

    // ---------------------------------------------------------- download, argument assembly

    #[test]
    fn download_passes_accession_and_destination() {
        let args = download_args("PXD000001", "out", &DownloadOptions::default());
        assert_eq!(args[0..2], ["pride".to_owned(), "download".to_owned()]);
        assert!(args.contains(&"PXD000001".to_owned()));
        assert!(args.contains(&"out".to_owned()));
    }

    #[test]
    fn download_omits_filters_when_not_asked_for() {
        let args = download_args("PXD000001", "out", &DownloadOptions::default());
        assert!(!args.contains(&"--category".to_owned()));
        assert!(!args.contains(&"--ext".to_owned()));
    }

    #[test]
    fn download_passes_category() {
        let options = DownloadOptions {
            category: Some("RAW".to_owned()),
            ..Default::default()
        };
        let args = download_args("PXD000001", "out", &options);
        let index = args.iter().position(|a| a == "--category").unwrap();
        assert_eq!(args[index + 1], "RAW");
    }

    #[test]
    fn download_joins_extensions_with_commas() {
        let options = DownloadOptions {
            extensions: vec![".raw".to_owned(), ".mzML".to_owned()],
            ..Default::default()
        };
        let args = download_args("PXD000001", "out", &options);
        let index = args.iter().position(|a| a == "--ext").unwrap();
        assert_eq!(args[index + 1], ".raw,.mzML");
    }

    #[test]
    fn overwrite_flag_is_not_inverted() {
        // The bridge takes --no-overwrite, the API takes overwrite. Getting the polarity wrong
        // re-downloads a project someone was resuming, or skips files they wanted replaced.
        let default_args = download_args("PXD000001", "out", &DownloadOptions::default());
        assert!(!default_args.contains(&"--no-overwrite".to_owned()));

        let options = DownloadOptions {
            overwrite: Some(false),
            ..Default::default()
        };
        assert!(download_args("PXD000001", "out", &options).contains(&"--no-overwrite".to_owned()));
    }

    #[test]
    fn download_defaults_to_no_timeout() {
        // Multi-gigabyte projects legitimately take a while; a default timeout would kill them.
        assert_eq!(DownloadOptions::default().timeout, None);
    }

    #[test]
    fn download_rejects_blank_accession_before_touching_the_network() {
        for accession in ["", "   "] {
            assert!(
                build_download_args(accession, Path::new("out"), &DownloadOptions::default())
                    .is_err()
            );
        }
    }

    #[test]
    fn blank_destination_is_refused_instead_of_writing_to_the_cwd() {
        for destination in ["", "   "] {
            let error = build_download_args(
                "PXD000001",
                Path::new(destination),
                &DownloadOptions::default(),
            )
            .unwrap_err();
            assert!(error
                .to_string()
                .contains("destination directory is required"));
        }
    }

    #[test]
    fn a_blank_category_is_refused_rather_than_selecting_everything() {
        let options = DownloadOptions {
            category: Some("  ".to_owned()),
            ..Default::default()
        };
        let error = build_download_args("PXD000001", Path::new("out"), &options).unwrap_err();
        assert!(error
            .to_string()
            .contains("Omit it to download every category"));
    }

    #[test]
    fn flag_like_filter_values_are_refused() {
        for value in ["--no-overwrite", "-x"] {
            let options = DownloadOptions {
                category: Some(value.to_owned()),
                ..Default::default()
            };
            let error = build_download_args("PXD000001", Path::new("out"), &options).unwrap_err();
            assert!(error.to_string().contains("may not begin with"), "{value}");
        }
    }

    #[test]
    fn an_extension_list_that_names_nothing_is_refused() {
        let error = normalise_extensions(&["  ".to_owned(), String::new()]).unwrap_err();
        assert!(error.to_string().contains("names no extensions"));
    }

    #[test]
    fn an_extension_containing_a_comma_is_refused() {
        let error = normalise_extensions(&[".raw,.mzML".to_owned()]).unwrap_err();
        assert!(error.to_string().contains("may not contain a comma"));
    }

    #[test]
    fn download_raises_when_a_filter_matches_nothing() {
        let options = DownloadOptions {
            category: Some("RAW".to_owned()),
            ..Default::default()
        };
        let error = check_filter_matched(&[], "PXD000001", &options).unwrap_err();
        assert!(error.to_string().contains("No file in PXD000001 matched"));
        assert!(error.to_string().contains(".gz"));
    }

    #[test]
    fn download_without_a_filter_may_legitimately_write_nothing() {
        // With --no-overwrite and everything already present, zero written files is the correct
        // answer, not an error.
        assert!(check_filter_matched(&[], "PXD000001", &DownloadOptions::default()).is_ok());
    }

    // ---------------------------------------------------------- download_files

    #[test]
    fn download_files_sends_the_selection_on_stdin() {
        // argv has a ~32 KB ceiling; a few thousand names exceed it.
        let files = recorded_files();
        let selection: Vec<PrideFile> = files
            .iter()
            .filter(|f| f.downloadable())
            .take(2)
            .cloned()
            .collect();
        let (args, stdin) = build_download_files_args(&selection, Path::new("out"), true).unwrap();

        assert!(args.contains(&"--names-from-stdin".to_owned()));
        for file in &selection {
            assert!(stdin.contains(&file.file_name), "{}", file.file_name);
        }
        assert_eq!(stdin.lines().count(), selection.len());
    }

    #[test]
    fn download_files_refuses_an_empty_selection() {
        let error = build_download_files_args(&[], Path::new("out"), true).unwrap_err();
        assert!(error.to_string().contains("No files selected"));
    }

    #[test]
    fn download_files_refuses_files_that_cannot_be_fetched() {
        let unreachable = PrideFile {
            https_url: None,
            ..recorded_files()[0].clone()
        };
        let error = build_download_files_args(&[unreachable], Path::new("out"), true).unwrap_err();
        assert!(error.to_string().contains("no HTTPS location"));
    }

    #[test]
    fn download_files_refuses_a_mixed_project_selection() {
        let files = recorded_files();
        let mut other = files[1].clone();
        other.project_accession = "PXD000002".to_owned();
        let error = build_download_files_args(&[files[0].clone(), other], Path::new("out"), true)
            .unwrap_err();
        assert!(error.to_string().contains("must come from one project"));
    }

    #[test]
    fn download_files_refuses_files_with_no_project_accession() {
        let orphan = PrideFile {
            project_accession: String::new(),
            ..recorded_files()[0].clone()
        };
        let error = build_download_files_args(&[orphan], Path::new("out"), true).unwrap_err();
        assert!(error.to_string().contains("Obtain them from list_files()"));
    }

    #[test]
    fn a_file_name_containing_a_newline_is_refused() {
        // Newline-delimited framing cannot represent it, and mis-parsing would silently select the
        // wrong files.
        let awkward = PrideFile {
            file_name: "two\nlines.raw".to_owned(),
            ..recorded_files()[0].clone()
        };
        let error = build_download_files_args(&[awkward], Path::new("out"), true).unwrap_err();
        assert!(error.to_string().contains("line break"));
    }

    #[test]
    fn download_files_honours_no_overwrite() {
        let files = recorded_files();
        let selection: Vec<PrideFile> = files
            .iter()
            .filter(|f| f.downloadable())
            .take(1)
            .cloned()
            .collect();
        let (args, _) = build_download_files_args(&selection, Path::new("out"), false).unwrap();
        assert!(args.contains(&"--no-overwrite".to_owned()));
    }

    #[test]
    fn download_files_refuses_a_blank_destination() {
        let files = recorded_files();
        assert!(build_download_files_args(&files[0..1], Path::new(""), true).is_err());
    }

    #[test]
    fn parse_paths_reads_what_was_written() {
        let data = serde_json::json!({"paths": ["out/a.raw", "out/b.raw"]});
        assert_eq!(
            parse_paths(&data),
            vec![PathBuf::from("out/a.raw"), PathBuf::from("out/b.raw")]
        );
        assert!(parse_paths(&serde_json::json!({})).is_empty());
    }
}
