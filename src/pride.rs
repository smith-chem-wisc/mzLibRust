//! PRIDE Archive access, backed by mzLib's `PrideArchiveClient`.
//!
//! | You want to | Call | mzLib type |
//! |---|---|---|
//! | Find datasets about a subject | [`search`] | `PrideArchiveClient.SearchProjectsAsync` |
//! | See a project's files, with categories and checksums | [`list_files`] | `PrideArchiveClient.GetProjectFilesAsync` |
//! | See **every** file, including ones the manifest omits | [`list_ftp_files`] | `PrideArchiveClient.GetProjectFilesFromFtpAsync` |
//! | Size a project before downloading it | [`total_size_bytes`], [`approximate_total_size_bytes`] | `TotalApproximateSizeBytes` |
//! | Download the files you picked | [`download_files`] | `PrideArchiveClient.DownloadFilesAsync` |
//! | Download by category or extension | [`download`] | `WhereCategory`, `WhereExtension` |
//! | Check a download against the MD5 and exact size PRIDE recorded | [`file_checksums`] | `PrideArchiveClient.GetFileChecksumsAsync` |
//!
//! Every example in this guide runs in CI against answers recorded from the live PRIDE Archive,
//! replayed by a stand-in bridge; the one that downloads says so.
//!
//! The [PRIDE Archive](https://www.ebi.ac.uk/pride/archive/) is EBI's public proteomics data
//! repository. This module lets a Rust user list what is in a project and pull files down, using
//! the same paging, URL-resolution, and safe-download logic that mzLib uses in C#.
//!
//! ```
//! # mzlib_replay::activate();
//! let files = mzlib::pride::list_files("PXD000001")?;
//! assert_eq!(files.len(), 8);
//! // 0.51 GB is PRIDE's REST manifest, which omits files, the two largest among them; see
//! // list_files.
//! assert_eq!(mzlib::pride::total_size_bytes(&files), 514_278_049);
//!
//! // Filter however you like — the full expressiveness of Rust — then fetch exactly that.
//! let small: Vec<_> = files
//!     .iter()
//!     .filter(|f| f.size_mb() < 5.0 && f.downloadable())
//!     .cloned()
//!     .collect();
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Don't have an accession yet? [`search`] is the discovery step that produces them — every page
//! fetched, no accession repeated, and each hit saying which fields matched:
//!
//! ```
//! # mzlib_replay::activate();
//! let hits = mzlib::pride::search("plasmodium falciparum schizont")?;
//! assert_eq!(hits[0].accession, "PXD070842");
//! assert_eq!(hits[0].matched_fields(), ["references", "title"]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! A hit is PRIDE's search projection, not the project's metadata: controlled vocabulary is
//! flattened to display strings, and a zero or an empty list means "not reported". Follow the
//! accession to [`list_files`] for what you can act on.
//!
//! An accession that does not exist is an error, not an empty list, so a typo cannot pass for a
//! project with no files:
//!
//! ```
//! # mzlib_replay::activate();
//! let error = mzlib::pride::list_files("PXD999999999").unwrap_err();
//! assert!(matches!(error, mzlib::MzLibError::ProjectNotFound(_)));
//! ```
//!
//! The REST manifest is not the whole project. The FTP tree is, subdirectories included, and for
//! PXD000001 it holds an mzML the manifest leaves out:
//!
//! ```
//! # mzlib_replay::activate();
//! let manifest = mzlib::pride::list_files("PXD000001")?;
//! let tree = mzlib::pride::list_ftp_files("PXD000001")?;
//! assert_eq!((tree.len(), manifest.len()), (14, 8));
//! let hidden: Vec<&str> = tree
//!     .iter()
//!     .filter(|f| f.extension() == ".mzml" && !manifest.iter().any(|m| m.file_name == f.file_name))
//!     .map(|f| f.relative_path.as_str())
//!     .collect();
//! assert_eq!(hidden, ["TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML"]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Fetching files is one more call.
//!
//! Not run: it downloads from EBI.
//!
//! ```no_run
//! # let small: Vec<mzlib::pride::PrideFile> = Vec::new();
//! mzlib::pride::download_files(&small, "downloads", &Default::default())?;
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Each file streams to a temporary `.partial` name and moves into place only once complete, so a
//! cancelled or crashed download leaves no truncated file under its real name for a later run to
//! mistake for a good one. A transfer that fails after mzLib's retries leaves `<name>.partial` and
//! `<name>.partial.validator` in the folder on purpose: the next call for that file resumes from
//! them instead of starting again (mzLib#1412). Leave them in place, or delete both to start over.
//!
//! **Check large files yourself.** mzLib can check a download against PRIDE's checksums, and it
//! also catches an EBI error page served at the right size (mzLib#1427, #1449). That check is not
//! yet on the path [`download`] and [`download_files`] use, so a file can be reported as
//! downloaded when its bytes are wrong. For files that matter, check each one against
//! [`file_checksums`] before you search it, as the next section shows. With
//! [`DownloadOptions::overwrite`] set to `false`, a file already present is skipped and not
//! checked either, so a damaged copy stays damaged.
//!
//! ## Checking a download against PRIDE's checksums
//!
//! **Is the file on your disk the file PRIDE holds?** [`download`] cannot tell you yet. PRIDE
//! records an MD5 and the exact size of every file a submitter deposits, and [`file_checksums`]
//! returns that record, keyed by file name:
//!
//! ```
//! # mzlib_replay::activate();
//! let checksums = mzlib::pride::file_checksums("PXD000001")?;
//! assert_eq!(checksums.len(), 12);
//! let raw = &checksums["TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw"];
//! assert_eq!(raw.md5, "c567cecbb9866796e690a4958746f7b6");
//! assert_eq!(raw.size_bytes, 220_475_548);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Look a file up by the name it was downloaded under. Each value is a [`PrideFileChecksum`]:
//!
//! - `md5` is 32 hexadecimal characters, as PRIDE wrote them. Compare it without regard to case.
//! - `size_bytes` is the exact length a download of the file comes to. Neither other listing gives
//!   you that: [`list_ftp_files`] rounds its sizes, and [`list_files`] can be wrong. Because each
//!   size is exact, so is their sum: what the deposited files weigh, which the other two listings
//!   can only estimate. It leaves out the `generated/` files, which have no row (below).
//!
//! ```
//! # mzlib_replay::activate();
//! let files = mzlib::pride::list_files("PXD000001")?;
//! let checksums = mzlib::pride::file_checksums("PXD000001")?;
//! let rest = files
//!     .iter()
//!     .find(|f| f.file_name == "PRIDE_Exp_Complete_Ac_22134.xml.gz")
//!     .unwrap();
//! // The download is 10,668,000 bytes: the checksum list's number, not the manifest's.
//! assert_eq!(rest.file_size_bytes, 10_677_205);
//! assert_eq!(checksums[&rest.file_name].size_bytes, 10_668_000);
//!
//! let deposited: u64 = checksums.values().map(|c| c.size_bytes).sum();
//! assert_eq!(format!("{:.2} GB", deposited as f64 / 1e9), "1.43 GB");
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! ### The check
//!
//! A file is PRIDE's when its length and its MD5 both match. The length costs nothing to check.
//! The MD5 reads the whole file, which takes a few seconds per gigabyte. Check both, because the
//! length alone does not prove the bytes: on 7 and 8 October 2026, EBI served about 23 of about 70
//! large transfers at the right size with its own error text written over part of them
//! (mzLib#1449). This crate has no digest of its own; the example uses the
//! [`md-5`](https://crates.io/crates/md-5) crate. Here it is on two damaged copies of the project's
//! 1,645-byte `README.txt`. The first was cut short. The second is the right size with EBI's error
//! text inside, and only the MD5 catches it:
//!
//! ```
//! # mzlib_replay::activate();
//! use std::path::Path;
//!
//! use md5::{Digest, Md5};
//! use mzlib::pride::PrideFileChecksum;
//!
//! /// None when the file is the one PRIDE recorded, otherwise what is wrong with it.
//! fn check(path: &Path, expected: &PrideFileChecksum) -> std::io::Result<Option<String>> {
//!     let size = std::fs::metadata(path)?.len();
//!     if size != expected.size_bytes {
//!         return Ok(Some(format!("{size} bytes where PRIDE records {}", expected.size_bytes)));
//!     }
//!     let mut md5 = Md5::new();
//!     std::io::copy(&mut std::fs::File::open(path)?, &mut md5)?;
//!     let actual = format!("{:x}", md5.finalize());
//!     if !actual.eq_ignore_ascii_case(&expected.md5) {
//!         return Ok(Some(format!("MD5 {actual} where PRIDE records {}", expected.md5)));
//!     }
//!     Ok(None)
//! }
//!
//! let checksums = mzlib::pride::file_checksums("PXD000001")?;
//! let readme = &checksums["README.txt"];
//! let damaged = Path::new("README.txt");
//!
//! std::fs::write(damaged, [b'x'; 1000])?;
//! assert_eq!(
//!     check(damaged, readme)?.as_deref(),
//!     Some("1000 bytes where PRIDE records 1645")
//! );
//!
//! let mut error_page = b"<Error><Code>ConnectionClosedException</Code></Error>".to_vec();
//! error_page.resize(1645, b'.');
//! std::fs::write(damaged, error_page)?;
//! assert_eq!(
//!     check(damaged, readme)?.as_deref(),
//!     Some(
//!         "MD5 6a8f1a7d39089d34bf9bc9b19e553dd4 where PRIDE records \
//!          492753137013825025e6729c8798664f"
//!     )
//! );
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! After a real download, check every file it wrote.
//!
//! Not run: it downloads from PRIDE.
//!
//! ```no_run
//! # use std::path::Path;
//! # fn check(_: &Path, _: &mzlib::pride::PrideFileChecksum) -> std::io::Result<Option<String>> {
//! #     Ok(None)
//! # }
//! use mzlib::pride::{download, file_checksums, DownloadOptions};
//!
//! let options = DownloadOptions {
//!     category: Some("RAW".to_owned()),
//!     ..Default::default()
//! };
//! let paths = download("PXD000001", "downloads", &options)?;
//! let checksums = file_checksums("PXD000001")?;
//! for path in &paths {
//!     let name = path.file_name().unwrap().to_string_lossy();
//!     match checksums.get(name.as_ref()) {
//!         None => println!("{name}: PRIDE has no checksum for it, so it cannot be checked"),
//!         Some(expected) => {
//!             if let Some(problem) = check(path, expected)? {
//!                 println!("{name}: {problem}. Download it again.");
//!             }
//!         }
//!     }
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! To download a damaged file again, pass it to [`download_files`], or call [`download`] with
//! [`DownloadOptions::overwrite`] left at its default. `Some(false)` would skip it, because it is
//! already there.
//!
//! ### What the list covers
//!
//! It covers the files the submitter deposited at the project's root. It is a third list, not a
//! view of either of the others, and three things follow from that.
//!
//! **PRIDE's own `generated/` files have no row.** A file without a row cannot be checked; that
//! does not make it invalid. In PXD000001 those are the two files under `generated/`, and
//! [`list_files`] names the same two without their folder:
//!
//! ```
//! # mzlib_replay::activate();
//! let checksums = mzlib::pride::file_checksums("PXD000001")?;
//! let tree = mzlib::pride::list_ftp_files("PXD000001")?;
//! let unchecked: Vec<&str> = tree
//!     .iter()
//!     .map(|f| f.relative_path.as_str())
//!     .filter(|path| !checksums.contains_key(*path))
//!     .collect();
//! assert_eq!(
//!     unchecked,
//!     [
//!         "generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz",
//!         "generated/PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz",
//!     ]
//! );
//!
//! let manifest = mzlib::pride::list_files("PXD000001")?;
//! let unlisted: Vec<&str> = manifest
//!     .iter()
//!     .map(|f| f.file_name.as_str())
//!     .filter(|name| !checksums.contains_key(*name))
//!     .collect();
//! assert_eq!(
//!     unlisted,
//!     [
//!         "PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz",
//!         "PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz",
//!     ]
//! );
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **It names files the REST manifest omits,** such as the 450 MB mzML that only the FTP listing
//! showed above:
//!
//! ```
//! # mzlib_replay::activate();
//! let checksums = mzlib::pride::file_checksums("PXD000001")?;
//! assert!(checksums
//!     .contains_key("TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **It can name files that are not deposits at all.** PXD015239's list includes a 2.3 GB NFS
//! temp file, `.nfs80360b91005fae5300002478`. Look up the files you have by name; do not read the
//! list as the project's contents.
//!
//! ### Two traps
//!
//! **[`PrideFile::checksum`] is not this MD5.** Where PRIDE fills it in, it is a 40-character
//! SHA-1, a different digest of the same bytes. For every PXD000001 file it is empty:
//!
//! ```
//! # mzlib_replay::activate();
//! let files = mzlib::pride::list_files("PXD000001")?;
//! assert!(files.iter().all(|f| f.checksum.is_empty()));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **An empty map is an answer.** It means the project exists but PRIDE holds no checksums for
//! it, so nothing can be checked. An accession PRIDE does not know is
//! [`MzLibError::ProjectNotFound`], and so is a project not yet published: PRIDE holds no
//! checksums before publication. A list mzLib cannot read (a new header, a malformed row, a name
//! listed twice) is [`MzLibError::Bridge`] with `error_type` `"MzLibException"`. It is not your
//! bug, and mzLib does not retry it, so report it.
//!
//! ## Retries and timeouts
//!
//! **mzLib retries for you.** A request that fails in a way that can pass (no answer, an HTTP 408,
//! 429 or 5xx, or a 403 from EBI's FTP host, its rate limit) is retried three times, after 5, 20
//! and 60 seconds, and that request alone (mzLib#1412). A paged manifest retries only the page
//! that failed, and a download resumes where it stopped.
//!
//! **This crate adds no timeout of its own by default**, on any PRIDE function. Each attempt is
//! already limited to 100 seconds, so an outage that hangs is reported as
//! [`MzLibError::ServiceUnavailable`] after about 8 minutes per request. Set a timeout
//! ([`ListOptions::timeout`], [`SearchOptions::timeout`], [`DownloadOptions::timeout`], or the
//! argument of [`list_ftp_files_with`] and [`file_checksums_with`]) to cut a call off sooner; it then fails with
//! [`MzLibError::Timeout`], which is not `ServiceUnavailable`.
//!
//! By the time `ServiceUnavailable` reaches you, mzLib has already retried, so retrying again at
//! once rarely helps. Wait minutes, not seconds. Two failures arrive as [`MzLibError::Bridge`]
//! with `error_type` `"HttpRequestException"` although they are EBI's: an FTP listing that came
//! back empty, and a 403 from EBI's FTP host, each still failing after mzLib's retries. Treat
//! those as an outage too.
//!
//! Not run: shows how to handle a PRIDE outage, which a recording cannot produce.
//!
//! ```no_run
//! use std::time::Duration;
//!
//! let files = match mzlib::pride::list_files("PXD000001") {
//!     Err(mzlib::MzLibError::ServiceUnavailable { .. }) => {
//!         // mzLib retried already; give EBI ten minutes.
//!         std::thread::sleep(Duration::from_secs(600));
//!         mzlib::pride::list_files("PXD000001")?
//!     }
//!     other => other?,
//! };
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! ## Cite
//!
//! Cite mzLib, this crate, and the archive:
//!
#![doc = include_str!("../docs/reference/cite.pride.md")]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Duration;

use chrono::{DateTime, FixedOffset, NaiveDate};
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
    /// Size in bytes as reported by PRIDE. For compressed files PRIDE often reports the
    /// **decompressed** size (PXD000001's `.mgf.gz`: 16,448,103 reported, 5,984,662 transferred).
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
/// PRIDE's REST manifest ([`list_files`]) is knowingly incomplete: for PXD000001 on 2026-10-21 it
/// listed 9 of the 14 files in the FTP tree, omitting five including the two largest. When
/// completeness or a true project size matters, [`list_ftp_files`] walks the directory
/// (subdirectories included) and returns everything (mzLib #1121). The trade-off is the size: [`PrideFtpFile::approximate_size_bytes`]
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

/// The MD5 and exact size PRIDE recorded for one file deposited in a project.
///
/// [`file_checksums`] returns these, keyed by file name. They are what a downloaded file is
/// checked against: it is the file PRIDE holds when its length is [`Self::size_bytes`] and its
/// MD5 is [`Self::md5`]. The [module guide](crate::pride#checking-a-download-against-prides-checksums)
/// shows the check.
///
/// Read strictly, unlike the other PRIDE rows: a row missing its MD5 or its size fails the call
/// with [`MzLibError::Protocol`]. Defaulting it to `""` or `0` would turn every later check of
/// that file into a mismatch that blames the download.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct PrideFileChecksum {
    /// The bare file name, spelled as [`PrideFile::file_name`] spells it. Never a path.
    /// Case-sensitive: `README.txt` and `readme.txt` would be two files.
    pub file_name: String,
    /// The file's MD5 as 32 hexadecimal characters, as PRIDE wrote it. Lower case in every list
    /// seen, but not rewritten, so compare it without regard to case, as mzLib's own checks do.
    /// Not [`PrideFile::checksum`], which is a different digest.
    pub md5: String,
    /// The file's exact size in bytes: the length a download of it comes to. Compare a
    /// downloaded file's length with this, not with [`PrideFile::file_size_bytes`], which can
    /// differ. For PXD000001's `PRIDE_Exp_Complete_Ac_22134.xml.gz` that one is 10,677,205 bytes
    /// and this one is 10,668,000, the length the download has.
    pub size_bytes: u64,
    /// The accession the list was fetched for.
    ///
    /// Not on the wire: stamped on by [`file_checksums`], as [`PrideFile::project_accession`] is.
    #[serde(skip)]
    pub project_accession: String,
}

/// How a manifest is fetched. By default the page size is mzLib's, and so is the time limit.
#[derive(Debug, Clone, Default)]
pub struct ListOptions {
    /// How many files to request per underlying API call. Only affects how the manifest is
    /// fetched, never what you get back. `None` sends none, so mzLib's default applies; a value
    /// must be at least 1.
    pub page_size: Option<u32>,
    /// Time to allow for the whole fetch.
    ///
    /// `None` (the default) adds no limit of this crate's own, because mzLib already bounds every
    /// request: it retries a transient failure three times, after 5, 20 and 60 s, and each attempt
    /// can take up to 100 s (mzLib#1412). An outage that hangs is therefore reported as
    /// [`MzLibError::ServiceUnavailable`] after about 8 minutes per request. A value cuts the whole
    /// call off after that long with [`MzLibError::Timeout`].
    pub timeout: Option<Duration>,
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

/// The wire's `--page-size` pair, or nothing for `None` so mzLib's own default applies.
fn page_size_args(page_size: Option<u32>) -> Result<Vec<String>> {
    let Some(page_size) = page_size else {
        return Ok(Vec::new());
    };
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
    Ok(vec!["--page-size".to_owned(), page_size.to_string()])
}

/// The argv for `pride files`.
fn build_list_args(accession: &str, page_size: Option<u32>) -> Result<Vec<String>> {
    let canonical = normalise_accession(accession)?;
    let paging = page_size_args(page_size)?;

    let mut args = vec![
        "pride".to_owned(),
        "files".to_owned(),
        "--accession".to_owned(),
        canonical,
    ];
    args.extend(paging);
    Ok(args)
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

/// Turn the `pride checksums` payload into rows keyed by file name, stamping the accession.
fn parse_checksums(
    data: &serde_json::Value,
    accession: &str,
) -> Result<BTreeMap<String, PrideFileChecksum>> {
    let raw = data
        .get("files")
        .cloned()
        .unwrap_or(serde_json::Value::Null);
    let rows: Vec<PrideFileChecksum> = if raw.is_null() {
        Vec::new()
    } else {
        serde_json::from_value(raw).map_err(|error| {
            MzLibError::Protocol(format!("PRIDE checksum list could not be read: {error}"))
        })?
    };
    Ok(rows
        .into_iter()
        .map(|mut row| {
            row.project_accession = accession.to_owned();
            (row.file_name.clone(), row)
        })
        .collect())
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

/// The file manifest of a PRIDE Archive project, with mzLib's default paging.
///
/// Paging is handled for you: however many pages the project spans, you get one list.
///
/// **This is what PRIDE's REST API publishes, which is not always everything in the project.**
/// On 2026-10-21, for PXD000001, the API listed **9** of the **14** files in the FTP tree, and the
/// five it omits include the two largest: `…60min_01-20141210.mzML` (450 MB) and the matching
/// `.mzXML` (472 MB), which are exactly the modern open-format conversions most people want. The
/// recording this crate's examples replay is from 2026-10-03, when the API listed 8. The omission is
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

/// List every file of a PRIDE project from PRIDE's REST manifest, with the paging and timeout
/// stated explicitly.
///
/// The manifest is knowingly incomplete for some projects; see [`list_files`] and
/// [`list_ftp_files`]. [`ListOptions::page_size`] is sent only when set; otherwise mzLib's own
/// default applies. A malformed accession, or a `page_size` of zero or larger than the API allows,
/// is refused before anything is spawned.
#[doc = include_str!("../docs/reference/pride.files.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::ProjectNotFound`] when PRIDE returns an empty manifest: an unknown accession and a
/// private project both come back empty rather than as an error, and an empty list would be
/// indistinguishable from a project that genuinely has nothing matching.
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::pride::{list_files_with, ListOptions};
///
/// let files = list_files_with("PXD000001", &ListOptions::default())?;
/// let raw: Vec<_> = files.iter().filter(|f| f.category == "RAW").collect();
/// assert!(!raw.is_empty());
/// assert!(files.iter().all(|f| f.project_accession == "PXD000001"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/pride.files.see-also.md")]
pub fn list_files_with(accession: &str, options: &ListOptions) -> Result<Vec<PrideFile>> {
    let args = build_list_args(accession, options.page_size)?;
    let canonical = normalise_accession(accession)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    parse_manifest(&data, &canonical)
}

/// Return the COMPLETE file list of a PRIDE project, read from its FTP directory tree.
///
/// The authoritative counterpart to [`list_files`]: where that returns PRIDE's REST manifest —
/// incomplete for some projects, omitting for PXD000001 five of the project's 14 files on
/// 2026-10-21, including the two largest — this walks
/// the FTP directory (subdirectories included) and returns everything the project holds (mzLib
/// #1121). Sizes are PRIDE's rounded index sizes, so [`approximate_total_size_bytes`] is an
/// estimate over the *whole* project, unlike [`total_size_bytes`].
///
/// This is a **listing** surface only. [`download`] and [`download_files`] operate on the REST
/// manifest, so a file that appears *only* here — the whole point of this function — is fetched
/// directly from its [`PrideFtpFile::url`] with an ordinary HTTPS client.
///
/// It adds no time limit of this crate's own, because mzLib already bounds every request (see
/// [`ListOptions::timeout`]); [`list_ftp_files_with`] takes one.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the accession is malformed; [`MzLibError::ProjectNotFound`] if no
/// project has that accession (or it lacks the publication date that locates its FTP directory), or
/// the walk found no files; [`MzLibError::ServiceUnavailable`] if EBI was unreachable.
///
/// An FTP listing that comes back empty used to reach this crate as no files, and so as
/// `ProjectNotFound`. mzLib now retries it and then fails (mzLib#1412), so it arrives as
/// [`MzLibError::Bridge`] with `error_type` `"HttpRequestException"`, as does a 403 from EBI's FTP
/// host that outlasts the retries. Both are EBI's failures: treat them as an outage.
pub fn list_ftp_files(accession: &str) -> Result<Vec<PrideFtpFile>> {
    list_ftp_files_with(accession, None)
}

/// Return the COMPLETE file list of a PRIDE project from its FTP directory tree, with the timeout
/// stated explicitly.
///
/// See [`list_ftp_files`] for when to prefer this to [`list_files`]. `timeout` is the time to allow
/// for the whole walk, which spans one request per directory. `None` adds no limit of this crate's
/// own, because mzLib already bounds every request: it retries a transient failure three times,
/// after 5, 20 and 60 s, and each attempt can take up to 100 s (mzLib#1412). An outage that hangs
/// is therefore reported as [`MzLibError::ServiceUnavailable`] after about 8 minutes per request. A
/// value cuts the whole walk off after that long with [`MzLibError::Timeout`].
#[doc = include_str!("../docs/reference/pride.ftp-files.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::ProjectNotFound`] if no project has that accession (or it lacks the publication
/// date that locates its FTP directory), or the walk found no files.
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use std::time::Duration;
///
/// // Give up after a minute, with MzLibError::Timeout, rather than wait out mzLib's retries.
/// let files = mzlib::pride::list_ftp_files_with("PXD000001", Some(Duration::from_secs(60)))?;
/// assert_eq!(files.len(), 14);
/// assert_eq!(files[0].relative_path, "F063721.dat");
/// assert!(files[0].url.starts_with("https://ftp.pride.ebi.ac.uk/"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/pride.ftp-files.see-also.md")]
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

/// The first pyMzLib whose bridge has `pride checksums`: its spec's `since.pymzlib`.
const CHECKSUMS_SINCE: &str = "0.6.0";

/// mzLib's wording for an unknown accession (`PrideArchiveClient.GetProjectAsync` at the pin).
/// Every verb that resolves the project first throws it, as an `MzLibException`.
const UNKNOWN_PROJECT_MARKER: &str = "has no project with accession";

/// The MD5 and exact size PRIDE recorded for each file deposited in a project, keyed by file name.
///
/// This is what a download is checked against: a file is the one PRIDE holds when its length is
/// [`PrideFileChecksum::size_bytes`] and its MD5 is [`PrideFileChecksum::md5`]. [`download`] and
/// [`download_files`] check neither yet, so for files that matter, check them yourself; the
/// [module guide](crate::pride#checking-a-download-against-prides-checksums) shows how.
///
/// **This is a third list, not a view of either of the others.** It covers the files the
/// submitter deposited at the project's root:
///
/// - PRIDE's own `generated/` files have no row. A file without a row cannot be checked; that
///   does not make it invalid.
/// - It names files [`list_files`] omits: 5 of PXD000001's 12 rows on 2026-10-21.
/// - It can name files that are not deposits at all. PXD015239's lists an NFS temp file,
///   `.nfs80360b91005fae5300002478`. Look up the files you have by name; do not read the list as
///   the project's contents.
///
/// [`PrideFile::checksum`] is not this MD5. Where PRIDE fills it in, it is a 40-character SHA-1,
/// and for every PXD000001 file it is empty.
///
/// The map is in file-name order, so iterating it is deterministic. **It is empty when the
/// project exists but PRIDE holds no checksums for it**: that is an answer, not an error, though
/// nothing can then be checked.
///
/// It adds no time limit of this crate's own, because mzLib already bounds every request (see
/// [`ListOptions::timeout`]); an empty list costs a second request, the project lookup.
/// [`file_checksums_with`] takes a limit.
///
/// ```
/// # mzlib_replay::activate();
/// let checksums = mzlib::pride::file_checksums("PXD000001")?;
/// assert_eq!(checksums.len(), 12);
/// let readme = &checksums["README.txt"];
/// assert_eq!(readme.md5, "492753137013825025e6729c8798664f");
/// assert_eq!(readme.size_bytes, 1645);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors
///
/// [`MzLibError::Usage`] if the accession is blank or malformed, or the bridge is older than
/// pyMzLib 0.6.0's; [`MzLibError::ProjectNotFound`] if PRIDE has no public project with that
/// accession (a project not yet published is reported the same way, since PRIDE holds no
/// checksums before publication); [`MzLibError::ServiceUnavailable`] if PRIDE was down,
/// rate-limiting or timed out, after mzLib's retries; [`MzLibError::Bridge`] with `error_type`
/// `"MzLibException"` if PRIDE sent a list mzLib cannot read (a new header, a malformed row, a
/// name listed twice), or another HTTP error.
pub fn file_checksums(accession: &str) -> Result<BTreeMap<String, PrideFileChecksum>> {
    file_checksums_with(accession, None)
}

/// The MD5 and exact size PRIDE recorded for each file of a project, with the timeout stated
/// explicitly.
///
/// See [`file_checksums`] for what the list covers and what it does not. The accession is trimmed
/// and upper-cased first, because PRIDE's checksum route is case-sensitive (`pxd000001` is an
/// unknown accession there). `timeout` is the time to allow for the whole call. `None` adds no
/// limit of this crate's own, because mzLib already bounds every request: it retries a transient
/// failure three times, after 5, 20 and 60 s, and each attempt can take up to 100 s (mzLib#1412).
/// An outage that hangs is therefore reported as [`MzLibError::ServiceUnavailable`] after about 8
/// minutes per request, and an empty list costs a second request, the project lookup. A value cuts
/// the whole call off after that long with [`MzLibError::Timeout`].
#[doc = include_str!("../docs/reference/pride.checksums.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::ProjectNotFound`] for the one `MzLibException` that means PRIDE has no public
/// project with that accession. Any other `MzLibException` (a list mzLib cannot read) stays
/// [`MzLibError::Bridge`]: reporting it as a typo would hide a broken contract.
///
/// [`MzLibError::Usage`] before anything is spawned for a malformed accession, and when the bridge
/// in use is older than pyMzLib 0.6.0's, which has no `pride checksums`.
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use std::time::Duration;
///
/// // Give up after ten minutes, with MzLibError::Timeout, rather than wait out every retry.
/// let checksums =
///     mzlib::pride::file_checksums_with("PXD000001", Some(Duration::from_secs(600)))?;
/// let raw = &checksums["TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw"];
/// assert_eq!(raw.md5, "c567cecbb9866796e690a4958746f7b6");
/// assert_eq!(raw.size_bytes, 220_475_548);
/// assert_eq!(raw.project_accession, "PXD000001");
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/pride.checksums.see-also.md")]
pub fn file_checksums_with(
    accession: &str,
    timeout: Option<Duration>,
) -> Result<BTreeMap<String, PrideFileChecksum>> {
    let canonical = normalise_accession(accession)?;
    bridge::require_verb("pride checksums", CHECKSUMS_SINCE)?;
    fetch_checksums(&canonical, timeout)
}

/// [`file_checksums_with`] after the accession and the bridge have been checked: the call itself.
fn fetch_checksums(
    canonical: &str,
    timeout: Option<Duration>,
) -> Result<BTreeMap<String, PrideFileChecksum>> {
    let data = bridge::invoke(&checksums_args(canonical), None, timeout)
        .map_err(|error| checksums_error(error, canonical))?;
    parse_checksums(&data, canonical)
}

/// The `pride checksums` call for a canonical accession: nothing but the accession is sent.
fn checksums_args(canonical: &str) -> Vec<String> {
    vec![
        "pride".to_owned(),
        "checksums".to_owned(),
        "--accession".to_owned(),
        canonical.to_owned(),
    ]
}

/// The error [`file_checksums_with`] reports for one the bridge returned.
///
/// PRIDE answers an unknown accession with an empty list, so mzLib asks for the project and throws
/// `MzLibException` when there is none. Only that wording is re-mapped: a list mzLib cannot read is
/// an `MzLibException` too, and reporting it as a typo would hide a broken contract.
/// `ServiceUnavailable` has its own variant and passes through unchanged.
fn checksums_error(error: MzLibError, canonical: &str) -> MzLibError {
    match error {
        MzLibError::Bridge {
            error_type,
            message,
        } if error_type == "MzLibException" && message.contains(UNKNOWN_PROJECT_MARKER) => {
            MzLibError::ProjectNotFound(format!(
                "PRIDE has no public project '{canonical}'. Check for a typo; a project not yet \
                 published looks the same, and PRIDE holds no checksums before publication."
            ))
        }
        other => other,
    }
}

/// Download a project's files, optionally filtered, and return where they landed.
///
/// Files are streamed to a temporary name and moved into place only once complete, so an
/// interrupted download never leaves a truncated file under its real name; one that fails after
/// mzLib's retries leaves a `.partial` pair the next call resumes from (see the
/// [module guide](crate::pride)). To fetch a selection you made yourself, use [`download_files`].
#[doc = include_str!("../docs/reference/pride.download.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] if a `category` or `extensions` filter was asked for and matched
/// nothing: an empty download is reported as an error here, not an empty success.
///
/// # Examples
///
/// Not run: it downloads from EBI.
///
/// ```no_run
/// use mzlib::pride::{download, DownloadOptions};
///
/// let written = download(
///     "PXD000001",
///     "downloads",
///     &DownloadOptions {
///         category: Some("PEAK".to_owned()),
///         extensions: vec![".gz".to_owned()],
///         ..Default::default()
///     },
/// )?;
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/pride.download.see-also.md")]
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
/// **It is also a sum over an incomplete manifest.** For PXD000001's recorded manifest (8 files,
/// 2026-10-03) this returns 0.51 GB; the project on disk is **1.44 GB**, because PRIDE's API omits
/// files including the two largest (see [`list_files`]). The two errors happen to run in opposite directions and do **not** cancel:
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

// ------------------------------------------------------------------ search

/// The longest keyword [`search_with`] sends. PRIDE answers a longer one with HTTP 500, which
/// cannot be told apart from an outage, so mzLib refuses it first and so does this crate. Mirrors
/// mzLib's `PrideArchiveClient.MaxKeywordLength`.
pub const MAX_KEYWORD_LENGTH: usize = 1000;

/// Downloads of a project in one year, as PRIDE reports them.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct YearlyDownloads {
    /// The year, as PRIDE writes it, e.g. `"2025"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub year: String,
    /// Downloads in that year.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub count: u64,
}

/// One hit from [`search`].
///
/// **This is not a project's full metadata, and the two are not interchangeable.** PRIDE serves
/// search from an Elasticsearch projection in which every controlled-vocabulary field is flattened
/// to a display string: the same project reports its instruments as `["Q Exactive"]` here and as
/// terms with accessions elsewhere, contacts collapse to display names, and publications to one
/// pre-formatted citation string. That is PRIDE's wire, not a simplification chosen here — the
/// accessions are simply not sent. Follow [`Self::accession`] when you need the vocabulary.
///
/// **A zero or an empty list means "not reported", never a measured zero.** PRIDE omits nothing as
/// null, so absence arrives as `0`, `""` or `[]`, and several fields are genuinely sparse: sampled
/// over 1,600 hits, `project_tags` was filled on 2.6%, `sdrf` on 2.4%, `other_omics_links` on 18%,
/// and the bot/hub/organic trio on under half. A `download_count` of 0 does not mean nobody
/// downloaded it.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct PrideProjectSearchResult {
    /// The project accession: usually `PXD…`, but search also returns legacy `PRD` and affinity
    /// `PAD` projects. Treat it as an opaque key, and pass it to [`list_files`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub accession: String,
    /// The project title.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub title: String,
    /// The submitter's free-text description.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub project_description: String,
    /// How the samples were prepared, as free text.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_processing_protocol: String,
    /// How the data were searched and processed, as free text.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub data_processing_protocol: String,
    /// The dataset DOI, or `""` when PRIDE has minted none.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub doi: String,
    /// `"COMPLETE"`, `"PARTIAL"`, `"AFFINITY"`, or `"PRIDE"` for legacy submissions. Not a closed
    /// set: PRIDE has added values before.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub submission_type: String,
    /// The project's SDRF term **values**, flattened by the search index into one space-joined
    /// string. Not a file, name or URL — nothing can be fetched with it, and the row and column
    /// structure is gone; `""` on most hits. For a real SDRF see [`crate::sdrf`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sdrf: String,
    /// A bare calendar date: no time and no offset, because PRIDE sent none, so this is never
    /// widened to a timestamp (unlike [`PrideFile::submission_date`]). `None`: PRIDE reported no
    /// date.
    #[serde(default, deserialize_with = "lenient_date")]
    pub submission_date: Option<NaiveDate>,
    /// As [`Self::submission_date`]. `None`: PRIDE reported no date.
    #[serde(default, deserialize_with = "lenient_date")]
    pub publication_date: Option<NaiveDate>,
    /// As [`Self::submission_date`]. `None`: PRIDE reported no date.
    #[serde(default, deserialize_with = "lenient_date")]
    pub updated_date: Option<NaiveDate>,
    /// PRIDE's coarse classification tags. Sparse.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub project_tags: Vec<String>,
    /// The submitter's keywords. **May contain empty and whitespace-only strings** — about 9% of
    /// hits carry one. Passed through rather than filtered, so every binding reports the same
    /// keywords; filter before joining.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub keywords: Vec<String>,
    /// Submitters' display names, flattened from structured contacts.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub submitters: Vec<String>,
    /// Lab heads' display names.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub lab_pis: Vec<String>,
    /// Affiliations, as display strings.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub affiliations: Vec<String>,
    /// Instrument display names; this endpoint sends no CV accessions.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub instruments: Vec<String>,
    /// Software display names.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub softwares: Vec<String>,
    /// Quantification methods, e.g. `"TMT"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub quantification_methods: Vec<String>,
    /// Sample characteristic **values**, e.g. `"liver"`. Which characteristic each value
    /// describes is not recoverable from a hit.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_attributes: Vec<String>,
    /// Organism display names, e.g. `"Homo sapiens (human)"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub organisms: Vec<String>,
    /// Organism-part display names.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub organism_parts: Vec<String>,
    /// Disease display names.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub diseases: Vec<String>,
    /// Publications, each one pre-formatted citation string: no separate PubMed id or DOI.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub references: Vec<String>,
    /// Experiment types, e.g. `"Data-independent acquisition"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub experiment_types: Vec<String>,
    /// File **names** only — a search convenience, not the manifest: no sizes, categories or
    /// locations. Use [`list_files`] or [`list_ftp_files`] to act on the files.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub project_file_names: Vec<String>,
    /// Links to related datasets in other omics repositories.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub other_omics_links: Vec<String>,
    /// Why this project matched: PRIDE field name → matched snippets, the matched terms wrapped in
    /// `<em>` markup. Keys vary per hit and per query, and cross unchanged. The one thing search
    /// returns that project metadata cannot.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub highlights: BTreeMap<String, Vec<String>>,
    /// Downloads per year. Empty when not reported.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub yearly_downloads: Vec<YearlyDownloads>,
    /// Total downloads. 0 means not reported, never a measured zero.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub download_count: u64,
    /// Mean downloads per file. 0 means not reported.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub avg_downloads_per_file: f64,
    /// The project's download-popularity percentile within PRIDE. 0 means not reported.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub percentile: u32,
    /// Downloads attributed to crawlers. 0 means not reported (filled on under half of hits).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub bot_count: u64,
    /// Downloads attributed to institutional or aggregating hubs. 0 means not reported.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub hub_count: u64,
    /// Downloads attributed to ordinary human traffic. 0 means not reported.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub organic_count: u64,
}

impl PrideProjectSearchResult {
    /// Which PRIDE fields the query hit, from [`Self::highlights`], sorted. Empty if PRIDE
    /// reported none.
    #[must_use]
    pub fn matched_fields(&self) -> Vec<&str> {
        self.highlights.keys().map(String::as_str).collect()
    }
}

/// How a search is fetched. By default the page size is mzLib's, and so is the time limit.
#[derive(Debug, Clone, Default)]
pub struct SearchOptions {
    /// How many projects to request per underlying API call. Changes how many requests the fetch
    /// takes, never what comes back. `None` sends none, so mzLib's default applies; a value must
    /// be at least 1.
    pub page_size: Option<u32>,
    /// Time to allow for the whole fetch, every page included.
    ///
    /// `None` (the default) adds no limit of this crate's own, because mzLib already bounds every
    /// request: it retries a transient failure three times, after 5, 20 and 60 s, and each attempt
    /// can take up to 100 s (mzLib#1412). An outage that hangs is therefore reported as
    /// [`MzLibError::ServiceUnavailable`] after about 8 minutes per request. A value cuts the whole
    /// call off after that long with [`MzLibError::Timeout`].
    pub timeout: Option<Duration>,
}

/// A bare calendar date from the bridge, treating anything unreadable as absent.
///
/// Deliberately **not** [`lenient_timestamp`]: the search endpoint sends `"2026-08-15"` with no
/// time and no offset, and widening it to midnight would invent a time PRIDE never reported.
fn lenient_date<'de, D>(deserializer: D) -> std::result::Result<Option<NaiveDate>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let raw: Option<String> = Option::deserialize(deserializer)?;
    Ok(raw.as_deref().and_then(parse_date))
}

fn parse_date(value: &str) -> Option<NaiveDate> {
    let head = value.trim().get(..10)?;
    NaiveDate::parse_from_str(head, "%Y-%m-%d").ok()
}

/// The argv for `pride search`, validated before anything is spawned.
fn build_search_args(keyword: &str, page_size: Option<u32>) -> Result<Vec<String>> {
    let keyword = keyword.trim();
    if keyword.is_empty() {
        return Err(MzLibError::Usage(
            "A search keyword is required, e.g. 'phosphoproteome'.".to_owned(),
        ));
    }
    let keyword = reject_flag_like("keyword", keyword)?;
    let length = keyword.chars().count();
    if length > MAX_KEYWORD_LENGTH {
        return Err(MzLibError::Usage(format!(
            "keyword may be at most {MAX_KEYWORD_LENGTH} characters; got {length}. PRIDE answers a \
             longer keyword with HTTP 500, which cannot be told apart from the service being down."
        )));
    }
    let paging = page_size_args(page_size)?;
    let mut args = vec![
        "pride".to_owned(),
        "search".to_owned(),
        "--keyword".to_owned(),
        keyword,
    ];
    args.extend(paging);
    Ok(args)
}

fn parse_search(data: serde_json::Value) -> Result<Vec<PrideProjectSearchResult>> {
    #[derive(Deserialize)]
    struct Payload {
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        results: Vec<PrideProjectSearchResult>,
    }
    let payload: Payload = serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!("PRIDE search payload could not be read: {error}"))
    })?;
    Ok(payload.results)
}

/// Find PRIDE projects by keyword, with mzLib's default paging.
///
/// **The discovery entry point**: every other function here takes an accession you already have,
/// and this is the one that produces them. See [`search_with`] for the reference.
///
/// # Errors
///
/// As [`search_with`].
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let hits = mzlib::pride::search("plasmodium falciparum schizont")?;
/// for hit in &hits {
///     println!("{}\t{}", hit.accession, hit.title);
/// }
/// # assert_eq!(hits.len(), 6);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
pub fn search(keyword: &str) -> Result<Vec<PrideProjectSearchResult>> {
    search_with(keyword, &SearchOptions::default())
}

/// Find PRIDE Archive projects by keyword, with every page fetched and no accession repeated.
///
/// Paging is handled for you: however many pages the result set spans, you get one list, in
/// PRIDE's ranking order. **An empty list is a real answer** — PRIDE reports no hits as an empty
/// result rather than an error, and unlike [`list_files`] there is no accession here that could
/// have been a typo, so it is not [`MzLibError::ProjectNotFound`]. [`SearchOptions::page_size`]
/// is sent only when set; otherwise mzLib's own default applies. A blank keyword, one longer than
/// [`MAX_KEYWORD_LENGTH`] characters or beginning with `-`, and a `page_size` of zero or larger
/// than the API allows are refused before anything is spawned.
#[doc = include_str!("../docs/reference/pride.search.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::pride::{search_with, SearchOptions};
///
/// let hits = search_with("plasmodium falciparum schizont", &SearchOptions::default())?;
/// assert_eq!(hits[0].accession, "PXD070842");
/// assert_eq!(hits[0].matched_fields(), ["references", "title"]);
/// // A calendar date, not a timestamp: PRIDE sent no time.
/// assert_eq!(hits[0].submission_date.map(|d| d.to_string()).as_deref(), Some("2025-11-17"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/pride.search.see-also.md")]
pub fn search_with(
    keyword: &str,
    options: &SearchOptions,
) -> Result<Vec<PrideProjectSearchResult>> {
    let args = build_search_args(keyword, options.page_size)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    parse_search(data)
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

    const FTP_FIXTURE: &str = include_str!("../tests/fixtures/pride_PXD000001_ftp_files.json");

    fn recorded_ftp() -> Vec<PrideFtpFile> {
        let data: serde_json::Value =
            serde_json::from_str(FTP_FIXTURE).expect("ftp fixture should be valid JSON");
        parse_ftp_files(&data).expect("ftp fixture should parse")
    }

    #[test]
    fn ftp_files_parse_every_file_including_the_rest_hidden_one() {
        let files = recorded_ftp();
        assert_eq!(files.len(), 14);
        // The whole point of the verb: the mzML the REST manifest omits is present.
        assert!(files.iter().any(|f| f.relative_path
            == "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML"));
    }

    #[test]
    fn ftp_nested_file_keeps_its_path_but_a_bare_leaf_name() {
        let nested = recorded_ftp()
            .into_iter()
            .find(|f| f.relative_path.contains('/'))
            .expect("a nested file");
        assert_eq!(
            nested.relative_path,
            "generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz"
        );
        assert_eq!(nested.file_name, "PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz");
        assert_eq!(nested.extension(), ".gz");
        assert!(nested.url.starts_with("https://"));
    }

    #[test]
    fn approximate_total_size_sums_the_complete_listing() {
        let files = recorded_ftp();
        let expected: u64 = files.iter().map(|f| f.approximate_size_bytes).sum();
        assert_eq!(approximate_total_size_bytes(&files), expected);

        let run = files
            .iter()
            .find(|f| f.file_name.ends_with(".raw"))
            .expect("the raw file");
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
        let error = build_list_args("PXD000001", Some(0)).unwrap_err();
        assert!(error.to_string().contains("must be positive"));
    }

    #[test]
    fn a_page_size_larger_than_the_api_allows_is_refused() {
        let error = build_list_args("PXD000001", Some(2_147_483_648)).unwrap_err();
        assert!(error.to_string().contains("larger than the API allows"));
    }

    #[test]
    fn list_args_carry_the_canonical_accession_and_page_size() {
        let args = build_list_args("  pxd000001 ", Some(100)).unwrap();
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
    fn list_files_without_a_page_size_sends_none_so_mzlibs_default_applies() {
        // The crate used to send 100 of its own. mzLib's GetProjectFilesAsync has its default.
        let args = build_list_args("PXD000001", ListOptions::default().page_size).unwrap();
        assert_eq!(args, ["pride", "files", "--accession", "PXD000001"]);
    }

    #[test]
    fn list_defaults_match_pymzlib() {
        let defaults = ListOptions::default();
        assert_eq!(defaults.page_size, None);
        assert_eq!(defaults.timeout, None);
    }

    #[test]
    fn search_defaults_match_pymzlib() {
        let defaults = SearchOptions::default();
        assert_eq!(defaults.page_size, None);
        assert_eq!(defaults.timeout, None);
    }

    // ---------------------------------------------------------- timeouts

    /// A PRIDE call reduced to whether it succeeded, so several fit in one table.
    type PrideCall = fn() -> Result<usize>;

    /// The time limit a caller sets for itself in these tests.
    const OWN_LIMIT: Option<Duration> = Some(Duration::from_secs(60));

    #[test]
    fn pride_metadata_calls_add_no_timeout_of_their_own() {
        // mzLib bounds each PRIDE request itself (#1412: retries after 5, 20 and 60 s, up to 100 s
        // per attempt). The 300 s limit this crate used to add turned a hanging outage into
        // MzLibError::Timeout instead of ServiceUnavailable, as pyMzLib's did (pyMzLib #87).
        let calls: [(&str, PrideCall); 5] = [
            ("list_files", || list_files("PXD000001").map(|f| f.len())),
            ("list_files_with", || {
                list_files_with("PXD000001", &ListOptions::default()).map(|f| f.len())
            }),
            ("list_ftp_files", || {
                list_ftp_files("PXD000001").map(|f| f.len())
            }),
            ("search", || search("phospho").map(|h| h.len())),
            ("search_with", || {
                search_with("phospho", &SearchOptions::default()).map(|h| h.len())
            }),
        ];
        for (name, call) in calls {
            let mut outcome = None;
            let seen = bridge::intercept::timeouts(|| outcome = Some(call()));
            assert_eq!(seen, [None], "{name} should leave the time limit to mzLib");
            assert!(
                matches!(outcome, Some(Err(MzLibError::ServiceUnavailable { .. }))),
                "{name}: an outage should stay an outage, got {outcome:?}"
            );
        }
    }

    #[test]
    fn a_callers_own_pride_timeout_still_reaches_the_bridge() {
        let calls: [(&str, PrideCall); 3] = [
            ("list_files_with", || {
                let options = ListOptions {
                    timeout: OWN_LIMIT,
                    ..Default::default()
                };
                list_files_with("PXD000001", &options).map(|f| f.len())
            }),
            ("list_ftp_files_with", || {
                list_ftp_files_with("PXD000001", OWN_LIMIT).map(|f| f.len())
            }),
            ("search_with", || {
                let options = SearchOptions {
                    timeout: OWN_LIMIT,
                    ..Default::default()
                };
                search_with("phospho", &options).map(|h| h.len())
            }),
        ];
        for (name, call) in calls {
            let seen = bridge::intercept::timeouts(|| {
                let _ = call();
            });
            assert_eq!(
                seen,
                [OWN_LIMIT],
                "{name} should pass the caller's timeout on"
            );
        }
    }

    // ---------------------------------------------------------- checksums (mzLib #1427)

    /// PRIDE's checksum list for PXD000001, recorded live through the bridge, shared with pyMzLib.
    const CHECKSUMS_FIXTURE: &str =
        include_str!("../tests/fixtures/pride_PXD000001_checksums.json");

    fn recorded_checksums() -> BTreeMap<String, PrideFileChecksum> {
        let data: serde_json::Value =
            serde_json::from_str(CHECKSUMS_FIXTURE).expect("checksums fixture should be JSON");
        parse_checksums(&data, "PXD000001").expect("checksums fixture should parse")
    }

    #[test]
    fn file_checksums_keys_every_row_by_its_file_name() {
        let data: serde_json::Value = serde_json::from_str(CHECKSUMS_FIXTURE).unwrap();
        let checksums = recorded_checksums();
        assert_eq!(checksums.len(), 12);
        assert_eq!(data["file_count"], 12);
        assert!(checksums.iter().all(|(name, row)| *name == row.file_name));
        // The wire's order (sorted by name, ordinal) is the map's order, so iterating is
        // deterministic and matches what pyMzLib's dict keeps.
        let wire: Vec<&str> = data["files"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| row["file_name"].as_str().unwrap())
            .collect();
        assert_eq!(
            checksums.keys().map(String::as_str).collect::<Vec<_>>(),
            wire
        );
    }

    #[test]
    fn file_checksum_fields_are_typed_and_carry_the_project() {
        let checksums = recorded_checksums();
        let raw = &checksums["TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw"];
        assert_eq!(raw.md5, "c567cecbb9866796e690a4958746f7b6");
        assert_eq!(raw.size_bytes, 220_475_548);
        assert_eq!(raw.project_accession, "PXD000001");
    }

    #[test]
    fn the_generated_files_have_no_checksum_row() {
        // The list covers deposits only: the FTP tree's generated/ files have no row (mzLib #1427).
        let checksums = recorded_checksums();
        let unchecked: Vec<String> = recorded_ftp()
            .into_iter()
            .map(|f| f.relative_path)
            .filter(|path| !checksums.contains_key(path))
            .collect();
        assert_eq!(
            unchecked,
            [
                "generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz",
                "generated/PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz",
            ]
        );
    }

    #[test]
    fn the_checksum_size_is_the_download_length_not_the_manifests() {
        // PRIDE's REST manifest over-reports this .xml.gz; the checksum list has the length a
        // download comes to (an HTTP HEAD's Content-Length on 2026-10-21).
        let name = "PRIDE_Exp_Complete_Ac_22134.xml.gz";
        let rest = recorded_files()
            .into_iter()
            .find(|f| f.file_name == name)
            .unwrap();
        assert_eq!(
            (rest.file_size_bytes, recorded_checksums()[name].size_bytes),
            (10_677_205, 10_668_000)
        );
    }

    #[test]
    fn file_checksums_sends_the_verb_and_the_normalised_accession() {
        // PRIDE's checksum route is case-sensitive: pxd000001 is an unknown accession there.
        let canonical = normalise_accession("  pxd000001 ").unwrap();
        assert_eq!(
            checksums_args(&canonical),
            ["pride", "checksums", "--accession", "PXD000001"],
            "nothing but the accession is sent"
        );
    }

    #[test]
    fn an_existing_project_with_no_checksums_is_an_empty_map_not_an_error() {
        // mzLib has already asked whether the project exists, so an empty list is an answer.
        let empty = serde_json::json!({"accession": "PXD000001", "file_count": 0, "files": []});
        assert!(parse_checksums(&empty, "PXD000001").unwrap().is_empty());
    }

    #[test]
    fn a_row_without_its_md5_fails_rather_than_defaulting() {
        // Defaulting it to "" would make every later check of that file blame the download.
        let broken = serde_json::json!({
            "accession": "PXD000001", "file_count": 1,
            "files": [{"file_name": "run1.raw", "size_bytes": 10}],
        });
        let error = parse_checksums(&broken, "PXD000001").unwrap_err();
        assert!(matches!(error, MzLibError::Protocol(_)), "{error:?}");
    }

    #[test]
    fn file_checksums_unknown_project_is_project_not_found() {
        let error = checksums_error(
            MzLibError::Bridge {
                error_type: "MzLibException".to_owned(),
                message: "PRIDE Archive has no project with accession 'PXD999999999'.".to_owned(),
            },
            "PXD999999999",
        );
        assert!(
            matches!(&error, MzLibError::ProjectNotFound(m) if m.contains("no public project 'PXD999999999'"))
        );
    }

    #[test]
    fn a_list_mzlib_cannot_read_is_not_reported_as_a_typo() {
        // Also an MzLibException, but the project exists and its list is broken. "Check for a
        // typo" would send the caller looking for the wrong mistake.
        for message in [
            "PRIDE Archive's checksum list for 'PXD000001' starts with 'File-Name,MD5', not the \
             expected 'File-Name,File-MD5Checksum,File-Size'.",
            "PRIDE Archive's checksum list for 'PXD000001' lists 'run1.raw' twice.",
        ] {
            let error = checksums_error(
                MzLibError::Bridge {
                    error_type: "MzLibException".to_owned(),
                    message: message.to_owned(),
                },
                "PXD000001",
            );
            assert!(
                matches!(&error, MzLibError::Bridge { error_type, .. } if error_type == "MzLibException"),
                "{error:?}"
            );
        }
    }

    #[test]
    fn a_checksum_outage_stays_service_unavailable_and_adds_no_timeout() {
        let mut outcome = None;
        let seen =
            bridge::intercept::timeouts(|| outcome = Some(fetch_checksums("PXD000001", None)));
        assert_eq!(
            seen,
            [None],
            "file_checksums should leave the time limit to mzLib"
        );
        assert!(matches!(
            outcome,
            Some(Err(MzLibError::ServiceUnavailable { .. }))
        ));
        let seen = bridge::intercept::timeouts(|| {
            let _ = fetch_checksums("PXD000001", OWN_LIMIT);
        });
        assert_eq!(seen, [OWN_LIMIT], "a caller's own limit reaches the bridge");
    }

    #[test]
    fn file_checksums_rejects_a_bad_accession_before_any_work() {
        for accession in ["", "   ", "not-an-accession"] {
            let seen = bridge::intercept::timeouts(|| {
                let error = file_checksums(accession).unwrap_err();
                assert!(matches!(error, MzLibError::Usage(_)), "{accession:?}");
            });
            assert!(seen.is_empty(), "{accession:?} reached the bridge");
        }
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

    // ---------------------------------------------------------- search

    const SEARCH_FIXTURE: &str = include_str!("../tests/fixtures/pride_search_plasmodium.json");

    fn recorded_hits() -> Vec<PrideProjectSearchResult> {
        parse_search(serde_json::from_str(SEARCH_FIXTURE).unwrap()).unwrap()
    }

    #[test]
    fn search_hits_parse_with_their_flattened_fields() {
        let hits = recorded_hits();
        assert_eq!(hits.len(), 6);
        let first = &hits[0];
        assert_eq!(first.accession, "PXD070842");
        assert_eq!(first.instruments, ["Orbitrap Fusion Lumos"]);
        assert_eq!(first.matched_fields(), ["references", "title"]);
        assert!(first.highlights["title"][0].contains("<em>"));
    }

    #[test]
    fn search_dates_are_calendar_dates_not_midnights() {
        let first = &recorded_hits()[0];
        assert_eq!(first.submission_date, NaiveDate::from_ymd_opt(2025, 11, 17));
        assert_eq!(parse_date(""), None);
        assert_eq!(parse_date("not a date"), None);
        assert_eq!(
            parse_date("2026-08-15T00:00:00"),
            NaiveDate::from_ymd_opt(2026, 8, 15)
        );
    }

    #[test]
    fn search_passes_empty_keywords_through_rather_than_filtering() {
        // PRIDE ships blank keywords on about 9% of hits; every binding reports them as sent.
        assert_eq!(recorded_hits()[0].keywords, [""]);
    }

    #[test]
    fn a_zero_download_count_is_not_reported_rather_than_measured() {
        // The fixture's first hit has no counts; they arrive as 0, never null.
        let first = &recorded_hits()[0];
        assert_eq!(first.download_count, 0);
        assert_eq!(first.avg_downloads_per_file, 0.0);
    }

    #[test]
    fn an_empty_search_is_an_empty_list_not_an_error() {
        let hits =
            parse_search(serde_json::json!({"keyword": "x", "result_count": 0, "results": []}))
                .unwrap();
        assert!(hits.is_empty());
    }

    #[test]
    fn search_args_are_assembled_in_the_documented_order() {
        assert_eq!(
            build_search_args("  phosphoproteome ", Some(50)).unwrap(),
            [
                "pride",
                "search",
                "--keyword",
                "phosphoproteome",
                "--page-size",
                "50"
            ]
        );
    }

    #[test]
    fn search_without_a_page_size_sends_none_so_mzlibs_default_applies() {
        // The crate used to send 100 of its own. mzLib's SearchProjectsAsync has its default.
        let args = build_search_args("phospho", SearchOptions::default().page_size).unwrap();
        assert_eq!(args, ["pride", "search", "--keyword", "phospho"]);
    }

    #[test]
    fn a_bad_search_is_refused_before_anything_is_spawned() {
        for (keyword, page_size) in [
            ("", None),
            ("   ", Some(100)),
            ("-x", None),
            ("--no-overwrite", Some(100)),
            ("ok", Some(0)),
            ("ok", Some(u32::MAX)),
        ] {
            let error = build_search_args(keyword, page_size).unwrap_err();
            assert!(
                matches!(error, MzLibError::Usage(_)),
                "{keyword:?} {page_size:?}"
            );
        }
        let long = "a".repeat(MAX_KEYWORD_LENGTH + 1);
        let error = build_search_args(&long, None).unwrap_err();
        assert!(error.to_string().contains("HTTP 500"), "{error}");
        assert!(build_search_args(&"a".repeat(MAX_KEYWORD_LENGTH), None).is_ok());
    }
}
