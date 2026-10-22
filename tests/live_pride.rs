//! Live canaries against the real PRIDE Archive, through the real bridge.
//!
//! These are the tests that would catch mzLib or PRIDE changing under us — the offline suite only
//! ever sees a recorded manifest. They **skip** rather than fail when EBI is unavailable, because
//! a red build that might just mean "EBI is having a bad morning" is an ambiguous red build, and
//! an ambiguous red build gets ignored.
//!
//! Run with `cargo test --features live`.

#![cfg(feature = "live")]

mod support;

use support::{external_service, require_bridge};

/// PXD000001 is the canonical PRIDE example project: small, stable, and public since 2012.
const PROJECT: &str = "PXD000001";

#[test]
fn the_api_still_answers_and_the_manifest_still_parses() {
    let Some(()) = require_bridge() else { return };

    let Some(files) = external_service("PRIDE Archive", mzlib::pride::list_files(PROJECT)) else {
        return;
    };

    assert!(!files.is_empty(), "PXD000001 has published files");
    assert!(
        mzlib::pride::total_size_bytes(&files) > 0,
        "the manifest should carry real sizes"
    );
}

#[test]
fn the_ftp_listing_is_more_complete_than_the_rest_manifest() {
    let Some(()) = require_bridge() else { return };

    let Some(rest) = external_service("PRIDE Archive", mzlib::pride::list_files(PROJECT)) else {
        return;
    };
    let Some(ftp) = external_service("PRIDE Archive", mzlib::pride::list_ftp_files(PROJECT)) else {
        return;
    };

    // The whole reason list_ftp_files exists. A comparison, not a hard count, so it survives PRIDE
    // re-curating the project - while still catching the walk finding nothing, or the REST manifest
    // quietly becoming complete (which would leave the extra verb without a purpose).
    assert!(
        !ftp.is_empty(),
        "the FTP walk returned nothing for a project known to have files"
    );
    assert!(
        ftp.len() > rest.len(),
        "the FTP listing ({}) should exceed the REST manifest ({})",
        ftp.len(),
        rest.len()
    );

    let first = &ftp[0];
    assert!(!first.relative_path.is_empty(), "relative_path");
    assert!(
        first.url.starts_with("https://"),
        "url must be the HTTPS location"
    );

    // Not FTP-total > REST-total: PRIDE over-reports decompressed .gz sizes, so that inequality is
    // not durable. What is: the total is positive and at least its largest single file.
    let largest = ftp
        .iter()
        .map(|f| f.approximate_size_bytes)
        .max()
        .unwrap_or(0);
    assert!(
        largest > 0 && mzlib::pride::approximate_total_size_bytes(&ftp) >= largest,
        "the approximate total must be positive and cover its largest file"
    );
}

#[test]
fn the_fields_the_rust_layer_reads_are_still_populated() {
    let Some(()) = require_bridge() else { return };

    let Some(files) = external_service("PRIDE Archive", mzlib::pride::list_files(PROJECT)) else {
        return;
    };

    let first = &files[0];
    assert!(!first.file_name.is_empty(), "file_name");
    assert!(first.file_size_bytes > 0, "file_size_bytes");
    assert!(!first.category.is_empty(), "category");
    assert_eq!(first.project_accession, PROJECT);
    assert!(
        files.iter().any(|file| file.submission_date.is_some()),
        "at least one file should carry a parseable timestamp"
    );
}

#[test]
fn at_least_one_file_is_still_reachable_over_https() {
    // If PRIDE ever stopped publishing HTTPS locations, `download` would break for everyone and
    // the offline suite would never notice.
    let Some(()) = require_bridge() else { return };

    let Some(files) = external_service("PRIDE Archive", mzlib::pride::list_files(PROJECT)) else {
        return;
    };

    assert!(
        files.iter().any(mzlib::pride::PrideFile::downloadable),
        "no file in {PROJECT} has an HTTPS location any more"
    );
}

#[test]
fn an_unknown_accession_raises_rather_than_reporting_an_empty_project() {
    // PRIDE answers an unknown accession with an empty result, not a 404. A binding that passes
    // that through lets a typo produce "0 files, done" and a green exit.
    let Some(()) = require_bridge() else { return };

    match mzlib::pride::list_files("PXD999999999") {
        Err(mzlib::MzLibError::ProjectNotFound(_)) => {}
        Err(mzlib::MzLibError::ServiceUnavailable { message, .. }) => {
            support::skip(&format!("PRIDE Archive unavailable ({message})"));
        }
        Err(other) => panic!("expected ProjectNotFound, got {other:?}"),
        Ok(files) => panic!("an unknown accession returned {} files", files.len()),
    }
}

#[test]
fn a_real_selection_can_be_downloaded_directly() {
    let Some(()) = require_bridge() else { return };

    let Some(files) = external_service("PRIDE Archive", mzlib::pride::list_files(PROJECT)) else {
        return;
    };

    // The smallest downloadable file, so the canary costs seconds rather than gigabytes.
    let mut candidates: Vec<_> = files
        .iter()
        .filter(|file| file.downloadable())
        .cloned()
        .collect();
    candidates.sort_by_key(|file| file.file_size_bytes);
    let Some(smallest) = candidates.first().cloned() else {
        support::skip("no downloadable file in the manifest");
        return;
    };

    let destination = std::env::temp_dir().join("mzlibrust-live-pride-selection");
    let _ = std::fs::remove_dir_all(&destination);

    let Some(written) = external_service(
        "PRIDE Archive",
        mzlib::pride::download_files(
            std::slice::from_ref(&smallest),
            &destination,
            &Default::default(),
        ),
    ) else {
        return;
    };

    assert_eq!(written.len(), 1, "exactly the selected file");
    let landed = &written[0];
    assert!(landed.is_file(), "{} should exist", landed.display());
    assert_eq!(
        std::fs::metadata(landed).unwrap().len(),
        smallest.file_size_bytes,
        "the downloaded file should be the size the manifest promised"
    );

    let _ = std::fs::remove_dir_all(&destination);
}

#[test]
fn a_real_download_still_works_end_to_end() {
    // The filtered path, as opposed to the explicit-selection path above. Both reach the same
    // bridge verb by different argument routes, and only this one exercises the filter.
    let Some(()) = require_bridge() else { return };

    let Some(files) = external_service("PRIDE Archive", mzlib::pride::list_files(PROJECT)) else {
        return;
    };

    // Pick a category that exists and whose files are small, so the canary stays cheap.
    let Some(smallest) = files
        .iter()
        .filter(|file| file.downloadable())
        .min_by_key(|file| file.file_size_bytes)
    else {
        support::skip("no downloadable file in the manifest");
        return;
    };
    let extension = smallest.extension();
    if extension.is_empty() {
        support::skip("smallest file has no extension to filter on");
        return;
    }

    let destination = std::env::temp_dir().join("mzlibrust-live-pride-download");
    let _ = std::fs::remove_dir_all(&destination);

    let options = mzlib::pride::DownloadOptions {
        extensions: vec![extension.clone()],
        ..Default::default()
    };
    let Some(written) = external_service(
        "PRIDE Archive",
        mzlib::pride::download(PROJECT, &destination, &options),
    ) else {
        return;
    };

    assert!(
        !written.is_empty(),
        "filtering on '{extension}' matched nothing, but the manifest said otherwise"
    );
    assert!(written.iter().all(|path| path.is_file()));

    let _ = std::fs::remove_dir_all(&destination);
}

#[test]
fn a_filter_that_matches_nothing_is_an_error_not_a_green_no_op() {
    // The doctrine the offline suite asserts, proven against the real repository: asking for a
    // filter and getting silence must not look like success.
    let Some(()) = require_bridge() else { return };

    let destination = std::env::temp_dir().join("mzlibrust-live-pride-nomatch");
    let options = mzlib::pride::DownloadOptions {
        extensions: vec![".definitely-not-a-real-extension".to_owned()],
        ..Default::default()
    };

    match mzlib::pride::download(PROJECT, &destination, &options) {
        Err(mzlib::MzLibError::Usage(message)) => {
            assert!(message.contains("matched"), "{message}");
        }
        Err(mzlib::MzLibError::ServiceUnavailable { message, .. }) => {
            support::skip(&format!("PRIDE Archive unavailable ({message})"));
        }
        Err(other) => panic!("expected a usage error, got {other:?}"),
        Ok(paths) => panic!("a nonsense filter wrote {} files", paths.len()),
    }

    let _ = std::fs::remove_dir_all(&destination);
}

#[test]
fn search_still_finds_the_first_project_and_says_why() {
    // pride search was already on the wire in the pinned bridge; this canary is what notices PRIDE
    // changing its search projection under the typed hit.
    let Some(()) = require_bridge() else { return };

    let Some(hits) = external_service("PRIDE Archive", mzlib::pride::search(PROJECT)) else {
        return;
    };

    let first = hits
        .iter()
        .find(|hit| hit.accession == PROJECT)
        .expect("searching for PXD000001 should find PXD000001");
    assert!(!first.title.is_empty(), "a hit carries its title");
    assert!(
        first.submission_date.is_some(),
        "the first public project has a submission date"
    );
    let mut accessions: Vec<&str> = hits.iter().map(|hit| hit.accession.as_str()).collect();
    accessions.sort_unstable();
    accessions.dedup();
    assert_eq!(accessions.len(), hits.len(), "no accession is repeated");
}

// ------------------------------------------------------------------- checksums (mzLib #1427)

#[test]
fn the_checksum_list_still_parses_and_covers_the_deposited_files() {
    // Asserted as relations to the FTP tree rather than counts, so it survives PRIDE re-curating
    // the project: every row names a file at the project's root, and none names one of PRIDE's
    // generated/ files.
    let Some(()) = support::require_verb("pride checksums") else {
        return;
    };
    let Some(checksums) = external_service("PRIDE Archive", mzlib::pride::file_checksums(PROJECT))
    else {
        return;
    };
    let Some(ftp) = external_service("PRIDE Archive", mzlib::pride::list_ftp_files(PROJECT)) else {
        return;
    };

    assert!(
        !checksums.is_empty(),
        "PRIDE has always listed checksums for the canary; an empty list means the route changed"
    );
    let root: Vec<&str> = ftp
        .iter()
        .map(|f| f.relative_path.as_str())
        .filter(|path| !path.contains('/'))
        .collect();
    let outside: Vec<&String> = checksums
        .keys()
        .filter(|name| !root.contains(&name.as_str()))
        .collect();
    assert!(outside.is_empty(), "rows outside the FTP root: {outside:?}");
    assert!(
        ftp.iter()
            .any(|f| f.relative_path.starts_with("generated/")),
        "the canary should still hold generated/ files, which have no checksum row"
    );
    assert!(checksums
        .values()
        .all(|c| c.md5.len() == 32 && c.size_bytes > 0 && c.project_accession == PROJECT));
}

#[test]
fn a_checksum_row_matches_the_bytes_a_download_delivers() {
    // The point of the list: its size and MD5 are those of the file a download writes. A small
    // file in both the REST manifest and the checksum list keeps this fast.
    use md5::{Digest, Md5};

    let Some(()) = support::require_verb("pride checksums") else {
        return;
    };
    let Some(checksums) = external_service("PRIDE Archive", mzlib::pride::file_checksums(PROJECT))
    else {
        return;
    };
    let Some(files) = external_service("PRIDE Archive", mzlib::pride::list_files(PROJECT)) else {
        return;
    };
    let Some(chosen) = files
        .iter()
        .filter(|f| f.downloadable() && checksums.contains_key(&f.file_name))
        .min_by_key(|f| checksums[&f.file_name].size_bytes)
        .cloned()
    else {
        support::skip("no downloadable file of the manifest has a checksum row");
        return;
    };

    let destination = std::env::temp_dir().join("mzlibrust-live-pride-checksum");
    let _ = std::fs::remove_dir_all(&destination);
    let Some(written) = external_service(
        "PRIDE Archive",
        mzlib::pride::download_files(
            std::slice::from_ref(&chosen),
            &destination,
            &Default::default(),
        ),
    ) else {
        return;
    };

    let expected = &checksums[&chosen.file_name];
    let bytes = std::fs::read(&written[0]).unwrap();
    assert_eq!(
        bytes.len() as u64,
        expected.size_bytes,
        "{}",
        chosen.file_name
    );
    let mut md5 = Md5::new();
    md5.update(&bytes);
    assert_eq!(format!("{:x}", md5.finalize()), expected.md5.to_lowercase());

    let _ = std::fs::remove_dir_all(&destination);
}

#[test]
fn an_unknown_accession_raises_rather_than_reporting_no_checksums() {
    // PRIDE answers an unknown accession with an empty list, as it does a project with no
    // checksums. mzLib asks for the project to tell them apart, so a typo is not "0 files to check".
    let Some(()) = support::require_verb("pride checksums") else {
        return;
    };
    match mzlib::pride::file_checksums("PXD999999999") {
        Err(mzlib::MzLibError::ProjectNotFound(_)) => {}
        Err(mzlib::MzLibError::ServiceUnavailable { message, .. }) => {
            support::skip(&format!("PRIDE Archive unavailable ({message})"));
        }
        Err(other) => panic!("expected ProjectNotFound, got {other:?}"),
        Ok(rows) => panic!("an unknown accession returned {} checksums", rows.len()),
    }
}
