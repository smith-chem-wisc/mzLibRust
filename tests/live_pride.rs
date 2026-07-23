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
