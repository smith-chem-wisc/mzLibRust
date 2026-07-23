//! Live canaries against the real UniProt, through the real bridge.
//!
//! These are the tests that would catch mzLib or UniProt changing under us. They **skip** rather
//! than fail when UniProt is unavailable.
//!
//! Run with `cargo test --features live`. The two histone tests are genuinely slow (modification
//! isoforms are enumerated combinatorially) and are marked `#[ignore]`; run them with
//! `cargo test --features live -- --ignored`.

#![cfg(feature = "live")]

mod support;

use mzlib::peptidoform::{fragments, fragments_with, FragmentOptions};
use support::{external_service, require_bridge};

/// Human serum albumin: large, heavily annotated, and mostly annotated with glycosylation sites
/// that have no defined mass — which is what makes it the right protein for the census.
const ALBUMIN: &str = "P02768";

/// Histone H3.1: where modification-isoform combinatorics actually bite.
const HISTONE: &str = "P68431";

/// The hydrogen atom, in daltons. `c_i + z•_(n-i)` closes on the peptide mass plus one of these.
const HYDROGEN_MASS: f64 = 1.007_825_032;

fn bare(max_modifications: u32, min_length: u32) -> FragmentOptions {
    FragmentOptions {
        max_modifications,
        min_length,
        ..Default::default()
    }
}

#[test]
fn the_workflow_still_answers_end_to_end() {
    let Some(()) = require_bridge() else { return };

    let Some(digest) = external_service("UniProt", fragments(ALBUMIN)) else {
        return;
    };

    assert_eq!(digest.accession, ALBUMIN);
    assert_eq!(digest.sequence_length, 609, "albumin's length changed");
    assert!(!digest.peptides.is_empty());
    assert!(digest.fragment_count() > 0);
}

#[test]
fn etd_produces_c_and_z_ions() {
    // The dissociation type must reach mzLib, or the caller silently gets the wrong chemistry.
    // On the y ions ETD also emits, see smith-chem-wisc/mzLib#1109 and the offline test that pins
    // it; this canary asserts only what ETD genuinely should produce.
    let Some(()) = require_bridge() else { return };

    let Some(digest) = external_service("UniProt", fragments_with(ALBUMIN, &bare(0, 20))) else {
        return;
    };

    let kinds: std::collections::HashSet<String> = digest
        .peptides
        .iter()
        .flat_map(|p| p.fragments.iter())
        .map(|f| f.product_type.clone())
        .collect();

    assert!(
        kinds.iter().any(|k| k.starts_with('c')),
        "expected c ions from ETD, got {kinds:?}"
    );
    assert!(
        kinds.iter().any(|k| k.to_lowercase().starts_with('z')),
        "expected z ions from ETD, got {kinds:?}"
    );
}

#[test]
fn fragment_series_close_on_the_precursor_mass() {
    // c_i + z•_(n-i) must equal the peptide mass plus one hydrogen, for every i.
    //
    // This is the invariant that catches a whole class of silent error: a wrong terminal group, a
    // wrong ion definition, a modification applied to the wrong terminus. Any of those produce a
    // fragment table with sensible spacings and monotonic series that is nonetheless wrong.
    let Some(()) = require_bridge() else { return };

    let Some(digest) = external_service("UniProt", fragments_with(ALBUMIN, &bare(0, 15))) else {
        return;
    };

    let mut checked = 0_u32;
    for peptide in digest.peptides.iter().take(5) {
        let c_ions: std::collections::HashMap<i32, f64> = peptide
            .fragments
            .iter()
            .filter(|f| f.product_type == "c" && f.neutral_loss == 0.0)
            .map(|f| (f.fragment_number, f.neutral_mass))
            .collect();
        let z_ions: std::collections::HashMap<i32, f64> = peptide
            .fragments
            .iter()
            .filter(|f| f.product_type.to_lowercase().starts_with('z') && f.neutral_loss == 0.0)
            .map(|f| (f.fragment_number, f.neutral_mass))
            .collect();

        for (&index, &c_mass) in &c_ions {
            let Some(&z_mass) = z_ions.get(&(peptide.length - index)) else {
                continue;
            };
            let closure = c_mass + z_mass - peptide.monoisotopic_mass;
            assert!(
                (closure - HYDROGEN_MASS).abs() < 5e-4,
                "{}: c{index} + z{} - M = {closure:.6}, expected {HYDROGEN_MASS:.6}",
                peptide.base_sequence,
                peptide.length - index
            );
            checked += 1;
        }
    }

    assert!(
        checked > 10,
        "expected many closure pairs to check, got {checked}"
    );
}

#[test]
fn the_annotation_census_reports_what_was_excluded() {
    // Albumin's annotations are mostly glycosylation sites, which have no defined mass.
    let Some(()) = require_bridge() else { return };

    let Some(digest) = external_service("UniProt", fragments_with(ALBUMIN, &bare(0, 7))) else {
        return;
    };
    let census = &digest.modification_census;

    assert!(census.annotated > census.applied);
    assert!(census.excluded() > 0);
    assert!(
        census.explain().contains("glycosylation"),
        "{}",
        census.explain()
    );
}

#[test]
fn modifications_change_the_answer_substantially() {
    // The control that shows the annotations are doing real work.
    let Some(()) = require_bridge() else { return };

    let Some(with_mods) = external_service("UniProt", fragments_with(ALBUMIN, &bare(1, 7))) else {
        return;
    };
    let without_options = FragmentOptions {
        modifications: false,
        ..Default::default()
    };
    let Some(without) = external_service("UniProt", fragments_with(ALBUMIN, &without_options))
    else {
        return;
    };

    assert!(with_mods.peptides.len() > without.peptides.len());
    assert!(!with_mods.modified_peptides().is_empty());
    assert!(without.modified_peptides().is_empty());
}

#[test]
#[ignore = "slow: histone modification isoforms are enumerated combinatorially"]
fn modification_isoforms_are_enumerated_combinatorially() {
    // Histones are where this matters: alternatives at one residue multiply across residues.
    let Some(()) = require_bridge() else { return };

    let Some(one) = external_service("UniProt", fragments_with(HISTONE, &bare(1, 7))) else {
        return;
    };
    let Some(two) = external_service("UniProt", fragments_with(HISTONE, &bare(2, 7))) else {
        return;
    };

    assert!(
        two.peptides.len() > one.peptides.len() * 2,
        "modification isoforms should multiply, not add: {} vs {}",
        two.peptides.len(),
        one.peptides.len()
    );
}

#[test]
#[ignore = "slow: raising the isoform cap on a histone enumerates tens of thousands of forms"]
fn the_isoform_cap_truncates_and_says_so() {
    // mzLib's default of 1024 isoforms per peptide truncates silently. It must not be silent here:
    // a truncated peptidoform list is indistinguishable from a short one.
    let Some(()) = require_bridge() else { return };

    let capped_options = FragmentOptions {
        max_modifications: 4,
        max_isoforms: 1024,
        timeout: None,
        ..Default::default()
    };
    let raised_options = FragmentOptions {
        max_modifications: 4,
        max_isoforms: 100_000,
        timeout: None,
        ..Default::default()
    };

    let Some(capped) = external_service("UniProt", fragments_with(HISTONE, &capped_options)) else {
        return;
    };
    let Some(raised) = external_service("UniProt", fragments_with(HISTONE, &raised_options)) else {
        return;
    };

    assert!(
        capped.truncated(),
        "the default cap binds on a histone at four modifications"
    );
    assert!(capped.peptides_at_cap > 0);
    assert!(!raised.truncated());
    assert!(
        raised.peptides.len() > capped.peptides.len(),
        "raising the cap must recover peptidoforms the default discarded"
    );
}

#[test]
fn an_unknown_accession_is_a_usage_error_not_an_empty_result() {
    let Some(()) = require_bridge() else { return };

    // Well-formed by UniProt's grammar, but no such entry — so it reaches UniProt and comes back
    // as a 404, which must surface as a usage error rather than an outage or an empty digest.
    match fragments("Q6ZZZ9") {
        Err(mzlib::MzLibError::Usage(message)) => {
            assert!(message.contains("Q6ZZZ9"), "{message}");
        }
        Err(mzlib::MzLibError::ServiceUnavailable { message, .. }) => {
            support::skip(&format!("UniProt unavailable ({message})"));
        }
        Err(other) => panic!("expected a usage error, got {other:?}"),
        Ok(digest) => panic!(
            "an unknown accession returned {} peptides",
            digest.peptides.len()
        ),
    }
}

#[test]
fn an_unknown_protease_names_the_alternatives() {
    // A rejection that does not say what IS allowed sends the user to the source.
    let Some(()) = require_bridge() else { return };

    let options = FragmentOptions {
        protease: "definitely-not-a-protease".to_owned(),
        ..Default::default()
    };

    match fragments_with(ALBUMIN, &options) {
        Err(mzlib::MzLibError::Usage(message)) => {
            assert!(message.contains("Unknown protease"), "{message}");
            assert!(message.contains("trypsin"), "{message}");
        }
        Err(mzlib::MzLibError::ServiceUnavailable { message, .. }) => {
            support::skip(&format!("UniProt unavailable ({message})"));
        }
        Err(other) => panic!("expected a usage error, got {other:?}"),
        Ok(_) => panic!("an unknown protease was accepted"),
    }
}
