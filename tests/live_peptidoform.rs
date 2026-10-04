//! Live canaries against the real UniProt, through the real bridge, and the sequence conversion
//! checked against the live bridge.
//!
//! These are the tests that would catch mzLib or UniProt changing under us. They **skip** rather
//! than fail when UniProt is unavailable. The `convert` tests need no network, only a bridge that
//! dispatches `peptidoform convert` (pyMzLib 0.4.0's, the pinned one, or later), and **skip** with
//! an older one.
//!
//! Run with `cargo test --features live`. The two histone tests are genuinely slow (modification
//! isoforms are enumerated combinatorially) and are marked `#[ignore]`; run them with
//! `cargo test --features live -- --ignored`.

#![cfg(feature = "live")]

mod support;

use mzlib::peptidoform::{
    convert, convert_with, fragments, fragments_with, ConversionMode, ConvertOptions,
    FragmentOptions, SequenceConversions,
};
use support::{external_service, require_bridge, require_verb};

/// Human serum albumin: large, heavily annotated, and mostly annotated with glycosylation sites,
/// which mzLib excludes on feature type — which is what makes it the right protein for the census.
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
fn the_documented_ground_truth_digest_counts_reproduce() {
    // Pins the bake-off ground truth (design/bakeoff/DESIGN.md, Task 2) against the live bridge so a
    // documented count cannot quietly rot again — the exact failure smith-chem-wisc/mzLibRust#1
    // caught, where the table read 254/269 while the real digest yields 243/257. Distinct base
    // sequences, ETD, both termini, 2 missed cleavages, up to 2 modifications applied — the defaults.
    let Some(()) = require_bridge() else { return };

    // (protease, min_length, expected distinct base sequences)
    let cases = [
        ("trypsin|P", 7, 195),
        ("trypsin|P", 1, 243),
        ("trypsin", 7, 202),
        ("trypsin", 1, 257),
    ];
    for (protease, min_length, expected) in cases {
        let options = FragmentOptions {
            protease: protease.to_owned(),
            min_length,
            ..Default::default()
        };
        let Some(digest) = external_service("UniProt", fragments_with(ALBUMIN, &options)) else {
            return;
        };
        assert_eq!(
            digest.distinct_base_sequences(),
            expected,
            "{protease} min_length={min_length}: albumin ground truth is {expected} distinct base \
             sequences"
        );
    }
}

#[test]
fn etd_produces_c_and_z_ions() {
    // The dissociation type must reach mzLib, or the caller silently gets the wrong chemistry.
    // ETD used to emit y ions as well (smith-chem-wisc/mzLib#1109, fixed by #1114); the offline
    // etd_produces_c_and_zdot_but_no_y_ions pins their absence. This canary asserts only what ETD
    // genuinely should produce.
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
    // Albumin's annotations are mostly glycosylation sites, which mzLib excludes on feature type.
    let Some(()) = require_bridge() else { return };

    let Some(digest) = external_service("UniProt", fragments_with(ALBUMIN, &bare(0, 7))) else {
        return;
    };
    let census = &digest.modification_census;

    assert!(census.annotated > census.applied);
    assert!(census.excluded() > 0);
    // The exclusion is reported at feature-type granularity, not modification-name level.
    assert!(
        census.explain().contains("glycosylation"),
        "{}",
        census.explain()
    );
    assert!(
        census.explain().contains("feature type"),
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

// ------------------------------------------------------------------ convert

/// The four probe sequences pyMzLib's tests and the bridge's C# tests assert.
const PROBE: [&str; 4] = [
    "[UniProt:N-acetylserine on S]SEQK",
    "PEPK[UniProt:N6,N6-dimethyllysine on K]R",
    "PEPM[Common Variable:Oxidation on M]K",
    "PEPK[Made Up:Not a modification on K]R",
];

fn recorded_conversion(name: &str) -> SequenceConversions {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name);
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// Everything a recording asserts that the live bridge must still say: the envelope and every row.
fn assert_same_conversion(live: &SequenceConversions, recorded: &SequenceConversions) {
    assert_eq!(live.source_format, recorded.source_format);
    assert_eq!(live.target_format, recorded.target_format);
    assert_eq!(live.mode, recorded.mode);
    assert_eq!(live.source_formats, recorded.source_formats);
    assert_eq!(live.target_formats, recorded.target_formats);
    assert_eq!(
        (
            live.record_count,
            live.converted_count,
            live.warned_count,
            live.failed_count
        ),
        (
            recorded.record_count,
            recorded.converted_count,
            recorded.warned_count,
            recorded.failed_count
        )
    );
    assert_eq!(live.columns, recorded.columns);
    assert_eq!(live.caveats, recorded.caveats);
}

#[test]
fn the_conversion_recordings_still_match_the_live_bridge() {
    let Some(()) = require_verb("peptidoform convert") else {
        return;
    };
    let probe = recorded_conversion("peptidoform_convert_unimod.json");
    assert_same_conversion(&convert(&PROBE).expect("the probe converts"), &probe);

    let psmtsv = recorded_conversion("peptidoform_convert_psmtsv.json");
    let inputs: Vec<String> = psmtsv
        .sequences()
        .unwrap()
        .into_iter()
        .map(|row| row.input)
        .collect();
    assert_same_conversion(&convert(&inputs).unwrap(), &psmtsv);

    let proforma = ConvertOptions {
        target: "ProForma".into(),
        ..Default::default()
    };
    assert_same_conversion(
        &convert_with(&inputs, &proforma).unwrap(),
        &recorded_conversion("peptidoform_convert_psmtsv_proforma.json"),
    );
}

#[test]
fn an_incompatible_modification_is_dropped_with_a_warning_when_asked() {
    let Some(()) = require_verb("peptidoform convert") else {
        return;
    };
    for mode in [
        ConversionMode::RemoveIncompatibleElements,
        ConversionMode::UsePrimarySequence,
    ] {
        let options = ConvertOptions {
            mode,
            ..Default::default()
        };
        let result = convert_with(&PROBE[3..], &options).unwrap();
        assert_eq!(result.mode, mode.as_str());
        let row = &result.sequences().unwrap()[0];
        assert_eq!(row.status, "converted_with_warnings", "{mode}");
        assert_eq!(row.output.as_deref(), Some("PEPKR"), "{mode}");
        assert_eq!(
            row.incompatible_items,
            ["Made Up:Not a modification on K @3(K)"],
            "{mode}"
        );
    }
}

#[test]
fn throw_exception_fails_the_call_naming_the_first_input_it_could_not_convert() {
    let Some(()) = require_verb("peptidoform convert") else {
        return;
    };
    let options = ConvertOptions {
        mode: ConversionMode::ThrowException,
        ..Default::default()
    };
    match convert_with(&PROBE, &options) {
        Err(mzlib::MzLibError::Usage(message)) => {
            assert!(
                message.contains("Made Up:Not a modification on K"),
                "{message}"
            );
        }
        other => panic!("expected a usage error, got {other:?}"),
    }
    // Nothing to throw on: the same rows as ReturnNull.
    let clean = convert_with(&PROBE[..3], &options).unwrap();
    assert_eq!(clean.converted_count, 3);
}

#[test]
fn an_unregistered_format_is_a_usage_error_listing_the_registered_ones() {
    let Some(()) = require_verb("peptidoform convert") else {
        return;
    };
    let options = ConvertOptions {
        target: "Mascot".into(),
        ..Default::default()
    };
    match convert_with(&PROBE, &options) {
        Err(mzlib::MzLibError::Usage(message)) => {
            assert!(
                message.contains("Unimod") && message.contains("ProForma"),
                "{message}"
            );
        }
        other => panic!("expected a usage error, got {other:?}"),
    }
    // Format names match case-insensitively, and the envelope spells them as mzLib registered them.
    let options = ConvertOptions {
        source: "MZLIB".into(),
        target: "unimod".into(),
        ..Default::default()
    };
    let result = convert_with(&PROBE[..1], &options).unwrap();
    assert_eq!(
        (result.source_format.as_str(), result.target_format.as_str()),
        ("mzLib", "Unimod")
    );
}

#[test]
fn the_rows_are_the_same_in_input_order_at_any_thread_count() {
    let Some(()) = require_verb("peptidoform convert") else {
        return;
    };
    let mut many: Vec<&str> = Vec::new();
    for _ in 0..50 {
        many.extend(PROBE);
    }
    let one = convert(&many).unwrap();
    let every_core = convert_with(
        &many,
        &ConvertOptions {
            threads: -1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(one.columns, every_core.columns);
    assert_eq!(one.record_count, 200);
    assert_eq!(one.failed_count, 50);
}
