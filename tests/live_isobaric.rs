//! Live checks of the isobaric module against the real bridge.
//!
//! The offline suite in `src/isobaric.rs` pins the projection against pyMzLib's recordings; these
//! pin that the recordings still match what the bridge emits. They need only the bridge, and
//! **skip** without one. Run with `cargo test --features live`.

#![cfg(feature = "live")]

mod support;

use mzlib::isobaric;

fn recorded(name: &str) -> isobaric::IsobaricKits {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name);
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

#[test]
fn the_kit_recordings_still_match_the_live_bridge() {
    let Some(()) = support::require_verb("isobaric kits") else {
        return;
    };
    let all = isobaric::kits(None).expect("every kit");
    let recorded_all = recorded("isobaric_kits.json");
    assert_eq!(all.kits, recorded_all.kits);
    assert_eq!(all.columns, recorded_all.columns);
    assert_eq!(all.caveats, recorded_all.caveats);

    let tmtpro = isobaric::kits(Some("TMT18")).expect("one kit");
    assert_eq!(tmtpro.columns, recorded("isobaric_kits_TMT18.json").columns);
}

#[test]
fn a_metamorpheus_modification_name_resolves_and_a_partial_name_is_refused() {
    let Some(()) = support::require_verb("isobaric kits") else {
        return;
    };
    let itraq = isobaric::kits(Some("iTRAQ-4plex on K")).expect("MetaMorpheus's name resolves");
    assert_eq!(itraq.kits[0].kit, "iTRAQ4");

    let error = isobaric::kits(Some("TMT10plex")).unwrap_err();
    assert!(
        matches!(error, mzlib::MzLibError::Usage(ref m) if m.contains("TMT10plex")),
        "{error}"
    );
}

#[test]
fn tmt16_is_the_first_sixteen_channels_of_tmt18() {
    let Some(()) = support::require_verb("isobaric kits") else {
        return;
    };
    let all = isobaric::kits(None).unwrap();
    let tmt16 = all.kit("TMT16").unwrap().reporter_ion_mzs();
    let tmt18 = all.kit("TMT18").unwrap().reporter_ion_mzs();
    assert_eq!(tmt16, tmt18[..16]);
}
