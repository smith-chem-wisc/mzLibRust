//! Live checks of the prediction module against the real bridge, and — for a few — the real Koina
//! server.
//!
//! The catalogue tests need only the bridge: `predict models` reflects over the loaded assembly and
//! touches no network, which makes them the cheapest possible guard on the fact this crate most
//! depends on — that the constraints it describes are the ones mzLib declares.
//!
//! The prediction tests do reach Koina, a public community server. They **skip** rather than fail
//! when it is unavailable, because a red build that might just mean "someone else's GPU is busy" is
//! an ambiguous red build, and an ambiguous red build gets ignored.
//!
//! Run with `cargo test --features live`.

#![cfg(feature = "live")]

mod support;

use support::{external_service, require_bridge};

// ---- the catalogue (bridge only, no network) --------------------------------------------------

#[test]
fn every_model_mzlib_can_call_is_constructible_and_described() {
    let Some(()) = require_bridge() else { return };

    let models = mzlib::prediction::models().expect("the models verb should answer");

    assert_eq!(
        models.len(),
        37,
        "mzLib ships 37 Koina models; a change here means this crate's documented count is stale"
    );
    for model in &models {
        assert!(
            model.error.is_none(),
            "{} could not be constructed: {:?}",
            model.r#type.as_deref().unwrap_or("?"),
            model.error
        );
        assert!(!model.model.is_empty());
        assert!(
            model.verb.is_some(),
            "every model names the verb that calls it"
        );
    }
}

#[test]
fn the_five_families_are_all_represented() {
    let Some(()) = require_bridge() else { return };

    let models = mzlib::prediction::models().expect("the models verb should answer");
    let families: std::collections::BTreeSet<_> =
        models.iter().map(|m| m.family.as_str()).collect();

    assert_eq!(
        families,
        [
            "collisional_cross_section",
            "crosslink_intensity",
            "detectability",
            "fragment_intensity",
            "retention_time",
        ]
        .into_iter()
        .collect()
    );
}

#[test]
fn the_retention_time_unit_is_per_model_not_per_family() {
    // The distinction that makes mzLib's bare IsIndexed boolean insufficient on the wire.
    let Some(()) = require_bridge() else { return };

    let models = mzlib::prediction::models().expect("the models verb should answer");
    let named = |name: &str| {
        models
            .iter()
            .find(|m| m.model == name)
            .unwrap_or_else(|| panic!("{name} should be in the catalogue"))
            .clone()
    };

    assert_eq!(
        named("Prosit_2019_irt").retention_time_unit.as_deref(),
        Some("indexed_retention_time")
    );
    assert_eq!(
        named("Chronologer_RT").retention_time_unit.as_deref(),
        Some("minutes")
    );
    // ...and it is meaningless outside that family.
    assert_eq!(named("IM2Deep").retention_time_unit, None);
}

#[test]
fn a_constraint_distinguishes_not_applicable_from_required_any() {
    // Reading mzLib's raw nullable set makes CID look permissive about collision energy when it
    // accepts none, and HCD look like it accepts none when it requires one.
    let Some(()) = require_bridge() else { return };

    let models = mzlib::prediction::models().expect("the models verb should answer");
    let energy = |name: &str| {
        models
            .iter()
            .find(|m| m.model == name)
            .unwrap_or_else(|| panic!("{name} should be in the catalogue"))
            .collision_energy
            .clone()
    };

    assert!(matches!(
        energy("Prosit_2020_intensity_HCD"),
        mzlib::prediction::Constraint::AnyValueRequired
    ));
    assert!(matches!(
        energy("Prosit_2020_intensity_CID"),
        mzlib::prediction::Constraint::NotApplicable
    ));
    assert!(matches!(
        energy("Altimeter_2024_intensities"),
        mzlib::prediction::Constraint::OneOf { .. }
    ));
}

#[test]
fn an_unknown_model_names_the_ones_this_verb_has() {
    let Some(()) = require_bridge() else { return };

    let error = mzlib::prediction::retention_time("not_a_model", &["PEPTIDEK"])
        .expect_err("an unknown model must not silently pick one");

    let message = error.to_string();
    assert!(message.contains("Available:"), "{message}");
}

#[test]
fn a_model_from_the_wrong_family_is_refused_by_name() {
    // Prosit_2019_irt is a real model and not a fragment-intensity one. Accepting it would fail
    // much later, inside a request, with an error about the payload rather than the name.
    let Some(()) = require_bridge() else { return };

    let error = mzlib::prediction::fragments("Prosit_2019_irt", &["PEPTIDEK"])
        .expect_err("a model from another family must be refused");

    assert!(error.to_string().contains("No model named"), "{error}");
}

// ---- predictions (Koina) ----------------------------------------------------------------------

#[test]
fn koina_still_answers_and_an_irt_model_is_labelled_as_one() {
    let Some(()) = require_bridge() else { return };

    let Some(result) = external_service(
        "Koina",
        mzlib::prediction::retention_time("Prosit_2019_irt", &["PEPTIDEK", "ELVISLIVESK"]),
    ) else {
        return;
    };

    assert_eq!(result.row_count, 2);
    assert_eq!(result.failed_row_count, 0);
    assert_eq!(result.retention_time_unit, "indexed_retention_time");
    assert!(result
        .columns
        .floats("retention_time")
        .expect("a numeric column")
        .iter()
        .all(Option::is_some));
}

#[test]
fn fragment_arrays_are_ragged_and_index_aligned_within_a_row() {
    let Some(()) = require_bridge() else { return };

    let peptides = [
        mzlib::prediction::Peptide::from("PEPTIDEK")
            .charge(2)
            .collision_energy(28),
        mzlib::prediction::Peptide::from("ELVISLIVESK")
            .charge(2)
            .collision_energy(28),
    ];

    let Some(result) = external_service(
        "Koina",
        mzlib::prediction::fragments("Prosit_2020_intensity_HCD", &peptides),
    ) else {
        return;
    };

    let mz = result.columns.float_arrays("fragment_mz").expect("arrays");
    let intensity = result
        .columns
        .float_arrays("fragment_intensity")
        .expect("arrays");
    let annotations = result
        .columns
        .string_arrays("fragment_annotations")
        .expect("arrays");

    let lengths: Vec<usize> = mz
        .iter()
        .map(|row| row.as_ref().map_or(0, Vec::len))
        .collect();
    // Koina returns a fixed-width grid with -1 for impossible ions and mzLib drops those, so each
    // row is as long as ITS peptide's possible ions. Indexing these as a rectangle is wrong.
    assert!(
        lengths.windows(2).any(|pair| pair[0] != pair[1]),
        "two peptides of different lengths must give arrays of different lengths: {lengths:?}"
    );
    // The model's published count is 174; a short tryptic peptide gets a fraction of it.
    assert!(lengths.iter().all(|&length| length < 174));

    for row in 0..mz.len() {
        let expected = mz[row].as_ref().map_or(0, Vec::len);
        assert_eq!(intensity[row].as_ref().map_or(0, Vec::len), expected);
        assert_eq!(annotations[row].as_ref().map_or(0, Vec::len), expected);
    }
    assert_eq!(result.intensity_scale, "relative");
}

#[test]
fn a_peptide_that_cannot_be_predicted_still_gets_a_row() {
    // Prosit_2020_intensity_HCD requires a collision energy. Omitting it must not lose the row, or
    // predictions would no longer line up with the peptides that were sent.
    let Some(()) = require_bridge() else { return };

    let peptides = [mzlib::prediction::Peptide::from("PEPTIDEK").charge(2)];
    let Some(result) = external_service(
        "Koina",
        mzlib::prediction::fragments("Prosit_2020_intensity_HCD", &peptides),
    ) else {
        return;
    };

    assert_eq!(result.row_count, 1);
    assert_eq!(result.failed_row_count, 1);
    let warnings = result.warnings().expect("a warning column");
    assert!(
        warnings
            .iter()
            .any(|(_, message)| message.contains("CollisionEnergy")),
        "{warnings:?}"
    );
}

#[test]
fn ccs_is_reported_in_square_angstroms() {
    let Some(()) = require_bridge() else { return };

    let peptides = [mzlib::prediction::Peptide::from("PEPTIDEK").charge(2)];
    let Some(result) = external_service("Koina", mzlib::prediction::ccs("IM2Deep", &peptides))
    else {
        return;
    };

    assert_eq!(result.collisional_cross_section_unit, "square_angstroms");
    assert!(result
        .columns
        .floats("collisional_cross_section")
        .expect("a numeric column")[0]
        .is_some_and(|value| value > 0.0));
}

#[test]
fn detectability_returns_four_classes_that_sum_to_one() {
    let Some(()) = require_bridge() else { return };

    let Some(result) = external_service(
        "Koina",
        mzlib::prediction::detectability("pfly_2024_fine_tuned", &["PEPTIDEK"]),
    ) else {
        return;
    };

    let total: f64 = [
        "not_detectable",
        "low_detectability",
        "intermediate_detectability",
        "high_detectability",
    ]
    .into_iter()
    .map(|name| result.columns.floats(name).expect("a numeric column")[0].unwrap_or(0.0))
    .sum();

    assert!(
        (total - 1.0).abs() < 1e-5,
        "the four classes sum to {total}"
    );
}
