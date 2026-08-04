//! Predict peptide properties: retention time, fragment intensities, CCS, detectability.
//!
//! mzLib ships clients for **37 published models** on the
//! [Koina](https://koina.wilhelmlab.org/) inference server, across five families. This calls them.
//!
//! ```no_run
//! # fn main() -> Result<(), mzlib::MzLibError> {
//! let r = mzlib::prediction::retention_time("Prosit_2019_irt", &["PEPTIDEK", "ELVISLIVESK"])?;
//! println!("{:?} in {}", r.columns.floats("retention_time")?, r.retention_time_unit);
//! //  [Some(5.5165), Some(129.7723)] in indexed_retention_time
//! # Ok(())
//! # }
//! ```
//!
//! Start with [`models`]. It is enumerated from mzLib rather than transcribed, and each entry
//! carries the constraints that decide whether a peptide can be sent at all — length bounds
//! (frequently 30, which excludes a lot of real tryptic peptides), which UNIMOD modifications the
//! model was trained on, and whether a collision energy or instrument type is required.
//!
//! ## Five things worth knowing before the first call
//!
//! **Koina is someone else's GPU** — public, shared, community-run, free and unauthenticated.
//! [`Politeness`] exists for a genuinely large job; its defaults are mzLib's and are not raised
//! here, because a binding that maximised throughput out of the box would be spending capacity
//! nobody here pays for.
//!
//! **A prediction is an opinion, not a measurement.** Nothing here has been matched against a
//! spectrum and no output is FDR-anything.
//!
//! **Retention time is not always in minutes.** Most of these models return *indexed* retention
//! time — iRT, a dimensionless scale anchored to standard peptides — and only `Chronologer_RT`
//! returns absolute minutes. [`Predictions::retention_time_unit`] says which. An iRT of 130 looks
//! exactly like a plausible 130-minute gradient, which is why the unit is a value rather than
//! prose.
//!
//! **Fragment arrays are ragged.** Koina returns a fixed-width grid with `-1` marking ions that
//! cannot exist for a peptide, and mzLib drops those, so each row's arrays are as long as *that*
//! peptide's possible ions — 28 for `PEPTIDEK` against a model whose nominal count is 174. Read
//! them with [`Table::float_arrays`], never as a rectangle.
//!
//! **A peptide that cannot be predicted still gets a row**, with `None` and a warning saying why,
//! so predictions always line up with the peptides you sent.
//!
//! ## Not exposed
//!
//! The **local TorchSharp Chronologer**: x64-only, extracting hundreds of megabytes of weights to a
//! shared temp path, and racing any concurrent process doing the same. `Chronologer_RT` reaches the
//! same model over Koina — but the two report *different units*, absolute retention time over the
//! network and % acetonitrile locally.

use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};
use crate::table::Table;

/// What a model requires of one optional input parameter.
///
/// A **tri-state**, and mzLib expresses it as a nullable set whose emptiness means the opposite of
/// what it looks like: `null` means *not applicable* and an *empty* set means *required, any
/// value*. Reading the raw collection is how you conclude that `Prosit_2020_intensity_CID` accepts
/// any collision energy — it accepts none, being fixed at NCE 35 — and that
/// `Prosit_2020_intensity_HCD` accepts none, when it requires one.
#[derive(Debug, Clone, Default, PartialEq, Eq, Deserialize)]
#[serde(tag = "requirement", rename_all = "snake_case")]
pub enum Constraint {
    /// The model has no such input. Do not send it.
    ///
    /// The default, because a payload that omits a constraint is describing a family that has no
    /// such parameter — not one that will accept anything.
    #[default]
    NotApplicable,
    /// You must send one; any value is accepted.
    AnyValueRequired,
    /// You must send one of these.
    OneOf {
        /// The accepted values, as JSON — numbers for charge and collision energy, strings for
        /// instrument and fragmentation type.
        #[serde(default)]
        values: Vec<serde_json::Value>,
    },
}

impl Constraint {
    /// Whether this parameter should be sent at all.
    #[must_use]
    pub fn applicable(&self) -> bool {
        !matches!(self, Self::NotApplicable)
    }
}

/// One Koina model mzLib can call, with the constraints that decide what you may send it.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct Model {
    /// The model's published Koina name — what you pass as `model`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub model: String,
    /// `retention_time`, `fragment_intensity`, `collisional_cross_section`, `detectability` or
    /// `crosslink_intensity`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub family: String,
    /// The bridge verb this family is called through, for cross-referencing.
    #[serde(default)]
    pub verb: Option<String>,
    /// The mzLib class name, for cross-referencing the mzLib source.
    #[serde(default)]
    pub r#type: Option<String>,
    /// Shortest base sequence the model accepts.
    #[serde(default)]
    pub min_peptide_length: Option<i64>,
    /// Longest. **Frequently 30**, which excludes a lot of real tryptic peptides.
    #[serde(default)]
    pub max_peptide_length: Option<i64>,
    /// Sequences per request as the server accepts them; mzLib batches for you.
    #[serde(default)]
    pub max_batch_size: Option<i64>,
    /// UNIMOD accessions the model was trained on. **Empty means the model accepts no
    /// modifications at all**, which is a real answer rather than a missing one.
    #[serde(default)]
    pub allowed_unimod_ids: Vec<i64>,
    /// Whether a precursor charge is required, and which.
    #[serde(default)]
    pub precursor_charge: Constraint,
    /// Whether a collision energy is required, and which.
    #[serde(default)]
    pub collision_energy: Constraint,
    /// Whether an instrument type is required, and which.
    #[serde(default)]
    pub instrument_type: Constraint,
    /// Whether a fragmentation type is required, and which.
    #[serde(default)]
    pub fragmentation_type: Constraint,
    /// `indexed_retention_time` or `minutes` for the retention-time family; `None` for the rest.
    #[serde(default)]
    pub retention_time_unit: Option<String>,
    /// The model's nominal ion count, or `None` when it is dynamic. **Not** the length of any
    /// row's fragment arrays — see the module docs.
    #[serde(default)]
    pub number_of_predicted_fragment_ions: Option<i64>,
    /// Set when the model could not be constructed, which means a broken mzLib build.
    #[serde(default)]
    pub error: Option<String>,
}

impl Model {
    /// Whether the model was trained on any modifications at all.
    #[must_use]
    pub fn accepts_modifications(&self) -> bool {
        !self.allowed_unimod_ids.is_empty()
    }
}

/// A table of predictions, one row per peptide sent.
#[derive(Debug, Clone)]
pub struct Predictions {
    /// The model that produced them.
    pub model: String,
    /// Rows returned — always equal to the number of peptides sent, so predictions line up with
    /// inputs even where some could not be predicted.
    pub row_count: u64,
    /// Rows whose prediction is `null` with a warning explaining why. **Not an error**: too long,
    /// an untrained modification, or a missing required parameter are normal outcomes.
    pub failed_row_count: u64,
    /// The predictions. Empty when the table was written to disk.
    pub columns: Table,
    /// `indexed_retention_time` or `minutes` for [`retention_time`]; empty for other verbs.
    pub retention_time_unit: String,
    /// `square_angstroms` for [`ccs`]; empty for other verbs.
    pub collisional_cross_section_unit: String,
    /// `relative` for the fragment verbs; empty for other verbs.
    pub intensity_scale: String,
    /// What these numbers cannot be trusted to mean.
    pub caveats: Vec<String>,
    /// Where the table was written, or `None` if it came back inline.
    pub output: Option<WrittenTable>,
}

impl Predictions {
    /// `(row index, message)` for every row that could not be predicted.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if the `warning` column is not strings.
    pub fn warnings(&self) -> Result<Vec<(usize, String)>> {
        if !self.columns.has("warning") {
            return Ok(Vec::new());
        }

        Ok(self
            .columns
            .strings("warning")?
            .into_iter()
            .enumerate()
            .filter_map(|(row, warning)| warning.map(|message| (row, message)))
            .collect())
    }
}

/// Where a prediction wrote its table.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct WrittenTable {
    /// The absolute path written.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Always `"tsv"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub format: String,
    /// Rows written, excluding the header.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
}

/// How hard to lean on a shared community server, and where to put the answer.
///
/// The defaults are mzLib's own and are deliberately not raised. Koina is public, free and
/// unauthenticated; these knobs exist for a genuinely large job.
#[derive(Debug, Clone, Default)]
pub struct Politeness {
    /// In-flight requests. `None` uses mzLib's default.
    pub max_batches: Option<u32>,
    /// Delay between request chunks, in milliseconds. `None` uses mzLib's default.
    pub throttle_ms: Option<u32>,
    /// Write the table here as a tab-separated file and return only a summary.
    pub out: Option<String>,
    /// Seconds to allow. `None` waits indefinitely, which a large batch legitimately needs.
    ///
    /// Note this bounds the **wait**, not the work: mzLib's public prediction API threads no
    /// cancellation token, so a timed-out batch is still running on the server.
    pub timeout: Option<Duration>,
}

/// One peptide's inputs.
///
/// Built from a bare sequence with [`From`], so the common case stays a one-liner, and extended
/// with the builder methods where a model needs more.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Peptide {
    /// The sequence, in mzLib `FullSequence` notation — **except** for the crosslink family, which
    /// requires raw UNIMOD brackets. See [`crosslink_fragments`].
    pub sequence: String,
    /// The partner sequence, for the crosslink family only.
    pub beta_sequence: Option<String>,
    /// Required by the CCS and fragment families.
    pub precursor_charge: Option<i32>,
    /// Required by many fragment models — check [`Model::collision_energy`].
    pub collision_energy: Option<i32>,
    /// Required by a few fragment models, e.g. `"QE"`, `"LUMOS"`.
    pub instrument_type: Option<String>,
    /// Required by a few fragment models, e.g. `"HCD"`, `"CID"`.
    pub fragmentation_type: Option<String>,
}

impl<T: Into<String>> From<T> for Peptide {
    fn from(sequence: T) -> Self {
        Self {
            sequence: sequence.into(),
            ..Self::default()
        }
    }
}

impl Peptide {
    /// Sets the precursor charge.
    #[must_use]
    pub fn charge(mut self, charge: i32) -> Self {
        self.precursor_charge = Some(charge);
        self
    }

    /// Sets the collision energy.
    #[must_use]
    pub fn collision_energy(mut self, energy: i32) -> Self {
        self.collision_energy = Some(energy);
        self
    }

    /// Sets the instrument type.
    #[must_use]
    pub fn instrument(mut self, instrument: impl Into<String>) -> Self {
        self.instrument_type = Some(instrument.into());
        self
    }

    /// Sets the fragmentation type.
    #[must_use]
    pub fn fragmentation(mut self, fragmentation: impl Into<String>) -> Self {
        self.fragmentation_type = Some(fragmentation.into());
        self
    }

    /// Sets the crosslink partner sequence.
    #[must_use]
    pub fn beta(mut self, beta: impl Into<String>) -> Self {
        self.beta_sequence = Some(beta.into());
        self
    }

    fn cell(&self, column: &str) -> String {
        let value = match column {
            "sequence" | "alpha_sequence" => Some(self.sequence.clone()),
            "beta_sequence" => self.beta_sequence.clone(),
            "precursor_charge" => self.precursor_charge.map(|v| v.to_string()),
            "collision_energy" => self.collision_energy.map(|v| v.to_string()),
            "instrument_type" => self.instrument_type.clone(),
            "fragmentation_type" => self.fragmentation_type.clone(),
            _ => None,
        };
        value.unwrap_or_default()
    }
}

/// Every Koina model mzLib can call.
///
/// Enumerated from mzLib itself, so it reflects the installed version and cannot go stale.
/// **This does not touch the network** — it describes what mzLib can call, not what the server
/// currently answers.
///
/// # Errors
///
/// [`MzLibError::Bridge`] if mzLib itself failed; [`MzLibError::Protocol`] if the payload cannot
/// be interpreted.
pub fn models() -> Result<Vec<Model>> {
    models_in(None)
}

/// [`models`], restricted to one family, e.g. `"retention_time"`.
///
/// # Errors
///
/// As [`models`], plus [`MzLibError::Usage`] for a blank family.
pub fn models_in(family: Option<&str>) -> Result<Vec<Model>> {
    #[derive(Deserialize)]
    struct Payload {
        #[serde(default)]
        models: Vec<Model>,
    }

    let mut args = vec!["predict".to_owned(), "models".to_owned()];
    if let Some(family) = family {
        if family.trim().is_empty() {
            return Err(MzLibError::Usage(
                "family must be a non-empty name, or None for every family.".to_owned(),
            ));
        }
        args.push("--family".to_owned());
        args.push(family.trim().to_owned());
    }

    let data = bridge::invoke(&args, None, Some(Duration::from_secs(60)))?;
    let payload: Payload = serde_json::from_value(data).map_err(protocol)?;
    Ok(payload.models)
}

/// Predict elution, one row per peptide.
///
/// **Check [`Predictions::retention_time_unit`]** — most of these models return iRT, not minutes.
///
/// # Errors
///
/// [`MzLibError::Usage`] for a blank model or an empty peptide list, or when the model is not in
/// this family — the message names the ones that are.
pub fn retention_time(model: &str, peptides: &[impl Clone + Into<Peptide>]) -> Result<Predictions> {
    retention_time_with(model, peptides, &Politeness::default())
}

/// [`retention_time`], with explicit throttling and destination.
///
/// # Errors
///
/// As [`retention_time`].
pub fn retention_time_with(
    model: &str,
    peptides: &[impl Clone + Into<Peptide>],
    politeness: &Politeness,
) -> Result<Predictions> {
    predict("retention-time", model, peptides, &["sequence"], politeness)
}

/// Predict MS2 fragment m/z and relative intensity.
///
/// The `fragment_annotations`, `fragment_mz` and `fragment_intensity` columns each hold **one array
/// per row**, of differing lengths — read them with [`Table::float_arrays`] and
/// [`Table::string_arrays`].
///
/// # Errors
///
/// As [`retention_time`].
pub fn fragments(model: &str, peptides: &[impl Clone + Into<Peptide>]) -> Result<Predictions> {
    fragments_with(model, peptides, &Politeness::default())
}

/// [`fragments`], with explicit throttling and destination.
///
/// # Errors
///
/// As [`fragments`].
pub fn fragments_with(
    model: &str,
    peptides: &[impl Clone + Into<Peptide>],
    politeness: &Politeness,
) -> Result<Predictions> {
    predict(
        "fragments",
        model,
        peptides,
        &[
            "sequence",
            "precursor_charge",
            "collision_energy",
            "instrument_type",
            "fragmentation_type",
        ],
        politeness,
    )
}

/// Predict collisional cross-section, in **square angstroms — never 1/K0**.
///
/// Converting to the reduced mobility a timsTOF reports needs drift-gas temperature and pressure,
/// which mzLib does not carry, so no conversion is offered rather than a guessed one.
///
/// # Errors
///
/// As [`retention_time`].
pub fn ccs(model: &str, peptides: &[impl Clone + Into<Peptide>]) -> Result<Predictions> {
    ccs_with(model, peptides, &Politeness::default())
}

/// [`ccs`], with explicit throttling and destination.
///
/// # Errors
///
/// As [`ccs`].
pub fn ccs_with(
    model: &str,
    peptides: &[impl Clone + Into<Peptide>],
    politeness: &Politeness,
) -> Result<Predictions> {
    predict(
        "ccs",
        model,
        peptides,
        &["sequence", "precursor_charge"],
        politeness,
    )
}

/// Predict flyability, as four class probabilities that sum to 1.
///
/// They are a distribution over classes, **not** an expected intensity and not a probability of
/// detection.
///
/// # Errors
///
/// As [`retention_time`].
pub fn detectability(model: &str, peptides: &[impl Clone + Into<Peptide>]) -> Result<Predictions> {
    detectability_with(model, peptides, &Politeness::default())
}

/// [`detectability`], with explicit throttling and destination.
///
/// # Errors
///
/// As [`detectability`].
pub fn detectability_with(
    model: &str,
    peptides: &[impl Clone + Into<Peptide>],
    politeness: &Politeness,
) -> Result<Predictions> {
    predict("detectability", model, peptides, &["sequence"], politeness)
}

/// Predict MS2 intensities for a crosslinked peptide pair.
///
/// # Errors
///
/// As [`retention_time`].
///
/// # Sequence notation
///
/// **This family takes a different sequence language from every other function here.** The others
/// accept mzLib's `FullSequence` notation and convert it; the crosslink models reject it and
/// require raw UNIMOD brackets — `K[UNIMOD:1896]`. That is mzLib's constraint, not a choice made
/// here, and it is repeated in [`Predictions::caveats`].
pub fn crosslink_fragments(
    model: &str,
    pairs: &[impl Clone + Into<Peptide>],
) -> Result<Predictions> {
    crosslink_fragments_with(model, pairs, &Politeness::default())
}

/// [`crosslink_fragments`], with explicit throttling and destination.
///
/// # Errors
///
/// As [`crosslink_fragments`].
pub fn crosslink_fragments_with(
    model: &str,
    pairs: &[impl Clone + Into<Peptide>],
    politeness: &Politeness,
) -> Result<Predictions> {
    predict(
        "crosslink-fragments",
        model,
        pairs,
        &[
            "alpha_sequence",
            "beta_sequence",
            "precursor_charge",
            "collision_energy",
        ],
        politeness,
    )
}

/// The shared body of every predict verb: validate, send, parse.
fn predict(
    verb: &str,
    model: &str,
    peptides: &[impl Clone + Into<Peptide>],
    columns: &[&str],
    politeness: &Politeness,
) -> Result<Predictions> {
    if model.trim().is_empty() {
        return Err(MzLibError::Usage(
            "A model name is required. prediction::models() lists them with their constraints."
                .to_owned(),
        ));
    }
    if peptides.is_empty() {
        return Err(MzLibError::Usage(
            "At least one peptide is required.".to_owned(),
        ));
    }

    let mut args = vec![
        "predict".to_owned(),
        verb.to_owned(),
        "--model".to_owned(),
        model.trim().to_owned(),
    ];

    if let Some(out) = &politeness.out {
        if out.trim().is_empty() {
            return Err(MzLibError::Usage(
                "out must be a non-empty path, or None to return the table.".to_owned(),
            ));
        }
        args.push("--out".to_owned());
        args.push(out.trim().to_owned());
    }
    if let Some(max_batches) = politeness.max_batches {
        if max_batches == 0 {
            return Err(MzLibError::Usage(
                "max_batches must be 1 or greater; None uses mzLib's default.".to_owned(),
            ));
        }
        args.push("--max-batches".to_owned());
        args.push(max_batches.to_string());
    }
    if let Some(throttle) = politeness.throttle_ms {
        args.push("--throttle-ms".to_owned());
        args.push(throttle.to_string());
    }

    // On stdin rather than argv: a real prediction run is thousands of peptides and argv has a
    // hard ceiling of roughly 32 KB.
    let mut stdin = columns.join("\t");
    for peptide in peptides {
        let peptide: Peptide = peptide.clone().into();
        stdin.push('\n');
        stdin.push_str(
            &columns
                .iter()
                .map(|column| peptide.cell(column))
                .collect::<Vec<_>>()
                .join("\t"),
        );
    }
    stdin.push('\n');

    let data = bridge::invoke(&args, Some(&stdin), politeness.timeout)?;

    #[derive(Deserialize)]
    struct Wire {
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        model: String,
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        row_count: u64,
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        failed_row_count: u64,
        #[serde(default)]
        column_names: Vec<String>,
        #[serde(default)]
        columns: Option<std::collections::BTreeMap<String, Vec<serde_json::Value>>>,
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        retention_time_unit: String,
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        collisional_cross_section_unit: String,
        #[serde(default, deserialize_with = "bridge::null_to_default")]
        intensity_scale: String,
        #[serde(default)]
        caveats: Vec<String>,
        #[serde(default)]
        output: Option<WrittenTable>,
    }

    let wire: Wire = serde_json::from_value(data).map_err(protocol)?;
    Ok(Predictions {
        model: wire.model,
        row_count: wire.row_count,
        failed_row_count: wire.failed_row_count,
        columns: Table::from_wire(wire.column_names, wire.columns),
        retention_time_unit: wire.retention_time_unit,
        collisional_cross_section_unit: wire.collisional_cross_section_unit,
        intensity_scale: wire.intensity_scale,
        caveats: wire.caveats,
        output: wire.output,
    })
}

fn protocol(error: serde_json::Error) -> MzLibError {
    MzLibError::Protocol(format!(
        "prediction payload could not be interpreted: {error}"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_bare_sequence_becomes_a_peptide() {
        let peptide: Peptide = "PEPTIDEK".into();
        assert_eq!(peptide.sequence, "PEPTIDEK");
        assert_eq!(peptide.precursor_charge, None);
    }

    #[test]
    fn the_builder_fills_the_optional_columns() {
        let peptide = Peptide::from("PEPTIDEK")
            .charge(2)
            .collision_energy(28)
            .instrument("QE");

        assert_eq!(peptide.cell("precursor_charge"), "2");
        assert_eq!(peptide.cell("collision_energy"), "28");
        assert_eq!(peptide.cell("instrument_type"), "QE");
        // An unset optional is an EMPTY cell, not the text "None": the latter would reach the
        // bridge as a value and fail to parse as an integer, reporting a confusing error about a
        // column the caller never set.
        assert_eq!(peptide.cell("fragmentation_type"), "");
    }

    #[test]
    fn the_alpha_column_reads_the_same_field_as_sequence() {
        // The crosslink family names its first column alpha_sequence, and a Peptide built from a
        // bare string must still fill it.
        let peptide: Peptide = "PEPTIDEK".into();
        assert_eq!(peptide.cell("alpha_sequence"), "PEPTIDEK");
    }

    #[test]
    fn a_blank_model_is_refused_before_anything_is_spawned() {
        let error = predict(
            "retention-time",
            "   ",
            &["PEPTIDEK"],
            &["sequence"],
            &Politeness::default(),
        )
        .expect_err("a blank model must not reach the bridge");
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    #[test]
    fn an_empty_peptide_list_is_refused() {
        let empty: [&str; 0] = [];
        let error = predict(
            "retention-time",
            "Prosit_2019_irt",
            &empty,
            &["sequence"],
            &Politeness::default(),
        )
        .expect_err("an empty list must not reach the bridge");
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    #[test]
    fn a_zero_max_batches_is_refused() {
        let error = predict(
            "retention-time",
            "Prosit_2019_irt",
            &["PEPTIDEK"],
            &["sequence"],
            &Politeness {
                max_batches: Some(0),
                ..Politeness::default()
            },
        )
        .expect_err("zero in-flight requests would never finish");
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    #[test]
    fn a_constraint_distinguishes_not_applicable_from_required_any() {
        // The trap this enum exists for: mzLib expresses both as a nullable set, where null means
        // "no such input" and EMPTY means "required, any value".
        let not_applicable: Constraint =
            serde_json::from_str(r#"{"requirement":"not_applicable","values":null}"#).unwrap();
        let required: Constraint =
            serde_json::from_str(r#"{"requirement":"any_value_required","values":null}"#).unwrap();
        let one_of: Constraint =
            serde_json::from_str(r#"{"requirement":"one_of","values":[20,21]}"#).unwrap();

        assert!(!not_applicable.applicable());
        assert!(required.applicable());
        assert!(matches!(one_of, Constraint::OneOf { ref values } if values.len() == 2));
    }

    #[test]
    fn a_model_with_no_allowed_modifications_says_so() {
        let model: Model = serde_json::from_str(
            r#"{"model":"pfly_2024_fine_tuned","family":"detectability","allowed_unimod_ids":[]}"#,
        )
        .unwrap();
        assert!(!model.accepts_modifications());
    }
}
