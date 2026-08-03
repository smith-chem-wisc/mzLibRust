//! Peptidoform-level questions: digest an annotated protein and fragment its peptides.
//!
//! The question this answers is the one a mass spectrometrist actually asks — *what fragments
//! would I see for this protein's peptides?* — in one call:
//!
//! ```no_run
//! # fn main() -> Result<(), mzlib::MzLibError> {
//! let digest = mzlib::peptidoform::fragments("P02768")?;
//! println!("{} peptides", digest.peptides.len());
//! println!("{}", digest.modification_census.explain());
//! # Ok(())
//! # }
//! ```
//!
//! The defaults are opinions, not placeholders. Tryptic with the proline rule, two missed
//! cleavages, ETD, both termini, UniProt's annotated modifications applied. They are the choices
//! this lab makes when it does not have a reason to choose otherwise, so the common question needs
//! no parameters — and every one of them is reachable, because the point is to open the doors, not
//! to hide them.

use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};

/// The proton mass, in daltons.
///
/// Stated here rather than buried, because which of the two nearby constants a library used is
/// invisible in its output and changes an answer by about 1 ppm. This is the proton
/// (1.007276), **not** the hydrogen atom (1.007825).
pub const PROTON_MASS: f64 = 1.00727646677;

/// One backbone fragment ion.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Fragment {
    /// The ion series, e.g. `"c"` or `"zDot"` for ETD, `"b"`/`"y"` for CID.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub product_type: String,
    /// Position in the series — `c3` is the third from the N-terminus.
    ///
    /// **The `c` and `y` series run `1..=length-1`; the `zDot` series runs `1..=length`.** The extra
    /// `z•` numbered `length` is the whole peptide minus NH₂ (`monoisotopic_mass − 16.01872`), from
    /// the N–Cα cleavage at the *first* residue — real ETD chemistry, deliberate in mzLib, and not a
    /// backbone fragment between two residues. Exclude it if you are counting cleavage sites. It is
    /// absent when the peptide begins with proline.
    ///
    /// **z• ions are suppressed N-terminal to proline; the complementary c ions are not.** ETD
    /// cleaves the N–Cα bond, and at a proline the ring keeps the halves tethered, so neither
    /// fragment should appear — but mzLib only drops the z•. On albumin that leaves 138 c ions
    /// (~4% of the c series) that cannot occur in a real spectrum. See
    /// [smith-chem-wisc/mzLib#1110](https://github.com/smith-chem-wisc/mzLib/issues/1110).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub fragment_number: i32,
    /// Monoisotopic neutral mass in daltons.
    ///
    /// **Not** an m/z: no proton has been added and no charge assumed. Fragments deliberately
    /// expose no `mz()`. Converting one correctly requires the fixed charge *within this
    /// fragment's span* — a c or z ion carries only the permanently charged modifications on the
    /// residues it contains, not the whole peptide's [`Peptide::fixed_charges`]. For an unmodified
    /// or neutrally-modified peptide, `(neutral_mass + z * PROTON_MASS) / z` is correct; for a
    /// fragment bearing a quaternary-ammonium modification it is not, and per-fragment charge
    /// accounting is not yet provided.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub neutral_mass: f64,
    /// Neutral loss in daltons, `0.0` when there is none.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub neutral_loss: f64,
    /// One-based residue position in the peptide.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub residue_position: i32,
}

/// One modification applied to a peptide.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Modification {
    /// One-based position within **the peptide's own residues**, or `None` for a terminal
    /// modification, which carries [`Modification::terminus`] instead.
    ///
    /// mzLib's internal dictionary reserves slot 1 for the N-terminus, so its keys are one past the
    /// residue they modify. That is corrected in the bridge rather than passed on: exposing the raw
    /// key as a "one-based position" was a lie about what the number is — a peptide `MAR` reported
    /// position 4 for its arginine, position 4 of a 3-mer.
    #[serde(default)]
    pub one_based_residue: Option<i32>,
    /// `"N"` or `"C"` for a terminal modification, `None` for one on a residue.
    #[serde(default)]
    pub terminus: Option<String>,
    /// The modification's identifier with its motif, e.g. `"N6-succinyllysine on K"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub id: String,
    /// The monoisotopic mass delta, if the modification has one.
    #[serde(default)]
    pub mass: Option<f64>,
    /// Charges the modification leaves permanently on the residue. See [`Peptide::fixed_charges`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub formal_charge: i32,
}

/// One digested peptide, with its modifications and fragment ions.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Peptide {
    /// The bare amino-acid sequence.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub base_sequence: String,
    /// The sequence with modifications written inline, as mzLib renders them.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub full_sequence: String,
    /// The neutral monoisotopic mass, modifications included.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub monoisotopic_mass: f64,
    /// The peptide's length in residues.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub length: i32,
    /// Start position within the parent protein.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub one_based_start: i32,
    /// End position within the parent protein.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub one_based_end: i32,
    /// How many cleavage sites the peptide spans.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub missed_cleavages: i32,
    /// Charges the peptide carries **before any protonation**, from modifications that leave a
    /// permanently charged residue.
    ///
    /// Trimethylation of a lysine ε-amine gives a quaternary ammonium, and UniProt records the
    /// delta as 43.054227 — C₃H₇ *minus an electron* — rather than the neutral 43.054775. So
    /// [`Peptide::monoisotopic_mass`] already carries that charge, and [`Peptide::mz`] accounts for
    /// it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub fixed_charges: i32,
    /// Each applied modification.
    #[serde(default)]
    pub modifications: Vec<Modification>,
    /// The fragment ions for the requested dissociation type.
    #[serde(default)]
    pub fragments: Vec<Fragment>,
}

impl Peptide {
    /// Whether this peptide carries at least one modification.
    #[must_use]
    pub fn is_modified(&self) -> bool {
        !self.modifications.is_empty()
    }

    /// The m/z of the intact peptide at a given total charge.
    ///
    /// Two conventions are handled explicitly here, because getting either wrong is invisible in
    /// the answer.
    ///
    /// **The proton mass, not the hydrogen atom.** The difference is 0.55 mDa — 1.1 ppm at m/z 500,
    /// which on an Orbitrap is a match versus a miss. Libraries differ on this and rarely say which
    /// they used. See [`PROTON_MASS`].
    ///
    /// **Fixed charges are not double-counted.** A mass that already carries a charge needs fewer
    /// protons added, not the same number — so only `charge - fixed_charges` protons are added.
    /// Adding a full complement would put a 2+ trimethylated peptide half a Thomson high, on the
    /// most important histone modification there is.
    ///
    /// A peptide with a fixed charge is therefore observable at that charge with no protonation at
    /// all, which is why `charge` may not be below [`Peptide::fixed_charges`].
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `charge` is below 1, or below this peptide's fixed charge.
    pub fn mz(&self, charge: i32) -> Result<f64> {
        if charge < 1 {
            return Err(MzLibError::Usage(format!(
                "charge must be a positive whole number; got {charge}."
            )));
        }
        if charge < self.fixed_charges {
            return Err(MzLibError::Usage(format!(
                "This peptide already carries {} fixed charge(s) from its modifications, so it \
                 cannot be observed at charge {charge}.",
                self.fixed_charges
            )));
        }
        Ok(
            (self.monoisotopic_mass + f64::from(charge - self.fixed_charges) * PROTON_MASS)
                / f64::from(charge),
        )
    }
}

/// One UniProt feature type, and whether mzLib could use it.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct FeatureType {
    /// The UniProt feature type, e.g. `"modified residue"`, `"glycosylation site"`.
    #[serde(rename = "type", default, deserialize_with = "bridge::null_to_default")]
    pub feature_type: String,
    /// How many features of this type the entry annotates.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub count: u32,
    /// Whether mzLib loads this type at all.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub loaded: bool,
}

/// What UniProt annotates, and what could actually be used.
///
/// mzLib loads only `modified residue` and `lipid moiety-binding region` annotations; every other
/// feature type is dropped **on feature type alone**, before any mass lookup. So the census sees the
/// world at feature-*type* granularity — one entry per type in [`ModificationCensus::by_type`], never
/// per modification *name*. On serum albumin the 24 excluded features all sit under the single type
/// `glycosylation site`; at UniProt's finer name level 22 of those 24 are specifically
/// `N-linked (Glc) (glycation) lysine`, which *does* have a defined mass — but the census never
/// surfaces that name, so "22" is a fact you confirm by reading the UniProt entry, not a number this
/// type reports. Read the exclusion as "wrong feature type", **not** "no defined mass".
///
/// The exclusion is still correct: glycation and glycosylation are labile, heterogeneous adducts, so
/// assigning one an exact mass and a clean fragment ladder would describe a species you cannot
/// observe. What this type exists for is that you should not have to *guess* it happened: for serum
/// albumin, 14 modifications are applied out of 38 annotated, and without this the 14 arrives with no
/// indication that a rule was ever applied. See smith-chem-wisc/mzLib#1112.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ModificationCensus {
    /// Distinct residue positions carrying at least one modification.
    ///
    /// A histone lists several alternatives at one residue — K9me1, K9me2, K9me3, K9ac are four
    /// modifications at one site — so this is always the smaller number and is **not** a
    /// modification count.
    pub sites: u32,
    /// Modifications actually placed on the protein.
    pub applied: u32,
    /// Modification-like features UniProt lists.
    pub annotated: u32,
    /// Modification names UniProt annotated that could not be resolved to a mass — usually because
    /// the name is absent from UniProt's own ptmlist.
    ///
    /// These vanish silently otherwise: on histone H3.1, seven N6-lactoyllysine sites were dropped
    /// while the type summary still reported "modified residue … loaded".
    pub unresolved: Vec<String>,
    /// One entry per feature type, with its count and whether it was loaded.
    pub by_type: Vec<FeatureType>,
}

impl ModificationCensus {
    /// Annotated features that could not be used because they have no defined mass.
    #[must_use]
    pub fn excluded(&self) -> u32 {
        self.annotated.saturating_sub(self.applied)
    }

    /// A one-paragraph, human-readable account of what was used and what was not.
    ///
    /// It names the excluded feature *types* and their counts — the only granularity the census
    /// has. It never reports a modification-*name*-level breakdown (e.g. "22 of 24 are glycation"),
    /// because [`ModificationCensus::by_type`] does not carry names; such a figure comes from
    /// reading the UniProt entry, not from this census.
    #[must_use]
    pub fn explain(&self) -> String {
        if self.excluded() == 0 {
            return format!(
                "All {} annotated modifications were applied, across {} residue positions.",
                self.annotated, self.sites
            );
        }

        let mut sentences = vec![format!(
            "{} of {} annotated modifications were applied, across {} residue positions.",
            self.applied, self.annotated, self.sites
        )];

        let excluded_types: Vec<&FeatureType> =
            self.by_type.iter().filter(|entry| !entry.loaded).collect();
        if !excluded_types.is_empty() {
            let named = excluded_types
                .iter()
                .map(|entry| format!("{} × {}", entry.count, entry.feature_type))
                .collect::<Vec<_>>()
                .join(", ");
            sentences.push(format!(
                "Excluded by type: {named} — mzLib loads only 'modified residue' and 'lipid \
                 moiety-binding region' annotations, so these were dropped on feature type alone. \
                 The exclusion is usually right: a glycation or glycosylation annotation describes \
                 a labile, heterogeneous adduct, so assigning it one exact mass and a clean \
                 fragment ladder would invent a species you cannot observe. But the reason is not \
                 reported, and the qualifier is not read — some annotations are marked 'in vitro' \
                 and some exist only in disease variants, which are different grounds for exclusion \
                 needing different judgements from you. Read the annotations on the UniProt entry \
                 before concluding anything about a specific site; this census can only tell you \
                 the count (smith-chem-wisc/mzLib#1112)."
            ));
        }

        if !self.unresolved.is_empty() {
            sentences.push(format!(
                "Could not be resolved to a mass: {} — annotated by UniProt but absent from its \
                 own modification list, so they were dropped.",
                self.unresolved.join(", ")
            ));
        }

        sentences.join(" ")
    }
}

/// The result of digesting a protein and fragmenting its peptides.
#[derive(Debug, Clone, PartialEq)]
pub struct Digest {
    /// The UniProt accession that was fetched.
    pub accession: String,
    /// The entry's short name, e.g. `"ALBU_HUMAN"`.
    pub name: String,
    /// The entry's full protein name, e.g. `"Albumin"`.
    pub full_name: String,
    /// The source organism.
    pub organism: String,
    /// The protein's length in residues.
    pub sequence_length: u32,
    /// The protease used, in mzLib's naming.
    pub protease: String,
    /// The dissociation type used.
    pub dissociation: String,
    /// Which terminus (or both) was fragmented.
    pub terminus: String,
    /// Whether UniProt's annotated modifications were applied.
    pub modifications_applied: bool,
    /// The maximum modifications considered per peptide.
    pub max_modifications: u32,
    /// The maximum modification isoforms allowed per peptide position.
    pub max_isoforms: u32,
    /// How many peptides hit that isoform cap.
    pub peptides_at_cap: u32,
    /// What UniProt annotated versus what was applied.
    pub modification_census: ModificationCensus,
    /// The digested peptides, each with its fragments.
    pub peptides: Vec<Peptide>,
}

impl Digest {
    /// Whether any peptide hit the isoform cap, meaning the result is incomplete.
    ///
    /// A short answer and a truncated answer look identical from the outside. Check this before
    /// treating a peptide list as exhaustive.
    #[must_use]
    pub fn truncated(&self) -> bool {
        self.peptides_at_cap > 0
    }

    /// Only the peptides carrying at least one modification.
    #[must_use]
    pub fn modified_peptides(&self) -> Vec<&Peptide> {
        self.peptides.iter().filter(|p| p.is_modified()).collect()
    }

    /// How many peptides there are, counted as **distinct base sequences**.
    ///
    /// [`Self::peptides`] holds *peptidoforms* — one entry per sequence-and-modification-placement
    /// — so `peptides.len()` is much the larger number whenever modifications are applied. On
    /// albumin at two modifications it is 303 peptidoforms over 195 distinct sequences. Both are
    /// legitimate answers to "how many peptides"; they are not interchangeable, and quoting one for
    /// the other is a very large error rather than a rounding one.
    #[must_use]
    pub fn distinct_base_sequences(&self) -> usize {
        self.peptides
            .iter()
            .map(|peptide| peptide.base_sequence.as_str())
            .collect::<std::collections::HashSet<_>>()
            .len()
    }

    /// Fragment ions per product type, e.g. `{"c": 3253, "y": 3253, "zDot": 3308}`.
    ///
    /// Prefer this to [`Self::fragment_count`] whenever the ion series matter, which for ETD is
    /// always: mzLib emits a spurious `y` series for ETD (about a third of the total, see
    /// [`FragmentOptions::dissociation`]), and `zDot` carries one extra full-length ion per peptide
    /// (see [`Fragment::fragment_number`]). A single total silently folds both in.
    #[must_use]
    pub fn fragments_by_series(&self) -> std::collections::BTreeMap<String, usize> {
        let mut counts = std::collections::BTreeMap::new();
        for fragment in self.peptides.iter().flat_map(|p| p.fragments.iter()) {
            *counts.entry(fragment.product_type.clone()).or_insert(0) += 1;
        }
        counts
    }

    /// Total fragment ions across every peptide.
    ///
    /// A bare total. For ETD this includes the spurious `y` series and the full-length `z•` per
    /// peptide, so it is **not** the number of ions you would look for in a spectrum — see
    /// [`Self::fragments_by_series`].
    #[must_use]
    pub fn fragment_count(&self) -> usize {
        self.peptides.iter().map(|p| p.fragments.len()).sum()
    }
}

/// How a protein is digested and fragmented. Every default is the lab's opinion, and reachable.
#[derive(Debug, Clone)]
pub struct FragmentOptions {
    /// **Read this if you are coming from MaxQuant or Mascot.**
    ///
    /// mzLib's `"trypsin|P"` applies the classic Keil rule — cleave after K/R *except* before
    /// proline — and is the default here because it is what a mass spectrometrist usually means.
    /// mzLib's plain `"trypsin"` cleaves before proline too. That is the **reverse** of the MaxQuant
    /// and Mascot convention, where `Trypsin/P` denotes *ignoring* the proline rule. On serum
    /// albumin the two differ by 7 peptides out of about 200 (195 vs 202, tryptic, 2 missed
    /// cleavages, min length 7) — a small count hiding a large semantic difference, since which
    /// peptides you get changes wherever a K/R precedes a proline.
    pub protease: String,
    /// `"HCD"`/`"CID"` (b and y ions), `"ETD"`, and the rest of mzLib's dissociation types.
    ///
    /// **`"ETD"` returns three series — `c`, `zDot` **and `y`** — not two.** The y ions are
    /// spurious: ETD cleaves N–Cα and yields c/z•, while b/y come from amide cleavage under
    /// vibrational activation, and mzLib's `EThcD` row correctly pairs y *with* b. ETD's does not.
    /// They are about **a third** of every ETD fragment list, so
    /// [`Digest::fragment_count`] over-counts real ETD ions by that much — use
    /// [`Digest::fragments_by_series`] to see the split. Tracked as
    /// [smith-chem-wisc/mzLib#1109](https://github.com/smith-chem-wisc/mzLib/issues/1109).
    pub dissociation: String,
    /// Apply UniProt's annotated modifications.
    ///
    /// **`false` is not a clean control.** It discards UniProt's whole feature table, which also
    /// carries the signal-peptide and propeptide boundaries mzLib digests at — so the *peptide
    /// list* changes too, not only the modifications on it. On albumin, `false` loses
    /// `MKWVTFISLLFLFSSAYS` (1–18) and `WVTFISLLFLFSSAYS` (3–18), both unmodified, both ending
    /// exactly at the signal-peptide cleavage site. Comparing the two runs therefore varies two
    /// things at once. Tracked as
    /// [smith-chem-wisc/pyMzLib#8](https://github.com/smith-chem-wisc/pyMzLib/issues/8).
    ///
    /// Historic wording, kept because it is still true of the modifications themselves: `false`
    /// gives the bare sequence — useful as a
    /// control, and the difference is usually large.
    pub modifications: bool,
    /// Maximum missed cleavage sites per peptide.
    pub missed_cleavages: u32,
    /// Shortest peptide to keep.
    ///
    /// The default of 7 silently discards shorter peptides — roughly a third of a histone digest —
    /// so pass 1 when you mean *every* peptide.
    pub min_length: u32,
    /// Longest peptide to keep. `None` means unbounded.
    pub max_length: Option<u32>,
    /// Maximum modifications considered per peptide.
    ///
    /// Modification isoforms are enumerated combinatorially: histone H3.1 yields 49 bare tryptic
    /// peptides, 2,563 at two modifications and 7,040 at three.
    pub max_modifications: u32,
    /// Maximum modification isoforms per peptide position.
    ///
    /// mzLib's default of 1024 **truncates silently** when it binds — on H3.1 at four modifications
    /// it discards about 30% of the peptidoforms. [`Digest::peptides_at_cap`] reports how many
    /// peptides hit it, so a truncated answer is visible rather than merely short.
    pub max_isoforms: u32,
    /// `"Both"`, `"N"` or `"C"`.
    pub terminus: String,
    /// Seconds to allow. Large proteins with many modification isoforms take longer.
    pub timeout: Option<Duration>,
}

impl Default for FragmentOptions {
    fn default() -> Self {
        Self {
            protease: "trypsin|P".to_owned(),
            dissociation: "ETD".to_owned(),
            modifications: true,
            missed_cleavages: 2,
            min_length: 7,
            max_length: None,
            max_modifications: 2,
            max_isoforms: 1024,
            terminus: "Both".to_owned(),
            timeout: Some(Duration::from_secs(300)),
        }
    }
}

// ------------------------------------------------------------------ validation

/// Whether a string matches UniProtKB's accession grammar.
///
/// <https://www.uniprot.org/help/accession_numbers>, hand-rolled rather than pulled from `regex`:
/// `[OPQ][0-9][A-Z0-9]{3}[0-9]` or `[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}`. Checking it here
/// means a typo costs nothing instead of a network round trip and a puzzling HTTP 400.
fn is_uniprot_accession(candidate: &str) -> bool {
    let bytes = candidate.as_bytes();
    let alnum_upper = |b: u8| b.is_ascii_uppercase() || b.is_ascii_digit();

    // [OPQ][0-9][A-Z0-9]{3}[0-9]
    if bytes.len() == 6
        && matches!(bytes[0], b'O' | b'P' | b'Q')
        && bytes[1].is_ascii_digit()
        && bytes[2..5].iter().copied().all(alnum_upper)
        && bytes[5].is_ascii_digit()
    {
        return true;
    }

    // [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
    if (bytes.len() == 6 || bytes.len() == 10)
        && bytes[0].is_ascii_uppercase()
        && !matches!(bytes[0], b'O' | b'P' | b'Q')
        && bytes[1].is_ascii_digit()
    {
        let groups = (bytes.len() - 2) / 4;
        return (1..=2).contains(&groups)
            && (0..groups).all(|group| {
                let start = 2 + group * 4;
                bytes[start].is_ascii_uppercase()
                    && bytes[start + 1..start + 3].iter().copied().all(alnum_upper)
                    && bytes[start + 3].is_ascii_digit()
            });
    }

    false
}

/// Validate and canonicalise a UniProt accession.
///
/// UniProt's accessions are upper-case and its API is case-sensitive, so the validated canonical
/// form is what gets sent, not the caller's original casing.
fn normalise_accession(accession: &str) -> Result<String> {
    let trimmed = accession.trim();
    if trimmed.is_empty() {
        return Err(MzLibError::Usage(
            "A UniProt accession is required, e.g. 'P02768'.".to_owned(),
        ));
    }

    let canonical = trimmed.to_ascii_uppercase();
    if !is_uniprot_accession(&canonical) {
        return Err(MzLibError::Usage(format!(
            "'{accession}' is not a valid UniProtKB accession. They look like 'P02768' or \
             'A0A0B4J2D5' — see https://www.uniprot.org/help/accession_numbers."
        )));
    }
    Ok(canonical)
}

// ------------------------------------------------------------------ argument assembly

/// The argv for `peptidoform fragments`.
fn build_args(accession: &str, options: &FragmentOptions) -> Result<Vec<String>> {
    let canonical = normalise_accession(accession)?;

    if options.max_isoforms < 1 {
        return Err(MzLibError::Usage(format!(
            "max_isoforms must be at least 1; got {}.",
            options.max_isoforms
        )));
    }

    let mut args = vec![
        "peptidoform".to_owned(),
        "fragments".to_owned(),
        "--accession".to_owned(),
        canonical,
        "--protease".to_owned(),
        options.protease.clone(),
        "--dissociation".to_owned(),
        options.dissociation.clone(),
        "--terminus".to_owned(),
        options.terminus.clone(),
        "--missed-cleavages".to_owned(),
        options.missed_cleavages.to_string(),
        "--min-length".to_owned(),
        options.min_length.to_string(),
        "--max-length".to_owned(),
        options.max_length.unwrap_or(0).to_string(),
        "--max-mods".to_owned(),
        options.max_modifications.to_string(),
        "--max-isoforms".to_owned(),
        options.max_isoforms.to_string(),
    ];

    if !options.modifications {
        args.push("--no-modifications".to_owned());
    }

    Ok(args)
}

// ------------------------------------------------------------------ parsing

/// The wire shape, before the flat census fields are gathered into a [`ModificationCensus`].
#[derive(Debug, Deserialize)]
struct WireDigest {
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    accession: String,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    name: String,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    full_name: String,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    organism: String,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    sequence_length: u32,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    modifications_applied: bool,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    annotated_modification_sites: u32,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    annotated_modifications_loaded: u32,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    uniprot_annotated_features: u32,
    #[serde(default)]
    unresolved_modifications: Vec<String>,
    #[serde(default)]
    uniprot_features_by_type: Vec<FeatureType>,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    protease: String,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    dissociation: String,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    terminus: String,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    max_modifications: u32,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    max_modification_isoforms: u32,
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    peptides_at_isoform_cap: u32,
    #[serde(default)]
    peptides: Vec<Peptide>,
}

/// Turn the `peptidoform fragments` payload into a [`Digest`].
fn parse(data: serde_json::Value) -> Result<Digest> {
    let wire: WireDigest = serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!("peptidoform payload could not be read: {error}"))
    })?;

    Ok(Digest {
        accession: wire.accession,
        name: wire.name,
        full_name: wire.full_name,
        organism: wire.organism,
        sequence_length: wire.sequence_length,
        protease: wire.protease,
        dissociation: wire.dissociation,
        terminus: wire.terminus,
        modifications_applied: wire.modifications_applied,
        max_modifications: wire.max_modifications,
        max_isoforms: wire.max_modification_isoforms,
        peptides_at_cap: wire.peptides_at_isoform_cap,
        modification_census: ModificationCensus {
            sites: wire.annotated_modification_sites,
            applied: wire.annotated_modifications_loaded,
            annotated: wire.uniprot_annotated_features,
            unresolved: wire.unresolved_modifications,
            by_type: wire.uniprot_features_by_type,
        },
        peptides: wire.peptides,
    })
}

// ------------------------------------------------------------------ the public surface

/// Fetch a UniProt entry, digest it, and fragment every peptide, with the lab's defaults.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the accession is not recognised; [`MzLibError::ServiceUnavailable`] if
/// UniProt was unreachable.
pub fn fragments(accession: &str) -> Result<Digest> {
    fragments_with(accession, &FragmentOptions::default())
}

/// [`fragments`], with the digestion and fragmentation stated explicitly.
///
/// Check [`Digest::modification_census`] before trusting a modification count — it reports what was
/// annotated as well as what was applied — and [`Digest::truncated`] before treating the peptide
/// list as exhaustive.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the accession, protease, dissociation type or terminus is not
/// recognised; [`MzLibError::ServiceUnavailable`] if UniProt was unreachable.
pub fn fragments_with(accession: &str, options: &FragmentOptions) -> Result<Digest> {
    let args = build_args(accession, options)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    parse(data)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The recorded albumin digest, the same fixture pyMzLib's offline suite uses.
    const FIXTURE: &str = include_str!("../tests/fixtures/peptidoform_P02768_small.json");

    fn recorded_digest() -> Digest {
        parse(serde_json::from_str(FIXTURE).expect("fixture should be valid JSON"))
            .expect("fixture should parse")
    }

    fn args(options: &FragmentOptions) -> Vec<String> {
        build_args("P02768", options).expect("should assemble")
    }

    fn value_after(args: &[String], flag: &str) -> String {
        let index = args
            .iter()
            .position(|a| a == flag)
            .unwrap_or_else(|| panic!("{flag} should be present in {args:?}"));
        args[index + 1].clone()
    }

    // ---------------------------------------------------------- parsing

    #[test]
    fn digest_and_peptides_parse() {
        let digest = recorded_digest();
        assert_eq!(digest.accession, "P02768");
        assert_eq!(digest.name, "ALBU_HUMAN");
        assert_eq!(digest.organism, "Homo sapiens");
        assert_eq!(digest.sequence_length, 609);
        assert!(!digest.peptides.is_empty());
        assert!(digest.peptides.iter().all(|p| !p.base_sequence.is_empty()));
    }

    #[test]
    fn fragments_parse_into_typed_objects() {
        let digest = recorded_digest();
        let fragment = &digest.peptides[0].fragments[0];
        assert!(!fragment.product_type.is_empty());
        assert!(fragment.fragment_number >= 1);
        assert!(fragment.neutral_mass > 0.0);
    }

    #[test]
    fn fragment_count_sums_across_peptides() {
        let digest = recorded_digest();
        let expected: usize = digest.peptides.iter().map(|p| p.fragments.len()).sum();
        assert_eq!(digest.fragment_count(), expected);
        assert!(digest.fragment_count() > 0);
    }

    #[test]
    fn etd_produces_c_and_z_ions_and_also_y_which_it_should_not() {
        // The dissociation type actually reached the engine — the defaults are opinions, and this
        // is the one that changes which ions come back.
        //
        // It also pins a live mzLib defect, smith-chem-wisc/mzLib#1109: ETD (and ECD) are mapped
        // to { c, y, zDot }, so a third of every ETD fragment list is y ions. ETD cleaves the
        // N–Cα bond and yields c and z• ions; b and y come from amide cleavage under vibrational
        // activation. If the y ions modelled residual activation, b ions would accompany them —
        // which is exactly what mzLib's EThcD row does. y without b has no mechanism.
        //
        // The assertion therefore records what mzLib DOES, not what it should do, and fails the
        // moment the upstream fix lands — which is the point. A test that quietly asserted the
        // wrong behaviour would let the fix look like a regression.
        let digest = recorded_digest();
        assert_eq!(digest.dissociation, "ETD");
        let series: std::collections::HashSet<&str> = digest
            .peptides
            .iter()
            .flat_map(|p| p.fragments.iter())
            .map(|f| f.product_type.as_str())
            .collect();

        assert!(
            series.contains("c"),
            "ETD should produce c ions: {series:?}"
        );
        assert!(
            series.contains("zDot"),
            "ETD should produce z• ions: {series:?}"
        );
        assert!(
            !series.contains("b"),
            "b ions would mean amide cleavage: {series:?}"
        );
        assert!(
            series.contains("y"),
            "mzLib#1109 has ETD emitting y ions; if this now fails, the upstream fix has landed \
             and this test plus the `spurious_etd_y_ions` note should be retired: {series:?}"
        );
    }

    #[test]
    fn the_scale_of_the_spurious_etd_y_ions_is_recorded() {
        // Quantifying it is what makes mzLib#1109 actionable rather than a shrug: a third of the
        // theoretical ion list for ETD data is ions ETD does not make.
        let digest = recorded_digest();
        let total = digest.fragment_count();
        let y_ions = digest
            .peptides
            .iter()
            .flat_map(|p| p.fragments.iter())
            .filter(|f| f.product_type == "y")
            .count();

        assert!(y_ions > 0);
        let fraction = y_ions as f64 / total as f64;
        assert!(
            (0.30..0.40).contains(&fraction),
            "expected roughly a third of ETD products to be y ions, got {fraction:.3} \
             ({y_ions}/{total})"
        );
    }

    // ---------------------------------------------------------- the m/z conversion

    #[test]
    fn mz_uses_the_proton_mass_not_hydrogen() {
        // 0.55 mDa apart — 1.1 ppm at m/z 500, which on an Orbitrap is a match versus a miss.
        let digest = recorded_digest();
        let peptide = &digest.peptides[0];
        let expected = (peptide.monoisotopic_mass + 2.0 * PROTON_MASS) / 2.0;
        assert!((peptide.mz(2).unwrap() - expected).abs() < 1e-9);

        let with_hydrogen = (peptide.monoisotopic_mass + 2.0 * 1.007_825) / 2.0;
        assert!(
            (peptide.mz(2).unwrap() - with_hydrogen).abs() > 1e-4,
            "the hydrogen-atom mass must not have been used"
        );
    }

    #[test]
    fn mz_does_not_double_count_a_fixed_charge() {
        // Trimethyllysine: the recorded mass already carries the charge, so only charge − fixed
        // protons are added. Adding a full complement puts a 2+ peptide half a Thomson high.
        let peptide = Peptide {
            monoisotopic_mass: 1000.0,
            fixed_charges: 1,
            ..recorded_digest().peptides[0].clone()
        };
        let expected = (1000.0 + 1.0 * PROTON_MASS) / 2.0;
        assert!((peptide.mz(2).unwrap() - expected).abs() < 1e-12);

        let naive = (1000.0 + 2.0 * PROTON_MASS) / 2.0;
        assert!((peptide.mz(2).unwrap() - naive).abs() > 0.4);
    }

    #[test]
    fn a_fixed_charge_peptide_is_observable_without_protonation() {
        let peptide = Peptide {
            monoisotopic_mass: 1000.0,
            fixed_charges: 1,
            ..recorded_digest().peptides[0].clone()
        };
        assert!((peptide.mz(1).unwrap() - 1000.0).abs() < 1e-12);
    }

    #[test]
    fn a_charge_below_the_fixed_charge_is_refused() {
        let peptide = Peptide {
            monoisotopic_mass: 1000.0,
            fixed_charges: 2,
            ..recorded_digest().peptides[0].clone()
        };
        let error = peptide.mz(1).unwrap_err();
        assert!(error.to_string().contains("already carries 2 fixed charge"));
    }

    #[test]
    fn mz_rejects_charges_that_are_not_positive() {
        // Python must also reject 1.5, "2", True and None; in Rust those are compile errors, so
        // only the reachable cases remain.
        let peptide = recorded_digest().peptides[0].clone();
        for charge in [0, -1] {
            assert!(peptide.mz(charge).is_err(), "charge {charge}");
        }
    }

    // ---------------------------------------------------------- the census

    #[test]
    fn census_separates_sites_from_modifications() {
        // A histone carries several alternatives at one residue — K9me1/2/3, K9ac are four
        // modifications at one site. Conflating the two made H3.1 look as though 93 annotations
        // had been dropped when they had all been loaded and merely shared residues.
        let census = recorded_digest().modification_census;
        assert_eq!(census.sites, 14);
        assert_eq!(census.applied, 14);
        assert_eq!(census.annotated, 38);
    }

    #[test]
    fn sites_and_modifications_are_not_the_same_number() {
        let census = ModificationCensus {
            sites: 30,
            applied: 93,
            annotated: 93,
            ..Default::default()
        };
        assert_ne!(census.sites, census.applied);
        assert_eq!(census.excluded(), 0);
        assert!(census.explain().contains("All 93"));
        assert!(census.explain().contains("30 residue positions"));
    }

    #[test]
    fn modification_positions_index_the_peptides_own_residues() {
        // mzLib's dictionary reserves slot 1 for the N-terminus, so its keys are one past the
        // residue. A peptide MAR reported position 4 for its arginine — position 4 of a 3-mer.
        let digest = recorded_digest();
        for peptide in &digest.peptides {
            for modification in &peptide.modifications {
                match modification.one_based_residue {
                    Some(position) => {
                        assert!(
                            position >= 1 && position <= peptide.length,
                            "{} on a {}-mer at {position}",
                            modification.id,
                            peptide.length
                        );
                        assert_eq!(modification.terminus, None);
                    }
                    // A terminal modification has no residue index and says which terminus.
                    None => assert!(modification.terminus.is_some()),
                }
            }
        }
    }

    #[test]
    fn unresolved_modification_names_are_reported() {
        // These vanish silently otherwise: on H3.1, seven N6-lactoyllysine sites were dropped
        // while the type summary still said "modified residue … loaded".
        let census = ModificationCensus {
            sites: 10,
            applied: 10,
            annotated: 17,
            unresolved: vec!["N6-lactoyllysine on K".to_owned()],
            by_type: vec![FeatureType {
                feature_type: "modified residue".to_owned(),
                count: 17,
                loaded: true,
            }],
        };
        let explanation = census.explain();
        assert!(
            explanation.contains("N6-lactoyllysine on K"),
            "{explanation}"
        );
        assert!(explanation.contains("absent from its"), "{explanation}");
    }

    #[test]
    fn the_two_exclusion_reasons_are_not_conflated() {
        // Excluded-by-type (no defined mass) and unresolved-by-name (absent from ptmlist) are
        // different failures. Reporting only the first accounted for 3 of 10 exclusions on H3.1
        // and said nothing about the other 7 — worse than saying nothing, because it creates the
        // feeling of having been told.
        let census = ModificationCensus {
            sites: 5,
            applied: 5,
            annotated: 15,
            unresolved: vec!["Something odd on K".to_owned()],
            by_type: vec![
                FeatureType {
                    feature_type: "glycosylation site".to_owned(),
                    count: 7,
                    loaded: false,
                },
                FeatureType {
                    feature_type: "modified residue".to_owned(),
                    count: 8,
                    loaded: true,
                },
            ],
        };
        let explanation = census.explain();
        assert!(
            explanation.contains("7 × glycosylation site"),
            "{explanation}"
        );
        assert!(explanation.contains("Something odd on K"), "{explanation}");
    }

    #[test]
    fn census_explains_what_was_excluded_and_why() {
        let census = recorded_digest().modification_census;
        let explanation = census.explain();
        assert!(explanation.contains("14 of 38"), "{explanation}");
        assert!(
            explanation.contains("24 × glycosylation site"),
            "{explanation}"
        );
        // The reason must not overclaim. 22 of albumin's 24 glycosylation sites are
        // `N-linked (Glc) (glycation) lysine`, which UniProt's own ptmlist defines with
        // CF C6H10O5 and MM 162.052823 — they are excluded by *feature type*, not for want of a
        // mass. The original wording said they had no defined composition, which made the
        // trap-disclosure itself the misinformation. See smith-chem-wisc/mzLib#1112.
        assert!(explanation.contains("feature type"), "{explanation}");
        assert!(explanation.contains("mzLib#1112"), "{explanation}");
        assert_eq!(census.excluded(), 24);
    }

    #[test]
    fn census_explain_stays_at_feature_type_granularity() {
        // The census sees feature *types*, not modification *names*: by_type carries
        // "glycosylation site", never "N-linked (Glc) (glycation) lysine". So explain() must
        // attribute the exclusion to feature type and never present a name-level count — e.g.
        // "22 of 24 are glycation" — or a protein-specific qualifier tally as a census fact.
        // Pins smith-chem-wisc/mzLibRust#2: a documented number at the wrong granularity.
        let census = recorded_digest().modification_census; // albumin: 24 × glycosylation site
        let explanation = census.explain();

        assert!(
            explanation.contains("24 × glycosylation site"),
            "{explanation}"
        );
        assert!(explanation.contains("feature type"), "{explanation}");
        // Facts the census cannot produce must not appear as its output.
        assert!(
            !explanation.contains("N-linked"),
            "name-level modification detail leaked into census output: {explanation}"
        );
        assert!(
            !explanation.contains("22 of"),
            "name-level '22 of 24' count leaked into census output: {explanation}"
        );
        assert!(
            !explanation.contains("14 of the 24"),
            "protein-specific qualifier tally baked into general census output: {explanation}"
        );
    }

    #[test]
    fn census_says_so_when_nothing_was_excluded() {
        let census = ModificationCensus {
            sites: 3,
            applied: 5,
            annotated: 5,
            ..Default::default()
        };
        assert_eq!(census.excluded(), 0);
        assert!(census.explain().starts_with("All 5 annotated"));
    }

    // ---------------------------------------------------------- the silent cap

    #[test]
    fn truncation_is_visible() {
        // mzLib truncates at maxModificationIsoforms and says nothing. A short answer and a
        // truncated answer must not look alike.
        let digest = Digest {
            peptides_at_cap: 12,
            ..recorded_digest()
        };
        assert!(digest.truncated());
    }

    #[test]
    fn no_truncation_when_nothing_hit_the_cap() {
        let digest = recorded_digest();
        assert_eq!(digest.peptides_at_cap, 0);
        assert!(!digest.truncated());
    }

    // ---------------------------------------------------------- the call itself

    #[test]
    fn defaults_are_the_labs_opinion_and_reach_the_bridge() {
        let assembled = args(&FragmentOptions::default());
        assert_eq!(
            assembled[0..2],
            ["peptidoform".to_owned(), "fragments".to_owned()]
        );
        assert_eq!(value_after(&assembled, "--accession"), "P02768");
        assert_eq!(value_after(&assembled, "--protease"), "trypsin|P");
        assert_eq!(value_after(&assembled, "--dissociation"), "ETD");
        assert_eq!(value_after(&assembled, "--terminus"), "Both");
        assert_eq!(value_after(&assembled, "--missed-cleavages"), "2");
        assert_eq!(value_after(&assembled, "--min-length"), "7");
        assert_eq!(value_after(&assembled, "--max-length"), "0");
        assert_eq!(value_after(&assembled, "--max-mods"), "2");
        assert_eq!(value_after(&assembled, "--max-isoforms"), "1024");
        assert!(!assembled.contains(&"--no-modifications".to_owned()));
    }

    #[test]
    fn every_default_is_reachable() {
        // The point is to open the doors, not to hide them.
        let options = FragmentOptions {
            protease: "chymotrypsin (don't cleave before proline)".to_owned(),
            dissociation: "HCD".to_owned(),
            modifications: false,
            missed_cleavages: 0,
            min_length: 1,
            max_length: Some(30),
            max_modifications: 4,
            max_isoforms: 65536,
            terminus: "N".to_owned(),
            timeout: None,
        };
        let assembled = args(&options);
        assert_eq!(
            value_after(&assembled, "--protease"),
            "chymotrypsin (don't cleave before proline)"
        );
        assert_eq!(value_after(&assembled, "--dissociation"), "HCD");
        assert_eq!(value_after(&assembled, "--terminus"), "N");
        assert_eq!(value_after(&assembled, "--missed-cleavages"), "0");
        assert_eq!(value_after(&assembled, "--min-length"), "1");
        assert_eq!(value_after(&assembled, "--max-length"), "30");
        assert_eq!(value_after(&assembled, "--max-mods"), "4");
        assert_eq!(value_after(&assembled, "--max-isoforms"), "65536");
        assert!(assembled.contains(&"--no-modifications".to_owned()));
    }

    #[test]
    fn an_unbounded_max_length_is_sent_as_zero() {
        // The bridge reads 0 as int.MaxValue; None must not become the string "None".
        let assembled = args(&FragmentOptions {
            max_length: None,
            ..Default::default()
        });
        assert_eq!(value_after(&assembled, "--max-length"), "0");
    }

    #[test]
    fn blank_accession_is_refused_before_any_work() {
        for accession in ["", "   "] {
            let error = build_args(accession, &FragmentOptions::default()).unwrap_err();
            assert!(error.to_string().contains("UniProt accession is required"));
        }
    }

    #[test]
    fn malformed_accession_is_refused_without_touching_the_network() {
        for accession in [
            "banana",
            "P0",
            "12345",
            "notanaccession",
            "P0276",
            "PP02768",
        ] {
            assert!(
                build_args(accession, &FragmentOptions::default()).is_err(),
                "{accession} should be rejected"
            );
        }
    }

    #[test]
    fn real_uniprot_accessions_are_accepted_in_any_casing() {
        for accession in ["P02768", "p02768", "  P02768 ", "A0A0B4J2D5", "Q9Y6K9"] {
            assert!(
                normalise_accession(accession).is_ok(),
                "{accession} should be accepted"
            );
        }
        assert_eq!(normalise_accession(" p02768 ").unwrap(), "P02768");
    }

    #[test]
    fn max_isoforms_below_one_is_refused() {
        let error = build_args(
            "P02768",
            &FragmentOptions {
                max_isoforms: 0,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(error.to_string().contains("must be at least 1"));
    }

    #[test]
    fn is_modified_and_modified_peptides_reflect_the_modifications() {
        let digest = recorded_digest();
        let modified = digest.modified_peptides();
        assert!(!modified.is_empty(), "the fixture has a modified peptide");
        assert!(modified.iter().all(|p| p.is_modified()));
        assert_eq!(
            modified.len(),
            digest
                .peptides
                .iter()
                .filter(|p| !p.modifications.is_empty())
                .count()
        );
    }

    #[test]
    fn a_modification_carries_its_mass_and_formal_charge() {
        let digest = recorded_digest();
        let modification = digest
            .modified_peptides()
            .first()
            .and_then(|p| p.modifications.first())
            .expect("the fixture has a modification")
            .clone();
        assert!(!modification.id.is_empty());
        assert!(modification.mass.is_some_and(|mass| mass > 0.0));
        assert_eq!(modification.formal_charge, 0);
    }
}
