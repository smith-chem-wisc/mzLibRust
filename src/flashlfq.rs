//! Label-free quantification with FlashLFQ: quantify a search's peptides across mzML runs.
//!
//! The question this answers is the one a quant workflow actually asks — *given these
//! identifications and these runs, how much of each peptide and protein is in each run?* — in one
//! call:
//!
//! ```no_run
//! # fn main() -> Result<(), mzlib::MzLibError> {
//! use mzlib::flashlfq::{quantify_with, QuantifyOptions, SpectraFile};
//!
//! let result = quantify_with(
//!     "AllPSMs.psmtsv",
//!     &[SpectraFile::from("run_3.mzML"), SpectraFile::from("run_4.mzML")],
//!     &QuantifyOptions { match_between_runs: true, ..Default::default() },
//! )?;
//! println!("{} peptides, {} proteins", result.peptides.len(), result.proteins.len());
//! println!("{} peptides rescued by MBR", result.mbr_rescued_peptide_count());
//! # Ok(())
//! # }
//! ```
//!
//! The whole pipeline is mzLib's: the result file is read by mzLib's `Readers`, turned into
//! FlashLFQ identifications by mzLib's own converter, and quantified by `FlashLfqEngine`.
//! MetaMorpheus is not involved — mzLib does it alone.
//!
//! # Three limits worth knowing before you trust a number
//!
//! **mzML only, for now.** Convert `.raw`/`.d` to mzML first; a non-mzML path is rejected up front.
//!
//! **A protein intensity can be `None`.** FlashLFQ's median-polish protein quant marks a protein
//! NaN when its peptide matrix is degenerate — too few peptides per run, or identical intensities
//! across runs, a real artifact documented in mzLib's own tests. NaN is not valid JSON, so it
//! arrives as `None` — "could not be quantified" — rather than a silently wrong number. A
//! *peptide* intensity, by contrast, is [`f64`] and `0.0` when missing. That asymmetry is the whole
//! reason [`ProteinGroup::intensities`] is `Option<f64>` and [`Peptide::intensities`] is not: the
//! compiler will not let you forget which one you are holding.
//!
//! **For match-between-runs, read the peaks, not the peptides.** The peptide roll-up
//! ([`FlashLfqResults::peptides`], mirroring `QuantifiedPeptides.tsv`) reports far fewer MBR
//! transfers than actually happened — a whole run's transfers can be absent. [`FlashLfqResults::peaks`]
//! is the complete surface, and [`FlashLfqResults::mbr_rescued_peptide_count`] is the number to
//! trust.
//!
//! MBR also **requires a complete, balanced design**: every condition and biological replicate must
//! carry the same set of fractions. And [`QuantifyOptions::mbr_q_value_threshold`] is the FDR
//! control that makes transfers trustworthy — without it a bake-off arm measured roughly 80% false
//! transfers.

use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};

/// One quantified run, mirroring mzLib's `MassSpectrometry.SpectraFileInfo`.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct SpectraFileInfo {
    /// The run's base name (no directory, no extension) — the key used everywhere else here to
    /// look up this run's intensity.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_name: String,
    /// The mzML path as provided.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub full_path: String,
    /// The sample-group label, or `""` if none was given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub condition: String,
    /// The experimental-design biological replicate.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub biological_replicate: u32,
    /// The experimental-design technical replicate.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub technical_replicate: u32,
    /// The experimental-design fraction.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub fraction: u32,
    /// Chromatographic peaks quantified in this run.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub peak_count: u32,
    /// Of those, how many were transferred by match-between-runs — peaks quantified in this run for
    /// a peptide that was never identified *in* it. Zero unless
    /// [`QuantifyOptions::match_between_runs`] was set.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub mbr_peak_count: u32,
}

/// A quantified peptide, mirroring FlashLFQ's `Peptide`.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Peptide {
    /// The full (modified) sequence, as FlashLFQ renders it — the identity FlashLFQ quantifies.
    /// Two different modification states of one base sequence are two peptides.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sequence: String,
    /// The bare amino-acid sequence.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub base_sequence: String,
    /// The protein group(s) this peptide belongs to, `;`-joined.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub protein_groups: String,
    /// Run base name → intensity in that run.
    ///
    /// Missing is `0.0`, **never `None`** — unlike [`ProteinGroup::intensities`]. This roll-up
    /// mirrors FlashLFQ's `QuantifiedPeptides.tsv` and does **not** fully reflect
    /// match-between-runs: many transferred peptides read `0.0` / `"NotDetected"` here. For
    /// MBR-inclusive quantities use [`FlashLfqResults::peaks`].
    ///
    /// **The two surfaces also aggregate differently, so their magnitudes disagree.** Where a
    /// peptide has several peaks in one run, this roll-up reports **one** of them rather than their
    /// sum — `LKEYEAAVEQLK` has peaks of 722,818 and 18,613 in K562_4, and the roll-up says
    /// 722,818. If you follow the advice above and pivot [`FlashLfqResults::peaks`] yourself, decide
    /// deliberately whether to sum or take the maximum, and know that neither will reproduce
    /// `QuantifiedPeptides.tsv` exactly. Counting *presence* is unaffected; summing *intensity* is
    /// not.
    #[serde(default, deserialize_with = "bridge::deserialize_intensities")]
    pub intensities: HashMap<String, f64>,
    /// Run base name → how it was quantified there.
    ///
    /// The values FlashLFQ emits: `"MSMS"` (identified and quantified in this run), `"MBR"`
    /// (transferred by match-between-runs), `"MSMSIdentifiedButNotQuantified"` (an ID here but no
    /// usable peak), `"MSMSAmbiguousPeakfinding"` (more than one peptide fits the peak), and
    /// `"NotDetected"`.
    #[serde(default)]
    pub detection_types: HashMap<String, String>,
}

impl Peptide {
    /// This peptide's intensity in the named run.
    ///
    /// **`0.0` — not `None` — means "not quantified here."** Only *protein* intensities are ever
    /// `None`. Treat `0.0` as missing, not as a measured absence: log-transforming it will mislead.
    /// And note a peptide that FlashLFQ *transferred* into this run by match-between-runs may still
    /// read `0.0` here — see [`FlashLfqResults::peaks`]. Returns `0.0` for a run never provided.
    #[must_use]
    pub fn intensity(&self, file_name: &str) -> f64 {
        self.intensities.get(file_name).copied().unwrap_or(0.0)
    }

    /// How this peptide was quantified in the named run (`"NotDetected"` if it was not).
    #[must_use]
    pub fn detection_type(&self, file_name: &str) -> &str {
        self.detection_types
            .get(file_name)
            .map_or("NotDetected", String::as_str)
    }
}

/// A quantified protein group, mirroring FlashLFQ's `ProteinGroup`.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct ProteinGroup {
    /// The protein group name (accession, or `;`-joined accessions).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub protein_group: String,
    /// The gene name, when the result file carried one.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub gene_name: String,
    /// The organism, when the result file carried one.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub organism: String,
    /// Sample label → protein intensity in that sample.
    ///
    /// **The key is a *sample*, not a file.** FlashLFQ rolls peptides up to proteins per sample,
    /// grouping runs by condition and biological replicate first. [`quantify`] and [`quantify_with`]
    /// give every run its own sample, so the label here is the run base name and the distinction
    /// costs you nothing. It would bite if runs were ever grouped into replicates: the key would
    /// then be the sample (`"condition_replicate"`), and several runs would share one entry. Compare
    /// [`Peptide::intensities`], which is keyed per run in every case — peptides are measured in
    /// files, proteins are resolved across samples.
    ///
    /// **May be `None`**: FlashLFQ's median-polish protein quant emits NaN (returned here as
    /// `None`) when the peptide matrix for the protein is degenerate — too few peptides per run to
    /// resolve, or several runs reporting the same intensity — protecting you from a
    /// fabricated-looking number. Missing, as opposed to *un-resolvable*, is `Some(0.0)`.
    ///
    /// **Calibration: `Some(0.0)` is the overwhelmingly common outcome, not `None`.** On the K562
    /// pair, 847 of 943 protein groups are `0.0` in both runs and only **2** are `None`. Most of
    /// the zeros are the default [`QuantifyOptions::use_shared_peptides_for_protein_quant`] doing
    /// its job: 842 of those 847 groups have no *unique* peptide, so nothing contributes to their
    /// quant. Setting that option to `true` drops the all-zero count from 847 to 34.
    ///
    /// So `Option<f64>` draws the line where the *arithmetic* is degenerate, not where *you* are
    /// stuck. If your question is "which proteins do I have no usable number for", that is
    /// `None` **plus** the zeros — 849 here, not 2.
    #[serde(default)]
    pub intensities: HashMap<String, Option<f64>>,
}

impl ProteinGroup {
    /// This protein's intensity in the named sample.
    ///
    /// `None` means FlashLFQ could not resolve a number (a degenerate peptide matrix);
    /// `Some(0.0)` means simply not measured in this sample.
    ///
    /// The name is a *sample* label — which, for results from [`quantify`] and [`quantify_with`],
    /// is the run base name, since each run is its own sample. See [`ProteinGroup::intensities`].
    #[must_use]
    pub fn intensity(&self, file_name: &str) -> Option<f64> {
        self.intensities
            .get(file_name)
            .copied()
            .unwrap_or(Some(0.0))
    }
}

/// One quantified chromatographic peak, mirroring FlashLFQ's `ChromatographicPeak`.
///
/// This is the surface to use for **match-between-runs**. Unlike the peptide roll-up
/// ([`Peptide::intensities`], which mirrors `QuantifiedPeptides.tsv` and drops most MBR transfers),
/// the peaks fully represent every quantified peak, transferred or not. To build an MBR-inclusive
/// peptide × run matrix, pivot these on `(sequence, file_name)`.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Peak {
    /// The run this peak was measured in (base name).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_name: String,
    /// The full (modified) sequence of the peptide the peak was assigned to.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sequence: String,
    /// The bare amino-acid sequence.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub base_sequence: String,
    /// The peak's intensity (`None` only if FlashLFQ could not integrate it).
    #[serde(default)]
    pub intensity: Option<f64>,
    /// `"MSMS"`, `"MBR"`, `"MSMSIdentifiedButNotQuantified"`, `"MSMSAmbiguousPeakfinding"` — the
    /// same vocabulary as [`Peptide::detection_types`]. Filter on `"MBR"` to see exactly the
    /// transferred peaks.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub detection_type: String,
    /// Apex retention time in minutes, or `None` if the peak has no apex.
    #[serde(default)]
    pub retention_time: Option<f64>,
    /// How many peptides could explain this peak. Greater than 1 means it is ambiguous and its
    /// intensity should be treated with care.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub num_identifications: u32,
    /// The protein group(s) the assigned identification(s) belong to, `;`-joined.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub protein_groups: String,
}

impl Peak {
    /// Whether this peak was transferred by match-between-runs.
    #[must_use]
    pub fn is_mbr(&self) -> bool {
        self.detection_type == "MBR"
    }
}

/// The FlashLFQ parameters actually used, echoed back with their mzLib names.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct FlashLfqParameters {
    /// Whether intensities were normalized across runs.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub normalize: bool,
    /// Mass tolerance for peak-finding, in ppm.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub ppm_tolerance: f64,
    /// Mass tolerance for isotope-envelope matching, in ppm.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub isotope_ppm_tolerance: f64,
    /// Whether peak intensities were integrated rather than taken at the apex.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub integrate: bool,
    /// Whether match-between-runs was on.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub match_between_runs: bool,
    /// Mass tolerance for MBR transfers, in ppm.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub mbr_ppm_tolerance: f64,
    /// The q-value cutoff below which an MBR transfer was accepted.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub mbr_q_value_threshold: f64,
    /// Whether peptides shared between protein groups contributed to protein quant.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub use_shared_peptides_for_protein_quant: bool,
    /// Whether FlashLFQ's Bayesian protein-fold-change engine was run.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub bayesian_protein_quant: bool,
    /// Worker threads used; `-1` lets FlashLFQ choose.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub max_threads: i32,
}

/// The result of a quantification run, mirroring mzLib's `FlashLfqResults`.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct FlashLfqResults {
    /// The absolute path of the PSM result file that was quantified.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub psm_file: String,
    /// How many identifications were read from it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub identification_count: u32,
    /// The FlashLFQ parameters actually used.
    pub parameters: FlashLfqParameters,
    /// One entry per quantified run.
    #[serde(default)]
    pub spectra_files: Vec<SpectraFileInfo>,
    /// One entry per quantified modified sequence.
    ///
    /// Mirrors `QuantifiedPeptides.tsv`; does **not** fully carry MBR — use [`Self::peaks`].
    #[serde(default)]
    pub peptides: Vec<Peptide>,
    /// One entry per quantified protein group.
    #[serde(default)]
    pub proteins: Vec<ProteinGroup>,
    /// Every quantified peak across all runs (mirrors `QuantifiedPeaks.tsv`) — the complete
    /// surface, and the one to use for match-between-runs.
    #[serde(default)]
    pub peaks: Vec<Peak>,
    /// Where the FlashLFQ TSVs were written, or `None` if none were.
    #[serde(default)]
    pub output_directory: Option<String>,
}

impl FlashLfqResults {
    /// Total match-between-runs peaks (transfers) across every run.
    ///
    /// This counts transferred **peaks**, not distinct peptides: one peptide rescued in two runs is
    /// two peaks here. For "how many peptides did MBR rescue", use
    /// [`Self::mbr_rescued_peptide_count`]. Either way, do not count MBR from the peptide roll-up
    /// ([`Peptide::detection_types`]) — it under-counts. Zero unless MBR was on.
    #[must_use]
    pub fn mbr_peak_count(&self) -> u32 {
        self.spectra_files.iter().map(|f| f.mbr_peak_count).sum()
    }

    /// Exactly the peaks transferred by match-between-runs.
    #[must_use]
    pub fn mbr_peaks(&self) -> Vec<&Peak> {
        self.peaks.iter().filter(|peak| peak.is_mbr()).collect()
    }

    /// Distinct peptides quantified in at least one run only by match-between-runs.
    ///
    /// Exactly: **the number of distinct `sequence` values among peaks whose `detection_type` is
    /// `"MBR"`.** Stated in code terms because the prose version — "peptides quantified in at least
    /// one run *only* by match-between-runs" — is subtly different, and on real data the two
    /// diverge.
    ///
    /// The gap is peptides that have **both** an MBR peak *and* a zero-intensity `MSMS` peak in the
    /// same run: FlashLFQ identified the peptide there but could not integrate a peak from the
    /// identification, so it transferred one as well. Those were identified in that run, so under
    /// the strict "never identified there" reading they are not rescues. On the K562 pair that is 5
    /// peptides — this method returns 140, the strict count is 135.
    ///
    /// Neither number is wrong; they answer different questions. If you want the strict one:
    ///
    /// ```no_run
    /// # use mzlib::flashlfq::FlashLfqResults;
    /// # fn strict(result: &FlashLfqResults) -> usize {
    /// use std::collections::HashSet;
    /// let identified: HashSet<(&str, &str)> = result
    ///     .peaks
    ///     .iter()
    ///     .filter(|p| p.detection_type == "MSMS")
    ///     .map(|p| (p.file_name.as_str(), p.sequence.as_str()))
    ///     .collect();
    /// result
    ///     .mbr_peaks()
    ///     .iter()
    ///     .filter(|p| !identified.contains(&(p.file_name.as_str(), p.sequence.as_str())))
    ///     .map(|p| p.sequence.as_str())
    ///     .collect::<HashSet<_>>()
    ///     .len()
    /// # }
    /// ```
    ///
    /// Do not read `mbr_rescued_peptide_count() == mbr_peak_count()` as reassurance that nothing
    /// was double-counted. On this data both are 140, and they coincide only because every MBR peak
    /// happened to carry a distinct sequence.
    #[must_use]
    pub fn mbr_rescued_peptide_count(&self) -> usize {
        self.mbr_peaks()
            .iter()
            .map(|peak| peak.sequence.as_str())
            .collect::<HashSet<_>>()
            .len()
    }
}

/// One run, with as much of the experimental design as you care to state.
///
/// A bare path is the common case; the design fields are what MBR needs.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SpectraFile {
    /// The mzML path.
    pub path: String,
    /// The sample-group label. `None` becomes a blank condition.
    pub condition: Option<String>,
    /// The biological replicate. `None` lets the bridge default it, as MetaMorpheus does with no
    /// experimental-design file: each file its own biological replicate.
    pub biological_replicate: Option<u32>,
    /// The technical replicate. `None` defaults to 0.
    pub technical_replicate: Option<u32>,
    /// The fraction. `None` defaults to 0.
    pub fraction: Option<u32>,
}

impl From<&str> for SpectraFile {
    fn from(path: &str) -> Self {
        Self {
            path: path.to_owned(),
            ..Self::default()
        }
    }
}

impl From<String> for SpectraFile {
    fn from(path: String) -> Self {
        Self {
            path,
            ..Self::default()
        }
    }
}

/// How a quantification is run. Every default is FlashLFQ's own.
#[derive(Debug, Clone)]
pub struct QuantifyOptions {
    /// Normalize intensities across runs (FlashLFQ `Normalize`).
    pub normalize: bool,
    /// Mass tolerance for peak-finding, in ppm (`PpmTolerance`).
    pub ppm_tolerance: f64,
    /// Mass tolerance for isotope-envelope matching, in ppm.
    pub isotope_ppm_tolerance: f64,
    /// Integrate peak intensities rather than taking the apex. FlashLFQ recommends leaving this
    /// off.
    pub integrate: bool,
    /// Quantify a peptide in a run where it was not identified, by transferring the identification
    /// from a run where it was (`MatchBetweenRuns`).
    ///
    /// The reason to reach for FlashLFQ; off by default because it makes assumptions worth opting
    /// into. **Requires a complete, balanced design**: every condition and biological replicate
    /// must carry the same set of fractions — a missing replicate or fraction breaks the
    /// complementarity MBR relies on. Read the transferred peaks from [`FlashLfqResults::peaks`],
    /// not the peptide table.
    pub match_between_runs: bool,
    /// Mass tolerance for MBR transfers, in ppm.
    pub mbr_ppm_tolerance: f64,
    /// The q-value cutoff below which an MBR transfer is accepted.
    ///
    /// This is the FDR control that makes transfers trustworthy. Without it, a bake-off arm
    /// measured roughly 80% false transfers.
    pub mbr_q_value_threshold: f64,
    /// Let peptides shared between protein groups contribute to protein quant.
    pub use_shared_peptides_for_protein_quant: bool,
    /// Run FlashLFQ's Bayesian protein-fold-change engine.
    pub bayesian_protein_quant: bool,
    /// Filter identifications on PEP q-value rather than q-value.
    pub use_pep_q_value: bool,
    /// Worker threads; `-1` lets FlashLFQ choose.
    ///
    /// **This is not only a performance knob — it changes results.** With `-1`, FlashLFQ's peptide
    /// roll-up nondeterministically drops some MBR intensities, so peptide and protein numbers vary
    /// between runs on byte-identical inputs. On the K562 pair, 6 peptides flip between `0.0` and a
    /// real intensity, which flips a borderline protein group between `None` and a number: it came
    /// back unquantifiable in 5 of 6 runs and quantified in the 6th. The `peaks` are stable
    /// throughout — only the roll-up wobbles.
    ///
    /// **Set `max_threads: 1` for anything you intend to publish or reproduce.** Tracked as
    /// [smith-chem-wisc/mzLib#1111](https://github.com/smith-chem-wisc/mzLib/issues/1111).
    pub max_threads: i32,
    /// If given, FlashLFQ also writes `QuantifiedPeaks.tsv`, `QuantifiedPeptides.tsv` and
    /// `QuantifiedProteins.tsv` there.
    pub output_directory: Option<String>,
    /// Seconds to allow. Large experiments legitimately take a while; `None` waits indefinitely.
    pub timeout: Option<Duration>,
}

impl Default for QuantifyOptions {
    fn default() -> Self {
        Self {
            normalize: false,
            ppm_tolerance: 10.0,
            isotope_ppm_tolerance: 5.0,
            integrate: false,
            match_between_runs: false,
            mbr_ppm_tolerance: 10.0,
            mbr_q_value_threshold: 0.05,
            use_shared_peptides_for_protein_quant: false,
            bayesian_protein_quant: false,
            use_pep_q_value: false,
            max_threads: -1,
            output_directory: None,
            timeout: None,
        }
    }
}

// ------------------------------------------------------------------ stdin rendering

/// Render the spectra argument as the bridge's tab-separated stdin lines.
///
/// The wire format is one run per line: `path[⇥condition[⇥biorep[⇥techrep[⇥fraction]]]]`. Omitted
/// trailing design fields are simply not written; the bridge applies the same defaults MetaMorpheus
/// does with no experimental-design file (blank condition, each run its own biological replicate,
/// fraction 0, technical replicate 0).
///
/// stdin rather than argv because a real experiment has many runs and argv has a hard size ceiling.
fn spectra_stdin(spectra: &[SpectraFile]) -> Result<String> {
    if spectra.is_empty() {
        return Err(MzLibError::Usage(
            "At least one spectra file is required.".to_owned(),
        ));
    }

    let mut lines = Vec::with_capacity(spectra.len());
    for (index, file) in spectra.iter().enumerate() {
        if file.path.trim().is_empty() {
            return Err(MzLibError::Usage(format!("spectra[{index}] has no path.")));
        }
        if file.path.contains('\t') || file.path.contains('\n') {
            return Err(MzLibError::Usage(format!(
                "spectra[{index}] path may not contain a tab or newline."
            )));
        }

        let mut fields = vec![
            file.path.clone(),
            file.condition.clone().unwrap_or_default(),
            file.biological_replicate
                .map(|v| v.to_string())
                .unwrap_or_default(),
            file.technical_replicate
                .map(|v| v.to_string())
                .unwrap_or_default(),
            file.fraction.map(|v| v.to_string()).unwrap_or_default(),
        ];

        // Drop trailing empties so a bare path stays a bare path on the wire. A blank field in the
        // MIDDLE is kept, because dropping it would shift every later field left by one and
        // silently reassign a fraction to a replicate.
        while fields.len() > 1 && fields.last().is_some_and(String::is_empty) {
            fields.pop();
        }
        lines.push(fields.join("\t"));
    }

    Ok(format!("{}\n", lines.join("\n")))
}

// ------------------------------------------------------------------ argument assembly

/// The argv for `quant flashlfq`.
fn build_args(psms: &str, options: &QuantifyOptions) -> Result<Vec<String>> {
    if psms.trim().is_empty() {
        return Err(MzLibError::Usage(
            "A PSM result file path is required, e.g. 'AllPSMs.psmtsv'.".to_owned(),
        ));
    }

    let mut args = vec![
        "quant".to_owned(),
        "flashlfq".to_owned(),
        "--psms".to_owned(),
        psms.to_owned(),
    ];

    if options.normalize {
        args.push("--normalize".to_owned());
    }
    args.push("--ppm".to_owned());
    args.push(bridge::format_number(
        options.ppm_tolerance,
        "ppm_tolerance",
    )?);
    args.push("--isotope-ppm".to_owned());
    args.push(bridge::format_number(
        options.isotope_ppm_tolerance,
        "isotope_ppm_tolerance",
    )?);
    if options.integrate {
        args.push("--integrate".to_owned());
    }
    if options.match_between_runs {
        args.push("--mbr".to_owned());
    }
    args.push("--mbr-ppm".to_owned());
    args.push(bridge::format_number(
        options.mbr_ppm_tolerance,
        "mbr_ppm_tolerance",
    )?);
    args.push("--mbr-q".to_owned());
    args.push(bridge::format_number(
        options.mbr_q_value_threshold,
        "mbr_q_value_threshold",
    )?);
    if options.use_shared_peptides_for_protein_quant {
        args.push("--shared-peptides".to_owned());
    }
    if options.bayesian_protein_quant {
        args.push("--bayesian".to_owned());
    }
    if options.use_pep_q_value {
        args.push("--use-pep-q".to_owned());
    }
    args.push("--threads".to_owned());
    args.push(options.max_threads.to_string());

    if let Some(directory) = &options.output_directory {
        if directory.trim().is_empty() {
            return Err(MzLibError::Usage(
                "output_directory must be a non-empty path or None.".to_owned(),
            ));
        }
        args.push("--out".to_owned());
        args.push(directory.clone());
    }

    Ok(args)
}

/// Turn the `quant flashlfq` payload into typed results.
fn parse(data: serde_json::Value) -> Result<FlashLfqResults> {
    serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!("FlashLFQ payload could not be read: {error}"))
    })
}

// ------------------------------------------------------------------ the public surface

/// Quantify a search's peptides across mzML runs with FlashLFQ, with FlashLFQ's own defaults.
///
/// `psms` is a PSM result file — a MetaMorpheus `.psmtsv` gives the full field set (q-values,
/// scores); an MSFragger result file also works. Every run named in it must have a matching mzML in
/// `spectra`; FlashLFQ matches identifications to runs by base file name, so base names must be
/// unique.
///
/// # Errors
///
/// [`MzLibError::Usage`] if an argument is malformed, a run is not mzML, an mzML is missing, or the
/// PSM file names a run with no mzML provided; [`MzLibError::Bridge`] if FlashLFQ itself failed.
pub fn quantify(psms: impl AsRef<Path>, spectra: &[SpectraFile]) -> Result<FlashLfqResults> {
    quantify_with(psms, spectra, &QuantifyOptions::default())
}

/// [`quantify`], with the FlashLFQ parameters stated explicitly.
///
/// # Errors
///
/// As [`quantify`].
pub fn quantify_with(
    psms: impl AsRef<Path>,
    spectra: &[SpectraFile],
    options: &QuantifyOptions,
) -> Result<FlashLfqResults> {
    let psms = psms.as_ref().to_string_lossy().into_owned();
    // The design is rendered BEFORE the arguments are assembled, so an empty or malformed spectra
    // list fails without the caller's other arguments having to be valid first.
    let stdin = spectra_stdin(spectra)?;
    let args = build_args(&psms, options)?;
    let data = bridge::invoke(&args, Some(&stdin), options.timeout)?;
    parse(data)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The recorded quantification payload, the same fixture pyMzLib's offline suite uses.
    const FIXTURE: &str = include_str!("../tests/fixtures/flashlfq_small.json");

    fn recorded() -> FlashLfqResults {
        parse(serde_json::from_str(FIXTURE).expect("fixture should be valid JSON"))
            .expect("fixture should parse")
    }

    fn two_runs() -> Vec<SpectraFile> {
        vec![
            SpectraFile::from("run_3.mzML"),
            SpectraFile::from("run_4.mzML"),
        ]
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
    fn results_parse_into_typed_objects() {
        let result = recorded();
        assert_eq!(result.identification_count, 4);
        assert_eq!(result.peptides.len(), 2);
        assert_eq!(result.proteins.len(), 2);
        assert_eq!(result.spectra_files.len(), 2);
        assert_eq!(result.parameters.ppm_tolerance, 10.0);
        assert!(result.parameters.match_between_runs);
        assert_eq!(result.output_directory, None);
    }

    #[test]
    fn peptide_intensity_and_detection_lookup() {
        let peptide = &recorded().peptides[0];
        assert_eq!(peptide.sequence, "PEPTIDEK");
        assert_eq!(peptide.protein_groups, "P12345");
        assert_eq!(peptide.intensity("run_3"), 1000.0);
        assert_eq!(peptide.detection_type("run_4"), "MBR");
        // A run that was never provided reads as not detected, not a panic.
        assert_eq!(peptide.detection_type("never_seen"), "NotDetected");
        // Missing peptide intensity is 0.0, never None (None is proteins-only).
        assert_eq!(peptide.intensity("never_seen"), 0.0);
    }

    #[test]
    fn a_null_peptide_intensity_reads_as_zero_not_as_unresolvable() {
        // The fixture carries `"run_4": null` for ACDEFR. pyMzLib documents peptide intensities as
        // "0.0 when missing, never None" but stores the wire value verbatim, so its own fixture
        // makes `intensity("run_4")` return None — the documented invariant is not enforced and no
        // Python test covers this case.
        //
        // Rust cannot be ambiguous here: `intensities` is HashMap<String, f64>, so the null is
        // resolved at the boundary and the documented invariant is true by construction. The
        // distinction that matters — "could not be resolved" — belongs to proteins alone, and is
        // Option<f64> there.
        let acdefr = recorded()
            .peptides
            .into_iter()
            .find(|p| p.base_sequence == "ACDEFR")
            .expect("the fixture has ACDEFR");
        assert_eq!(acdefr.intensity("run_4"), 0.0);
        assert_eq!(acdefr.detection_type("run_4"), "NotDetected");
    }

    #[test]
    fn nan_protein_intensity_arrives_as_none() {
        // FlashLFQ's median-polish protein quant emits NaN for a degenerate peptide matrix. NaN is
        // not valid JSON, so it crosses as null and must stay distinguishable from a real zero.
        let result = recorded();
        let unquantified = result
            .proteins
            .iter()
            .find(|g| g.protein_group == "P67890")
            .expect("the fixture has P67890");
        assert_eq!(unquantified.intensity("run_3"), None);
        assert_eq!(unquantified.intensity("run_4"), None);

        let quantified = result
            .proteins
            .iter()
            .find(|g| g.protein_group == "P12345")
            .expect("the fixture has P12345");
        assert_eq!(quantified.intensity("run_3"), Some(1000.0));
        // Not measured in a run that was never provided is zero, not unresolvable.
        assert_eq!(quantified.intensity("never_seen"), Some(0.0));
    }

    #[test]
    fn the_none_and_zero_distinction_is_carried_by_the_types() {
        // The selling point, asserted rather than assumed: you cannot read a protein intensity
        // without confronting the possibility that FlashLFQ could not resolve one.
        let result = recorded();
        let peptide_total: f64 = result.peptides.iter().map(|p| p.intensity("run_3")).sum();
        assert!(peptide_total > 0.0);

        let resolved: Vec<f64> = result
            .proteins
            .iter()
            .filter_map(|g| g.intensity("run_3"))
            .collect();
        assert_eq!(resolved.len(), 1, "one of two proteins resolved in run_3");
    }

    #[test]
    fn mbr_peak_count_sums_across_runs() {
        assert_eq!(recorded().mbr_peak_count(), 1);
    }

    #[test]
    fn spectra_file_fields_parse() {
        let result = recorded();
        let run_4 = &result.spectra_files[1];
        assert_eq!(run_4.file_name, "run_4");
        assert_eq!(run_4.condition, "treated");
        assert_eq!(run_4.biological_replicate, 1);
        assert_eq!(run_4.peak_count, 2);
    }

    // ---------------------------------------------------------- peaks (the MBR surface)

    #[test]
    fn peaks_parse_into_typed_objects() {
        let result = recorded();
        assert_eq!(result.peaks.len(), 3);
        let peak = &result.peaks[0];
        assert_eq!(peak.file_name, "run_3");
        assert_eq!(peak.sequence, "PEPTIDEK");
        assert_eq!(peak.intensity, Some(1000.0));
        assert_eq!(peak.detection_type, "MSMS");
        assert_eq!(peak.retention_time, Some(30.1));
        assert!(!peak.is_mbr());
    }

    #[test]
    fn mbr_peaks_are_the_transferred_ones() {
        let result = recorded();
        let mbr = result.mbr_peaks();
        assert_eq!(mbr.len(), 1);
        assert!(mbr[0].is_mbr());
        assert_eq!(mbr[0].file_name, "run_4");
        assert_eq!(mbr[0].sequence, "PEPTIDEK");
    }

    #[test]
    fn peaks_carry_mbr_the_peptide_table_drops() {
        // The finding the pyMzLib bake-off surfaced: an MBR transfer visible in the peaks need not
        // appear as "MBR" in the peptide roll-up. The peaks are the authoritative surface.
        let result = recorded();
        assert!(result.peaks.iter().any(Peak::is_mbr));
    }

    #[test]
    fn mbr_rescued_peptide_count_is_distinct_sequences() {
        // One MBR peak in the fixture (PEPTIDEK in run_4) -> exactly one distinct rescued peptide.
        assert_eq!(recorded().mbr_rescued_peptide_count(), 1);
    }

    #[test]
    fn rescued_peptide_count_counts_sequences_not_peaks() {
        // One peptide rescued in two runs is two peaks but one rescued peptide.
        let mut result = recorded();
        let mut second = result.peaks[1].clone();
        second.file_name = "run_5".to_owned();
        result.peaks.push(second);
        assert_eq!(result.mbr_peaks().len(), 2);
        assert_eq!(result.mbr_rescued_peptide_count(), 1);
    }

    // ---------------------------------------------------------- argument assembly

    #[test]
    fn quantify_passes_expected_args_and_stdin() {
        let options = QuantifyOptions {
            match_between_runs: true,
            ppm_tolerance: 7.5,
            ..Default::default()
        };
        let args = build_args("AllPSMs.psmtsv", &options).unwrap();

        assert_eq!(args[0..2], ["quant".to_owned(), "flashlfq".to_owned()]);
        assert_eq!(value_after(&args, "--psms"), "AllPSMs.psmtsv");
        assert_eq!(value_after(&args, "--ppm"), "7.5");
        assert_eq!(value_after(&args, "--isotope-ppm"), "5.0");
        assert_eq!(value_after(&args, "--mbr-ppm"), "10.0");
        assert_eq!(value_after(&args, "--mbr-q"), "0.05");
        assert_eq!(value_after(&args, "--threads"), "-1");
        assert!(args.contains(&"--mbr".to_owned()));

        let stdin = spectra_stdin(&two_runs()).unwrap();
        assert_eq!(stdin, "run_3.mzML\nrun_4.mzML\n");
    }

    #[test]
    fn flags_absent_when_false() {
        let args = build_args("AllPSMs.psmtsv", &QuantifyOptions::default()).unwrap();
        for flag in [
            "--normalize",
            "--integrate",
            "--mbr",
            "--shared-peptides",
            "--bayesian",
            "--use-pep-q",
            "--out",
        ] {
            assert!(!args.contains(&flag.to_owned()), "{flag} should be absent");
        }
    }

    #[test]
    fn every_flag_is_reachable() {
        let options = QuantifyOptions {
            normalize: true,
            integrate: true,
            match_between_runs: true,
            use_shared_peptides_for_protein_quant: true,
            bayesian_protein_quant: true,
            use_pep_q_value: true,
            max_threads: 4,
            ..Default::default()
        };
        let args = build_args("AllPSMs.psmtsv", &options).unwrap();
        for flag in [
            "--normalize",
            "--integrate",
            "--mbr",
            "--shared-peptides",
            "--bayesian",
            "--use-pep-q",
        ] {
            assert!(args.contains(&flag.to_owned()), "{flag} should be present");
        }
        assert_eq!(value_after(&args, "--threads"), "4");
    }

    #[test]
    fn output_directory_becomes_out_flag() {
        let options = QuantifyOptions {
            output_directory: Some("quant-out".to_owned()),
            ..Default::default()
        };
        let args = build_args("AllPSMs.psmtsv", &options).unwrap();
        assert_eq!(value_after(&args, "--out"), "quant-out");
    }

    #[test]
    fn a_blank_output_directory_is_refused() {
        let options = QuantifyOptions {
            output_directory: Some("   ".to_owned()),
            ..Default::default()
        };
        assert!(build_args("AllPSMs.psmtsv", &options).is_err());
    }

    // ---------------------------------------------------------- stdin rendering

    #[test]
    fn spectra_stdin_bare_paths_are_bare() {
        assert_eq!(
            spectra_stdin(&two_runs()).unwrap(),
            "run_3.mzML\nrun_4.mzML\n"
        );
    }

    #[test]
    fn spectra_stdin_renders_design_fields() {
        let files = vec![SpectraFile {
            path: "a.mzML".to_owned(),
            condition: Some("treated".to_owned()),
            biological_replicate: Some(1),
            technical_replicate: Some(0),
            fraction: Some(2),
        }];
        assert_eq!(spectra_stdin(&files).unwrap(), "a.mzML\ttreated\t1\t0\t2\n");
    }

    #[test]
    fn spectra_stdin_trims_trailing_empty_fields() {
        let files = vec![SpectraFile {
            path: "a.mzML".to_owned(),
            condition: Some("treated".to_owned()),
            ..Default::default()
        }];
        assert_eq!(spectra_stdin(&files).unwrap(), "a.mzML\ttreated\n");
    }

    #[test]
    fn spectra_stdin_keeps_empty_middle_field() {
        // Dropping a blank field in the middle would shift every later field left by one, silently
        // reading a fraction as a replicate.
        let files = vec![SpectraFile {
            path: "a.mzML".to_owned(),
            condition: None,
            biological_replicate: Some(3),
            ..Default::default()
        }];
        assert_eq!(spectra_stdin(&files).unwrap(), "a.mzML\t\t3\n");
    }

    // ---------------------------------------------------------- validation

    #[test]
    fn empty_spectra_raises_before_invoke() {
        let error = spectra_stdin(&[]).unwrap_err();
        assert!(error.to_string().contains("At least one spectra file"));
    }

    #[test]
    fn blank_psms_raises() {
        for psms in ["", "   "] {
            let error = build_args(psms, &QuantifyOptions::default()).unwrap_err();
            assert!(error
                .to_string()
                .contains("PSM result file path is required"));
        }
    }

    #[test]
    fn a_spectra_entry_without_a_path_raises() {
        let files = vec![SpectraFile {
            path: "  ".to_owned(),
            ..Default::default()
        }];
        assert!(spectra_stdin(&files).is_err());
    }

    #[test]
    fn tab_in_path_raises() {
        // A tab in a path would be read as a field separator and silently become a condition.
        let files = vec![SpectraFile {
            path: "has\ttab.mzML".to_owned(),
            ..Default::default()
        }];
        let error = spectra_stdin(&files).unwrap_err();
        assert!(error.to_string().contains("tab or newline"));
    }

    #[test]
    fn newline_in_path_raises() {
        let files = vec![SpectraFile {
            path: "has\nnewline.mzML".to_owned(),
            ..Default::default()
        }];
        assert!(spectra_stdin(&files).is_err());
    }

    #[test]
    fn non_finite_tolerance_raises() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let options = QuantifyOptions {
                ppm_tolerance: bad,
                ..Default::default()
            };
            let error = build_args("AllPSMs.psmtsv", &options).unwrap_err();
            assert!(error.to_string().contains("finite number"), "{bad}");
        }
    }

    #[test]
    fn tolerances_are_formatted_invariantly() {
        // The bridge parses with InvariantCulture; a comma-decimal locale must never reach it.
        let options = QuantifyOptions {
            ppm_tolerance: 12.5,
            mbr_q_value_threshold: 0.01,
            ..Default::default()
        };
        let args = build_args("AllPSMs.psmtsv", &options).unwrap();
        assert_eq!(value_after(&args, "--ppm"), "12.5");
        assert_eq!(value_after(&args, "--mbr-q"), "0.01");
    }

    #[test]
    fn a_negative_replicate_is_impossible_by_construction() {
        // Python must reject these at runtime; here they are compile errors, because the design
        // fields are u32. The test documents the guarantee rather than defending it.
        let files = vec![SpectraFile {
            path: "a.mzML".to_owned(),
            biological_replicate: Some(0),
            ..Default::default()
        }];
        assert_eq!(spectra_stdin(&files).unwrap(), "a.mzML\t\t0\n");
    }
}
