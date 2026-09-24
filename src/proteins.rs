//! Protein databases: what each protein *is*, which gene it belongs to, and whether a peptide
//! identifies it.
//!
//! Three questions a search result cannot answer on its own, each answered by mzLib from the
//! protein database you searched.
//!
//! **What is this accession?** [`read_with`] loads UniProt XML or FASTA and returns one row per
//! protein — organism, NCBI taxonomy id, gene names, length, monoisotopic mass — and, on request,
//! its Gene Ontology terms and Ensembl gene links as long tables:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::proteins::{read_with, ProteinReadOptions, ProteinTable};
//!
//! let db = read_with(
//!     &["human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"],
//!     &ProteinReadOptions {
//!         contaminants: vec!["contaminants.fasta".into()],
//!         tables: ProteinTable::ALL.to_vec(),
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!(db.taxonomy()?["Q9Z0X1"].as_deref(), Some("10090"));
//! assert_eq!(db.organisms()?["P02769"].as_deref(), Some("Bos taurus"));
//! let go = db.go_terms.as_ref().expect("go_terms was asked for");
//! assert_eq!(go.columns.strings("term_name")?[0].as_deref(), Some("cytoplasm"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Which gene is it, reproducibly?** [`resolve_genes_with`] resolves every protein to a stable
//! Ensembl gene id, counted against **a gene set you supply and pin** — an Ensembl GTF for one
//! release — and says, per protein, *how* it resolved (`resolved`, `multi_gene`,
//! `off_primary_only`, `not_in_source`, …). Every row carries the sha256 of the database, the GTF
//! and the optional cross-reference table it was computed from:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::proteins::{resolve_genes_with, GeneResolveOptions};
//!
//! let genes = resolve_genes_with(
//!     &["human_subset.xml", "human_extra.fasta"],
//!     &GeneResolveOptions {
//!         gtf: Some("Homo_sapiens.GRCh38.116.gtf".into()),
//!         xref: Some("Homo_sapiens.GRCh38.116.uniprot.tsv".into()),
//!         contaminants: vec!["contaminants.fasta".into()],
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!(genes.gene_set.release.as_deref(), Some("116"));
//! assert_eq!(genes.gene_set.genome_build.as_deref(), Some("GRCh38.p14"));
//! assert_eq!(genes.outcome_counts["resolved"], 1);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Does this peptide identify one protein?** [`classify_peptides_with`] classifies each peptide
//! as `Unique`, `SharedWithinGene`, `SharedAcrossGenes` or `NotInDatabase`, treating **I and L as
//! the same residue**, because a mass spectrometer cannot tell them apart:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::proteins::{classify_peptides_with, ClassifyOptions};
//!
//! let calls = classify_peptides_with(
//!     &["VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"],
//!     &["human_subset.xml", "human_extra.fasta"],
//!     &ClassifyOptions { contaminants: vec!["contaminants.fasta".into()], ..Default::default() },
//! )?;
//! let sharing: Vec<Option<String>> = calls.columns.strings("sharing")?;
//! assert_eq!(sharing[3].as_deref(), Some("SharedAcrossGenes"));   // human and bovine albumin
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! These examples replay the recordings in `tests/fixtures/`, made over the small databases in
//! `tests/fixtures/proteins/` — the same ones pyMzLib's examples replay.
//!
//! All three take **one database or many**. A list is read in one bridge call, in order;
//! `threads` says how many are read at once (default 1), and the answer is identical at any value.
//! Mark contaminant databases with `contaminants` — it changes answers: a contaminant is never
//! mapped to a gene, and a peptide it shares with a target is shared.
//!
//! # What these functions will not tell you, and say so instead
//!
//! - **A FASTA knows no GO terms and no Ensembl genes.** Its headers carry organism (`OS=`),
//!   taxonomy (`OX=`) and gene name (`GN=`) and nothing else, so it contributes no `go_terms` or
//!   `ensembl_genes` rows, and [`resolve_genes_with`] calls every FASTA protein `not_in_source`.
//!   That is the format's silence, not a biological absence, and each result says so in
//!   [`DatabaseFile::absent_fields`] and its caveats. Read the UniProt XML of the same proteome.
//! - **No decoys are generated and no variants are expanded.** A database is read as written. A
//!   decoy already in the file (an accession starting `DECOY`) is kept and flagged `is_decoy`. A
//!   UniProt XML that records a **genotype** still has those variants applied by mzLib, which
//!   renames the accession (`P38936_C117Y`); the file's caveats say how many.
//! - **Masses are of the unmodified sequence as written**, plus one water — the precursor,
//!   initiator methionine and signal peptide included — so they are not the mass of the mature
//!   protein.
//!
//! Wire verbs: `proteins read`, `genes resolve` and `proteins classify-peptides`.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::Deserialize;
use serde_json::Value;

use crate::bridge::{self, MzLibError, OnError, Result};
use crate::readers::Table;

/// Every table [`read_with`] can return, in the wire's order.
pub const TABLES: [&str; 3] = ["proteins", "go_terms", "ensembl_genes"];

/// Every protein containing the peptide has the same sequence. Identical sequences under two
/// accessions count as one, because no peptide can tell them apart.
pub const UNIQUE: &str = "Unique";
/// The peptide is in more than one distinct sequence, and all of them share a gene: it supports
/// the gene, not any one isoform.
pub const SHARED_WITHIN_GENE: &str = "SharedWithinGene";
/// The peptide is in sequences with no gene in common.
pub const SHARED_ACROSS_GENES: &str = "SharedAcrossGenes";
/// No target protein in the databases contains the peptide.
pub const NOT_IN_DATABASE: &str = "NotInDatabase";

/// Every outcome [`resolve_genes_with`] can report, as written on the wire (mzLib's
/// `GeneResolutionTsv.OutcomeName`). Every protein gets exactly one.
pub const OUTCOMES: [&str; 6] = [
    "resolved",
    "multi_gene",
    "off_primary_only",
    "not_in_source",
    "unrecognized_accession",
    "contaminant_not_mapped",
];

/// The line that separates the database paths from a second list when both travel on stdin.
const SECTION: &str = "--";

// ---------------------------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------------------------

/// A table [`read_with`] can return.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ProteinTable {
    /// One row per protein: the result's own [`ProteinDatabase::columns`].
    Proteins,
    /// One row per (protein, GO term): [`ProteinDatabase::go_terms`].
    GoTerms,
    /// One row per (protein, Ensembl transcript): [`ProteinDatabase::ensembl_genes`].
    EnsemblGenes,
}

impl ProteinTable {
    /// All three, in the wire's order.
    pub const ALL: [ProteinTable; 3] = [Self::Proteins, Self::GoTerms, Self::EnsemblGenes];

    /// The wire name: `"proteins"`, `"go_terms"` or `"ensembl_genes"`.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Proteins => "proteins",
            Self::GoTerms => "go_terms",
            Self::EnsemblGenes => "ensembl_genes",
        }
    }
}

/// How [`read_with`] reads its databases, and what it returns.
#[derive(Debug, Clone)]
pub struct ProteinReadOptions {
    /// Databases to load as contaminants (mzLib's `isContaminant`): `is_contaminant` is then true
    /// on their rows. They come after the main databases in `source_index` order.
    pub contaminants: Vec<PathBuf>,
    /// Which tables to return. Default: [`ProteinTable::Proteins`] alone. Ask for GO terms and
    /// Ensembl genes when you want them: a whole human proteome has about twenty GO rows per
    /// protein. A table not asked for comes back `None` (or, for the proteins table, empty).
    pub tables: Vec<ProteinTable>,
    /// Keep only proteins whose accession **exactly equals** one of these, in every table. Misses
    /// are listed in [`ProteinDatabase::accessions_not_found`]; `P04406-1` does not find `P04406`.
    /// Sent on stdin, so tens of thousands are fine. `None` reads every protein.
    pub accessions: Option<Vec<String>>,
    /// Add a `sequence` column to the proteins table.
    pub sequences: bool,
    /// Databases loaded at once: 1 or more, or -1 for every core. Default 1. A resource choice
    /// only: the output is byte-identical at any value (mzLib's own per-load thread count is
    /// pinned to 1, so this is the only degree of parallelism).
    pub threads: i32,
    /// [`OnError::Fail`] (default) returns the failure of the lowest failing input;
    /// [`OnError::Skip`] records it in that file's [`DatabaseFile::error`] and reads the rest.
    pub on_error: OnError,
    /// Time to allow. `None` (default) waits: a proteome XML takes a while.
    pub timeout: Option<Duration>,
}

impl Default for ProteinReadOptions {
    fn default() -> Self {
        Self {
            contaminants: Vec::new(),
            tables: vec![ProteinTable::Proteins],
            accessions: None,
            sequences: false,
            threads: 1,
            on_error: OnError::Fail,
            timeout: None,
        }
    }
}

/// How [`resolve_genes_with`] resolves: the gene set to count against, and the databases' roles.
///
/// Give **exactly one** of `gtf` and `gene_set`; there is no default gene set, because the
/// release it pins is part of the answer.
#[derive(Debug, Clone)]
pub struct GeneResolveOptions {
    /// An Ensembl GTF, plain or `.gz`; only its `gene` rows are read. Use the **primary-assembly**
    /// GTF (`Species.Assembly.Release.gtf.gz`, not `chr_patch_hapl_scaff`) and keep Ensembl's file
    /// name, which carries the release.
    pub gtf: Option<PathBuf>,
    /// Instead of `gtf`: a compact gene table written by mzLib's `EnsemblGeneSetWriter` (about
    /// 0.5 MB for human, against a 141 MB GTF). It carries the GTF's provenance, so rows are keyed
    /// exactly as against the GTF.
    pub gene_set: Option<PathBuf>,
    /// Optional: Ensembl's `Species.Assembly.Release.uniprot.tsv(.gz)`, for a second opinion. Each
    /// gene row then says whether Ensembl agrees, and a gene only Ensembl links gets its own row
    /// (`source == "ensembl_xref"`). Without it, agreement is unknown, not false.
    pub xref: Option<PathBuf>,
    /// Databases to load as contaminants. Their proteins are `contaminant_not_mapped` — never
    /// mapped, and never silently dropped.
    pub contaminants: Vec<PathBuf>,
    /// Databases loaded (and hashed) at once: 1 or more, or -1 for every core. Default 1; the
    /// output is identical at any value.
    pub threads: i32,
    /// [`OnError::Fail`] (default) or [`OnError::Skip`], as for [`read_with`].
    pub on_error: OnError,
    /// Time to allow. `None` (default) waits: a gzipped GTF takes tens of seconds.
    pub timeout: Option<Duration>,
}

impl Default for GeneResolveOptions {
    fn default() -> Self {
        Self {
            gtf: None,
            gene_set: None,
            xref: None,
            contaminants: Vec::new(),
            threads: 1,
            on_error: OnError::Fail,
            timeout: None,
        }
    }
}

/// How [`classify_peptides_with`] builds its search space.
///
/// There is no `on_error`: every database is part of one search space, so the wire accepts only
/// "fail" — dropping a database that failed to load would report the peptides it contains as
/// unique.
#[derive(Debug, Clone)]
pub struct ClassifyOptions {
    /// Contaminant databases, searched alongside: they are real sequences in the search space.
    pub contaminants: Vec<PathBuf>,
    /// Databases loaded at once: 1 or more, or -1 for every core. Default 1. The classification
    /// itself is one pass; the output is identical at any value.
    pub threads: i32,
    /// Time to allow. `None` (default) waits.
    pub timeout: Option<Duration>,
}

impl Default for ClassifyOptions {
    fn default() -> Self {
        Self {
            contaminants: Vec::new(),
            threads: 1,
            timeout: None,
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------------------------

/// Why one database could not be read, under [`OnError::Skip`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct FileError {
    /// `"usage"` (the file is missing, or not a database this reads) or `"correctness"` (mzLib
    /// could not parse it).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub kind: String,
    /// The error type as the bridge classified it: `"usage"`, or the .NET exception name.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub r#type: String,
    /// What went wrong, naming the file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub message: String,
}

/// What happened to one input database: one per input, in input order, whether it was read or not.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct DatabaseFile {
    /// Position in the input: the databases first, then the contaminants. Every table's
    /// `source_index` column points back here.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub source_index: u64,
    /// The absolute path.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `"UniProtXml"` or `"Fasta"`, chosen from the extension (`.xml`, or
    /// `.fasta`/`.fa`/`.faa`/`.fas`, each optionally `.gz`). `None`: the file failed before its type
    /// was known.
    #[serde(default)]
    pub file_type: Option<String>,
    /// The mzLib loader used, e.g. `"ProteinDbLoader.LoadProteinXML"`. `None`: the file failed first.
    #[serde(default)]
    pub reader: Option<String>,
    /// Whether this database was loaded as contaminants.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub contaminant: bool,
    /// Proteins mzLib loaded from it, before any accession filter.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub protein_count: u64,
    /// Of those, how many proteins the file itself marks as decoys (accession starts `DECOY`).
    /// None are generated.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub decoy_count: u64,
    /// Rows this database contributed to the main table: for [`read_with`], proteins that passed
    /// the accession filter; for [`resolve_genes_with`], resolution rows; for
    /// [`classify_peptides_with`], target proteins searched.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// For [`resolve_genes_with`] only: the lower-case hex sha256 of the database's
    /// **decompressed** bytes, so `.xml` and `.xml.gz` agree. `None` from the other two.
    #[serde(default)]
    pub search_database_sha256: Option<String>,
    /// What this file cannot tell you, or what mzLib did to it: a FASTA has no GO or Ensembl
    /// source; genotype variants were applied; FASTA lines the loader skipped.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Tables and columns this **format** has no source for, so an empty value reads as "this file
    /// cannot say" and never as "there is none". For a FASTA:
    /// `["go_terms", "ensembl_genes", "ensembl_gene_ids"]`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absent_fields: Vec<String>,
    /// Why the file could not be read ([`OnError::Skip`] only); `None` when it was.
    #[serde(default)]
    pub error: Option<FileError>,
}

impl DatabaseFile {
    /// Whether this database was read (its [`Self::error`] is `None`).
    #[must_use]
    pub fn ok(&self) -> bool {
        self.error.is_none()
    }
}

/// A long table that travels beside a result's main one, e.g. [`ProteinDatabase::go_terms`].
#[derive(Debug, Clone, Deserialize)]
pub struct LongTable {
    /// Rows in the table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table itself; its column order is [`Table::names`]. Cells that are lists
    /// (`evidence_codes`, `projects`) read with [`string_lists`].
    #[serde(flatten)]
    pub columns: Table,
}

/// What [`read_with`] returns: the proteins, and the long tables asked for.
///
/// [`Self::columns`] — the proteins table — has one row per protein, in database order:
/// `source_index`, `source_path`, `accession` (as the database wrote it), `name` (entry name,
/// `G3P_HUMAN`), `full_name`, `organism`, `ncbi_taxonomy_id` (a string: it is an identifier),
/// `primary_gene_name`, `gene_names` (a list), `length` (residues), `monoisotopic_mass` (Da;
/// `None` when the sequence holds a letter with no defined mass: X, B, Z, J), `is_contaminant`,
/// `is_decoy`, `is_entrapment`, `ensembl_gene_ids` (a list; empty for a FASTA), and `sequence` with
/// [`ProteinReadOptions::sequences`].
#[derive(Debug, Clone, Deserialize)]
pub struct ProteinDatabase {
    /// Database files given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Database files read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Database files that failed ([`OnError::Skip`] only; otherwise the call returns the error).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// Proteins that passed the accession filter: the rows of the proteins table when it was
    /// requested.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// The tables returned, in the order `proteins`, `go_terms`, `ensembl_genes`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub tables: Vec<String>,
    /// Distinct accessions in the filter; `None` when no accession filter was given.
    #[serde(default)]
    pub accession_filter_count: Option<u64>,
    /// Filter accessions no database contained, in the order given; `None` when no accession
    /// filter was given. **Check this**: matching is exact, so `P04406-1` does not find `P04406`.
    #[serde(default)]
    pub accessions_not_found: Option<Vec<String>>,
    /// Traps in this result as a whole, e.g. proteins with no computable mass.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// One [`DatabaseFile`] per input, in input order — present for one database too.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<DatabaseFile>,
    /// The proteins table, one row per protein; see the type's docs for its columns. **Empty when
    /// the proteins table was not requested** (the wire sends null), which [`Self::tables`] says.
    #[serde(flatten)]
    pub columns: Table,
    /// One row per (protein, GO id): `source_index`, `source_path`, `accession`, `go_id`
    /// (`GO:0005737`), `aspect` (`BiologicalProcess`, `CellularComponent`, `MolecularFunction`, or
    /// `Unknown` when UniProt gave no `C:`/`F:`/`P:` prefix — never guessed), `term_name`,
    /// `evidence_codes` and `projects` (lists). A GO id UniProt repeats with several lines of
    /// evidence is **one** row with the evidence unioned. `None` when go_terms was not requested.
    #[serde(default)]
    pub go_terms: Option<LongTable>,
    /// One row per (protein, Ensembl transcript): `source_index`, `source_path`, `accession`,
    /// `transcript_id` (versioned, as UniProt wrote it), `protein_id`, `gene_id` (stable — **join
    /// on this**), `versioned_gene_id` and `gene_version` (`None` when unversioned, which is not
    /// version 0). `None` when ensembl_genes was not requested.
    #[serde(default)]
    pub ensembl_genes: Option<LongTable>,
}

impl ProteinDatabase {
    /// Accession → NCBI taxonomy id, for every protein in the table.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if the proteins table was not requested.
    pub fn taxonomy(&self) -> Result<BTreeMap<String, Option<String>>> {
        self.by_accession("ncbi_taxonomy_id")
    }

    /// Accession → organism name, for every protein in the table.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if the proteins table was not requested.
    pub fn organisms(&self) -> Result<BTreeMap<String, Option<String>>> {
        self.by_accession("organism")
    }

    fn by_accession(&self, column: &str) -> Result<BTreeMap<String, Option<String>>> {
        if !self.tables.iter().any(|t| t == "proteins") {
            return Err(MzLibError::Usage(
                "The proteins table was not requested; include ProteinTable::Proteins in \
                 ProteinReadOptions::tables."
                    .to_owned(),
            ));
        }
        let accessions = self.columns.strings("accession")?;
        let values = self.columns.strings(column)?;
        Ok(accessions
            .into_iter()
            .zip(values)
            .filter_map(|(accession, value)| accession.map(|a| (a, value)))
            .collect())
    }
}

/// The Ensembl gene set every resolution was counted against, and where it came from.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct GeneSet {
    /// The GTF's file name. With [`GeneResolveOptions::gene_set`], still the **GTF's** name: a
    /// compact table carries the provenance of the GTF it was made from.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub source_file_name: String,
    /// Lower-case hex sha256 of the GTF's bytes as read — the compressed bytes for a `.gz`, which is
    /// what Ensembl publishes and checksums.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sha256: String,
    /// The Ensembl release parsed from the file name (`"116"` from
    /// `Homo_sapiens.GRCh38.116.gtf.gz`); `None` when the name carries none. Keep Ensembl's name.
    #[serde(default)]
    pub release: Option<String>,
    /// The `#!genome-build` header, e.g. `"GRCh38.p14"`; `None` if absent.
    #[serde(default)]
    pub genome_build: Option<String>,
    /// The `#!genebuild-last-updated` header; `None` if absent.
    #[serde(default)]
    pub genebuild_last_updated: Option<String>,
    /// Genes in the set.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub gene_count: u64,
}

/// Ensembl's own accession-to-gene cross-references, when [`GeneResolveOptions::xref`] was given.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct XrefTable {
    /// The file's name, e.g. `Homo_sapiens.GRCh38.116.uniprot.tsv.gz`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub source_file_name: String,
    /// Lower-case hex sha256 of the file's bytes as read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sha256: String,
    /// The Ensembl release from the file name; `None` if absent.
    #[serde(default)]
    pub release: Option<String>,
    /// Distinct accessions the table links.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub accession_count: u64,
}

/// What [`resolve_genes_with`] returns: one row per (protein, gene), or one outcome row per protein
/// with no gene. Never a `|`-joined cell, and never a pick among several genes.
///
/// The columns after `source_index` and `source_path` are **mzLib's own** `GeneResolutionTsv`
/// schema, so this table and one written by mzLib's TSV writer mean the same: `accession`,
/// `entry_accession`, `isoform`, `namespace`, `outcome` (one of [`OUTCOMES`]), `n_genes`,
/// `gene_id` (stable, `ENSG…`), `versioned_gene_id`, `gene_symbol` (the **release's** symbol —
/// display only), `gene_biotype`, `off_primary_genes` (linked genes outside the set: ALT
/// haplotypes, patches), `uniprot_gene_name`, `source` (`search_database_dbreference`, or
/// `ensembl_xref` for a gene only Ensembl links), `search_database_sha256`, `gene_set_release`,
/// `gene_set_sha256`, `ensembl_xref_agrees` (`None` with no gene on the row or no xref: **unknown,
/// not false**), `ensembl_xref_info_type` and `ensembl_xref_sha256`.
#[derive(Debug, Clone, Deserialize)]
pub struct GeneResolutions {
    /// Database files given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Database files read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Database files that failed ([`OnError::Skip`] only).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// Proteins resolved.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub protein_count: u64,
    /// Rows in the table: at least one per protein, one per gene for `multi_gene`, plus any
    /// `ensembl_xref` rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// The gene set and its provenance.
    pub gene_set: GeneSet,
    /// The xref table and its provenance; `None` when no xref was given.
    #[serde(default)]
    pub xref: Option<XrefTable>,
    /// Outcome → number of **proteins** with it. Every key in [`OUTCOMES`] is present, zeros
    /// included.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub outcome_counts: BTreeMap<String, u64>,
    /// Traps in this result: FASTA proteins (all `not_in_source`), a GTF with no release or genome
    /// build, no xref (agreement unknown, not false).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// One [`DatabaseFile`] per input, each with its `search_database_sha256`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<DatabaseFile>,
    /// The table; see the type's docs for its columns.
    #[serde(flatten)]
    pub columns: Table,
}

/// What [`classify_peptides_with`] returns: one row per peptide, in the order given.
///
/// Columns: `peptide` (exactly as given — **not** I/L-folded), `sharing` ([`UNIQUE`],
/// [`SHARED_WITHIN_GENE`], [`SHARED_ACROSS_GENES`] or [`NOT_IN_DATABASE`]), `accession_count`
/// (target proteins containing it), `accessions` (a list: those proteins, distinct, ordinal order)
/// and `shared_gene_keys` (a list: the gene keys common to all of them, e.g.
/// `ensembl:ENSG00000111640`, `gene:Homo sapiens:GAPDH`, `entry:P04406`; empty for
/// `NotInDatabase` and `SharedAcrossGenes`). Read the list columns with [`string_lists`].
#[derive(Debug, Clone, Deserialize)]
pub struct PeptideClassification {
    /// Database files searched.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Rows: peptides, one per input peptide, duplicates included.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub peptide_count: u64,
    /// Target proteins searched. Contaminants count: they are real sequences in the search space.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub target_protein_count: u64,
    /// Decoy proteins in the databases, which are never searched.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub decoy_proteins_ignored: u64,
    /// Always `true`: I and L are one residue for matching. On the wire so that the rule travels
    /// with every result rather than living only in documentation.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub i_and_l_equivalent: bool,
    /// Class → peptides in it; every class present, zeros included.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sharing_counts: BTreeMap<String, u64>,
    /// One [`DatabaseFile`] per input; its `record_count` is the target proteins searched from it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<DatabaseFile>,
    /// The table; see the type's docs for its columns.
    #[serde(flatten)]
    pub columns: Table,
}

impl PeptideClassification {
    /// Peptide → sharing class. A peptide given twice appears once: it classifies the same.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if the payload has no `peptide` or `sharing` column of strings.
    pub fn sharing_of(&self) -> Result<BTreeMap<String, String>> {
        let peptides = self.columns.strings("peptide")?;
        let sharing = self.columns.strings("sharing")?;
        Ok(peptides
            .into_iter()
            .zip(sharing)
            .filter_map(|(peptide, class)| Some((peptide?, class?)))
            .collect())
    }
}

/// A column whose every cell is a list of strings — `gene_names`, `ensembl_gene_ids`,
/// `accessions`, `shared_gene_keys`, `evidence_codes`, `projects` — with a wire `null` as `None`.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the column is absent (the message names the columns there are),
/// [`MzLibError::Protocol`] if a cell is not a list of strings.
pub fn string_lists(table: &Table, column: &str) -> Result<Vec<Option<Vec<String>>>> {
    let cells = table.raw(column).ok_or_else(|| {
        MzLibError::Usage(format!(
            "No column '{column}' in this table. Its columns are: {}.",
            table.names().join(", ")
        ))
    })?;
    cells
        .iter()
        .enumerate()
        .map(|(row, cell)| match cell {
            Value::Null => Ok(None),
            Value::Array(items) => items
                .iter()
                .map(|item| item.as_str().map(str::to_owned))
                .collect::<Option<Vec<_>>>()
                .map(Some)
                .ok_or_else(|| list_error(column, row, cell)),
            other => Err(list_error(column, row, other)),
        })
        .collect()
}

fn list_error(column: &str, row: usize, cell: &Value) -> MzLibError {
    MzLibError::Protocol(format!(
        "Column '{column}' row {row} is not a list of strings: {cell}"
    ))
}

// ---------------------------------------------------------------------------------------------
// Argument assembly
// ---------------------------------------------------------------------------------------------

/// The database options, and the stdin path lines when there are several databases.
///
/// One database travels as `--path` (plus `--contaminant`); several go on stdin, one per line, a
/// contaminant's line ending in a tab and `contaminant`.
fn database_args<P: AsRef<Path>>(
    databases: &[P],
    contaminants: &[PathBuf],
) -> Result<(Vec<String>, Option<Vec<String>>)> {
    let targets = if databases.is_empty() {
        Vec::new()
    } else {
        bridge::path_lines(databases, "database")?
    };
    let extra = if contaminants.is_empty() {
        Vec::new()
    } else {
        bridge::path_lines(contaminants, "contaminant database")?
    };
    if targets.is_empty() && extra.is_empty() {
        return Err(MzLibError::Usage(
            "At least one protein database is required, e.g. 'human.xml'.".to_owned(),
        ));
    }
    let everything: Vec<&String> = targets.iter().chain(&extra).collect();
    for (index, path) in everything.iter().enumerate() {
        if path.as_str() == SECTION {
            return Err(MzLibError::Usage(format!(
                "'{SECTION}' is not a usable database path."
            )));
        }
        if everything[..index].contains(path) {
            return Err(MzLibError::Usage(format!(
                "The database '{path}' is listed twice."
            )));
        }
    }

    if everything.len() == 1 {
        let mut args = vec!["--path".to_owned(), everything[0].clone()];
        if !extra.is_empty() {
            args.push("--contaminant".to_owned());
        }
        return Ok((args, None));
    }

    let mut lines = targets;
    lines.extend(extra.into_iter().map(|p| format!("{p}\tcontaminant")));
    Ok((vec!["--paths-stdin".to_owned()], Some(lines)))
}

/// A list of accessions or peptides for stdin, one per line.
fn entry_lines<S: AsRef<str>>(values: &[S], what: &str) -> Result<Vec<String>> {
    if values.is_empty() {
        return Err(MzLibError::Usage(format!("{what} is empty.")));
    }
    values
        .iter()
        .map(|value| {
            let text = value.as_ref().trim();
            if text.is_empty() {
                return Err(MzLibError::Usage(format!("{what} contains a blank entry.")));
            }
            if text.contains(['\t', '\n', '\r']) {
                return Err(MzLibError::Usage(format!(
                    "An entry in {what} contains a tab or line break: {text:?}."
                )));
            }
            Ok(text.to_owned())
        })
        .collect()
}

/// The path lines and a second list, separated by a `--` line when both are present.
fn stdin(path_lines: Option<Vec<String>>, second: Option<Vec<String>>) -> Option<String> {
    let mut lines = match (path_lines, second) {
        (None, None) => return None,
        (Some(paths), None) => paths,
        (None, Some(second)) => second,
        (Some(mut paths), Some(second)) => {
            paths.push(SECTION.to_owned());
            paths.extend(second);
            paths
        }
    };
    lines.push(String::new());
    Some(lines.join("\n"))
}

fn optional_path(path: &Path, what: &str) -> Result<String> {
    let text = path.to_str().ok_or_else(|| {
        MzLibError::Usage(format!(
            "{what} is not valid UTF-8, which the bridge requires."
        ))
    })?;
    let text = text.trim();
    if text.is_empty() {
        return Err(MzLibError::Usage(format!(
            "{what} must be a non-empty path."
        )));
    }
    Ok(text.to_owned())
}

fn read_request<P: AsRef<Path>>(
    databases: &[P],
    options: &ProteinReadOptions,
) -> Result<(Vec<String>, Option<String>)> {
    if options.tables.is_empty() {
        return Err(MzLibError::Usage(format!(
            "tables must name at least one of {}.",
            TABLES.join(", ")
        )));
    }
    let (db_args, path_lines) = database_args(databases, &options.contaminants)?;
    let mut tables: Vec<&str> = Vec::new();
    for table in ProteinTable::ALL {
        if options.tables.contains(&table) {
            tables.push(table.as_str());
        }
    }
    let mut args = vec!["proteins".to_owned(), "read".to_owned()];
    args.extend(db_args);
    args.extend([
        "--tables".to_owned(),
        tables.join(","),
        "--threads".to_owned(),
        bridge::threads_arg(options.threads)?,
        "--on-error".to_owned(),
        options.on_error.as_str().to_owned(),
    ]);
    if options.sequences {
        args.push("--sequences".to_owned());
    }
    let filter = match &options.accessions {
        None => None,
        Some(accessions) => {
            let lines = entry_lines(accessions, "accessions").map_err(|error| match error {
                MzLibError::Usage(message) if accessions.is_empty() => {
                    MzLibError::Usage(format!("{message} Pass None to read every protein."))
                }
                other => other,
            })?;
            args.push("--accessions-stdin".to_owned());
            Some(lines)
        }
    };
    Ok((args, stdin(path_lines, filter)))
}

fn resolve_request<P: AsRef<Path>>(
    databases: &[P],
    options: &GeneResolveOptions,
) -> Result<(Vec<String>, Option<String>)> {
    let gene_source =
        match (&options.gtf, &options.gene_set) {
            (Some(gtf), None) => ("--gtf", optional_path(gtf, "gtf")?),
            (None, Some(set)) => ("--gene-set", optional_path(set, "gene_set")?),
            _ => return Err(MzLibError::Usage(
                "Give exactly one of gtf (an Ensembl GTF, e.g. 'Homo_sapiens.GRCh38.116.gtf.gz') \
                 or gene_set (a compact table made from one). There is no default gene set: the \
                 release it pins is part of the answer."
                    .to_owned(),
            )),
        };
    let (db_args, path_lines) = database_args(databases, &options.contaminants)?;
    let mut args = vec!["genes".to_owned(), "resolve".to_owned()];
    args.extend(db_args);
    args.extend([
        "--threads".to_owned(),
        bridge::threads_arg(options.threads)?,
        "--on-error".to_owned(),
        options.on_error.as_str().to_owned(),
        gene_source.0.to_owned(),
        gene_source.1,
    ]);
    if let Some(xref) = &options.xref {
        args.push("--xref".to_owned());
        args.push(optional_path(xref, "xref")?);
    }
    Ok((args, stdin(path_lines, None)))
}

fn classify_request<S: AsRef<str>, P: AsRef<Path>>(
    peptides: &[S],
    databases: &[P],
    options: &ClassifyOptions,
) -> Result<(Vec<String>, Option<String>)> {
    let peptide_lines = entry_lines(peptides, "peptides").map_err(|error| match error {
        MzLibError::Usage(_) if peptides.is_empty() => {
            MzLibError::Usage("At least one peptide is required, e.g. [\"PEPTIDEK\"].".to_owned())
        }
        other => other,
    })?;
    let (db_args, path_lines) = database_args(databases, &options.contaminants)?;
    let mut args = vec!["proteins".to_owned(), "classify-peptides".to_owned()];
    args.extend(db_args);
    args.extend([
        "--threads".to_owned(),
        bridge::threads_arg(options.threads)?,
    ]);
    Ok((args, stdin(path_lines, Some(peptide_lines))))
}

fn call<T: serde::de::DeserializeOwned>(
    (args, stdin): (Vec<String>, Option<String>),
    timeout: Option<Duration>,
) -> Result<T> {
    let data = bridge::invoke(&args, stdin.as_deref(), timeout)?;
    serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!(
            "proteins payload could not be interpreted: {error}"
        ))
    })
}

// ---------------------------------------------------------------------------------------------
// The public surface
// ---------------------------------------------------------------------------------------------

/// Read protein databases into one row per protein, with every default: the proteins table
/// alone, no contaminants, no accession filter.
///
/// See [`read_with`] for the reference.
///
/// # Errors
///
/// As [`read_with`].
pub fn read<P: AsRef<Path>>(databases: &[P]) -> Result<ProteinDatabase> {
    read_with(databases, &ProteinReadOptions::default())
}

/// Read protein databases (UniProt XML or FASTA) into one row per protein, with GO terms and
/// Ensembl gene links as long tables on request.
///
/// Loads each database with mzLib's `ProteinDbLoader`, exactly as a search would but with no decoys
/// generated. A UniProt XML carries everything; a FASTA carries organism, taxon (`OX=`) and gene
/// name (`GN=`) only — see [`DatabaseFile::absent_fields`]. One database is sent as a path, several
/// on stdin, in order, databases first and then [`ProteinReadOptions::contaminants`].
///
/// No database, a blank or repeated path, an empty table selection, an empty accession list, or a
/// `threads` of 0 or below -1 is refused before anything is spawned.
#[doc = include_str!("../docs/reference/proteins.read.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::proteins::{read_with, string_lists, ProteinReadOptions, ProteinTable};
///
/// let db = read_with(
///     &["human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"],
///     &ProteinReadOptions {
///         contaminants: vec!["contaminants.fasta".into()],
///         tables: ProteinTable::ALL.to_vec(),
///         ..Default::default()
///     },
/// )?;
/// assert_eq!(db.record_count, 8);
/// assert_eq!(db.go_terms.as_ref().map(|t| t.row_count), Some(47));
/// assert_eq!(db.ensembl_genes.as_ref().map(|t| t.row_count), Some(6));
/// // A FASTA is silent about GO and Ensembl, and says so rather than looking empty.
/// assert_eq!(db.files[1].file_type.as_deref(), Some("Fasta"));
/// assert_eq!(db.files[1].absent_fields, ["go_terms", "ensembl_genes", "ensembl_gene_ids"]);
/// let genes = string_lists(&db.columns, "gene_names")?;
/// assert_eq!(genes[0].as_ref().unwrap()[0], "GAPDH");
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/proteins.read.see-also.md")]
pub fn read_with<P: AsRef<Path>>(
    databases: &[P],
    options: &ProteinReadOptions,
) -> Result<ProteinDatabase> {
    call(read_request(databases, options)?, options.timeout)
}

/// Resolve every protein to stable Ensembl gene ids against an Ensembl GTF, with every other
/// default.
///
/// See [`resolve_genes_with`] for the reference.
///
/// # Errors
///
/// As [`resolve_genes_with`].
pub fn resolve_genes<P: AsRef<Path>>(
    databases: &[P],
    gtf: impl AsRef<Path>,
) -> Result<GeneResolutions> {
    resolve_genes_with(
        databases,
        &GeneResolveOptions {
            gtf: Some(gtf.as_ref().to_path_buf()),
            ..GeneResolveOptions::default()
        },
    )
}

/// Resolve every protein in the given databases to stable Ensembl gene ids, counted against a
/// caller-supplied, release-pinned gene set, with one outcome per protein and the hash of every
/// input.
///
/// Wraps mzLib's `EnsemblGeneResolver`. The gene links come from the database itself (UniProt's
/// `<dbReference type="Ensembl">`), and each is **counted against the gene set you pass**: a gene in
/// the set counts; a gene outside it (an ALT haplotype, a patch) is dropped from `n_genes` but
/// counted in `off_primary_genes`, and a protein with only such genes is `off_primary_only`.
///
/// **You supply the gene set, and you should pin it.** Gene ids, versions, symbols and the very set
/// of genes on the primary assembly change between Ensembl releases, so a resolution means something
/// only relative to one release. Record [`GeneSet::sha256`] with your results. Nothing is downloaded
/// or defaulted here.
///
/// Neither or both of `gtf` and `gene_set`, and anything [`read_with`] would refuse, are refused
/// before anything is spawned.
#[doc = include_str!("../docs/reference/genes.resolve.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::proteins::{resolve_genes_with, GeneResolveOptions};
///
/// let genes = resolve_genes_with(
///     &["human_subset.xml", "human_extra.fasta"],
///     &GeneResolveOptions {
///         gtf: Some("Homo_sapiens.GRCh38.116.gtf".into()),
///         xref: Some("Homo_sapiens.GRCh38.116.uniprot.tsv".into()),
///         contaminants: vec!["contaminants.fasta".into()],
///         ..Default::default()
///     },
/// )?;
/// assert_eq!(genes.outcome_counts["contaminant_not_mapped"], 1);
/// // Every FASTA protein is not_in_source: the format carries no Ensembl links.
/// assert_eq!(genes.outcome_counts["not_in_source"], 4);
/// assert!(genes.files.iter().all(|f| f.search_database_sha256.is_some()));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/genes.resolve.see-also.md")]
pub fn resolve_genes_with<P: AsRef<Path>>(
    databases: &[P],
    options: &GeneResolveOptions,
) -> Result<GeneResolutions> {
    call(resolve_request(databases, options)?, options.timeout)
}

/// Classify peptides by how widely they are shared across the databases, with every other
/// default: no contaminants, one database loaded at a time.
///
/// See [`classify_peptides_with`] for the reference.
///
/// # Errors
///
/// As [`classify_peptides_with`].
pub fn classify_peptides<S: AsRef<str>, P: AsRef<Path>>(
    peptides: &[S],
    databases: &[P],
) -> Result<PeptideClassification> {
    classify_peptides_with(peptides, databases, &ClassifyOptions::default())
}

/// Classify peptides as unique to one sequence, shared within a gene, shared across genes, or not
/// in the databases, treating I and L as the same residue.
///
/// Wraps mzLib's `PeptideUniquenessClassifier.Classify`. A protein contains a peptide when its
/// sequence contains it **anywhere, whatever the protease** — deliberately conservative, so a
/// peptide called `Unique` cannot be explained by another entry at a site the search's cleavage
/// rules happened to skip. **I and L are the same residue** throughout: `LVLNGNPLTLFQER` finds
/// GAPDH's `LVINGNPITIFQER`.
///
/// Two proteins share a gene when they share any key: an Ensembl gene id, organism plus primary gene
/// name (so human and bovine ALB stay apart), or the UniProt entry (so `P12345-2` meets `P12345`).
/// Decoys are ignored; contaminants are included, because they are real sequences in the search
/// space.
///
/// `peptides` are unmodified base sequences, upper case, one result row each, in order, duplicates
/// included; they travel on stdin. Strip modifications first: mzLib refuses `"PEPT[Phospho]IDEK"`
/// and `"peptidek"` rather than guess. No peptide, a blank one, or one with a tab or line break is
/// refused before anything is spawned.
#[doc = include_str!("../docs/reference/proteins.classify-peptides.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::proteins::{classify_peptides_with, string_lists, ClassifyOptions, SHARED_ACROSS_GENES};
///
/// let calls = classify_peptides_with(
///     &["VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"],
///     &["human_subset.xml", "human_extra.fasta"],
///     &ClassifyOptions { contaminants: vec!["contaminants.fasta".into()], ..Default::default() },
/// )?;
/// assert!(calls.i_and_l_equivalent);
/// assert_eq!(calls.sharing_of()?["YLYEIAR"], SHARED_ACROSS_GENES);
/// let accessions = string_lists(&calls.columns, "accessions")?;
/// assert_eq!(accessions[3].as_deref(), Some(&["P02768".to_owned(), "P02769".to_owned()][..]));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/proteins.classify-peptides.see-also.md")]
pub fn classify_peptides_with<S: AsRef<str>, P: AsRef<Path>>(
    peptides: &[S],
    databases: &[P],
    options: &ClassifyOptions,
) -> Result<PeptideClassification> {
    call(
        classify_request(peptides, databases, options)?,
        options.timeout,
    )
}

#[cfg(test)]
mod tests {
    //! Offline. The recordings are pyMzLib's, byte for byte, made against the small databases in
    //! `tests/fixtures/proteins/`.

    use super::*;

    const READ: &str = include_str!("../tests/fixtures/proteins_read_human.json");
    const FILTERED: &str = include_str!("../tests/fixtures/proteins_read_filtered.json");
    const GENES: &str = include_str!("../tests/fixtures/genes_resolve_human.json");
    const CLASSIFY: &str = include_str!("../tests/fixtures/proteins_classify_peptides.json");

    fn contaminants() -> Vec<PathBuf> {
        vec!["contaminants.fasta".into()]
    }

    // ---- read -----------------------------------------------------------------------------------

    #[test]
    fn a_read_parses_the_proteins_and_both_long_tables() {
        let db: ProteinDatabase = serde_json::from_str(READ).unwrap();
        assert_eq!((db.file_count, db.read_count, db.failed_count), (4, 4, 0));
        assert_eq!(db.record_count, 8);
        assert_eq!(db.columns.rows(), 8);
        assert_eq!(db.tables, TABLES);
        assert_eq!(db.accession_filter_count, None);
        assert_eq!(db.accessions_not_found, None);
        assert_eq!(db.go_terms.as_ref().unwrap().row_count, 47);
        assert_eq!(db.go_terms.as_ref().unwrap().columns.rows(), 47);
        assert_eq!(db.ensembl_genes.as_ref().unwrap().row_count, 6);
        assert_eq!(db.taxonomy().unwrap()["Q9Z0X1"].as_deref(), Some("10090"));
        assert_eq!(
            db.organisms().unwrap()["P02769"].as_deref(),
            Some("Bos taurus")
        );
    }

    #[test]
    fn a_fasta_says_it_is_silent_rather_than_looking_empty() {
        let db: ProteinDatabase = serde_json::from_str(READ).unwrap();
        let fasta = &db.files[1];
        assert_eq!(fasta.file_type.as_deref(), Some("Fasta"));
        assert_eq!(
            fasta.absent_fields,
            ["go_terms", "ensembl_genes", "ensembl_gene_ids"]
        );
        assert!(!fasta.caveats.is_empty());
        assert!(db.files[3].contaminant);
        assert!(db.files.iter().all(DatabaseFile::ok));
        // Its proteins have an empty Ensembl list: absent, per absent_fields, not "no gene".
        let ids = string_lists(&db.columns, "ensembl_gene_ids").unwrap();
        assert_eq!(ids[0].as_deref(), Some(&["ENSG00000111640".to_owned()][..]));
    }

    #[test]
    fn an_unversioned_gene_is_none_not_version_zero() {
        let db: ProteinDatabase = serde_json::from_str(READ).unwrap();
        let versions = db
            .ensembl_genes
            .unwrap()
            .columns
            .integers("gene_version")
            .unwrap();
        assert_eq!(versions.last(), Some(&None));
        assert_eq!(versions[0], Some(15));
    }

    #[test]
    fn a_filtered_read_names_the_accessions_it_did_not_find() {
        let db: ProteinDatabase = serde_json::from_str(FILTERED).unwrap();
        assert_eq!(db.accession_filter_count, Some(2));
        assert_eq!(
            db.accessions_not_found.as_deref(),
            Some(&["P04406-1".to_owned()][..])
        );
        assert!(db.ensembl_genes.is_none(), "a table not asked for is None");
        assert_eq!(
            db.files.len(),
            1,
            "files[] is present with one database too"
        );
    }

    #[test]
    fn the_proteins_table_not_requested_is_refused_by_the_helpers() {
        let mut payload: Value = serde_json::from_str(FILTERED).unwrap();
        payload["tables"] = serde_json::json!(["go_terms"]);
        payload["column_names"] = Value::Null;
        payload["columns"] = Value::Null;
        let db: ProteinDatabase = serde_json::from_value(payload).unwrap();
        assert!(db.columns.is_empty());
        assert!(matches!(db.taxonomy(), Err(MzLibError::Usage(_))));
    }

    #[test]
    fn several_databases_travel_on_stdin_with_their_roles_and_the_filter_after_a_separator() {
        let (args, stdin) = read_request(
            &["a.xml", "b.fasta"],
            &ProteinReadOptions {
                contaminants: contaminants(),
                tables: vec![ProteinTable::EnsemblGenes, ProteinTable::Proteins],
                accessions: Some(vec!["P04406".to_owned(), " Q13409 ".to_owned()]),
                sequences: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "proteins",
                "read",
                "--paths-stdin",
                "--tables",
                "proteins,ensembl_genes",
                "--threads",
                "1",
                "--on-error",
                "fail",
                "--sequences",
                "--accessions-stdin"
            ]
        );
        assert_eq!(
            stdin.unwrap(),
            "a.xml\nb.fasta\ncontaminants.fasta\tcontaminant\n--\nP04406\nQ13409\n"
        );
    }

    #[test]
    fn one_database_is_a_path_and_a_filter_alone_fills_stdin() {
        let (args, stdin) = read_request(
            &["a.xml"],
            &ProteinReadOptions {
                accessions: Some(vec!["P04406".to_owned()]),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(&args[2..4], ["--path", "a.xml"]);
        assert!(!args.contains(&"--contaminant".to_owned()));
        assert_eq!(stdin.as_deref(), Some("P04406\n"));

        let (args, stdin) = read_request(
            &[] as &[&str],
            &ProteinReadOptions {
                contaminants: contaminants(),
                on_error: OnError::Skip,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            &args[2..5],
            ["--path", "contaminants.fasta", "--contaminant"]
        );
        assert!(args.ends_with(&["--on-error".to_owned(), "skip".to_owned()]));
        assert_eq!(stdin, None);
    }

    #[test]
    fn a_read_is_refused_before_anything_is_spawned() {
        let none: &[&str] = &[];
        for (databases, options) in [
            (none, ProteinReadOptions::default()),
            (&["a.xml", "a.xml"][..], ProteinReadOptions::default()),
            (
                &["a.xml"][..],
                ProteinReadOptions {
                    tables: vec![],
                    ..Default::default()
                },
            ),
            (
                &["a.xml"][..],
                ProteinReadOptions {
                    accessions: Some(vec![]),
                    ..Default::default()
                },
            ),
            (
                &["a.xml"][..],
                ProteinReadOptions {
                    accessions: Some(vec!["  ".to_owned()]),
                    ..Default::default()
                },
            ),
            (
                &["a.xml"][..],
                ProteinReadOptions {
                    threads: 0,
                    ..Default::default()
                },
            ),
            (&["--"][..], ProteinReadOptions::default()),
            (&["a\tb.xml"][..], ProteinReadOptions::default()),
        ] {
            let error = read_request(databases, &options).unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{databases:?}");
        }
    }

    // ---- genes resolve --------------------------------------------------------------------------

    #[test]
    fn a_resolution_carries_its_provenance_and_one_outcome_per_protein() {
        let genes: GeneResolutions = serde_json::from_str(GENES).unwrap();
        assert_eq!((genes.protein_count, genes.record_count), (7, 7));
        assert_eq!(genes.gene_set.release.as_deref(), Some("116"));
        assert_eq!(genes.gene_set.genome_build.as_deref(), Some("GRCh38.p14"));
        assert_eq!(genes.gene_set.gene_count, 3);
        assert_eq!(genes.xref.as_ref().unwrap().accession_count, 2);
        for outcome in OUTCOMES {
            assert!(genes.outcome_counts.contains_key(outcome), "{outcome}");
        }
        assert_eq!(genes.outcome_counts.values().sum::<u64>(), 7);
        assert_eq!(
            genes.files[0].search_database_sha256.as_deref(),
            Some("6717e4ebfe0b02023a54da34a605a1d23e855bb9c27cac1f8fed94160759ffe8")
        );
        let outcomes = genes.columns.strings("outcome").unwrap();
        assert_eq!(outcomes[0].as_deref(), Some("resolved"));
        assert_eq!(outcomes[6].as_deref(), Some("contaminant_not_mapped"));
    }

    #[test]
    fn no_xref_is_none_not_an_empty_table() {
        let mut payload: Value = serde_json::from_str(GENES).unwrap();
        payload["xref"] = Value::Null;
        let genes: GeneResolutions = serde_json::from_value(payload).unwrap();
        assert!(genes.xref.is_none());
    }

    #[test]
    fn exactly_one_gene_source_is_required() {
        let both = GeneResolveOptions {
            gtf: Some("a.gtf".into()),
            gene_set: Some("a.tsv".into()),
            ..Default::default()
        };
        for options in [GeneResolveOptions::default(), both] {
            let error = resolve_request(&["a.xml"], &options).unwrap_err();
            assert!(error.to_string().contains("exactly one of gtf"), "{error}");
        }
    }

    #[test]
    fn a_resolution_names_its_gene_source_and_xref() {
        let (args, stdin) = resolve_request(
            &["a.xml", "b.fasta"],
            &GeneResolveOptions {
                gene_set: Some("genes.tsv".into()),
                xref: Some("x.tsv".into()),
                contaminants: contaminants(),
                threads: -1,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "genes",
                "resolve",
                "--paths-stdin",
                "--threads",
                "-1",
                "--on-error",
                "fail",
                "--gene-set",
                "genes.tsv",
                "--xref",
                "x.tsv"
            ]
        );
        assert_eq!(
            stdin.unwrap(),
            "a.xml\nb.fasta\ncontaminants.fasta\tcontaminant\n"
        );
    }

    // ---- classify-peptides ----------------------------------------------------------------------

    #[test]
    fn a_classification_parses_with_i_and_l_on_the_wire() {
        let calls: PeptideClassification = serde_json::from_str(CLASSIFY).unwrap();
        assert!(calls.i_and_l_equivalent);
        assert_eq!(calls.peptide_count, 6);
        assert_eq!(calls.target_protein_count, 7);
        assert_eq!(calls.sharing_counts.values().sum::<u64>(), 6);
        let sharing = calls.sharing_of().unwrap();
        assert_eq!(sharing["LVLNGNPLTLFQER"], UNIQUE); // GAPDH's LVINGNPITIFQER, I = L
        assert_eq!(sharing["ALSEQINIFFDYSGR"], SHARED_WITHIN_GENE);
        assert_eq!(sharing["PEPTIDEK"], NOT_IN_DATABASE);
        let keys = string_lists(&calls.columns, "shared_gene_keys").unwrap();
        assert_eq!(
            keys[3].as_deref(),
            Some(&[][..]),
            "SharedAcrossGenes has no common key"
        );
    }

    #[test]
    fn peptides_follow_the_paths_after_a_separator_or_fill_stdin_alone() {
        let (args, stdin) = classify_request(
            &["PEPTIDEK", "AEFVEVTK"],
            &["a.xml", "b.fasta"],
            &ClassifyOptions::default(),
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "proteins",
                "classify-peptides",
                "--paths-stdin",
                "--threads",
                "1"
            ]
        );
        assert_eq!(stdin.unwrap(), "a.xml\nb.fasta\n--\nPEPTIDEK\nAEFVEVTK\n");

        let (args, stdin) =
            classify_request(&["PEPTIDEK"], &["a.xml"], &ClassifyOptions::default()).unwrap();
        assert_eq!(&args[2..4], ["--path", "a.xml"]);
        assert_eq!(stdin.as_deref(), Some("PEPTIDEK\n"));
    }

    #[test]
    fn no_peptides_or_a_blank_one_is_refused() {
        let error =
            classify_request(&[] as &[&str], &["a.xml"], &ClassifyOptions::default()).unwrap_err();
        assert!(
            error.to_string().contains("At least one peptide"),
            "{error}"
        );
        let error = classify_request(&["PEPTIDEK", " "], &["a.xml"], &ClassifyOptions::default())
            .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    #[test]
    fn a_list_column_that_is_not_a_list_is_a_protocol_error() {
        let calls: PeptideClassification = serde_json::from_str(CLASSIFY).unwrap();
        assert!(matches!(
            string_lists(&calls.columns, "peptide"),
            Err(MzLibError::Protocol(_))
        ));
        assert!(matches!(
            string_lists(&calls.columns, "nonesuch"),
            Err(MzLibError::Usage(_))
        ));
    }
}
