//! Live checks of the proteins module against the real bridge and the small databases the
//! recordings were made from (`tests/fixtures/proteins/`, shared byte for byte with pyMzLib).
//!
//! The offline suite in `src/proteins.rs` pins the projection against those recordings. These pin
//! that the recordings still match what the bridge emits. They need a bridge that dispatches the
//! mzLib 1.0.592 verbs (pyMzLib 0.2.0 or later) and **skip** rather than fail without one.
//!
//! Run with `cargo test --features live`.

#![cfg(feature = "live")]

mod support;

use std::path::PathBuf;

use mzlib::proteins::{
    self, string_lists, ClassifyOptions, GeneResolutions, GeneResolveOptions,
    PeptideClassification, ProteinDatabase, ProteinReadOptions, ProteinTable,
};
use support::require_verb;

fn db(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("proteins")
        .join(name)
}

fn recorded<T: serde::de::DeserializeOwned>(json: &str) -> T {
    serde_json::from_str(json).unwrap()
}

#[test]
fn a_read_matches_its_recording() {
    let Some(()) = require_verb("proteins read") else {
        return;
    };
    let live = proteins::read_with(
        &[
            db("human_subset.xml"),
            db("human_extra.fasta"),
            db("mouse_aifm1.fasta"),
        ],
        &ProteinReadOptions {
            contaminants: vec![db("contaminants.fasta")],
            tables: ProteinTable::ALL.to_vec(),
            ..Default::default()
        },
    )
    .expect("the recorded read should succeed live");
    let want: ProteinDatabase = recorded(include_str!("fixtures/proteins_read_human.json"));

    assert_eq!(live.record_count, want.record_count);
    assert_eq!(live.tables, want.tables);
    assert_eq!(live.columns.names(), want.columns.names());
    assert_eq!(
        live.columns.strings("accession").unwrap(),
        want.columns.strings("accession").unwrap()
    );
    assert_eq!(live.taxonomy().unwrap(), want.taxonomy().unwrap());
    assert_eq!(
        string_lists(&live.columns, "gene_names").unwrap(),
        string_lists(&want.columns, "gene_names").unwrap()
    );
    assert_eq!(
        live.go_terms.as_ref().unwrap().row_count,
        want.go_terms.as_ref().unwrap().row_count
    );
    assert_eq!(
        live.ensembl_genes.as_ref().unwrap().row_count,
        want.ensembl_genes.as_ref().unwrap().row_count
    );
    let absent: Vec<_> = live.files.iter().map(|f| f.absent_fields.clone()).collect();
    let recorded_absent: Vec<_> = want.files.iter().map(|f| f.absent_fields.clone()).collect();
    assert_eq!(absent, recorded_absent);
}

#[test]
fn an_accession_filter_is_exact_and_says_what_it_missed() {
    let Some(()) = require_verb("proteins read") else {
        return;
    };
    let live = proteins::read_with(
        &[db("human_subset.xml")],
        &ProteinReadOptions {
            tables: vec![ProteinTable::Proteins, ProteinTable::GoTerms],
            accessions: Some(vec!["P04406".to_owned(), "P04406-1".to_owned()]),
            ..Default::default()
        },
    )
    .unwrap();
    let want: ProteinDatabase = recorded(include_str!("fixtures/proteins_read_filtered.json"));
    assert_eq!(live.accession_filter_count, want.accession_filter_count);
    assert_eq!(live.accessions_not_found, want.accessions_not_found);
    assert_eq!(live.record_count, want.record_count);
    assert!(live.ensembl_genes.is_none());
}

#[test]
fn a_resolution_matches_its_recording_and_its_database_hashes() {
    let Some(()) = require_verb("genes resolve") else {
        return;
    };
    let live = proteins::resolve_genes_with(
        &[db("human_subset.xml"), db("human_extra.fasta")],
        &GeneResolveOptions {
            gtf: Some(db("Homo_sapiens.GRCh38.116.gtf")),
            xref: Some(db("Homo_sapiens.GRCh38.116.uniprot.tsv")),
            contaminants: vec![db("contaminants.fasta")],
            ..Default::default()
        },
    )
    .unwrap();
    let want: GeneResolutions = recorded(include_str!("fixtures/genes_resolve_human.json"));

    assert_eq!(live.outcome_counts, want.outcome_counts);
    assert_eq!(
        live.columns.strings("outcome").unwrap(),
        want.columns.strings("outcome").unwrap()
    );
    // The databases are byte-exact (-text), so their decompressed sha256 is the recorded one.
    let hashes = |r: &GeneResolutions| -> Vec<Option<String>> {
        r.files
            .iter()
            .map(|f| f.search_database_sha256.clone())
            .collect()
    };
    assert_eq!(hashes(&live), hashes(&want));
    assert_eq!(live.gene_set.release, want.gene_set.release);
    assert_eq!(live.gene_set.genome_build, want.gene_set.genome_build);
    assert_eq!(live.gene_set.gene_count, want.gene_set.gene_count);
    // Not compared with the recording: pyMzLib's GTF and xref files carry CRLF line endings, while
    // the recording hashed their LF form, so the recorded sha256 of those two is of other bytes.
    assert_eq!(live.gene_set.sha256.len(), 64);
    assert_eq!(live.xref.as_ref().unwrap().sha256.len(), 64);
}

#[test]
fn a_classification_matches_its_recording() {
    let Some(()) = require_verb("proteins classify-peptides") else {
        return;
    };
    let live = proteins::classify_peptides_with(
        &[
            "VGVNGFGR",
            "LVLNGNPLTLFQER",
            "ALSEQINIFFDYSGR",
            "YLYEIAR",
            "AEFVEVTK",
            "PEPTIDEK",
        ],
        &[db("human_subset.xml"), db("human_extra.fasta")],
        &ClassifyOptions {
            contaminants: vec![db("contaminants.fasta")],
            ..Default::default()
        },
    )
    .unwrap();
    let want: PeptideClassification =
        recorded(include_str!("fixtures/proteins_classify_peptides.json"));
    assert_eq!(live.sharing_of().unwrap(), want.sharing_of().unwrap());
    assert_eq!(
        string_lists(&live.columns, "accessions").unwrap(),
        string_lists(&want.columns, "accessions").unwrap()
    );
    assert_eq!(live.sharing_counts, want.sharing_counts);
}

#[test]
fn a_modified_peptide_is_refused_by_mzlib_as_a_usage_error() {
    let Some(()) = require_verb("proteins classify-peptides") else {
        return;
    };
    let error =
        proteins::classify_peptides(&["PEPT[Phospho]IDEK"], &[db("human_subset.xml")]).unwrap_err();
    assert!(
        matches!(error, mzlib::MzLibError::Usage(_)),
        "mzLib refuses a modified peptide rather than guess: {error:?}"
    );
}
