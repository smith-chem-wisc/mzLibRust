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

// ---- Gene Ontology (mzLib 1.0.593) ----------------------------------------------------------
//
// Local files only: the committed PXD036557 table, its five UniProt entries and the trimmed GO
// release. update_go downloads from GO's PURL, so it skips on an outage.

fn go_options() -> proteins::GoAnnotateOptions {
    proteins::GoAnnotateOptions {
        go_obo: db("go-pxd036557.obo"),
        ..Default::default()
    }
}

#[test]
fn a_go_annotation_matches_its_recording() {
    let Some(()) = require_verb("proteins annotate-go") else {
        return;
    };
    let live = proteins::annotate_go_with(
        db("PXD036557_AllQuantifiedProteinGroups.tsv"),
        db("pxd036557_proteins.xml"),
        &proteins::GoAnnotateOptions {
            category_map: Some(db("organelle_map.tsv")),
            ..go_options()
        },
    )
    .expect("the recorded annotation should succeed live");
    let want: proteins::GoAnnotations =
        recorded(include_str!("fixtures/proteins_annotate_go_pxd036557.json"));

    assert_eq!(
        (live.group_count, live.row_count, live.go.term_count),
        (5, 563, 412)
    );
    assert_eq!(live.groups_file_sha256, want.groups_file_sha256);
    assert_eq!(live.go.sha256, want.go.sha256);
    assert_eq!(
        live.annotation_database.sha256,
        want.annotation_database.sha256
    );
    assert_eq!(live.columns.names(), want.columns.names());
    for column in ["protein_group", "go_id", "annotation_status"] {
        assert_eq!(
            live.columns.strings(column).unwrap(),
            want.columns.strings(column).unwrap(),
            "{column}"
        );
    }
    assert_eq!(
        live.columns.integers("n_with").unwrap(),
        want.columns.integers("n_with").unwrap()
    );
    assert_eq!(live.header["status_annotated"], "4");
    assert_eq!(live.header["status_contaminant"], "1");
    assert_eq!(live.categories.as_ref().map(|c| c.row_count), Some(30));
}

#[test]
fn out_holds_the_whole_table_while_the_wire_holds_none() {
    let Some(()) = require_verb("proteins annotate-go") else {
        return;
    };
    let scratch = std::env::temp_dir().join(format!("mzlib-live-go-{}", std::process::id()));
    std::fs::create_dir_all(&scratch).unwrap();
    let out = scratch.join("go.tsv");
    let live = proteins::annotate_go_with(
        db("PXD036557_AllQuantifiedProteinGroups.tsv"),
        db("pxd036557_proteins.xml"),
        &proteins::GoAnnotateOptions {
            out: Some(out.clone()),
            limit: Some(0),
            ..go_options()
        },
    )
    .expect("writing the table should succeed live");

    let text = std::fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "#!go_annotation_format 1");
    assert_eq!(lines.iter().filter(|l| !l.starts_with("#!")).count(), 564); // header + 563 rows
    assert_eq!(live.written.as_ref().and_then(|w| w.row_count), Some(563));
    assert_eq!((live.returned_count, live.truncated), (0, true));
    let _ = std::fs::remove_dir_all(&scratch);
}

#[test]
fn a_wrong_extension_or_a_missing_go_obo_is_a_usage_error() {
    let Some(()) = require_verb("proteins annotate-go") else {
        return;
    };
    let scratch = std::env::temp_dir().join(format!("mzlib-live-go-usage-{}", std::process::id()));
    std::fs::create_dir_all(&scratch).unwrap();
    let groups = db("PXD036557_AllQuantifiedProteinGroups.tsv");
    let xml = db("pxd036557_proteins.xml");

    let error = proteins::annotate_go_with(
        &groups,
        &xml,
        &proteins::GoAnnotateOptions {
            out: Some(scratch.join("go.csv")),
            ..go_options()
        },
    )
    .unwrap_err();
    assert!(
        matches!(&error, mzlib::MzLibError::Usage(m) if m.contains(".tsv")),
        "{error}"
    );

    let missing = scratch.join("go.obo");
    let error = proteins::annotate_go(&groups, &xml, &missing).unwrap_err();
    assert!(
        matches!(&error, mzlib::MzLibError::Usage(m) if m.contains("update-go")),
        "{error}"
    );
    assert!(
        !missing.exists(),
        "annotate-go must never download a go.obo"
    );
    let _ = std::fs::remove_dir_all(&scratch);
}

#[test]
fn update_go_fetches_a_release_that_loads() {
    let Some(()) = require_verb("proteins update-go") else {
        return;
    };
    let scratch = std::env::temp_dir().join(format!("mzlib-live-update-go-{}", std::process::id()));
    std::fs::create_dir_all(&scratch).unwrap();
    let target = scratch.join("go.obo");
    let Some(update) = support::external_service("GO's PURL", proteins::update_go(&target)) else {
        let _ = std::fs::remove_dir_all(&scratch);
        return;
    };
    assert!(!update.existed_before && update.changed);
    assert_eq!(update.url, "https://purl.obolibrary.org/obo/go.obo");
    assert!(update.go.term_count > 40_000, "{}", update.go.term_count);
    assert!(target.exists());
    let _ = std::fs::remove_dir_all(&scratch);
}
