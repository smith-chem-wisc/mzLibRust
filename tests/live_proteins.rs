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
    self, string_lists, AlignOptions, Alignments, ClassifyOptions, EndGaps, GeneResolutions,
    GeneResolveOptions, OrthologyOptions, OrthologySnapshot, PeptideClassification,
    ProteinDatabase, ProteinReadOptions, ProteinTable,
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

// ---- pyMzLib 0.6.0: align, parse-accession, genes orthology ---------------------------------
//
// They need the bridge from pyMzLib 0.6.0, and skip on an older one.

fn sequences() -> std::collections::HashMap<String, String> {
    let read = proteins::read_with(
        &[db("human_extra.fasta")],
        &ProteinReadOptions {
            contaminants: vec![db("contaminants.fasta")],
            sequences: true,
            ..Default::default()
        },
    )
    .expect("the human FASTA reads with its sequences");
    read.columns
        .strings("accession")
        .unwrap()
        .into_iter()
        .zip(read.columns.strings("sequence").unwrap())
        .filter_map(|(accession, sequence)| Some((accession?, sequence?)))
        .collect()
}

#[test]
fn align_leaves_every_option_to_mzlib_and_says_what_it_used() {
    let Some(()) = require_verb("proteins align") else {
        return;
    };
    let aln =
        proteins::align_with(&[("MLLLPEPTIDE", "PEPTIDE")], &AlignOptions::default()).unwrap();
    assert_eq!(
        (
            aln.matrix.as_str(),
            aln.gap_open,
            aln.gap_extend,
            aln.end_gaps.as_str(),
            aln.max_cells
        ),
        ("BLOSUM62", 11, 1, "free", 200_000_000)
    );
    assert_eq!(
        aln.aligner_id,
        "global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"
    );
    assert_eq!(aln.columns.integers("score").unwrap(), [Some(39)]);
}

#[test]
fn align_maps_canonical_dync1i2_onto_its_isoforms_as_recorded() {
    let Some(()) = require_verb("proteins align") else {
        return;
    };
    let seq = sequences();
    let live = proteins::align_with(
        &[
            (&seq["Q13409"], &seq["Q13409-2"]),
            (&seq["Q13409"], &seq["Q13409-3"]),
        ],
        &AlignOptions::default(),
    )
    .unwrap();
    let want: Alignments = recorded(include_str!("fixtures/proteins_align_dync1i2.json"));
    assert_eq!(live.columns, want.columns);
    assert_eq!(
        live.residue_map.as_ref().unwrap().columns,
        want.residue_map.as_ref().unwrap().columns
    );
    assert_eq!(
        live.alignments().unwrap()[0]
            .map_source_to_target(82)
            .unwrap(),
        Some(76)
    );
}

#[test]
fn align_options_reach_mzlib() {
    let Some(()) = require_verb("proteins align") else {
        return;
    };
    let aln = proteins::align_with(
        &[("MLLLPEPTIDE", "PEPTIDE")],
        &AlignOptions {
            gap_open: Some(5),
            gap_extend: Some(2),
            end_gaps: Some(EndGaps::Charged),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        aln.aligner_id,
        "global-affine;BLOSUM62;open=5;extend=2;end-gaps=charged"
    );
    assert_eq!(aln.alignments().unwrap()[0].score, 39 - (5 + 4 * 2));
    let one = proteins::align("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ").unwrap();
    assert_eq!(one.map_source_to_target(12).unwrap(), Some(9));
}

#[test]
fn mzlibs_alignment_refusals_are_usage_errors() {
    let Some(()) = require_verb("proteins align") else {
        return;
    };
    /// The pairs, the options, and what mzLib's refusal says.
    type Refused<'a> = (&'a [(&'a str, &'a str)], AlignOptions, &'a str);
    let cases: [Refused; 4] = [
        (
            &[("PEP-TIDE", "PEPTIDE")],
            AlignOptions::default(),
            "gap character",
        ),
        (
            &[("PEPTIDE", "PEPTIDE"), ("AAAA", "AAAC")],
            AlignOptions {
                max_cells: Some(1),
                ..Default::default()
            },
            "Pair 1 was refused",
        ),
        (
            &[("PEPTIDE", "PEPTIDE")],
            AlignOptions {
                gap_open: Some(-1),
                ..Default::default()
            },
            "--gap-open was refused by mzLib",
        ),
        (
            &[("PEPTIDE", "PEPTIDE")],
            AlignOptions {
                gap_extend: Some(0),
                ..Default::default()
            },
            "--gap-extend was refused by mzLib",
        ),
    ];
    for (pairs, options, message) in cases {
        let error = proteins::align_with(pairs, &options).unwrap_err();
        assert!(
            matches!(&error, mzlib::MzLibError::Usage(m) if m.contains(message)),
            "{message}: {error:?}"
        );
    }
}

#[test]
fn parse_accessions_never_repairs_what_it_is_given() {
    let Some(()) = require_verb("proteins parse-accession") else {
        return;
    };
    let given = [
        "P04637-2",
        " P04637",
        "P04637 ",
        "NP_000537.3_R72P",
        "P04637_2",
    ];
    let parsed = proteins::parse_accessions(&given).unwrap();
    let text = |name: &str| parsed.columns.strings(name).unwrap();
    assert_eq!(text("input"), given.map(|s| Some(s.to_owned())));
    assert_eq!(
        text("namespace"),
        [
            "uniprot",
            "unrecognized",
            "unrecognized",
            "refseq",
            "unrecognized"
        ]
        .map(|s| Some(s.to_owned()))
    );
    assert_eq!(
        text("applied_variants"),
        [None, None, None, Some("R72P".to_owned()), None]
    );
}

#[test]
fn parse_accessions_and_resolve_genes_share_one_namespace_vocabulary() {
    let Some(()) = require_verb("proteins parse-accession") else {
        return;
    };
    let genes = proteins::resolve_genes_with(
        &[db("human_subset.xml")],
        &GeneResolveOptions {
            gtf: Some(db("Homo_sapiens.GRCh38.116.gtf")),
            ..Default::default()
        },
    )
    .unwrap();
    let accessions: Vec<String> = genes
        .columns
        .strings("accession")
        .unwrap()
        .into_iter()
        .flatten()
        .collect();
    let parsed = proteins::parse_accessions(&accessions).unwrap();
    for column in ["namespace", "entry_accession"] {
        assert_eq!(
            parsed.columns.strings(column).unwrap(),
            genes.columns.strings(column).unwrap(),
            "{column}"
        );
    }
    assert_eq!(
        parsed.columns.integers("isoform").unwrap(),
        genes.columns.integers("isoform").unwrap()
    );
}

fn orthology_input(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("orthology")
        .join(name)
}

const SPECIES: [&str; 3] = ["homo_sapiens", "mus_musculus", "rattus_norvegicus"];

/// mzLib's three-species test world, as committed (byte-exact, -text).
fn orthology() -> OrthologyOptions {
    OrthologyOptions {
        gtfs: vec![
            (
                SPECIES[0].to_owned(),
                orthology_input("Homo_sapiens.GRCh38.116.gtf"),
            ),
            (
                SPECIES[1].to_owned(),
                orthology_input("Mus_musculus.GRCm39.116.gtf"),
            ),
            (
                SPECIES[2].to_owned(),
                orthology_input("Rattus_norvegicus.GRCr8.116.gtf"),
            ),
        ],
        dumps: SPECIES
            .iter()
            .map(|s| {
                orthology_input("homologies")
                    .join(s)
                    .join("Compara.116.protein_default.homologies.tsv.gz")
            })
            .collect(),
        gene_trees: orthology_input("vertebrates.GeneTree_content.default.e116.txt.gz"),
        ..Default::default()
    }
}

#[test]
fn orthology_matches_its_recording_over_the_committed_inputs() {
    let Some(()) = require_verb("genes orthology") else {
        return;
    };
    let live = proteins::build_orthology_with(
        "116",
        &OrthologyOptions {
            pair_status: Some((SPECIES[0].to_owned(), SPECIES[1].to_owned())),
            ..orthology()
        },
    )
    .unwrap();
    let want: OrthologySnapshot =
        recorded(include_str!("fixtures/genes_orthology_pair_status.json"));
    assert_eq!(
        live.snapshot_id, want.snapshot_id,
        "the inputs are byte-exact (-text), so the id is too"
    );
    assert_eq!(
        live.pair_status.as_ref().unwrap().columns,
        want.pair_status.as_ref().unwrap().columns
    );
    assert_eq!(live.pairs, want.pairs);
    assert_eq!(live.caveats, want.caveats);
}

#[test]
fn orthology_writes_parquet_once_and_refuses_to_overwrite() {
    let Some(()) = require_verb("genes orthology") else {
        return;
    };
    let scratch = std::env::temp_dir().join(format!("mzlib-live-orthology-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&scratch);
    let out = scratch.join("compara-116");
    let snap = proteins::build_orthology_with(
        "116",
        &OrthologyOptions {
            restrict_gene_trees: true,
            out: Some(out.clone()),
            ..orthology()
        },
    )
    .unwrap();
    let written = snap.written.as_ref().expect("out was given");
    let manifest: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(written.manifest()).unwrap()).unwrap();
    assert_eq!(manifest["snapshot_id"], snap.snapshot_id.as_str());
    let listed: Vec<&str> = manifest["files"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f["path"].as_str().unwrap())
        .collect();
    assert_eq!(
        listed,
        written
            .files
            .iter()
            .map(|f| f.path.as_str())
            .collect::<Vec<_>>()
    );
    for f in &written.files {
        assert_eq!(
            std::fs::metadata(out.join(&f.path)).unwrap().len(),
            f.bytes,
            "{}",
            f.path
        );
    }

    let again = proteins::build_orthology_with(
        "116",
        &OrthologyOptions {
            out: Some(out.clone()),
            ..orthology()
        },
    )
    .unwrap_err();
    assert!(
        matches!(&again, mzlib::MzLibError::Usage(m) if m.contains("is not empty")),
        "{again:?}"
    );
    let _ = std::fs::remove_dir_all(&scratch);
}

#[test]
fn orthology_without_every_dump_is_refused_by_mzlib() {
    let Some(()) = require_verb("genes orthology") else {
        return;
    };
    let mut options = orthology();
    options.dumps.remove(1);
    let error = proteins::build_orthology_with("116", &options).unwrap_err();
    assert!(
        matches!(&error, mzlib::MzLibError::Usage(m) if m.contains("No homology dump for mus_musculus")),
        "{error:?}"
    );
}

#[test]
fn orthology_pair_status_outside_the_species_is_refused_before_reading() {
    let Some(()) = require_verb("genes orthology") else {
        return;
    };
    let error = proteins::build_orthology_with(
        "116",
        &OrthologyOptions {
            pair_status: Some((SPECIES[0].to_owned(), "danio_rerio".to_owned())),
            ..orthology()
        },
    )
    .unwrap_err();
    assert!(
        matches!(&error, mzlib::MzLibError::Usage(m) if m.contains("not one of the species on stdin")),
        "{error:?}"
    );
}

#[test]
fn orthology_a_malformed_dump_is_mzlibs_correctness_failure() {
    let Some(()) = require_verb("genes orthology") else {
        return;
    };
    let scratch = std::env::temp_dir().join(format!("mzlib-live-bad-dump-{}", std::process::id()));
    std::fs::create_dir_all(&scratch).unwrap();
    let bad = scratch.join("Compara.116.protein_default.homologies.tsv");
    std::fs::write(&bad, "gene_stable_id\tspecies\n").unwrap();
    let mut options = orthology();
    options.dumps[2] = bad;
    match proteins::build_orthology_with("116", &options).unwrap_err() {
        mzlib::MzLibError::Bridge {
            error_type,
            message,
        } => {
            assert_eq!(error_type, "InvalidDataException");
            assert!(message.contains("unexpected columns"), "{message}");
        }
        other => panic!("expected mzLib's InvalidDataException, got {other:?}"),
    }
    let _ = std::fs::remove_dir_all(&scratch);
}

// ---- genes map-residues (pyMzLib 0.7.0) -----------------------------------------------------------
//
// test_proteins_map_residues_live.py, test for test, over the committed cut of Ensembl release 116
// in tests/fixtures/map_residues/ (byte for byte: -text). Skips on a bridge without the verb.

use mzlib::proteins::{map_residues_with, ComparaAlignment, GeneTrees, MapResiduesOptions};

const TP53: &str = "ENSP00000269305";
const TRP53: &str = "ENSMUSP00000104298";

fn map_input(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("map_residues")
        .join(name)
}

fn map_options() -> MapResiduesOptions {
    MapResiduesOptions {
        alignment: map_input("Compara.116.protein_default.aa.fasta.gz"),
        gene_trees: map_input("vertebrates.GeneTree_content.default.e116.txt.gz"),
        gtfs: vec![
            (
                "homo_sapiens".into(),
                map_input("Homo_sapiens.GRCh38.116.gtf"),
            ),
            (
                "mus_musculus".into(),
                map_input("Mus_musculus.GRCm39.116.gtf"),
            ),
            (
                "rattus_norvegicus".into(),
                map_input("Rattus_norvegicus.GRCr8.116.gtf"),
            ),
        ],
        ..Default::default()
    }
}

fn map_recording(name: &str) -> serde_json::Value {
    serde_json::from_str(
        &std::fs::read_to_string(
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("fixtures")
                .join(name),
        )
        .unwrap(),
    )
    .unwrap()
}

fn map_usage(result: mzlib::Result<proteins::ResidueMappings>) -> String {
    match result {
        Err(mzlib::MzLibError::Usage(message)) => message,
        other => panic!("expected a usage error, got {other:?}"),
    }
}

#[test]
fn each_map_residues_recording_is_what_the_bridge_answers_now() {
    let Some(()) = require_verb("genes map-residues") else {
        return;
    };
    for (fixture, with_proteins) in [
        ("genes_map_residues_tp53.json", true),
        ("genes_map_residues_chain.json", false),
    ] {
        let recorded = map_recording(fixture);
        let columns = &recorded["columns"];
        let text = |name: &str, i: usize| columns[name][i].as_str().unwrap().to_owned();
        let queries: Vec<(String, i64, String)> = (0..recorded["row_count"].as_u64().unwrap()
            as usize)
            .map(|i| {
                (
                    text("source_protein_id", i),
                    columns["position"][i].as_i64().unwrap(),
                    text("target_protein_id", i),
                )
            })
            .collect();
        let live = map_residues_with(
            &queries,
            &MapResiduesOptions {
                proteins: with_proteins,
                ..map_options()
            },
        )
        .unwrap();
        for name in recorded["column_names"].as_array().unwrap() {
            let name = name.as_str().unwrap();
            assert_eq!(
                serde_json::Value::Array(live.columns.raw(name).unwrap().to_vec()),
                columns[name],
                "{fixture} {name}"
            );
        }
        assert_eq!(
            serde_json::to_value(&live.outcome_counts).unwrap(),
            recorded["outcome_counts"],
            "{fixture}"
        );
        // Byte-exact inputs (-text), so every sha256 is the recording's.
        let alignment: ComparaAlignment =
            serde_json::from_value(recorded["alignment"].clone()).unwrap();
        assert_eq!(live.alignment, alignment, "{fixture}");
        let trees: GeneTrees = serde_json::from_value(recorded["gene_trees"].clone()).unwrap();
        assert_eq!(live.gene_trees, trees, "{fixture}");
        assert_eq!(
            serde_json::to_value(&live.caveats).unwrap(),
            recorded["caveats"],
            "{fixture}"
        );
        match (&live.proteins, recorded["proteins"].is_null()) {
            (None, true) => {}
            (Some(kept), false) => {
                for name in recorded["proteins"]["column_names"].as_array().unwrap() {
                    let name = name.as_str().unwrap();
                    assert_eq!(
                        serde_json::Value::Array(kept.columns.raw(name).unwrap().to_vec()),
                        recorded["proteins"]["columns"][name],
                        "{fixture} proteins {name}"
                    );
                }
            }
            (live, _) => panic!("{fixture}: proteins {live:?}"),
        }
    }
}

#[test]
fn a_position_past_the_source_fails_the_whole_call_with_mzlibs_message() {
    let Some(()) = require_verb("genes map-residues") else {
        return;
    };
    // mzLib's MapResidue throws for it (the bridge checks no length), so the first query's answer
    // is lost too: no row comes back.
    let message = map_usage(map_residues_with(
        &[(TP53, 15, TRP53), (TP53, 394, TRP53)],
        &map_options(),
    ));
    assert!(
        message.contains("Query 2 (ENSP00000269305, position 394, ENSMUSP00000104298)"),
        "{message}"
    );
    assert!(
        message.contains("ENSP00000269305 has 393 residues."),
        "mzLib's own words: {message}"
    );
    assert!(
        message.contains("no query of this call was answered"),
        "{message}"
    );
    let message = map_usage(map_residues_with(&[(TP53, 0, TRP53)], &map_options()));
    assert!(message.contains("has 393 residues"), "{message}");
}

#[test]
fn a_source_not_kept_is_an_outcome_before_its_position_is_looked_at() {
    let Some(()) = require_verb("genes map-residues") else {
        return;
    };
    // MapResidue checks the source before the position, so an unkept source with any position is
    // SourceNotInAlignment, not a refusal.
    let mapped = map_residues_with(&[("ENSP00000352610", 100_000, TRP53)], &map_options()).unwrap();
    assert_eq!(mapped.rows().unwrap()[0].outcome, "SourceNotInAlignment");
}

#[test]
fn without_gene_sets_mzlibs_default_keeps_every_gene_and_a_caveat_says_so() {
    let Some(()) = require_verb("genes map-residues") else {
        return;
    };
    let mapped = map_residues_with(
        &[(TP53, 15, TRP53)],
        &MapResiduesOptions {
            gtfs: Vec::new(),
            ..map_options()
        },
    )
    .unwrap();
    assert!(mapped.gene_sets.is_empty());
    assert_eq!(mapped.gene_trees.restricted_to_gene_sets, None);
    assert_eq!(
        mapped.rows().unwrap()[0].target_position,
        Some(18),
        "the same answer for the same proteins"
    );
    assert!(mapped
        .caveats
        .iter()
        .any(|c| c.contains("No gene sets were given, so mzLib kept every gene")));
}

#[test]
fn a_ragged_alignment_is_mzlibs_correctness_failure_and_a_missing_one_a_usage_error() {
    let Some(()) = require_verb("genes map-residues") else {
        return;
    };
    let scratch =
        std::env::temp_dir().join(format!("mzlib-live-map-residues-{}", std::process::id()));
    std::fs::create_dir_all(&scratch).unwrap();
    // Plain text: mzLib reads the alignment plain or .gz, and this crate carries no gzip writer.
    let bad = scratch.join("Compara.116.protein_default.aa.fasta");
    std::fs::write(
        &bad,
        ">ENSDARP00000000001\nPEPT\n>ENSDARP00000000002\nPEP\n//\n",
    )
    .unwrap();
    match map_residues_with(
        &[(TP53, 15, TRP53)],
        &MapResiduesOptions {
            alignment: bad,
            ..map_options()
        },
    )
    .unwrap_err()
    {
        mzlib::MzLibError::Bridge {
            error_type,
            message,
        } => {
            assert_eq!(error_type, "InvalidDataException");
            assert!(
                message.contains("is 3 columns wide, the alignment 4"),
                "{message}"
            );
        }
        other => panic!("expected mzLib's InvalidDataException, got {other:?}"),
    }
    let message = map_usage(map_residues_with(
        &[(TP53, 15, TRP53)],
        &MapResiduesOptions {
            alignment: scratch.join("absent.aa.fasta.gz"),
            ..map_options()
        },
    ));
    assert!(
        message.contains("Compara gene-tree alignment not found"),
        "{message}"
    );
    let _ = std::fs::remove_dir_all(&scratch);
}
