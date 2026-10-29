//! Protein databases: what each protein *is*, which gene it belongs to, whether a peptide
//! identifies it, and what each protein group does.
//!
//! | You want to … | Call | mzLib type |
//! |---|---|---|
//! | know a protein's organism, taxon, genes and mass, and its GO terms and Ensembl links | [`read_with`] | `ProteinDbLoader` |
//! | map proteins to stable Ensembl gene ids against a release you pin | [`resolve_genes_with`] | `EnsemblGeneResolver` |
//! | know whether a peptide is unique, shared within a gene, or shared across genes | [`classify_peptides_with`] | `PeptideUniquenessClassifier` |
//! | annotate MetaMorpheus protein groups with GO terms, every member kept | [`annotate_go_with`] | `GoGroupAnnotator` |
//! | fetch the current GO release, on purpose | [`update_go`] | `Loaders.UpdateGeneOntology` |
//! | carry a position to an isoform, an ortholog or a newer release | [`align_with`], [`align`] | `PairwiseAligner` |
//! | read an accession back into entry, isoform or version, and variants | [`parse_accessions`] | `VariantApplication.ParseAccession` |
//! | know which genes are orthologous across species, as Compara states it | [`build_orthology_with`] | `OrthologySnapshot` |
//! | know which mouse residue is human p53's Ser15, by Compara's gene-tree alignment | [`map_residues_with`] | `ComparaGeneTreeAlignment` |
//!
//! Questions a search result cannot answer on its own, each answered by mzLib from the protein
//! database you searched.
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
//! **What does each protein group do, every member kept?** [`annotate_go_with`] reads the
//! protein-group table MetaMorpheus wrote, the UniProt XML it searched and a go.obo you pinned, and
//! returns one row per (group, GO term) that **any** member holds — directly or through an ancestor
//! in the ontology — naming the members that carry it. MetaMorpheus picks no leading protein, so
//! neither does this: consensus and direct-only views are filters on the rows. Here the table is a
//! real MetaMorpheus search (PXD036557):
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::proteins::{annotate_go_with, GoAnnotateOptions};
//!
//! let go = annotate_go_with(
//!     "PXD036557_AllQuantifiedProteinGroups.tsv",
//!     "pxd036557_proteins.xml",
//!     &GoAnnotateOptions {
//!         go_obo: "go-pxd036557.obo".into(),
//!         category_map: Some("organelle_map.tsv".into()),
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!((go.group_count, go.row_count), (5, 563));
//! assert_eq!(go.header["status_contaminant"], "1");   // bovine albumin: one row, saying why
//!
//! // Consensus is a filter: of the histone group's terms, the ones both members carry.
//! let groups = go.columns.strings("protein_group")?;
//! let (n_with, n_members) = (go.columns.integers("n_with")?, go.columns.integers("n_members")?);
//! let consensus = (0..groups.len())
//!     .filter(|&i| groups[i].as_deref() == Some("P0C0S5|Q71UI9") && n_with[i] == n_members[i])
//!     .count();
//! assert_eq!(consensus, 57);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Pin the ontology.** Terms and their ancestors change between GO releases. [`annotate_go_with`]
//! reads only a go.obo you name; [`update_go`] is the one function that fetches one. Record
//! [`GoAnnotations::go`] (its sha256) with your results.
//!
//! These examples replay the recordings in `tests/fixtures/`, made over the small databases in
//! `tests/fixtures/proteins/` (and, for orthology, mzLib's own three-species test world in
//! `tests/fixtures/orthology/`) — the same ones pyMzLib's examples replay.
//!
//! [`read_with`], [`resolve_genes_with`] and [`classify_peptides_with`] take **one database or
//! many**. A list is read in one bridge call, in order;
//! `threads` says how many are read at once (default 1), and the answer is identical at any value.
//! Mark contaminant databases with `contaminants` — it changes answers: a contaminant is never
//! mapped to a gene, and a peptide it shares with a target is shared.
//!
//! MetaMorpheus writes the organism *name* in its protein-group table, not the NCBI taxon, and a
//! mixed-species search (a host and a pathogen, a xenograft, a spike-in) needs the taxon to join
//! against anything taxonomic. The database you searched has it; ask for your accessions:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::proteins::{read_with, ProteinReadOptions};
//!
//! let found = read_with(
//!     &["human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"],
//!     &ProteinReadOptions {
//!         contaminants: vec!["contaminants.fasta".into()],
//!         accessions: Some(vec!["P04406".into(), "Q9Z0X1".into(), "P02769".into(), "P04406-1".into()]),
//!         ..Default::default()
//!     },
//! )?;
//! let taxa = found.taxonomy()?;
//! assert_eq!(taxa["P04406"].as_deref(), Some("9606"));    // human GAPDH
//! assert_eq!(taxa["Q9Z0X1"].as_deref(), Some("10090"));   // mouse AIFM1
//! assert_eq!(found.organisms()?["P02769"].as_deref(), Some("Bos taurus"));   // the contaminant
//! // An isoform suffix is a different accession: matching is exact, and a miss is named.
//! assert_eq!(found.accessions_not_found.as_deref(), Some(&["P04406-1".to_owned()][..]));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! ## Map a residue from one sequence onto another
//!
//! **The problem.** A position means something only in the sequence it was counted in. A
//! phosphosite reported on the canonical isoform, a variant called against last year's UniProt
//! release, a cysteine you know in the human protein: before you look for it in an isoform, an
//! ortholog or a newer entry, you need the residue that *corresponds* to it there. Counting the
//! same number along the other sequence is wrong as soon as there is one insertion or deletion
//! before it.
//!
//! [`align_with`] aligns pairs of sequences with mzLib's `PairwiseAligner` (mzLib #1428) and maps
//! every residue to its counterpart in the other sequence, or to `None` when it has none. Here the
//! canonical DYNC1I2 sequence (Q13409) is aligned with its isoforms 2B and 2C, read from the
//! database with [`ProteinReadOptions::sequences`]:
//!
//! ```
//! # mzlib_replay::activate();
//! use std::collections::HashMap;
//!
//! use mzlib::proteins::{align_with, read_with, AlignOptions, ProteinReadOptions};
//!
//! let db = read_with(
//!     &["human_extra.fasta"],
//!     &ProteinReadOptions {
//!         contaminants: vec!["contaminants.fasta".into()],
//!         accessions: Some(
//!             ["Q13409", "Q13409-2", "Q13409-3", "P02768", "P02769"].map(String::from).to_vec(),
//!         ),
//!         sequences: true,
//!         ..Default::default()
//!     },
//! )?;
//! let seq: HashMap<String, String> = db
//!     .columns
//!     .strings("accession")?
//!     .into_iter()
//!     .zip(db.columns.strings("sequence")?)
//!     .filter_map(|(accession, sequence)| Some((accession?, sequence?)))
//!     .collect();
//!
//! let isoforms = align_with(
//!     &[(&seq["Q13409"], &seq["Q13409-2"]), (&seq["Q13409"], &seq["Q13409-3"])],
//!     &AlignOptions::default(),
//! )?;
//! let each = isoforms.alignments()?;
//! let (two_b, two_c) = (&each[0], &each[1]);
//! for p in [75, 76, 82, 100, 120, 133] {
//!     println!("{p} {:?} {:?}", two_b.map_source_to_target(p)?, two_c.map_source_to_target(p)?);
//! }
//! //   75 Some(75) Some(75)
//! //   76 None None
//! //   82 Some(76) Some(76)
//! //   100 Some(94) Some(94)
//! //   120 Some(114) None
//! //   133 Some(127) Some(107)
//! assert_eq!(two_b.map_source_to_target(76)?, None);
//! assert_eq!(two_c.map_source_to_target(133)?, Some(107));
//!
//! // Where the gaps fall: '-' wherever a residue faces one (characters 60 to 139 of 2C's rows).
//! assert_eq!(
//!     &two_c.aligned_source[60..140],
//!     "EALLQSMGLTPESPIVFSEYWVPPPMSPSSKSVSTPSEAGSQDSGDGAVGSRTLHWDTDPSVLQLHSDSDLGRGPIKLGM"
//! );
//! assert_eq!(
//!     &two_c.aligned_target[60..140],
//!     "EALLQSMGLTPESPI------VPPPMSPSSKSVSTPSEAGSQDSGDGAVGSR--------------------RGPIKLGM"
//! );
//! assert_eq!((two_b.identical, two_b.substituted, two_b.source_residues_facing_gap), (632, 0, 6));
//! assert_eq!(isoforms.aligner_id, "global-affine;BLOSUM62;open=11;extend=1;end-gaps=free");
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Reading it row by row:
//!
//! - **Up to residue 75 the numbers agree.** Both isoforms begin with the same 75 residues.
//! - **Residues 76 to 81 have no counterpart in either isoform**: `None`, never the nearest
//!   residue. Both isoforms lack those six residues of the canonical sequence.
//! - **After the gap the numbers are shifted**: canonical residue 100 is residue 94 of each
//!   isoform.
//! - **Isoform 2C lacks 20 more** (canonical 113 to 132), so canonical 120 has no counterpart in
//!   2C, and canonical 133 is 2C's residue 107.
//!
//! [`Alignment::map_target_to_source`] goes the other way. [`Alignment::aligned_source`] and
//! [`Alignment::aligned_target`] hold `-` wherever a residue faces a gap: print the two one above
//! the other to see the alignment, as the last asserts above do for isoform 2C.
//!
//! **What each number means.** `identical` and
//! `substituted` count the paired residues that are the same letter and that differ;
//! `source_residues_facing_gap` and `target_residues_facing_gap` count the residues with no
//! counterpart. `score` is the alignment's BLOSUM62 score in half-bits (see the traps).
//! [`Alignments::aligner_id`] names the method and every parameter, so store it beside any
//! position you map.
//!
//! **Orthologs: the counterpart can be a different amino acid.** Across species most residues
//! still pair, but many pair with a *different* residue. Human serum albumin's free thiol, Cys34
//! of the mature protein, is residue 58 of the UniProt sequence, which counts the 24-residue
//! signal peptide and propeptide. In bovine albumin it is a cysteine at 58 too. Residue 26 pairs
//! as well, but it is A in human and T in bovine albumin:
//!
//! ```
//! # mzlib_replay::activate();
//! # use std::collections::HashMap;
//! # use mzlib::proteins::{align_with, read_with, AlignOptions, ProteinReadOptions};
//! # let db = read_with(
//! #     &["human_extra.fasta"],
//! #     &ProteinReadOptions {
//! #         contaminants: vec!["contaminants.fasta".into()],
//! #         accessions: Some(
//! #             ["Q13409", "Q13409-2", "Q13409-3", "P02768", "P02769"].map(String::from).to_vec(),
//! #         ),
//! #         sequences: true,
//! #         ..Default::default()
//! #     },
//! # )?;
//! # let seq: HashMap<String, String> = db.columns.strings("accession")?.into_iter()
//! #     .zip(db.columns.strings("sequence")?)
//! #     .filter_map(|(a, s)| Some((a?, s?))).collect();
//! let albumin = mzlib::proteins::align(&seq["P02768"], &seq["P02769"])?;
//! assert_eq!(albumin.map_source_to_target(58)?, Some(58));
//! assert_eq!(&seq["P02769"][57..58], "C");
//! assert_eq!(
//!     (albumin.identical, albumin.substituted, albumin.source_residues_facing_gap),
//!     (465, 142, 2)
//! );
//! assert_eq!(albumin.map_source_to_target(26)?, Some(26));
//! assert_eq!((&seq["P02768"][25..26], &seq["P02769"][25..26]), ("A", "T"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Always check the residue at the mapped position** before you carry a modification across: a
//! phosphosite on S is not a site on A. [`Alignments::residue_map`] has the whole map as a table,
//! one row per residue of both sequences, with both letters side by side in its `residue` and
//! `mapped_residue` columns.
//!
//! **Many pairs, one call.** [`align_with`] takes a list of `(source, target)` pairs and aligns
//! them in one bridge process, in order, so a few thousand isoform pairs cost one process start.
//! The residue map has a `pair_index` column, and [`Alignments::alignments`] gives one
//! [`Alignment`] per pair. For a large batch where you only want the scores and counts, set
//! [`AlignOptions::residue_map`] to `false`: the map has one row per residue of both sequences,
//! and leaving it out keeps the answer small.
//!
//! **The options, and whose defaults they are.** Every option left `None` is **mzLib's to set**.
//! Nothing is sent for it, and the result says what the aligner used: the isoform alignment above
//! passed no options, so its [`Alignments::matrix`], [`Alignments::gap_open`],
//! [`Alignments::gap_extend`], [`Alignments::end_gaps`] and [`Alignments::max_cells`] are mzLib's
//! own defaults in this build.
//!
//! | Option | Meaning |
//! |---|---|
//! | [`AlignOptions::gap_open`] | cost of opening a gap, in BLOSUM62's half-bits; a gap of *k* residues costs `gap_open + k * gap_extend` (BLAST's convention; mzLib's defaults are BLASTP's) |
//! | [`AlignOptions::gap_extend`] | cost of each residue in a gap, the first included |
//! | [`AlignOptions::end_gaps`] | [`EndGaps::Free`]: an overhang at either end costs nothing (EMBOSS needle's default), so a signal peptide or a longer tail is not penalised. [`EndGaps::Charged`]: every gap costs |
//! | [`AlignOptions::max_cells`] | the largest table mzLib will build, `(source length + 1) * (target length + 1)` cells of one byte; a larger pair is refused rather than allowed to exhaust memory |
//!
//! The matrix is always BLOSUM62, mzLib's default and the only one offered.
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::proteins::{align_with, AlignOptions, EndGaps};
//!
//! // mzLib's deletion and overhang cases, with every option left to mzLib.
//! let free = align_with(
//!     &[("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ"), ("MLLLPEPTIDE", "PEPTIDE")],
//!     &AlignOptions::default(),
//! )?;
//! assert_eq!(
//!     (free.matrix.as_str(), free.gap_open, free.gap_extend, free.end_gaps.as_str(), free.max_cells),
//!     ("BLOSUM62", 11, 1, "free", 200_000_000)
//! );
//! let charged = align_with(
//!     &[("MLLLPEPTIDE", "PEPTIDE")],
//!     &AlignOptions { end_gaps: Some(EndGaps::Charged), ..Default::default() },
//! )?;
//! // The four-residue overhang costs 11 + 4 * 1 only when it is charged.
//! assert_eq!(free.columns.integers("score")?[1], Some(39));
//! assert_eq!(charged.columns.integers("score")?, [Some(24)]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **The traps.**
//!
//! - **`score` is not an identity, a bit score or an E-value.** It is a raw sum of BLOSUM62 scores
//!   in half-bits, less the gaps, and it grows with length. Compare scores only between
//!   alignments with the same `aligner_id`. This crate computes no percent identity, because
//!   there are several definitions: divide `identical` by the length you mean.
//! - **It is a global alignment.** Every residue of both sequences is placed, paired or facing a
//!   gap. With free end gaps a short sequence is placed on the best-scoring stretch of a long one,
//!   even when it has no real match there.
//! - **One optimal alignment, chosen the same way every time.** Where two alignments score the
//!   same, mzLib's traceback picks one, always the same one. Isoform 2B above is a real case:
//!   canonical V76 and V82 are both V, either could pair with the isoform's V76, and mzLib pairs
//!   V82, leaving V76 to W81 facing the gap. Another aligner may give the equally good other
//!   answer.
//! - **Pass bare one-letter sequences.** mzLib refuses only `-` and `.`, as
//!   [`MzLibError::Usage`]. Any other character BLOSUM62 does not name (selenocysteine `U`,
//!   pyrrolysine `O`, a digit, a space, a modification bracket) is scored as an unknown residue
//!   `X`, silently. Strip modifications first.
//! - **Long pairs need a larger `max_cells`.** Time and memory grow with the product of the two
//!   lengths. Two titin isoforms need about 1.2 billion cells, past the default; mzLib accepts up
//!   to 2,147,483,591. Two identical sequences are aligned residue for residue without the table,
//!   so the limit does not apply to them.
//!
//! ## Read an accession into entry, isoform and variants
//!
//! **The problem.** The accessions a search reports are not all plain entries. UniProt isoforms
//! carry `-2`; RefSeq accessions carry a `.3` version; and mzLib names a protein with applied
//! sequence variants after its entry, `P04637_S70N_A80T` or `NP_000537.3_R72P`. To join results to
//! anything you keep per entry (genes, orthologs, GO terms) you need the entry back, without
//! losing the isoform or the variants. Cutting at the first `_` turns `NP_000537` into `NP`.
//!
//! [`parse_accessions`] reads each accession with mzLib's own grammars
//! (`VariantApplication.ParseAccession`, mzLib #1382):
//!
//! ```
//! # mzlib_replay::activate();
//! let parsed = mzlib::proteins::parse_accessions(&[
//!     "P04637", "P04637-2", "Q13409-2", "NP_000537.3_R72P", "P04637_S70N_A80T",
//!     "P04637_2", "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768",
//! ])?;
//! let c = &parsed.columns;
//! let (input, namespace, entry) = (c.strings("input")?, c.strings("namespace")?, c.strings("entry_accession")?);
//! let (isoform, version, variants) = (c.integers("isoform")?, c.integers("version")?, c.strings("applied_variants")?);
//! for i in 0..input.len() {
//!     println!("{:?} {:?} {:?} {:?} {:?} {:?}", input[i], namespace[i], entry[i], isoform[i], version[i], variants[i]);
//! }
//! //   "P04637"              uniprot       P04637               -  -  -
//! //   "P04637-2"            uniprot       P04637               2  -  -
//! //   "Q13409-2"            uniprot       Q13409               2  -  -
//! //   "NP_000537.3_R72P"    refseq        NP_000537            -  3  R72P
//! //   "P04637_S70N_A80T"    uniprot       P04637               -  -  S70N_A80T
//! //   "P04637_2"            unrecognized  P04637_2             -  -  -
//! //   "DECOY_P02768"        unrecognized  DECOY_P02768         -  -  -
//! //   "sp|P04637|P53_HUMAN" unrecognized  sp|P04637|P53_HUMAN  -  -  -
//! //   "P04637|P02768"       unrecognized  P04637|P02768        -  -  -
//! assert_eq!(namespace[5].as_deref(), Some("unrecognized"));
//! assert_eq!((version[3], variants[3].as_deref()), (Some(3), Some("R72P")));
//! assert_eq!(
//!     [parsed.namespace_counts["unrecognized"], parsed.namespace_counts["uniprot"], parsed.namespace_counts["refseq"]],
//!     [4, 4, 1]
//! );
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! What each column means:
//!
//! - **`entry_accession`** is the UniProt entry or the unversioned RefSeq accession. It answers
//!   "same entry?", **not** "same sequence?": `P04637` and `P04637-2` share it, but they are
//!   different sequences. Join per-entry data on it; join sequence-level data on `input`.
//! - **`isoform`** and **`version`** are the `-n` and `.n` numbers, `None` when there is none.
//! - **`parent_accession`** is the entry the variants were applied to, isoform or version kept:
//!   `NP_000537.3` for `NP_000537.3_R72P`. Look the parent sequence up under it.
//! - **`applied_variants`** is the variant suffix exactly as written, `None` when there is none.
//!   It is not checked against any sequence.
//! - **`namespace`** is one of [`NAMESPACES`], spelled as [`resolve_genes_with`] spells it, so the
//!   two tables join.
//!
//! [`ParsedAccessions::entry_of`] is the `input` to `entry_accession` map, for grouping.
//!
//! **It parses and never repairs.** Everything outside the two grammars comes back
//! `unrecognized`, with `entry_accession` equal to the input, so a data problem stays visible
//! instead of being fixed silently:
//!
//! - **`P04637_2` is not a variant of P04637.** It is mzLib's loader counter for a *different*
//!   entry whose accession collided with P04637's.
//! - **Decoy, contaminant and entrapment prefixes** (`DECOY_P02768`) and **FASTA headers**
//!   (`sp|P04637|P53_HUMAN`) are left as they are.
//! - **A protein group** (`P04637|P02768`) is not split. Split groups on `|` first, if that is
//!   what you mean.
//! - **Lower case and spaces count.** `" P04637"` is unrecognized: every accession is sent exactly
//!   as you give it. Only RefSeq *protein* accessions (`NP_`, `XP_`, `YP_`, `WP_`, `AP_`) are
//!   recognized, so a transcript (`NM_`) or an Ensembl id (`ENSP...`) is unrecognized too.
//! - **An isoform or version past 2,147,483,647** reads as `None`, while the namespace still says
//!   `uniprot` or `refseq`: mzLib keeps the entry and drops the number.
//!
//! ## Orthologs across species: an Ensembl Compara snapshot
//!
//! **The problem.** You have a mouse result and a human question, or identifications from three
//! species to line up: which human gene is this mouse gene's ortholog? Ensembl Compara answers
//! that, per species pair, in bulk files with three traps:
//!
//! 1. **The dumps partition.** Each species has its own homology dump, and Compara puts each of a
//!    pair's relationships in *one* of the pair's two dumps, arbitrarily. In release 116 every
//!    human-mouse orthology is in the **mouse** dump and none is in the human one. Reading the
//!    dump named after your species gives a complete-looking, badly incomplete answer.
//! 2. **Orthology is not transitive.** Human A ~ mouse B and mouse B ~ rat C does not make A ~ C.
//!    Chaining pairs into groups invents tuples Compara never asserted.
//! 3. **"No ortholog" means three things.** The gene's tree may hold a gene of the other species
//!    with no ortholog called, hold none of that species, or the gene may be in no tree at all.
//!
//! [`build_orthology_with`] builds mzLib's orthology snapshot (`OrthologySnapshot`, mzLib #1381),
//! which makes the first trap structural (**mzLib refuses a call without every species' dump**)
//! and keeps the other two visible. This example is mzLib's own test world, written to
//! `tests/fixtures/orthology/`: four human genes, three mouse and two rat, with every human-mouse
//! row in the mouse dump:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::proteins::{build_orthology_with, OrthologyOptions};
//!
//! let species = ["homo_sapiens", "mus_musculus", "rattus_norvegicus"];
//! let options = OrthologyOptions {
//!     gtfs: vec![
//!         ("homo_sapiens".into(), "Homo_sapiens.GRCh38.116.gtf".into()),
//!         ("mus_musculus".into(), "Mus_musculus.GRCm39.116.gtf".into()),
//!         ("rattus_norvegicus".into(), "Rattus_norvegicus.GRCr8.116.gtf".into()),
//!     ],
//!     dumps: species
//!         .iter()
//!         .map(|s| format!("homologies/{s}/Compara.116.protein_default.homologies.tsv.gz").into())
//!         .collect(),
//!     gene_trees: "vertebrates.GeneTree_content.default.e116.txt.gz".into(),
//!     ..Default::default()
//! };
//! let snap = build_orthology_with(
//!     "116",
//!     &OrthologyOptions {
//!         // Every human gene, in mouse.
//!         pair_status: Some(("homo_sapiens".into(), "mus_musculus".into())),
//!         ..options.clone()
//!     },
//! )?;
//! let status = snap.pair_status.as_ref().expect("pair_status was asked for");
//! let (ids, statuses) = (status.columns.strings("gene_id")?, status.columns.strings("status")?);
//! for (id, status) in ids.iter().zip(&statuses) {
//!     println!("{} {}", id.as_deref().unwrap_or(""), status.as_deref().unwrap_or(""));
//! }
//! //   ENSG00000000001 has_ortholog
//! //   ENSG00000000002 has_ortholog
//! //   ENSG00000000003 tree_lacks_target_species
//! //   ENSG00000000004 not_in_any_tree
//! assert_eq!(statuses[3].as_deref(), Some("not_in_any_tree"));
//! assert_eq!(&snap.snapshot_id[..12], "5487ebcabe48");
//!
//! // Every species' dump, always: the human dump holds one row (human-rat); all three human-mouse
//! // rows came from the mouse dump.
//! let kept: Vec<(&str, u64)> = snap.dumps.iter().map(|(s, d)| (s.as_str(), d.rows_kept)).collect();
//! assert_eq!(kept, [("homo_sapiens", 1), ("mus_musculus", 5), ("rattus_norvegicus", 1)]);
//! for p in &snap.pairs {
//!     println!("{} {} {} {}", p.species_a, p.species_b, p.ortholog_count, p.paralog_count);
//! }
//! //   homo_sapiens homo_sapiens 0 0
//! //   homo_sapiens mus_musculus 3 0
//! //   homo_sapiens rattus_norvegicus 1 0
//! //   mus_musculus mus_musculus 0 1
//! //   mus_musculus rattus_norvegicus 2 0
//! //   rattus_norvegicus rattus_norvegicus 0 0
//! assert_eq!(snap.pair("mus_musculus", "homo_sapiens")?.row_count, 3);   // either order
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Record [`OrthologySnapshot::snapshot_id`] with any result you draw from the snapshot. It is
//! mzLib's sha256 over the format, the release, the species and every input's sha256, so the same
//! inputs always give the same id, and a different id means different inputs.
//!
//! Leave the mouse dump out and mzLib refuses: *"No homology dump for mus_musculus. Each pair's
//! relationships are split between its two species' dumps, so every species' dump is required."*
//! That is [`MzLibError::Usage`]. A pair with nothing in it, such as the human paralogs here,
//! still has its entry.
//!
//! **What to download, once per release.** Nothing is downloaded for you. For release 116, from
//! `https://ftp.ensembl.org/pub/release-116/`:
//!
//! | Input | Path under the release | Size (human, 116) | Option |
//! |---|---|---|---|
//! | Each species' primary-assembly GTF | `gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz` | 141 MB | [`OrthologyOptions::gtfs`] |
//! | Each species' homology dump | `tsv/ensembl-compara/homologies/homo_sapiens/Compara.116.protein_default.homologies.tsv.gz` | 109 MB | [`OrthologyOptions::dumps`] |
//! | The gene trees, one file for all species | `compara/vertebrates.GeneTree_content.default.e116.txt.gz` | 79 MB | [`OrthologyOptions::gene_trees`] |
//!
//! - **The species names are Ensembl production names**: `homo_sapiens`, `mus_musculus`,
//!   `rattus_norvegicus`, exactly as the dumps' `species` column spells them.
//! - **Take the dump in each species' own folder.** Every species' dump has the same file name, so
//!   keep them in their folders. One folder up, `homologies/` also holds a 10 GB file of the same
//!   name; it is not one species' dump.
//! - **Keep Ensembl's file names.** mzLib holds the gene trees and every dump to your `release`
//!   only when their names carry it (`.e116.txt.gz`, `Compara.116.`). A renamed file is not
//!   checked, and a caveat says which.
//! - **A GTF's own number is not compared.** Ensembl publishes some species with their own
//!   annotation release. mzLib holds a gene set to the dumps by content instead: a homology that
//!   names a gene the gene set lacks is refused.
//! - **A compact gene table works too.** If you made one for [`resolve_genes_with`] with mzLib's
//!   `EnsemblGeneSetWriter`, pass it in [`OrthologyOptions::gene_sets`] instead. It carries the
//!   GTF's provenance, so the snapshot and its id are the same.
//!
//! **What each status means.** [`OrthologyOptions::pair_status`] gives every gene of
//! `from_species`'s gene set exactly one status in `to_species`. It says how far Compara's
//! inference got, not what biology holds:
//!
//! | status | what it means | read it as |
//! |---|---|---|
//! | `has_ortholog` | at least one ortholog in `to_species` | an answer; the Parquet has the rows |
//! | `no_edge_in_shared_tree` | the gene's tree holds a `to_species` gene, and Compara called no ortholog to it | the nearest the data comes to "no ortholog exists" |
//! | `tree_lacks_target_species` | the gene's tree holds no `to_species` gene | nothing was compared |
//! | `not_in_any_tree` | the gene is in the gene set but in no gene tree | no call was attempted |
//!
//! [`PairStatus::status_counts`] counts genes per status, every status present, and
//! [`PairStatus::genes`] lists one status's genes. **A status is not symmetric**: ask for
//! `("mus_musculus", "homo_sapiens")` separately. **Every biotype is a row**, so choose your
//! denominator by filtering on `gene_biotype`.
//!
//! **Write the snapshot, then query it.** The relationships themselves go to disk. Give
//! [`OrthologyOptions::out`] a directory that does not exist yet:
//!
//! ```
//! # mzlib_replay::activate();
//! # use mzlib::proteins::{build_orthology_with, OrthologyOptions};
//! # let species = ["homo_sapiens", "mus_musculus", "rattus_norvegicus"];
//! # let options = OrthologyOptions {
//! #     gtfs: vec![
//! #         ("homo_sapiens".into(), "Homo_sapiens.GRCh38.116.gtf".into()),
//! #         ("mus_musculus".into(), "Mus_musculus.GRCm39.116.gtf".into()),
//! #         ("rattus_norvegicus".into(), "Rattus_norvegicus.GRCr8.116.gtf".into()),
//! #     ],
//! #     dumps: species.iter()
//! #         .map(|s| format!("homologies/{s}/Compara.116.protein_default.homologies.tsv.gz").into())
//! #         .collect(),
//! #     gene_trees: "vertebrates.GeneTree_content.default.e116.txt.gz".into(),
//! #     ..Default::default()
//! # };
//! let written = build_orthology_with(
//!     "116",
//!     &OrthologyOptions {
//!         restrict_gene_trees: true,
//!         out: Some("compara-116".into()),
//!         ..options
//!     },
//! )?;
//! assert!(written.snapshot_id.starts_with("5487ebcabe48"));   // the same id as above
//! let files = &written.written.as_ref().expect("out was given").files;
//! for f in &files[..7] {
//!     println!("{} {:?}", f.path, f.rows);
//! }
//! //   genes/homo_sapiens.parquet Some(4)
//! //   genes/mus_musculus.parquet Some(3)
//! //   genes/rattus_norvegicus.parquet Some(2)
//! //   members/homo_sapiens.parquet Some(3)
//! //   members/mus_musculus.parquet Some(3)
//! //   members/rattus_norvegicus.parquet Some(2)
//! //   pairs/homo_sapiens__homo_sapiens.parquet Some(0)
//! let last = files.last().unwrap();
//! assert_eq!((last.path.as_str(), last.rows), ("views.sql", None));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! The directory holds `genes/SPECIES.parquet` (the gene set), `members/SPECIES.parquet` (gene to
//! Compara tree), `pairs/A__B.parquet` (every relationship between A and B, A sorting first;
//! A == B holds the paralogs), `views.sql` and `manifest.json` (every input and file with its
//! sha256, at [`WrittenSnapshot::manifest`]). The same inputs always give the same bytes. **A
//! directory that is not empty is refused**: a snapshot is written once, and next release's goes
//! in a new directory to compare with this one. [`OrthologyOptions::restrict_gene_trees`] keeps
//! only your species' genes from the gene-tree file; without it, mzLib keeps every gene in it and
//! a caveat says how many. The snapshot is the same either way, as the equal ids show.
//!
//! **mzLib writes Parquet and ships no reader, and this crate reads none.** Query the files with
//! DuckDB (its command-line shell, or the `duckdb` crate), which runs the macros mzLib writes in
//! `views.sql`. These are the answers for the test world above:
//!
//! ```sql
//! .read compara-116/views.sql
//!
//! -- Every human protein-coding gene's status in mouse: pair_status, from the files.
//! SELECT status, count(*) AS genes
//! FROM pair_status('compara-116', 'homo_sapiens', 'mus_musculus')
//! WHERE gene_biotype = 'protein_coding'
//! GROUP BY status ORDER BY status;
//! -- has_ortholog 2 | not_in_any_tree 1 | tree_lacks_target_species 1
//!
//! -- The orthologs, oriented from human: one row per relationship, one-to-many kept as several.
//! SELECT gene, partner_gene, relationship_type, identity, partner_identity
//! FROM orthologs('compara-116', 'homo_sapiens', 'mus_musculus') ORDER BY gene, partner_gene;
//! -- ENSG00000000001 ENSMUSG00000000001 ortholog_one2one  91.0 90.0
//! -- ENSG00000000002 ENSMUSG00000000002 ortholog_one2many 71.0 70.0
//! -- ENSG00000000002 ENSMUSG00000000003 ortholog_one2many 61.0 60.0
//!
//! -- One gene per species where all three pairs are orthologs; never chained.
//! SELECT * FROM species_set3('compara-116', 'homo_sapiens', 'mus_musculus', 'rattus_norvegicus');
//! -- ENSG00000000001 ENSMUSG00000000001 ENSRNOG00000000001 true
//! ```
//!
//! H2 is orthologous to mouse M2, and M2 to rat R2, but H2 ~ R2 was never called, so
//! `species_set3` does not return H2. That is trap 2, handled.
//!
//! **What the Parquet holds.** `pairs/A__B.parquet` has one row per relationship, ordered by
//! Compara's `homology_id`. Values are Compara's, verbatim, and **a NULL stays null, never 0 or
//! false**: `identity_a` and `identity_b` (percent; the share of **that side's** protein identical
//! to the other's, so the two differ), `goc_score` (0 to 100), `wga_coverage` (percent), `dn` and
//! `ds` (substitutions per site; null for most comparisons), `is_high_confidence` (null: Compara
//! declines to say, which is not false), and `source_dump` (the **file name** of the dump the row
//! came from: provenance, not meaning). `members/SPECIES.parquet`'s `group_id` is `compara-116:`
//! plus the tree id, so trees are never compared across releases. The full column lists, with
//! types and units, are on [`build_orthology_with`].
//!
//! **The traps.**
//!
//! - **Never chain pairs.** Use `species_set3` (or mzLib's `OrthologySnapshot.SpeciesSet` for more
//!   species), which requires every two members of a tuple to be orthologs.
//! - **One-to-many is several rows.** Nothing is picked: a human gene with two mouse orthologs has
//!   two rows, both `ortholog_one2many`. Decide per question.
//! - **The gene set and the dumps must be one release.** A homology naming a gene its gene set
//!   lacks is refused as [`MzLibError::Bridge`] (`InvalidDataException`): that is what a GTF and
//!   dumps from different releases look like. Use the primary-assembly GTF, as for
//!   [`resolve_genes_with`].
//! - **`source_dump` cannot tell species apart.** Ensembl gives every species' dump the same file
//!   name, so the column holds that one name everywhere, and a caveat says so. `manifest.json`'s
//!   `inputs` list each dump with its species and sha256.
//! - **Join with [`resolve_genes_with`] on the gene set.** Its `gene_set_sha256` column and
//!   [`OrthologySnapshot::gene_sets`]`[species].sha256` are the same value when both were counted
//!   against the same GTF. Check that before joining a resolution to a snapshot.
//! - **It is slow and large.** Every input is read whole. For human, mouse and rat at release 116,
//!   mzLib reports about 20 s to build and 40 s more to write 13 files (9.0 MB). Leave the timeout
//!   at `None` and set `restrict_gene_trees` when memory is short.
//!
//! ## Map a residue across species: Compara's gene-tree alignment
//!
//! **The problem.** You have a site in a human protein, human p53's phosphoserine 15, and a mouse
//! experiment. Which mouse residue is it? Counting along the two sequences fails as soon as one has
//! an insertion. Aligning just the two proteins works ([`align_with`]), but it is your alignment,
//! made with your parameters. Ensembl Compara has already aligned every protein of each gene tree,
//! across all its species, to build the tree. A residue's column in that alignment is the published
//! answer.
//!
//! [`map_residues_with`] reads Compara's peptide alignments (`ComparaGeneTreeAlignment.Load`,
//! mzLib #1465), keeping only the canonical proteins of your species' genes, and answers each query
//! with mzLib's `MapResidue`. This example is a real cut of Ensembl release 116
//! (`tests/fixtures/map_residues/`): the TP53 and GAPDH gene trees' alignments, with the human,
//! mouse and rat rows exactly as Ensembl publishes them. The same queries on the full files give the
//! same rows.
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::proteins::{map_residues_with, MapResiduesOptions, COLUMN_OUTCOMES};
//!
//! let gtfs = vec![
//!     ("homo_sapiens".to_owned(), "Homo_sapiens.GRCh38.116.gtf".into()),
//!     ("mus_musculus".to_owned(), "Mus_musculus.GRCm39.116.gtf".into()),
//!     ("rattus_norvegicus".to_owned(), "Rattus_norvegicus.GRCr8.116.gtf".into()),
//! ];
//! let (tp53, trp53, rat_tp53) = ("ENSP00000269305", "ENSMUSP00000104298", "ENSRNOP00000074031");
//! let queries = [
//!     (tp53, 15, trp53), (tp53, 20, trp53), (tp53, 46, trp53), (tp53, 392, trp53),
//!     (trp53, 18, tp53),                     // and back
//!     (tp53, 15, rat_tp53),
//!     (tp53, 15, "ENSMUSP00000113213"),      // mouse GAPDH: another gene tree
//!     ("ENSP00000352610", 15, trp53),        // a TP53 isoform Compara did not align
//!     (tp53, 15, "ENSMUSP00000005371"),      // a Trp53 isoform, likewise
//! ];
//! let mapped = map_residues_with(
//!     &queries,
//!     &MapResiduesOptions {
//!         alignment: "Compara.116.protein_default.aa.fasta.gz".into(),
//!         gene_trees: "vertebrates.GeneTree_content.default.e116.txt.gz".into(),
//!         gtfs,
//!         proteins: true,
//!         ..Default::default()
//!     },
//! )?;
//! for r in mapped.rows()? {
//!     println!("{} {} {} {} {:?} {:?} {:?}", r.source_protein_id, r.position, r.target_protein_id,
//!              r.outcome, r.column, r.target_position, r.target_residue);
//! }
//! //   ENSP00000269305 15 ENSMUSP00000104298 Aligned Some(476) Some(18) Some("S")
//! //   ENSP00000269305 20 ENSMUSP00000104298 Aligned Some(547) Some(23) Some("S")
//! //   ENSP00000269305 46 ENSMUSP00000104298 GapInTarget Some(863) None None
//! //   ENSP00000269305 392 ENSMUSP00000104298 Aligned Some(3036) Some(389) Some("S")
//! //   ENSMUSP00000104298 18 ENSP00000269305 Aligned Some(476) Some(15) Some("S")
//! //   ENSP00000269305 15 ENSRNOP00000074031 Aligned Some(476) Some(15) Some("S")
//! //   ENSP00000269305 15 ENSMUSP00000113213 DifferentAlignments None None None
//! //   ENSP00000352610 15 ENSMUSP00000104298 SourceNotInAlignment None None None
//! //   ENSP00000269305 15 ENSMUSP00000005371 TargetNotInAlignment None None None
//! let rows = mapped.rows()?;
//! assert_eq!(
//!     rows.iter().map(|r| r.target_position).collect::<Vec<_>>(),
//!     [Some(18), Some(23), None, Some(389), Some(15), Some(15), None, None, None]
//! );
//! assert_eq!((rows[2].outcome.as_str(), rows[2].column), ("GapInTarget", Some(863)));
//!
//! // Every outcome is counted, zeros included.
//! let counts: Vec<u64> = COLUMN_OUTCOMES.iter().map(|o| mapped.outcome_counts[*o]).collect();
//! assert_eq!(counts, [5, 1, 1, 1, 1]);
//! assert_eq!(
//!     (mapped.alignment.release.as_deref(), mapped.alignment.protein_count, mapped.gene_trees.member_count),
//!     (Some("116"), 5, 5)
//! );
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Reading it row by row:
//!
//! - **Human Ser15 is mouse Ser18**, Ser20 is Ser23, and Ser392 is Ser389. The mouse protein is
//!   three residues shorter near its N-terminus, so the numbers shift, and the column says by how
//!   much.
//! - **Human Ser46 faces a gap in mouse.** No mouse residue sits in column 863, and the outcome is
//!   `GapInTarget`. It never offers the nearest residue instead.
//! - **Mapping back gives the residue you started from.** Mouse residue 18 is in the same column,
//!   476.
//! - **Rat is a separate query** against the same alignment: human Ser15 is rat Ser15.
//! - **GAPDH is in another gene tree**, so the alignment says nothing about it:
//!   `DifferentAlignments`.
//! - **Only each gene's canonical protein is aligned.** `ENSP00000352610` is a real TP53 protein,
//!   but not the one Compara used, so it is `SourceNotInAlignment`. It is not mapped through the
//!   canonical protein. The gene-tree file names each gene's canonical protein, in the row flagged
//!   `Y`.
//!
//! | outcome | what it means | `column` | `target_position`, `target_residue` |
//! |---|---|---|---|
//! | `Aligned` | the target has a residue in the source residue's column | set | set |
//! | `GapInTarget` | the target has a gap in that column: no residue corresponds | set | `None` |
//! | `DifferentAlignments` | the two proteins are in different gene trees | `None` | `None` |
//! | `SourceNotInAlignment` | the source protein was not kept: not in the file, not canonical, or outside your gene sets | `None` | `None` |
//! | `TargetNotInAlignment` | the target protein was not kept | `None` | `None` |
//!
//! The outcome names are mzLib's own (`ComparaColumnOutcome`), and [`COLUMN_OUTCOMES`] lists them in
//! mzLib's order.
//!
//! ### From a UniProt position to a UniProt position
//!
//! Compara's positions are positions in **Ensembl** proteins. A site from a search against UniProt
//! needs two more steps around [`map_residues_with`], each an alignment:
//!
//! 1. Align the UniProt sequence to the canonical Ensembl protein, and map the position across.
//! 2. Map the Ensembl position through the gene-tree column with [`map_residues_with`].
//! 3. Align the target's canonical Ensembl protein to its UniProt sequence, and map across again.
//!
//! Any of the three can meet a gap. [`MapResiduesOptions::proteins`] returned the Ensembl
//! sequences, so no second file is needed for them. The UniProt sequences come from UniProt's FASTA
//! for human p53 (P04637) and mouse p53 (P02340):
//!
//! ```
//! # mzlib_replay::activate();
//! # use mzlib::proteins::{map_residues_with, MapResiduesOptions};
//! # let gtfs = vec![
//! #     ("homo_sapiens".to_owned(), "Homo_sapiens.GRCh38.116.gtf".into()),
//! #     ("mus_musculus".to_owned(), "Mus_musculus.GRCm39.116.gtf".into()),
//! #     ("rattus_norvegicus".to_owned(), "Rattus_norvegicus.GRCr8.116.gtf".into()),
//! # ];
//! # let (tp53, trp53, rat_tp53) = ("ENSP00000269305", "ENSMUSP00000104298", "ENSRNOP00000074031");
//! # let queries = [
//! #     (tp53, 15, trp53), (tp53, 20, trp53), (tp53, 46, trp53), (tp53, 392, trp53), (trp53, 18, tp53),
//! #     (tp53, 15, rat_tp53), (tp53, 15, "ENSMUSP00000113213"), ("ENSP00000352610", 15, trp53),
//! #     (tp53, 15, "ENSMUSP00000005371"),
//! # ];
//! # let options = MapResiduesOptions {
//! #     alignment: "Compara.116.protein_default.aa.fasta.gz".into(),
//! #     gene_trees: "vertebrates.GeneTree_content.default.e116.txt.gz".into(),
//! #     gtfs,
//! #     ..Default::default()
//! # };
//! # let mapped = map_residues_with(&queries, &MapResiduesOptions { proteins: true, ..options.clone() })?;
//! use mzlib::proteins::{align_with, read_with, AlignOptions, ProteinReadOptions};
//!
//! let uniprot = read_with(
//!     &["P04637.fasta", "P02340.fasta"],
//!     &ProteinReadOptions { sequences: true, ..Default::default() },
//! )?;
//! let sequences = uniprot.columns.strings("sequence")?;
//! let (p04637, p02340) = (sequences[0].clone().unwrap(), sequences[1].clone().unwrap());
//! let ensembl = mapped.proteins.as_ref().expect("proteins: true").sequences()?;
//! let legs = align_with(
//!     &[
//!         (p04637.as_str(), ensembl[tp53].as_str()),  // leg 1: UniProt human -> Ensembl human
//!         (ensembl[trp53].as_str(), p02340.as_str()), // leg 3: Ensembl mouse -> UniProt mouse
//!     ],
//!     &AlignOptions::default(),
//! )?
//! .alignments()?;
//! let (leg1, leg3) = (&legs[0], &legs[1]);
//! assert_eq!((leg1.identical, leg1.source_length, leg3.identical, leg3.target_length), (393, 393, 390, 390));
//!
//! let ensembl_position = leg1.map_source_to_target(15)?.expect("no gap in leg 1");
//! let site = map_residues_with(&[(tp53, ensembl_position as i64, trp53)], &options)?;
//! let mouse_position = site.rows()?[0].target_position.expect("Aligned") as usize;
//! let uniprot_mouse = leg3.map_source_to_target(mouse_position)?.expect("no gap in leg 3");
//! assert_eq!((uniprot_mouse, &p02340[uniprot_mouse - 1..uniprot_mouse]), (18, "S"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! For p53 the UniProt and Ensembl canonical sequences are identical in both species, so legs 1 and
//! 3 map every residue to itself. That is common for canonical sequences, but check it every time
//! (`identical` equal to both lengths) rather than assume it. For an isoform, or an entry whose
//! sequence UniProt and Ensembl annotate differently, those legs shift positions or meet gaps, as
//! the [alignment section](#map-a-residue-from-one-sequence-onto-another) shows. **Check the residue
//! at every step.** A column pairs positions the alignment puts together, which is homology, not
//! identity. Here it is serine on both sides. Where the target residue is a different amino acid,
//! the site does not carry over: a phosphosite on S is not a site on A.
//!
//! **What to download, once per release.** Nothing is downloaded for you. For release 116, from
//! `https://ftp.ensembl.org/pub/release-116/`:
//!
//! | Input | Path under the release | Option |
//! |---|---|---|
//! | The peptide alignments of every gene tree | `emf/ensembl-compara/homologies/Compara.116.protein_default.aa.fasta.gz` | [`MapResiduesOptions::alignment`] |
//! | The gene trees, one file for all species | `compara/vertebrates.GeneTree_content.default.e116.txt.gz` | [`MapResiduesOptions::gene_trees`] |
//! | Each species' primary-assembly GTF | `gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz` | [`MapResiduesOptions::gtfs`] |
//!
//! - **Use the alignment and the gene trees of one release.** mzLib does not compare them, and a
//!   gene's canonical protein and tree can change between releases. Keep Ensembl's file names: the
//!   result reports the release each name carries, and a caveat says when they differ or one is
//!   missing.
//! - **Give your species' gene sets.** They are used only to choose what is kept. Without them,
//!   mzLib's default keeps every gene in the gene-tree file, millions for the vertebrates, and every
//!   canonical protein in the alignment, and a caveat says so. The answers for your species'
//!   proteins are the same either way. A compact gene table ([`MapResiduesOptions::gene_sets`])
//!   works as it does for [`build_orthology_with`].
//! - **Batch your queries.** Reading the alignment is the cost: for human, mouse and rat at release
//!   116, the whole call took about a minute. Each query after that is a lookup, so send them all
//!   in one call.
//!
//! **The traps.**
//!
//! - **A position past the end of the source protein fails the whole call.** mzLib's `MapResidue`
//!   throws for it rather than returning an outcome, so the call returns [`MzLibError::Usage`] with
//!   mzLib's message (*"ENSP00000269305 has 393 residues."*) and no row, not even for the queries
//!   that were fine. This crate does not check lengths itself. Keep each position within its
//!   protein: with [`MapResiduesOptions::proteins`], a sequence's length is the last valid position.
//! - **Use Ensembl ids and Ensembl positions.** A UniProt accession is `SourceNotInAlignment`, never
//!   an error. Cross over with [`align_with`] first, as above.
//! - **Sharing a gene tree is not orthology.** A tree holds a gene family, paralogs included, so an
//!   `Aligned` answer between two proteins of one tree does not make them orthologs. Ask
//!   [`build_orthology_with`] which genes are orthologs, then map their residues here.
//! - **`DifferentAlignments` is not "unrelated".** It says only that this alignment cannot answer.
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
//! - **GO annotation takes exactly one database**, because every row carries its sha256. A
//!   contaminant member whose entry is not in that XML reads `no_entry`, not `no_go_terms`.
//!
//! Wire verbs: `proteins read`, `genes resolve`, `genes orthology`, `genes map-residues`,
//! `proteins classify-peptides`, `proteins annotate-go`, `proteins update-go`, `proteins align` and
//! `proteins parse-accession`.
//!
//! ## Cite
//!
#![doc = include_str!("../docs/reference/cite.proteins.md")]

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

// ---------------------------------------------------------------------------------------------
// Gene Ontology: annotate stored protein groups, and fetch a release on purpose
// ---------------------------------------------------------------------------------------------

/// Every `annotation_status` [`annotate_go_with`] can report, as mzLib writes it
/// (`GoAnnotationTsv.StatusName`). A row with a GO term is always `annotated`; the other three are
/// the single term-less row of a group that has no term, saying why.
pub const ANNOTATION_STATUSES: [&str; 4] = ["annotated", "no_go_terms", "no_entry", "contaminant"];

/// How [`annotate_go_with`] annotates: the ontology to propagate over, and what to return.
///
/// `go_obo` is required and has no default: a GO annotation means something only relative to one
/// release, so the release is always a file you chose to keep. [`update_go`] fetches one.
#[derive(Debug, Clone, Default)]
pub struct GoAnnotateOptions {
    /// A go.obo file. It must exist; nothing is downloaded. Fetch a release on purpose with
    /// [`update_go`], keep the file, and record [`GoAnnotations::go`] with your results.
    pub go_obo: PathBuf,
    /// Your own term-to-category map, in mzLib's format (`#!category_map_format 1`, `#!map_name`,
    /// `#!map_version`, then `category`, `subcategory` and `anchor_go_id` columns). mzLib ships no
    /// vocabulary. Adds [`GoAnnotations::categories`]. `None` (default): no categories.
    pub category_map: Option<PathBuf>,
    /// `false` (default) fails when the database cites a GO id the release lacks — usually a
    /// UniProt release newer than the go.obo — naming every missing id. `true` drops each such id
    /// and lists it in [`GoAnnotations::unresolved_go_ids`], so one new term does not cost the run.
    pub skip_unknown_go_ids: bool,
    /// Write the **whole** table here with mzLib's own `GoAnnotationTsv` writer, provenance header
    /// included, whatever `offset` and `limit` are. Must end in `.tsv`, must not name an input,
    /// and its folder must exist. For a large run, pair it with `limit: Some(0)`.
    pub out: Option<PathBuf>,
    /// Write the category table here with mzLib's `GoCategoryTsv` writer. Must end in `.tsv`, and
    /// needs `category_map`.
    pub categories_out: Option<PathBuf>,
    /// Rows to skip before the returned window, in rows (not groups). Default 0.
    pub offset: u64,
    /// Return at most this many rows (rows, not groups) in [`GoAnnotations::columns`]. `None`
    /// (default) returns every row; `Some(0)` returns only the summary. Never shortens `out`.
    pub limit: Option<u64>,
    /// Time to allow. `None` (default) waits: a whole proteome XML takes a while.
    pub timeout: Option<Duration>,
}

/// The go.obo a result was computed against. Record it with your results.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct GoRelease {
    /// The file name read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub source_file_name: String,
    /// Lower-case hex sha256 of the file's bytes. This, not the release, is what proves two runs
    /// used the same ontology.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sha256: String,
    /// The file's `data-version`, e.g. `"releases/2026-07-26"`. `None`: the file has no
    /// `data-version` header.
    #[serde(default)]
    pub release: Option<String>,
    /// Terms in the file, obsolete terms included (alternative ids are not counted).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub term_count: u64,
}

/// The protein database whose GO terms annotated the groups.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct AnnotationDatabase {
    /// The path read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `"UniProtXml"`, the only type [`annotate_go_with`] accepts.
    #[serde(default)]
    pub file_type: Option<String>,
    /// The mzLib loader used, `"ProteinDbLoader.LoadProteinXML"`.
    #[serde(default)]
    pub reader: Option<String>,
    /// Entries loaded, as written: no decoys generated, no variants expanded.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub protein_count: u64,
    /// Lower-case hex sha256 of the **decompressed** database, so `.xml` and `.xml.gz` agree.
    /// Every row's `annotation_db_sha256` is this.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sha256: String,
    /// What mzLib did to the file while loading it, e.g. genotype variants applied.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// A file [`annotate_go_with`] wrote, when asked to.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct WrittenTable {
    /// The path written.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Rows written: the **whole** table, never the `limit`/`offset` window. `None` for the
    /// category table, whose rows are counted in [`GoCategories::row_count`].
    #[serde(default)]
    pub row_count: Option<u64>,
}

/// Your category map applied to the annotated terms: one row per (term, category, subcategory).
///
/// A term belongs to a category when one of that category's anchors is the term itself or one of
/// its ancestors. A term under no anchor has **no row**, so absence means "outside your map". Join
/// to [`GoAnnotations`] on `go_id`: both tables come from the same go.obo release.
///
/// Columns: `go_id` (an annotated term at or below one of the map's anchors), `category` (the
/// label, as your map writes it) and `subcategory` (`"category:subcategory"`, the most specific
/// anchor within the category; `None` where the term reaches only the category's own anchor).
#[derive(Debug, Clone, Deserialize)]
pub struct GoCategories {
    /// The map's declared `map_name`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub map_name: String,
    /// The map's declared `map_version`. Yours to bump; mzLib only records it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub map_version: String,
    /// The map file read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub source_file_name: String,
    /// Lower-case hex sha256 of the map file. A version is a claim; the hash is proof.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sha256: String,
    /// Rows (anchors) in the map.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub anchor_count: u64,
    /// Rows in the table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table; see the type's docs for its columns.
    #[serde(flatten)]
    pub columns: Table,
}

/// What [`annotate_go_with`] returns: one row per (protein group, GO term) that **any** member of
/// the group holds, directly or through an ancestor.
///
/// No member of a group is privileged — MetaMorpheus never picks a leading protein — so the table
/// is the union, and every row says which members carry its term. The views people usually want
/// are filters you apply: **consensus** is `n_with == n_members`; **direct annotations only** is
/// `propagated == Some(false)`; **leave out isoform inheritance** is `inherited == Some(false)`.
///
/// The column names are **mzLib's own** `GoAnnotationTsv` schema, in its order, so this table and
/// the file `out` writes (or mzLib writes inside a search) mean the same thing:
///
/// | column | read with | meaning; `None` means |
/// |---|---|---|
/// | `protein_group` | [`Table::strings`] | the group, as MetaMorpheus named it (`P0C0S5\|Q71UI9`) |
/// | `accession_used` | [`string_lists`] | members carrying the term, directly or by propagation |
/// | `accession_direct` | [`string_lists`] | of those, members annotated to this exact term |
/// | `accession_inherited` | [`string_lists`] | of those, members whose terms were borrowed: an isoform (`P04406-2`) or sequence variant (`P04406_A20T`) absent from the database takes its entry's terms |
/// | `go_id` | [`Table::strings`] | the term's primary id; `None`: a term-less row |
/// | `go_name` | [`Table::strings`] | the term's name in this release; `None`: a term-less row |
/// | `aspect` | [`Table::strings`] | `biological_process`, `cellular_component`, `molecular_function` or `unknown`; `None`: a term-less row |
/// | `evidence` | [`string_lists`] | evidence codes (ECO) pooled over the carrying members |
/// | `evidence_by_member` | [`string_list_maps`] | member → its own evidence codes for this term; empty on a term-less row |
/// | `inherited` | [`Table::booleans`] | every carrying member's terms were borrowed; `None`: a term-less row |
/// | `propagated` | [`Table::booleans`] | no member is annotated to this exact term — it is implied by a more specific one; `None`: a term-less row |
/// | `n_members` | [`Table::integers`] | members (proteins) in the group |
/// | `n_with` | [`Table::integers`] | members (proteins) carrying the term; 0 on a term-less row |
/// | `entrapment_members` | [`string_lists`] | the group's entrapment members, on every row of the group |
/// | `annotation_status` | [`Table::strings`] | one of [`ANNOTATION_STATUSES`] |
/// | `q_value` | [`Table::floats`] | the group's q-value (a fraction, 0 to 1) from the table; never filtered on |
/// | `go_release` | [`Table::strings`] | the go.obo `data-version`; `None`: the file has none |
/// | `go_obo_sha256` | [`Table::strings`] | sha256 of the go.obo |
/// | `annotation_db_sha256` | [`Table::strings`] | sha256 of the decompressed annotation database |
#[derive(Debug, Clone, Deserialize)]
pub struct GoAnnotations {
    /// The protein-group table read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub groups_file: String,
    /// Lower-case hex sha256 of that table's bytes; also [`Self::header`]'s `source_file_sha256`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub groups_file_sha256: String,
    /// Rows in the protein-group table, decoys included.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub table_row_count: u64,
    /// Of those, decoy groups, which are skipped: decoys carry no GO. A group is a decoy when any
    /// `|`-joined parent label is `D` or the entrapment decoy `ED` (mzLib #1371).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub decoy_group_count: u64,
    /// Groups annotated: every non-decoy group, contaminants included, each once even when
    /// MetaMorpheus wrote it twice (mzLib #1366).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub group_count: u64,
    /// The database the terms came from.
    pub annotation_database: AnnotationDatabase,
    /// The ontology release the terms were propagated over.
    pub go: GoRelease,
    /// Whether [`GoAnnotateOptions::skip_unknown_go_ids`] was set.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub skip_unknown_go_ids: bool,
    /// GO ids the database cites that the release lacks, dropped because
    /// [`GoAnnotateOptions::skip_unknown_go_ids`] was set, in ordinal order. Always empty otherwise
    /// (the call fails instead).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub unresolved_go_ids: Vec<String>,
    /// mzLib's `#!key value` header, as its writer produced it: `go_annotation_format`,
    /// `mzlib_version`, `mzlib_release`, `go_release`, every sha256, and the counters
    /// `n_multi_member_groups` and one `status_<status>` per status. **The counters count groups,
    /// not rows, and only groups at q <= `counter_q_value_max` (0.01)**; the rows themselves are not
    /// filtered. Every value is a string, as in the file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub header: BTreeMap<String, String>,
    /// Rows in the whole table.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Rows in [`Self::columns`]: the `limit`/`offset` window.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// The offset applied, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// True whenever rows were left out of [`Self::columns`] by `limit` or `offset`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// What this result cannot tell you: ids dropped from the database, a go.obo with no release,
    /// what mzLib did loading the database.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
    /// Where the whole table was written; `None`: `out` was not given.
    #[serde(default)]
    pub written: Option<WrittenTable>,
    /// Your category map applied to these terms; `None`: `category_map` was not given.
    #[serde(default)]
    pub categories: Option<GoCategories>,
    /// Where the category table was written; `None`: `categories_out` was not given.
    #[serde(default)]
    pub categories_written: Option<WrittenTable>,
    /// The table, for the returned window; see the type's docs for its columns.
    #[serde(flatten)]
    pub columns: Table,
}

/// What [`update_go`] returns: which go.obo is now on disk, and whether it changed.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct GoUpdate {
    /// The path written.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub go_obo_file: String,
    /// Where it was fetched from: GO's PURL, `https://purl.obolibrary.org/obo/go.obo`, which
    /// always serves the **current** release.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub url: String,
    /// Whether a file was already at the path.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub existed_before: bool,
    /// That file's sha256; `None`: no file was there.
    #[serde(default)]
    pub previous_sha256: Option<String>,
    /// Whether the file on disk now differs from what was there (true for a first download). When
    /// it does and a file existed, the old one was kept beside it as
    /// `<name>.<yyyyMMdd-HHmmss-fff>`; see [`Self::caveats`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub changed: bool,
    /// The release now at the path, read back with mzLib's `GeneOntologyGraph.Load`, so a
    /// truncated download fails here rather than in a later [`annotate_go_with`].
    pub go: GoRelease,
    /// The backup kept when a different file was replaced; a file with no `data-version`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One cell of a map-valued column such as `evidence_by_member`: member accession → its strings.
pub type StringListMap = BTreeMap<String, Vec<String>>;

/// A column whose every cell maps a string to a list of strings — `evidence_by_member` — with a
/// wire `null` as `None`.
///
/// # Errors
///
/// [`MzLibError::Usage`] if the column is absent (the message names the columns there are),
/// [`MzLibError::Protocol`] if a cell is not such a map.
pub fn string_list_maps(table: &Table, column: &str) -> Result<Vec<Option<StringListMap>>> {
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
            Value::Object(_) => serde_json::from_value(cell.clone())
                .map(Some)
                .map_err(|_| map_error(column, row, cell)),
            other => Err(map_error(column, row, other)),
        })
        .collect()
}

fn map_error(column: &str, row: usize, cell: &Value) -> MzLibError {
    MzLibError::Protocol(format!(
        "Column '{column}' row {row} is not a map of string lists: {cell}"
    ))
}

fn annotate_go_request(
    groups: &Path,
    database: &Path,
    options: &GoAnnotateOptions,
) -> Result<Vec<String>> {
    let mut args = vec![
        "proteins".to_owned(),
        "annotate-go".to_owned(),
        "--groups".to_owned(),
        optional_path(groups, "groups")?,
        "--database".to_owned(),
        optional_path(database, "database")?,
        "--go-obo".to_owned(),
        optional_path(&options.go_obo, "go_obo").map_err(|_| {
            MzLibError::Usage(
                "go_obo is required: a go.obo you keep, e.g. 'go.obo'. Nothing is downloaded; \
                 fetch a release on purpose with mzlib::proteins::update_go."
                    .to_owned(),
            )
        })?,
    ];
    for (option, value, what) in [
        ("--category-map", &options.category_map, "category_map"),
        ("--out", &options.out, "out"),
        (
            "--categories-out",
            &options.categories_out,
            "categories_out",
        ),
    ] {
        if let Some(path) = value {
            args.push(option.to_owned());
            args.push(optional_path(path, what)?);
        }
    }
    if options.skip_unknown_go_ids {
        args.push("--skip-unknown-go-ids".to_owned());
    }
    if let Some(limit) = options.limit {
        args.push("--limit".to_owned());
        args.push(limit.to_string());
    }
    args.push("--offset".to_owned());
    args.push(options.offset.to_string());
    Ok(args)
}

fn update_go_request(go_obo: &Path) -> Result<Vec<String>> {
    Ok(vec![
        "proteins".to_owned(),
        "update-go".to_owned(),
        "--go-obo".to_owned(),
        optional_path(go_obo, "go_obo")?,
    ])
}

/// Annotate MetaMorpheus protein groups with Gene Ontology terms against a go.obo you keep, with
/// every other default.
///
/// See [`annotate_go_with`] for the reference.
///
/// # Errors
///
/// As [`annotate_go_with`].
pub fn annotate_go(
    groups: impl AsRef<Path>,
    database: impl AsRef<Path>,
    go_obo: impl AsRef<Path>,
) -> Result<GoAnnotations> {
    annotate_go_with(
        groups,
        database,
        &GoAnnotateOptions {
            go_obo: go_obo.as_ref().to_path_buf(),
            ..GoAnnotateOptions::default()
        },
    )
}

/// Annotate a stored MetaMorpheus protein-group table with Gene Ontology terms, keeping every
/// member: one row per (group, term) that any member holds, directly or through an ancestor.
///
/// Wraps mzLib's `GoGroupAnnotator` over the table MetaMorpheus wrote
/// (`ProteinGroupFromTsv.ToGoAnnotationGroups`). For each group it returns one row per GO term that
/// **any** member holds — directly, or by propagation up `is_a` and `part_of` — and each row names
/// the members that carry it. Nothing is collapsed: the union, the consensus and the direct-only
/// views are filters on the rows (see [`GoAnnotations`]).
///
/// **Every non-decoy group gets at least one row.** A group with no term gets a single row whose
/// `annotation_status` says why: `no_go_terms` (its members have none), `no_entry` (a member is not
/// in the database — annotate against the database the search used), or `contaminant`.
///
/// **Pin the ontology.** Terms and their ancestors change between GO releases, so a result means
/// something only relative to one go.obo. This reads the file you name and never downloads one;
/// fetch a release on purpose with [`update_go`], keep the file, and record [`GoAnnotations::go`]
/// (its `sha256`) with your results.
///
/// `groups` is a MetaMorpheus protein-group table — `AllQuantifiedProteinGroups.tsv`, or, from a
/// search without quantification, `AllProteinGroups.tsv` or a file's `<file>_ProteinGroups.tsv`.
/// `database` is the UniProt XML (`.xml` or `.xml.gz`) the search used; a FASTA is refused, because
/// it carries no GO and every group would read `no_go_terms` whatever the proteins are.
///
/// Needs the bridge from pyMzLib 0.3.0 or later: an older one is refused before anything is read.
#[doc = include_str!("../docs/reference/proteins.annotate-go.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for an empty `groups`, `database` or
/// `go_obo` path, or a bridge too old to dispatch the verb.
///
/// # Examples
///
/// A real MetaMorpheus search (PXD036557): five protein groups and a decoy, the five UniProt
/// entries they name, and the GO release mzLib propagated over:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::proteins::{annotate_go_with, string_list_maps, GoAnnotateOptions};
///
/// let go = annotate_go_with(
///     "PXD036557_AllQuantifiedProteinGroups.tsv",
///     "pxd036557_proteins.xml",
///     &GoAnnotateOptions {
///         go_obo: "go-pxd036557.obo".into(),
///         category_map: Some("organelle_map.tsv".into()),
///         ..Default::default()
///     },
/// )?;
/// assert_eq!((go.group_count, go.row_count), (5, 563));
/// assert_eq!(go.go.release.as_deref(), Some("releases/2026-07-26"));
/// assert_eq!(go.header["status_annotated"], "4");     // groups at q <= 0.01, not rows
/// assert_eq!(go.header["status_contaminant"], "1");
///
/// // The histone group: the union holds 106 terms, of which both members carry 57.
/// let groups = go.columns.strings("protein_group")?;
/// let n_with = go.columns.integers("n_with")?;
/// let n_members = go.columns.integers("n_members")?;
/// let histones: Vec<usize> = (0..groups.len())
///     .filter(|&i| groups[i].as_deref() == Some("P0C0S5|Q71UI9"))
///     .collect();
/// let consensus = histones.iter().filter(|&&i| n_with[i] == n_members[i]).count();
/// assert_eq!((histones.len(), consensus), (106, 57));
///
/// // Evidence stays per member.
/// let evidence = string_list_maps(&go.columns, "evidence_by_member")?;
/// assert!(evidence[histones[0]].as_ref().is_some_and(|by| by.contains_key("P0C0S5")));
///
/// let categories = go.categories.as_ref().expect("a category map was given");
/// assert_eq!((categories.map_name.as_str(), categories.row_count), ("organelle", 30));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// For a whole proteome, write the table with mzLib's own writer and take back only the summary:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::proteins::{annotate_go_with, GoAnnotateOptions};
///
/// let go = annotate_go_with(
///     "PXD036557_AllQuantifiedProteinGroups.tsv",
///     "pxd036557_proteins.xml",
///     &GoAnnotateOptions {
///         go_obo: "go-pxd036557.obo".into(),
///         out: Some("go_annotations.tsv".into()),
///         limit: Some(0),
///         ..Default::default()
///     },
/// )?;
/// let written = go.written.as_ref().expect("out was given");
/// assert_eq!(written.row_count, Some(563));            // the whole table, on disk
/// assert_eq!((go.returned_count, go.truncated), (0, true));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/proteins.annotate-go.see-also.md")]
pub fn annotate_go_with(
    groups: impl AsRef<Path>,
    database: impl AsRef<Path>,
    options: &GoAnnotateOptions,
) -> Result<GoAnnotations> {
    let args = annotate_go_request(groups.as_ref(), database.as_ref(), options)?;
    bridge::require_verb("proteins annotate-go", bridge::MZLIB_1_0_593_BRIDGE)?;
    call((args, None), options.timeout)
}

/// Download the current Gene Ontology release (go.obo) to a path, on purpose, keeping any
/// different file already there as a timestamped backup.
///
/// Wraps mzLib's `Loaders.UpdateGeneOntology`. The whole go.obo (tens of megabytes) is streamed from
/// GO's PURL, which always serves the **current** release. When a file is already at the path, it
/// is kept beside the new one as `<name>.<yyyyMMdd-HHmmss-fff>` if the download differs, and left
/// alone if it is the same, so earlier runs stay reproducible. A failed download leaves any
/// existing file untouched.
///
/// This is the only function in this crate that fetches a go.obo. [`annotate_go_with`] never
/// does, so the release a result was computed against is always a file you chose to keep. mzLib
/// itself gives up after two minutes without data; it sets no overall timeout, and neither does
/// this.
///
/// Needs the bridge from pyMzLib 0.3.0 or later: an older one is refused before anything is
/// fetched.
#[doc = include_str!("../docs/reference/proteins.update-go.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for an empty path, or a bridge too old to
/// dispatch the verb.
///
/// # Examples
///
/// The recording was made into an empty folder, so this was a first download:
///
/// ```
/// # mzlib_replay::activate();
/// let update = mzlib::proteins::update_go("go.obo")?;
/// assert_eq!(update.go.release.as_deref(), Some("releases/2026-07-26"));
/// assert_eq!(update.go.term_count, 48340);
/// assert_eq!((update.existed_before, update.changed), (false, true));
/// assert_eq!(update.previous_sha256, None);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/proteins.update-go.see-also.md")]
pub fn update_go(go_obo: impl AsRef<Path>) -> Result<GoUpdate> {
    let args = update_go_request(go_obo.as_ref())?;
    bridge::require_verb("proteins update-go", bridge::MZLIB_1_0_593_BRIDGE)?;
    call((args, None), None)
}

// ---------------------------------------------------------------------------------------------
// Alignment (mzLib 1.0.595, #1428)
// ---------------------------------------------------------------------------------------------

/// The pyMzLib release whose bridge first dispatches `proteins align` and
/// `proteins parse-accession`: their specs' `since.pymzlib`.
const ALIGN_SINCE: &str = "0.6.0";

/// The end-gap rule [`align_with`] uses, as mzLib's aligner id spells it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EndGaps {
    /// A gap at either end of either sequence costs nothing, so a signal peptide or a longer tail
    /// on one side is not penalised (EMBOSS needle's default). mzLib's default.
    Free,
    /// Every gap costs, at the ends too.
    Charged,
}

impl EndGaps {
    /// The wire spelling: `"free"` or `"charged"`.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Free => "free",
            Self::Charged => "charged",
        }
    }
}

/// How [`align_with`] aligns. **Every option left `None` is mzLib's to set**: nothing is sent for
/// it, and the result says what the aligner used ([`Alignments::gap_open`] and the rest).
#[derive(Debug, Clone)]
pub struct AlignOptions {
    /// Cost of opening a gap, in half-bits (BLOSUM62 score units), charged once per gap: a gap of
    /// `k` residues costs `gap_open + k * gap_extend`. `None`: mzLib's, 11 at the pin (BLASTP's
    /// for BLOSUM62). mzLib refuses a negative value.
    pub gap_open: Option<i64>,
    /// Cost of each residue in a gap, the first included, in half-bits (BLOSUM62 score units).
    /// `None`: mzLib's, 1 at the pin. mzLib refuses less than 1.
    pub gap_extend: Option<i64>,
    /// The end-gap rule. `None`: mzLib's, [`EndGaps::Free`].
    pub end_gaps: Option<EndGaps>,
    /// The largest table mzLib will build, `(source length + 1) * (target length + 1)` cells of
    /// one byte each; a larger pair is refused rather than allowed to exhaust memory. `None`:
    /// mzLib's, 200,000,000 cells (about 200 MB) at the pin. Two titin isoforms need about 1.2
    /// billion; mzLib accepts up to 2,147,483,591.
    pub max_cells: Option<u64>,
    /// Return [`Alignments::residue_map`], one row per residue of both sequences. `true` (the
    /// default) sends nothing; `false` sends `--no-residue-map` and returns only the pairs table,
    /// for a large batch where you want the scores and counts.
    pub residue_map: bool,
    /// Time to allow. `None` (the default) waits.
    pub timeout: Option<Duration>,
}

impl Default for AlignOptions {
    fn default() -> Self {
        Self {
            gap_open: None,
            gap_extend: None,
            end_gaps: None,
            max_cells: None,
            residue_map: true,
            timeout: None,
        }
    }
}

/// What [`align_with`] returns: one row per pair, and every residue's counterpart.
///
/// [`Self::columns`] is the pairs table, one row per input pair, in input order: `pair_index`
/// (from 0), `source_length` and `target_length` (residues), `score` (half-bits), `identical` and
/// `substituted` (paired residues that are the same letter, and that differ),
/// `source_residues_facing_gap` and `target_residues_facing_gap` (residues with no counterpart),
/// and `aligned_source` and `aligned_target` (each sequence with `-` where it faces a gap).
/// [`Self::alignments`] gives the same rows as [`Alignment`] structs, each with its residue map.
#[derive(Debug, Clone, Deserialize)]
pub struct Alignments {
    /// mzLib's name for the method and every parameter, e.g.
    /// `"global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"`. One aligner per call. Store it
    /// with any position you map, so the mapping says how it was made.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub aligner_id: String,
    /// The substitution matrix, always `"BLOSUM62"`: mzLib's default and the only one offered.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub matrix: String,
    /// The gap-open cost used, in half-bits (BLOSUM62 score units): yours, or mzLib's default.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub gap_open: i64,
    /// The cost per gap residue used, in half-bits (BLOSUM62 score units): yours, or mzLib's.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub gap_extend: i64,
    /// `"free"` or `"charged"`: the end-gap rule used, yours or mzLib's.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub end_gaps: String,
    /// The table limit used, in cells: yours, or mzLib's default.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub max_cells: u64,
    /// Pairs aligned, in pairs: the rows of [`Self::columns`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// The pairs table, one row per pair; its column order is [`Table::names`].
    #[serde(flatten)]
    pub columns: Table,
    /// One row per residue of each sequence of each pair, pair by pair: a pair's source residues
    /// first, then its target residues. Columns `pair_index`, `sequence` (`"source"` or
    /// `"target"`), `position` (1-based), `residue`, `mapped_position` (1-based position of its
    /// counterpart in the other sequence; null when it faces a gap, **never** the nearest residue)
    /// and `mapped_residue` (the counterpart; null when it faces a gap). `None` when
    /// [`AlignOptions::residue_map`] was `false`.
    #[serde(default)]
    pub residue_map: Option<LongTable>,
    /// What the score and the mapping do and do not mean.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One pair's alignment, from [`align`] or [`Alignments::alignments`].
///
/// Positions are **1-based**, as in mzLib and UniProt. A residue facing a gap has no counterpart,
/// and [`Self::map_source_to_target`] returns `None` for it, never the nearest residue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// Position of the pair in the input, from 0.
    pub pair_index: u64,
    /// mzLib's name for the method and every parameter (the call's [`Alignments::aligner_id`]).
    pub aligner_id: String,
    /// Substitution scores summed over the paired residues, less every charged gap, in BLOSUM62's
    /// half-bits. A raw score: it grows with length, and is not a bit score, an E-value or a
    /// percent identity.
    pub score: i64,
    /// Paired residues that are the same letter (case-insensitive).
    pub identical: u64,
    /// Paired residues that differ.
    pub substituted: u64,
    /// Source residues with no counterpart in the target.
    pub source_residues_facing_gap: u64,
    /// Target residues with no counterpart in the source.
    pub target_residues_facing_gap: u64,
    /// The source's length, in residues.
    pub source_length: usize,
    /// The target's length, in residues.
    pub target_length: usize,
    /// The source with `-` wherever it faces a target residue across a gap.
    pub aligned_source: String,
    /// The target with `-` wherever it faces a source residue across a gap. As long as
    /// [`Self::aligned_source`]; print the two one above the other to see the alignment.
    pub aligned_target: String,
    /// For each source position `p`, item `p - 1` is the paired target position, or `None` when
    /// that residue faces a gap. `None` as a whole when the call asked for no residue map.
    pub source_to_target: Option<Vec<Option<usize>>>,
    /// The same, from each target position to the source.
    pub target_to_source: Option<Vec<Option<usize>>>,
    /// What the score and the mapping do and do not mean (the call's caveats).
    pub caveats: Vec<String>,
}

impl Alignment {
    /// The target position paired with source `position` (1-based), or `None` if it faces a gap.
    ///
    /// It reads the residue map mzLib's `TryMapSourceToTarget` filled; nothing is computed here.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `position` is outside 1 to [`Self::source_length`], or the call
    /// asked for no residue map.
    pub fn map_source_to_target(&self, position: usize) -> Result<Option<usize>> {
        lookup(
            self.source_to_target.as_deref(),
            position,
            self.source_length,
            "source",
        )
    }

    /// The source position paired with target `position` (1-based), or `None` if it faces a gap.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `position` is outside 1 to [`Self::target_length`], or the call
    /// asked for no residue map.
    pub fn map_target_to_source(&self, position: usize) -> Result<Option<usize>> {
        lookup(
            self.target_to_source.as_deref(),
            position,
            self.target_length,
            "target",
        )
    }
}

fn lookup(
    table: Option<&[Option<usize>]>,
    position: usize,
    length: usize,
    what: &str,
) -> Result<Option<usize>> {
    let table = table.ok_or_else(|| {
        MzLibError::Usage(
            "This alignment has no residue map: it came from a call with residue_map: false."
                .to_owned(),
        )
    })?;
    if !(1..=length).contains(&position) {
        return Err(MzLibError::Usage(format!(
            "position {position} is outside the {what}, which has positions 1 to {length} \
             (positions are 1-based, as in mzLib and UniProt)."
        )));
    }
    Ok(table[position - 1])
}

impl Alignments {
    /// One [`Alignment`] per pair, in input order, each with its residue map.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is missing or not the type the wire contract says.
    pub fn alignments(&self) -> Result<Vec<Alignment>> {
        let count = |name: &str| -> Result<Vec<u64>> {
            self.columns
                .integers(name)?
                .into_iter()
                .map(|value| {
                    value.and_then(|v| u64::try_from(v).ok()).ok_or_else(|| {
                        MzLibError::Protocol(format!(
                            "proteins align returned a {name} that is not a count: {value:?}"
                        ))
                    })
                })
                .collect()
        };
        let pair_index = count("pair_index")?;
        let source_length = count("source_length")?;
        let target_length = count("target_length")?;
        let identical = count("identical")?;
        let substituted = count("substituted")?;
        let source_gaps = count("source_residues_facing_gap")?;
        let target_gaps = count("target_residues_facing_gap")?;
        let score = self.columns.integers("score")?;
        let aligned_source = self.columns.strings("aligned_source")?;
        let aligned_target = self.columns.strings("aligned_target")?;

        // (pair, sequence) -> position -> mapped position, from the residue map.
        let mut maps: BTreeMap<(u64, String), BTreeMap<usize, Option<usize>>> = BTreeMap::new();
        if let Some(map) = &self.residue_map {
            let pairs = map.columns.integers("pair_index")?;
            let sequences = map.columns.strings("sequence")?;
            let positions = map.columns.integers("position")?;
            let mapped = map.columns.integers("mapped_position")?;
            for i in 0..pairs.len() {
                let (Some(pair), Some(sequence), Some(position)) =
                    (pairs[i], sequences[i].clone(), positions[i])
                else {
                    return Err(MzLibError::Protocol(format!(
                        "proteins align returned a residue-map row with no pair, sequence or \
                         position (row {i})"
                    )));
                };
                maps.entry((pair as u64, sequence))
                    .or_default()
                    .insert(position as usize, mapped[i].map(|m| m as usize));
            }
        }
        let table = |pair: u64, sequence: &str, length: usize| -> Option<Vec<Option<usize>>> {
            self.residue_map.as_ref()?;
            let found = maps.get(&(pair, sequence.to_owned()));
            Some(
                (1..=length)
                    .map(|p| found.and_then(|f| f.get(&p).copied().flatten()))
                    .collect(),
            )
        };

        (0..pair_index.len())
            .map(|i| {
                let pair = pair_index[i];
                let (source, target) = (source_length[i] as usize, target_length[i] as usize);
                Ok(Alignment {
                    pair_index: pair,
                    aligner_id: self.aligner_id.clone(),
                    score: score[i].ok_or_else(|| {
                        MzLibError::Protocol(format!("proteins align returned no score (row {i})"))
                    })?,
                    identical: identical[i],
                    substituted: substituted[i],
                    source_residues_facing_gap: source_gaps[i],
                    target_residues_facing_gap: target_gaps[i],
                    source_length: source,
                    target_length: target,
                    aligned_source: aligned_source[i].clone().unwrap_or_default(),
                    aligned_target: aligned_target[i].clone().unwrap_or_default(),
                    source_to_target: table(pair, "source", source),
                    target_to_source: table(pair, "target", target),
                    caveats: self.caveats.clone(),
                })
            })
            .collect()
    }
}

/// The pairs as stdin lines, `source<TAB>target`, unchanged, refusing what cannot travel.
fn pair_lines<S: AsRef<str>, T: AsRef<str>>(pairs: &[(S, T)]) -> Result<Vec<String>> {
    if pairs.is_empty() {
        return Err(MzLibError::Usage(
            "At least one pair is required, e.g. [(\"MKTAYIAK\", \"MKTSYIAK\")].".to_owned(),
        ));
    }
    pairs
        .iter()
        .enumerate()
        .map(|(index, (source, target))| {
            let (source, target) = (source.as_ref(), target.as_ref());
            for (name, value) in [("source", source), ("target", target)] {
                // One pair is one stdin line, its two sequences split on the tab, so neither may
                // hold a tab or a line break: either would shift or split every later pair.
                if value.contains(['\t', '\n', '\r']) {
                    return Err(MzLibError::Usage(format!(
                        "Pair {index}'s {name} contains a tab or line break: {value:?}."
                    )));
                }
            }
            // The bridge skips a blank line, which would put every later result on the wrong pair.
            if format!("{source}{target}").trim().is_empty() {
                return Err(MzLibError::Usage(format!(
                    "Pair {index} has two empty sequences; there is nothing to align."
                )));
            }
            Ok(format!("{source}\t{target}"))
        })
        .collect()
}

/// The arguments and stdin of a `proteins align` call: only the options that were given.
fn align_request<S: AsRef<str>, T: AsRef<str>>(
    pairs: &[(S, T)],
    options: &AlignOptions,
) -> Result<(Vec<String>, Option<String>)> {
    let lines = pair_lines(pairs)?;
    let mut args = vec!["proteins".to_owned(), "align".to_owned()];
    for (name, value) in [
        ("--gap-open", options.gap_open.map(|v| v.to_string())),
        ("--gap-extend", options.gap_extend.map(|v| v.to_string())),
        ("--max-cells", options.max_cells.map(|v| v.to_string())),
    ] {
        if let Some(value) = value {
            args.push(name.to_owned());
            args.push(value);
        }
    }
    if let Some(end_gaps) = options.end_gaps {
        args.push("--end-gaps".to_owned());
        args.push(end_gaps.as_str().to_owned());
    }
    if !options.residue_map {
        args.push("--no-residue-map".to_owned());
    }
    Ok((args, Some(lines.join("\n") + "\n")))
}

/// Align two protein sequences with mzLib, with every default, and map residues from one onto
/// the other.
///
/// The one-pair form of [`align_with`]: the same bridge call and the same mzLib defaults,
/// returning one [`Alignment`] with its residue map. To set an option for one pair, pass
/// `&[(source, target)]` to [`align_with`].
///
/// ```
/// # mzlib_replay::activate();
/// // mzLib's own deletion case: QRQ (source 9 to 11) is missing from the target.
/// let aln = mzlib::proteins::align("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")?;
/// assert_eq!(aln.aligned_source, "MKTAYIAKQRQISFVKSHFSRQ");
/// assert_eq!(aln.aligned_target, "MKTAYIAK---ISFVKSHFSRQ");
/// assert_eq!((aln.score, aln.identical, aln.source_residues_facing_gap), (80, 19, 3));
/// let mapped: Vec<Option<usize>> = [8, 9, 10, 11, 12]
///     .iter()
///     .map(|&p| aln.map_source_to_target(p))
///     .collect::<Result<_, _>>()?;
/// assert_eq!(mapped, [Some(8), None, None, None, Some(9)]);
/// assert_eq!(aln.map_target_to_source(9)?, Some(12));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors
///
/// As [`align_with`].
pub fn align(source: &str, target: &str) -> Result<Alignment> {
    let result = align_with(&[(source, target)], &AlignOptions::default())?;
    result.alignments()?.into_iter().next().ok_or_else(|| {
        MzLibError::Protocol("proteins align returned no row for the pair it was given".to_owned())
    })
}

/// Align pairs of protein sequences in one call, and map every residue across each pair.
///
/// Wraps mzLib's `PairwiseAligner.Align` (mzLib #1428): a **global** alignment with affine gaps
/// (Gotoh's form of Needleman-Wunsch) scored with **BLOSUM62**, and the residue map its
/// `TryMapSourceToTarget` and `TryMapTargetToSource` give. Use it to carry a position (a
/// modification site, a variant, a cleavage site) from one sequence to another: an isoform, an
/// ortholog, an older release of the same entry. Pairs are aligned in one bridge process, in
/// order, so a few thousand isoform pairs cost one process start.
///
/// Each sequence is sent exactly as given: case is kept (and compared case-insensitively), and
/// either sequence may be empty, not both. mzLib refuses `-` and `.`; it scores any other letter
/// BLOSUM62 does not name (`U`, `O`, a digit, a space) as `X`, so strip modifications first. The
/// [module guide](crate::proteins#map-a-residue-from-one-sequence-onto-another) carries a phosphosite
/// across isoforms and an ortholog.
#[doc = include_str!("../docs/reference/proteins.align.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for no pairs, two empty sequences, or a
/// sequence holding a tab or line break; and when the bridge in use is older than pyMzLib 0.6.0's,
/// which has no `proteins align`.
///
/// # Examples
///
/// A deletion and an extra N-terminus, mzLib's own test cases, with mzLib's defaults:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::proteins::{align_with, AlignOptions};
///
/// let alns = align_with(
///     &[("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ"), ("MLLLPEPTIDE", "PEPTIDE")],
///     &AlignOptions::default(),
/// )?;
/// assert_eq!(alns.aligner_id, "global-affine;BLOSUM62;open=11;extend=1;end-gaps=free");
/// assert_eq!(alns.columns.integers("score")?, [Some(80), Some(39)]);
/// assert_eq!(alns.columns.integers("source_residues_facing_gap")?, [Some(3), Some(4)]);
/// let first_residues: Vec<Option<usize>> = alns
///     .alignments()?
///     .iter()
///     .map(|a| a.map_target_to_source(1))
///     .collect::<Result<_, _>>()?;
/// assert_eq!(first_residues, [Some(1), Some(5)]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// Charging the end gaps changes the score of an overhang, and the aligner id says so:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::proteins::{align_with, AlignOptions, EndGaps};
///
/// let charged = align_with(
///     &[("MLLLPEPTIDE", "PEPTIDE")],
///     &AlignOptions { end_gaps: Some(EndGaps::Charged), ..Default::default() },
/// )?;
/// let only = &charged.alignments()?[0];
/// assert_eq!((only.aligned_target.as_str(), only.score), ("----PEPTIDE", 24));
/// assert_eq!(charged.aligner_id, "global-affine;BLOSUM62;open=11;extend=1;end-gaps=charged");
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/proteins.align.see-also.md")]
pub fn align_with<S: AsRef<str>, T: AsRef<str>>(
    pairs: &[(S, T)],
    options: &AlignOptions,
) -> Result<Alignments> {
    let request = align_request(pairs, options)?;
    bridge::require_verb("proteins align", ALIGN_SINCE)?;
    call(request, options.timeout)
}

// ---------------------------------------------------------------------------------------------
// Accession parsing (mzLib 1.0.595, #1382)
// ---------------------------------------------------------------------------------------------

/// Every `namespace` [`parse_accessions`] can report: mzLib's `AccessionNamespace`, lower case,
/// exactly as [`resolve_genes_with`] writes it.
pub const NAMESPACES: [&str; 3] = ["uniprot", "refseq", "unrecognized"];

/// What [`parse_accessions`] returns: one row per accession, in the order given.
///
/// [`Self::columns`] has `input` (exactly as given), `entry_accession` (the UniProt entry or
/// unversioned RefSeq accession, `P04637` for `P04637-2_S70N`; `input` when unrecognized: an
/// **entry**, not a sequence), `isoform` (UniProt `-n`; null when there is no suffix, always for
/// RefSeq, or a number past 2,147,483,647), `version` (RefSeq `.n`; null likewise),
/// `namespace` (one of [`NAMESPACES`]), `parent_accession` (the entry the variants were applied
/// to, isoform or version kept; `input` when there are no applied variants), `applied_variants`
/// (the variant suffix as written, `S70N_A80T`; null when none) and `has_applied_variants`.
#[derive(Debug, Clone, Deserialize)]
pub struct ParsedAccessions {
    /// Accessions parsed, in accessions: the rows, duplicates included.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Namespace to the accessions in it. Every key in [`NAMESPACES`] is present, zeros included.
    /// A large `"unrecognized"` count is worth reading: see [`Self::caveats`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub namespace_counts: BTreeMap<String, u64>,
    /// Accessions that name applied variants, in accessions.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub with_applied_variants_count: u64,
    /// The table, one row per accession; its column order is [`Table::names`].
    #[serde(flatten)]
    pub columns: Table,
    /// What mzLib leaves unrecognized rather than repairs, and other traps.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl ParsedAccessions {
    /// Accession to entry accession, to group isoforms and variants under their entry. An
    /// accession given twice appears once (it parses the same); an unrecognized one maps to
    /// itself.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is missing or not the type the wire contract says.
    pub fn entry_of(&self) -> Result<BTreeMap<String, String>> {
        let inputs = self.columns.strings("input")?;
        let entries = self.columns.strings("entry_accession")?;
        Ok(inputs
            .into_iter()
            .zip(entries)
            .filter_map(|(input, entry)| Some((input?, entry?)))
            .collect())
    }
}

/// The accessions as stdin lines, unchanged: mzLib parses and never repairs, so neither does this
/// crate. A blank entry or a line break would shift every later row, so both are refused.
fn accession_lines<S: AsRef<str>>(accessions: &[S]) -> Result<Vec<String>> {
    if accessions.is_empty() {
        return Err(MzLibError::Usage(
            "At least one accession is required, e.g. [\"P04637-2\"].".to_owned(),
        ));
    }
    accessions
        .iter()
        .map(|value| {
            let value = value.as_ref();
            // A space around an accession is NOT stripped: mzLib reads it as unrecognized, and
            // that is the answer.
            if value.trim().is_empty() {
                return Err(MzLibError::Usage(
                    "accessions contains a blank entry; every entry is one result row.".to_owned(),
                ));
            }
            if value.contains(['\n', '\r']) {
                return Err(MzLibError::Usage(format!(
                    "An accession contains a line break: {value:?}."
                )));
            }
            Ok(value.to_owned())
        })
        .collect()
}

/// Read accessions into entry, isoform or version, and applied variants, with mzLib.
///
/// Wraps mzLib's `VariantApplication.ParseAccession` (mzLib #1382), and through it
/// `ProteinAccession.Parse`. It is the inverse of the names mzLib gives variant proteoforms:
/// `P04637_S70N` is entry `P04637` with `S70N` applied, and `NP_000537.3_R72P` is RefSeq
/// `NP_000537`, version 3, with `R72P`. Use it to link a search result's accession to anything you
/// keep per entry (genes, orthologs, GO terms) without losing the isoform or the variants.
///
/// **It parses and never repairs.** An accession outside UniProt's and RefSeq's protein grammars
/// comes back `"unrecognized"` and verbatim: lower case, a stray space, a FASTA header, a `DECOY_`
/// prefix, a `|`-joined protein group, and `P12345_2`, mzLib's loader counter for a *different*
/// entry whose accession collided. Each accession is sent exactly as given, one stdin line each,
/// duplicates included. The [module
/// guide](crate::proteins#read-an-accession-into-entry-isoform-and-variants) reads every case.
#[doc = include_str!("../docs/reference/proteins.parse-accession.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for no accessions, or one that is blank or
/// holds a line break (either would shift every later row); and when the bridge in use is older
/// than pyMzLib 0.6.0's, which has no `proteins parse-accession`.
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let parsed = mzlib::proteins::parse_accessions(&[
///     "P04637", "P04637-2", "Q13409-2", "NP_000537.3_R72P", "P04637_S70N_A80T",
///     "P04637_2", "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768",
/// ])?;
/// let entries = parsed.columns.strings("entry_accession")?;
/// assert_eq!(entries[3].as_deref(), Some("NP_000537"));
/// assert_eq!(parsed.columns.integers("isoform")?[..3], [None, Some(2), Some(2)]);
/// assert_eq!(parsed.columns.integers("version")?[3], Some(3));
/// let variants = parsed.columns.strings("applied_variants")?;
/// assert_eq!((variants[3].as_deref(), variants[4].as_deref()), (Some("R72P"), Some("S70N_A80T")));
/// assert_eq!(parsed.namespace_counts["unrecognized"], 4);
/// assert_eq!(parsed.entry_of()?["P04637_S70N_A80T"], "P04637");
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/proteins.parse-accession.see-also.md")]
pub fn parse_accessions<S: AsRef<str>>(accessions: &[S]) -> Result<ParsedAccessions> {
    let lines = accession_lines(accessions)?;
    bridge::require_verb("proteins parse-accession", ALIGN_SINCE)?;
    call(
        (
            vec!["proteins".to_owned(), "parse-accession".to_owned()],
            Some(lines.join("\n") + "\n"),
        ),
        None,
    )
}

// ---------------------------------------------------------------------------------------------
// Orthology (mzLib 1.0.595, #1381)
// ---------------------------------------------------------------------------------------------

/// The pyMzLib release whose bridge first dispatches `genes orthology`.
const ORTHOLOGY_SINCE: &str = "0.6.0";

/// Every status [`build_orthology_with`] can report in [`OrthologySnapshot::pair_status`], as mzLib
/// writes it (`OrthologySnapshot.StatusName`), in mzLib's order. Every gene gets exactly one, so
/// "Compara found nothing", "there was nothing to compare against" and "no call was attempted"
/// never share an empty result. `not_in_gene_set` is for a gene id from elsewhere, and never occurs
/// in a pair-status table, whose rows are the gene set itself.
pub const ORTHOLOGY_STATUSES: [&str; 5] = [
    "has_ortholog",
    "no_edge_in_shared_tree",
    "tree_lacks_target_species",
    "not_in_any_tree",
    "not_in_gene_set",
];

/// What [`build_orthology_with`] reads, and what it writes. Every default is mzLib's.
///
/// `gene_trees`, `dumps` and at least one species in `gtfs` or `gene_sets` are required; the
/// default leaves them empty, which is refused before anything is spawned.
#[derive(Debug, Clone, Default)]
pub struct OrthologyOptions {
    /// Each species, as its Ensembl production name (`"homo_sapiens"`, as the dumps' `species`
    /// column spells it), with its Ensembl primary-assembly GTF, plain or `.gz`. These are the
    /// snapshot's species, sent on stdin one line each, in this order.
    pub gtfs: Vec<(String, PathBuf)>,
    /// Instead of a GTF, a compact gene table mzLib's `EnsemblGeneSetWriter` made from one (it keeps
    /// the GTF's provenance, so the snapshot is the same). A species may be in `gtfs` or
    /// `gene_sets`, not both.
    pub gene_sets: Vec<(String, PathBuf)>,
    /// One Compara homology dump per species
    /// (`homologies/<species>/Compara.<release>.protein_default.homologies.tsv.gz`), in any order:
    /// each dump's species is read from its rows. **Every species' dump is required**; mzLib
    /// refuses a call without one.
    pub dumps: Vec<PathBuf>,
    /// Compara's gene-tree file, `<collection>.GeneTree_content.default.e<release>.txt.gz`.
    /// Required.
    pub gene_trees: PathBuf,
    /// Keep only the gene-tree members in your gene sets. `false` (the default) is mzLib's default
    /// and keeps every gene in the file: 4.2 million for the vertebrates of release 116. The
    /// snapshot is the same either way; only memory differs.
    pub restrict_gene_trees: bool,
    /// A directory to write the snapshot to: Parquet tables (`genes/`, `members/`, `pairs/`),
    /// `views.sql` and `manifest.json`. It must not exist yet or be empty: mzLib writes a snapshot
    /// once and never updates one in place. `None` (the default): nothing is written.
    pub out: Option<PathBuf>,
    /// `(from_species, to_species)`: return every gene of `from_species` with its status in
    /// `to_species` ([`PairStatus`]). Two different species of the snapshot, in that order: a
    /// status is not symmetric. `None` (the default): not computed.
    pub pair_status: Option<(String, String)>,
    /// Time to allow. `None` (the default) waits: for human, mouse and rat at release 116, mzLib
    /// reports about 20 s to build and 40 s more to write.
    pub timeout: Option<Duration>,
}

/// Compara's gene-tree file a snapshot was built with, and what was kept from it.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct GeneTrees {
    /// The file's name, e.g. `vertebrates.GeneTree_content.default.e116.txt.gz`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub source_file_name: String,
    /// Lower-case hex sha256 of the file's bytes as read (the compressed bytes for a `.gz`).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sha256: String,
    /// The release in the file name (`.e116.txt`); `None` when the name carries none, and then
    /// mzLib could not hold the file to the release you gave (a caveat says so).
    #[serde(default)]
    pub release: Option<String>,
    /// The Compara collection in the file name, e.g. `"vertebrates"`; `None` when the name does
    /// not start `<collection>.GeneTree_content.`.
    #[serde(default)]
    pub collection: Option<String>,
    /// Genes held in memory from the file, in genes: every gene in it, or with
    /// [`OrthologyOptions::restrict_gene_trees`] only those in your gene sets.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub member_count: u64,
    /// The sha256s of the gene sets the members were restricted to; `None` when every gene in the
    /// file was kept.
    #[serde(default)]
    pub restricted_to_gene_sets: Option<Vec<String>>,
}

/// One species' Compara homology dump, as read.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct HomologyDump {
    /// The file's name. Ensembl names every species' dump
    /// `Compara.<release>.protein_default.homologies.tsv.gz`, so this alone does not say whose it
    /// is: the key of [`OrthologySnapshot::dumps`] does, read from the dump's rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub source_file_name: String,
    /// Lower-case hex sha256 of the file's bytes as read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sha256: String,
    /// The release in the file name; `None` when it carries none (then not checked).
    #[serde(default)]
    pub release: Option<String>,
    /// Every data row read and checked, kept or not, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub rows_read: u64,
    /// Rows whose two species are both snapshot species, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub rows_kept: u64,
}

/// How many relationships one species pair holds: one row each in `pairs/<a>__<b>.parquet`.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct SpeciesPair {
    /// The species that sorts first: side A of every row of the pair.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub species_a: String,
    /// The other; the same species for its paralogs, the pair `(s, s)`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub species_b: String,
    /// Relationships between the two, in rows. One row each, never merged or picked.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Of those, orthologs (`ortholog_one2one`, `ortholog_one2many`, `ortholog_many2many`), in rows.
    /// Always 0 for `(s, s)`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub ortholog_count: u64,
    /// Paralogs (`within_species_paralog`, `other_paralog`, `gene_split`), in rows. Always 0
    /// between two species.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub paralog_count: u64,
    /// Homoeologs, between subgenomes of one polyploid species, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub homoeolog_count: u64,
}

/// One file of a written snapshot.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct SnapshotFile {
    /// Relative to [`WrittenSnapshot::directory`], with forward slashes, e.g.
    /// `"pairs/homo_sapiens__mus_musculus.parquet"`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Rows in the Parquet table, in rows; `None` for a file that is not a table (`views.sql`). A
    /// pair with no relationship is still written, with 0 rows.
    #[serde(default)]
    pub rows: Option<u64>,
    /// The file's size in bytes.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub bytes: u64,
    /// Lower-case hex sha256 of the file, as `manifest.json` lists it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sha256: String,
}

/// Where [`build_orthology_with`] wrote the snapshot, with [`OrthologyOptions::out`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct WrittenSnapshot {
    /// The `out` directory, absolute.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub directory: String,
    /// Every file written but `manifest.json`, in order of path: the Parquet tables
    /// `genes/<species>.parquet`, `members/<species>.parquet` and `pairs/<a>__<b>.parquet`, and
    /// `views.sql`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<SnapshotFile>,
}

impl WrittenSnapshot {
    /// The path of `manifest.json`, which lists every input and file with its sha256.
    #[must_use]
    pub fn manifest(&self) -> PathBuf {
        Path::new(&self.directory).join("manifest.json")
    }
}

/// Every gene of one species with how far Compara's inference got in another.
///
/// [`Self::columns`] has one row per gene of [`Self::from_species`]'s gene set, every biotype, in
/// order of gene id: `gene_id` (stable Ensembl id, unversioned), `gene_name` (the release's gene
/// name, a label; null when the GTF gives none), `gene_biotype` (`unknown` when the GTF omits it)
/// and `status` (one of [`ORTHOLOGY_STATUSES`] but `not_in_gene_set`).
#[derive(Debug, Clone, Deserialize)]
pub struct PairStatus {
    /// Whose genes the rows are.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub from_species: String,
    /// The species each status is about.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub to_species: String,
    /// Genes: every gene of `from_species`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Status to the genes with it. Every key in [`ORTHOLOGY_STATUSES`] is present, zeros
    /// included; `not_in_gene_set` is always 0 here.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub status_counts: BTreeMap<String, u64>,
    /// The table, one row per gene; its column order is [`Table::names`].
    #[serde(flatten)]
    pub columns: Table,
}

impl PairStatus {
    /// The gene ids with one status, in order of gene id.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `status` is not one of [`ORTHOLOGY_STATUSES`];
    /// [`MzLibError::Protocol`] if a column is missing or not a string column.
    pub fn genes(&self, status: &str) -> Result<Vec<String>> {
        if !ORTHOLOGY_STATUSES.contains(&status) {
            return Err(MzLibError::Usage(format!(
                "status must be one of {}; got {status:?}.",
                ORTHOLOGY_STATUSES.join(", ")
            )));
        }
        let ids = self.columns.strings("gene_id")?;
        let statuses = self.columns.strings("status")?;
        Ok(ids
            .into_iter()
            .zip(statuses)
            .filter(|(_, s)| s.as_deref() == Some(status))
            .filter_map(|(id, _)| id)
            .collect())
    }
}

/// What [`build_orthology_with`] returns: an Ensembl Compara orthology snapshot, as built.
///
/// The relationships themselves are in the Parquet files [`OrthologyOptions::out`] writes, not
/// here: this says what went in, what came out, and, with [`OrthologyOptions::pair_status`],
/// every gene's status in one other species.
#[derive(Debug, Clone, Deserialize)]
pub struct OrthologySnapshot {
    /// The Ensembl release you gave.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub release: String,
    /// The snapshot's species, Ensembl production names, in sorted order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub species: Vec<String>,
    /// `"ensembl-orthology-snapshot"`, as `manifest.json` names the format.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub format: String,
    /// The format's version.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub format_version: u32,
    /// mzLib's id for the snapshot: a sha256 over the format, the release, the species and every
    /// input's sha256. The same inputs give the same id, written or not, and it is
    /// `manifest.json`'s `snapshot_id`. Record it with any result drawn from the snapshot.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub snapshot_id: String,
    /// Homology rows found in both of a pair's dumps, identical, and counted once, in rows. 0 in
    /// release 116, whose dumps partition.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub identical_duplicates_dropped: u64,
    /// Species to the [`GeneSet`] it was counted against, the same type [`resolve_genes_with`]
    /// reports. Its `sha256` is what `manifest.json`'s `gene_set_sha256` and
    /// [`resolve_genes_with`]'s `gene_set_sha256` column carry, so a resolution and a snapshot can
    /// be checked against one gene set.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub gene_sets: BTreeMap<String, GeneSet>,
    /// The gene-tree file.
    pub gene_trees: GeneTrees,
    /// Species to its [`HomologyDump`]. The species is read from the dump's rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub dumps: BTreeMap<String, HomologyDump>,
    /// One [`SpeciesPair`] per species pair `(a, b)`, `a <= b`, `(s, s)` included.
    /// [`Self::pair`] finds one in either order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub pairs: Vec<SpeciesPair>,
    /// Where the snapshot was written; `None` without [`OrthologyOptions::out`], when nothing was
    /// written.
    #[serde(default)]
    pub written: Option<WrittenSnapshot>,
    /// Every gene's status in one other species; `None` without [`OrthologyOptions::pair_status`].
    #[serde(default)]
    pub pair_status: Option<PairStatus>,
    /// This call's traps: a file name with no release (not checked against it), dumps sharing one
    /// file name, a GTF named for another release, duplicate rows, unrestricted gene trees.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl OrthologySnapshot {
    /// The counts for species `a` and `b`, given in either order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if `a` or `b` is not one of [`Self::species`].
    pub fn pair(&self, a: &str, b: &str) -> Result<&SpeciesPair> {
        let (first, second) = if a <= b { (a, b) } else { (b, a) };
        self.pairs
            .iter()
            .find(|pair| pair.species_a == first && pair.species_b == second)
            .ok_or_else(|| {
                MzLibError::Usage(format!(
                    "({a:?}, {b:?}) is not a pair of this snapshot's species: {}.",
                    self.species.join(", ")
                ))
            })
    }
}

/// A species name for stdin or `--pair-status`: no tab, line break or comma.
fn species_name(name: &str, what: &str) -> Result<String> {
    let name = name.trim();
    if name.is_empty() {
        return Err(MzLibError::Usage(format!(
            "Every species in {what} must be a non-empty Ensembl production name."
        )));
    }
    if name.contains(['\t', '\r', '\n', ',']) {
        return Err(MzLibError::Usage(format!(
            "A species name cannot hold a tab, line break or comma: {name:?}."
        )));
    }
    Ok(name.to_owned())
}

/// The stdin species lines, `species<TAB>gtf|gene-set<TAB>path`, in the order given: the GTFs,
/// then the gene sets. `genes orthology` and `genes map-residues` read them the same way.
fn species_lines(
    gtfs: &[(String, PathBuf)],
    gene_sets: &[(String, PathBuf)],
) -> Result<Vec<String>> {
    let mut lines = Vec::new();
    let mut seen: Vec<String> = Vec::new();
    for (kind, list, what) in [("gtf", gtfs, "gtfs"), ("gene-set", gene_sets, "gene_sets")] {
        for (species, path) in list {
            let name = species_name(species, what)?;
            if seen.contains(&name) {
                return Err(MzLibError::Usage(format!(
                    "'{name}' is given twice in gtfs and gene_sets; give each species one gene set."
                )));
            }
            let path = optional_path(path, &format!("{what} {name}"))?;
            if path.contains(['\t', '\n', '\r']) {
                return Err(MzLibError::Usage(format!(
                    "The path for {name} contains a tab or line break: {path:?}."
                )));
            }
            lines.push(format!("{name}\t{kind}\t{path}"));
            seen.push(name);
        }
    }
    Ok(lines)
}

/// The arguments and stdin of a `genes orthology` call: the species lines
/// (`species<TAB>gtf|gene-set<TAB>path`), a `--` line, then one dump per line.
fn orthology_request(
    release: &str,
    options: &OrthologyOptions,
) -> Result<(Vec<String>, Option<String>)> {
    let release = release.trim();
    if release.is_empty() || release.contains(char::is_whitespace) {
        return Err(MzLibError::Usage(format!(
            "release must be an Ensembl release, e.g. \"116\"; got {release:?}."
        )));
    }
    let mut lines = species_lines(&options.gtfs, &options.gene_sets)?;
    if lines.is_empty() {
        return Err(MzLibError::Usage(
            "At least one species is required: gtfs: vec![(\"homo_sapiens\".into(), \
             \"Homo_sapiens.GRCh38.116.gtf.gz\".into())] (or gene_sets for compact gene tables)."
                .to_owned(),
        ));
    }
    if options.dumps.is_empty() {
        return Err(MzLibError::Usage(
            "dumps is empty. Give every species' Compara homology dump: each pair's \
             relationships are split between the pair's two dumps, so mzLib requires them all."
                .to_owned(),
        ));
    }
    let dumps = bridge::path_lines(&options.dumps, "homology dump")?;
    if let Some(bad) = dumps.iter().find(|line| line.as_str() == SECTION) {
        return Err(MzLibError::Usage(format!(
            "{bad:?} is not a usable homology dump path: it is the line that separates the species \
             from the dumps."
        )));
    }

    if options.gene_trees.as_os_str().is_empty() {
        return Err(MzLibError::Usage(
            "gene_trees is required: Compara's \
             <collection>.GeneTree_content.default.e<release>.txt.gz."
                .to_owned(),
        ));
    }
    let mut args = vec![
        "genes".to_owned(),
        "orthology".to_owned(),
        "--release".to_owned(),
        release.to_owned(),
        "--gene-trees".to_owned(),
        optional_path(&options.gene_trees, "gene_trees")?,
    ];
    if options.restrict_gene_trees {
        args.push("--restrict-gene-trees".to_owned());
    }
    if let Some(out) = &options.out {
        args.push("--out".to_owned());
        args.push(optional_path(out, "out")?);
    }
    if let Some((from, to)) = &options.pair_status {
        let from = species_name(from, "pair_status")?;
        let to = species_name(to, "pair_status")?;
        args.push("--pair-status".to_owned());
        args.push(format!("{from},{to}"));
    }

    lines.push(SECTION.to_owned());
    lines.extend(dumps);
    Ok((args, Some(lines.join("\n") + "\n")))
}

/// Build an Ensembl Compara orthology snapshot for your species, and write it as Parquet.
///
/// Wraps mzLib's `OrthologySnapshot.Build` and `OrthologySnapshotWriter.WriteAsync` (mzLib
/// #1381). Which genes of one species are orthologous to which genes of another, kept exactly as
/// Compara states it: one-to-many and many-to-many orthologies are several rows, paralogs are kept,
/// and nothing is picked, merged or collapsed. Local files only; nothing is downloaded.
///
/// **Every species' homology dump is required, and mzLib refuses a call without one.** Compara
/// puts each of a pair's relationships in one of the pair's two dumps, arbitrarily, and in practice
/// they partition: in release 116 every human-mouse orthology is in the *mouse* dump and none is in
/// the human one. Reading the dump named after your species gives a complete-looking, badly
/// incomplete answer.
///
/// **Orthology is pairwise and not transitive.** A~B and B~C do not make A~C, so never chain
/// pairs into groups. The written `views.sql` has a `species_set3` macro that forms a tuple only
/// when every two members are orthologs.
///
/// The relationships are in the Parquet files [`OrthologyOptions::out`] writes, which mzLib cannot
/// read back and this crate does not read (it carries no Parquet reader): query them with DuckDB,
/// as the [module guide](crate::proteins#orthologs-across-species-an-ensembl-compara-snapshot)
/// shows. Without `out`, the call checks that the inputs fit together and, with
/// [`OrthologyOptions::pair_status`], says gene by gene how far Compara's inference got in one
/// other species.
///
/// `release` is the Ensembl release of the Compara inputs, e.g. `"116"`. Required: mzLib holds the
/// gene trees and every dump to it when their file names carry a release. A GTF's own number is not
/// compared with it (Ensembl publishes some species with their own annotation release).
#[doc = include_str!("../docs/reference/genes.orthology.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for a blank release or one holding a space;
/// no species, a species given twice or a name holding a tab, line break or comma; no dumps, or a
/// blank one; no gene trees; a blank `out`; and when the bridge in use is older than pyMzLib
/// 0.6.0's, which has no `genes orthology`.
///
/// # Examples
///
/// mzLib's own test world: four human genes, three mouse and two rat, with every human-mouse row
/// in the mouse dump, as in release 116.
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::proteins::{build_orthology_with, OrthologyOptions};
///
/// let species = ["homo_sapiens", "mus_musculus", "rattus_norvegicus"];
/// let gtfs = vec![
///     ("homo_sapiens".to_owned(), "Homo_sapiens.GRCh38.116.gtf".into()),
///     ("mus_musculus".to_owned(), "Mus_musculus.GRCm39.116.gtf".into()),
///     ("rattus_norvegicus".to_owned(), "Rattus_norvegicus.GRCr8.116.gtf".into()),
/// ];
/// let dumps = species
///     .iter()
///     .map(|s| format!("homologies/{s}/Compara.116.protein_default.homologies.tsv.gz").into())
///     .collect();
/// let snap = build_orthology_with(
///     "116",
///     &OrthologyOptions {
///         gtfs,
///         dumps,
///         gene_trees: "vertebrates.GeneTree_content.default.e116.txt.gz".into(),
///         pair_status: Some(("homo_sapiens".into(), "mus_musculus".into())),
///         ..Default::default()
///     },
/// )?;
/// assert_eq!(snap.dumps["mus_musculus"].rows_kept, 5);
/// assert_eq!(snap.pair("mus_musculus", "homo_sapiens")?.row_count, 3);   // either order
/// let status = snap.pair_status.as_ref().expect("pair_status was asked for");
/// assert_eq!(
///     status.columns.strings("status")?,
///     ["has_ortholog", "has_ortholog", "tree_lacks_target_species", "not_in_any_tree"]
///         .map(|s| Some(s.to_owned()))
/// );
/// assert_eq!(status.genes("not_in_any_tree")?, ["ENSG00000000004"]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/genes.orthology.see-also.md")]
pub fn build_orthology_with(
    release: &str,
    options: &OrthologyOptions,
) -> Result<OrthologySnapshot> {
    let request = orthology_request(release, options)?;
    bridge::require_verb("genes orthology", ORTHOLOGY_SINCE)?;
    call(request, options.timeout)
}

// ---------------------------------------------------------------------------------------------
// Residues across species: Compara's gene-tree alignment (mzLib 1.0.596, #1465)
// ---------------------------------------------------------------------------------------------

/// The pyMzLib release whose bridge first dispatches `genes map-residues`.
const MAP_RESIDUES_SINCE: &str = "0.7.0";

/// Every outcome [`map_residues_with`] can report: mzLib's `ComparaColumnOutcome` names, verbatim
/// and in mzLib's order. Every query gets exactly one: a residue (`Aligned`), a gap in the target's
/// row (`GapInTarget`), proteins in different gene trees (`DifferentAlignments`), or a protein the
/// alignment does not hold (`SourceNotInAlignment`, `TargetNotInAlignment`).
pub const COLUMN_OUTCOMES: [&str; 5] = [
    "Aligned",
    "GapInTarget",
    "DifferentAlignments",
    "SourceNotInAlignment",
    "TargetNotInAlignment",
];

/// What [`map_residues_with`] reads. Every default is mzLib's.
///
/// `alignment` and `gene_trees` are required; the default leaves them empty, which is refused
/// before anything is spawned. Give your species' gene sets in `gtfs` or `gene_sets`: they choose
/// what is kept, and without them mzLib keeps every gene.
#[derive(Debug, Clone, Default)]
pub struct MapResiduesOptions {
    /// Compara's peptide alignments of its gene trees,
    /// `Compara.<release>.protein_<collection>.aa.fasta.gz` (plain or `.gz`), from
    /// `emf/ensembl-compara/homologies/`: one gapped multiple alignment per tree, read whole by
    /// mzLib's `ComparaGeneTreeAlignment.Load`, which keeps only the canonical proteins of the genes
    /// kept. 908 MB compressed for `protein_default` in release 116. Required.
    pub alignment: PathBuf,
    /// Compara's gene-tree file of the same release,
    /// `<collection>.GeneTree_content.default.e<release>.txt.gz` (plain or `.gz`). It names each
    /// gene's tree and its canonical protein (the row flagged `Y`), the only protein of the gene the
    /// alignment holds. Required.
    pub gene_trees: PathBuf,
    /// Each species, as its Ensembl production name (`"homo_sapiens"`), with its Ensembl
    /// primary-assembly GTF, as for [`build_orthology_with`]. Only the genes in these are kept, and
    /// so only their canonical proteins. Empty with no `gene_sets` (the default): mzLib's default
    /// keeps every gene in the gene-tree file, millions for the vertebrates, and a caveat says so.
    pub gtfs: Vec<(String, PathBuf)>,
    /// Instead of a GTF, a compact gene table mzLib's `EnsemblGeneSetWriter` made from one. A
    /// species may be in `gtfs` or `gene_sets`, not both.
    pub gene_sets: Vec<(String, PathBuf)>,
    /// Also return the kept proteins with their ungapped sequences
    /// ([`ResidueMappings::proteins`]), to align a UniProt sequence to an Ensembl one with
    /// [`align_with`] without a second file. `false` (the default): `proteins` is `None`.
    pub proteins: bool,
    /// Time to allow. `None` (the default) waits: for human, mouse and rat at release 116 the whole
    /// call takes about a minute, nearly all of it reading the alignment.
    pub timeout: Option<Duration>,
}

/// Compara's gene-tree alignment file a [`map_residues_with`] call read, and what was kept.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct ComparaAlignment {
    /// The file's name, e.g. `Compara.116.protein_default.aa.fasta.gz`
    /// (`ComparaGeneTreeAlignment.SourceFileName`).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub source_file_name: String,
    /// Lower-case hex sha256 of the file's bytes as read (the compressed bytes for a `.gz`).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sha256: String,
    /// The Ensembl release in the file name (`Compara.116.`); `None` when the name does not start
    /// `Compara.<release>.`.
    #[serde(default)]
    pub release: Option<String>,
    /// The collection in the file name, e.g. `"protein_default"`; `None` when the name does not
    /// match `Compara.<release>.<collection>.aa.fasta`.
    #[serde(default)]
    pub collection: Option<String>,
    /// Alignments (gene trees) in the file, kept proteins or not, in alignments.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub alignment_count: u64,
    /// Proteins kept, in proteins: the canonical proteins of the genes kept from the gene-tree
    /// file that the alignment holds.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub protein_count: u64,
    /// The sha256 of the gene-tree file whose canonical proteins were kept; the same as
    /// [`ResidueMappings::gene_trees`]' `sha256`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub restricted_to_gene_tree_content: String,
}

/// The proteins a [`map_residues_with`] call kept from the alignment, with
/// [`MapResiduesOptions::proteins`].
///
/// [`Self::columns`] has one row per kept protein, in ordinal order of protein id: `protein_id`
/// (the Ensembl protein stable id, as the alignment's header gives it), `alignment_index` (which
/// alignment holds it, from 0 in file order; two proteins share a column only when they share
/// this), `tree_id` (its gene's Compara gene tree, e.g. `ENSGT00950000183153`) and `sequence` (its
/// residues with the alignment's gaps removed: the sequence every position counts in).
#[derive(Debug, Clone, Deserialize)]
pub struct AlignedProteins {
    /// Proteins kept, in proteins.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The table, one row per protein; its column order is [`Table::names`].
    #[serde(flatten)]
    pub columns: Table,
}

impl AlignedProteins {
    /// Protein id to its ungapped sequence, to align with a UniProt sequence by [`align_with`]. A
    /// sequence's length is the last position a query may name in it.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is missing or not a string column.
    pub fn sequences(&self) -> Result<BTreeMap<String, String>> {
        let ids = self.columns.strings("protein_id")?;
        let sequences = self.columns.strings("sequence")?;
        Ok(ids
            .into_iter()
            .zip(sequences)
            .filter_map(|(id, sequence)| Some((id?, sequence?)))
            .collect())
    }
}

/// What [`map_residues_with`] returns: each query's residue in the other protein, or why there is
/// none.
///
/// [`Self::columns`] has one row per query, in the order given; [`Self::rows`] gives them as typed
/// [`ResidueMapping`]s.
#[derive(Debug, Clone, Deserialize)]
pub struct ResidueMappings {
    /// The alignment file, and what was kept from it.
    pub alignment: ComparaAlignment,
    /// The gene-tree file, and what was kept from it: [`GeneTrees::restricted_to_gene_sets`] is
    /// `None` when no gene set was given, and every gene was kept.
    pub gene_trees: GeneTrees,
    /// Species to the [`GeneSet`] its genes were kept from, the type [`resolve_genes_with`]
    /// reports; empty when no gene set was given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub gene_sets: BTreeMap<String, GeneSet>,
    /// Queries, one row each, in queries.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Outcome to the queries with it, in queries. Every key in [`COLUMN_OUTCOMES`] is present,
    /// zeros included.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub outcome_counts: BTreeMap<String, u64>,
    /// The mapping table, one row per query; its column order is [`Table::names`].
    #[serde(flatten)]
    pub columns: Table,
    /// Whether the call asked for [`Self::proteins`] ([`MapResiduesOptions::proteins`]).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub proteins_included: bool,
    /// The kept proteins with their sequences; `None` without [`MapResiduesOptions::proteins`].
    #[serde(default)]
    pub proteins: Option<AlignedProteins>,
    /// This call's traps: gene trees read unrestricted (every gene kept), an alignment and
    /// gene-tree file named for different releases or for none, queries naming a protein the
    /// alignment did not keep.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// One row of [`ResidueMappings`]: one query, and its answer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResidueMapping {
    /// The query's source protein id, as given.
    pub source_protein_id: String,
    /// The query's 1-based residue position in the source protein, as given.
    pub position: i64,
    /// The query's target protein id, as given.
    pub target_protein_id: String,
    /// One of [`COLUMN_OUTCOMES`], mzLib's `ComparaColumnOutcome` name, verbatim.
    pub outcome: String,
    /// The source residue's 1-based column in the gene tree's alignment. `None`: the two proteins do
    /// not share an alignment (`DifferentAlignments`, `SourceNotInAlignment` or
    /// `TargetNotInAlignment`).
    pub column: Option<i64>,
    /// The 1-based position of the target's residue in that column, in the target's ungapped
    /// sequence. `None`: the outcome is not `Aligned` (a gap is `None`, **never** the nearest
    /// residue).
    pub target_position: Option<i64>,
    /// The target's residue there, one letter as the alignment writes it. It may differ from the
    /// source's residue: a column pairs homologous positions, not identical ones. `None`: the
    /// outcome is not `Aligned`.
    pub target_residue: Option<String>,
}

impl ResidueMappings {
    /// Every row, typed, in the order the queries were given.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn rows(&self) -> Result<Vec<ResidueMapping>> {
        let t = &self.columns;
        if t.is_empty() {
            return Ok(Vec::new());
        }
        let source = t.strings("source_protein_id")?;
        let position = t.integers("position")?;
        let target = t.strings("target_protein_id")?;
        let outcome = t.strings("outcome")?;
        let column = t.integers("column")?;
        let target_position = t.integers("target_position")?;
        let target_residue = t.strings("target_residue")?;
        Ok((0..source.len())
            .map(|i| ResidueMapping {
                source_protein_id: source[i].clone().unwrap_or_default(),
                position: position[i].unwrap_or_default(),
                target_protein_id: target[i].clone().unwrap_or_default(),
                outcome: outcome[i].clone().unwrap_or_default(),
                column: column[i],
                target_position: target_position[i],
                target_residue: target_residue[i].clone(),
            })
            .collect())
    }
}

/// The queries as stdin lines, `source<TAB>position<TAB>target`, refusing what cannot travel. The
/// position is sent as given: mzLib owns the rule for one outside the source.
fn query_lines<S: AsRef<str>, T: AsRef<str>>(queries: &[(S, i64, T)]) -> Result<Vec<String>> {
    if queries.is_empty() {
        return Err(MzLibError::Usage(
            "queries is empty. Give at least one (source_protein_id, position, target_protein_id)."
                .to_owned(),
        ));
    }
    queries
        .iter()
        .enumerate()
        .map(|(index, (source, position, target))| {
            let mut ids = Vec::with_capacity(2);
            for (name, value) in [("source", source.as_ref()), ("target", target.as_ref())] {
                if value.trim().is_empty() {
                    return Err(MzLibError::Usage(format!(
                        "Query {index}'s {name} must be a protein id; got {value:?}."
                    )));
                }
                if value.contains(['\t', '\n', '\r']) {
                    return Err(MzLibError::Usage(format!(
                        "Query {index}'s {name} contains a tab or line break: {value:?}."
                    )));
                }
                ids.push(value.trim());
            }
            Ok(format!("{}\t{position}\t{}", ids[0], ids[1]))
        })
        .collect()
}

/// The arguments and stdin of a `genes map-residues` call: the species lines (none is mzLib's
/// default, every gene kept), a `--` line, then one query per line.
fn map_residues_request<S: AsRef<str>, T: AsRef<str>>(
    queries: &[(S, i64, T)],
    options: &MapResiduesOptions,
) -> Result<(Vec<String>, Option<String>)> {
    let queries = query_lines(queries)?;
    let mut lines = species_lines(&options.gtfs, &options.gene_sets)?;
    for (path, what, example) in [
        (
            &options.alignment,
            "alignment",
            "Compara's Compara.<release>.protein_<collection>.aa.fasta.gz",
        ),
        (
            &options.gene_trees,
            "gene_trees",
            "Compara's <collection>.GeneTree_content.default.e<release>.txt.gz",
        ),
    ] {
        if path.as_os_str().is_empty() {
            return Err(MzLibError::Usage(format!("{what} is required: {example}.")));
        }
    }
    let mut args = vec![
        "genes".to_owned(),
        "map-residues".to_owned(),
        "--alignment".to_owned(),
        optional_path(&options.alignment, "alignment")?,
        "--gene-trees".to_owned(),
        optional_path(&options.gene_trees, "gene_trees")?,
    ];
    if options.proteins {
        args.push("--proteins".to_owned());
    }
    lines.push(SECTION.to_owned());
    lines.extend(queries);
    Ok((args, Some(lines.join("\n") + "\n")))
}

/// Map residues across species through Ensembl Compara's gene-tree alignments.
///
/// Wraps mzLib's `ComparaGeneTreeAlignment` (mzLib #1465): `Load` reads Compara's peptide alignment
/// of every gene tree, keeping the canonical proteins of your species' genes, and `MapResidue`
/// answers, per query, which residue of the target protein sits in the same alignment column as
/// residue `position` of the source protein. Local files only; nothing is downloaded.
///
/// Every query gets exactly one outcome ([`COLUMN_OUTCOMES`]): `Aligned` (with the column, the
/// target position and its residue), `GapInTarget` (no target residue in that column, **never** the
/// nearest one), `DifferentAlignments` (the proteins are in different gene trees), or
/// `SourceNotInAlignment` / `TargetNotInAlignment` (that protein was not kept: only each gene's
/// canonical protein is aligned).
///
/// `queries` are `(source_protein_id, position, target_protein_id)`: Ensembl protein stable ids and
/// a 1-based residue position in the source. One row each, in this order, duplicates included.
/// They travel on stdin, after the species lines.
///
/// **Positions are Ensembl's.** To start from a UniProt position, align the UniProt sequence to the
/// canonical Ensembl protein with [`align_with`] first, and align the target's Ensembl protein to
/// its UniProt sequence after; [`MapResiduesOptions::proteins`] returns the Ensembl sequences for
/// that. **A position outside the source protein fails the whole call**: mzLib's `MapResidue`
/// throws for it rather than returning an outcome, and this crate passes positions through
/// unchecked, so keep them within each source protein.
///
/// Needs the bridge from pyMzLib 0.7.0 or later.
#[doc = include_str!("../docs/reference/genes.map-residues.md")]
///
/// # Examples
///
/// Human p53's Ser15 (Ensembl's canonical TP53 protein) in mouse p53, through release 116's TP53
/// gene-tree alignment:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::proteins::{map_residues_with, MapResiduesOptions};
///
/// let gtfs = vec![
///     ("homo_sapiens".to_owned(), "Homo_sapiens.GRCh38.116.gtf".into()),
///     ("mus_musculus".to_owned(), "Mus_musculus.GRCm39.116.gtf".into()),
///     ("rattus_norvegicus".to_owned(), "Rattus_norvegicus.GRCr8.116.gtf".into()),
/// ];
/// let site = map_residues_with(
///     &[("ENSP00000269305", 15, "ENSMUSP00000104298")],
///     &MapResiduesOptions {
///         alignment: "Compara.116.protein_default.aa.fasta.gz".into(),
///         gene_trees: "vertebrates.GeneTree_content.default.e116.txt.gz".into(),
///         gtfs,
///         ..Default::default()
///     },
/// )?;
/// let row = &site.rows()?[0];
/// assert_eq!(
///     (row.outcome.as_str(), row.column, row.target_position, row.target_residue.as_deref()),
///     ("Aligned", Some(476), Some(18), Some("S"))
/// );
/// assert_eq!(
///     (site.alignment.release.as_deref(), site.gene_trees.collection.as_deref()),
///     (Some("116"), Some("vertebrates"))
/// );
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for no queries, a query whose source or target
/// is blank or holds a tab or line break; a blank `alignment` or `gene_trees`; a species given
/// twice, or a name holding a tab, line break or comma; and when the bridge in use is older than
/// pyMzLib 0.7.0's, which has no `genes map-residues`. From mzLib: a position outside its source
/// protein (0, negative, or past its last residue) is [`MzLibError::Usage`] quoting mzLib, and no
/// row is returned; a file mzLib's reader refuses is [`MzLibError::Bridge`]
/// (`InvalidDataException`, naming the file and line).
#[doc = include_str!("../docs/reference/genes.map-residues.see-also.md")]
pub fn map_residues_with<S: AsRef<str>, T: AsRef<str>>(
    queries: &[(S, i64, T)],
    options: &MapResiduesOptions,
) -> Result<ResidueMappings> {
    let request = map_residues_request(queries, options)?;
    bridge::require_verb("genes map-residues", MAP_RESIDUES_SINCE)?;
    call(request, options.timeout)
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

    // ---- annotate-go and update-go ----------------------------------------------------------------

    const GO: &str = include_str!("../tests/fixtures/proteins_annotate_go_pxd036557.json");
    const GO_OUT: &str = include_str!("../tests/fixtures/proteins_annotate_go_out.json");
    const GO_UPDATE: &str = include_str!("../tests/fixtures/proteins_update_go.json");

    fn go_options() -> GoAnnotateOptions {
        GoAnnotateOptions {
            go_obo: "go-pxd036557.obo".into(),
            ..Default::default()
        }
    }

    /// Row indices whose `protein_group` is `group`.
    fn rows_of(go: &GoAnnotations, group: &str) -> Vec<usize> {
        let groups = go.columns.strings("protein_group").unwrap();
        (0..groups.len())
            .filter(|&i| groups[i].as_deref() == Some(group))
            .collect()
    }

    #[test]
    fn every_non_decoy_group_is_annotated_and_the_decoy_skipped() {
        let go: GoAnnotations = serde_json::from_str(GO).unwrap();
        assert_eq!(
            (go.table_row_count, go.decoy_group_count, go.group_count),
            (6, 1, 5)
        );
        assert_eq!((go.row_count, go.returned_count), (563, 563));
        assert_eq!(go.columns.rows(), 563);
        assert!(!go.truncated);
        let mut groups: Vec<String> = Vec::new();
        for group in go
            .columns
            .strings("protein_group")
            .unwrap()
            .into_iter()
            .flatten()
        {
            if !groups.contains(&group) {
                groups.push(group);
            }
        }
        assert_eq!(
            groups,
            ["P68363", "P05141", "P0C0S5|Q71UI9", "P02769", "P63104"]
        );
    }

    #[test]
    fn a_contaminant_is_one_term_less_row_that_says_why() {
        let go: GoAnnotations = serde_json::from_str(GO).unwrap();
        let albumin = rows_of(&go, "P02769");
        assert_eq!(albumin.len(), 1);
        let i = albumin[0];
        assert_eq!(
            go.columns.strings("annotation_status").unwrap()[i].as_deref(),
            Some("contaminant")
        );
        assert_eq!(go.columns.strings("go_id").unwrap()[i], None);
        assert_eq!(go.columns.strings("aspect").unwrap()[i], None);
        assert_eq!(go.columns.booleans("inherited").unwrap()[i], None);
        assert_eq!(go.columns.booleans("propagated").unwrap()[i], None);
        assert_eq!(go.columns.integers("n_with").unwrap()[i], Some(0));
        assert_eq!(
            string_lists(&go.columns, "accession_used").unwrap()[i],
            Some(Vec::new())
        );
    }

    #[test]
    fn the_union_keeps_both_histones_and_the_consensus_is_a_filter() {
        let go: GoAnnotations = serde_json::from_str(GO).unwrap();
        let histones = rows_of(&go, "P0C0S5|Q71UI9");
        let n_with = go.columns.integers("n_with").unwrap();
        let n_members = go.columns.integers("n_members").unwrap();
        let consensus = histones
            .iter()
            .filter(|&&i| n_with[i] == n_members[i])
            .count();
        assert_eq!((histones.len(), consensus), (106, 57));
    }

    #[test]
    fn evidence_is_kept_per_member() {
        let go: GoAnnotations = serde_json::from_str(GO).unwrap();
        let names = go.columns.strings("go_name").unwrap();
        let evidence = string_list_maps(&go.columns, "evidence_by_member").unwrap();
        let nucleosome = rows_of(&go, "P0C0S5|Q71UI9")
            .into_iter()
            .find(|&i| names[i].as_deref() == Some("nucleosome"))
            .unwrap();
        let mut want = BTreeMap::new();
        want.insert("P0C0S5".to_owned(), vec!["ECO:0000353".to_owned()]);
        want.insert("Q71UI9".to_owned(), vec!["ECO:0000353".to_owned()]);
        assert_eq!(evidence[nucleosome], Some(want));
        assert!(matches!(
            string_list_maps(&go.columns, "go_id"),
            Err(MzLibError::Protocol(_))
        ));
        assert!(matches!(
            string_list_maps(&go.columns, "nonesuch"),
            Err(MzLibError::Usage(_))
        ));
    }

    #[test]
    fn the_column_names_are_mzlibs_schema() {
        let go: GoAnnotations = serde_json::from_str(GO).unwrap();
        let names = go.columns.names();
        assert_eq!(names.len(), 19);
        assert_eq!(names[0], "protein_group");
        assert_eq!(
            names[16..],
            ["go_release", "go_obo_sha256", "annotation_db_sha256"]
        );
    }

    #[test]
    fn provenance_pins_every_input() {
        let go: GoAnnotations = serde_json::from_str(GO).unwrap();
        assert_eq!(go.go.source_file_name, "go-pxd036557.obo");
        assert_eq!(go.go.release.as_deref(), Some("releases/2026-07-26"));
        assert_eq!(go.go.term_count, 412);
        assert_eq!(
            go.annotation_database.file_type.as_deref(),
            Some("UniProtXml")
        );
        assert_eq!(go.annotation_database.protein_count, 5);
        assert_eq!(go.header["source_file_sha256"], go.groups_file_sha256);
        assert_eq!(
            go.header["annotation_db_sha256"],
            go.annotation_database.sha256
        );
        assert_eq!(go.header["counter_q_value_max"], "0.01");
        assert_eq!(go.header["n_multi_member_groups"], "1");
        assert!(go
            .columns
            .strings("go_obo_sha256")
            .unwrap()
            .iter()
            .all(|sha| sha.as_deref() == Some(go.go.sha256.as_str())));
    }

    #[test]
    fn categories_join_on_go_id_and_say_the_subcategory() {
        let go: GoAnnotations = serde_json::from_str(GO).unwrap();
        let cats = go.categories.as_ref().unwrap();
        assert_eq!(
            (
                cats.map_name.as_str(),
                cats.map_version.as_str(),
                cats.anchor_count,
                cats.row_count
            ),
            ("organelle", "1", 9, 30)
        );
        let ids = cats.columns.strings("go_id").unwrap();
        let subs = cats.columns.strings("subcategory").unwrap();
        let sub_of = |id: &str| {
            let i = ids.iter().position(|x| x.as_deref() == Some(id)).unwrap();
            subs[i].clone()
        };
        assert_eq!(
            sub_of("GO:0005743").as_deref(),
            Some("mitochondrion:inner_membrane")
        );
        assert_eq!(sub_of("GO:0005739"), None);
        assert!(go.written.is_none() && go.categories_written.is_none());
        assert!(!go.skip_unknown_go_ids && go.unresolved_go_ids.is_empty());
    }

    #[test]
    fn out_with_limit_zero_returns_only_the_summary() {
        let go: GoAnnotations = serde_json::from_str(GO_OUT).unwrap();
        assert_eq!(
            go.written,
            Some(WrittenTable {
                path: "go_annotations.tsv".to_owned(),
                row_count: Some(563)
            })
        );
        assert_eq!((go.returned_count, go.columns.rows()), (0, 0));
        assert!(go.truncated && go.categories.is_none());
    }

    #[test]
    fn the_required_inputs_go_on_argv_and_the_optional_ones_follow() {
        let args = annotate_go_request(Path::new("groups.tsv"), Path::new("db.xml"), &go_options())
            .unwrap();
        assert_eq!(
            args,
            [
                "proteins",
                "annotate-go",
                "--groups",
                "groups.tsv",
                "--database",
                "db.xml",
                "--go-obo",
                "go-pxd036557.obo",
                "--offset",
                "0"
            ]
        );
        let args = annotate_go_request(
            Path::new("groups.tsv"),
            Path::new("db.xml"),
            &GoAnnotateOptions {
                category_map: Some("map.tsv".into()),
                skip_unknown_go_ids: true,
                out: Some("go.tsv".into()),
                categories_out: Some("cats.tsv".into()),
                limit: Some(0),
                offset: 5,
                ..go_options()
            },
        )
        .unwrap();
        for (option, value) in [
            ("--category-map", "map.tsv"),
            ("--out", "go.tsv"),
            ("--categories-out", "cats.tsv"),
            ("--limit", "0"),
            ("--offset", "5"),
        ] {
            let at = args.iter().position(|a| a == option).unwrap();
            assert_eq!(args[at + 1], value);
        }
        assert!(args.contains(&"--skip-unknown-go-ids".to_owned()));
    }

    #[test]
    fn a_missing_go_obo_or_input_is_refused_before_anything_is_spawned() {
        let error = annotate_go_request(
            Path::new("groups.tsv"),
            Path::new("db.xml"),
            &GoAnnotateOptions::default(),
        )
        .unwrap_err();
        assert!(error.to_string().contains("update_go"), "{error}");
        for (groups, database) in [("", "db.xml"), ("groups.tsv", " ")] {
            let error = annotate_go_request(Path::new(groups), Path::new(database), &go_options())
                .unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)), "{error}");
        }
    }

    #[test]
    fn update_go_reports_the_release_now_on_disk() {
        let update: GoUpdate = serde_json::from_str(GO_UPDATE).unwrap();
        assert_eq!(
            update_go_request(Path::new("go.obo")).unwrap(),
            ["proteins", "update-go", "--go-obo", "go.obo"]
        );
        assert_eq!(
            (
                update.existed_before,
                update.previous_sha256.as_deref(),
                update.changed
            ),
            (false, None, true)
        );
        assert_eq!(update.go.release.as_deref(), Some("releases/2026-07-26"));
        assert_eq!(update.go.term_count, 48340);
        assert_eq!(update.url, "https://purl.obolibrary.org/obo/go.obo");
        assert!(matches!(
            update_go_request(Path::new("")),
            Err(MzLibError::Usage(_))
        ));
    }

    #[test]
    fn the_statuses_are_mzlibs_four() {
        assert_eq!(
            ANNOTATION_STATUSES,
            ["annotated", "no_go_terms", "no_entry", "contaminant"]
        );
    }

    // ---- align (mzLib #1428) --------------------------------------------------------------------

    fn recorded<T: serde::de::DeserializeOwned>(name: &str) -> T {
        let text = std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("fixtures")
                .join(name),
        )
        .unwrap();
        serde_json::from_str(&text).unwrap()
    }

    const DELETION: (&str, &str) = ("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ");

    #[test]
    fn align_sends_the_pair_on_stdin_and_no_option_it_was_not_given() {
        let (args, stdin) = align_request(&[DELETION], &AlignOptions::default()).unwrap();
        assert_eq!(args, ["proteins", "align"]);
        assert_eq!(
            stdin.as_deref(),
            Some("MKTAYIAKQRQISFVKSHFSRQ\tMKTAYIAKISFVKSHFSRQ\n")
        );
        assert_eq!(AlignOptions::default().timeout, None);
    }

    #[test]
    fn align_projects_mzlibs_counts_and_rows() {
        let alns: Alignments = recorded("proteins_align_deletion.json");
        let aln = &alns.alignments().unwrap()[0];
        assert_eq!(aln.aligned_source, "MKTAYIAKQRQISFVKSHFSRQ");
        assert_eq!(aln.aligned_target, "MKTAYIAK---ISFVKSHFSRQ");
        assert_eq!((aln.score, aln.identical, aln.substituted), (80, 19, 0));
        assert_eq!(
            (
                aln.source_residues_facing_gap,
                aln.target_residues_facing_gap
            ),
            (3, 0)
        );
        assert_eq!((aln.source_length, aln.target_length), (22, 19));
        assert_eq!(aln.aligner_id, alns.aligner_id);
        assert_eq!(aln.caveats, alns.caveats);
        assert_eq!(alns.residue_map.as_ref().unwrap().row_count, 41);
    }

    #[test]
    fn a_residue_facing_a_gap_maps_to_none_never_a_neighbour() {
        let alns: Alignments = recorded("proteins_align_deletion.json");
        let aln = &alns.alignments().unwrap()[0];
        let mapped: Vec<Option<usize>> = (8..=12)
            .map(|p| aln.map_source_to_target(p).unwrap())
            .collect();
        assert_eq!(mapped, [Some(8), None, None, None, Some(9)]);
        assert_eq!(aln.map_target_to_source(9).unwrap(), Some(12));
    }

    #[test]
    fn a_position_outside_the_sequence_is_refused() {
        let alns: Alignments = recorded("proteins_align_deletion.json");
        let aln = &alns.alignments().unwrap()[0];
        for position in [0, 23] {
            let error = aln.map_source_to_target(position).unwrap_err();
            assert!(
                matches!(&error, MzLibError::Usage(m) if m.contains("1-based")),
                "{position}: {error:?}"
            );
        }
        assert!(aln.map_target_to_source(20).is_err());
    }

    #[test]
    fn each_align_option_given_is_sent_and_nothing_else() {
        let cases: [(AlignOptions, &[&str]); 4] = [
            (
                AlignOptions {
                    gap_open: Some(5),
                    ..Default::default()
                },
                &["--gap-open", "5"],
            ),
            (
                AlignOptions {
                    gap_extend: Some(2),
                    ..Default::default()
                },
                &["--gap-extend", "2"],
            ),
            (
                AlignOptions {
                    end_gaps: Some(EndGaps::Charged),
                    ..Default::default()
                },
                &["--end-gaps", "charged"],
            ),
            (
                AlignOptions {
                    max_cells: Some(1000),
                    ..Default::default()
                },
                &["--max-cells", "1000"],
            ),
        ];
        for (options, sent) in cases {
            let (args, _) = align_request(&[DELETION], &options).unwrap();
            assert_eq!(args[2..], *sent);
        }
    }

    #[test]
    fn charged_end_gaps_come_back_in_the_envelope() {
        let alns: Alignments = recorded("proteins_align_charged.json");
        assert_eq!(alns.end_gaps, "charged");
        assert_eq!(
            alns.aligner_id,
            "global-affine;BLOSUM62;open=11;extend=1;end-gaps=charged"
        );
        let aln = &alns.alignments().unwrap()[0];
        assert_eq!(
            (aln.aligned_target.as_str(), aln.score),
            ("----PEPTIDE", 24)
        );
    }

    #[test]
    fn align_with_keeps_input_order_and_one_alignment_per_pair() {
        let (_, stdin) = align_request(
            &[DELETION, ("MLLLPEPTIDE", "PEPTIDE")],
            &AlignOptions::default(),
        )
        .unwrap();
        assert_eq!(
            stdin.as_deref(),
            Some("MKTAYIAKQRQISFVKSHFSRQ\tMKTAYIAKISFVKSHFSRQ\nMLLLPEPTIDE\tPEPTIDE\n")
        );
        let alns: Alignments = recorded("proteins_align_many.json");
        assert_eq!(alns.record_count, 2);
        let each = alns.alignments().unwrap();
        assert_eq!(
            each.iter().map(|a| a.pair_index).collect::<Vec<_>>(),
            [0, 1]
        );
        assert_eq!(
            alns.columns.integers("score").unwrap(),
            [Some(80), Some(39)]
        );
        assert_eq!(each[1].map_target_to_source(1).unwrap(), Some(5));
        assert_eq!(
            (
                alns.matrix.as_str(),
                alns.gap_open,
                alns.gap_extend,
                alns.max_cells
            ),
            ("BLOSUM62", 11, 1, 200_000_000)
        );
    }

    #[test]
    fn residue_map_false_sends_the_flag_and_leaves_no_map() {
        let options = AlignOptions {
            residue_map: false,
            ..Default::default()
        };
        let (args, _) = align_request(&[DELETION], &options).unwrap();
        assert_eq!(args, ["proteins", "align", "--no-residue-map"]);
        let alns: Alignments = recorded("proteins_align_no_map.json");
        assert!(alns.residue_map.is_none());
        let aln = &alns.alignments().unwrap()[0];
        assert_eq!(aln.source_to_target, None);
        let error = aln.map_source_to_target(1).unwrap_err();
        assert!(error.to_string().contains("no residue map"));
    }

    #[test]
    fn real_isoforms_shift_every_position_after_the_missing_exon() {
        let alns: Alignments = recorded("proteins_align_dync1i2.json");
        let each = alns.alignments().unwrap();
        let (two_b, two_c) = (&each[0], &each[1]);
        let mapped: Vec<(Option<usize>, Option<usize>)> = [75, 76, 82, 100, 120, 133]
            .iter()
            .map(|&p| {
                (
                    two_b.map_source_to_target(p).unwrap(),
                    two_c.map_source_to_target(p).unwrap(),
                )
            })
            .collect();
        assert_eq!(
            mapped,
            [
                (Some(75), Some(75)),
                (None, None),
                (Some(76), Some(76)),
                (Some(94), Some(94)),
                (Some(114), None),
                (Some(127), Some(107))
            ]
        );
        assert_eq!(
            (
                two_b.identical,
                two_b.substituted,
                two_b.source_residues_facing_gap
            ),
            (632, 0, 6)
        );
    }

    #[test]
    fn real_orthologs_substituted_residues_keep_their_own_letters() {
        let alns: Alignments = recorded("proteins_align_albumin.json");
        let albumin = &alns.alignments().unwrap()[0];
        assert_eq!(albumin.map_source_to_target(58).unwrap(), Some(58));
        assert_eq!(
            (
                albumin.identical,
                albumin.substituted,
                albumin.source_residues_facing_gap
            ),
            (465, 142, 2)
        );
        let map = alns.residue_map.as_ref().unwrap();
        let residues = map.columns.strings("residue").unwrap();
        let mapped = map.columns.strings("mapped_residue").unwrap();
        let positions = map.columns.integers("position").unwrap();
        let sequences = map.columns.strings("sequence").unwrap();
        let row = (0..positions.len())
            .find(|&i| positions[i] == Some(26) && sequences[i].as_deref() == Some("source"))
            .unwrap();
        assert_eq!(
            (residues[row].as_deref(), mapped[row].as_deref()),
            (Some("A"), Some("T")),
            "a pair can be two different amino acids"
        );
    }

    #[test]
    fn bad_pairs_are_refused_before_the_bridge() {
        let empty: [(&str, &str); 0] = [];
        let cases: Vec<(Vec<(&str, &str)>, &str)> = vec![
            (empty.to_vec(), "At least one pair"),
            (vec![("PEP\tTIDE", "PEPTIDE")], "tab or line break"),
            (vec![("PEPTIDE", "PEP\nTIDE")], "tab or line break"),
            (vec![("", "")], "two empty sequences"),
            (vec![(" ", "")], "two empty sequences"),
        ];
        for (pairs, message) in cases {
            let error = align_request(&pairs, &AlignOptions::default()).unwrap_err();
            assert!(
                matches!(&error, MzLibError::Usage(m) if m.contains(message)),
                "{message}: {error:?}"
            );
        }
    }

    #[test]
    fn an_empty_sequence_on_one_side_is_sent() {
        let (_, stdin) = align_request(&[("", "PEP")], &AlignOptions::default()).unwrap();
        assert_eq!(stdin.as_deref(), Some("\tPEP\n"));
    }

    // ---- parse-accession (mzLib #1382) ----------------------------------------------------------

    #[test]
    fn parse_accessions_sends_each_accession_unchanged() {
        let lines = accession_lines(&[" P04637", "P04637 ", "p04637", "P04637|P02768"]).unwrap();
        assert_eq!(
            lines.join("\n") + "\n",
            " P04637\nP04637 \np04637\nP04637|P02768\n",
            "not stripped, not split"
        );
    }

    #[test]
    fn parse_accessions_projects_mzlibs_reading() {
        let parsed: ParsedAccessions = recorded("proteins_parse_accession.json");
        assert_eq!(parsed.record_count, 9);
        assert_eq!(
            parsed.namespace_counts,
            BTreeMap::from([
                ("refseq".to_owned(), 1),
                ("uniprot".to_owned(), 4),
                ("unrecognized".to_owned(), 4)
            ])
        );
        assert_eq!(parsed.with_applied_variants_count, 2);
        let entries = parsed.columns.strings("entry_accession").unwrap();
        assert_eq!(
            entries[..5],
            ["P04637", "P04637", "Q13409", "NP_000537", "P04637"].map(|s| Some(s.to_owned()))
        );
        let namespaces = parsed.columns.strings("namespace").unwrap();
        assert_eq!(
            namespaces[4..7],
            ["uniprot", "unrecognized", "unrecognized"].map(|s| Some(s.to_owned()))
        );
        assert!(namespaces
            .iter()
            .all(|n| NAMESPACES.contains(&n.as_deref().unwrap())));
        assert_eq!(
            parsed.columns.strings("parent_accession").unwrap()[3].as_deref(),
            Some("NP_000537.3")
        );
        assert_eq!(parsed.entry_of().unwrap()["NP_000537.3_R72P"], "NP_000537");
        assert_eq!(
            parsed.entry_of().unwrap()["P04637_2"],
            "P04637_2",
            "unrecognized maps to itself"
        );
    }

    #[test]
    fn bad_accession_lists_are_refused_before_the_bridge() {
        let empty: [&str; 0] = [];
        for (accessions, message) in [
            (&empty[..], "At least one accession"),
            (&["P04637", " "][..], "blank entry"),
            (&["P04637\nP02768"][..], "line break"),
        ] {
            let error = accession_lines(accessions).unwrap_err();
            assert!(
                matches!(&error, MzLibError::Usage(m) if m.contains(message)),
                "{message}: {error:?}"
            );
        }
    }

    // ---- genes orthology (mzLib #1381) -----------------------------------------------------------

    const HUMAN: &str = "homo_sapiens";
    const MOUSE: &str = "mus_musculus";
    const RAT: &str = "rattus_norvegicus";
    const TREES: &str = "vertebrates.GeneTree_content.default.e116.txt.gz";

    fn orthology_options() -> OrthologyOptions {
        OrthologyOptions {
            gtfs: vec![
                (HUMAN.to_owned(), "Homo_sapiens.GRCh38.116.gtf".into()),
                (MOUSE.to_owned(), "Mus_musculus.GRCm39.116.gtf".into()),
                (RAT.to_owned(), "Rattus_norvegicus.GRCr8.116.gtf".into()),
            ],
            dumps: [HUMAN, MOUSE, RAT]
                .iter()
                .map(|s| {
                    format!("homologies/{s}/Compara.116.protein_default.homologies.tsv.gz").into()
                })
                .collect(),
            gene_trees: TREES.into(),
            ..Default::default()
        }
    }

    #[test]
    fn orthology_sends_the_species_lines_then_the_dumps_on_stdin() {
        let (args, stdin) = orthology_request("116", &orthology_options()).unwrap();
        assert_eq!(
            args,
            [
                "genes",
                "orthology",
                "--release",
                "116",
                "--gene-trees",
                TREES
            ]
        );
        assert_eq!(
            stdin.unwrap(),
            "homo_sapiens\tgtf\tHomo_sapiens.GRCh38.116.gtf\n\
             mus_musculus\tgtf\tMus_musculus.GRCm39.116.gtf\n\
             rattus_norvegicus\tgtf\tRattus_norvegicus.GRCr8.116.gtf\n\
             --\n\
             homologies/homo_sapiens/Compara.116.protein_default.homologies.tsv.gz\n\
             homologies/mus_musculus/Compara.116.protein_default.homologies.tsv.gz\n\
             homologies/rattus_norvegicus/Compara.116.protein_default.homologies.tsv.gz\n"
        );
        assert_eq!(OrthologyOptions::default().timeout, None);
    }

    #[test]
    fn every_orthology_option_given_is_sent_in_its_wire_spelling() {
        let options = OrthologyOptions {
            restrict_gene_trees: true,
            out: Some("compara-116".into()),
            pair_status: Some((MOUSE.to_owned(), HUMAN.to_owned())),
            ..orthology_options()
        };
        let (args, _) = orthology_request("116", &options).unwrap();
        assert!(args.contains(&"--restrict-gene-trees".to_owned()));
        let after = |flag: &str| args[args.iter().position(|a| a == flag).unwrap() + 1].clone();
        assert_eq!(after("--out"), "compara-116");
        assert_eq!(
            after("--pair-status"),
            "mus_musculus,homo_sapiens",
            "from, then to"
        );
    }

    #[test]
    fn gene_sets_are_sent_as_gene_set_lines_beside_gtf_lines() {
        let options = OrthologyOptions {
            gtfs: vec![(HUMAN.to_owned(), "human.gtf.gz".into())],
            gene_sets: vec![
                (MOUSE.to_owned(), "mouse.genes.tsv".into()),
                (RAT.to_owned(), "rat.genes.tsv".into()),
            ],
            ..orthology_options()
        };
        let (_, stdin) = orthology_request("116", &options).unwrap();
        let stdin = stdin.unwrap();
        let species: Vec<&str> = stdin.split("--\n").next().unwrap().lines().collect();
        assert_eq!(
            species,
            [
                "homo_sapiens\tgtf\thuman.gtf.gz",
                "mus_musculus\tgene-set\tmouse.genes.tsv",
                "rattus_norvegicus\tgene-set\trat.genes.tsv"
            ]
        );
    }

    #[test]
    fn malformed_orthology_arguments_are_usage_errors_before_the_bridge_runs() {
        let base = orthology_options;
        let cases: Vec<(&str, OrthologyOptions, &str)> = vec![
            (
                "116",
                OrthologyOptions {
                    gtfs: Vec::new(),
                    ..base()
                },
                "At least one species",
            ),
            (
                "116",
                OrthologyOptions {
                    gene_sets: vec![(HUMAN.to_owned(), "h.tsv".into())],
                    ..base()
                },
                "'homo_sapiens' is given twice",
            ),
            (
                "116",
                OrthologyOptions {
                    gtfs: vec![("homo\tsapiens".to_owned(), "x.gtf".into())],
                    ..base()
                },
                "cannot hold a tab",
            ),
            (
                "116",
                OrthologyOptions {
                    gtfs: vec![(String::new(), "x.gtf".into())],
                    ..base()
                },
                "non-empty Ensembl production name",
            ),
            (
                "116",
                OrthologyOptions {
                    gtfs: vec![(HUMAN.to_owned(), "".into())],
                    ..base()
                },
                "must be a non-empty path",
            ),
            (
                "116",
                OrthologyOptions {
                    dumps: Vec::new(),
                    ..base()
                },
                "dumps is empty",
            ),
            (
                "116",
                OrthologyOptions {
                    dumps: vec!["--".into()],
                    ..base()
                },
                "is not a usable homology dump path",
            ),
            (
                "116",
                OrthologyOptions {
                    gene_trees: PathBuf::new(),
                    ..base()
                },
                "gene_trees is required",
            ),
            (
                "116",
                OrthologyOptions {
                    pair_status: Some((HUMAN.to_owned(), "mus,musculus".to_owned())),
                    ..base()
                },
                "cannot hold a tab, line break or comma",
            ),
            (
                "116",
                OrthologyOptions {
                    out: Some(PathBuf::from(" ")),
                    ..base()
                },
                "out must be a non-empty path",
            ),
            ("", base(), "release must be an Ensembl release"),
            ("11 6", base(), "release must be an Ensembl release"),
        ];
        for (release, options, message) in cases {
            let error = orthology_request(release, &options).unwrap_err();
            assert!(
                matches!(&error, MzLibError::Usage(m) if m.contains(message)),
                "{message}: {error:?}"
            );
        }
    }

    #[test]
    fn the_snapshot_says_what_went_in() {
        let snap: OrthologySnapshot = recorded("genes_orthology_pair_status.json");
        assert_eq!(snap.release, "116");
        assert_eq!(snap.species, [HUMAN, MOUSE, RAT]);
        assert_eq!(snap.format, "ensembl-orthology-snapshot");
        assert_eq!(snap.format_version, 1);
        assert_eq!(snap.snapshot_id.len(), 64);
        assert!(snap.snapshot_id.starts_with("5487ebcabe48"));
        assert_eq!(snap.identical_duplicates_dropped, 0);
        assert_eq!(snap.gene_sets[HUMAN].gene_count, 4);
        assert_eq!(snap.gene_trees.source_file_name, TREES);
        assert_eq!(snap.gene_trees.release.as_deref(), Some("116"));
        assert_eq!(
            snap.gene_sets[HUMAN].genome_build, None,
            "the test GTFs carry no #!genome-build"
        );
        assert_eq!(snap.gene_trees.collection.as_deref(), Some("vertebrates"));
        assert_eq!(snap.gene_trees.member_count, 9);
        assert_eq!(snap.gene_trees.restricted_to_gene_sets, None);
        assert!(snap.written.is_none());
    }

    #[test]
    fn every_human_mouse_row_came_from_the_mouse_dump_and_none_was_lost() {
        let snap: OrthologySnapshot = recorded("genes_orthology_pair_status.json");
        let kept: Vec<(&str, u64)> = snap
            .dumps
            .iter()
            .map(|(s, d)| (s.as_str(), d.rows_kept))
            .collect();
        assert_eq!(kept, [(HUMAN, 1), (MOUSE, 5), (RAT, 1)]);
        let counts: Vec<(&str, &str, u64, u64)> = snap
            .pairs
            .iter()
            .map(|p| {
                (
                    p.species_a.as_str(),
                    p.species_b.as_str(),
                    p.ortholog_count,
                    p.paralog_count,
                )
            })
            .collect();
        assert_eq!(
            counts,
            [
                (HUMAN, HUMAN, 0, 0),
                (HUMAN, MOUSE, 3, 0),
                (HUMAN, RAT, 1, 0),
                (MOUSE, MOUSE, 0, 1),
                (MOUSE, RAT, 2, 0),
                (RAT, RAT, 0, 0)
            ]
        );
        assert_eq!(snap.pair(MOUSE, HUMAN).unwrap().row_count, 3);
        assert_eq!(snap.pair(HUMAN, MOUSE).unwrap().row_count, 3);
    }

    #[test]
    fn a_pair_of_a_species_not_in_the_snapshot_is_a_usage_error() {
        let snap: OrthologySnapshot = recorded("genes_orthology_pair_status.json");
        let error = snap.pair(HUMAN, "danio_rerio").unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    #[test]
    fn pair_status_is_mzlibs_status_for_every_gene() {
        let snap: OrthologySnapshot = recorded("genes_orthology_pair_status.json");
        let status = snap.pair_status.unwrap();
        assert_eq!(
            (status.from_species.as_str(), status.to_species.as_str()),
            (HUMAN, MOUSE)
        );
        assert_eq!(status.row_count, 4);
        assert_eq!(status.status_counts["has_ortholog"], 2);
        assert_eq!(status.status_counts["not_in_gene_set"], 0);
        assert_eq!(
            status.status_counts.keys().count(),
            ORTHOLOGY_STATUSES.len(),
            "every status present, zeros included"
        );
        assert_eq!(
            status.genes("has_ortholog").unwrap(),
            ["ENSG00000000001", "ENSG00000000002"]
        );
        assert_eq!(
            status.columns.names(),
            ["gene_id", "gene_name", "gene_biotype", "status"]
        );
        assert_eq!(
            status.columns.strings("gene_name").unwrap(),
            [None, None, None, None],
            "no gene_name in the GTF: None, not an empty string"
        );
        let error = status.genes("no_ortholog").unwrap_err();
        assert!(matches!(&error, MzLibError::Usage(m) if m.contains("status must be one of")));
    }

    #[test]
    fn the_orthology_caveats_name_the_shared_dump_name_and_the_unrestricted_trees() {
        let snap: OrthologySnapshot = recorded("genes_orthology_pair_status.json");
        assert!(snap.caveats.iter().any(
            |c| c.contains("share the file name Compara.116.protein_default.homologies.tsv.gz")
        ));
        assert!(snap
            .caveats
            .iter()
            .any(|c| c.contains("read unrestricted: 9 genes")));
        let names: std::collections::BTreeSet<&str> = snap
            .dumps
            .values()
            .map(|d| d.source_file_name.as_str())
            .collect();
        assert_eq!(
            names,
            ["Compara.116.protein_default.homologies.tsv.gz"].into()
        );
    }

    #[test]
    fn a_written_snapshot_lists_every_file_but_the_manifest() {
        let snap: OrthologySnapshot = recorded("genes_orthology_written.json");
        let mut sets: Vec<String> = snap.gene_sets.values().map(|g| g.sha256.clone()).collect();
        sets.sort();
        assert_eq!(snap.gene_trees.restricted_to_gene_sets, Some(sets));
        let written = snap.written.unwrap();
        assert!(written.manifest().ends_with("manifest.json"));
        assert!(written.files.iter().all(|f| f.path != "manifest.json"));
        assert!(written
            .files
            .iter()
            .all(|f| f.sha256.len() == 64 && f.bytes > 0));
        let pairs: Vec<Option<u64>> = written
            .files
            .iter()
            .filter(|f| f.path.starts_with("pairs/"))
            .map(|f| f.rows)
            .collect();
        assert_eq!(pairs, [0, 3, 1, 1, 2, 0].map(Some));
        assert_eq!(
            written.files.last().map(|f| (f.path.as_str(), f.rows)),
            Some(("views.sql", None))
        );
        assert_eq!(written.files.len(), 13);
        let members = written
            .files
            .iter()
            .find(|f| f.path == "members/homo_sapiens.parquet")
            .unwrap();
        assert_eq!(members.rows, Some(3), "ENSG00000000004 is in no tree");
        assert!(written.manifest().starts_with(&written.directory));
    }

    #[test]
    fn restricting_the_gene_trees_does_not_change_the_snapshot() {
        let plain: OrthologySnapshot = recorded("genes_orthology_pair_status.json");
        let restricted: OrthologySnapshot = recorded("genes_orthology_written.json");
        assert_eq!(plain.snapshot_id, restricted.snapshot_id);
        assert_eq!(plain.pairs, restricted.pairs);
        assert_eq!(
            (
                plain.gene_trees.member_count,
                restricted.gene_trees.member_count
            ),
            (9, 8)
        );
        assert!(!restricted
            .caveats
            .iter()
            .any(|c| c.contains("unrestricted")));
    }

    #[test]
    fn an_orthology_payload_without_optional_blocks_projects_to_none() {
        let mut data: serde_json::Value = recorded("genes_orthology_written.json");
        data["written"] = serde_json::Value::Null;
        let snap: OrthologySnapshot = serde_json::from_value(data).unwrap();
        assert!(snap.written.is_none() && snap.pair_status.is_none());
    }

    // ---- genes map-residues ---------------------------------------------------------------------

    const TP53: &str = "ENSP00000269305";
    const TRP53: &str = "ENSMUSP00000104298";

    fn map_options() -> MapResiduesOptions {
        MapResiduesOptions {
            alignment: "Compara.116.protein_default.aa.fasta.gz".into(),
            gene_trees: "vertebrates.GeneTree_content.default.e116.txt.gz".into(),
            gtfs: vec![
                ("homo_sapiens".into(), "Homo_sapiens.GRCh38.116.gtf".into()),
                ("mus_musculus".into(), "Mus_musculus.GRCm39.116.gtf".into()),
                (
                    "rattus_norvegicus".into(),
                    "Rattus_norvegicus.GRCr8.116.gtf".into(),
                ),
            ],
            ..Default::default()
        }
    }

    #[test]
    fn map_residues_sends_the_species_lines_then_the_queries_and_no_option_it_was_not_given() {
        let (args, stdin) = map_residues_request(&[(TP53, 15, TRP53)], &map_options()).unwrap();
        assert_eq!(
            args,
            [
                "genes",
                "map-residues",
                "--alignment",
                "Compara.116.protein_default.aa.fasta.gz",
                "--gene-trees",
                "vertebrates.GeneTree_content.default.e116.txt.gz"
            ]
        );
        assert_eq!(
            stdin.unwrap(),
            "homo_sapiens\tgtf\tHomo_sapiens.GRCh38.116.gtf\n\
             mus_musculus\tgtf\tMus_musculus.GRCm39.116.gtf\n\
             rattus_norvegicus\tgtf\tRattus_norvegicus.GRCr8.116.gtf\n\
             --\nENSP00000269305\t15\tENSMUSP00000104298\n"
        );
        assert_eq!(map_options().timeout, None);

        let (args, _) = map_residues_request(
            &[(TP53, 15, TRP53)],
            &MapResiduesOptions {
                proteins: true,
                ..map_options()
            },
        )
        .unwrap();
        assert_eq!(args.last().map(String::as_str), Some("--proteins"));
    }

    #[test]
    fn without_gene_sets_no_species_line_is_sent_and_mzlib_keeps_every_gene() {
        let (_, stdin) = map_residues_request(
            &[(TP53, 15, TRP53)],
            &MapResiduesOptions {
                gtfs: Vec::new(),
                ..map_options()
            },
        )
        .unwrap();
        assert_eq!(
            stdin.unwrap(),
            "--\nENSP00000269305\t15\tENSMUSP00000104298\n",
            "stdin starts with the separator"
        );
    }

    #[test]
    fn a_position_is_sent_as_given_never_checked_against_a_length() {
        // mzLib owns the rule (MapResidue throws for a position outside the source).
        let (_, stdin) = map_residues_request(
            &[
                (TP53, 0, TRP53),
                (TP53, -3, TRP53),
                (TP53, 1_000_000, TRP53),
            ],
            &map_options(),
        )
        .unwrap();
        let stdin = stdin.unwrap();
        let queries: Vec<&str> = stdin.split("--\n").nth(1).unwrap().lines().collect();
        assert_eq!(
            queries,
            [
                format!("{TP53}\t0\t{TRP53}"),
                format!("{TP53}\t-3\t{TRP53}"),
                format!("{TP53}\t1000000\t{TRP53}")
            ]
        );
    }

    #[test]
    fn gene_sets_are_sent_as_gene_set_lines() {
        let (_, stdin) = map_residues_request(
            &[(TP53, 15, TRP53)],
            &MapResiduesOptions {
                gtfs: vec![("homo_sapiens".into(), "human.gtf.gz".into())],
                gene_sets: vec![("mus_musculus".into(), "mouse.genes.tsv".into())],
                ..map_options()
            },
        )
        .unwrap();
        let stdin = stdin.unwrap();
        let species: Vec<&str> = stdin.split("--\n").next().unwrap().lines().collect();
        assert_eq!(
            species,
            [
                "homo_sapiens\tgtf\thuman.gtf.gz",
                "mus_musculus\tgene-set\tmouse.genes.tsv"
            ]
        );
    }

    #[test]
    fn malformed_map_residues_arguments_are_usage_errors_before_the_bridge_runs() {
        type Query<'a> = (&'a str, i64, &'a str);
        let none: [Query; 0] = [];
        let cases: Vec<(Vec<Query>, MapResiduesOptions, &str)> = vec![
            (none.to_vec(), map_options(), "queries is empty"),
            (
                vec![("", 15, TRP53)],
                map_options(),
                "source must be a protein id",
            ),
            (
                vec![(TP53, 15, " ")],
                map_options(),
                "target must be a protein id",
            ),
            (
                vec![(TP53, 15, "ENSMUSP\t1")],
                map_options(),
                "contains a tab or line break",
            ),
            (
                vec![(TP53, 15, TRP53)],
                MapResiduesOptions {
                    alignment: PathBuf::new(),
                    ..map_options()
                },
                "alignment is required",
            ),
            (
                vec![(TP53, 15, TRP53)],
                MapResiduesOptions {
                    gene_trees: PathBuf::new(),
                    ..map_options()
                },
                "gene_trees is required",
            ),
            (
                vec![(TP53, 15, TRP53)],
                MapResiduesOptions {
                    gene_sets: vec![("homo_sapiens".into(), "h.tsv".into())],
                    ..map_options()
                },
                "'homo_sapiens' is given twice in gtfs and gene_sets",
            ),
        ];
        for (queries, options, needle) in cases {
            match map_residues_request(&queries, &options) {
                Err(MzLibError::Usage(message)) => assert!(message.contains(needle), "{message}"),
                other => panic!("{needle}: {other:?}"),
            }
        }
    }

    #[test]
    fn every_outcome_is_mzlibs_name_and_a_gap_is_never_the_nearest_residue() {
        let mapped: ResidueMappings = recorded("genes_map_residues_tp53.json");
        assert_eq!(mapped.row_count, 9);
        assert_eq!(
            mapped.columns.names(),
            [
                "source_protein_id",
                "position",
                "target_protein_id",
                "outcome",
                "column",
                "target_position",
                "target_residue"
            ]
        );
        let rows = mapped.rows().unwrap();
        assert_eq!(
            rows.iter().map(|r| r.outcome.as_str()).collect::<Vec<_>>(),
            [
                "Aligned",
                "Aligned",
                "GapInTarget",
                "Aligned",
                "Aligned",
                "Aligned",
                "DifferentAlignments",
                "SourceNotInAlignment",
                "TargetNotInAlignment"
            ]
        );
        let mut keys: Vec<&str> = mapped.outcome_counts.keys().map(String::as_str).collect();
        let mut outcomes = COLUMN_OUTCOMES.to_vec();
        keys.sort_unstable();
        outcomes.sort_unstable();
        assert_eq!(keys, outcomes, "every outcome");
        assert_eq!(
            COLUMN_OUTCOMES.map(|o| mapped.outcome_counts[o]),
            [5, 1, 1, 1, 1]
        );
        // Human S15, S20 and S392 are mouse S18, S23 and S389; S46 faces a gap in mouse.
        assert_eq!(
            rows.iter().map(|r| r.target_position).collect::<Vec<_>>(),
            [
                Some(18),
                Some(23),
                None,
                Some(389),
                Some(15),
                Some(15),
                None,
                None,
                None
            ]
        );
        assert_eq!(
            rows.iter().map(|r| r.column).collect::<Vec<_>>(),
            [
                Some(476),
                Some(547),
                Some(863),
                Some(3036),
                Some(476),
                Some(476),
                None,
                None,
                None
            ]
        );
        assert_eq!(
            rows[2],
            ResidueMapping {
                source_protein_id: TP53.into(),
                position: 46,
                target_protein_id: TRP53.into(),
                outcome: "GapInTarget".into(),
                column: Some(863),
                target_position: None,
                target_residue: None,
            }
        );
    }

    #[test]
    fn the_map_residues_envelope_says_what_went_in() {
        let mapped: ResidueMappings = recorded("genes_map_residues_tp53.json");
        let a = &mapped.alignment;
        assert_eq!(
            (
                a.source_file_name.as_str(),
                a.release.as_deref(),
                a.collection.as_deref()
            ),
            (
                "Compara.116.protein_default.aa.fasta.gz",
                Some("116"),
                Some("protein_default")
            )
        );
        assert_eq!((a.alignment_count, a.protein_count), (2, 5));
        assert_eq!(a.restricted_to_gene_tree_content, mapped.gene_trees.sha256);
        assert_eq!(mapped.gene_trees.member_count, 5);
        let mut set_hashes: Vec<String> = mapped
            .gene_sets
            .values()
            .map(|g| g.sha256.clone())
            .collect();
        set_hashes.sort();
        assert_eq!(mapped.gene_trees.restricted_to_gene_sets, Some(set_hashes));
        assert_eq!(
            mapped.gene_sets["homo_sapiens"].genome_build.as_deref(),
            Some("GRCh38.p14")
        );
        assert_eq!(
            mapped.gene_sets.keys().collect::<Vec<_>>(),
            ["homo_sapiens", "mus_musculus", "rattus_norvegicus"]
        );
        assert!(mapped
            .caveats
            .iter()
            .any(|c| c.contains("2 queries name a protein the alignment did not keep")));
    }

    #[test]
    fn proteins_hold_the_kept_canonical_proteins_with_their_ungapped_sequences() {
        let mapped: ResidueMappings = recorded("genes_map_residues_tp53.json");
        let kept = mapped.proteins.as_ref().unwrap();
        assert_eq!(kept.row_count, 5);
        let ids = kept.columns.strings("protein_id").unwrap();
        let mut sorted = ids.clone();
        sorted.sort();
        assert_eq!(ids, sorted, "mzLib's ordinal order");
        let sequences = kept.sequences().unwrap();
        assert_eq!((sequences[TP53].len(), sequences[TRP53].len()), (393, 390));
        assert!(sequences.values().all(|s| !s.contains('-')));
        assert_eq!(
            (&sequences[TP53][14..15], &sequences[TRP53][17..18]),
            ("S", "S")
        );
        let trees: BTreeMap<String, String> = ids
            .into_iter()
            .zip(kept.columns.strings("tree_id").unwrap())
            .filter_map(|(id, tree)| Some((id?, tree?)))
            .collect();
        assert_eq!(trees[TP53], "ENSGT00950000183153");
        assert_eq!(trees[TRP53], trees[TP53]);
        assert!(
            !sequences.contains_key("ENSP00000352610"),
            "a TP53 isoform Compara did not align"
        );
    }

    #[test]
    fn without_proteins_the_block_is_none() {
        let site: ResidueMappings = recorded("genes_map_residues_chain.json");
        assert!(!site.proteins_included && site.proteins.is_none());
        assert_eq!(
            site.rows().unwrap(),
            [ResidueMapping {
                source_protein_id: TP53.into(),
                position: 15,
                target_protein_id: TRP53.into(),
                outcome: "Aligned".into(),
                column: Some(476),
                target_position: Some(18),
                target_residue: Some("S".into()),
            }]
        );
        assert!(site.caveats.is_empty());
    }
}
