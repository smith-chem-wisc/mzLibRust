//! SDRF-Proteomics experimental-design files: read one, pool several, and ask mzLib whether
//! they are well-formed, consistent and informative.
//!
//! Every other reader in this crate answers *what did the search find*. SDRF answers *what was
//! searched* — which sample, which organism part, which replicate, which instrument settings — and
//! that is the half you need to group results across experiments.
//!
//! | You want to | Call | mzLib type |
//! |---|---|---|
//! | Read what the file says, cell for cell | [`read`], [`read_with`] | `SdrfDocument` |
//! | Put several experiments in one table | [`pool_labelled`], [`pool`], [`pool_with`] | `SdrfCollection.Merge` |
//! | Know whether a file is well-formed | [`validate`], [`validate_many`] | `SdrfValidator` |
//! | Know whether several files write the same thing the same way | [`lint_labelled`], [`lint`] | `SdrfDriftLint` |
//! | Know whether a file describes its samples at all | [`assess`], [`assess_many`] | `SdrfSampleInformativeness` |
//! | Get what it says about each sample | [`samples`], [`samples_many`] | `SdrfSampleBlock` |
//! | Read ages in years | [`parse_ages`] | `SdrfAge.TryParse` |
//! | Drive a quantification from it, or learn why it cannot | [`design_with`], [`design`] | `SdrfLabelFreeDesign` |
//! | The same for a TMT, iTRAQ or DiLeu experiment: which sample is in which channel | [`isobaric_design_with`] | `SdrfIsobaricDesign` |
//! | Know which rows describe the files a search read | [`restrict_with`], [`restrict`] | `SdrfSearchScope` |
//!
//! ```
//! # mzlib_replay::activate();
//! let doc = mzlib::sdrf::read("PXD000070.sdrf.tsv")?;
//! assert_eq!((doc.row_count, doc.columns.len()), (6, 31));
//! let organisms = doc.value("characteristics[organism]");              // Vec<Option<&str>>
//! assert_eq!(organisms[0], Some("plasmodium falciparum"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Pool several experiments into one analysis table, naming each one yourself:
//!
//! ```
//! # mzlib_replay::activate();
//! use std::path::PathBuf;
//! use mzlib::sdrf::{pool_with, PoolInput, PoolOptions, ReadOptions};
//!
//! let documents = PoolInput::Labelled(vec![
//!     (PathBuf::from("PXD000070.sdrf.tsv"), "malaria".to_owned()),
//!     (PathBuf::from("PXD026824.sdrf.tsv"), "colon".to_owned()),
//! ]);
//! let options = PoolOptions {
//!     read: ReadOptions { limit: Some(4), ..Default::default() },
//!     ..Default::default()
//! };
//! let pooled = pool_with(&documents, &options)?;
//! assert_eq!((pooled.document_count, pooled.document.row_count), (2, 24));
//! assert_eq!(pooled.labels, ["malaria", "colon"]);
//! assert!(pooled.document.truncated);                 // 4 of the 24 rows came back
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Use this, not [`crate::readers::read_records`], for SDRF.** `read_records` recognises SDRF,
//! but it joins each row's cells into one semicolon-separated string, and SDRF's own key=value
//! grammar (`NT=Oxidation;AC=UNIMOD:35`) puts semicolons *inside* cells. The joined string cannot
//! be split back apart.
//!
//! **This module is row-major, and every other reader here is columnar.** That is not a style
//! choice. [`crate::readers::Table`] is keyed by column name because those names are a schema.
//! SDRF's are *data*, and they **repeat**: 649 files in the curated corpus carry
//! `comment[modification parameters]` more than once, up to eight times in one file, and one file
//! repeats an *empty* name 23 times. A map keyed by name would silently keep one occurrence and
//! drop the rest. So [`SdrfDocument::columns`] is a `Vec` that may contain duplicates,
//! [`SdrfDocument::rows`] is a `Vec` of cell `Vec`s, and position is what links them. Use
//! [`SdrfDocument::value`] for the first cell under a name and [`SdrfDocument::all`] for every one.
//!
//! **Three things worth knowing before you index anything:**
//!
//! * *Rows are ragged.* A row may be shorter than `columns`. Real files are like this — PXD059974
//!   in mzLib's own fixtures has a 46-column header with 17 of its 22 rows carrying 42 cells — and
//!   mzLib preserves it rather than padding, so the file round-trips byte for byte.
//!   [`SdrfDocument::value`] returns `None` for a position a row does not reach.
//! * *Cells are raw strings, never interpreted.* The key=value grammar arrives exactly as written.
//!   It is not decoded, because it cannot be told apart from a cell that merely contains `=` and
//!   `;` — `comment[file uri]` routinely carries pre-signed download URLs whose query strings
//!   contain `Signature=` and `Expires=`.
//! * *A reserved word is a real value.* `"not available"` and `"not applicable"` mean *the
//!   experiment stated an absence*, which is not the same as a column the document does not have.
//!   `None` means the latter. Do not collapse the two.
//!
//! ## Three questions about a document, and one about ages
//!
//! Each is mzLib's own answer, projected once in the bridge for all three bindings, never
//! reimplemented here — the mzLib types are in the table at the top of this page.
//!
//! **The first three are blind in different places, which is why there are three.** A file of
//! `"not available"` validates cleanly and lints clean, and only [`assess`] sees that it says
//! nothing. Twenty individually valid files can still disagree, and only [`lint_labelled`] looks
//! across them. Each result's `caveats` names the blind spot it leaves.
//!
//! ```
//! # mzlib_replay::activate();
//! let cohort = "sdrf_cohort.sdrf.tsv";
//! assert!(mzlib::sdrf::validate(cohort)?.is_valid);                 // well-formed...
//! assert_eq!(mzlib::sdrf::assess(cohort)?.verdict, "Informative");  // ...and says something
//! let drift = mzlib::sdrf::lint_labelled(&[
//!     (cohort, "cohort"),
//!     ("sdrf_cohort_partner.sdrf.tsv", "partner"),
//! ])?;
//! assert_eq!(drift.finding_count, 4);                                // ...but not like its partner
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! [`assess`] asks whether a file describes its samples at all. PXD000070's sample columns are
//! filled and its replicates differ, but its factor value never varies, so it is `Partial`; a
//! file whose characteristics are all blank is a `Skeleton`:
//!
//! ```
//! # mzlib_replay::activate();
//! let a = mzlib::sdrf::assess("PXD000070.sdrf.tsv")?;
//! assert_eq!(a.verdict, "Partial");
//! assert_eq!(
//!     (a.factor_value_varies, a.sample_is_described, a.biological_replicate_varies),
//!     (false, true, true)
//! );
//! assert_eq!(mzlib::sdrf::assess("sdrf_skeleton.sdrf.tsv")?.verdict, "Skeleton");
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! [`validate_many`], [`assess_many`] and [`samples_many`] read a whole corpus in **one** bridge
//! call, [`BulkOptions::threads`] documents at a time, and return one long table whose first two
//! columns say which document each row came from. The answer is identical at any thread count;
//! there is deliberately no Rust-side loop, because parallelism is the bridge's to express once.
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::sdrf::{assess_many, BulkOptions};
//!
//! let corpus = ["sdrf_cohort.sdrf.tsv", "sdrf_skeleton.sdrf.tsv", "PXD000070.sdrf.tsv"];
//! let batch = assess_many(&corpus, &BulkOptions::default())?;
//! assert_eq!(batch.paths_with(&["Skeleton"])?, ["sdrf_skeleton.sdrf.tsv"]);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Ages are read into **years** — a month is 1/12, a week 7/365.25, a day 1/365.25 — and a cell
//! that would need a guess is refused with its reason, never assumed. A bare `63` is refused,
//! because 63 years and 63 days are both plausible in one study.
//!
//! ## Turn an SDRF into a quantification design
//!
//! FlashLFQ and MetaMorpheus need an experimental design: which runs are the same condition,
//! which are biological replicates, which are fractions of one sample. The SDRF already says all
//! of that, and an invalid hand-written `ExperimentalDesign.tsv` is worse than none: MetaMorpheus
//! skips quantification with one warning and no error. [`design_with`] reads the design out of the
//! SDRF with mzLib's `SdrfLabelFreeDesign` (mzLib #1363), or refuses and lists **every** reason at
//! once, so you fix the file once rather than once per failure.
//!
//! PXD067622, a TurboID proximity-labelling study — two genotypes by four treatments, three
//! biological replicates each. Name the factor columns that make up the condition; their values
//! are joined with `_`, in that order:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::sdrf::{design_with, DesignOptions};
//!
//! let both = DesignOptions {
//!     condition_columns: vec!["factor value[genotype]".into(), "factor value[treatment]".into()],
//!     ..Default::default()
//! };
//! let d = design_with("PXD067622.sdrf.tsv", &both)?;
//! assert_eq!((d.is_valid, d.file_count), (true, 24));
//! for run in &d.files()?[..3] {
//!     println!("{} | {} | {}", run.file_name, run.condition, run.biological_replicate);
//! }
//! assert_eq!(d.files()?[2].biological_replicate, 2);  // 0-based: the third replicate
//!
//! let spectra = d.spectra()?;                          // what flashlfq::quantify_with takes
//! assert_eq!(spectra[0].condition.as_deref(), Some("SPRTN-TurboID CA_DMSO (vehicle)"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **The design is 0-based; `ExperimentalDesign.tsv` is 1-based.** `biological_replicate`,
//! `technical_replicate` and `fraction` are mzLib's `SpectraFileInfo` coordinates, exactly what
//! [`crate::flashlfq::quantify_with`] and [`crate::flashlfq::median_polish_with`] take.
//! [`DesignOptions::out`] writes MetaMorpheus's file, and mzLib adds the one when it writes and
//! nowhere else. Never add it yourself.
//!
//! **A refusal is a result, not an error.** PXD049018's depositors did not say which pulldown was
//! treated: its `factor value[treatment]` is `not available` on every row, and a condition built
//! from an unknown factor would pair samples nobody said were alike:
//!
//! ```
//! # mzlib_replay::activate();
//! # use mzlib::sdrf::{design_with, DesignOptions};
//! # let both = DesignOptions {
//! #     condition_columns: vec!["factor value[genotype]".into(), "factor value[treatment]".into()],
//! #     ..Default::default()
//! # };
//! let refused = design_with("PXD049018.sdrf.tsv", &both)?;
//! assert_eq!((refused.is_valid, refused.refusals.len()), (false, 20)); // one per run
//! println!("{}", refused.report);                                       // mzLib's whole account
//! assert!(refused.spectra().is_err());       // so it cannot reach FlashLFQ by accident
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Replicate numbers are ranked within each condition**, and [`SdrfDesign::notes`] is the only
//! record of the renumbering. A drafted SDRF often copies the study-wide run index from the file
//! names (WT_DMSO1-3, CA_DMSO4-6, … CA_FA22-24); mzLib ranks them back, and the result is the
//! hand-numbered design exactly:
//!
//! ```
//! # mzlib_replay::activate();
//! # use mzlib::sdrf::{design_with, DesignOptions};
//! # let both = DesignOptions {
//! #     condition_columns: vec!["factor value[genotype]".into(), "factor value[treatment]".into()],
//! #     ..Default::default()
//! # };
//! let ranked = design_with("PXD067622_studywide.sdrf.tsv", &both)?;
//! assert_eq!(ranked.notes.len(), 7);
//! assert!(ranked.notes[1].ends_with("biological replicates renumbered 22 -> 1, 23 -> 2, 24 -> 3."));
//! assert_eq!(ranked.columns, design_with("PXD067622.sdrf.tsv", &both)?.columns);
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Without condition columns, mzLib uses the document's only factor column; PXD067622 has two, so
//! leaving them undeclared is itself a refusal. When the SDRF names files without the directory
//! you keep them in, or the search read only some of them, pass [`DesignOptions::searched_files`]:
//! each must be named exactly (case and extension) by one row, and rows for other files are
//! dropped and listed in `notes`. **Label-free only**: an isobaric SDRF needs a channel design;
//! read it with [`isobaric_design_with`], below. A search that read `X-calib.mzML` matches no row
//! naming `X.raw`: [`restrict_with`] the SDRF to the search first.
//!
//! ## Read an isobaric (TMT) design
//!
//! **The question.** In a TMT, TMTpro, iTRAQ or DiLeu experiment, many samples share one run:
//! each sits in a channel, and its reporter ion says how much of each peptide came from it. To
//! compare samples, MetaMorpheus needs to know which sample is in which channel of which file, and
//! it reads that only from a `TmtDesign.txt` beside the spectra. Without one, it writes reporter
//! intensities per PSM and quantifies nothing per channel. The SDRF already holds the design: one
//! row per file and channel.
//!
//! **[`isobaric_design_with`] reads it out, or refuses and says why.** It calls mzLib's
//! `SdrfIsobaricDesign` (mzLib #1380), which runs MetaMorpheus's own checks on a TMT design and
//! reports every failure at once: one plex, fraction and technical replicate per file; each channel
//! described one way within a plex; no sample twice at one fraction of a plex. It needs two facts
//! that mzLib will not guess from the SDRF:
//!
//! | You must say | How | Why mzLib will not guess it |
//! |---|---|---|
//! | The kit the search used | [`IsobaricDesignOptions::kit`] | One kit's labels can be a subset of another's: a TMT11 plex that left 131C empty reads exactly like TMT10 |
//! | Where each file's plex is written | exactly one of `plex_column`, `plex_file_name_pattern`, `single_plex` | The column that looks like a plex often is not one (see PXD061609 below) |
//!
//! The kit is named as [`crate::isobaric::kits`] takes it, by the same mzLib code: `"TMT10"`,
//! `"tmt18"`, `"TMT6-plex"`, `"iTRAQ4"`, `"TMT11 on K"`. A name mzLib does not know is
//! [`MzLibError::Usage`], listing the known kits.
//!
//! The examples use mzLib's own isobaric design files, cut from
//! [bigbio/sdrf-annotated-datasets](https://github.com/bigbio/sdrf-annotated-datasets): PXD008841
//! (breast tumours, TMT10; two fractions each of two plexes of a nine-plex study) and PXD061609
//! (eighteen mouse samples, ten tissues of both sexes, TMTpro 18; two fractions of each of three
//! fractionation arms of one plex).
//!
//! PXD008841 has no plex column. The plex survives only in the file names, spelled two ways
//! (`TMTpool1`, `TMT_pool5`), so one regular expression reads both. Its group named `plex` is the
//! plex:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::sdrf::{isobaric_design_with, IsobaricDesignOptions};
//!
//! let tmt = isobaric_design_with(
//!     "PXD008841.subset.sdrf.tsv",
//!     &IsobaricDesignOptions {
//!         kit: Some("TMT10".into()),
//!         plex_file_name_pattern: Some(r"TMT_?pool(?<plex>\d+)".into()),
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!((tmt.is_valid, tmt.file_count, tmt.record_count), (true, 4, 40));
//! assert_eq!(tmt.plex_ids, [("1".to_owned(), 1), ("5".to_owned(), 2)].into());
//! let channels = tmt.channels()?;
//! for c in &channels[..3] {
//!     println!("{} | {} | {} | {} | {}", c.plex, c.channel, c.sample_name, c.condition, c.fraction);
//! }
//! //   1 | 126 | OSL.53E | basal-like breast carcinoma | 1
//! //   1 | 127N | OSL.567 | luminal A breast carcinoma | 1
//! //   1 | 127C | OSL.3FF | basal-like breast carcinoma | 1
//! assert_eq!(channels[1].sample_name, "OSL.567");
//! assert_eq!((channels[9].channel.as_str(), channels[9].sample_name.as_str()), ("131N", "pool"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **What each row means.** One row per file and annotated channel: files in SDRF order, each
//! file's channels in reporter m/z order ([`SdrfIsobaricDesign::kit_channels`]).
//!
//! | Field | Meaning |
//! |---|---|
//! | `full_path` | The file, as the SDRF names it, or the searched path you passed for it |
//! | `plex`, `plex_id` | The plex, as your source gives it, and its number, 1..N in name order (MetaMorpheus numbers plexes the same way) |
//! | `sample_name`, `channel` | The SDRF `source name` in that channel, and the channel as the kit spells it (`127N`), whatever form `comment[label]` used |
//! | `condition` | The condition columns' values joined with `_`, exactly as [`design_with`] builds it |
//! | `biological_replicate`, `fraction`, `technical_replicate` | The SDRF's own numbers, **1-based**, never ranked or shifted |
//!
//! **This design is 1-based; [`design_with`] is 0-based.** [`isobaric_design_with`] returns the
//! SDRF's numbers as written, starting at 1, which is also what `TmtDesign.txt` holds.
//! [`design_with`] returns mzLib's `SpectraFileInfo` coordinates, starting at 0. The split is
//! mzLib's own, and neither the bridge nor this crate converts either table. **Never join or
//! compare the two without adding or removing 1.**
//!
//! Channel 131N is the reference pool of PXD008841's bridge design. mzLib does not read sample
//! type yet, so every annotated channel is a study sample, the pool included: leave it out of a
//! comparison yourself. A kit channel the SDRF does not annotate has no row. Read PXD008841 as
//! TMT11 and `kit_channels` gains `131C`; the table does not.
//!
//! ### Nothing is guessed
//!
//! Declare no plex source, and mzLib refuses rather than infer one:
//!
//! ```
//! # mzlib_replay::activate();
//! # use mzlib::sdrf::{isobaric_design_with, IsobaricDesignOptions};
//! let unknown = isobaric_design_with(
//!     "PXD008841.subset.sdrf.tsv",
//!     &IsobaricDesignOptions { kit: Some("TMT10".into()), ..Default::default() },
//! )?;
//! assert_eq!((unknown.is_valid, unknown.plex_source.as_str()), (false, "(not declared)"));
//! assert!(unknown.refusals[0].starts_with(
//!     "No plex source was declared. Name the column that holds the plex, give a file-name \
//!      pattern whose group captures it, or name the single plex."
//! ));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! PXD061609 is the reason. Its `comment[sample preparation batch]` holds `NoSAX`, `HighSalt` and
//! `LowSalt`: three fractionation arms of **one** plex, every arm carrying all eighteen samples
//! and numbering its fractions from 1. Read as the one plex it is, the arms' fraction 1s collide,
//! and MetaMorpheus could not tell them apart, so mzLib refuses, once for each sample and
//! fraction. Declare the batch column as the plex, and each arm becomes a plex of its own, which
//! MetaMorpheus accepts:
//!
//! ```
//! # mzlib_replay::activate();
//! # use mzlib::sdrf::{isobaric_design_with, IsobaricDesignOptions};
//! let one = isobaric_design_with(
//!     "PXD061609.subset.sdrf.tsv",
//!     &IsobaricDesignOptions {
//!         kit: Some("TMT18".into()),
//!         single_plex: Some("tube1".into()),
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!((one.is_valid, one.refusals.len()), (false, 36));
//! assert!(one.refusals[0].starts_with(
//!     "Plex 'tube1' has sample 'brain_female' biorep 1 fraction 1 techrep 1 more than once:"
//! ));
//!
//! let arms = isobaric_design_with(
//!     "PXD061609.subset.sdrf.tsv",
//!     &IsobaricDesignOptions {
//!         kit: Some("TMT18".into()),
//!         plex_column: Some("comment[sample preparation batch]".into()),
//!         ..Default::default()
//!     },
//! )?;
//! assert!(arms.is_valid);
//! assert_eq!(arms.plex_ids.keys().collect::<Vec<_>>(), ["HighSalt", "LowSalt", "NoSAX"]);
//! assert_eq!(arms.record_count, 108);
//! let second = &arms.channels()?[1];
//! assert_eq!(
//!     (second.sample_name.as_str(), second.biological_replicate, second.condition.as_str()),
//!     ("brain_male", 2, "NT=brain;AC=UBERON:0000955")
//! );
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Whether a fractionation arm should be quantified as a plex is a decision about your
//! experiment, not something mzLib can read from the file; that is why it makes you declare it.
//! Note two more things about the rows. `condition` is the factor cell as written, key=value
//! grammar and all. And the biological replicate is the SDRF's own number (here, a tissue's female
//! and male samples are replicates 1 and 2): an isobaric design never ranks replicates the way
//! [`design_with`] does.
//!
//! **The traps.**
//!
//! - **1-based against 0-based**, above. It is the one to remember.
//! - **A refusal is a result.** `is_valid` `false` returns `Ok`, with `refusals` filled and the
//!   table empty. A missing kit, no plex source, two plex sources and a pattern that does not parse
//!   are all refusals, not errors. Print `report` for mzLib's whole account.
//! - **Searched files match by exact name.** [`IsobaricDesignOptions::searched_files`] works as
//!   for [`design_with`]: each must be named exactly by some row. A search that read `X-calib.mzML`
//!   matches no row naming `X.raw`; restrict the SDRF to the search first (next section).
//! - **`out` names a folder.** mzLib writes MetaMorpheus's `TmtDesign.txt` inside it (only for a
//!   valid design), replacing one already there. MetaMorpheus reads it only from the folder of the
//!   first spectra file, and refuses a folder that also holds an `ExperimentalDesign.tsv`.
//! - **TMT16 is TMTpro 16-plex.** MetaMorpheus searches it with its TMT18 label; a TMT16 design is
//!   read as is by such a run, with 134C and 135N empty.
//!
//! ## Restrict an SDRF to one search's files
//!
//! **The question.** A deposited SDRF describes every file the depositors acquired. Your search
//! read some of them, often under other names: MetaMorpheus calibrates `X.raw` and searches
//! `X-calib.mzML`. Which rows describe the files your search read, and what does the SDRF fail to
//! describe?
//!
//! **[`restrict_with`] answers both.** It calls mzLib's `SdrfSearchScope.Restrict` (mzLib #1451),
//! which joins each searched file to its acquisition **by stem**, ignoring case, after removing
//! the derivative suffixes MetaMorpheus writes (`-calib`, `-averaged`, in any order). Each kept row
//! is carried cell for cell, with the searched file's name written into a new
//! [`SEARCHED_DATA_FILE_COLUMN`], directly after `comment[data file]`. Nothing is renumbered.
//!
//! Here is a search over PXD008841 that hits every problem but one. It has two calibrated runs,
//! the raw file of one of them as well, a renamed file whose acquisition you know, and a file the
//! SDRF does not describe:
//!
//! ```
//! # mzlib_replay::activate();
//! use mzlib::sdrf::{restrict_with, RestrictOptions, SearchedFile};
//!
//! let pool1 = "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0";
//! let pool5 = "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0";
//! let scoped = restrict_with(
//!     "PXD008841.subset.sdrf.tsv",
//!     &[
//!         SearchedFile::from(format!(r"D:\search\{pool1}1-calib.mzML")),
//!         SearchedFile::from(format!("{pool1}2-calib.mzML")),
//!         SearchedFile::from(format!("{pool1}2.raw")),
//!         SearchedFile::acquired_from("sample_A.mzML", format!("{pool5}1.raw")),
//!         SearchedFile::from("not_in_sdrf.mzML"),
//!     ],
//!     &RestrictOptions { limit: Some(2), ..Default::default() },
//! )?;
//! assert_eq!((scoped.input_row_count, scoped.document.row_count), (40, 30));
//! assert_eq!(
//!     scoped.document.value("comment[data file]")[0],
//!     Some("HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01.raw")
//! );
//! assert_eq!(
//!     scoped.searched_data_files()[0],
//!     Some("HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01-calib.mzML")
//! );
//! assert_eq!(scoped.searched_without_row, ["not_in_sdrf.mzML"]);
//! assert_eq!(
//!     scoped.dropped_data_files,
//!     ["HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr02.raw"]
//! );
//! // Two searched files claimed fraction 2's acquisition; the first in name order kept its rows.
//! assert_eq!(
//!     scoped.ambiguous,
//!     [format!("{pool1}2: {pool1}2-calib.mzML, {pool1}2.raw")]
//! );
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Pass [`SearchedFile::from`] a path or name when you know only the searched name, and
//! [`SearchedFile::acquired_from`] when you also know which acquired file a searched one came from
//! (MetaMorpheus reads it from an mzML's own `sourceFile`). A name you give wins over the suffix
//! rule, which is how `sample_A.mzML` found pool 5's first fraction above.
//!
//! **What each report means.**
//!
//! | Report | Meaning | What to do |
//! |---|---|---|
//! | `searched_without_row` | Searched files the SDRF has no row for | They have no design. Fix the SDRF, or leave them out of the quantification |
//! | `dropped_data_files` | Acquired files the SDRF describes that the search did not read | Usually expected; check that none was meant to be searched |
//! | `ambiguous` | Two or more searched files claimed one acquisition, as `stem: searched, searched` | The first in name order kept the rows. The others have none, and are not in `searched_without_row` |
//! | `shared_stem` | One searched file joined rows naming two **different** acquired files (`X.raw` and `X.mzML`), as `searched: acquired, acquired` | Every such row is kept, so the document describes one run twice. Only the depositor knows which row is right |
//!
//! [`RestrictedSdrf::document`] is an [`SdrfDocument`]: [`SdrfDocument::value`],
//! [`SdrfDocument::all`] and [`SdrfDocument::records`] work on it. A row shorter than the header is
//! padded with empty cells, so unlike [`read`], its rows are never ragged.
//!
//! ### Restrict, then read the design
//!
//! [`design_with`] and [`isobaric_design_with`] match `searched_files` by **exact** name. That is
//! deliberate: mzLib leaves the stem join to [`restrict_with`]. Both readers key on
//! `comment[searched data file]` when a document has it, so restrict first, write the result, and
//! read the design from it:
//!
//! ```
//! # mzlib_replay::activate();
//! use std::path::PathBuf;
//!
//! use mzlib::sdrf::{
//!     isobaric_design_with, restrict_with, IsobaricDesignOptions, RestrictOptions, SearchedFile,
//!     WrittenSdrf,
//! };
//!
//! let runs: Vec<PathBuf> = [
//!     "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0",
//!     "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0",
//! ]
//! .iter()
//! .flat_map(|pool| (1..=2).map(move |i| PathBuf::from(format!(r"D:\search\{pool}{i}-calib.mzML"))))
//! .collect();
//! let searched: Vec<SearchedFile> = runs.iter().cloned().map(SearchedFile::from).collect();
//!
//! let kept = restrict_with(
//!     "PXD008841.subset.sdrf.tsv",
//!     &searched,
//!     &RestrictOptions {
//!         out: Some("PXD008841_searched.sdrf.tsv".into()),
//!         limit: Some(2),
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!(
//!     kept.written,
//!     Some(WrittenSdrf { path: "PXD008841_searched.sdrf.tsv".to_owned(), row_count: 40 })
//! );
//!
//! let tmt = isobaric_design_with(
//!     "PXD008841_searched.sdrf.tsv",
//!     &IsobaricDesignOptions {
//!         kit: Some("TMT10".into()),
//!         plex_file_name_pattern: Some(r"TMT_?pool(?<plex>\d+)".into()),
//!         searched_files: Some(runs),
//!         out: Some("design".into()),
//!         ..Default::default()
//!     },
//! )?;
//! assert_eq!(tmt.file_key_column.as_deref(), Some("comment[searched data file]"));
//! // Recorded on Windows. Read a path from these as text: a reader on another platform does not
//! // split it on '\'.
//! assert_eq!(
//!     tmt.channels()?[0].full_path,
//!     r"D:\search\HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01-calib.mzML"
//! );
//! let written = tmt.written.as_ref().expect("a valid design with out writes TmtDesign.txt");
//! assert_eq!(written.record_count, 40);
//! assert!(written.path.ends_with("TmtDesign.txt"));
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! Give the same searched files to [`isobaric_design_with`] with the unrestricted SDRF and every
//! one is refused: the SDRF names `X.raw`, and the search read `X-calib.mzML`.
//!
//! **The traps.**
//!
//! - **`out` must end in `.tsv`; name it `.sdrf.tsv`.** mzLib adds `.sdrf.tsv` to a name that does
//!   not already end in it, so `searched.tsv` is written as `searched.tsv.sdrf.tsv`.
//!   [`WrittenSdrf::path`] always names the file actually written, for [`pool_with`] as well.
//! - **`out` writes the whole document**, whatever the limit and offset are, with CRLF line
//!   endings whatever the input used.
//! - **`comment[searched data file]` holds the file name, not the folder.** The design readers
//!   take the folder from their `searched_files`.
//! - **mzLib's design readers split a path on this platform's separator only** (the bridge's
//!   thread, SDRF-B5). `restrict_with` reduces either separator on any platform, but a design read
//!   on Linux from `D:\search\X-calib.mzML` keys on the whole string. Pass searched paths in the
//!   form the platform you run on writes them.
//! - **An SDRF with no `comment[data file]` column cannot be restricted.** mzLib raises, and this
//!   crate returns [`MzLibError::Bridge`] with `error_type` `"ArgumentException"`.
//!
//! Ported from pyMzLib's `pymzlib.sdrf`, which decided the verbs, the wire fields and the caveats.
//!
//! ## Cite
//!
#![doc = include_str!("../docs/reference/cite.sdrf.md")]

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};
use crate::readers::Table;

pub use crate::bridge::OnError;

/// The column [`pool`] adds to record which document each row came from.
pub const SOURCE_DOCUMENT_COLUMN: &str = "comment[source document]";

/// The default timeout, matching pyMzLib's `sdrf` module. SDRF files are small.
const DEFAULT_TIMEOUT: Duration = Duration::from_secs(60);

// ---------------------------------------------------------------------------------------------
// The document
// ---------------------------------------------------------------------------------------------

/// One SDRF-Proteomics document: an ordered header, and rows of raw cells.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct SdrfDocument {
    /// The path that was read. Empty for a pooled table, which is not a file that was read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// The column names, **verbatim and in document order**. Names may repeat, and the order is
    /// part of the document. Names are never case-normalised: the corpus contains
    /// `"comment[MS min charge]"` and `"Material Type"`, and rewriting them would silently alter
    /// someone's file.
    #[serde(
        rename = "column_names",
        default,
        deserialize_with = "bridge::null_to_default"
    )]
    pub columns: Vec<String>,
    /// One `Vec` of cells per row. **Ragged**: a row may be shorter than [`Self::columns`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub rows: Vec<Vec<String>>,
    /// Data rows in the **whole document** (the header is not a row), regardless of any limit or
    /// offset.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Data rows actually carried back in [`Self::rows`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub returned_count: u64,
    /// The offset that was applied, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub offset: u64,
    /// **Whether rows were left behind**, by either the limit or the offset. A short answer and a
    /// complete one must never look alike, so check this rather than comparing counts yourself.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub truncated: bool,
    /// What this document's data cannot tell you about itself — raggedness, repeated names,
    /// reserved words. Worth printing the first time you read an unfamiliar file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfDocument {
    /// The position of the first column with this name, or `None` if the document has none.
    ///
    /// Comparison is exact and case-**sensitive**, matching the SDRF specification and mzLib.
    pub fn index_of(&self, column: &str) -> Option<usize> {
        self.columns.iter().position(|name| name == column)
    }

    /// Every position carrying this name, in document order. Empty when the column is absent.
    pub fn indexes_of(&self, column: &str) -> Vec<usize> {
        self.columns
            .iter()
            .enumerate()
            .filter(|(_, name)| *name == column)
            .map(|(i, _)| i)
            .collect()
    }

    /// The first cell under `column`, one entry per returned row.
    ///
    /// `None` means **the document does not have this column, or the row is too short to reach
    /// it**. It does not mean "empty": the SDRF reserved words `"not available"` and
    /// `"not applicable"` are real values that an experiment chose to write, and they come back as
    /// themselves.
    pub fn value(&self, column: &str) -> Vec<Option<&str>> {
        match self.index_of(column) {
            None => vec![None; self.rows.len()],
            Some(i) => self
                .rows
                .iter()
                .map(|row| row.get(i).map(String::as_str))
                .collect(),
        }
    }

    /// Every cell under `column`, one `Vec` per returned row.
    ///
    /// The accessor for a multi-cardinality column such as `comment[modification parameters]`,
    /// which legitimately repeats — up to eight times in one corpus file. Positions a row is too
    /// short to reach are skipped rather than reported, so each inner `Vec` holds only cells that
    /// exist.
    pub fn all(&self, column: &str) -> Vec<Vec<&str>> {
        let indexes = self.indexes_of(column);
        self.rows
            .iter()
            .map(|row| {
                indexes
                    .iter()
                    .filter_map(|&i| row.get(i).map(String::as_str))
                    .collect()
            })
            .collect()
    }

    /// The rows as name-to-cell maps, for the common case of a document with no repeated names.
    ///
    /// **Lossy when a name repeats** — a later position overwrites an earlier one — which is
    /// exactly why it is not the primary shape. [`Self::has_repeated_columns`] says whether that
    /// applies to this document; for a repeating document use [`Self::all`] instead. A position a
    /// row is too short to reach is absent from that row's map.
    pub fn records(&self) -> Vec<HashMap<&str, &str>> {
        self.rows
            .iter()
            .map(|row| {
                self.columns
                    .iter()
                    .zip(row)
                    .map(|(name, cell)| (name.as_str(), cell.as_str()))
                    .collect()
            })
            .collect()
    }

    /// Whether any column name appears more than once — see [`Self::records`].
    pub fn has_repeated_columns(&self) -> bool {
        let mut seen = std::collections::HashSet::with_capacity(self.columns.len());
        !self.columns.iter().all(|name| seen.insert(name))
    }

    /// How many returned rows carry fewer cells than there are columns.
    pub fn ragged_row_count(&self) -> usize {
        self.rows
            .iter()
            .filter(|row| row.len() < self.columns.len())
            .count()
    }
}

/// Where [`pool_with`] or [`restrict_with`] wrote its document, when asked to write one.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct WrittenSdrf {
    /// The file mzLib actually wrote, which is not always the name you gave: mzLib adds
    /// `.sdrf.tsv` to a name that does not already end in it, so an `out` of `x.tsv` writes
    /// `x.tsv.sdrf.tsv`. Bridges before pyMzLib 0.6.0 reported the name as given for [`pool_with`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Rows written — the **whole** document, not the windowed slice. The limit and offset shape
    /// what comes back over the wire; they never shorten the file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
}

/// Several SDRF documents merged into one table.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct PooledSdrf {
    /// The merged table. Its `path` is empty: a pooled table is not a file that was read.
    #[serde(flatten)]
    pub document: SdrfDocument,
    /// How many documents were pooled.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub document_count: u64,
    /// The paths pooled, in the order given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub paths: Vec<String>,
    /// The provenance label used for each, in the same order — either what you supplied or
    /// mzLib's `containing-folder/file-stem` default.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub labels: Vec<String>,
    /// Where the merged document was written, or `None` if [`PoolOptions::out`] was not set.
    #[serde(default)]
    pub written: Option<WrittenSdrf>,
}

impl PooledSdrf {
    /// The provenance label of each returned row — which document it came from.
    pub fn source_documents(&self) -> Vec<Option<&str>> {
        self.document.value(SOURCE_DOCUMENT_COLUMN)
    }
}

// ---------------------------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------------------------

/// How much of one document to return.
#[derive(Debug, Clone)]
pub struct ReadOptions {
    /// Return at most this many rows. `None` returns all of them. `Some(0)` returns the header
    /// alone.
    pub limit: Option<u64>,
    /// Skip this many rows.
    pub offset: u64,
    /// Time to allow. `None` waits indefinitely. Defaults to 60 seconds.
    pub timeout: Option<Duration>,
}

impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            limit: None,
            offset: 0,
            timeout: Some(DEFAULT_TIMEOUT),
        }
    }
}

/// How much of a pooled table to return, and where to write the whole of it.
#[derive(Debug, Clone, Default)]
pub struct PoolOptions {
    /// The window and the timeout, as for [`read_with`].
    pub read: ReadOptions,
    /// Write the merged document here as SDRF. The **whole** document is written regardless of
    /// the limit and offset. A name that does not end in `.sdrf.tsv` gets it added, so `x.tsv`
    /// writes `x.tsv.sdrf.tsv`; [`WrittenSdrf::path`] says which file was written.
    pub out: Option<String>,
}

/// The documents to pool: bare paths, or each path with the provenance label to stamp on its rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PoolInput {
    /// Paths only. mzLib falls back to `containing-folder/file-stem` for provenance, which depends
    /// on where the files happen to sit, so the result carries a caveat saying it is not
    /// reproducible elsewhere. Prefer [`PoolInput::Labelled`].
    Paths(Vec<PathBuf>),
    /// Each path with the label you chose for it. Every label must be non-blank.
    Labelled(Vec<(PathBuf, String)>),
}

// ---------------------------------------------------------------------------------------------
// Argument assembly
// ---------------------------------------------------------------------------------------------

fn path_text<'a>(path: &'a Path, what: &str) -> Result<&'a str> {
    let text = path.to_str().ok_or_else(|| {
        MzLibError::Usage(format!(
            "The {what} is not valid UTF-8, which the bridge requires."
        ))
    })?;
    Ok(text.trim())
}

fn window_args(args: &mut Vec<String>, options: &ReadOptions) {
    if let Some(limit) = options.limit {
        args.push("--limit".to_owned());
        args.push(limit.to_string());
    }
    if options.offset > 0 {
        args.push("--offset".to_owned());
        args.push(options.offset.to_string());
    }
}

fn read_args(path: &Path, options: &ReadOptions) -> Result<Vec<String>> {
    let path = path_text(path, "file path")?;
    if path.is_empty() {
        return Err(MzLibError::Usage(
            "A file path is required, e.g. 'PXD000070.sdrf.tsv'.".to_owned(),
        ));
    }
    let mut args = vec![
        "sdrf".to_owned(),
        "read".to_owned(),
        "--path".to_owned(),
        path.to_owned(),
    ];
    window_args(&mut args, options);
    Ok(args)
}

/// The arguments and the stdin for a pool, validated before anything is spawned.
///
/// Documents cross on stdin, one per line, as `path[\tlabel]` — the bridge keeps named options in
/// a dictionary, so a repeated `--path` would silently keep only the last.
fn pool_request(documents: &PoolInput, options: &PoolOptions) -> Result<(Vec<String>, String)> {
    let lines = document_lines(documents)?;

    let mut args = vec!["sdrf".to_owned(), "pool".to_owned()];
    window_args(&mut args, &options.read);
    if let Some(out) = &options.out {
        if out.trim().is_empty() {
            return Err(MzLibError::Usage(
                "out must be a non-empty path, or None to skip writing.".to_owned(),
            ));
        }
        args.push("--out".to_owned());
        args.push(out.trim().to_owned());
    }

    Ok((args, lines.join("\n")))
}

/// The `path[\tlabel]` stdin lines `sdrf pool` and `sdrf lint` share, validated.
fn document_lines(documents: &PoolInput) -> Result<Vec<String>> {
    let pairs: Vec<(&Path, Option<&str>)> = match documents {
        PoolInput::Paths(paths) => paths.iter().map(|p| (p.as_path(), None)).collect(),
        PoolInput::Labelled(pairs) => pairs
            .iter()
            .map(|(p, label)| (p.as_path(), Some(label.as_str())))
            .collect(),
    };
    if pairs.is_empty() {
        return Err(MzLibError::Usage(
            "At least one SDRF document is required.".to_owned(),
        ));
    }

    let mut lines = Vec::with_capacity(pairs.len());
    for (path, label) in pairs {
        let path = path_text(path, "document path")?;
        if path.is_empty() {
            return Err(MzLibError::Usage(
                "A document path may not be blank.".to_owned(),
            ));
        }
        let label = label.map(str::trim);
        if label.is_some_and(str::is_empty) {
            return Err(MzLibError::Usage(format!(
                "The label for '{path}' is blank. Give every document a label, or use \
                 PoolInput::Paths to accept mzLib's path-derived default for all of them."
            )));
        }
        let separator = |text: &str| text.contains(['\t', '\n', '\r']);
        if separator(path) || label.is_some_and(separator) {
            return Err(MzLibError::Usage(format!(
                "A path or label contains a tab or newline, which the bridge uses to separate \
                 them: {path:?} / {:?}.",
                label.unwrap_or("")
            )));
        }
        lines.push(match label {
            Some(label) => format!("{path}\t{label}"),
            None => path.to_owned(),
        });
    }
    Ok(lines)
}

// ---------------------------------------------------------------------------------------------
// The public surface
// ---------------------------------------------------------------------------------------------

/// Read one SDRF-Proteomics file, every row.
///
/// See [`read_with`] for the reference.
///
/// # Errors
///
/// As [`read_with`].
pub fn read(path: impl AsRef<Path>) -> Result<SdrfDocument> {
    read_with(path, &ReadOptions::default())
}

/// Read one SDRF-Proteomics file with every cell intact, in a row-major shape that keeps the
/// document's repeated column names.
///
/// A blank path is refused before anything is spawned.
#[doc = include_str!("../docs/reference/sdrf.read.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// // PXD059974: a 46-column header over rows of 42 cells. mzLib keeps the raggedness.
/// let doc = mzlib::sdrf::read("PXD059974.sdrf.tsv")?;
/// assert_eq!(doc.ragged_row_count(), 17);
/// let last = doc.columns.last().unwrap();
/// assert!(doc.value(last).iter().any(Option::is_none)); // a short row reaches no cell there
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.read.see-also.md")]
pub fn read_with(path: impl AsRef<Path>, options: &ReadOptions) -> Result<SdrfDocument> {
    let args = read_args(path.as_ref(), options)?;
    let data = bridge::invoke(&args, None, options.timeout)?;
    serde_json::from_value(data).map_err(protocol)
}

/// Merge several SDRF documents into one analysis table, with provenance from their paths.
///
/// Prefer [`pool_labelled`]: see [`PoolInput::Paths`] for why.
///
/// # Errors
///
/// As [`pool_with`].
pub fn pool<P: AsRef<Path>>(paths: &[P]) -> Result<PooledSdrf> {
    pool_with(
        &PoolInput::Paths(paths.iter().map(|p| p.as_ref().to_path_buf()).collect()),
        &PoolOptions::default(),
    )
}

/// Merge several SDRF documents into one analysis table, naming each one yourself.
///
/// # Errors
///
/// As [`pool_with`].
pub fn pool_labelled<P: AsRef<Path>, L: AsRef<str>>(documents: &[(P, L)]) -> Result<PooledSdrf> {
    pool_with(
        &PoolInput::Labelled(
            documents
                .iter()
                .map(|(p, l)| (p.as_ref().to_path_buf(), l.as_ref().to_owned()))
                .collect(),
        ),
        &PoolOptions::default(),
    )
}

/// Merge several SDRF documents into one analysis table, with provenance.
///
/// Columns are the union of every document's, ordered by SDRF's own block structure, and a name
/// that repeats is carried at the highest multiplicity any single document used, so nothing is
/// dropped. A cell a document did not have is filled with the reserved word `"not available"`,
/// and a [`SOURCE_DOCUMENT_COLUMN`] records which document each row came from.
///
/// **The result is an analysis table, not something to deposit.** `source name` + `assay name` +
/// `comment[label]` is unique within one document, but two experiments may both have a
/// `"Sample 1"`, so a pooled table will usually violate SDRF's uniqueness rule. Use the
/// source-document column as part of any key.
///
/// No documents, a blank path or label, or a path or label containing a tab or newline, is refused
/// before anything is spawned.
#[doc = include_str!("../docs/reference/sdrf.pool.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{pool_with, PoolInput, PoolOptions, ReadOptions, SOURCE_DOCUMENT_COLUMN};
///
/// let pooled = pool_with(
///     &PoolInput::Labelled(vec![
///         ("PXD000070.sdrf.tsv".into(), "malaria".to_owned()),
///         ("PXD026824.sdrf.tsv".into(), "colon".to_owned()),
///     ]),
///     &PoolOptions { read: ReadOptions { limit: Some(4), ..Default::default() }, out: None },
/// )?;
/// assert!(pooled.document.columns.iter().any(|c| c == SOURCE_DOCUMENT_COLUMN));
/// assert_eq!(pooled.source_documents()[0], Some("malaria"));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.pool.see-also.md")]
pub fn pool_with(documents: &PoolInput, options: &PoolOptions) -> Result<PooledSdrf> {
    let (args, stdin) = pool_request(documents, options)?;
    let data = bridge::invoke(&args, Some(&stdin), options.read.timeout)?;
    serde_json::from_value(data).map_err(protocol)
}

// ---------------------------------------------------------------------------------------------
// Validate, lint, assess, samples and ages (mzLib 1.0.592)
// ---------------------------------------------------------------------------------------------

/// The three verdicts [`assess`] can return, best first, as mzLib spells them
/// (`SdrfSampleVerdict`).
pub const VERDICTS: [&str; 3] = ["Informative", "Partial", "Skeleton"];

/// How much of an age a cell pins down, as mzLib spells it (`SdrfAgePrecision`). See
/// [`parse_ages`].
pub const AGE_PRECISIONS: [&str; 4] = ["Exact", "Range", "LowerBound", "UpperBound"];

/// Why an age cell was refused. mzLib's `SdrfAge.TryParse` answers only yes or no; the bridge
/// names the reason from the cell: `empty`; `reserved_word` (`"not available"` and friends, any
/// case); `no_unit` (a bare number such as `63`, where years and days cannot be told apart);
/// `unreadable`.
pub const AGE_REFUSALS: [&str; 4] = ["empty", "reserved_word", "no_unit", "unreadable"];

/// How a many-document call reads its corpus: how many documents at once, and what one
/// unreadable document does to the rest (BULK.md §1).
///
/// The whole list goes to **one** bridge process, which reads `threads` documents at once and
/// returns one table in input order. There is deliberately no Rust-side loop or thread pool:
/// parallelism is the bridge's to express, once, where it can be counted.
#[derive(Debug, Clone)]
pub struct BulkOptions {
    /// Documents read at once: 1 or more, or `-1` for one per core. The default is 1. The answer is
    /// byte-identical at any value (tested at 1, 4 and -1), so this is a resource choice only.
    pub threads: i32,
    /// [`OnError::Fail`] (the default) returns the first unreadable document's error, in input
    /// order. [`OnError::Skip`] records the failure in that document's `files` entry and reads
    /// the rest.
    pub on_error: OnError,
    /// Time to allow for the whole batch. `None`, the default, waits indefinitely.
    pub timeout: Option<Duration>,
}

impl Default for BulkOptions {
    fn default() -> Self {
        Self {
            threads: 1,
            on_error: OnError::Fail,
            timeout: None,
        }
    }
}

/// Why one document of a `*_many` call produced no result, under [`OnError::Skip`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct FileError {
    /// `"usage"` (the file is missing) or `"correctness"` (mzLib could not read it; the message
    /// then starts with the .NET exception type, e.g. `"MzLibException: ..."`).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub kind: String,
    /// What went wrong, naming the file.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub message: String,
}

// ---- validate -------------------------------------------------------------------------------

/// One finding from mzLib's `SdrfValidator`, located as precisely as its rule allows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidationMessage {
    /// `"Error"` (the document cannot be reliably consumed: a ragged row, a missing required
    /// column, two indistinguishable rows) or `"Warning"` (it deviates from the specification but
    /// can still be joined: casing, an empty cell, a missing recommended column).
    pub severity: String,
    /// The stable rule id, e.g. `"RequiredColumn"`, `"RowWidth"`, `"RowKeyUniqueness"`,
    /// `"ReservedWordCase"`. Stable, so you can count or suppress by it.
    pub rule: String,
    /// Human-readable, including the offending value.
    pub message: String,
    /// 0-based index into the document's data rows; `None` for a finding about the whole
    /// document (a missing column, an empty header, no rows).
    pub row_index: Option<i64>,
    /// The 1-based line in the file, `row_index + 2` because the header is line 1; `None` exactly
    /// when `row_index` is.
    pub line_number: Option<i64>,
    /// The column involved; `None` when the finding is not about one column.
    pub column_name: Option<String>,
}

fn messages(table: &Table) -> Result<Vec<ValidationMessage>> {
    let severity = table.strings("severity")?;
    let rule = table.strings("rule")?;
    let message = table.strings("message")?;
    let row_index = table.integers("row_index")?;
    let line_number = table.integers("line_number")?;
    let column_name = table.strings("column_name")?;
    Ok((0..rule.len())
        .map(|i| ValidationMessage {
            severity: severity[i].clone().unwrap_or_default(),
            rule: rule[i].clone().unwrap_or_default(),
            message: message[i].clone().unwrap_or_default(),
            row_index: row_index[i],
            line_number: line_number[i],
            column_name: column_name[i].clone(),
        })
        .collect())
}

/// The structural findings for one SDRF document — what [`validate`] returns.
///
/// The findings are a table, one row per finding: `severity`, `rule`, `message`, `row_index`,
/// `line_number` and `column_name`, in mzLib's order. [`SdrfValidation::messages`] gives the same
/// rows as [`ValidationMessage`] values.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfValidation {
    /// The absolute path validated.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `true` when there is no `Error` finding. **Warnings never make a document invalid**: mzLib
    /// calibrated every severity against the 1,236-file curated corpus, and a rule that fired on
    /// most curated files was judged wrong rather than the files.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub is_valid: bool,
    /// Findings of severity `Error`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub error_count: u64,
    /// Findings of severity `Warning`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub warning_count: u64,
    /// All findings: the rows of [`Self::columns`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub message_count: u64,
    /// Data rows in the document (the header is not a row).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The findings table, one row per finding, in mzLib's order.
    #[serde(flatten)]
    pub columns: Table,
    /// What a clean result does *not* mean. Read these once.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfValidation {
    /// Every finding as a [`ValidationMessage`], in mzLib's order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn messages(&self) -> Result<Vec<ValidationMessage>> {
        messages(&self.columns)
    }

    /// The findings that make the document invalid.
    ///
    /// # Errors
    ///
    /// As [`Self::messages`].
    pub fn errors(&self) -> Result<Vec<ValidationMessage>> {
        Ok(self
            .messages()?
            .into_iter()
            .filter(|m| m.severity == "Error")
            .collect())
    }

    /// The findings worth fixing that do not make the document invalid.
    ///
    /// # Errors
    ///
    /// As [`Self::messages`].
    pub fn warnings(&self) -> Result<Vec<ValidationMessage>> {
        Ok(self
            .messages()?
            .into_iter()
            .filter(|m| m.severity == "Warning")
            .collect())
    }
}

/// One document's summary in a [`validate_many`] result. Every count is `None` when the document
/// could not be read, and [`Self::error`] then says why.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct ValidatedFile {
    /// The path as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// No `Error` findings; `None` if unread.
    #[serde(default)]
    pub is_valid: Option<bool>,
    /// `Error` findings; `None` if unread.
    #[serde(default)]
    pub error_count: Option<u64>,
    /// `Warning` findings; `None` if unread.
    #[serde(default)]
    pub warning_count: Option<u64>,
    /// All findings; `None` if unread.
    #[serde(default)]
    pub message_count: Option<u64>,
    /// Data rows in the document; `None` if unread.
    #[serde(default)]
    pub row_count: Option<u64>,
    /// Why it was not read, or `None` when it was.
    #[serde(default)]
    pub error: Option<FileError>,
}

/// The findings for many SDRF documents as one long table — what [`validate_many`] returns.
///
/// The table has the columns of [`SdrfValidation`], preceded by `source_index` (the document's
/// 0-based position in the list you passed) and `source_path`. Rows are in input order, then
/// mzLib's order within a document, whatever [`BulkOptions::threads`] was.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfValidationBatch {
    /// Documents given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Documents read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Documents that could not be read; non-zero only under [`OnError::Skip`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// Rows in [`Self::columns`] over every document read — the count every bulk verb reports
    /// (BULK.md §2). Equal to [`Self::message_count`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Documents read with no `Error`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub valid_count: u64,
    /// Findings over every document read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub message_count: u64,
    /// One entry per input, in input order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<ValidatedFile>,
    /// The findings table: `source_index`, `source_path`, then the single-document columns.
    #[serde(flatten)]
    pub columns: Table,
    /// As for one document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfValidationBatch {
    /// Every finding as a [`ValidationMessage`]; read [`Self::columns`]' `source_index` to know
    /// which document each came from.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn messages(&self) -> Result<Vec<ValidationMessage>> {
        messages(&self.columns)
    }

    /// The inputs that could not be read, under [`OnError::Skip`].
    #[must_use]
    pub fn failed_files(&self) -> Vec<&ValidatedFile> {
        self.files.iter().filter(|f| f.error.is_some()).collect()
    }
}

// ---- lint -----------------------------------------------------------------------------------

/// One spelling of one inconsistently written concept: a row of an [`SdrfDrift`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DriftVariant {
    /// Which finding the row belongs to, 0-based, most impactful first.
    pub finding_index: i64,
    /// `AccessionNameConflict`, `NameAccessionConflict`, `MixedTermAndFreeText`,
    /// `ColumnNameVariant` or `ValueCaseVariant` (mzLib's `SdrfDriftKind`).
    pub kind: String,
    /// What was written inconsistently: an accession, a name, a normalised value, or a normalised
    /// column name, depending on `kind`.
    pub concept: String,
    /// The column whose values drift; `None` for `ColumnNameVariant`, a finding about a name.
    pub column_name: Option<String>,
    /// 0 for the majority spelling, then by frequency. Descriptive, not advice.
    pub variant_rank: i64,
    /// This spelling exactly as written.
    pub value: String,
    /// How widely this spelling was used, in documents.
    pub occurrences: i64,
    /// The labels of the documents that used it.
    pub documents: Vec<String>,
}

/// Concepts a set of SDRF documents wrote inconsistently, from mzLib's `SdrfDriftLint` — what
/// [`lint_labelled`] and [`lint`] return.
///
/// The table is long, **one row per finding × spelling**: a finding with three spellings is three
/// rows sharing a `finding_index`. Its columns are `finding_index`, `kind`, `concept`,
/// `column_name`, `variant_rank`, `value`, `occurrences` and `documents` — the last a **list** per
/// row, since labels are your text and no delimiter could be undone.
/// [`SdrfDrift::findings`] groups the rows.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfDrift {
    /// Documents linted.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub document_count: u64,
    /// The paths, in the order given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub paths: Vec<String>,
    /// The label each document is named by in the `documents` column: yours, or mzLib's
    /// `containing-folder/file-stem` default.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub labels: Vec<String>,
    /// Distinct findings (not rows).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub finding_count: u64,
    /// The findings table, one row per finding × spelling, findings in mzLib's impact order and
    /// spellings majority first.
    #[serde(flatten)]
    pub columns: Table,
    /// Including why the majority spelling is not advice.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfDrift {
    /// Every row as a [`DriftVariant`], in the table's order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn variants(&self) -> Result<Vec<DriftVariant>> {
        let t = &self.columns;
        let finding_index = t.integers("finding_index")?;
        let kind = t.strings("kind")?;
        let concept = t.strings("concept")?;
        let column_name = t.strings("column_name")?;
        let variant_rank = t.integers("variant_rank")?;
        let value = t.strings("value")?;
        let occurrences = t.integers("occurrences")?;
        let documents = t.raw("documents").unwrap_or(&[]);
        (0..finding_index.len())
            .map(|i| {
                let labels = match documents.get(i) {
                    Some(serde_json::Value::Array(items)) => items
                        .iter()
                        .map(|item| item.as_str().map(str::to_owned))
                        .collect::<Option<Vec<_>>>(),
                    Some(serde_json::Value::Null) | None => Some(Vec::new()),
                    Some(_) => None,
                }
                .ok_or_else(|| {
                    MzLibError::Protocol(format!(
                        "Column 'documents' row {i} is not a list of labels."
                    ))
                })?;
                Ok(DriftVariant {
                    finding_index: finding_index[i].unwrap_or_default(),
                    kind: kind[i].clone().unwrap_or_default(),
                    concept: concept[i].clone().unwrap_or_default(),
                    column_name: column_name[i].clone(),
                    variant_rank: variant_rank[i].unwrap_or_default(),
                    value: value[i].clone().unwrap_or_default(),
                    occurrences: occurrences[i].unwrap_or_default(),
                    documents: labels,
                })
            })
            .collect()
    }

    /// The rows grouped by finding: one `Vec` of spellings per finding, majority first.
    ///
    /// # Errors
    ///
    /// As [`Self::variants`].
    pub fn findings(&self) -> Result<Vec<Vec<DriftVariant>>> {
        let mut grouped: Vec<Vec<DriftVariant>> = Vec::new();
        for variant in self.variants()? {
            match grouped.last_mut() {
                Some(group) if group[0].finding_index == variant.finding_index => {
                    group.push(variant);
                }
                _ => grouped.push(vec![variant]),
            }
        }
        Ok(grouped)
    }
}

// ---- assess ---------------------------------------------------------------------------------

/// Whether one SDRF document's sample half describes an experimental design — what [`assess`]
/// returns.
///
/// mzLib's `SdrfSampleInformativeness` asks three questions, and the verdict is how many pass: all
/// three is `Informative`, none is `Skeleton`, anything between is `Partial`. The table is the
/// evidence, one row per column a check read: `role` (`factor_value`, `sample_characteristic` or
/// `biological_replicate`), `column_name`, `rows`, `filled` (rows with a real answer), `absent`
/// (empty or a reserved word), `distinct_values` (compared ignoring case and surrounding space)
/// and `fill_rate` (`filled / rows`, a fraction from 0 to 1, which may cross as a whole number:
/// read it with [`Table::floats`]).
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfAssessment {
    /// The absolute path assessed.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `"Informative"`, `"Partial"` or `"Skeleton"` — see [`VERDICTS`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub verdict: String,
    /// Some `factor value[...]` column holds two or more distinct real answers.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub factor_value_varies: bool,
    /// Some `characteristics[...]` column other than organism and biological replicate holds a
    /// real answer. Organism is left out because a search can fill it without a human.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_is_described: bool,
    /// `characteristics[biological replicate]` holds two or more distinct real answers; `false`
    /// also when the column is missing (the table then has no `biological_replicate` row).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub biological_replicate_varies: bool,
    /// Data rows in the document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// The evidence table: factor values, then sample characteristics, then biological
    /// replicate, each in mzLib's order.
    #[serde(flatten)]
    pub columns: Table,
    /// Including when `Partial` is perfectly legitimate.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

/// How many documents of an [`assess_many`] batch got each verdict.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Deserialize)]
pub struct VerdictCounts {
    /// Documents judged `Informative`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub informative: u64,
    /// Documents judged `Partial`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub partial: u64,
    /// Documents judged `Skeleton`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub skeleton: u64,
}

/// One document's verdict in an [`assess_many`] result; `None` fields mean it was not read.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct AssessedFile {
    /// The path as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `"Informative"`, `"Partial"` or `"Skeleton"`; `None` if unread.
    #[serde(default)]
    pub verdict: Option<String>,
    /// As [`SdrfAssessment::factor_value_varies`]; `None` if unread.
    #[serde(default)]
    pub factor_value_varies: Option<bool>,
    /// As [`SdrfAssessment::sample_is_described`]; `None` if unread.
    #[serde(default)]
    pub sample_is_described: Option<bool>,
    /// As [`SdrfAssessment::biological_replicate_varies`]; `None` if unread.
    #[serde(default)]
    pub biological_replicate_varies: Option<bool>,
    /// Data rows; `None` if unread.
    #[serde(default)]
    pub row_count: Option<u64>,
    /// Why it was not read, or `None`.
    #[serde(default)]
    pub error: Option<FileError>,
}

/// Verdicts for many SDRF documents, with the evidence as one long table — what [`assess_many`]
/// returns: the gate for a corpus.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfAssessmentBatch {
    /// Documents given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Documents assessed.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Documents that could not be read; non-zero only under [`OnError::Skip`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// Rows in [`Self::columns`] over every document read — the count every bulk verb reports
    /// (BULK.md §2).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Documents with each verdict, over the documents read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub verdict_counts: VerdictCounts,
    /// One entry per input, in input order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<AssessedFile>,
    /// The evidence table: `source_index`, `source_path`, then the columns of [`SdrfAssessment`].
    #[serde(flatten)]
    pub columns: Table,
    /// As for one document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfAssessmentBatch {
    /// The paths whose verdict is one of `verdicts`, in input order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] for a verdict that is not one of [`VERDICTS`]: a typo would otherwise
    /// silently select nothing.
    pub fn paths_with(&self, verdicts: &[&str]) -> Result<Vec<&str>> {
        if let Some(unknown) = verdicts.iter().find(|v| !VERDICTS.contains(v)) {
            return Err(MzLibError::Usage(format!(
                "Unknown verdict '{unknown}'; expected some of {VERDICTS:?}."
            )));
        }
        Ok(self
            .files
            .iter()
            .filter(|f| f.verdict.as_deref().is_some_and(|v| verdicts.contains(&v)))
            .map(|f| f.path.as_str())
            .collect())
    }
}

// ---- samples --------------------------------------------------------------------------------

/// One sample's `characteristics[age]`, read into years: a row of [`SdrfSamples::ages`].
#[derive(Debug, Clone, PartialEq)]
pub struct SampleAge {
    /// The sample, as its first row spelled it.
    pub source_name: String,
    /// The cell, verbatim.
    pub value: Option<String>,
    /// The age, a range's midpoint, or a bound, in years; `None` when refused.
    pub years: Option<f64>,
    /// The youngest the cell allows, in years; `None` when refused.
    pub min_years: Option<f64>,
    /// The oldest the cell allows, in years; `None` when refused **or unbounded** (a
    /// `LowerBound` such as `>=90Y`).
    pub max_years: Option<f64>,
    /// One of [`AGE_PRECISIONS`]; `None` when refused.
    pub precision: Option<String>,
    /// `true` for the specification's grammar, `false` for unambiguous words; `None` when refused.
    pub follows_specification: Option<bool>,
    /// `None` when the cell was read, else one of [`AGE_REFUSALS`].
    pub refusal: Option<String>,
}

/// Every sample in one SDRF document and what its rows agree on, as a long table — what
/// [`samples`] returns.
///
/// A sample is a `source name`, matched ignoring case and surrounding space, exactly as mzLib's
/// `SdrfSampleBlock.BySourceName` keys it: one sample measured in three fractions is one sample.
/// The table has one row per sample × column × position, samples in document order, with
/// `source_name`, `sample_row_count`, `column_name` (the header's own spelling), `column_kind`
/// (`source_name`, `characteristic` or `factor_value`), `position` (`None` on a conflicting row),
/// `status` (`agreed`, or `conflicting` when the sample's rows disagree and mzLib withheld the
/// column), `value` (verbatim; `None` when withheld), and the six `age_*` columns, filled only on
/// `characteristics[age]` rows with the meanings [`ParsedAges`] documents.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfSamples {
    /// The absolute path read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Distinct samples (source names).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// Data rows in the document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub row_count: u64,
    /// Sample × column pairs mzLib withheld because the sample's rows disagree.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub conflict_count: u64,
    /// Rows mzLib could not place: no `source name` column, or a blank source name. Empty for a
    /// well-formed document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub problems: Vec<String>,
    /// The samples table.
    #[serde(flatten)]
    pub columns: Table,
    /// Including how to key samples across documents.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

fn sample_ages(table: &Table) -> Result<Vec<SampleAge>> {
    let source_name = table.strings("source_name")?;
    let value = table.strings("value")?;
    let years = table.floats("age_years")?;
    let min_years = table.floats("age_min_years")?;
    let max_years = table.floats("age_max_years")?;
    let precision = table.strings("age_precision")?;
    let follows = table.booleans("age_follows_specification")?;
    let refusal = table.strings("age_refusal")?;
    Ok((0..source_name.len())
        .filter(|&i| years[i].is_some() || refusal[i].is_some())
        .map(|i| SampleAge {
            source_name: source_name[i].clone().unwrap_or_default(),
            value: value[i].clone(),
            years: years[i],
            min_years: min_years[i],
            max_years: max_years[i],
            precision: precision[i].clone(),
            follows_specification: follows[i],
            refusal: refusal[i].clone(),
        })
        .collect())
}

impl SdrfSamples {
    /// Only the parsed-age rows: one per sample whose `characteristics[age]` was agreed.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn ages(&self) -> Result<Vec<SampleAge>> {
        sample_ages(&self.columns)
    }

    /// `(source_name, column_name)` for every column mzLib withheld as conflicting.
    ///
    /// # Errors
    ///
    /// As [`Self::ages`].
    pub fn conflicts(&self) -> Result<Vec<(String, String)>> {
        let status = self.columns.strings("status")?;
        let source_name = self.columns.strings("source_name")?;
        let column_name = self.columns.strings("column_name")?;
        Ok((0..status.len())
            .filter(|&i| status[i].as_deref() == Some("conflicting"))
            .map(|i| {
                (
                    source_name[i].clone().unwrap_or_default(),
                    column_name[i].clone().unwrap_or_default(),
                )
            })
            .collect())
    }
}

/// One document's summary in a [`samples_many`] result; `None` fields mean it was not read.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct SampledFile {
    /// The path as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Distinct samples; `None` if unread.
    #[serde(default)]
    pub sample_count: Option<u64>,
    /// Data rows; `None` if unread.
    #[serde(default)]
    pub row_count: Option<u64>,
    /// Columns withheld as conflicting; `None` if unread.
    #[serde(default)]
    pub conflict_count: Option<u64>,
    /// Rows mzLib could not place; `None` if unread.
    #[serde(default)]
    pub problems: Option<Vec<String>>,
    /// Why it was not read, or `None`.
    #[serde(default)]
    pub error: Option<FileError>,
}

/// The samples of many SDRF documents, as one long table — what [`samples_many`] returns.
///
/// The table has the columns of [`SdrfSamples`], preceded by `source_index` and `source_path`.
/// **A source name is only unique within one document** — two studies may both have a
/// `"Sample 1"` — so key a sample on `(source_index, source_name)`.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfSamplesBatch {
    /// Documents given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Documents read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub read_count: u64,
    /// Documents that could not be read; non-zero only under [`OnError::Skip`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub failed_count: u64,
    /// Rows in [`Self::columns`] over every document read — the count every bulk verb reports
    /// (BULK.md §2).
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Samples over every document read.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub sample_count: u64,
    /// One entry per input, in input order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub files: Vec<SampledFile>,
    /// The samples table: `source_index`, `source_path`, then the columns of [`SdrfSamples`].
    #[serde(flatten)]
    pub columns: Table,
    /// As for one document.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

// ---- parse_ages -----------------------------------------------------------------------------

/// One age cell read into years: a row of [`ParsedAges::ages`].
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedAge {
    /// The cell exactly as given.
    pub cell: String,
    /// The single figure to place the sample on an age axis, in years: the age, a range's
    /// **midpoint**, or a bound. `None` when refused.
    pub years: Option<f64>,
    /// The youngest the cell allows, in years; `0` for an upper bound. `None` when refused.
    pub min_years: Option<f64>,
    /// The oldest the cell allows, in years. `None` when refused **or** when there is no upper
    /// bound (`precision == "LowerBound"`, e.g. `>=90Y`): mzLib's +infinity, which JSON cannot
    /// carry. [`Self::precision`] tells the two apart.
    pub max_years: Option<f64>,
    /// One of [`AGE_PRECISIONS`]; `None` when refused.
    pub precision: Option<String>,
    /// `true` for the specification's `nYnMnD` / `nW` grammar, `false` for unambiguous words
    /// (`"3 year"`, `"6-8 weeks"`); `None` when refused.
    pub follows_specification: Option<bool>,
    /// `None` when the cell was read, else one of [`AGE_REFUSALS`].
    pub refusal: Option<String>,
}

/// Age cells read into years by mzLib's `SdrfAge.TryParse`, one row per cell given, in input
/// order — what [`parse_ages`] returns.
///
/// The table's columns are `cell`, `years`, `min_years`, `max_years`, `precision`,
/// `follows_specification` and `refusal`, with the meanings of [`ParsedAge`]. Every age is in
/// **years**: a month is 1/12 year, a week 7/365.25, a day 1/365.25, an hour 1/(365.25 × 24).
#[derive(Debug, Clone, Deserialize)]
pub struct ParsedAges {
    /// Cells given: the rows of [`Self::columns`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub cell_count: u64,
    /// Cells read as an age.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub parsed_count: u64,
    /// The ages table, one row per input cell.
    #[serde(flatten)]
    pub columns: Table,
    /// The unit conventions and refusal rules.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl ParsedAges {
    /// Every row as a [`ParsedAge`], row `i` for cell `i`.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is not the type the wire contract says it is.
    pub fn ages(&self) -> Result<Vec<ParsedAge>> {
        let t = &self.columns;
        let cell = t.strings("cell")?;
        let years = t.floats("years")?;
        let min_years = t.floats("min_years")?;
        let max_years = t.floats("max_years")?;
        let precision = t.strings("precision")?;
        let follows = t.booleans("follows_specification")?;
        let refusal = t.strings("refusal")?;
        Ok((0..cell.len())
            .map(|i| ParsedAge {
                cell: cell[i].clone().unwrap_or_default(),
                years: years[i],
                min_years: min_years[i],
                max_years: max_years[i],
                precision: precision[i].clone(),
                follows_specification: follows[i],
                refusal: refusal[i].clone(),
            })
            .collect())
    }
}

// ---- argument assembly for the analysis verbs -----------------------------------------------

fn one_path_args(verb: &str, path: &Path) -> Result<Vec<String>> {
    let path = path_text(path, "file path")?;
    if path.is_empty() {
        return Err(MzLibError::Usage(
            "A file path is required, e.g. 'PXD000070.sdrf.tsv'.".to_owned(),
        ));
    }
    Ok(vec![
        "sdrf".to_owned(),
        verb.to_owned(),
        "--path".to_owned(),
        path.to_owned(),
    ])
}

/// The arguments and stdin of a `*_many` call: the list on stdin, one path per line, and the
/// bulk options as flags. One process reads the whole corpus; nothing here loops over files.
fn bulk_request<P: AsRef<Path>>(
    verb: &str,
    paths: &[P],
    options: &BulkOptions,
) -> Result<(Vec<String>, String)> {
    let lines = bridge::path_lines(paths, "SDRF document")?;
    let args = vec![
        "sdrf".to_owned(),
        verb.to_owned(),
        "--paths-stdin".to_owned(),
        "--threads".to_owned(),
        bridge::threads_arg(options.threads)?,
        "--on-error".to_owned(),
        options.on_error.as_str().to_owned(),
    ];
    Ok((args, lines.join("\n")))
}

/// The stdin of `parse_ages`: one cell per line, blank cells kept, and a trailing newline so a
/// final blank cell survives the trip.
fn age_lines<S: AsRef<str>>(cells: &[S]) -> Result<String> {
    if cells.is_empty() {
        return Err(MzLibError::Usage(
            "parse_ages needs at least one cell.".to_owned(),
        ));
    }
    let mut text = String::new();
    for (index, cell) in cells.iter().enumerate() {
        let cell = cell.as_ref();
        if cell.contains(['\n', '\r']) {
            return Err(MzLibError::Usage(format!(
                "Cell {index} contains a line break, which separates cells: {cell:?}."
            )));
        }
        text.push_str(cell);
        text.push('\n');
    }
    Ok(text)
}

fn call<T: serde::de::DeserializeOwned>(
    args: &[String],
    stdin: Option<&str>,
    timeout: Option<Duration>,
) -> Result<T> {
    let data = bridge::invoke(args, stdin, timeout)?;
    serde_json::from_value(data).map_err(protocol)
}

// ---- the verbs ------------------------------------------------------------------------------

/// Check one SDRF file against the SDRF-Proteomics specification's structural rules, one row per
/// finding.
///
/// Calls mzLib's `SdrfValidator.Validate`. **Structure only** — required and recommended columns,
/// column order, casing, malformed names, ragged rows, integer replicate and fraction columns,
/// reserved-word casing, and the one hard row rule: `source name` + `assay name` +
/// `comment[label]` must be unique. Controlled-vocabulary accessions are **not** resolved. A
/// document of reserved words validates cleanly: use [`assess`] for whether it says anything, and
/// [`lint_labelled`] for whether several agree.
#[doc = include_str!("../docs/reference/sdrf.validate.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let result = mzlib::sdrf::validate("sdrf_cohort.sdrf.tsv")?;
/// assert_eq!((result.is_valid, result.error_count, result.warning_count), (true, 0, 2));
/// for m in result.warnings()? {
///     println!("line {:?}: {} in {:?}", m.line_number, m.rule, m.column_name);
/// }
/// assert_eq!(result.warnings()?[0].line_number, Some(10));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// A file generated from a list of data files fails on structure too:
///
/// ```
/// # mzlib_replay::activate();
/// let skeleton = mzlib::sdrf::validate("sdrf_skeleton.sdrf.tsv")?;
/// let first = &skeleton.errors()?[0];
/// assert!(!skeleton.is_valid);
/// assert_eq!((first.rule.as_str(), first.line_number), ("RequiredColumn", None)); // whole-document
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.validate.see-also.md")]
pub fn validate(path: impl AsRef<Path>) -> Result<SdrfValidation> {
    let args = one_path_args("validate", path.as_ref())?;
    call(&args, None, Some(DEFAULT_TIMEOUT))
}

/// Validate many SDRF files in one bridge call.
///
/// One process for the whole corpus instead of one per file — roughly 120 ms saved per document —
/// with the parallelism inside the bridge. The result does not depend on
/// [`BulkOptions::threads`]: rows are always in input order.
#[doc = include_str!("../docs/reference/sdrf.validate.bulk.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{validate_many, BulkOptions, OnError};
///
/// let batch = validate_many(
///     &["sdrf_skeleton.sdrf.tsv", "sdrf_cohort.sdrf.tsv", "missing.sdrf.tsv"],
///     &BulkOptions { on_error: OnError::Skip, ..Default::default() },
/// )?;
/// assert_eq!((batch.file_count, batch.read_count, batch.valid_count), (3, 2, 1));
/// let valid: Vec<Option<bool>> = batch.files.iter().map(|f| f.is_valid).collect();
/// assert_eq!(valid, [Some(false), Some(true), None]);        // None: not read
/// assert_eq!(batch.failed_files()[0].error.as_ref().unwrap().kind, "usage");
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors
///
/// As [`validate`], per the bulk rules above: an empty list, a blank path or a bad `threads` is
/// refused before anything is spawned.
pub fn validate_many<P: AsRef<Path>>(
    paths: &[P],
    options: &BulkOptions,
) -> Result<SdrfValidationBatch> {
    let (args, stdin) = bulk_request("validate", paths, options)?;
    call(&args, Some(&stdin), options.timeout)
}

/// Find the concepts a set of SDRF documents wrote inconsistently, naming each document by
/// mzLib's `containing-folder/file-stem` default.
///
/// Prefer [`lint_labelled`]: the default label depends on where the files sit. See it for the
/// reference.
///
/// # Errors
///
/// As [`lint_labelled`].
pub fn lint<P: AsRef<Path>>(paths: &[P]) -> Result<SdrfDrift> {
    lint_documents(&PoolInput::Paths(
        paths.iter().map(|p| p.as_ref().to_path_buf()).collect(),
    ))
}

/// Find the concepts a set of SDRF documents wrote inconsistently, one row per finding and
/// spelling, naming each document yourself.
///
/// Validity is a property of one file; comparability is a property of the relationship between
/// files. Every document can pass [`validate`] and the pooled table still be unusable because one
/// file writes `"Homo sapiens"` and another `"homo sapiens"`. This is mzLib's `SdrfDriftLint`
/// check for that. **Label the documents**: the labels are how each row's `documents` names them.
#[doc = include_str!("../docs/reference/sdrf.lint.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let drift = mzlib::sdrf::lint_labelled(&[
///     ("sdrf_cohort.sdrf.tsv", "cohort"),
///     ("sdrf_cohort_partner.sdrf.tsv", "partner"),
/// ])?;
/// assert_eq!(drift.finding_count, 4);
/// for spellings in drift.findings()? {
///     let values: Vec<&str> = spellings.iter().map(|v| v.value.as_str()).collect();
///     println!("{}: {values:?}", spellings[0].kind);
/// }
/// let last = drift.findings()?.pop().unwrap();
/// assert_eq!(last[0].kind, "ValueCaseVariant");
/// assert_eq!((last[0].value.as_str(), last[1].value.as_str()), ("Homo sapiens", "homo sapiens"));
/// assert_eq!(last[1].documents, ["partner"]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.lint.see-also.md")]
pub fn lint_labelled<P: AsRef<Path>, L: AsRef<str>>(documents: &[(P, L)]) -> Result<SdrfDrift> {
    lint_documents(&PoolInput::Labelled(
        documents
            .iter()
            .map(|(p, l)| (p.as_ref().to_path_buf(), l.as_ref().to_owned()))
            .collect(),
    ))
}

fn lint_documents(documents: &PoolInput) -> Result<SdrfDrift> {
    let lines = document_lines(documents)?;
    let args = ["sdrf".to_owned(), "lint".to_owned()];
    call(&args, Some(&lines.join("\n")), Some(DEFAULT_TIMEOUT))
}

/// Decide whether one SDRF file's sample half describes an experimental design — `Informative`,
/// `Partial` or `Skeleton` — with the per-column counts behind the verdict.
///
/// A file generated from a list of data files — reserved words in every sample column, one
/// replicate number, no factor — passes [`validate`], because reserved words are the
/// specification's correct way to say nothing. For grouping results by biology it is the same as
/// having no SDRF at all. mzLib's `SdrfSampleInformativeness.Assess` is the gate for that.
#[doc = include_str!("../docs/reference/sdrf.assess.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let a = mzlib::sdrf::assess("sdrf_cohort.sdrf.tsv")?;
/// assert_eq!(a.verdict, "Informative");
/// assert!(a.factor_value_varies && a.sample_is_described && a.biological_replicate_varies);
/// assert_eq!(a.columns.strings("column_name")?[0].as_deref(), Some("factor value[disease]"));
/// assert_eq!(a.columns.integers("distinct_values")?[0], Some(2));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.assess.see-also.md")]
pub fn assess(path: impl AsRef<Path>) -> Result<SdrfAssessment> {
    let args = one_path_args("assess", path.as_ref())?;
    call(&args, None, Some(DEFAULT_TIMEOUT))
}

/// Assess many SDRF files in one bridge call — the gate for a corpus.
#[doc = include_str!("../docs/reference/sdrf.assess.bulk.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{assess_many, BulkOptions};
///
/// let batch = assess_many(
///     &["sdrf_cohort.sdrf.tsv", "sdrf_skeleton.sdrf.tsv", "PXD000070.sdrf.tsv"],
///     &BulkOptions::default(),
/// )?;
/// let counts = batch.verdict_counts;
/// assert_eq!((counts.informative, counts.partial, counts.skeleton), (1, 1, 1));
/// let usable = batch.paths_with(&["Informative", "Partial"])?;
/// assert_eq!(usable.len(), 2);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors
///
/// As [`assess`], with the bulk rules of [`validate_many`].
pub fn assess_many<P: AsRef<Path>>(
    paths: &[P],
    options: &BulkOptions,
) -> Result<SdrfAssessmentBatch> {
    let (args, stdin) = bulk_request("assess", paths, options)?;
    call(&args, Some(&stdin), options.timeout)
}

/// One row per sample, column and position: what each source name's rows agree on, the columns
/// they disagree about named and withheld, and `characteristics[age]` read into years.
///
/// Calls mzLib's `SdrfSampleBlock.BySourceName`: the sample half of the document (`source name`,
/// every `characteristics[...]`, every `factor value[...]`) per sample, merged over the sample's
/// rows. Where those rows **disagree** — one fraction says `normal` and another `COVID-19` —
/// mzLib names the column and withholds it rather than let the first row win.
#[doc = include_str!("../docs/reference/sdrf.samples.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let s = mzlib::sdrf::samples("sdrf_cohort.sdrf.tsv")?;
/// assert_eq!(s.sample_count, 6);
/// assert_eq!(s.conflicts()?, [("S6".to_owned(), "characteristics[disease]".to_owned())]);
/// for age in s.ages()? {
///     println!("{} {:?} {:?} {:?}", age.source_name, age.value, age.years, age.refusal);
/// }
/// let ages = s.ages()?;
/// assert_eq!((ages[1].years, ages[1].precision.as_deref()), (Some(62.5), Some("Range")));
/// assert_eq!(ages[3].refusal.as_deref(), Some("no_unit"));            // "63": years or days?
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.samples.see-also.md")]
pub fn samples(path: impl AsRef<Path>) -> Result<SdrfSamples> {
    let args = one_path_args("samples", path.as_ref())?;
    call(&args, None, Some(DEFAULT_TIMEOUT))
}

/// Lift the samples of many SDRF files in one bridge call.
#[doc = include_str!("../docs/reference/sdrf.samples.bulk.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{samples_many, BulkOptions};
///
/// let batch = samples_many(
///     &["sdrf_cohort.sdrf.tsv", "sdrf_cohort_partner.sdrf.tsv"],
///     &BulkOptions::default(),
/// )?;
/// let per_file: Vec<Option<u64>> = batch.files.iter().map(|f| f.sample_count).collect();
/// assert_eq!((batch.sample_count, per_file), (8, vec![Some(6), Some(2)]));
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// # Errors
///
/// As [`samples`], with the bulk rules of [`validate_many`].
pub fn samples_many<P: AsRef<Path>>(
    paths: &[P],
    options: &BulkOptions,
) -> Result<SdrfSamplesBatch> {
    let (args, stdin) = bulk_request("samples", paths, options)?;
    call(&args, Some(&stdin), options.timeout)
}

/// Read `characteristics[age]` cells into years with an honest precision, refusing any cell that
/// would need a guess.
///
/// Calls mzLib's `SdrfAge.TryParse` on each cell. It reads the specification's grammar (`58Y`,
/// `30Y6M`, `16W`), ranges (`40Y-85Y`, `6-8 weeks`), bounds (`>=90Y`, `<1Y`) and unambiguous
/// words (`3 year`, `4 hour`). It **refuses** a bare number — `63` is 11% of real age cells, and 63
/// years and 63 days are both plausible in one study — as well as reserved words and free text.
/// One bridge call for any number of cells, so pass them all at once rather than looping. An empty
/// string is sent as an empty cell and comes back refused as `"empty"`, so the result stays
/// aligned with the input.
#[doc = include_str!("../docs/reference/sdrf.parse-age.md")]
///
/// # Examples
///
/// ```
/// # mzlib_replay::activate();
/// let ages = mzlib::sdrf::parse_ages(&[
///     "58Y", "30Y6M", "40Y-85Y", "40Y-40Y", ">=90Y", "<1Y", "6-8 weeks",
///     "63", "not available", "", "about forty",
/// ])?;
/// assert_eq!((ages.parsed_count, ages.cell_count), (7, 11));
/// let rows = ages.ages()?;
/// assert_eq!((rows[2].years, rows[2].min_years, rows[2].max_years), (Some(62.5), Some(40.0), Some(85.0)));
/// assert_eq!((rows[4].max_years, rows[4].precision.as_deref()), (None, Some("LowerBound")));
/// let refusals: Vec<Option<&str>> = rows[7..].iter().map(|r| r.refusal.as_deref()).collect();
/// assert_eq!(refusals, [Some("no_unit"), Some("reserved_word"), Some("empty"), Some("unreadable")]);
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.parse-age.see-also.md")]
pub fn parse_ages<S: AsRef<str>>(cells: &[S]) -> Result<ParsedAges> {
    let stdin = age_lines(cells)?;
    let args = ["sdrf".to_owned(), "parse-age".to_owned()];
    call(&args, Some(&stdin), Some(DEFAULT_TIMEOUT))
}

// ---- design ---------------------------------------------------------------------------------

/// How [`design_with`] reads a design. Every default is mzLib's.
#[derive(Debug, Clone)]
pub struct DesignOptions {
    /// The `factor value[...]` columns the condition is built from, by exact (case-sensitive)
    /// name, in the order their values are joined with `_`.
    ///
    /// **Empty** (the default) sends nothing, and mzLib uses the document's only factor column —
    /// and refuses a document with several, because joining every factor could split a condition
    /// on a nuisance factor. A blank name, or one holding a tab or line break, is a usage error.
    pub condition_columns: Vec<String>,
    /// The files the search will read, as paths or bare names. Each must be named **exactly**
    /// (case and extension) by one SDRF row; rows for other files are dropped and reported in
    /// [`SdrfDesign::notes`], and these paths become `full_path`.
    ///
    /// `None` (the default) takes the SDRF's own file names. `Some(vec![])` is sent as an empty
    /// list, which mzLib refuses as a design, not as a usage error.
    pub searched_files: Option<Vec<PathBuf>>,
    /// Also write MetaMorpheus's `ExperimentalDesign.tsv` (1-based) here. Must end in `.tsv`,
    /// checked before the SDRF is read. Written only when the design is valid: a refused design
    /// writes nothing and [`SdrfDesign::written`] is `None`. MetaMorpheus finds the file only when
    /// it is named `ExperimentalDesign.tsv` and sits beside the spectra.
    pub out: Option<PathBuf>,
    /// The longest the call may take. `None` waits indefinitely.
    pub timeout: Option<Duration>,
}

impl Default for DesignOptions {
    fn default() -> Self {
        Self {
            condition_columns: Vec::new(),
            searched_files: None,
            out: None,
            timeout: Some(DEFAULT_TIMEOUT),
        }
    }
}

/// The `ExperimentalDesign.tsv` that [`design_with`] wrote under [`DesignOptions::out`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct DesignWritten {
    /// The path written, as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Runs written, in files.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
}

/// One run of a label-free design, in mzLib's `SpectraFileInfo` coordinates: a row of
/// [`SdrfDesign::files`].
///
/// The design coordinates are **0-based**, exactly as [`crate::flashlfq::quantify_with`] and
/// [`crate::flashlfq::median_polish_with`] take them. `ExperimentalDesign.tsv` writes the same
/// numbers plus one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DesignedFile {
    /// The file as the SDRF names it, or the searched path you passed for it.
    pub full_path: String,
    /// `full_path`'s file name without its extension — the run key FlashLFQ uses
    /// (`Intensity_<file_name>`).
    pub file_name: String,
    /// The condition columns' values joined with `_`.
    pub condition: String,
    /// 0-based, ranked within the condition.
    pub biological_replicate: u32,
    /// 0-based.
    pub technical_replicate: u32,
    /// 0-based.
    pub fraction: u32,
}

/// A label-free experimental design read from an SDRF, or every reason it was refused — what
/// [`design_with`] returns.
///
/// **A refusal is an answer, not an error.** When [`Self::is_valid`] is `false`,
/// [`Self::refusals`] lists every reason at once and the table is empty, and [`Self::spectra`]
/// and [`Self::run_design`] refuse, so a design MetaMorpheus would reject cannot reach FlashLFQ
/// by mistake. mzLib refuses rather than repairs because MetaMorpheus skips quantification with
/// only a warning when its design file is invalid.
///
/// The table has one row per run, in SDRF row order: `full_path`, `file_name`, `condition`,
/// `biological_replicate`, `technical_replicate` and `fraction`, with the meanings of
/// [`DesignedFile`]. [`Self::files`] gives the same rows as structs.
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfDesign {
    /// The SDRF path as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `true` when there are no refusals. `false` is an answer: the table is then empty.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub is_valid: bool,
    /// The column runs were keyed on — `comment[searched data file]` when the SDRF has it, else
    /// `comment[data file]`. `None` when it has neither, so no row names a file; that is itself a
    /// refusal.
    #[serde(default)]
    pub file_key_column: Option<String>,
    /// The factor columns the condition was built from, in join order; empty when none could be
    /// chosen.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub condition_columns: Vec<String>,
    /// Whether [`DesignOptions::condition_columns`] was given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub condition_columns_declared: bool,
    /// Whether [`DesignOptions::searched_files`] was given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub searched_files_given: bool,
    /// Runs in the design, in files; `0` when refused.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Every reason the design was refused, in mzLib's words. Empty when valid.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub refusals: Vec<String>,
    /// Every relabelling on the way to a valid design: biological replicates ranked within a
    /// condition, with the `old -> new` mapping, and rows dropped because the search does not
    /// read their file. **The only record of a renumbering — read it.**
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub notes: Vec<String>,
    /// mzLib's human-readable account of all of the above (`SdrfLabelFreeDesign.Report`). Print
    /// it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub report: String,
    /// The design table, one row per run.
    #[serde(flatten)]
    pub columns: Table,
    /// The `ExperimentalDesign.tsv` written; `None` when [`DesignOptions::out`] was not given,
    /// or the design was refused and nothing was written.
    #[serde(default)]
    pub written: Option<DesignWritten>,
    /// What the design does and does not cover — read these once.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfDesign {
    /// Every run as a [`DesignedFile`], in SDRF row order. Empty when refused.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is missing or not the type the wire contract says.
    pub fn files(&self) -> Result<Vec<DesignedFile>> {
        if self.columns.names().is_empty() {
            return Ok(Vec::new());
        }
        let full_path = self.columns.strings("full_path")?;
        let file_name = self.columns.strings("file_name")?;
        let condition = self.columns.strings("condition")?;
        let biological = self.columns.integers("biological_replicate")?;
        let technical = self.columns.integers("technical_replicate")?;
        let fraction = self.columns.integers("fraction")?;
        let coordinate = |value: Option<i64>, what: &str| -> Result<u32> {
            value.and_then(|v| u32::try_from(v).ok()).ok_or_else(|| {
                MzLibError::Protocol(format!(
                    "sdrf design returned a {what} that is not a 0-based integer: {value:?}"
                ))
            })
        };
        (0..full_path.len())
            .map(|i| {
                Ok(DesignedFile {
                    full_path: full_path[i].clone().unwrap_or_default(),
                    file_name: file_name[i].clone().unwrap_or_default(),
                    condition: condition[i].clone().unwrap_or_default(),
                    biological_replicate: coordinate(biological[i], "biological_replicate")?,
                    technical_replicate: coordinate(technical[i], "technical_replicate")?,
                    fraction: coordinate(fraction[i], "fraction")?,
                })
            })
            .collect()
    }

    fn require_valid(&self, what: &str) -> Result<()> {
        if self.is_valid {
            return Ok(());
        }
        Err(MzLibError::Usage(format!(
            "The design was refused, so there is no {what} to give. mzLib's reasons:\n  {}",
            self.refusals.join("\n  ")
        )))
    }

    /// The design as [`crate::flashlfq::quantify_with`] takes its spectra: one
    /// [`SpectraFile`](crate::flashlfq::SpectraFile) per run, `path` set to `full_path` and the
    /// four design fields set, so the SDRF drives FlashLFQ with no hand-written design.
    ///
    /// Pass [`DesignOptions::searched_files`] when the SDRF names files without the directory you
    /// keep them in.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] when the design was refused (mzLib's `ToExperimentalDesign` refuses
    /// too); otherwise as [`Self::files`].
    pub fn spectra(&self) -> Result<Vec<crate::flashlfq::SpectraFile>> {
        self.require_valid("spectra list")?;
        Ok(self
            .files()?
            .into_iter()
            .map(|f| crate::flashlfq::SpectraFile {
                path: f.full_path,
                condition: Some(f.condition),
                biological_replicate: Some(f.biological_replicate),
                technical_replicate: Some(f.technical_replicate),
                fraction: Some(f.fraction),
            })
            .collect())
    }

    /// The design as [`crate::flashlfq::MedianPolishOptions::design`] takes it: one
    /// [`DesignEntry`](crate::flashlfq::DesignEntry) per run, keyed by `file_name` (the
    /// `Intensity_<file_name>` column of a `QuantifiedPeptides.tsv`).
    ///
    /// # Errors
    ///
    /// As [`Self::spectra`].
    pub fn run_design(&self) -> Result<Vec<crate::flashlfq::DesignEntry>> {
        self.require_valid("run design")?;
        Ok(self
            .files()?
            .into_iter()
            .map(|f| crate::flashlfq::DesignEntry {
                file_name: f.file_name,
                condition: Some(f.condition),
                biological_replicate: Some(f.biological_replicate),
                technical_replicate: Some(f.technical_replicate),
                fraction: Some(f.fraction),
            })
            .collect())
    }
}

/// The arguments and stdin of a `sdrf design` call, validated before anything is spawned.
fn design_request(path: &Path, options: &DesignOptions) -> Result<(Vec<String>, Option<String>)> {
    let mut args = one_path_args("design", path)?;
    condition_columns_args(&mut args, &options.condition_columns)?;
    let stdin = searched_files_stdin(&mut args, options.searched_files.as_ref())?;

    if let Some(out) = &options.out {
        let out = path_text(out, "out path")?;
        if out.is_empty() {
            return Err(MzLibError::Usage(
                "out must be a file path ending in .tsv, e.g. 'ExperimentalDesign.tsv'.".to_owned(),
            ));
        }
        args.push("--out".to_owned());
        args.push(out.to_owned());
    }
    Ok((args, stdin))
}

/// Read a label-free experimental design out of an SDRF with every default: the document's only
/// factor column as the condition, and the SDRF's own file names.
///
/// See [`design_with`] for the reference.
///
/// # Errors
///
/// As [`design_with`].
pub fn design(path: impl AsRef<Path>) -> Result<SdrfDesign> {
    design_with(path, &DesignOptions::default())
}

/// Read a label-free experimental design out of an SDRF, in MetaMorpheus's terms, or every
/// reason it cannot be read.
///
/// Calls mzLib's `SdrfLabelFreeDesign.Read` (mzLib #1363). It runs every check MetaMorpheus's own
/// design validator runs and **refuses rather than repairs**: a factor written `not available`,
/// two files claiming the same replicate and fraction, a document with several factor columns and
/// none declared — each is a refusal, and all are reported at once. It relabels in only two ways,
/// both recorded in [`SdrfDesign::notes`]: study-wide biological replicate numbers are ranked
/// within each condition (`22 -> 1`), and rows for files the search does not read are dropped.
///
/// **Label-free only.** An isobaric (TMT, iTRAQ, DiLeu) SDRF needs a channel design: read it
/// with [`isobaric_design_with`]. A search that read `X-calib.mzML` for an SDRF naming `X.raw`
/// matches no row by exact name; [`restrict_with`] the SDRF first.
#[doc = include_str!("../docs/reference/sdrf.design.md")]
///
/// # Examples
///
/// A valid two-factor design — PXD067622, two genotypes by four treatments, three biological
/// replicates each:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{design_with, DesignOptions};
///
/// let both = DesignOptions {
///     condition_columns: vec!["factor value[genotype]".into(), "factor value[treatment]".into()],
///     ..Default::default()
/// };
/// let d = design_with("PXD067622.sdrf.tsv", &both)?;
/// assert!(d.is_valid);
/// assert_eq!(d.file_count, 24);
/// assert_eq!(d.file_key_column.as_deref(), Some("comment[data file]"));
///
/// let first = &d.files()?[0];
/// assert_eq!(first.file_name, "20240830_HF_LC3_MAA_RK_12032_CA_DMSO4");
/// assert_eq!(first.condition, "SPRTN-TurboID CA_DMSO (vehicle)");
/// assert_eq!(first.biological_replicate, 0);                 // 0-based, as FlashLFQ takes it
///
/// // Straight into FlashLFQ, with no hand-written design.
/// let spectra = d.spectra()?;
/// assert_eq!(spectra[0].path, "20240830_HF_LC3_MAA_RK_12032_CA_DMSO4.raw");
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// A refusal lists every reason at once — here one per run, because PXD049018's treatment is
/// `not available`:
///
/// ```
/// # mzlib_replay::activate();
/// # use mzlib::sdrf::{design_with, DesignOptions};
/// # let both = DesignOptions {
/// #     condition_columns: vec!["factor value[genotype]".into(), "factor value[treatment]".into()],
/// #     ..Default::default()
/// # };
/// let refused = design_with("PXD049018.sdrf.tsv", &both)?;
/// assert_eq!((refused.is_valid, refused.refusals.len(), refused.file_count), (false, 20, 0));
/// assert!(refused.refusals[0].starts_with(
///     "Line 2 (MSB67868ABand_01.raw): 'factor value[treatment]' is 'not available'."
/// ));
/// assert!(refused.spectra().is_err());                       // cannot reach FlashLFQ by mistake
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.design.see-also.md")]
pub fn design_with(path: impl AsRef<Path>, options: &DesignOptions) -> Result<SdrfDesign> {
    let (args, stdin) = design_request(path.as_ref(), options)?;
    call(&args, stdin.as_deref(), options.timeout)
}

// ---- shared by design and isobaric design ---------------------------------------------------

/// `--condition-columns` for [`design_with`] and [`isobaric_design_with`]: nothing for an empty
/// list, so mzLib uses the document's only factor column.
fn condition_columns_args(args: &mut Vec<String>, condition_columns: &[String]) -> Result<()> {
    if condition_columns.is_empty() {
        return Ok(());
    }
    for (index, name) in condition_columns.iter().enumerate() {
        if name.trim().is_empty() {
            return Err(MzLibError::Usage(format!(
                "Condition column {index} is blank; every condition column must name a \
                 factor value[...] column."
            )));
        }
        if name.contains(['\t', '\n', '\r']) {
            return Err(MzLibError::Usage(format!(
                "Condition column {index} contains a tab or line break, which the wire uses \
                 to separate names: {name:?}."
            )));
        }
    }
    args.push("--condition-columns".to_owned());
    args.push(condition_columns.join("\t"));
    Ok(())
}

/// `--searched-files-stdin` and its stdin for [`design_with`] and [`isobaric_design_with`]:
/// nothing for `None`. An empty list is sent, and mzLib refuses it as a design.
fn searched_files_stdin(
    args: &mut Vec<String>,
    searched_files: Option<&Vec<PathBuf>>,
) -> Result<Option<String>> {
    let Some(files) = searched_files else {
        return Ok(None);
    };
    let lines = if files.is_empty() {
        Vec::new()
    } else {
        bridge::path_lines(files, "searched file")?
    };
    args.push("--searched-files-stdin".to_owned());
    Ok(Some(lines.join("\n") + "\n"))
}

// ---- isobaric design (mzLib 1.0.595) ------------------------------------------------------------

/// The first pyMzLib release whose bridge has `sdrf isobaric-design` and `sdrf restrict`: their
/// specs' `since.pymzlib`.
const ISOBARIC_RESTRICT_SINCE: &str = "0.6.0";

/// How [`isobaric_design_with`] reads a design. Every default is mzLib's, and mzLib's default for
/// the kit and the plex source is to refuse: neither is ever guessed.
///
/// Declare the plex with **exactly one** of [`Self::plex_column`], [`Self::plex_file_name_pattern`]
/// and [`Self::single_plex`]. None, or several, is mzLib's refusal, not a usage error.
#[derive(Debug, Clone)]
pub struct IsobaricDesignOptions {
    /// The kit the search used, named as [`crate::isobaric::kits`] takes it (`"TMT10"`, `"tmt18"`,
    /// `"TMT6-plex"`, `"TMT11 on K"`), resolved by the same mzLib code. `None` (the default) sends
    /// none, and mzLib then refuses: the kit is never guessed from the labels, because one kit's
    /// labels can be a subset of another's. A blank name is refused here; an unknown one is the
    /// bridge's usage error listing the known kits.
    pub kit: Option<String>,
    /// The column each row's plex is read from, by its exact name.
    pub plex_column: Option<String>,
    /// A .NET regular expression matched, ignoring case, against each row's file name; the plex is
    /// its group named `plex`, else its first group, e.g. `r"TMT_?pool(?<plex>\d+)"`.
    pub plex_file_name_pattern: Option<String>,
    /// The name of the one plex every file belongs to.
    pub single_plex: Option<String>,
    /// The `factor value[...]` columns the condition is built from, exactly as
    /// [`DesignOptions::condition_columns`]: empty (the default) sends nothing, and mzLib uses the
    /// document's only factor column.
    pub condition_columns: Vec<String>,
    /// The files the search will read, exactly as [`DesignOptions::searched_files`]: each must be
    /// named **exactly** by some row. A search that read `X-calib.mzML` for an SDRF naming `X.raw`
    /// matches nothing; [`restrict_with`] the SDRF first.
    pub searched_files: Option<Vec<PathBuf>>,
    /// A **folder**, created when missing. mzLib writes MetaMorpheus's `TmtDesign.txt` (1-based)
    /// inside it, replacing one already there, but only when the design is valid. MetaMorpheus
    /// reads it only from the folder of the first spectra file. A path that names an existing
    /// file is the bridge's usage error, checked before the SDRF is read.
    pub out: Option<PathBuf>,
    /// The longest the call may take. `None` waits indefinitely. Defaults to 60 seconds, as for
    /// every function in this module.
    pub timeout: Option<Duration>,
}

impl Default for IsobaricDesignOptions {
    fn default() -> Self {
        Self {
            kit: None,
            plex_column: None,
            plex_file_name_pattern: None,
            single_plex: None,
            condition_columns: Vec::new(),
            searched_files: None,
            out: None,
            timeout: Some(DEFAULT_TIMEOUT),
        }
    }
}

/// The `TmtDesign.txt` that [`isobaric_design_with`] wrote under [`IsobaricDesignOptions::out`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct IsobaricDesignWritten {
    /// The file written: the `out` folder joined with `TmtDesign.txt`.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// Files written, in files.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Rows written, one per file and annotated channel, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
}

/// One annotated channel of one file in an isobaric design, as mzLib reads it: a row of
/// [`SdrfIsobaricDesign::channels`].
///
/// **Every number here is 1-based, as the SDRF writes it**, unlike [`DesignedFile`], whose
/// numbers are 0-based. That is mzLib's own split (`SdrfIsobaricDesign` against
/// `SpectraFileInfo`), and this crate converts neither.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DesignedChannel {
    /// The file: the searched path you passed for it, else the file name the SDRF gives.
    pub full_path: String,
    /// The plex, as the source you declared gives it (a column cell, a pattern's group, or the
    /// single plex).
    pub plex: String,
    /// The plex's number in [`SdrfIsobaricDesign::plex_ids`], 1-based.
    pub plex_id: u32,
    /// The SDRF `source name` of the sample in this channel.
    pub sample_name: String,
    /// The channel as the kit spells it, e.g. `127N`, whatever form `comment[label]` used.
    pub channel: String,
    /// The condition columns' values joined with `_`.
    pub condition: String,
    /// 1-based, exactly as written; never ranked.
    pub biological_replicate: u32,
    /// 1-based, as written; 1 when the SDRF has no fraction column.
    pub fraction: u32,
    /// 1-based, as written; 1 when the SDRF has no technical-replicate column.
    pub technical_replicate: u32,
}

/// An isobaric (TMT, TMTpro, iTRAQ, DiLeu) design read from an SDRF, or every reason it was
/// refused: what [`isobaric_design_with`] returns.
///
/// **A refusal is an answer, not an error.** When [`Self::is_valid`] is `false`,
/// [`Self::refusals`] lists every reason at once and the table is empty. A missing kit, no plex
/// source, two plex sources and a bad pattern are refusals too.
///
/// The table has one row per file and annotated channel: `full_path`, `plex`, `plex_id`,
/// `sample_name`, `channel`, `condition`, `biological_replicate`, `fraction` and
/// `technical_replicate`, with the meanings of [`DesignedChannel`]. [`Self::channels`] gives the
/// same rows as structs. **All of the numbers are 1-based.**
#[derive(Debug, Clone, Deserialize)]
pub struct SdrfIsobaricDesign {
    /// The SDRF path as given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub path: String,
    /// `true` when there are no refusals. `false` is an answer: the table is then empty.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub is_valid: bool,
    /// The kit name exactly as you gave it. `None` when you gave none, and mzLib then refused.
    #[serde(default)]
    pub kit: Option<String>,
    /// mzLib's name for the kit it resolved to (`TMT10`, `TMT18`, `iTRAQ4`, ...), as
    /// [`crate::isobaric::kits`] names it. `None` when no kit was given.
    #[serde(default)]
    pub kit_type: Option<String>,
    /// Every channel of the kit, in reporter m/z order. A channel here and absent from the table
    /// is one the SDRF does not annotate. `None` when no kit was given.
    #[serde(default)]
    pub kit_channels: Option<Vec<String>>,
    /// The plex column you gave, or `None`.
    #[serde(default)]
    pub plex_column: Option<String>,
    /// The plex pattern you gave, or `None`.
    #[serde(default)]
    pub plex_file_name_pattern: Option<String>,
    /// The single plex you gave, or `None`.
    #[serde(default)]
    pub single_plex: Option<String>,
    /// Where mzLib read the plex from, as its report prints it; `"(not declared)"` when no single
    /// source could be used.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub plex_source: String,
    /// Plex name to its number, 1..N in case-insensitive name order (the order MetaMorpheus
    /// numbers plexes in). The map iterates in this crate's name order; the number is mzLib's.
    /// Empty when refused.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub plex_ids: BTreeMap<String, u32>,
    /// Plexes in the design, a count of plexes; `0` when refused.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub plex_count: u64,
    /// The column files were keyed on: `comment[searched data file]` when the SDRF has it, else
    /// `comment[data file]`. `None` when it has neither, so no row names a file; that is itself a
    /// refusal.
    #[serde(default)]
    pub file_key_column: Option<String>,
    /// The `factor value[...]` columns the condition was built from, in join order; empty when
    /// none could be chosen.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub condition_columns: Vec<String>,
    /// Whether [`IsobaricDesignOptions::condition_columns`] was given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub condition_columns_declared: bool,
    /// Whether [`IsobaricDesignOptions::searched_files`] was given.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub searched_files_given: bool,
    /// Files in the design, in files; `0` when refused.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub file_count: u64,
    /// Rows in the table, one per file and annotated channel, in rows; `0` when refused.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// Every reason the design was refused, in mzLib's words. Empty when valid.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub refusals: Vec<String>,
    /// What was dropped or filled on the way: rows for files the search does not read, and a file
    /// with no row of its own for one of its plex's channels.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub notes: Vec<String>,
    /// mzLib's human-readable account of all of the above (`SdrfIsobaricDesign.Report`). Print it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub report: String,
    /// The design table: files in SDRF row order, then each file's channels in reporter m/z order.
    #[serde(flatten)]
    pub columns: Table,
    /// The `TmtDesign.txt` written; `None` when [`IsobaricDesignOptions::out`] was not given, or
    /// the design was refused and nothing was written.
    #[serde(default)]
    pub written: Option<IsobaricDesignWritten>,
    /// What the design does and does not cover. The first says it is 1-based; read them once.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl SdrfIsobaricDesign {
    /// Every row as a [`DesignedChannel`]: files in SDRF order, each file's channels in reporter
    /// m/z order. Empty when refused.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is missing or not the type the wire contract says.
    pub fn channels(&self) -> Result<Vec<DesignedChannel>> {
        if self.columns.is_empty() {
            return Ok(Vec::new());
        }
        let text = |name: &str| self.columns.strings(name);
        let (full_path, plex, sample_name, channel, condition) = (
            text("full_path")?,
            text("plex")?,
            text("sample_name")?,
            text("channel")?,
            text("condition")?,
        );
        let plex_id = self.columns.integers("plex_id")?;
        let biological = self.columns.integers("biological_replicate")?;
        let fraction = self.columns.integers("fraction")?;
        let technical = self.columns.integers("technical_replicate")?;
        let number = |value: Option<i64>, what: &str| -> Result<u32> {
            value.and_then(|v| u32::try_from(v).ok()).ok_or_else(|| {
                MzLibError::Protocol(format!(
                    "sdrf isobaric-design returned a {what} that is not a whole number: {value:?}"
                ))
            })
        };
        (0..full_path.len())
            .map(|i| {
                Ok(DesignedChannel {
                    full_path: full_path[i].clone().unwrap_or_default(),
                    plex: plex[i].clone().unwrap_or_default(),
                    plex_id: number(plex_id[i], "plex_id")?,
                    sample_name: sample_name[i].clone().unwrap_or_default(),
                    channel: channel[i].clone().unwrap_or_default(),
                    condition: condition[i].clone().unwrap_or_default(),
                    biological_replicate: number(biological[i], "biological_replicate")?,
                    fraction: number(fraction[i], "fraction")?,
                    technical_replicate: number(technical[i], "technical_replicate")?,
                })
            })
            .collect()
    }
}

/// `--<name> value` for a string option, nothing for `None`; a blank value or a line break is
/// refused before anything is spawned.
fn optional_text(
    args: &mut Vec<String>,
    value: Option<&String>,
    name: &str,
    example: &str,
) -> Result<()> {
    let Some(value) = value else {
        return Ok(());
    };
    if value.trim().is_empty() {
        return Err(MzLibError::Usage(format!(
            "{name} must be non-blank, or None, e.g. {example:?}."
        )));
    }
    if value.contains(['\n', '\r']) {
        return Err(MzLibError::Usage(format!(
            "{name} may not contain a line break: {value:?}."
        )));
    }
    args.push(format!("--{}", name.replace('_', "-")));
    args.push(value.clone());
    Ok(())
}

/// The arguments and stdin of a `sdrf isobaric-design` call, validated before anything is spawned.
fn isobaric_design_request(
    path: &Path,
    options: &IsobaricDesignOptions,
) -> Result<(Vec<String>, Option<String>)> {
    let mut args = one_path_args("isobaric-design", path)?;
    optional_text(&mut args, options.kit.as_ref(), "kit", "TMT10")?;
    optional_text(
        &mut args,
        options.plex_column.as_ref(),
        "plex_column",
        "comment[plex]",
    )?;
    optional_text(
        &mut args,
        options.plex_file_name_pattern.as_ref(),
        "plex_file_name_pattern",
        r"TMT_?pool(?<plex>\d+)",
    )?;
    optional_text(
        &mut args,
        options.single_plex.as_ref(),
        "single_plex",
        "plex1",
    )?;
    condition_columns_args(&mut args, &options.condition_columns)?;
    let stdin = searched_files_stdin(&mut args, options.searched_files.as_ref())?;
    if let Some(out) = &options.out {
        let out = path_text(out, "out path")?;
        if out.is_empty() {
            return Err(MzLibError::Usage(
                "out must name a folder; TmtDesign.txt is written inside it, e.g. 'design'."
                    .to_owned(),
            ));
        }
        args.push("--out".to_owned());
        args.push(out.to_owned());
    }
    Ok((args, stdin))
}

/// Read an isobaric (TMT, TMTpro, iTRAQ, DiLeu) design out of an SDRF, or every reason it cannot
/// be read.
///
/// Calls mzLib's `SdrfIsobaricDesign.Read` (mzLib #1380). It runs MetaMorpheus's own checks on a
/// TMT design (one plex, fraction and technical replicate per file; a channel described once per
/// plex; no sample twice at one fraction in a plex) and **refuses rather than repairs**, reporting
/// every failure at once. **It guesses nothing:** not the kit, because one kit's labels can be a
/// subset of another's, and not the plex, which you declare with exactly one of
/// [`IsobaricDesignOptions::plex_column`], [`IsobaricDesignOptions::plex_file_name_pattern`] and
/// [`IsobaricDesignOptions::single_plex`]. Giving none, or several, is mzLib's refusal.
///
/// **The numbers are 1-based, as the SDRF writes them.** [`design_with`]'s are 0-based. That split
/// is mzLib's, and this crate converts neither. The
/// [module guide](crate::sdrf#read-an-isobaric-tmt-design) walks through two real deposits.
#[doc = include_str!("../docs/reference/sdrf.isobaric-design.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for a blank path, a blank kit, plex option or
/// condition column, one holding a line break (a tab, for a condition column), a searched file
/// that is blank or holds a tab or line break, or a blank `out`; and when the bridge in use is
/// older than pyMzLib 0.6.0's, which has no `sdrf isobaric-design`.
///
/// # Examples
///
/// PXD008841 writes no plex column; the plex survives only in the file names, spelled two ways
/// (`TMTpool1`, `TMT_pool5`), so one pattern reads both:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{isobaric_design_with, IsobaricDesignOptions};
///
/// let d = isobaric_design_with(
///     "PXD008841.subset.sdrf.tsv",
///     &IsobaricDesignOptions {
///         kit: Some("TMT10".into()),
///         plex_file_name_pattern: Some(r"TMT_?pool(?<plex>\d+)".into()),
///         ..Default::default()
///     },
/// )?;
/// assert_eq!((d.is_valid, d.file_count, d.record_count), (true, 4, 40));
/// assert_eq!(d.plex_ids, [("1".to_owned(), 1), ("5".to_owned(), 2)].into());
/// let c = &d.channels()?[9];
/// assert_eq!((c.channel.as_str(), c.sample_name.as_str()), ("131N", "pool"));
/// assert_eq!((c.plex_id, c.fraction, c.biological_replicate), (1, 1, 1));   // 1-based
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
///
/// Declare no plex source and mzLib refuses rather than guess one:
///
/// ```
/// # mzlib_replay::activate();
/// # use mzlib::sdrf::{isobaric_design_with, IsobaricDesignOptions};
/// let refused = isobaric_design_with(
///     "PXD008841.subset.sdrf.tsv",
///     &IsobaricDesignOptions { kit: Some("TMT10".into()), ..Default::default() },
/// )?;
/// assert_eq!((refused.is_valid, refused.record_count), (false, 0));
/// assert!(refused.refusals[0].starts_with("No plex source was declared."));
/// assert!(refused.channels()?.is_empty());
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.isobaric-design.see-also.md")]
pub fn isobaric_design_with(
    path: impl AsRef<Path>,
    options: &IsobaricDesignOptions,
) -> Result<SdrfIsobaricDesign> {
    let (args, stdin) = isobaric_design_request(path.as_ref(), options)?;
    bridge::require_verb("sdrf isobaric-design", ISOBARIC_RESTRICT_SINCE)?;
    call(&args, stdin.as_deref(), options.timeout)
}

// ---- restrict (mzLib 1.0.595) -------------------------------------------------------------------

/// The column [`restrict_with`] writes each kept row's searched file name into, directly after
/// `comment[data file]`.
pub const SEARCHED_DATA_FILE_COLUMN: &str = "comment[searched data file]";

/// One file a search read, for [`restrict_with`]: its name, and the acquired file it came from
/// when you know it.
///
/// Build one from a path or name (`"X-calib.mzML".into()`), or with [`Self::acquired_from`] when
/// you know the acquisition, as MetaMorpheus does from an mzML's own `sourceFile`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SearchedFile {
    /// The file the search read, as a path or a bare name.
    pub searched: PathBuf,
    /// The acquired file it came from, when you know it. It overrides mzLib's suffix rule for this
    /// file. `None` (what `From` gives) leaves the join to the rule.
    pub acquired: Option<PathBuf>,
}

impl SearchedFile {
    /// A searched file whose acquired file you know.
    pub fn acquired_from(searched: impl Into<PathBuf>, acquired: impl Into<PathBuf>) -> Self {
        Self {
            searched: searched.into(),
            acquired: Some(acquired.into()),
        }
    }
}

impl From<&str> for SearchedFile {
    fn from(searched: &str) -> Self {
        Self {
            searched: searched.into(),
            acquired: None,
        }
    }
}

impl From<String> for SearchedFile {
    fn from(searched: String) -> Self {
        Self {
            searched: searched.into(),
            acquired: None,
        }
    }
}

impl From<PathBuf> for SearchedFile {
    fn from(searched: PathBuf) -> Self {
        Self {
            searched,
            acquired: None,
        }
    }
}

impl From<&Path> for SearchedFile {
    fn from(searched: &Path) -> Self {
        Self {
            searched: searched.to_path_buf(),
            acquired: None,
        }
    }
}

/// How much of a restricted document to return, and where to write the whole of it. Every
/// default is mzLib's.
#[derive(Debug, Clone)]
pub struct RestrictOptions {
    /// Return at most this many restricted rows. `None` (the default) returns every row.
    pub limit: Option<u64>,
    /// Skip this many restricted rows in the payload. The file [`Self::out`] writes is never
    /// windowed.
    pub offset: u64,
    /// Also write the **whole** restricted document here as SDRF, whatever the limit and offset
    /// are, with CRLF line endings. Must end in `.tsv`, checked by the bridge before the SDRF is
    /// read; a missing folder is created. **Name it `.sdrf.tsv`**: mzLib adds `.sdrf.tsv` to a
    /// name that does not already end in it, so `searched.tsv` is written as
    /// `searched.tsv.sdrf.tsv`. [`RestrictedSdrf::written`] names the file written.
    pub out: Option<PathBuf>,
    /// The longest the call may take. `None` waits indefinitely. Defaults to 60 seconds.
    pub timeout: Option<Duration>,
}

impl Default for RestrictOptions {
    fn default() -> Self {
        Self {
            limit: None,
            offset: 0,
            out: None,
            timeout: Some(DEFAULT_TIMEOUT),
        }
    }
}

/// An SDRF restricted to the files one search read, and what the restriction could not place:
/// what [`restrict_with`] returns.
///
/// [`Self::document`] is the **restricted** document, with every [`SdrfDocument`] accessor: the
/// input's rows for the acquisitions the search read, carried cell for cell, with
/// [`SEARCHED_DATA_FILE_COLUMN`] added directly after `comment[data file]`. A kept row shorter
/// than the header is padded with empty cells, so unlike [`read`] the rows are never ragged.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct RestrictedSdrf {
    /// The restricted document. Its `path` is the SDRF that was restricted, as given; its
    /// `row_count` counts the rows kept, before the window.
    #[serde(flatten)]
    pub document: SdrfDocument,
    /// The searched files as sent, in order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub searched_files: Vec<String>,
    /// How many searched files you named an acquired file for, in files. `0` when none was, and
    /// mzLib was then given no map.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub acquired_name_count: u64,
    /// Data rows in the SDRF before the restriction, in rows.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub input_row_count: u64,
    /// Searched files (file names) the SDRF has no row for. Reported, never invented.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub searched_without_row: Vec<String>,
    /// Acquired files the SDRF describes that this search did not read, in row order. A row
    /// naming no file is dropped too, and not listed.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub dropped_data_files: Vec<String>,
    /// Acquisitions two or more searched files claimed, as `"stem: searched, searched"`. The first
    /// in name order kept the rows; the others have none and are not in
    /// [`Self::searched_without_row`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub ambiguous: Vec<String>,
    /// Searched files whose rows name two or more **different** acquired files (one run listed as
    /// `X.raw` and `X.mzML`), as `"searched: acquired, acquired"`. Every such row is kept, so the
    /// document can describe one run twice.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub shared_stem: Vec<String>,
    /// Where the whole restricted document was written, or `None` if [`RestrictOptions::out`] was
    /// not set.
    #[serde(default)]
    pub written: Option<WrittenSdrf>,
}

impl RestrictedSdrf {
    /// The searched file each returned row was joined to: its name, no folder.
    pub fn searched_data_files(&self) -> Vec<Option<&str>> {
        self.document.value(SEARCHED_DATA_FILE_COLUMN)
    }
}

/// The arguments and stdin of a `sdrf restrict` call, validated before anything is spawned.
///
/// One stdin line per searched file: `searched`, or `searched<TAB>acquired` where an acquired file
/// is named.
fn restrict_request(
    path: &Path,
    searched_files: &[SearchedFile],
    options: &RestrictOptions,
) -> Result<(Vec<String>, String)> {
    if searched_files.is_empty() {
        return Err(MzLibError::Usage(
            "restrict needs at least one searched file.".to_owned(),
        ));
    }
    let mut lines = Vec::with_capacity(searched_files.len());
    for (index, file) in searched_files.iter().enumerate() {
        let searched = path_text(&file.searched, "searched file")?;
        if searched.is_empty() {
            return Err(MzLibError::Usage(format!(
                "Searched file {index} is blank; every searched file must be a path or a name."
            )));
        }
        let mut fields = vec![searched];
        if let Some(acquired) = &file.acquired {
            let acquired = path_text(acquired, "acquired file")?;
            if acquired.is_empty() {
                return Err(MzLibError::Usage(format!(
                    "The acquired file for '{searched}' is blank; give a name, or None."
                )));
            }
            fields.push(acquired);
        }
        if let Some(bad) = fields.iter().find(|text| text.contains(['\t', '\n', '\r'])) {
            return Err(MzLibError::Usage(format!(
                "A searched or acquired file contains a tab or newline, which the wire uses to \
                 separate them: {bad:?}."
            )));
        }
        lines.push(fields.join("\t"));
    }

    let mut args = one_path_args("restrict", path)?;
    window_args(
        &mut args,
        &ReadOptions {
            limit: options.limit,
            offset: options.offset,
            timeout: None,
        },
    );
    if let Some(out) = &options.out {
        let out = path_text(out, "out path")?;
        if out.is_empty() {
            return Err(MzLibError::Usage(
                "out must be a file path ending in .tsv, e.g. 'searched.sdrf.tsv'.".to_owned(),
            ));
        }
        args.push("--out".to_owned());
        args.push(out.to_owned());
    }
    Ok((args, lines.join("\n") + "\n"))
}

/// Restrict an SDRF to the files one search read, with every default.
///
/// See [`restrict_with`] for the reference.
///
/// # Errors
///
/// As [`restrict_with`].
pub fn restrict(path: impl AsRef<Path>, searched_files: &[SearchedFile]) -> Result<RestrictedSdrf> {
    restrict_with(path, searched_files, &RestrictOptions::default())
}

/// Restrict an SDRF to the files one search read, joining each searched file to the row of its
/// acquisition, and report what the join could not place.
///
/// Calls mzLib's `SdrfSearchScope.Restrict` (mzLib #1451). A search reads `X-calib.mzML` while the
/// SDRF names the acquired `X.raw`; mzLib joins them **by stem**, ignoring case, after removing
/// the derivative suffixes MetaMorpheus writes (`-calib`, `-averaged`, in any order and number).
/// Each kept row is carried cell for cell, with the searched file's name in a new
/// [`SEARCHED_DATA_FILE_COLUMN`]; nothing is renumbered. A [`SearchedFile::acquired`] you name
/// overrides the suffix rule for that file.
///
/// **Restrict first, then read a design.** [`design_with`] and [`isobaric_design_with`] match
/// `searched_files` by exact name, and key on `comment[searched data file]` when a document has
/// it, so a restricted document joins the search's own file names. The
/// [module guide](crate::sdrf#restrict-an-sdrf-to-one-searchs-files) shows both steps.
#[doc = include_str!("../docs/reference/sdrf.restrict.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] before anything is spawned for a blank path, no searched file, a blank
/// searched or acquired file, one holding a tab or line break, or a blank `out`; and when the
/// bridge in use is older than pyMzLib 0.6.0's, which has no `sdrf restrict`.
///
/// # Examples
///
/// Two calibrated runs, the raw file of one of them as well, a renamed file whose acquisition you
/// know, and one the SDRF does not describe:
///
/// ```
/// # mzlib_replay::activate();
/// use mzlib::sdrf::{restrict_with, RestrictOptions, SearchedFile};
///
/// let pool1 = "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0";
/// let pool5 = "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0";
/// let searched = [
///     SearchedFile::from(format!(r"D:\search\{pool1}1-calib.mzML")),
///     SearchedFile::from(format!("{pool1}2-calib.mzML")),
///     SearchedFile::from(format!("{pool1}2.raw")),
///     SearchedFile::acquired_from("sample_A.mzML", format!("{pool5}1.raw")),
///     SearchedFile::from("not_in_sdrf.mzML"),
/// ];
/// let scoped = restrict_with(
///     "PXD008841.subset.sdrf.tsv",
///     &searched,
///     &RestrictOptions { limit: Some(2), ..Default::default() },
/// )?;
/// let doc = &scoped.document;
/// assert_eq!((scoped.input_row_count, doc.row_count, doc.returned_count), (40, 30, 2));
/// let first = format!("{pool1}1-calib.mzML");
/// assert_eq!(scoped.searched_data_files(), [Some(first.as_str()); 2]);   // the name, no folder
/// assert_eq!(scoped.searched_without_row, ["not_in_sdrf.mzML"]);
/// assert_eq!(scoped.dropped_data_files, [format!("{pool5}2.raw")]);
/// assert_eq!(
///     scoped.ambiguous,
///     [format!("{pool1}2: {pool1}2-calib.mzML, {pool1}2.raw")]
/// );
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/sdrf.restrict.see-also.md")]
pub fn restrict_with(
    path: impl AsRef<Path>,
    searched_files: &[SearchedFile],
    options: &RestrictOptions,
) -> Result<RestrictedSdrf> {
    let (args, stdin) = restrict_request(path.as_ref(), searched_files, options)?;
    bridge::require_verb("sdrf restrict", ISOBARIC_RESTRICT_SINCE)?;
    call(&args, Some(&stdin), options.timeout)
}

fn protocol(error: serde_json::Error) -> MzLibError {
    MzLibError::Protocol(format!("sdrf payload could not be interpreted: {error}"))
}

#[cfg(test)]
mod tests {
    //! Offline. The fixtures are recorded from the real bridge by pyMzLib and shared verbatim, so
    //! they carry the shape mzLib actually produces — PXD059974's raggedness and PXD000070's eight
    //! repetitions of `comment[modification parameters]` are the point, not an accident.

    use super::*;

    const READ_FIXTURE: &str = include_str!("../tests/fixtures/sdrf_read_PXD000070.json");
    const RAGGED_FIXTURE: &str = include_str!("../tests/fixtures/sdrf_read_ragged.json");
    const POOL_FIXTURE: &str = include_str!("../tests/fixtures/sdrf_pool_two.json");

    fn document() -> SdrfDocument {
        serde_json::from_str(READ_FIXTURE).unwrap()
    }

    fn ragged() -> SdrfDocument {
        serde_json::from_str(RAGGED_FIXTURE).unwrap()
    }

    fn pooled() -> PooledSdrf {
        serde_json::from_str(POOL_FIXTURE).unwrap()
    }

    fn labelled(pairs: &[(&str, &str)]) -> PoolInput {
        PoolInput::Labelled(
            pairs
                .iter()
                .map(|(p, l)| (PathBuf::from(p), (*l).to_owned()))
                .collect(),
        )
    }

    // ---- the shape ----------------------------------------------------------------------------

    #[test]
    fn columns_is_a_vec_because_sdrf_names_repeat() {
        // The whole reason this module is row-major. A map keyed by name would keep one of eight.
        let doc = document();
        assert!(doc.has_repeated_columns());
        assert_eq!(doc.indexes_of("comment[modification parameters]").len(), 8);
    }

    #[test]
    fn all_returns_every_occurrence_where_value_returns_the_first() {
        let doc = document();
        let every = &doc.all("comment[modification parameters]")[0];
        let first = doc.value("comment[modification parameters]")[0].unwrap();

        assert_eq!(every.len(), 8);
        assert_eq!(every[0], first);
        assert!(first.starts_with("NT=Carbamidomethyl"));
    }

    #[test]
    fn a_cv_cell_crosses_intact_rather_than_being_split_on_its_semicolons() {
        // The defect this module exists to close: read_records joins an SdrfRow's cells with ";".
        let doc = document();
        assert_eq!(
            doc.value("comment[modification parameters]")[0],
            Some("NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed")
        );
    }

    #[test]
    fn an_absent_column_is_none_and_a_reserved_word_is_itself() {
        // None means "no such column"; "not applicable" is a real answer a curator wrote.
        let doc = document();
        assert_eq!(
            doc.value("characteristics[nonesuch]"),
            vec![None; doc.rows.len()]
        );
        assert_eq!(doc.index_of("characteristics[nonesuch]"), None);
        assert_eq!(
            doc.value("characteristics[disease]")[0],
            Some("not applicable")
        );
    }

    #[test]
    fn a_short_row_reports_none_rather_than_shifting_later_columns() {
        let doc = ragged();
        assert_eq!(doc.ragged_row_count(), 17);
        let short = doc
            .rows
            .iter()
            .position(|row| row.len() < doc.columns.len())
            .unwrap();
        let full = doc
            .rows
            .iter()
            .position(|row| row.len() == doc.columns.len())
            .unwrap();

        let last = doc.columns.last().unwrap();
        assert!(doc.value(last)[full].is_some());
        assert!(doc.value(last)[short].is_none());
        // A position both rows reach still lines up, which is what "no shifting" means.
        assert!(doc.value("characteristics[organism]")[short].is_some());
    }

    #[test]
    fn caveats_name_the_raggedness_rather_than_leaving_it_to_be_discovered() {
        assert!(ragged().caveats.iter().any(|c| c.contains("SHORT")));
    }

    #[test]
    fn truncated_says_rows_were_left_behind() {
        // Recorded with limit=4 against a 24-row merge.
        let pooled = pooled();
        assert!(pooled.document.truncated);
        assert_eq!(pooled.document.returned_count, 4);
        assert_eq!(pooled.document.row_count, 24);
    }

    #[test]
    fn a_complete_read_is_not_marked_truncated() {
        let doc = ragged();
        assert!(!doc.truncated);
        assert_eq!(doc.returned_count, doc.row_count);
    }

    #[test]
    fn records_is_offered_but_lossy_when_names_repeat() {
        let doc = document();
        assert!(doc.has_repeated_columns());
        assert!(doc.records()[0].len() < doc.columns.len());
    }

    #[test]
    fn a_null_list_on_the_wire_is_empty_rather_than_an_error() {
        let doc: SdrfDocument = serde_json::from_str(
            r#"{"path": "x", "column_names": null, "rows": null, "caveats": null}"#,
        )
        .unwrap();
        assert!(doc.columns.is_empty() && doc.rows.is_empty() && doc.caveats.is_empty());
    }

    // ---- pool ---------------------------------------------------------------------------------

    #[test]
    fn pool_carries_provenance_and_the_labels_it_was_given() {
        let pooled = pooled();
        assert_eq!(pooled.document_count, 2);
        assert_eq!(pooled.labels, ["malaria", "colon"]);
        assert!(pooled.written.is_none());
        assert!(pooled
            .document
            .columns
            .iter()
            .any(|c| c == SOURCE_DOCUMENT_COLUMN));
        assert!(pooled
            .source_documents()
            .iter()
            .all(|s| matches!(s, Some("malaria" | "colon"))));
    }

    #[test]
    fn pool_renders_one_tab_separated_stdin_line_per_document() {
        let (args, stdin) = pool_request(
            &labelled(&[("a.sdrf.tsv", "malaria"), ("b.sdrf.tsv", "colon")]),
            &PoolOptions::default(),
        )
        .unwrap();
        assert_eq!(stdin, "a.sdrf.tsv\tmalaria\nb.sdrf.tsv\tcolon");
        assert_eq!(args, ["sdrf", "pool"]);
    }

    #[test]
    fn pool_without_labels_sends_bare_paths() {
        let (_, stdin) = pool_request(
            &PoolInput::Paths(vec!["a.sdrf.tsv".into(), "b.sdrf.tsv".into()]),
            &PoolOptions::default(),
        )
        .unwrap();
        assert_eq!(stdin, "a.sdrf.tsv\nb.sdrf.tsv");
    }

    #[test]
    fn pool_passes_the_fallback_warning_through() {
        // The recorded fixture was pooled WITH labels, so it carries no warning; the bridge adds
        // one when labels are absent, and it must reach the caller.
        assert!(!pooled()
            .document
            .caveats
            .iter()
            .any(|c| c.contains("reproducible")));

        let mut payload: serde_json::Value = serde_json::from_str(POOL_FIXTURE).unwrap();
        payload["caveats"]
            .as_array_mut()
            .unwrap()
            .push("... not reproducible on another machine ...".into());
        let unlabelled: PooledSdrf = serde_json::from_value(payload).unwrap();
        assert!(unlabelled
            .document
            .caveats
            .iter()
            .any(|c| c.contains("reproducible")));
    }

    #[test]
    fn pool_writes_the_whole_document_and_reports_it() {
        let (args, _) = pool_request(
            &PoolInput::Paths(vec!["a.sdrf.tsv".into()]),
            &PoolOptions {
                read: ReadOptions {
                    limit: Some(4),
                    ..ReadOptions::default()
                },
                out: Some("merged.sdrf.tsv".to_owned()),
            },
        )
        .unwrap();
        assert_eq!(
            args,
            ["sdrf", "pool", "--limit", "4", "--out", "merged.sdrf.tsv"]
        );

        let mut payload: serde_json::Value = serde_json::from_str(POOL_FIXTURE).unwrap();
        payload["written"] = serde_json::json!({"path": "merged.sdrf.tsv", "row_count": 24});
        let written: PooledSdrf = serde_json::from_value(payload).unwrap();
        assert_eq!(
            written.written,
            Some(WrittenSdrf {
                path: "merged.sdrf.tsv".to_owned(),
                row_count: 24
            })
        );
    }

    // ---- argument validation, before anything is spawned --------------------------------------

    #[test]
    fn an_empty_selection_is_refused() {
        let error = pool_request(&PoolInput::Paths(vec![]), &PoolOptions::default()).unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(error.to_string().contains("At least one"), "{error}");
    }

    #[test]
    fn a_blank_label_is_refused_rather_than_silently_defaulted() {
        let error = pool_request(
            &labelled(&[("a.sdrf.tsv", "malaria"), ("b.sdrf.tsv", "   ")]),
            &PoolOptions::default(),
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert!(
            error.to_string().contains("label for 'b.sdrf.tsv'"),
            "{error}"
        );
    }

    #[test]
    fn a_tab_in_a_label_is_refused_because_it_is_the_field_separator() {
        let error = pool_request(
            &labelled(&[("a.sdrf.tsv", "mal\taria")]),
            &PoolOptions::default(),
        )
        .unwrap_err();
        assert!(error.to_string().contains("tab or newline"), "{error}");
    }

    #[test]
    fn a_blank_path_is_refused() {
        for path in ["", "   "] {
            let error = read_args(Path::new(path), &ReadOptions::default()).unwrap_err();
            assert!(matches!(error, MzLibError::Usage(_)));
            assert!(
                error.to_string().contains("file path is required"),
                "{error}"
            );
        }
        let error =
            pool_request(&PoolInput::Paths(vec![" ".into()]), &PoolOptions::default()).unwrap_err();
        assert!(error.to_string().contains("may not be blank"), "{error}");
    }

    #[test]
    fn an_empty_out_is_refused() {
        let error = pool_request(
            &PoolInput::Paths(vec!["a.sdrf.tsv".into()]),
            &PoolOptions {
                out: Some("  ".to_owned()),
                ..PoolOptions::default()
            },
        )
        .unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
    }

    #[test]
    fn the_window_options_reach_the_bridge() {
        let args = read_args(
            Path::new("a.sdrf.tsv"),
            &ReadOptions {
                limit: Some(2),
                offset: 1,
                ..ReadOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "sdrf",
                "read",
                "--path",
                "a.sdrf.tsv",
                "--limit",
                "2",
                "--offset",
                "1"
            ]
        );
    }

    #[test]
    fn a_zero_limit_is_sent_and_a_zero_offset_is_not() {
        // The bridge accepts --limit 0 (a header-only answer); a default offset sent explicitly is
        // a default the bridge could later disagree with.
        let args = read_args(
            Path::new("a.sdrf.tsv"),
            &ReadOptions {
                limit: Some(0),
                ..ReadOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            args,
            ["sdrf", "read", "--path", "a.sdrf.tsv", "--limit", "0"]
        );
    }

    #[test]
    fn the_default_timeout_matches_pymzlib() {
        assert_eq!(
            ReadOptions::default().timeout,
            Some(Duration::from_secs(60))
        );
    }

    // ---- validate, lint, assess, samples, parse_ages ------------------------------------------

    fn fixture<T: serde::de::DeserializeOwned>(text: &str) -> T {
        serde_json::from_str(text).unwrap()
    }

    #[test]
    fn a_validation_reads_into_typed_findings() {
        let cohort: SdrfValidation =
            fixture(include_str!("../tests/fixtures/sdrf_validate_cohort.json"));
        assert!(cohort.is_valid);
        assert_eq!((cohort.error_count, cohort.warning_count), (0, 2));
        assert_eq!(cohort.columns.rows() as u64, cohort.message_count);
        let warnings = cohort.warnings().unwrap();
        assert_eq!(warnings.len(), 2);
        assert!(warnings.iter().all(|m| m.rule == "ReservedWordCase"));
        // line_number is row_index + 2: the header is line 1.
        for m in &warnings {
            assert_eq!(m.line_number, m.row_index.map(|r| r + 2));
        }
        assert!(cohort.errors().unwrap().is_empty());
    }

    #[test]
    fn a_document_level_finding_has_no_row_and_no_line() {
        let skeleton: SdrfValidation = fixture(include_str!(
            "../tests/fixtures/sdrf_validate_skeleton.json"
        ));
        assert!(!skeleton.is_valid);
        assert_eq!(
            skeleton.errors().unwrap().len() as u64,
            skeleton.error_count
        );
        let first = &skeleton.errors().unwrap()[0];
        assert_eq!(first.rule, "RequiredColumn");
        assert_eq!((first.row_index, first.line_number), (None, None));
    }

    #[test]
    fn a_skipped_document_is_named_and_counts_nothing() {
        let batch: SdrfValidationBatch =
            fixture(include_str!("../tests/fixtures/sdrf_validate_bulk.json"));
        assert_eq!(
            (batch.file_count, batch.read_count, batch.failed_count),
            (3, 2, 1)
        );
        assert_eq!(batch.record_count, batch.message_count);
        let missing = &batch.files[2];
        assert_eq!(missing.is_valid, None); // not read, which is not "invalid"
        assert_eq!(missing.error.as_ref().unwrap().kind, "usage");
        assert_eq!(batch.failed_files().len(), 1);
        assert_eq!(batch.columns.names()[..2], ["source_index", "source_path"]);
        assert_eq!(batch.messages().unwrap().len() as u64, batch.message_count);
    }

    #[test]
    fn drift_groups_by_finding_and_keeps_the_label_lists() {
        let drift: SdrfDrift = fixture(include_str!("../tests/fixtures/sdrf_lint_cohort.json"));
        assert_eq!(drift.labels, ["cohort", "partner"]);
        let findings = drift.findings().unwrap();
        assert_eq!(findings.len() as u64, drift.finding_count);
        assert!(findings.iter().all(|f| f[0].variant_rank == 0));
        let names = &findings[2];
        assert_eq!(names[0].kind, "ColumnNameVariant");
        assert_eq!(names[0].column_name, None); // a finding about a name, not a column
        assert_eq!(findings[0][0].documents, ["partner"]);
    }

    #[test]
    fn an_assessment_reads_its_verdict_and_integer_fill_rates_as_floats() {
        let a: SdrfAssessment = fixture(include_str!("../tests/fixtures/sdrf_assess_cohort.json"));
        assert_eq!(a.verdict, "Informative");
        assert!(VERDICTS.contains(&a.verdict.as_str()));
        // fill_rate crosses as 1, not 1.0: a float accessor must still read it.
        let rates = a.columns.floats("fill_rate").unwrap();
        assert!(rates
            .iter()
            .all(|r| r.is_some_and(|r| (0.0..=1.0).contains(&r))));
    }

    #[test]
    fn paths_with_refuses_a_verdict_that_does_not_exist() {
        let batch: SdrfAssessmentBatch =
            fixture(include_str!("../tests/fixtures/sdrf_assess_bulk.json"));
        assert_eq!(
            batch.paths_with(&["Skeleton"]).unwrap(),
            ["sdrf_skeleton.sdrf.tsv"]
        );
        let error = batch.paths_with(&["Informativ"]).unwrap_err();
        assert!(matches!(error, MzLibError::Usage(_)));
        assert_eq!(
            batch.verdict_counts,
            VerdictCounts {
                informative: 1,
                partial: 1,
                skeleton: 1
            }
        );
    }

    #[test]
    fn samples_withhold_a_conflict_and_parse_ages() {
        let s: SdrfSamples = fixture(include_str!("../tests/fixtures/sdrf_samples_cohort.json"));
        assert_eq!(s.conflict_count as usize, s.conflicts().unwrap().len());
        let ages = s.ages().unwrap();
        assert_eq!(ages.len(), 6);
        assert_eq!(ages[0].years, Some(58.0));
        // >=90Y has no upper bound: None, told apart from a refusal by its precision.
        assert_eq!(
            (ages[2].max_years, ages[2].precision.as_deref()),
            (None, Some("LowerBound"))
        );
        assert_eq!(ages[4].refusal.as_deref(), Some("reserved_word"));
        let batch: SdrfSamplesBatch =
            fixture(include_str!("../tests/fixtures/sdrf_samples_bulk.json"));
        assert_eq!(batch.sample_count, 8);
        assert_eq!(batch.columns.rows() as u64, batch.record_count);
    }

    #[test]
    fn parsed_ages_stay_aligned_with_the_cells() {
        let ages: ParsedAges = fixture(include_str!("../tests/fixtures/sdrf_parse_age.json"));
        let rows = ages.ages().unwrap();
        assert_eq!(rows.len() as u64, ages.cell_count);
        assert_eq!(rows[9].cell, "");
        assert_eq!(rows[9].refusal.as_deref(), Some("empty"));
        for row in &rows {
            assert_eq!(
                row.refusal.is_none(),
                row.precision.is_some(),
                "{}",
                row.cell
            );
            if let Some(refusal) = &row.refusal {
                assert!(AGE_REFUSALS.contains(&refusal.as_str()));
            }
        }
    }

    #[test]
    fn a_bulk_call_sends_its_options_and_one_path_per_line() {
        let (args, stdin) = bulk_request(
            "validate",
            &["a.sdrf.tsv", "b.sdrf.tsv"],
            &BulkOptions {
                threads: -1,
                on_error: OnError::Skip,
                timeout: None,
            },
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "sdrf",
                "validate",
                "--paths-stdin",
                "--threads",
                "-1",
                "--on-error",
                "skip"
            ]
        );
        assert_eq!(stdin, "a.sdrf.tsv\nb.sdrf.tsv");
        assert_eq!(BulkOptions::default().threads, 1);
        assert_eq!(BulkOptions::default().on_error, OnError::Fail);
    }

    #[test]
    fn bulk_calls_refuse_bad_input_before_anything_is_spawned() {
        let options = BulkOptions::default();
        assert!(bulk_request::<&str>("assess", &[], &options).is_err());
        assert!(bulk_request("assess", &["a", "  "], &options).is_err());
        assert!(bulk_request("assess", &["a\nb"], &options).is_err());
        let zero = BulkOptions {
            threads: 0,
            ..BulkOptions::default()
        };
        assert!(matches!(
            bulk_request("assess", &["a"], &zero),
            Err(MzLibError::Usage(_))
        ));
        assert!(one_path_args("samples", Path::new(" ")).is_err());
        assert_eq!(
            one_path_args("samples", Path::new("x.sdrf.tsv")).unwrap(),
            ["sdrf", "samples", "--path", "x.sdrf.tsv"]
        );
    }

    #[test]
    fn age_cells_keep_blanks_and_end_with_a_newline() {
        assert_eq!(age_lines(&["58Y", "", "63"]).unwrap(), "58Y\n\n63\n");
        // A final blank cell must survive: without the trailing newline it would vanish.
        assert_eq!(age_lines(&["58Y", ""]).unwrap(), "58Y\n\n");
        assert!(age_lines::<&str>(&[]).is_err());
        assert!(age_lines(&["58Y\n63"]).is_err());
    }

    #[test]
    fn lint_shares_pools_document_lines() {
        let lines =
            document_lines(&labelled(&[("a.sdrf.tsv", "one"), ("b.sdrf.tsv", "two")])).unwrap();
        assert_eq!(lines, ["a.sdrf.tsv\tone", "b.sdrf.tsv\ttwo"]);
        assert!(document_lines(&labelled(&[("a.sdrf.tsv", " ")])).is_err());
    }

    // ---- design -------------------------------------------------------------------------------

    fn two_factors() -> DesignOptions {
        DesignOptions {
            condition_columns: vec![
                "factor value[genotype]".to_owned(),
                "factor value[treatment]".to_owned(),
            ],
            ..Default::default()
        }
    }

    #[test]
    fn a_valid_design_reads_into_zero_based_runs() {
        let d: SdrfDesign = fixture(include_str!("../tests/fixtures/sdrf_design_PXD067622.json"));
        assert!(d.is_valid);
        assert_eq!(d.file_count, 24);
        assert_eq!(d.written, None);
        let files = d.files().unwrap();
        assert_eq!(files.len(), 24);
        assert_eq!(
            files[0],
            DesignedFile {
                full_path: "20240830_HF_LC3_MAA_RK_12032_CA_DMSO4.raw".to_owned(),
                file_name: "20240830_HF_LC3_MAA_RK_12032_CA_DMSO4".to_owned(),
                condition: "SPRTN-TurboID CA_DMSO (vehicle)".to_owned(),
                biological_replicate: 0,
                technical_replicate: 0,
                fraction: 0,
            }
        );
        let conditions: std::collections::BTreeSet<&str> =
            files.iter().map(|f| f.condition.as_str()).collect();
        assert_eq!(conditions.len(), 8);
    }

    #[test]
    fn a_design_feeds_flashlfq_unchanged() {
        let d: SdrfDesign = fixture(include_str!("../tests/fixtures/sdrf_design_PXD067622.json"));
        let files = d.files().unwrap();
        let spectra = d.spectra().unwrap();
        let runs = d.run_design().unwrap();
        for ((file, spectrum), run) in files.iter().zip(&spectra).zip(&runs) {
            assert_eq!(spectrum.path, file.full_path);
            assert_eq!(
                spectrum.biological_replicate,
                Some(file.biological_replicate)
            );
            assert_eq!(run.file_name, file.file_name);
            assert_eq!(run.condition.as_deref(), Some(file.condition.as_str()));
            assert_eq!(run.fraction, Some(file.fraction));
        }
    }

    #[test]
    fn a_refusal_is_an_answer_and_refuses_to_feed_flashlfq() {
        let refused: SdrfDesign =
            fixture(include_str!("../tests/fixtures/sdrf_design_PXD049018.json"));
        assert!(!refused.is_valid);
        assert_eq!((refused.refusals.len(), refused.file_count), (20, 0));
        assert!(refused.files().unwrap().is_empty());
        let error = refused.spectra().unwrap_err();
        assert!(
            matches!(error, MzLibError::Usage(ref m) if m.contains("MSB67868ABand_01.raw")),
            "{error}"
        );
        assert!(refused.run_design().is_err());
    }

    #[test]
    fn a_renumbering_is_recorded_and_lands_on_the_hand_numbered_design() {
        let ranked: SdrfDesign =
            fixture(include_str!("../tests/fixtures/sdrf_design_studywide.json"));
        let hand: SdrfDesign =
            fixture(include_str!("../tests/fixtures/sdrf_design_PXD067622.json"));
        assert_eq!(ranked.notes.len(), 7);
        assert!(ranked.notes[1].contains("22 -> 1, 23 -> 2, 24 -> 3"));
        assert_eq!(ranked.columns, hand.columns);
    }

    #[test]
    fn design_sends_tab_joined_columns_searched_files_on_stdin_and_out() {
        let (args, stdin) =
            design_request(Path::new("a.sdrf.tsv"), &DesignOptions::default()).unwrap();
        assert_eq!(args, ["sdrf", "design", "--path", "a.sdrf.tsv"]);
        assert_eq!(stdin, None);

        let options = DesignOptions {
            searched_files: Some(vec![PathBuf::from("/data/r1.raw"), PathBuf::from("r2.raw")]),
            out: Some(PathBuf::from("ExperimentalDesign.tsv")),
            ..two_factors()
        };
        let (args, stdin) = design_request(Path::new("a.sdrf.tsv"), &options).unwrap();
        assert_eq!(
            args,
            [
                "sdrf",
                "design",
                "--path",
                "a.sdrf.tsv",
                "--condition-columns",
                "factor value[genotype]\tfactor value[treatment]",
                "--searched-files-stdin",
                "--out",
                "ExperimentalDesign.tsv"
            ]
        );
        assert_eq!(stdin.as_deref(), Some("/data/r1.raw\nr2.raw\n"));

        // An empty searched list is sent: mzLib's refusal, not a usage error here.
        let empty = DesignOptions {
            searched_files: Some(Vec::new()),
            ..Default::default()
        };
        let (args, stdin) = design_request(Path::new("a.sdrf.tsv"), &empty).unwrap();
        assert!(args.contains(&"--searched-files-stdin".to_owned()));
        assert_eq!(stdin.as_deref(), Some("\n"));
    }

    #[test]
    fn design_refuses_bad_input_before_anything_is_spawned() {
        let blank = DesignOptions {
            condition_columns: vec![" ".to_owned()],
            ..Default::default()
        };
        assert!(design_request(Path::new("a.sdrf.tsv"), &blank).is_err());
        let tabbed = DesignOptions {
            condition_columns: vec!["a\tb".to_owned()],
            ..Default::default()
        };
        assert!(design_request(Path::new("a.sdrf.tsv"), &tabbed).is_err());
        let bad_file = DesignOptions {
            searched_files: Some(vec![PathBuf::from("a\nb.raw")]),
            ..Default::default()
        };
        assert!(design_request(Path::new("a.sdrf.tsv"), &bad_file).is_err());
        let empty_out = DesignOptions {
            out: Some(PathBuf::from("")),
            ..Default::default()
        };
        assert!(design_request(Path::new("a.sdrf.tsv"), &empty_out).is_err());
        assert!(design_request(Path::new(""), &DesignOptions::default()).is_err());
        assert_eq!(DesignOptions::default().timeout, Some(DEFAULT_TIMEOUT));
    }

    // ---- isobaric design (mzLib #1380) ----------------------------------------------------------

    const POOL_PATTERN: &str = r"TMT_?pool(?<plex>\d+)";
    const BATCH: &str = "comment[sample preparation batch]";

    /// A recording in tests/fixtures/, read at run time: the isobaric and restrict ones are large.
    fn recorded_file<T: serde::de::DeserializeOwned>(name: &str) -> T {
        let text = std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("fixtures")
                .join(name),
        )
        .unwrap();
        serde_json::from_str(&text).unwrap()
    }

    fn isobaric(name: &str) -> SdrfIsobaricDesign {
        recorded_file(name)
    }

    fn tmt10(pattern: Option<&str>) -> IsobaricDesignOptions {
        IsobaricDesignOptions {
            kit: Some("TMT10".to_owned()),
            plex_file_name_pattern: pattern.map(str::to_owned),
            ..Default::default()
        }
    }

    #[test]
    fn a_valid_isobaric_design_projects_every_annotated_channel() {
        let (args, stdin) = isobaric_design_request(
            Path::new("PXD008841.subset.sdrf.tsv"),
            &tmt10(Some(POOL_PATTERN)),
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "sdrf",
                "isobaric-design",
                "--path",
                "PXD008841.subset.sdrf.tsv",
                "--kit",
                "TMT10",
                "--plex-file-name-pattern",
                POOL_PATTERN
            ]
        );
        assert_eq!(stdin, None);

        let d = isobaric("sdrf_isobaric_design_PXD008841.json");
        assert!(d.is_valid && d.refusals.is_empty() && d.notes.is_empty());
        assert_eq!(d.kit.as_deref(), Some("TMT10"));
        assert_eq!(d.kit_type.as_deref(), Some("TMT10"));
        assert_eq!(d.kit_channels.as_ref().map(Vec::len), Some(10));
        assert_eq!(d.plex_file_name_pattern.as_deref(), Some(POOL_PATTERN));
        assert_eq!(
            (d.plex_column.as_deref(), d.single_plex.as_deref()),
            (None, None)
        );
        assert_eq!(
            d.plex_ids,
            BTreeMap::from([("1".to_owned(), 1), ("5".to_owned(), 2)])
        );
        assert_eq!((d.plex_count, d.file_count, d.record_count), (2, 4, 40));
        assert_eq!(d.channels().unwrap().len(), 40);
        assert_eq!(d.file_key_column.as_deref(), Some("comment[data file]"));
        assert_eq!(d.condition_columns, ["factor value[disease]"]);
        assert!(!d.condition_columns_declared && !d.searched_files_given);
        assert_eq!(d.written, None);
        assert_eq!(
            d.columns.names(),
            [
                "full_path",
                "plex",
                "plex_id",
                "sample_name",
                "channel",
                "condition",
                "biological_replicate",
                "fraction",
                "technical_replicate"
            ]
        );
        assert_eq!(d.caveats.len(), 8);
        assert!(d.caveats[0].starts_with("1-BASED"));
        assert!(d.report.contains("4 file(s), 2 plex(es), kit TMT10"));
    }

    #[test]
    fn the_isobaric_numbers_stay_one_based() {
        let d = isobaric("sdrf_isobaric_design_PXD008841.json");
        let channels = d.channels().unwrap();
        let set = |f: fn(&DesignedChannel) -> u32| -> std::collections::BTreeSet<u32> {
            channels.iter().map(f).collect()
        };
        assert_eq!(set(|c| c.fraction), [1, 2].into());
        assert_eq!(set(|c| c.technical_replicate), [1].into());
        assert_eq!(set(|c| c.biological_replicate), [1].into());
        assert_eq!(set(|c| c.plex_id), [1, 2].into());
        let first_file: Vec<&DesignedChannel> = channels
            .iter()
            .filter(|c| c.full_path == channels[0].full_path)
            .collect();
        let order: Vec<&str> = first_file.iter().map(|c| c.channel.as_str()).collect();
        assert_eq!(
            Some(order),
            d.kit_channels
                .as_ref()
                .map(|k| k.iter().map(String::as_str).collect()),
            "reporter m/z order"
        );
        assert_eq!(first_file.last().unwrap().sample_name, "pool");
    }

    #[test]
    fn an_isobaric_refusal_is_a_result() {
        let (args, _) =
            isobaric_design_request(Path::new("PXD008841.subset.sdrf.tsv"), &tmt10(None)).unwrap();
        assert_eq!(
            args,
            [
                "sdrf",
                "isobaric-design",
                "--path",
                "PXD008841.subset.sdrf.tsv",
                "--kit",
                "TMT10"
            ]
        );
        let d = isobaric("sdrf_isobaric_design_PXD008841_noplex.json");
        assert!(!d.is_valid);
        assert!(d.refusals[0].starts_with("No plex source was declared"));
        assert_eq!(d.plex_source, "(not declared)");
        assert_eq!((d.file_count, d.record_count), (0, 0));
        assert!(d.channels().unwrap().is_empty() && d.plex_ids.is_empty());
    }

    #[test]
    fn a_declared_batch_column_gives_three_plexes() {
        let options = IsobaricDesignOptions {
            kit: Some("TMT18".to_owned()),
            plex_column: Some(BATCH.to_owned()),
            ..Default::default()
        };
        let (args, _) =
            isobaric_design_request(Path::new("PXD061609.subset.sdrf.tsv"), &options).unwrap();
        assert_eq!(args[4..], ["--kit", "TMT18", "--plex-column", BATCH]);
        let d = isobaric("sdrf_isobaric_design_PXD061609.json");
        assert!(d.is_valid);
        assert_eq!(
            d.plex_ids,
            BTreeMap::from([
                ("HighSalt".to_owned(), 1),
                ("LowSalt".to_owned(), 2),
                ("NoSAX".to_owned(), 3)
            ])
        );
        assert_eq!(d.record_count, 108);
        assert_eq!(d.kit_channels.map(|k| k.len()), Some(18));
    }

    #[test]
    fn one_plex_over_three_arms_is_refused_every_reason_at_once() {
        let d = isobaric("sdrf_isobaric_design_PXD061609_oneplex.json");
        assert!(!d.is_valid);
        assert_eq!(d.single_plex.as_deref(), Some("tube1"));
        assert_eq!(d.refusals.len(), 36);
        assert!(d.refusals.iter().all(|r| r.contains("more than once")));
    }

    #[test]
    fn isobaric_condition_columns_and_searched_files_cross_as_design_sends_them() {
        let options = IsobaricDesignOptions {
            condition_columns: vec!["factor value[disease]".to_owned()],
            searched_files: Some(vec![
                PathBuf::from("/data/a-calib.mzML"),
                PathBuf::from("b-calib.mzML"),
            ]),
            out: Some(PathBuf::from("design")),
            ..tmt10(Some(POOL_PATTERN))
        };
        let (args, stdin) =
            isobaric_design_request(Path::new("PXD008841_searched.sdrf.tsv"), &options).unwrap();
        assert_eq!(
            args[args.len() - 5..],
            [
                "--condition-columns",
                "factor value[disease]",
                "--searched-files-stdin",
                "--out",
                "design"
            ]
        );
        assert_eq!(stdin.as_deref(), Some("/data/a-calib.mzML\nb-calib.mzML\n"));
    }

    #[test]
    fn a_restricted_document_keys_on_the_searched_name() {
        let d = isobaric("sdrf_isobaric_design_searched.json");
        assert_eq!(
            d.file_key_column.as_deref(),
            Some(SEARCHED_DATA_FILE_COLUMN)
        );
        assert!(d.searched_files_given);
        let paths: std::collections::BTreeSet<String> = d
            .channels()
            .unwrap()
            .into_iter()
            .map(|c| c.full_path)
            .collect();
        // Recorded on Windows: compared as text, never split on a separator.
        let expected: std::collections::BTreeSet<String> = [
            r"D:\search\HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01-calib.mzML",
            r"D:\search\HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02-calib.mzML",
            r"D:\search\HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr01-calib.mzML",
            r"D:\search\HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr02-calib.mzML",
        ]
        .map(str::to_owned)
        .into();
        assert_eq!(paths, expected);
        let written = d.written.unwrap();
        assert_eq!((written.file_count, written.record_count), (4, 40));
        assert!(written.path.ends_with("TmtDesign.txt"));
    }

    #[test]
    fn no_isobaric_option_is_sent_unless_given() {
        let (args, stdin) = isobaric_design_request(
            Path::new("PXD008841.subset.sdrf.tsv"),
            &IsobaricDesignOptions::default(),
        )
        .unwrap();
        assert_eq!(
            args,
            [
                "sdrf",
                "isobaric-design",
                "--path",
                "PXD008841.subset.sdrf.tsv"
            ]
        );
        assert_eq!(stdin, None);
        assert_eq!(
            IsobaricDesignOptions::default().timeout,
            Some(DEFAULT_TIMEOUT)
        );
    }

    #[test]
    fn bad_isobaric_arguments_are_refused_before_the_bridge() {
        let cases: Vec<(IsobaricDesignOptions, &str)> = vec![
            (
                IsobaricDesignOptions {
                    kit: Some(String::new()),
                    ..Default::default()
                },
                "kit must be non-blank",
            ),
            (
                IsobaricDesignOptions {
                    plex_column: Some("  ".to_owned()),
                    ..Default::default()
                },
                "plex_column must be non-blank",
            ),
            (
                IsobaricDesignOptions {
                    plex_file_name_pattern: Some("a\nb".to_owned()),
                    ..Default::default()
                },
                "line break",
            ),
            (
                IsobaricDesignOptions {
                    single_plex: Some(String::new()),
                    ..Default::default()
                },
                "single_plex must be non-blank",
            ),
            (
                IsobaricDesignOptions {
                    condition_columns: vec![" ".to_owned()],
                    ..Default::default()
                },
                "blank",
            ),
            (
                IsobaricDesignOptions {
                    searched_files: Some(vec![PathBuf::from("")]),
                    ..Default::default()
                },
                "blank",
            ),
            (
                IsobaricDesignOptions {
                    out: Some(PathBuf::from("")),
                    ..Default::default()
                },
                "out must name a folder",
            ),
        ];
        for (options, message) in cases {
            let error = isobaric_design_request(Path::new("a.sdrf.tsv"), &options).unwrap_err();
            assert!(
                matches!(&error, MzLibError::Usage(m) if m.contains(message)),
                "{message}: {error:?}"
            );
        }
        assert!(isobaric_design_request(Path::new(""), &tmt10(None)).is_err());
    }

    // ---- restrict (mzLib #1451) -----------------------------------------------------------------

    const POOL1: &str = "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0";
    const POOL5: &str = "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0";

    fn messy() -> Vec<SearchedFile> {
        vec![
            SearchedFile::from(format!(r"D:\search\{POOL1}1-calib.mzML")),
            SearchedFile::from(format!("{POOL1}2-calib.mzML")),
            SearchedFile::from(format!("{POOL1}2.raw")),
            SearchedFile::acquired_from("sample_A.mzML", format!("{POOL5}1.raw")),
            SearchedFile::from("not_in_sdrf.mzML"),
        ]
    }

    fn restricted(name: &str) -> RestrictedSdrf {
        recorded_file(name)
    }

    #[test]
    fn restrict_sends_an_acquired_name_only_where_one_is_given() {
        let options = RestrictOptions {
            limit: Some(2),
            ..Default::default()
        };
        let (args, stdin) =
            restrict_request(Path::new("PXD008841.subset.sdrf.tsv"), &messy(), &options).unwrap();
        assert_eq!(
            args,
            [
                "sdrf",
                "restrict",
                "--path",
                "PXD008841.subset.sdrf.tsv",
                "--limit",
                "2"
            ]
        );
        assert_eq!(
            stdin,
            format!(
                "D:\\search\\{POOL1}1-calib.mzML\n{POOL1}2-calib.mzML\n{POOL1}2.raw\n\
                 sample_A.mzML\t{POOL5}1.raw\nnot_in_sdrf.mzML\n"
            )
        );
        let scoped = restricted("sdrf_restrict_PXD008841.json");
        let sent: Vec<String> = messy()
            .into_iter()
            .map(|f| f.searched.to_string_lossy().into_owned())
            .collect();
        assert_eq!(scoped.searched_files, sent);
        assert_eq!(scoped.acquired_name_count, 1);
    }

    #[test]
    fn the_restrict_reports_are_mzlibs() {
        let scoped = restricted("sdrf_restrict_PXD008841.json");
        let doc = &scoped.document;
        assert_eq!((scoped.input_row_count, doc.row_count), (40, 30));
        assert_eq!(
            (doc.returned_count, doc.truncated, doc.offset),
            (2, true, 0)
        );
        assert_eq!(scoped.searched_without_row, ["not_in_sdrf.mzML"]);
        assert_eq!(scoped.dropped_data_files, [format!("{POOL5}2.raw")]);
        assert_eq!(
            scoped.ambiguous,
            [format!("{POOL1}2: {POOL1}2-calib.mzML, {POOL1}2.raw")]
        );
        assert!(scoped.shared_stem.is_empty());
        assert_eq!(scoped.written, None);
        assert_eq!(doc.caveats.len(), 6);
    }

    #[test]
    fn the_restricted_document_is_an_sdrf_document() {
        let scoped = restricted("sdrf_restrict_PXD008841.json");
        let doc = &scoped.document;
        let data_file = doc.index_of("comment[data file]").unwrap();
        assert_eq!(doc.columns[data_file + 1], SEARCHED_DATA_FILE_COLUMN);
        let raw = format!("{POOL1}1.raw");
        assert_eq!(doc.value("comment[data file]"), [Some(raw.as_str()); 2]);
        let calib = format!("{POOL1}1-calib.mzML");
        assert_eq!(
            scoped.searched_data_files(),
            [Some(calib.as_str()); 2],
            "the name, no folder"
        );
        assert_eq!(doc.ragged_row_count(), 0, "mzLib pads a short row");
    }

    #[test]
    fn restrict_sends_one_name_per_line_and_passes_out_through() {
        let clean: Vec<SearchedFile> = [POOL1, POOL5]
            .iter()
            .flat_map(|p| {
                (1..=2).map(move |i| SearchedFile::from(format!(r"D:\search\{p}{i}-calib.mzML")))
            })
            .collect();
        let options = RestrictOptions {
            limit: Some(2),
            out: Some(PathBuf::from("PXD008841_searched.sdrf.tsv")),
            ..Default::default()
        };
        let (args, stdin) =
            restrict_request(Path::new("PXD008841.subset.sdrf.tsv"), &clean, &options).unwrap();
        assert_eq!(
            args[args.len() - 4..],
            ["--limit", "2", "--out", "PXD008841_searched.sdrf.tsv"]
        );
        assert_eq!(stdin.lines().count(), 4);
        assert!(
            !stdin.contains('\t'),
            "no acquired name is sent unless one is given"
        );

        let scoped = restricted("sdrf_restrict_PXD008841_out.json");
        assert_eq!(
            scoped.written,
            Some(WrittenSdrf {
                path: "PXD008841_searched.sdrf.tsv".to_owned(),
                row_count: 40
            })
        );
        assert_eq!(
            (scoped.document.row_count, scoped.document.returned_count),
            (40, 2)
        );
    }

    #[test]
    fn bad_searched_files_are_refused_before_the_bridge() {
        let path = Path::new("PXD008841.subset.sdrf.tsv");
        let options = RestrictOptions::default();
        let cases: Vec<(Vec<SearchedFile>, &str)> = vec![
            (Vec::new(), "at least one searched file"),
            (vec![SearchedFile::from("")], "blank"),
            (vec![SearchedFile::from("a\tb")], "tab or newline"),
            (
                vec![SearchedFile::acquired_from("a.mzML", "")],
                "acquired file for 'a.mzML'",
            ),
            (
                vec![SearchedFile::acquired_from("a.mzML", "b\n.raw")],
                "tab or newline",
            ),
        ];
        for (files, message) in cases {
            let error = restrict_request(path, &files, &options).unwrap_err();
            assert!(
                matches!(&error, MzLibError::Usage(m) if m.contains(message)),
                "{message}: {error:?}"
            );
        }
        let blank_out = RestrictOptions {
            out: Some(PathBuf::from(" ")),
            ..Default::default()
        };
        let error = restrict_request(path, &messy(), &blank_out).unwrap_err();
        assert!(error.to_string().contains("ending in .tsv"));
        assert!(restrict_request(Path::new(""), &messy(), &options).is_err());
        assert_eq!(RestrictOptions::default().timeout, Some(DEFAULT_TIMEOUT));
    }

    // ---- pool --out names the file mzLib wrote ----------------------------------------------------

    #[test]
    fn written_sdrf_reads_the_path_mzlib_reports() {
        // pyMzLib 0.6.0's bridge reports the file mzLib wrote, which is x.tsv.sdrf.tsv for an out of
        // x.tsv; this crate passes it through, whatever it is.
        let written: WrittenSdrf =
            serde_json::from_str(r#"{"path": "merged.tsv.sdrf.tsv", "row_count": 24}"#).unwrap();
        assert_eq!(written.path, "merged.tsv.sdrf.tsv");
        assert_eq!(written.row_count, 24);
    }
}
