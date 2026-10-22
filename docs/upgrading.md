# Upgrading: what changed in your results

mzLibRust drives one mzLib bridge, and every bridge it pins moves mzLib forward. Most of what a new
mzLib brings is new capability. This page is the rest: **changes that give you a different answer
from the same file and the same call**. Read the section for the bridge you are moving to before
you compare numbers across the upgrade.

Each entry names the mzLib pull request it comes from, so you can read the change itself.
[`CHANGELOG.md`](../CHANGELOG.md) has everything else. Which mzLib a bridge was built from is in
`mzlib::bridge_version()?.mzlib`. pyMzLib keeps the same page, entry for entry, because both
bindings drive the same bridge.

## The pyMzLib 0.6.0 bridge (still built from mzLib 1.0.595)

mzLib did not move. These are the bridge's own.

### Same file, different value

- **`stats::fit` and `fit_with` on fewer than 3 features report no shared prior variance.** With
  so few features mzLib fits no prior ([mzLib#1398](https://github.com/smith-chem-wisc/mzLib/pull/1398),
  as limma's `squeezeVar` does): `VariancePrior::df` is 0 and each feature keeps its own variance.
  `VariancePrior::scale` used to report the first feature's variance as though it were the prior;
  it is now `None`, and each feature's value is in `FitRow::prior_variance`. Tables of 3 or more
  features are unchanged.
- **`stats::fit_with`'s caveats name the variance-prior estimator.** The first caveat now begins
  `moments_legacy:` (or `marginal_likelihood:` when you choose it with `FitOptions::estimator`), and
  the `residual_df_differ` caveat names `marginal_likelihood` as mzLib's estimator for that case.
  The `bh_adjusted` caveat says "per coefficient and per contrast". Every number is unchanged;
  only code that matches caveat text needs a look.

## The pyMzLib 0.5.0 bridge (built from mzLib 1.0.595)

This bridge moves mzLib from 1.0.593 to 1.0.595, so it carries 1.0.594's changes too. Three of the
changes below are this crate's own: it stopped choosing defaults that are mzLib's to set, and it
stopped adding a PRIDE time limit on top of mzLib's.

### Files that used to be refused now read

- **DIA-NN 2.x `report.parquet`** ([mzLib#1384](https://github.com/smith-chem-wisc/mzLib/pull/1384)).
  A `.parquet` was `MzLibError::Usage` (`File type not supported`). One whose columns include
  `Precursor.Id`, `Stripped.Sequence` and `Run` now reads as `DiaNnReportParquet`, a `quantifiable`
  type like the DIA-NN `report.tsv`, with retention time in minutes. Any other `.parquet` is still
  `MzLibError::Usage`, now `Parquet file type not supported`. `readers::formats()` lists one more
  entry, and one more of them offers `quantifiable`: code that asserts either count needs updating.

### Calls that used to succeed are now refused

- **ProForma input that lost part of itself** ([mzLib#1462](https://github.com/smith-chem-wisc/mzLib/pull/1462)).
  `peptidoform::convert_with` with `source: "ProForma"` converted `PEPTIDE-[UNIMOD:2][UNIMOD:35]`
  with one C-terminal modification dropped, and read `PEPTIDE-[UNIMOD:2]K` as `PEPTIDEK`. Their
  rows are now `status` `"failed"`: stacked terminal modifications with `failure_reason`
  `"IncompatibleModifications"` and `incompatible_items` `["Multiple C-terminal modifications (2
  found)."]`, and a residue after a C-terminal modification, a leading `-` or a modification with
  no sequence with `failure_reason` `"InvalidSequence"`.
- **SUMO remnants named `... on D`** ([mzLib#1432](https://github.com/smith-chem-wisc/mzLib/pull/1432)).
  mzLib's human SUMO-1 and SUMO-2/3 remnants (`DVIEVYQEQTGG (SUMO-1 Site human)`,
  `DVFQQQTGG (SUMO-2/3 Site human)`) now target K, where they targeted D. In `peptidoform::convert`
  and `convert_with`, a full sequence that names the `on D` form no longer converts (its row is
  `status` `"failed"`); the `on K` form does.
- **TMTpro on S, T, Y or H** ([mzLib#1419](https://github.com/smith-chem-wisc/mzLib/pull/1419)). mzLib
  used to turn these into TMT, 224 Da, so `PEPS[Common Fixed:TMTpro on S]K` converted at the wrong
  mass. They now fail to convert: the row is `status` `"failed"`, with the modification in
  `incompatible_items`. TMTpro on K is unchanged (+304.2071 Da).

### Calls that used to fail now succeed

- **`flashlfq::quantify` with gapped replicate numbers** ([mzLib#1422](https://github.com/smith-chem-wisc/mzLib/pull/1422)).
  A design whose biological replicates skip a number, with fractions, returned an error; it now
  quantifies. So does the numbering the bridge gives runs whose `SpectraFile::biological_replicate`
  is `None`, which skips at every condition after the first. `QuantifyOptions::bayesian_protein_quant`
  on such a design now writes `BayesianProteinQuant.tsv` to `QuantifyOptions::output_directory`.
- **A condition whose name ends in digits** ([mzLib#1424](https://github.com/smith-chem-wisc/mzLib/pull/1424)).
  Condition `A1` biorep 0 beside condition `A` biorep 10 failed with `MzLibError::Bridge`
  (`IndexOutOfRangeException`) in `flashlfq::quantify` and `flashlfq::median_polish`, because the
  two spell the same text. It now quantifies.

### Same file, different value

- **Entrapment and contaminant flags read every joined label** ([mzLib#1371](https://github.com/smith-chem-wisc/mzLib/pull/1371)).
  In `readers::read_records` and `readers::read_protein_groups`, the `is_decoy`, `is_contaminant`
  and `is_entrapment` columns are true when any `|`-joined parent is D, C or E. A protein group
  labelled `T|C` is now a contaminant and `T|ET` entrapment; on a `.psmtsv` the "E" in
  MetaMorpheus's "Output too long for Excel" no longer reads as entrapment. Labels MetaMorpheus
  writes itself read as before.
- **`.psmtsv` and `.osmtsv` records gain an `entrapment_fraction` column** ([mzLib#1371](https://github.com/smith-chem-wisc/mzLib/pull/1371)),
  read with `columns.floats("entrapment_fraction")`: the share of each PSM that is an entrapment
  discovery. A shared `T|ET` counts half when the sequences differ and nothing when both proteins
  carry the same peptide. Use it, not `is_entrapment`, to count entrapment discoveries. It is
  `None` where mzLib cannot line the label up with the full sequence.
- **Lactylation on K weighs one hydrogen more** ([mzLib#1417](https://github.com/smith-chem-wisc/mzLib/pull/1417)).
  mzLib's formula is now C3H4O2, 72.021129 Da, where it was C3H3O2, 71.013304 Da.
  `peptidoform::convert_with` with `target: "MassShift"` turns `PEPK[Less Common:Lactylation on K]R`
  into `PEPK[+72.021129]R`, and the `"ProForma"` and `"Unimod"` targets write `UNIMOD:2114`. The
  `pro_forma` column of `read_records` on a `.psmtsv` changes the same way.
- **`flashlfq::quantify_with` with `normalize: true`** ([mzLib#1422](https://github.com/smith-chem-wisc/mzLib/pull/1422)).
  The reference is now the first condition's lowest biological replicate, and each fraction's
  lowest technical replicate, taken in numeric order. A design numbered from 1 is now normalized;
  before, biorep normalization did nothing for it. A gapless design numbered from 0 and listed in
  order is unchanged.
- **`stats::fit` and `fit_with` on fewer than 3 features** ([mzLib#1398](https://github.com/smith-chem-wisc/mzLib/pull/1398)).
  mzLib fits no prior, as limma does: `VariancePrior::df` is 0, each feature keeps its own
  variance, and `t` and `p_value` are ordinary per-feature t-tests. They used to borrow a prior
  fitted to the two features. **On this bridge `VariancePrior::scale` still reports the first
  feature's own variance there**, as though it were the prior; the pyMzLib 0.6.0 bridge reports
  `None` (above). Read each feature's value from `FitRow::prior_variance`.
- **`peptidoform::fragments_with` with a `|P` protease** ([mzLib#1320](https://github.com/smith-chem-wisc/mzLib/pull/1320)).
  `protease: Some("trypsin|P".into())` cut a protein that begins `RP` or `KP` after its first
  residue; it no longer does, and the other `|P` proteases change the same way. mzLib's default
  `"trypsin"`, and any protein that begins with M, are unaffected.
- **PRIDE: a failure that can pass is retried** ([mzLib#1412](https://github.com/smith-chem-wisc/mzLib/pull/1412)).
  Every `pride` function retries such a failure three times, after 5, 20 and 60 seconds, so an
  outage takes longer to report. An FTP listing that comes back empty, which `pride::list_ftp_files`
  used to return as no files (and so `MzLibError::ProjectNotFound`), is retried and then returned
  as `MzLibError::Bridge`. A failed download leaves `<name>.partial` and `<name>.partial.validator`
  behind, so the next call resumes.
- **ProForma names UniProt modifications by their UNIMOD accession**
  ([mzLib#1403](https://github.com/smith-chem-wisc/mzLib/pull/1403)). `peptidoform::convert_with`
  with `target: "ProForma"` turns `YPIEH[UniProt:Tele-methylhistidine on H]GIVTNWDDMEK` into
  `YPIEH[UNIMOD:34]GIVTNWDDMEK`, where it wrote the UniProt name back. The `pro_forma` column of
  `read_records` on a `.psmtsv` written without a ProForma column changes the same way. A string
  comparison of ProForma across the upgrade will see these peptides as different.
  `SequenceConversions::caveats` holds five caveats instead of six, because the one describing
  this gap is gone.
- **`peptidoform::fragments`, and `fragments_with` without `protease` or `dissociation`**
  (pyMzLib #83). Every digestion option in `FragmentOptions` is now an `Option`, `None` by default,
  and an option left `None` is not sent, so it takes mzLib's default. Two of those differ from what
  this crate chose before:
  - the protease is mzLib's `"trypsin"`, which also cuts before proline, where the crate sent
    `"trypsin|P"`, which does not;
  - the fragments are HCD's b and y ions, where the crate sent ETD for c and z• ions. mzLib names no
    default dissociation type, so the bridge uses HCD until it does
    ([mzLib#1488](https://github.com/smith-chem-wisc/mzLib/issues/1488)).

  Serum albumin's default digest goes from 303 to 310 peptidoforms. For the old answer, pass
  `FragmentOptions { protease: Some("trypsin|P".into()), dissociation: Some("ETD".into()), ..Default::default() }`.
  The other defaults have the same values as before, now set by mzLib.
- **`stats::fit_with` with `trend: true` on fewer than 30 features**
  ([mzLib#1398](https://github.com/smith-chem-wisc/mzLib/pull/1398), pyMzLib #82). Without
  `spline_basis`, the trend spline now has limma's number of basis functions,
  1 + (n ≥ 3) + (n ≥ 6) + (n ≥ 30) for n features, where the bridge always asked for 4. A small
  trend fit's prior and moderated statistics change to limma's. Thirty or more features give the
  same answer, and `spline_basis: Some(4)` gives the old answer at any size.

### Same failure, different error

- **A PRIDE outage that hangs is `MzLibError::ServiceUnavailable`, not `MzLibError::Timeout`.**
  This crate's own change. `pride::list_files`, `list_ftp_files` and `search`, and
  `ListOptions::default()` and `SearchOptions::default()`, no longer set a 300-second timeout of
  their own. mzLib bounds each request itself ([mzLib#1412](https://github.com/smith-chem-wisc/mzLib/pull/1412)),
  so a hang is reported as an outage after about 8 minutes per request. To keep the old limit, set
  `ListOptions::timeout` or `SearchOptions::timeout` to `Some(Duration::from_secs(300))`, or call
  `list_ftp_files_with(accession, Some(Duration::from_secs(300)))`.

## The pyMzLib 0.3.0 and 0.4.0 bridges (built from mzLib 1.0.593)

### Files that used to be refused now read

- **MetaMorpheus protein-group tables written without quantification**
  ([mzLib#1365](https://github.com/smith-chem-wisc/mzLib/pull/1365)). `AllProteinGroups.tsv` and
  each file's `<file>_ProteinGroups.tsv` failed with `Tsv file type not supported`. They now read
  through `read_protein_groups` as `MetaMorpheusQuantifiedProteinGroups`, with spectral counts and
  **no intensities**: `intensity` is in `absent_fields`, meaning no basis, not zero.
- **RNA quantification tables** ([mzLib#1388](https://github.com/smith-chem-wisc/mzLib/pull/1388)).
  `AllQuantifiedTranscriptGroups.tsv` and `AllQuantifiedOligos.tsv` are new file types, so
  `readers::formats()` lists more entries than before. Code that asserts a count of formats needs
  updating.

### Same file, different value

- **Four modifications write their Unimod accession in `pro_forma`**
  ([mzLib#1328](https://github.com/smith-chem-wisc/mzLib/pull/1328)). In `read_records` on a
  MetaMorpheus `.psmtsv` or `.osmtsv`, `GG (Ubiquitination Site)` is now `[UNIMOD:121]`, both
  `Myristoylation` entries `[UNIMOD:45]` and `EQIGG` `[UNIMOD:846]`, where they were written by name.
  **Masses are unchanged.** A string comparison of `pro_forma` across the upgrade will see these
  peptides as different; a pipeline keyed on Unimod will now see ubiquitination sites it silently
  dropped before.

### Same failure, different error

- **A read fault on an existing `.mzid`** ([mzLib#1362](https://github.com/smith-chem-wisc/mzLib/pull/1362))
  is now `MzLibError::Bridge` with `error_type` `MzLibException`, naming the file, where it was an
  `IOException`. A missing file is still `MzLibError::Usage`. A plain `.mzid` is now streamed from
  disk, so reading a large one needs far less memory.
- **PRIDE errors name the file and the host, never the URL**
  ([mzLib#1350](https://github.com/smith-chem-wisc/mzLib/pull/1350)), so a reviewer token in a
  query string cannot reach a log. Every transport failure is still
  `MzLibError::ServiceUnavailable`. Code that parsed a URL out of an error message will find none.

## The pyMzLib 0.2.0 bridge (built from mzLib 1.0.592)

### Files that used to be refused now read

- **mzIdentML** (`.mzid`, and `.mzid.gz` without unpacking it,
  [mzLib#1313](https://github.com/smith-chem-wisc/mzLib/pull/1313)), and MetaMorpheus's
  `AllQuantifiedProteinGroups.tsv` and `AllQuantifiedPeptides.tsv`
  ([mzLib#1347](https://github.com/smith-chem-wisc/mzLib/pull/1347)).
- **Pytheas match output** ([mzLib#1277](https://github.com/smith-chem-wisc/mzLib/pull/1277)). A
  `.txt` whose first five lines carry `#theoretical_digest` reads as Pytheas; any other `.txt` is
  still Crux.
- **Current FlashLFQ and MetaMorpheus `QuantifiedPeaks.tsv`**
  ([mzLib#1345](https://github.com/smith-chem-wisc/mzLib/pull/1345)), which have no `MBR Score`
  column, read instead of failing with `HeaderValidationException`.

### Same file, different value

- **`pro_forma` is filled on older MetaMorpheus files**
  ([mzLib#1346](https://github.com/smith-chem-wisc/mzLib/pull/1346)). MetaMorpheus 1.1.11 and
  earlier wrote no ProForma column, so `read_records` gave `None`; mzLib now converts
  `full_sequence`. Ambiguous (`|`-joined) or unconvertible rows are still `None`.
- **`mbr_score` is `None` where it was `0`** in a FlashLFQ `QuantifiedPeaks.tsv`
  ([mzLib#1345](https://github.com/smith-chem-wisc/mzLib/pull/1345)), for every peak without a
  score. A `0` used to say "scored, and scored nothing".
- **MGF `one_based_precursor_scan_number`** is no longer always `None`
  ([mzLib#1227](https://github.com/smith-chem-wisc/mzLib/pull/1227)): mzLib reads it from the
  `PRECURSORSCAN` header it writes. It is still `None` on MGF files mzLib did not write.
- **FlashLFQ with match-between-runs is reproducible across thread counts**
  ([mzLib#1155](https://github.com/smith-chem-wisc/mzLib/pull/1155)). Before, more than one thread
  could give different peptide and protein numbers for identical inputs.

## Earlier bridges

- **MSFragger retention times are minutes** from the pyMzLib 0.1.0 bridge
  ([mzLib#1116](https://github.com/smith-chem-wisc/mzLib/pull/1116)). Before, an MSFragger
  `psm.tsv` passed seconds through, and FlashLFQ quantified it in the wrong window.
