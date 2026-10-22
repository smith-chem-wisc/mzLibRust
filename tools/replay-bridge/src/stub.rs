//! The stand-in bridge executable. Compiled by `build.rs`, not by Cargo, with the replay table from
//! `table.rs` baked in; std only.
//!
//! It answers a call from a fixture recorded from the real bridge, and only when that recording is
//! **consistent with the call**, so an example cannot print output recorded for other arguments.
//! The rules are pyMzLib's `pkg/python/tests/replay_bridge.py`, line for line:
//!
//! * `--path` must name the file the recording was made from (compared by file name), and any
//!   other `--option value` whose snake_case name is a top-level key of the recording must equal it;
//! * `--limit`/`--offset` must reproduce the recording's `returned_count` from its `record_count`
//!   (or `row_count`), and a `--flag` with a `<flag>_included` key must match it;
//! * a `--paths-stdin` call fits only a bulk recording (`files[]` whose entries carry a `path`, and
//!   no top-level `path`, per BULK.md), and a one-path call only a one-document recording;
//! * a recording with a `written` key fits an `--out` call only if `written` is set, and a call
//!   without `--out` only if it is null;
//! * an input-file option the bridge echoes as `<name>_file` (`--responses` -> `responses_file`), or
//!   under another name (`ECHOED_AS`), must name the recording's file;
//! * a `proteins read` recording made with an accession filter answers only a filtered call, and
//!   the reverse; a recording for one isobaric `kit` answers only a call that names a kit;
//! * a `proteins read` recording with a `sequence` column answers only a `--sequences` call, and
//!   the reverse; a `proteins align` recording with a residue map answers only a call without
//!   `--no-residue-map`; a `genes orthology` recording answers `--pair-status` only for the
//!   direction it holds, and `--restrict-gene-trees` only when its gene trees were restricted;
//! * for `quant flashlfq` (`STDIN_ECHO`), the runs sent on stdin must be the recording's runs, for
//!   `peptidoform convert` and `proteins parse-accession` the lines sent must be the recording's
//!   `input` column, in order, for `proteins align` the pairs must be its aligned rows without
//!   their gaps, for `sdrf restrict` the searched files must be its `searched_files`, in order,
//!   with as many acquired names as its `acquired_name_count`, and for `genes orthology` the
//!   species must be its `species`, with as many dumps;
//! * a recording listed in `RECORDED_WITH` fits only a call that sends the options it was made with.
//!
//! No match, or more than one, is answered as a usage error naming the candidates, so the doctest
//! fails and says why.

use std::io::Write;

#[allow(dead_code)]
enum V {
    Null,
    Bool(bool),
    Text(&'static str),
    Compound,
}

struct Recording {
    fixture: &'static str,
    envelope: &'static str,
    fields: &'static [(&'static str, V)],
    bulk: bool,
    /// The file names of the recording's `spectra_files[].full_path`, for `STDIN_ECHO` verbs.
    runs: &'static [&'static str],
    /// What a `STDIN_ECHO` recording was made from, one entry per stdin line: the `input` column
    /// of `peptidoform convert` and `proteins parse-accession`, the `searched_files` of
    /// `sdrf restrict`, the gapless pairs of `proteins align`, and the species, `--` and one entry
    /// per dump of `genes orthology`.
    inputs: &'static [&'static str],
    /// The recording's `column_names`.
    column_names: &'static [&'static str],
    /// (record_count or row_count, returned_count, offset) when the recording has a window.
    window: Option<(u64, u64, u64)>,
}

include!(env!("MZLIB_REPLAY_TABLE"));

impl Recording {
    fn get(&self, key: &str) -> Option<&V> {
        self.fields
            .iter()
            .find(|(name, _)| *name == key)
            .map(|(_, value)| value)
    }
}

/// Python's `str()` of a recorded scalar, which is what the option value is compared with.
fn python_str(value: &V) -> String {
    match value {
        V::Null => "None".to_owned(),
        V::Bool(true) => "True".to_owned(),
        V::Bool(false) => "False".to_owned(),
        V::Text(text) => (*text).to_owned(),
        V::Compound => String::new(),
    }
}

fn base(path: &str) -> &str {
    path.rsplit(['/', '\\']).next().unwrap_or(path)
}

fn parse(argv: &[String]) -> (String, Vec<(String, Option<String>)>) {
    let mut i = 0;
    let mut verb = Vec::new();
    while i < argv.len() && !argv[i].starts_with("--") {
        verb.push(argv[i].clone());
        i += 1;
    }
    let mut options = Vec::new();
    while i < argv.len() {
        let name = argv[i][2..].to_owned();
        if i + 1 < argv.len() && !argv[i + 1].starts_with("--") {
            options.push((name, Some(argv[i + 1].clone())));
            i += 2;
        } else {
            options.push((name, None));
            i += 1;
        }
    }
    (verb.join(" "), options)
}

/// Wire options a verb echoes under another name, so a recording made with one value cannot answer
/// a call with another (`peptidoform fragments` echoes `--max-mods` as `max_modifications`).
const ECHOED_AS: &[(&str, &str)] = &[
    ("max-mods", "max_modifications"),
    ("max-isoforms", "max_modification_isoforms"),
    ("psms", "psm_file"),
    ("peptides", "peptides_file"),
    ("from", "source_format"),
    ("to", "target_format"),
];

/// Verbs whose input travels on stdin and is echoed in the recording, so a recording answers only a
/// call that sent the same input. stdin is read only for these.
const STDIN_ECHO: &[&str] = &[
    "quant flashlfq",
    "peptidoform convert",
    "stats combine-p",
    "stats fit",
    "stats fit-logistic",
    "stats fit-mixed",
    "diff design",
    "diff observations",
    "sdrf restrict",
    "proteins parse-accession",
    "proteins align",
    "genes orthology",
];

/// Recordings made with an option the call must send as well. The echo alone cannot say so: a call
/// that leaves the option out takes a default (mzLib's) that the replay does not know, and the echo
/// of a value that was given looks the same as the echo of a default. So a call without `--max-mods`
/// cannot be answered from the max-mods 1 albumin digest. pyMzLib's `RECORDED_WITH`, entry for
/// entry; `None` is a flag, sent without a value.
const RECORDED_WITH: &[(&str, &[(&str, Option<&str>)])] = &[
    ("peptidoform_P02768_small.json", &[("max-mods", Some("1"))]),
    (
        "proteins_align_charged.json",
        &[("end-gaps", Some("charged"))],
    ),
    // library read echoes both, but a call that leaves them out takes the defaults, and the echo
    // would not stop a recording made with them from answering it.
    ("library_read_fragments.json", &[("fragments", None)]),
    (
        "library_read_index_only.json",
        &[("index-only-threshold", Some("0"))],
    ),
    ("stats_fit_mixed_lmm_ml.json", &[("estimator", Some("ml"))]),
    (
        "stats_fit_malat_marginal.json",
        &[("estimator", Some("marginal_likelihood"))],
    ),
    (
        "diff_observations_spike_in_normalized.json",
        &[("min-peptides", Some("50")), ("normalize", None)],
    ),
    // ptmqtl: a flag or an option left out takes mzLib's default, which the echo cannot tell apart.
    (
        "ptmqtl_pairs_batch_run_level.json",
        &[("remove-run-level", None)],
    ),
    (
        "ptmqtl_trait_effects_strict.json",
        &[("min-replicates", Some("6"))],
    ),
    // sdrf isobaric-design echoes the plex source, but a call that declares none must get the
    // refusal recorded without one, not the design read with the pattern.
    (
        "sdrf_isobaric_design_PXD008841.json",
        &[("plex-file-name-pattern", Some(r"TMT_?pool(?<plex>\d+)"))],
    ),
];

/// Optional input files a verb echoes as `<name>_file`, null when the option was not given (stats
/// correlate --given, diff observations --protein-groups and --samples, ptmqtl --samples and
/// --run-level-groups, stats fit --contrasts). A recording made with one answers only a call that
/// gives one, and the reverse. pyMzLib's `OPTIONAL_FILES`.
const OPTIONAL_FILES: &[&str] = &[
    "given",
    "protein-groups",
    "samples",
    "run-level-groups",
    "contrasts",
];

/// Options a verb echoes inside an object rather than at the top level (stats fit --estimator as
/// prior.estimator, --confidence as contrasts.confidence_level), flattened by `build.rs` as
/// `outer.inner`. pyMzLib's `NESTED_ECHOES`.
const NESTED_ECHOES: &[(&str, &str)] = &[
    ("estimator", "prior.estimator"),
    ("confidence", "contrasts.confidence_level"),
];

/// Optional output files a verb reports in its own field, null when the option was not given
/// (diff observations --responses-out). pyMzLib's `OPTIONAL_OUTPUTS`.
const OPTIONAL_OUTPUTS: &[(&str, &str)] = &[("responses-out", "responses_written")];

/// Echoes that say whether an option was given, when its value echo cannot (stats entrapment-fdp
/// echoes the thresholds it used, which are mzLib's default when --thresholds was absent).
/// pyMzLib's `GIVEN_ECHOES`.
const GIVEN_ECHOES: &[(&str, &str)] = &[("thresholds_given", "thresholds")];

/// A numeric echo against an option's text: JSON writes 1.0 as 1, Rust sends `1.0`. pyMzLib's
/// `same_number`.
fn same_number(recorded: &str, value: &str) -> bool {
    match (recorded.parse::<f64>(), value.parse::<f64>()) {
        (Ok(a), Ok(b)) => a == b,
        _ => false,
    }
}

/// pyMzLib's `NESTED_ECHOES`, `GIVEN_ECHOES`, `OPTIONAL_FILES` and `OPTIONAL_OUTPUTS` checks: why
/// this recording cannot answer the call, or `None`.
fn echo_mismatch(recording: &Recording, options: &[(String, Option<String>)]) -> Option<String> {
    for (name, key) in NESTED_ECHOES {
        if let (Some(Some(value)), Some(V::Text(echo))) = (option(options, name), recording.get(key))
        {
            if echo != value && !same_number(echo, value) {
                return Some(format!("recorded with {key}='{echo}', not '{value}'"));
            }
        }
    }
    for (key, name) in GIVEN_ECHOES {
        if let Some(V::Bool(given)) = recording.get(key) {
            if option(options, name).is_some() != *given {
                let with = if *given { "with" } else { "without" };
                return Some(format!("recorded {with} --{name}"));
            }
        }
    }
    let optional_files = OPTIONAL_FILES
        .iter()
        .map(|name| (*name, format!("{}_file", name.replace('-', "_"))));
    let optional_outputs = OPTIONAL_OUTPUTS
        .iter()
        .map(|(name, key)| (*name, (*key).to_owned()));
    for (name, key) in optional_files.chain(optional_outputs) {
        if let Some(echo) = recording.get(&key) {
            let recorded_with = !matches!(echo, V::Null);
            if option(options, name).is_some() != recorded_with {
                let with = if recorded_with { "with" } else { "without" };
                return Some(format!("recorded {with} --{name}"));
            }
        }
    }
    None
}

/// Why this recording cannot answer a call that does not send what it was recorded with, or `None`.
fn recorded_with_mismatch(
    recording: &Recording,
    options: &[(String, Option<String>)],
) -> Option<String> {
    let (_, required) = RECORDED_WITH
        .iter()
        .find(|(fixture, _)| *fixture == recording.fixture)?;
    required.iter().find_map(|(name, value)| {
        let sent = option(options, name);
        let fits = match value {
            Some(value) => sent.and_then(|v| v.as_deref()) == Some(*value),
            None => matches!(sent, Some(None)),
        };
        (!fits).then(|| {
            let shown = value.map_or_else(String::new, |value| format!(" {value}"));
            format!("recorded with --{name}{shown}, which the call does not send")
        })
    })
}

fn option<'a>(options: &'a [(String, Option<String>)], name: &str) -> Option<&'a Option<String>> {
    options.iter().find(|(n, _)| n == name).map(|(_, v)| v)
}

/// Why this recording cannot be the answer to these options, or `None` if it can.
fn mismatch(recording: &Recording, options: &[(String, Option<String>)]) -> Option<String> {
    for (name, value) in options {
        let key = ECHOED_AS
            .iter()
            .find(|(wire, _)| wire == name)
            .map_or_else(|| name.replace('-', "_"), |(_, echoed)| (*echoed).to_owned());
        if matches!(name.as_str(), "limit" | "offset" | "out") {
            continue;
        }
        let Some(value) = value else {
            let flag = recording
                .get(&format!("{key}_included"))
                .or_else(|| recording.get(&key));
            match flag {
                None | Some(V::Null) | Some(V::Bool(true)) => {}
                Some(other) => {
                    return Some(format!(
                        "--{name} given, but the recording has {key}={}",
                        python_str(other)
                    ))
                }
            }
            continue;
        };
        // An input file option echoed back as <name>_file (--peptides -> peptides_file,
        // --responses -> responses_file): hold the call's file to the recording's, by file name.
        if recording.get(&key).is_none() {
            if let Some(V::Text(file)) = recording.get(&format!("{key}_file")) {
                if base(file) != base(value) {
                    return Some(format!(
                        "recorded from '{}', not '{}'",
                        base(file),
                        base(value)
                    ));
                }
                continue;
            }
        }
        let recorded = match recording.get(&key) {
            None | Some(V::Compound) => continue,
            Some(recorded) => python_str(recorded),
        };
        if key == "path" || key.ends_with("_file") {
            if base(&recorded) != base(value) {
                return Some(format!(
                    "recorded from '{}', not '{}'",
                    base(&recorded),
                    base(value)
                ));
            }
        } else if recorded != *value && !same_number(&recorded, value) {
            return Some(format!("recorded with {key}='{recorded}', not '{value}'"));
        }
    }

    if option(options, "paths-stdin").is_some() != recording.bulk {
        return Some(if recording.bulk {
            "a bulk (--paths-stdin) recording".to_owned()
        } else {
            "a one-document recording".to_owned()
        });
    }

    // A recording that wrote a file answers only a call that asked for one, and the reverse: the
    // payload's `written` block is the evidence, so an `out` example cannot print a recording that
    // wrote nothing.
    if let Some(written) = recording.get("written") {
        let wrote = !matches!(written, V::Null);
        if option(options, "out").is_some() != wrote {
            return Some(if wrote {
                "a recording that wrote out=".to_owned()
            } else {
                "a recording without out=".to_owned()
            });
        }
    }

    // proteins read: a filtered read (--accessions-stdin) and an unfiltered one never stand in for
    // each other, or an example would print rows its filter did not select.
    if let Some(filter) = recording.get("accession_filter_count") {
        if option(options, "accessions-stdin").is_some() == matches!(filter, V::Null) {
            return Some(
                "a recording with the other accession filter (filtered vs unfiltered)".to_owned(),
            );
        }
    }

    // --sequences adds a sequence column: a recording with it and one without never stand in for
    // each other, or an example would print sequences it did not ask for (or lack them).
    if recording.get("accession_filter_count").is_some()
        && option(options, "sequences").is_some() != recording.column_names.contains(&"sequence")
    {
        return Some("a recording with the other --sequences choice".to_owned());
    }

    // proteins align: --no-residue-map leaves residue_map null, and only then.
    if let Some(map) = recording.get("residue_map") {
        if option(options, "no-residue-map").is_some() != matches!(map, V::Null) {
            return Some("a recording with the other --no-residue-map choice".to_owned());
        }
    }

    // genes orthology: --pair-status FROM,TO fills pair_status for that direction, and only then;
    // and only a --restrict-gene-trees call says which gene sets its gene trees were restricted to.
    if let Some(status) = recording.get("pair_status") {
        let given = option(options, "pair-status");
        if given.is_some() == matches!(status, V::Null) {
            return Some("a recording with the other --pair-status choice".to_owned());
        }
        if let Some(Some(sent)) = given {
            let direction = format!(
                "{},{}",
                recording
                    .get("pair_status.from_species")
                    .map(python_str)
                    .unwrap_or_default(),
                recording
                    .get("pair_status.to_species")
                    .map(python_str)
                    .unwrap_or_default()
            );
            if *sent != direction {
                return Some(format!("recorded for --pair-status {direction}"));
            }
        }
    }
    if let Some(restricted) = recording.get("gene_trees.restricted_to_gene_sets") {
        if option(options, "restrict-gene-trees").is_some() == matches!(restricted, V::Null) {
            return Some("a recording with the other --restrict-gene-trees choice".to_owned());
        }
    }

    if let Some(kit) = recording.get("kit") {
        if !matches!(kit, V::Null) && option(options, "kit").is_none() {
            return Some(format!(
                "recorded for kit={}, but the call asks for every kit",
                python_str(kit)
            ));
        }
    }

    if let Some(ms_order) = recording.get("ms_order") {
        if !matches!(ms_order, V::Null) && option(options, "ms-order").is_none() {
            return Some(format!(
                "recorded with ms_order={}, but the call has no ms-order",
                python_str(ms_order)
            ));
        }
    }

    if let (Some((total, returned, recorded_offset)), None) =
        (recording.window, option(options, "out"))
    {
        let number = |name: &str| -> u64 {
            option(options, name)
                .and_then(|v| v.as_deref())
                .and_then(|v| v.parse().ok())
                .unwrap_or(0)
        };
        let offset = number("offset");
        let mut expected = total.saturating_sub(offset);
        if option(options, "limit").is_some() {
            expected = expected.min(number("limit"));
        }
        if recorded_offset != offset || returned != expected {
            return Some(format!(
                "recorded window offset={recorded_offset} returned={returned} of {total}, not \
                 what limit/offset ask for"
            ));
        }
    }

    if matches!(recording.get("peaks_included"), Some(V::Bool(true)))
        && option(options, "peaks").is_none()
    {
        return Some("recorded with peaks, but the call did not ask for them".to_owned());
    }
    None
}

fn usage(message: &str) -> String {
    format!(
        "{{\"ok\":false,\"data\":null,\"error\":{{\"type\":\"usage\",\"message\":{}}}}}",
        json_string(&format!("replay bridge: {message}"))
    )
}

fn json_string(text: &str) -> String {
    let mut out = String::from('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn describe(options: &[(String, Option<String>)]) -> String {
    options
        .iter()
        .map(|(name, value)| match value {
            Some(value) => format!("--{name} {value}"),
            None => format!("--{name}"),
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// For `STDIN_ECHO` verbs: the call's stdin must be the recording's input. `quant flashlfq`: the
/// same runs (the first tab-separated cell of each stdin line), compared by file name.
/// `peptidoform convert`: the same sequences, in the same order, as the recording's `input` column.
fn stdin_mismatch(recording: &Recording, stdin: &str, verb: &str) -> Option<String> {
    if let Some(why) = model_stdin_mismatch(recording, stdin, verb) {
        return why;
    }
    // sdrf restrict: the same searched files, in the same order, as the recording's
    // searched_files, and as many acquired names as its acquired_name_count.
    if verb == "sdrf restrict" {
        let lines: Vec<Vec<&str>> = stdin
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| line.split('\t').collect())
            .collect();
        let sent: Vec<&str> = lines.iter().map(|fields| fields[0].trim()).collect();
        if sent != recording.inputs {
            return Some(format!(
                "recorded for searched files {:?}, not {sent:?}",
                recording.inputs
            ));
        }
        let acquired = lines
            .iter()
            .filter(|fields| fields.len() > 1 && !fields[1].trim().is_empty())
            .count();
        let recorded = recording.get("acquired_name_count").map(python_str);
        return (recorded.as_deref() != Some(acquired.to_string().as_str())).then(|| {
            format!(
                "recorded with {} acquired names, not {acquired}",
                recorded.unwrap_or_default()
            )
        });
    }
    // proteins align: the same pairs, in the same order, as the recording's gapless aligned rows.
    if verb == "proteins align" {
        let sent: Vec<&str> = stdin.lines().filter(|line| !line.trim().is_empty()).collect();
        return (sent != recording.inputs).then(|| {
            let short: Vec<String> = recording
                .inputs
                .iter()
                .map(|pair| pair.chars().take(12).collect())
                .collect();
            format!("recorded for pairs starting {short:?}, not these")
        });
    }
    // genes orthology: the same species (the first field of each line before `--`) as the
    // recording's species, and as many dump paths after it as the recording has dumps.
    if verb == "genes orthology" {
        let lines: Vec<&str> = stdin
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .collect();
        let cut = lines.iter().position(|line| *line == "--").unwrap_or(lines.len());
        let mut species: Vec<&str> = lines[..cut]
            .iter()
            .map(|line| line.split('\t').next().unwrap_or("").trim())
            .collect();
        species.sort_unstable();
        let recorded_cut = recording
            .inputs
            .iter()
            .position(|line| *line == "--")
            .unwrap_or(recording.inputs.len());
        if species != recording.inputs[..recorded_cut] {
            return Some(format!(
                "recorded for species {:?}, not {species:?}",
                &recording.inputs[..recorded_cut]
            ));
        }
        let dumps = lines.len().saturating_sub(cut + 1);
        let recorded_dumps = recording.inputs.len().saturating_sub(recorded_cut + 1);
        return (dumps != recorded_dumps)
            .then(|| format!("recorded with {recorded_dumps} dumps, not {dumps}"));
    }
    if verb == "peptidoform convert" || verb == "proteins parse-accession" {
        let sent: Vec<&str> = stdin.lines().filter(|line| !line.trim().is_empty()).collect();
        return (sent != recording.inputs)
            .then(|| format!("recorded for sequences {:?}, not {sent:?}", recording.inputs));
    }
    let mut sent: Vec<&str> = stdin
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| base(line.split('\t').next().unwrap_or("").trim()))
        .collect();
    let mut recorded: Vec<&str> = recording.runs.to_vec();
    sent.sort_unstable();
    sent.dedup();
    recorded.sort_unstable();
    recorded.dedup();
    (sent != recorded).then(|| format!("recorded for runs {recorded:?}, not {sent:?}"))
}

/// The non-blank stdin lines, each split into trimmed tab-separated fields.
fn stdin_fields(stdin: &str) -> Vec<Vec<&str>> {
    stdin
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split('\t').map(str::trim).collect())
        .collect()
}

fn recorded_text<'a>(recording: &'a Recording, key: &str) -> &'a str {
    match recording.get(key) {
        Some(V::Text(text)) => text,
        _ => "",
    }
}

/// A `diff design` call's stdin in `build.rs`'s `recorded_model` form. pyMzLib's `diff_model`.
fn canonical_model(stdin: &str) -> String {
    let number = |text: &str| {
        text.parse::<f64>()
            .map_or_else(|_| String::new(), |n| format!("{n:?}"))
    };
    let field = |fields: &[&str], i: usize| fields.get(i).copied().unwrap_or("").to_owned();
    let mut kinds: [Vec<String>; 5] = Default::default();
    for f in stdin_fields(stdin) {
        let (kind, line) = match f[0] {
            "factor" => (0, format!("factor\t{}\t{}", field(&f, 1), field(&f, 2))),
            "covariate" => (
                1,
                format!(
                    "covariate\t{}\t{}\t{}\t{}\t{}",
                    field(&f, 1),
                    field(&f, 2),
                    number(&field(&f, 3)),
                    field(&f, 4),
                    number(&field(&f, 5))
                ),
            ),
            "interaction" => (2, format!("interaction\t{}\t{}", field(&f, 1), field(&f, 2))),
            "stratum" => (3, format!("stratum\t{}", field(&f, 1))),
            "batch" => (4, format!("batch\t{}", field(&f, 1))),
            _ => continue,
        };
        kinds[kind].push(line);
    }
    kinds.concat().join("\n")
}

/// The stdin echoes of the stats and diff verbs, or `None` for another verb. pyMzLib's
/// `stdin_mismatch` branches for them:
///
/// * stats fit, fit-logistic and fit-mixed: the same coefficients, in the same order, as the
///   recording's `tested_coefficients` (none, for a stats fit of contrasts alone);
/// * stats combine-p: the same features in the same first-appearance order, the same number of
///   tests, and weights and effects given or not as the recording's `weighted` and `signed` say;
/// * diff design: the same model;
/// * diff observations: the same runs, or, with `--samples`, runs among those sent.
fn model_stdin_mismatch(recording: &Recording, stdin: &str, verb: &str) -> Option<Option<String>> {
    let why = match verb {
        "stats fit" | "stats fit-logistic" | "stats fit-mixed" => {
            let sent: Vec<&str> = stdin
                .lines()
                .map(str::trim)
                .filter(|line| !line.is_empty())
                .collect();
            let recorded = recorded_text(recording, "@tested_coefficients");
            (sent.join("\n") != recorded)
                .then(|| format!("recorded for coefficients {recorded:?}, not {sent:?}"))
        }
        "stats combine-p" => {
            let lines = stdin_fields(stdin);
            let mut features: Vec<&str> = Vec::new();
            for fields in &lines {
                if !features.contains(&fields[0]) {
                    features.push(fields[0]);
                }
            }
            let weighted = lines.iter().any(|f| {
                f.len() >= 3 && !["", "na", "nan"].contains(&f[2].to_lowercase().as_str())
            });
            let signed = lines.iter().any(|f| f.len() == 4);
            let flag = |key: &str| matches!(recording.get(key), Some(V::Bool(true)));
            if features.join("\n") != recorded_text(recording, "@features") {
                Some(format!("recorded for other features, not {features:?}"))
            } else if lines.len().to_string() != recorded_text(recording, "line_count") {
                Some(format!(
                    "recorded for {} tests, not {}",
                    recorded_text(recording, "line_count"),
                    lines.len()
                ))
            } else if (weighted, signed) != (flag("weighted"), flag("signed")) {
                Some(format!(
                    "recorded with weighted={} signed={}, not weighted={weighted} signed={signed}",
                    flag("weighted"),
                    flag("signed")
                ))
            } else {
                None
            }
        }
        "diff design" => {
            let sent = canonical_model(stdin);
            let recorded = recorded_text(recording, "@model");
            (sent != recorded).then(|| format!("recorded for the model {recorded:?}, not {sent:?}"))
        }
        "diff observations" => {
            let mut sent: Vec<String> = stdin_fields(stdin)
                .iter()
                .map(|f| {
                    let whole = |i: usize| f.get(i).and_then(|v| v.parse::<i64>().ok()).unwrap_or(0);
                    format!(
                        "{}\t{}\t{}\t{}",
                        f[0],
                        f.get(1).copied().unwrap_or(""),
                        whole(2),
                        whole(3)
                    )
                })
                .collect();
            sent.sort();
            sent.dedup();
            let recorded: Vec<&str> = recorded_text(recording, "@runs")
                .split('\n')
                .filter(|l| !l.is_empty())
                .collect();
            // --samples narrows the recorded runs to one stratum's, so they need only be among
            // those sent.
            let narrowed = !matches!(recording.get("samples_file"), None | Some(V::Null));
            let fits = if narrowed {
                recorded.iter().all(|r| sent.iter().any(|s| s == r))
            } else {
                recorded == sent
            };
            (!fits).then(|| format!("recorded for other runs, not {:?}", &sent[..sent.len().min(3)]))
        }
        _ => return None,
    };
    Some(why)
}

fn answer(argv: &[String], stdin: &str) -> (bool, String) {
    let (verb, options) = parse(argv);
    let Some((_, candidates)) = TABLE.iter().find(|(name, _)| *name == verb) else {
        return (
            false,
            usage(&format!("no fixture is recorded for '{verb}'")),
        );
    };

    let mut fits = Vec::new();
    let mut reasons = Vec::new();
    for recording in *candidates {
        let why = echo_mismatch(recording, &options)
            .or_else(|| mismatch(recording, &options))
            .or_else(|| recorded_with_mismatch(recording, &options))
            .or_else(|| {
                if STDIN_ECHO.contains(&verb.as_str()) {
                    stdin_mismatch(recording, stdin, &verb)
                } else {
                    None
                }
            });
        match why {
            Some(why) => reasons.push(format!("{}: {why}", recording.fixture)),
            None => fits.push(recording),
        }
    }
    match fits.as_slice() {
        [one] => (true, one.envelope.to_owned()),
        [] => (
            false,
            usage(&format!(
                "no recording of '{verb}' fits [{}]: {}",
                describe(&options),
                reasons.join("; ")
            )),
        ),
        several => (
            false,
            usage(&format!(
                "'{verb}' [{}] fits several recordings: {}",
                describe(&options),
                several
                    .iter()
                    .map(|r| r.fixture)
                    .collect::<Vec<_>>()
                    .join(", ")
            )),
        ),
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let (verb, _) = parse(&argv);
    let mut stdin = String::new();
    if STDIN_ECHO.contains(&verb.as_str()) {
        use std::io::Read;
        let _ = std::io::stdin().read_to_string(&mut stdin);
    }
    let (ok, envelope) = answer(&argv, &stdin);
    let mut stdout = std::io::stdout();
    let _ = stdout.write_all(envelope.as_bytes());
    let _ = stdout.flush();
    std::process::exit(if ok { 0 } else { 2 });
}
