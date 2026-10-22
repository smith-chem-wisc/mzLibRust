//! Builds the stand-in bridge: reads the recordings, works out what call each one can answer, and
//! compiles `src/stub.rs` into an executable with that table baked in.
//!
//! Which recordings a wire verb may answer from comes from the `examples` of its spec in
//! `docs/specs/`, plus [`REPLAY_EXTRA`] for recordings no spec lists yet. That is the rule
//! pyMzLib's `pkg/src/conftest.py` follows, so a recording that is a spec's example answers the same
//! call in every binding.
//!
//! All JSON handling happens here, at build time, so the executable is std-only and compiles in
//! well under a second. It is compiled by invoking the same `rustc` that is building this crate,
//! because Cargo cannot build a dependency's binary target for a doctest to run.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

/// Recordings for verbs whose spec does not list them as an example, or that have no spec yet.
///
/// Mirrors `REPLAY_EXTRA` in pyMzLib's `pkg/python/src/conftest.py`, plus the small recordings
/// this crate's function examples have always replayed. When a verb's spec lands with one of these
/// among its `examples`, delete it here.
const REPLAY_EXTRA: &[(&str, &[&str])] = &[
    (
        "readers read-spectra",
        &["readers_spectra_ms2.json", "readers_spectra_peaks.json"],
    ),
    (
        "readers read-protein-groups",
        &["readers_protein_groups_rna.json"],
    ),
    ("readers identify", &["readers_identify_fragger.json"]),
    ("readers read-features", &["readers_features_topfd.json"]),
    (
        "readers read-matches",
        &[
            "readers_matches_casanovo.json",
            "readers_matches_mspathfinder.json",
            "readers_matches_mzid.json",
            // --scores: the bridge's checker cannot hold a present_when column to a spec example.
            "readers_matches_mzid_scores.json",
        ],
    ),
    // The identify bulk shape is not a spec example: check_verbs.py has no table-less bulk verb.
    ("readers identify", &["readers_many_identify.json"]),
    ("version", &["bridge_version.json"]),
    // The PRIDE guide: PRIDE's empty answer for an unknown accession. The complete FTP listing it
    // also reads is the `pride ftp-files` spec's example.
    (
        "pride files",
        &[
            "pride_PXD000001_files.json",
            "pride_PXD999999999_files.json",
        ],
    ),
    // The peptidoforms guide: MetaMorpheus's BottomUpExample.psmtsv (mzLib's test data), whose
    // full sequences it converts to Unimod and ProForma.
    ("readers read-results", &["readers_results_psmtsv.json"]),
    ("pride search", &["pride_search_plasmodium.json"]),
    // The FlashLFQ guide: mzLib's two K562 test runs with MBR, and median polish of the peptide
    // table that run wrote (renamed K562_QuantifiedPeptides.tsv so it is told apart).
    (
        "quant flashlfq",
        &["flashlfq_small.json", "flashlfq_k562_mbr.json"],
    ),
    (
        "quant median-polish",
        &["median_polish_small.json", "median_polish_k562.json"],
    ),
    // The full albumin digest the peptidoforms guide reads.
    (
        "peptidoform fragments",
        &["peptidoform_P02768_small.json", "peptidoform_P02768.json"],
    ),
    // The SDRF guide: assess on a Partial file and on a skeleton.
    (
        "sdrf assess",
        &["sdrf_assess_PXD000070.json", "sdrf_assess_skeleton.json"],
    ),
    // The proteins guide: three databases and a contaminant, filtered to four accessions; and the
    // DYNC1I2 isoforms and both albumins with their sequences, which the alignment section aligns.
    (
        "proteins read",
        &[
            "proteins_read_accessions.json",
            "proteins_read_sequences.json",
        ],
    ),
    // The diff and stats guides: stats fit on the design file diff design wrote for the MALAT1
    // dilution series, without and with its contrasts file; and that series' design with no
    // reference level (every pair of loads a contrast), which the contrasts example reads.
    (
        "stats fit",
        &[
            "stats_fit_malat_diff_design.json",
            "stats_fit_malat_diff_design_contrasts.json",
        ],
    ),
    ("diff design", &["diff_design_malat_pairs.json"]),
];

fn main() {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let root = manifest.join("..").join("..");
    let specs = root.join("docs").join("specs");
    let fixtures = root.join("tests").join("fixtures");
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=src/stub.rs");
    println!("cargo:rerun-if-changed={}", specs.display());
    println!("cargo:rerun-if-changed={}", fixtures.display());

    let table = replay_table(&specs);
    let recordings = out.join("recordings");
    std::fs::create_dir_all(&recordings).unwrap();

    let mut code = String::from(
        "// Generated by tools/replay-bridge/build.rs. Do not edit.\n\
         pub static TABLE: &[(&str, &[Recording])] = &[\n",
    );
    for (verb, names) in &table {
        writeln!(code, "    ({verb:?}, &[").unwrap();
        for name in names {
            let path = fixtures.join(name);
            let text = std::fs::read_to_string(&path).unwrap_or_else(|error| {
                panic!(
                    "the replay table names fixture {name} for '{verb}', which cannot be read at \
                     {}: {error}",
                    path.display()
                )
            });
            let parsed: Value = serde_json::from_str(&text)
                .unwrap_or_else(|error| panic!("{name} is not JSON: {error}"));
            // A fixture is either the bare `data` payload or a whole envelope.
            let data = match &parsed {
                Value::Object(map) if map.contains_key("ok") && map.contains_key("data") => {
                    map["data"].clone()
                }
                other => other.clone(),
            };
            let envelope = serde_json::json!({"ok": true, "data": data, "error": null});
            let envelope_path = recordings.join(name);
            std::fs::write(&envelope_path, envelope.to_string()).unwrap();
            code.push_str(&recording(name, &data, &envelope_path));
        }
        code.push_str("    ]),\n");
    }
    code.push_str("];\n");

    let table_path = out.join("table.rs");
    std::fs::write(&table_path, code).unwrap();

    compile_stub(&manifest.join("src").join("stub.rs"), &table_path, &out);
}

/// Wire verb -> the recordings it may answer from, in a stable order.
fn replay_table(specs: &Path) -> BTreeMap<String, Vec<String>> {
    let mut table: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut entries: Vec<PathBuf> = std::fs::read_dir(specs)
        .unwrap_or_else(|error| panic!("cannot read {}: {error}", specs.display()))
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "yaml"))
        .collect();
    entries.sort();

    for path in entries {
        // A line scan rather than a YAML parser: the two keys needed are one per line in every
        // spec (`verb:` at the top level, `fixture:` in each example), and a build dependency on a
        // YAML crate would be the only reason to have one.
        let text = std::fs::read_to_string(&path).unwrap();
        let mut verb = None;
        for line in text.lines() {
            if let Some(rest) = line.strip_prefix("verb:") {
                verb = Some(rest.split('#').next().unwrap().trim().to_owned());
            } else if let Some(rest) = ["- fixture:", "- {fixture:", "fixture:"]
                .iter()
                .find_map(|key| line.trim_start().strip_prefix(key))
            {
                let name = rest
                    .split([',', '}'])
                    .next()
                    .unwrap()
                    .trim()
                    .trim_matches(|c| c == '"' || c == '\'');
                let verb = verb
                    .clone()
                    .unwrap_or_else(|| panic!("{}: fixture before verb", path.display()));
                push_unique(table.entry(verb).or_default(), name);
            }
        }
    }
    for (verb, names) in REPLAY_EXTRA {
        for name in *names {
            push_unique(table.entry((*verb).to_owned()).or_default(), name);
        }
    }
    table
}

fn push_unique(list: &mut Vec<String>, name: &str) {
    if !list.iter().any(|existing| existing == name) {
        list.push(name.to_owned());
    }
}

/// One `Recording { .. }` literal: the envelope, plus the facts `stub.rs` checks a call against.
fn recording(name: &str, data: &Value, envelope: &Path) -> String {
    let empty = serde_json::Map::new();
    let map = data.as_object().unwrap_or(&empty);

    let literal = |value: &Value| match value {
        Value::Null => "V::Null".to_owned(),
        Value::Bool(flag) => format!("V::Bool({flag})"),
        Value::String(text) => format!("V::Text({text:?})"),
        Value::Number(number) => format!("V::Text({:?})", number.to_string()),
        Value::Array(_) | Value::Object(_) => "V::Compound".to_owned(),
    };
    let mut fields = String::new();
    for (key, value) in map {
        write!(fields, "({key:?}, {}), ", literal(value)).unwrap();
        // One level of an echoed object as `outer.inner`, for the rules that read inside one
        // (`genes orthology`'s `pair_status.from_species`). A dotted key names no option, so the
        // generic option check never meets it.
        if let Value::Object(inner) = value {
            for (sub, value) in inner {
                write!(fields, "(\"{key}.{sub}\", {}), ", literal(value)).unwrap();
            }
        }
    }
    for (key, text) in derived_fields(map) {
        write!(fields, "({key:?}, V::Text({text:?})), ").unwrap();
    }

    // BULK.md: a --paths-stdin call has its own envelope. A recording is bulk when it carries the
    // per-input files[] (entries with a path) and no top-level path. `pride files` also has a
    // files list, of PRIDE files, which carry no path.
    // An empty list says nothing either way (PRIDE's answer for an unknown accession is files=[]),
    // so it counts as bulk only when BULK.md's read_count says so.
    let bulk = match map.get("files") {
        Some(Value::Array(files)) => {
            (!files.is_empty() || map.contains_key("read_count"))
                && files
                    .iter()
                    .all(|f| f.as_object().is_some_and(|f| f.contains_key("path")))
                && !map.contains_key("path")
        }
        _ => false,
    };

    // Python's data.get("record_count", data.get("row_count")): a present-but-null record_count
    // wins over row_count.
    let total = match map.get("record_count") {
        Some(value) => value.as_u64(),
        None => map.get("row_count").and_then(Value::as_u64),
    };
    let window = match (total, map.get("returned_count").and_then(Value::as_u64)) {
        (Some(total), Some(returned)) => {
            let offset = map.get("offset").and_then(Value::as_u64).unwrap_or(0);
            format!("Some(({total}, {returned}, {offset}))")
        }
        _ => "None".to_owned(),
    };

    // The runs a `quant flashlfq` recording was made from, by file name (stub.rs STDIN_ECHO).
    let runs: Vec<String> = match map.get("spectra_files") {
        Some(Value::Array(files)) => files
            .iter()
            .filter_map(|f| f.get("full_path").and_then(Value::as_str))
            .map(|path| path.rsplit(['/', '\\']).next().unwrap_or(path).to_owned())
            .collect(),
        _ => Vec::new(),
    };

    // The input a STDIN_ECHO verb's recording was made from, one entry per stdin line (stub.rs):
    // the `input` column of a `peptidoform convert` or `proteins parse-accession` recording; the
    // searched files of an `sdrf restrict` one (its `searched_files`, as sent); the pairs of a
    // `proteins align` one, `source<TAB>target`, from its aligned rows with the gaps removed
    // (mzLib refuses a '-' in an input sequence, so removing them is exact); and for
    // `genes orthology`, its `species`, a `--`, then one entry per dump.
    let strings = |value: Option<&Value>| -> Vec<String> {
        match value {
            Some(Value::Array(rows)) => rows
                .iter()
                .filter_map(Value::as_str)
                .map(str::to_owned)
                .collect(),
            _ => Vec::new(),
        }
    };
    let column = |name: &str| strings(map.get("columns").and_then(|c| c.get(name)));
    let inputs: Vec<String> =
        if map.contains_key("source_format") || map.contains_key("namespace_counts") {
            column("input")
        } else if map.contains_key("searched_without_row") {
            strings(map.get("searched_files"))
        } else if map.contains_key("aligner_id") {
            column("aligned_source")
                .iter()
                .zip(column("aligned_target"))
                .map(|(s, t)| format!("{}\t{}", s.replace('-', ""), t.replace('-', "")))
                .collect()
        } else if map.contains_key("snapshot_id") {
            let mut lines = strings(map.get("species"));
            lines.push("--".to_owned());
            if let Some(Value::Object(dumps)) = map.get("dumps") {
                lines.extend(dumps.keys().cloned());
            }
            lines
        } else {
            Vec::new()
        };

    // The recording's column order, for the rules that ask whether a column is there
    // (`proteins read --sequences`).
    let column_names = strings(map.get("column_names"));

    format!(
        "        Recording {{ fixture: {name:?}, envelope: include_str!({path:?}), fields: &[{fields}], \
         bulk: {bulk}, runs: &{runs:?}, inputs: &{inputs:?}, column_names: &{column_names:?}, \
         window: {window} }},\n",
        path = envelope.display().to_string(),
    )
}

/// Facts `stub.rs` checks that are not top-level scalars, flattened into text under keys no wire
/// option can produce (they start with `@`). One level of an echoed object, `outer.inner`, is
/// written by `recording` itself: pyMzLib's `NESTED_ECHOES` reads stats fit's `--estimator` as
/// `prior.estimator` and `--confidence` as `contrasts.confidence_level` there.
///
/// * `@tested_coefficients`, one per line: the coefficients a `stats fit`, `fit-logistic` or
///   `fit-mixed` recording was made with (empty for a fit of contrasts alone);
/// * `@features`, one per line, for a `stats combine-p` recording: its features in order;
/// * `@model`, a `diff design` recording's model, in the canonical form `stub.rs` also builds from a
///   call's stdin;
/// * `@runs`, a `diff observations` recording's runs as sorted `file, sample, fraction, technical
///   replicate` lines.
fn derived_fields(map: &serde_json::Map<String, Value>) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let strings = |value: Option<&Value>| -> Option<Vec<String>> {
        value?
            .as_array()?
            .iter()
            .map(|v| v.as_str().map(str::to_owned))
            .collect()
    };
    if let Some(coefficients) = strings(map.get("tested_coefficients")) {
        out.push(("@tested_coefficients".to_owned(), coefficients.join("\n")));
    }
    if map.contains_key("weighted") && map.contains_key("signed") {
        let features =
            strings(map.get("columns").and_then(|c| c.get("feature"))).unwrap_or_default();
        out.push(("@features".to_owned(), features.join("\n")));
    }
    if map.contains_key("factors") && map.contains_key("stratum_columns") {
        out.push(("@model".to_owned(), recorded_model(map)));
    }
    if let Some(Value::Array(runs)) = map.get("runs") {
        let mut lines: Vec<String> = runs
            .iter()
            .filter_map(|run| {
                let text = |key: &str| {
                    run.get(key)
                        .and_then(Value::as_str)
                        .unwrap_or("")
                        .to_owned()
                };
                let whole = |key: &str| run.get(key).and_then(Value::as_i64).unwrap_or(0);
                run.get("file_name")?;
                Some(format!(
                    "{}\t{}\t{}\t{}",
                    text("file_name"),
                    text("sample_id"),
                    whole("fraction"),
                    whole("technical_replicate")
                ))
            })
            .collect();
        lines.sort();
        lines.dedup();
        out.push(("@runs".to_owned(), lines.join("\n")));
    }
    out
}

/// A `diff design` recording's model as `stub.rs`'s `canonical_model` writes a call's stdin: one
/// line per term, factors, then covariates, interactions, strata and batch, numbers as Rust's
/// shortest round-trip form. pyMzLib's `diff_model`, as text.
fn recorded_model(map: &serde_json::Map<String, Value>) -> String {
    let text = |item: &Value, key: &str| {
        item.get(key)
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_owned()
    };
    let number = |item: &Value, key: &str| {
        item.get(key)
            .and_then(Value::as_f64)
            .map_or_else(String::new, |n| format!("{n:?}"))
    };
    let list = |key: &str| {
        map.get(key)
            .and_then(Value::as_array)
            .cloned()
            .unwrap_or_default()
    };
    let mut lines = Vec::new();
    for f in list("factors") {
        lines.push(format!(
            "factor\t{}\t{}",
            text(&f, "name"),
            text(&f, "reference_level")
        ));
    }
    for c in list("covariates") {
        lines.push(format!(
            "covariate\t{}\t{}\t{}\t{}\t{}",
            text(&c, "name"),
            text(&c, "unit"),
            number(&c, "scale"),
            text(&c, "scaled_unit"),
            number(&c, "centre")
        ));
    }
    for i in list("interactions") {
        lines.push(format!(
            "interaction\t{}\t{}",
            text(&i, "first"),
            text(&i, "second")
        ));
    }
    for s in list("stratum_columns") {
        lines.push(format!("stratum\t{}", s.as_str().unwrap_or("")));
    }
    if let Some(batch) = map.get("batch_column").and_then(Value::as_str) {
        lines.push(format!("batch\t{batch}"));
    }
    lines.join("\n")
}

fn compile_stub(source: &Path, table: &Path, out: &Path) {
    let exe = out.join(if cfg!(windows) {
        "mzlib-replay-bridge.exe"
    } else {
        "mzlib-replay-bridge"
    });
    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".to_owned());
    let target = std::env::var("TARGET").unwrap();

    let status = Command::new(rustc)
        .args(["--edition=2021", "--crate-name", "mzlib_replay_bridge"])
        .args(["--crate-type=bin", "--cap-lints=allow"])
        .args(["-C", "opt-level=0", "-C", "debuginfo=0"])
        .args(["--target", &target])
        .arg(source)
        .arg("-o")
        .arg(&exe)
        .env("MZLIB_REPLAY_TABLE", table)
        .status()
        .expect("could not run rustc to build the replay bridge");
    assert!(status.success(), "rustc failed to build the replay bridge");

    println!("cargo:rustc-env=MZLIB_REPLAY_BRIDGE={}", exe.display());
}
