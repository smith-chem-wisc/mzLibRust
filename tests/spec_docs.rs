//! The reference facts, rendered from the bridge's per-verb specs, and the docs held to them.
//!
//! Every wire verb has one YAML spec in the bridge repository (`design/verbs/`), vendored here as
//! `docs/specs/` by `scripts/sync-specs.ps1`. A spec owns the **facts** about its verb: parameters
//! with their units and ranges, result fields with their units and what a null means, error kinds,
//! caveats, the mzLib code it wraps, citations, the spelling in each binding, and `since`. This
//! crate owns the Rust idiom and the prose. This file keeps the two honest:
//!
//! * [`rendered_fragments_are_current`] renders one set of Markdown fragments per spec into
//!   `docs/reference/`, which the functions include with `#[doc = include_str!(...)]`. It fails
//!   when a committed fragment differs from what its spec renders to. Regenerate with
//!   `MZLIB_RENDER_SPEC_DOCS=1 cargo test --test spec_docs`.
//! * [`every_param_and_field_is_documented_with_its_unit`] reads this crate's source with `syn`,
//!   and fails when a spec parameter is not a documented field of the options struct (or an
//!   argument of the function), when a result field is not a documented field of the returned
//!   struct, or when either is documented without its spec unit. So "Skip this many" cannot ship
//!   without saying *scans*.
//!
//! A deliberate difference from a spec - another name, or a fact carried somewhere other than a
//! struct field - is declared in [`DEVIATIONS`] with a reason, and nowhere else. A declaration that
//! names nothing in its spec fails [`every_deviation_names_a_real_param_or_field`], so the table
//! cannot rot. See `docs/reference-facts.md`.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::PathBuf;

use serde_json::Value;

// =================================================================================================
// Declarations: every place the Rust projection differs from a spec, and why
// =================================================================================================

/// One declared difference between a spec and this crate.
struct Deviation {
    /// The spec's `verb`, or `"*"` for a convention that applies to every verb.
    verb: &'static str,
    /// `param.<wire name>`, `field.<wire name>`, `bulk.param.<wire name>`, `bulk.field.<wire name>`,
    /// `table.<table>.<wire name>`, or `name` / `bulk` for the function itself.
    key: &'static str,
    /// The Rust spelling, or `None` when the fact is not a named field or argument at all.
    rust: Option<&'static str>,
    /// Why. Required.
    why: &'static str,
}

const fn dev(
    verb: &'static str,
    key: &'static str,
    rust: Option<&'static str>,
    why: &'static str,
) -> Deviation {
    Deviation {
        verb,
        key,
        rust,
        why,
    }
}

/// Every deliberate difference between the specs and this crate, one list per module so that a
/// change to one module's projection touches only its own list. Nothing else may be skipped.
const DEVIATIONS: &[&[Deviation]] = &[COMMON, READERS, SDRF, PROTEINS, PRIDE, PEPTIDOFORM, QUANT];

/// Deviations: conventions that hold for every verb.
const COMMON: &[Deviation] = &[
    dev(
        "*",
        "field.column_names",
        Some("columns"),
        "a columnar result carries its column order inside the Table: Table::names()",
    ),
    dev(
        "*",
        "bulk.field.column_names",
        Some("columns"),
        "a many-file result carries its column order inside the Table: Table::names()",
    ),
    dev(
        "*",
        "bulk.param.paths-stdin",
        None,
        "the `_many` function is the bulk form: it takes the list, and sends it on stdin",
    ),
    dev(
        "*",
        "field.error",
        None,
        "a single-path call that cannot read its file returns Err(MzLibError) instead; the wire \
         field is always null there, and only a many-file result's per-file entries carry one",
    ),
];

/// Deviations: readers.
const READERS: &[Deviation] = &[
    dev(
        "readers formats",
        "field.format_count",
        None,
        "formats() returns Vec<Format>; the count is its len()",
    ),
    dev(
        "readers formats",
        "field.formats",
        None,
        "formats() returns this list itself, as Format values",
    ),
];

/// Deviations: sdrf.
const SDRF: &[Deviation] = &[
    dev(
        "sdrf read",
        "field.column_names",
        Some("columns"),
        "SDRF column names repeat, so the header is a Vec<String> named columns, not a Table",
    ),
    dev(
        "sdrf pool",
        "field.column_names",
        Some("columns"),
        "SDRF column names repeat, so the header is a Vec<String> named columns, not a Table",
    ),
    dev(
        "sdrf validate",
        "param.paths-stdin",
        None,
        "validate_many(paths, options) is the many-document form; passing a list sets the flag",
    ),
    dev(
        "sdrf validate",
        "bulk.field.column_names",
        Some("columns"),
        "a columnar result carries its column order inside the Table: Table::names()",
    ),
    dev(
        "sdrf assess",
        "param.paths-stdin",
        None,
        "assess_many(paths, options) is the many-document form; passing a list sets the flag",
    ),
    dev(
        "sdrf assess",
        "bulk.field.column_names",
        Some("columns"),
        "a columnar result carries its column order inside the Table: Table::names()",
    ),
    dev(
        "sdrf samples",
        "param.paths-stdin",
        None,
        "samples_many(paths, options) is the many-document form; passing a list sets the flag",
    ),
    dev(
        "sdrf samples",
        "bulk.field.column_names",
        Some("columns"),
        "a columnar result carries its column order inside the Table: Table::names()",
    ),
    dev(
        "sdrf lint",
        "param.stdin",
        Some("documents"),
        "the stdin lines are rendered from lint_labelled's documents, one path[TAB label] per line",
    ),
    dev(
        "sdrf parse-age",
        "param.stdin",
        Some("cells"),
        "one stdin line per element of parse_ages' cells",
    ),
];

/// Deviations: proteins and genes.
const PROTEINS: &[Deviation] = &[
    // The three verbs take their databases the same way: a slice argument, plus contaminants in
    // the options, with --path / --paths-stdin / --contaminant chosen from the list's shape.
    dev(
        "proteins read",
        "param.path",
        Some("databases"),
        "one database or a list, as a slice argument; a list travels as --paths-stdin",
    ),
    dev(
        "proteins read",
        "param.contaminant",
        Some("contaminants"),
        "contaminant databases are their own list in the options, tagged on stdin or sent as --contaminant",
    ),
    dev(
        "proteins read",
        "param.paths-stdin",
        None,
        "chosen from the number of databases: more than one travels on stdin",
    ),
    dev(
        "genes resolve",
        "param.path",
        Some("databases"),
        "one database or a list, as a slice argument; a list travels as --paths-stdin",
    ),
    dev(
        "genes resolve",
        "param.contaminant",
        Some("contaminants"),
        "contaminant databases are their own list in the options, tagged on stdin or sent as --contaminant",
    ),
    dev(
        "genes resolve",
        "param.paths-stdin",
        None,
        "chosen from the number of databases: more than one travels on stdin",
    ),
    dev(
        "proteins classify-peptides",
        "param.path",
        Some("databases"),
        "one database or a list, as a slice argument; a list travels as --paths-stdin",
    ),
    dev(
        "proteins classify-peptides",
        "param.contaminant",
        Some("contaminants"),
        "contaminant databases are their own list in the options, tagged on stdin or sent as --contaminant",
    ),
    dev(
        "proteins classify-peptides",
        "param.paths-stdin",
        None,
        "chosen from the number of databases: more than one travels on stdin",
    ),
    dev(
        "proteins read",
        "param.accessions-stdin",
        Some("accessions"),
        "the accession list itself, an Option; Some sets the flag and fills stdin after the paths",
    ),
    dev(
        "proteins classify-peptides",
        "param.on-error",
        None,
        "the wire accepts only fail, the default, so there is nothing to choose",
    ),
];

/// Deviations: pride.
const PRIDE: &[Deviation] = &[
    // The listing functions return the entries: a Vec is what a caller iterates, and the envelope
    // around it is either the caller's own argument or derivable from the list.
    dev(
        "pride files",
        "field.accession",
        None,
        "the caller's own argument; each PrideFile also carries it as project_accession",
    ),
    dev(
        "pride files",
        "field.file_count",
        None,
        "list_files returns Vec<PrideFile>; this is its len()",
    ),
    dev(
        "pride files",
        "field.total_size_bytes",
        None,
        "pride::total_size_bytes(&files) computes it, and says why it is not a transfer size",
    ),
    dev(
        "pride files",
        "field.files",
        None,
        "list_files returns this list itself, as PrideFile values",
    ),
    dev(
        "pride ftp-files",
        "field.accession",
        None,
        "the caller's own argument",
    ),
    dev(
        "pride ftp-files",
        "field.file_count",
        None,
        "list_ftp_files returns Vec<PrideFtpFile>; this is its len()",
    ),
    dev(
        "pride ftp-files",
        "field.approximate_total_size_bytes",
        None,
        "pride::approximate_total_size_bytes(&files) computes it",
    ),
    dev(
        "pride ftp-files",
        "field.files",
        None,
        "list_ftp_files returns this list itself, as PrideFtpFile values",
    ),
    dev(
        "pride download",
        "param.dest",
        Some("destination"),
        "a Rust argument, spelled out; it takes any AsRef<Path>",
    ),
    dev(
        "pride download",
        "param.ext",
        Some("extensions"),
        "DownloadOptions::extensions, spelled out, because it takes several",
    ),
    dev(
        "pride download",
        "param.no-overwrite",
        Some("overwrite"),
        "stated positively: DownloadOptions::overwrite = Some(false) sends --no-overwrite",
    ),
    dev(
        "pride download",
        "param.names-from-stdin",
        None,
        "download_files is the selection form: passing it PrideFile values sets this flag",
    ),
    dev(
        "pride download",
        "param.stdin",
        None,
        "the file names of download_files' `files` argument, one per line",
    ),
    dev(
        "pride download",
        "field.accession",
        None,
        "download returns the written paths; the accession is the caller's own argument",
    ),
    dev(
        "pride download",
        "field.destination_directory",
        None,
        "the caller's own destination argument",
    ),
    dev(
        "pride download",
        "field.downloaded_count",
        None,
        "download returns Vec<PathBuf>; this is its len()",
    ),
    dev(
        "pride download",
        "field.paths",
        None,
        "download returns this list itself, as PathBuf values",
    ),
    dev(
        "pride search",
        "field.keyword",
        None,
        "search returns the hits themselves; the keyword is the caller's own argument",
    ),
    dev(
        "pride search",
        "field.result_count",
        None,
        "search returns Vec<PrideProjectSearchResult>; this is its len()",
    ),
    dev(
        "pride search",
        "field.results",
        None,
        "search returns this list itself, as PrideProjectSearchResult values",
    ),
];

/// Deviations: peptidoform.
const PEPTIDOFORM: &[Deviation] = &[
    dev(
        "peptidoform fragments",
        "param.no-modifications",
        Some("modifications"),
        "stated positively: FragmentOptions::modifications = false sends --no-modifications",
    ),
    dev(
        "peptidoform fragments",
        "param.max-mods",
        Some("max_modifications"),
        "mzLib's own name, as the result field that echoes it",
    ),
    dev(
        "peptidoform fragments",
        "field.annotated_modification_sites",
        Some("sites"),
        "grouped with its three siblings as Digest::modification_census, whose explain() reads them",
    ),
    dev(
        "peptidoform fragments",
        "field.annotated_modifications_loaded",
        Some("applied"),
        "ModificationCensus::applied: the modifications actually placed on the protein",
    ),
    dev(
        "peptidoform fragments",
        "field.uniprot_annotated_features",
        Some("annotated"),
        "ModificationCensus::annotated",
    ),
    dev(
        "peptidoform fragments",
        "field.unresolved_modifications",
        Some("unresolved"),
        "ModificationCensus::unresolved",
    ),
    dev(
        "peptidoform fragments",
        "field.uniprot_features_by_type",
        Some("by_type"),
        "ModificationCensus::by_type",
    ),
    dev(
        "peptidoform fragments",
        "field.max_modification_isoforms",
        Some("max_isoforms"),
        "the same name as the FragmentOptions field it echoes",
    ),
    dev(
        "peptidoform fragments",
        "field.peptides_at_isoform_cap",
        Some("peptides_at_cap"),
        "short, beside Digest::truncated(), which reads it",
    ),
    dev(
        "peptidoform fragments",
        "field.peptide_count",
        None,
        "Digest::peptides is a Vec; this is its len()",
    ),
    dev(
        "peptidoform fragments",
        "field.modification_count",
        None,
        "Peptide::modifications is a Vec; this is its len()",
    ),
];

/// Deviations: quant (flashlfq).
const QUANT: &[Deviation] = &[
    // QuantifyOptions uses FlashLFQ's own parameter names (the crate's "names follow mzLib"
    // convention), which are also the names the result's `parameters` echoes.
    dev(
        "quant flashlfq",
        "param.stdin",
        Some("spectra"),
        "the runs are quantify_with's `spectra` argument, one stdin line per SpectraFile",
    ),
    dev(
        "quant flashlfq",
        "param.ppm",
        Some("ppm_tolerance"),
        "FlashLFQ's name, echoed by parameters.ppm_tolerance",
    ),
    dev(
        "quant flashlfq",
        "param.isotope-ppm",
        Some("isotope_ppm_tolerance"),
        "FlashLFQ's name",
    ),
    dev(
        "quant flashlfq",
        "param.mbr",
        Some("match_between_runs"),
        "FlashLFQ's name",
    ),
    dev(
        "quant flashlfq",
        "param.mbr-ppm",
        Some("mbr_ppm_tolerance"),
        "FlashLFQ's name",
    ),
    dev(
        "quant flashlfq",
        "param.mbr-q",
        Some("mbr_q_value_threshold"),
        "FlashLFQ's name",
    ),
    dev(
        "quant flashlfq",
        "param.shared-peptides",
        Some("use_shared_peptides_for_protein_quant"),
        "FlashLFQ's name",
    ),
    dev(
        "quant flashlfq",
        "param.bayesian",
        Some("bayesian_protein_quant"),
        "FlashLFQ's name",
    ),
    dev(
        "quant flashlfq",
        "param.use-pep-q",
        Some("use_pep_q_value"),
        "FlashLFQ's name",
    ),
    dev(
        "quant flashlfq",
        "param.threads",
        Some("max_threads"),
        "FlashLFQ's MaxThreads",
    ),
    dev(
        "quant flashlfq",
        "param.out",
        Some("output_directory"),
        "it names a directory, and FlashLFQ writes several files into it",
    ),
    dev(
        "quant flashlfq",
        "field.peptide_count",
        None,
        "FlashLfqResults::peptides is a Vec; this is its len()",
    ),
    dev(
        "quant flashlfq",
        "field.protein_count",
        None,
        "FlashLfqResults::proteins is a Vec; this is its len()",
    ),
    dev(
        "sdrf pool",
        "param.stdin",
        Some("documents"),
        "the stdin lines are rendered from the PoolInput argument, one path[TAB label] per line",
    ),
];

/// Spec facts this crate does not project **yet**: `(verb, key)`. Each is a gap, not a choice,
/// and the lint skips it only until it is closed. The test fails on an entry that is no longer a
/// gap, so this list can only shrink.
const PENDING: &[(&str, &str)] = &[];

// =================================================================================================
// The vendored specs
// =================================================================================================

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn specs_dir() -> PathBuf {
    root().join("docs").join("specs")
}

fn reference_dir() -> PathBuf {
    root().join("docs").join("reference")
}

struct Spec {
    /// `readers.read-spectra.yaml`.
    file: String,
    body: Value,
}

impl Spec {
    /// `readers.read-spectra`.
    fn stem(&self) -> &str {
        self.file.trim_end_matches(".yaml")
    }

    fn verb(&self) -> &str {
        self.body["verb"].as_str().unwrap_or("")
    }

    fn list(&self, pointer: &str) -> Vec<Value> {
        self.body
            .pointer(pointer)
            .and_then(Value::as_array)
            .cloned()
            .unwrap_or_default()
    }

    fn rust(&self, key: &str) -> Option<String> {
        self.body
            .pointer(&format!("/bindings/rust/{key}"))
            .and_then(Value::as_str)
            .map(str::to_owned)
    }

    fn since(&self, binding: &str) -> Option<String> {
        match self.body.pointer(&format!("/since/{binding}")) {
            Some(Value::String(text)) => Some(text.clone()),
            Some(Value::Number(number)) => Some(number.to_string()),
            _ => None,
        }
    }
}

fn load_specs() -> Vec<Spec> {
    let mut paths: Vec<PathBuf> = std::fs::read_dir(specs_dir())
        .expect("docs/specs/ is missing; run scripts/sync-specs.ps1")
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "yaml"))
        .collect();
    paths.sort();
    paths
        .into_iter()
        .map(|path| {
            let text = std::fs::read_to_string(&path).unwrap();
            let body: Value = serde_yaml::from_str(&text)
                .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
            Spec {
                file: path.file_name().unwrap().to_string_lossy().into_owned(),
                body,
            }
        })
        .collect()
}

fn source_commit() -> String {
    std::fs::read_to_string(specs_dir().join("SOURCE"))
        .unwrap_or_default()
        .lines()
        .find_map(|line| line.strip_prefix("commit:"))
        .map_or_else(|| "unknown".to_owned(), |c| c.trim().to_owned())
}

fn text(value: &Value) -> Option<String> {
    match value {
        Value::Null => None,
        Value::String(text) => Some(text.clone()),
        Value::Bool(flag) => Some(flag.to_string()),
        Value::Number(number) => Some(number.to_string()),
        other => Some(other.to_string()),
    }
}

fn field(item: &Value, key: &str) -> Option<String> {
    item.get(key).and_then(text)
}

fn deviation(verb: &str, key: &str) -> Option<&'static Deviation> {
    let all = || DEVIATIONS.iter().flat_map(|module| module.iter());
    all()
        .find(|d| d.verb == verb && d.key == key)
        .or_else(|| all().find(|d| d.verb == "*" && d.key == key))
}

/// The Rust spelling of a spec param or field: its declared deviation, or the wire name in
/// `snake_case`. `None` when the deviation says it is not a named item.
fn rust_name(verb: &str, key: &str, wire: &str) -> Option<String> {
    match deviation(verb, key) {
        Some(declared) => declared.rust.map(str::to_owned),
        None => Some(wire.replace('-', "_")),
    }
}

// =================================================================================================
// Rendering
// =================================================================================================

const MZLIB_BLOB: &str = "https://github.com/smith-chem-wisc/mzLib/blob";
const FIXTURE_BLOB: &str = "https://github.com/smith-chem-wisc/mzLibRust/blob/main/tests/fixtures";
const DASH: &str = "—";

/// Prose from a spec, made safe for a Markdown table cell that rustdoc renders.
///
/// Rustdoc reads `[age]` in `characteristics[age]` as an intra-doc link and `<em>` as an HTML tag,
/// and underscores in `_ms1.feature` as emphasis, so every Markdown-active character is escaped.
fn prose(value: Option<String>) -> String {
    let Some(value) = value.filter(|v| !v.is_empty()) else {
        return DASH.to_owned();
    };
    let joined = value.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut out = String::with_capacity(joined.len());
    for c in joined.chars() {
        match c {
            '\\' | '*' | '_' | '[' | ']' | '|' | '`' | '#' => {
                out.push('\\');
                out.push(c);
            }
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '&' => out.push_str("&amp;"),
            c => out.push(c),
        }
    }
    out
}

/// A value shown as code. Inside a table only `|` needs escaping.
fn code(value: Option<String>) -> String {
    match value.filter(|v| !v.is_empty()) {
        None => DASH.to_owned(),
        Some(value) => format!("`{}`", value.replace('|', "\\|").replace('`', "'")),
    }
}

fn table(header: &[&str], rows: &[Vec<String>]) -> Vec<String> {
    let mut lines = vec![
        format!("| {} |", header.join(" | ")),
        format!("|{}|", vec!["---"; header.len()].join("|")),
    ];
    lines.extend(rows.iter().map(|row| format!("| {} |", row.join(" | "))));
    lines
}

fn param_rows(spec: &Spec, params: &[Value], prefix: &str) -> Vec<Vec<String>> {
    params
        .iter()
        .map(|p| {
            let wire = field(p, "name").unwrap_or_default();
            let shown = match rust_name(spec.verb(), &format!("{prefix}param.{wire}"), &wire) {
                Some(name) => format!("{} (wire `--{wire}`)", code(Some(name))),
                None => format!("wire `--{wire}`"),
            };
            let default = if p.get("required") == Some(&Value::Bool(true)) {
                "required".to_owned()
            } else {
                match p.get("default").and_then(text) {
                    None => "absent".to_owned(),
                    Some(value) => code(Some(value)),
                }
            };
            vec![
                shown,
                code(field(p, "type")),
                default,
                prose(field(p, "unit")),
                code(field(p, "range")),
                prose(field(p, "doc")),
            ]
        })
        .collect()
}

fn field_rows(spec: &Spec, fields: &[Value], key_prefix: &str) -> Vec<Vec<String>> {
    fields
        .iter()
        .map(|f| {
            let wire = field(f, "wire").unwrap_or_default();
            let shown = match rust_name(spec.verb(), &format!("{key_prefix}{wire}"), &wire) {
                Some(name) if name == wire => code(Some(name)),
                Some(name) => format!("{} (wire `{wire}`)", code(Some(name))),
                None => format!("wire `{wire}`"),
            };
            let null = if f.get("nullable") == Some(&Value::Bool(true)) {
                format!("yes: {}", prose(field(f, "null_means")))
            } else {
                "never".to_owned()
            };
            let mut doc = prose(field(f, "doc"));
            if let Some(when) = field(f, "present_when") {
                doc.push_str(&format!(" *Present only with `{when}`.*"));
            }
            vec![
                shown,
                code(field(f, "type")),
                prose(field(f, "unit")),
                null,
                doc,
            ]
        })
        .collect()
}

const FIELD_HEADER: &[&str] = &["Field", "Type", "Unit", "Null?", "Meaning"];
const PARAM_HEADER: &[&str] = &["Parameter", "Type", "Default", "Unit", "Range", "Meaning"];

/// `<stem>.md`: summary, wraps, parameters, returns, errors, caveats - the facts a caller needs
/// before the example, in the shared reference-page order.
fn render_facts(spec: &Spec, commit: &str) -> String {
    let verb = spec.verb();
    let mut lines = vec![
        generated_banner(spec, commit),
        String::new(),
        format!(
            "**Wire verb `{verb}`**: {}",
            prose(spec.body.get("summary").and_then(text))
        ),
        String::new(),
        "# Wraps".to_owned(),
        String::new(),
    ];
    for w in spec.list("/wraps") {
        let symbol = field(&w, "symbol").unwrap_or_default();
        let path = field(&w, "path").unwrap_or_default();
        let pin = field(&w, "pin").unwrap_or_default();
        let short: String = pin.chars().take(8).collect();
        lines.push(format!(
            "- [`{symbol}`]({MZLIB_BLOB}/{pin}/{path}) in `{path}` at mzLib `{short}`"
        ));
    }

    lines.extend(["".to_owned(), "# Parameters".to_owned(), String::new()]);
    let params = spec.list("/params");
    if params.is_empty() {
        lines.push("This verb takes no parameters.".to_owned());
    } else {
        lines.extend(table(PARAM_HEADER, &param_rows(spec, &params, "")));
    }

    lines.extend(["".to_owned(), "# Returns".to_owned(), String::new()]);
    let envelope = spec.list("/result/envelope_fields");
    if !envelope.is_empty() {
        lines.push("Top-level fields:".to_owned());
        lines.push(String::new());
        lines.extend(table(FIELD_HEADER, &field_rows(spec, &envelope, "field.")));
    }
    match spec.body.pointer("/result/columns") {
        Some(Value::Array(columns)) if !columns.is_empty() => {
            lines.extend([
                String::new(),
                "Per-row fields (one value per record, in `columns` or each list entry):"
                    .to_owned(),
                String::new(),
            ]);
            lines.extend(table(FIELD_HEADER, &field_rows(spec, columns, "field.")));
        }
        Some(Value::String(kind)) if kind == "per-format" => {
            lines.extend([
                String::new(),
                "**The columns are per-format**: they are the file's own record fields, listed in \
                 `column_names` on every result. Every cell follows these rules:"
                    .to_owned(),
                String::new(),
            ]);
            for rule in spec.list("/result/cell_rules") {
                lines.push(format!("- {}", prose(text(&rule))));
            }
        }
        _ => {}
    }
    if let Some(Value::Object(tables)) = spec.body.pointer("/result/tables") {
        for (name, fields) in tables {
            let fields = fields.as_array().cloned().unwrap_or_default();
            lines.extend([
                String::new(),
                format!("Fields of each entry of `{name}`:"),
                String::new(),
            ]);
            lines.extend(table(
                FIELD_HEADER,
                &field_rows(spec, &fields, &format!("table.{name}.")),
            ));
        }
    }

    lines.extend(["".to_owned(), "# Errors".to_owned(), String::new()]);
    let errors: Vec<Vec<String>> = spec
        .list("/errors")
        .iter()
        .map(|e| {
            let kind = field(e, "kind").unwrap_or_default();
            let variant = match kind.as_str() {
                "usage" => "[`MzLibError::Usage`](crate::MzLibError::Usage)",
                "service_unavailable" => {
                    "[`MzLibError::ServiceUnavailable`](crate::MzLibError::ServiceUnavailable)"
                }
                "correctness" => "[`MzLibError::Bridge`](crate::MzLibError::Bridge)",
                _ => DASH,
            };
            let when = match field(e, "when").as_deref() {
                Some("never") => "never returned by this verb".to_owned(),
                other => prose(other.map(str::to_owned)),
            };
            vec![code(Some(kind)), variant.to_owned(), when]
        })
        .collect();
    lines.extend(table(&["Kind", "mzLibRust returns", "When"], &errors));
    lines.push(String::new());
    lines.push(
        concat!(
            "Any call can also fail on the way to mzLib or back, whatever the verb: ",
            "[`MzLibError::BridgeNotFound`](crate::MzLibError::BridgeNotFound), ",
            "[`MzLibError::Timeout`](crate::MzLibError::Timeout), ",
            "[`MzLibError::Protocol`](crate::MzLibError::Protocol) or ",
            "[`MzLibError::Io`](crate::MzLibError::Io)."
        )
        .to_owned(),
    );

    lines.extend(["".to_owned(), "# Caveats".to_owned(), String::new()]);
    let caveats = spec.list("/caveats");
    if caveats.is_empty() {
        lines.push("None recorded.".to_owned());
    }
    for caveat in caveats {
        lines.push(format!("- {}", prose(text(&caveat))));
    }
    lines.join("\n") + "\n"
}

/// `<stem>.bulk.md`: the many-input form's parameters and envelope (BULK.md).
fn render_bulk(spec: &Spec, commit: &str) -> Option<String> {
    let bulk = spec.body.pointer("/result/bulk")?;
    let mut lines = vec![
        generated_banner(spec, commit),
        String::new(),
        format!(
            "**Wire verb `{}` with `--paths-stdin`**: many inputs, one bridge process, one long \
             table, per BULK.md.",
            spec.verb()
        ),
    ];
    let params = bulk
        .get("params")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    if !params.is_empty() {
        lines.extend(["".to_owned(), "# Parameters".to_owned(), String::new()]);
        lines.extend(table(PARAM_HEADER, &param_rows(spec, &params, "bulk.")));
    }
    let envelope = bulk
        .get("envelope_fields")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    if !envelope.is_empty() {
        lines.extend(["".to_owned(), "# Returns".to_owned(), String::new()]);
        lines.extend(table(
            FIELD_HEADER,
            &field_rows(spec, &envelope, "bulk.field."),
        ));
    }
    Some(lines.join("\n") + "\n")
}

/// `<stem>.see-also.md`: performance, the same verb elsewhere, references, since, and what the
/// spec has not verified - the facts that close a reference page.
fn render_tail(spec: &Spec, commit: &str, shipped_in_rust: bool) -> String {
    let mut lines = vec![generated_banner(spec, commit)];

    if let Some(performance) = spec.body.get("performance").and_then(text) {
        lines.extend(["".to_owned(), "# Performance".to_owned(), String::new()]);
        lines.push(prose(Some(performance)));
    }

    let examples = spec.list("/examples");
    if !examples.is_empty() {
        lines.extend([
            "".to_owned(),
            "# Recorded examples".to_owned(),
            String::new(),
        ]);
        lines.push(
            "The examples on this page run in CI against these recordings from the live bridge, \
             replayed by a stand-in bridge that answers only a call the recording fits:"
                .to_owned(),
        );
        lines.push(String::new());
        for example in examples {
            let fixture = field(&example, "fixture").unwrap_or_default();
            lines.push(format!("- [`{fixture}`]({FIXTURE_BLOB}/{fixture})"));
        }
    }

    let spelling = |lang: &str| -> String {
        let name = spec
            .body
            .pointer(&format!("/bindings/{lang}/name"))
            .and_then(Value::as_str);
        let Some(name) = name else {
            return DASH.to_owned();
        };
        let name = if lang == "rust" {
            deviation(spec.verb(), "name")
                .and_then(|d| d.rust)
                .unwrap_or(name)
        } else {
            name
        };
        let mut shown = code(Some(name.to_owned()));
        if let Some(options) = spec
            .body
            .pointer(&format!("/bindings/{lang}/options"))
            .and_then(Value::as_str)
        {
            shown.push_str(&format!(" with `{options}`"));
        }
        if let Some(bulk) = spec
            .body
            .pointer(&format!("/bindings/{lang}/bulk"))
            .and_then(Value::as_str)
        {
            let bulk = if lang == "rust" {
                deviation(spec.verb(), "bulk")
                    .and_then(|d| d.rust)
                    .unwrap_or(bulk)
            } else {
                bulk
            };
            shown.push_str(&format!("; many inputs: `{bulk}`"));
        }
        shown
    };
    lines.extend([
        "".to_owned(),
        "# Same verb in other bindings".to_owned(),
        String::new(),
    ]);
    lines.extend(table(
        &["Python (pyMzLib)", "Rust (mzLibRust)", "R (mzLibR)"],
        &[vec![spelling("py"), spelling("rust"), spelling("r")]],
    ));

    lines.extend(["".to_owned(), "# References".to_owned(), String::new()]);
    let cites = spec.list("/cite");
    if cites.is_empty() {
        lines.push("None: no publication is attached to this verb.".to_owned());
    }
    for cite in cites {
        let doi = field(&cite, "doi").unwrap_or_default();
        lines.push(format!(
            "- [doi:{doi}](https://doi.org/{doi}): {}",
            prose(field(&cite, "for"))
        ));
    }

    lines.extend(["".to_owned(), "# Since".to_owned(), String::new()]);
    let since = |binding: &str| match spec.since(binding) {
        Some(version) => version,
        None if binding == "mzlibrust" && shipped_in_rust => {
            "next release (the spec has not recorded it yet)".to_owned()
        }
        None => "not yet shipped".to_owned(),
    };
    lines.extend(table(
        &["Wire protocol", "pyMzLib", "mzLibRust", "mzLibR"],
        &[vec![
            since("protocol"),
            since("pymzlib"),
            since("mzlibrust"),
            since("mzlibr"),
        ]],
    ));

    let questions = spec.list("/open_questions");
    if !questions.is_empty() {
        lines.extend([
            "".to_owned(),
            "# Not yet verified".to_owned(),
            String::new(),
        ]);
        lines.push("The spec lists these as open; they are reported, never hidden:".to_owned());
        lines.push(String::new());
        for question in questions {
            lines.push(format!("- {}", prose(text(&question))));
        }
    }
    lines.join("\n") + "\n"
}

fn generated_banner(spec: &Spec, commit: &str) -> String {
    let short: String = commit.chars().take(12).collect();
    format!(
        "<!-- GENERATED by tests/spec_docs.rs from docs/specs/{} (bridge {short}). Do not edit: fix \
         the spec in the bridge, sync it, re-render. -->",
        spec.file
    )
}

fn expected_fragments(specs: &[Spec], rust: &RustSource) -> BTreeMap<PathBuf, String> {
    let commit = source_commit();
    let mut out = BTreeMap::new();
    for spec in specs {
        let dir = reference_dir();
        out.insert(
            dir.join(format!("{}.md", spec.stem())),
            render_facts(spec, &commit),
        );
        if let Some(bulk) = render_bulk(spec, &commit) {
            out.insert(dir.join(format!("{}.bulk.md", spec.stem())), bulk);
        }
        let shipped = binding_function(spec, rust).is_some();
        out.insert(
            dir.join(format!("{}.see-also.md", spec.stem())),
            render_tail(spec, &commit, shipped),
        );
    }
    out
}

// =================================================================================================
// This crate's own source, read with syn
// =================================================================================================

#[derive(Debug, Clone)]
struct FieldInfo {
    name: String,
    /// The wire name, from `#[serde(rename = "..")]`, when it differs.
    wire: Option<String>,
    doc: String,
    ty: syn::Type,
}

#[derive(Debug, Clone)]
struct StructInfo {
    fields: Vec<FieldInfo>,
}

#[derive(Debug, Clone)]
struct FnInfo {
    doc: String,
    /// `#[doc = include_str!("..")]` targets, relative to the source file.
    includes: Vec<String>,
    inputs: Vec<(String, syn::Type)>,
    output: Option<syn::Type>,
}

/// Every struct and free function, keyed by `(module, name)`; `lib` is the crate root.
struct RustSource {
    structs: BTreeMap<(String, String), StructInfo>,
    functions: BTreeMap<(String, String), FnInfo>,
}

fn doc_of(attrs: &[syn::Attribute]) -> (String, Vec<String>) {
    let mut doc = String::new();
    let mut includes = Vec::new();
    for attr in attrs {
        if !attr.path().is_ident("doc") {
            continue;
        }
        if let syn::Meta::NameValue(pair) = &attr.meta {
            match &pair.value {
                syn::Expr::Lit(syn::ExprLit {
                    lit: syn::Lit::Str(line),
                    ..
                }) => {
                    doc.push_str(&line.value());
                    doc.push('\n');
                }
                syn::Expr::Macro(mac) if mac.mac.path.is_ident("include_str") => {
                    if let Ok(target) = mac.mac.parse_body::<syn::LitStr>() {
                        includes.push(target.value());
                    }
                }
                _ => {}
            }
        }
    }
    (doc, includes)
}

fn serde_rename(attrs: &[syn::Attribute]) -> Option<String> {
    let mut rename = None;
    for attr in attrs.iter().filter(|a| a.path().is_ident("serde")) {
        let _ = attr.parse_nested_meta(|meta| {
            if meta.path.is_ident("rename") {
                let value: syn::LitStr = meta.value()?.parse()?;
                rename = Some(value.value());
            } else if meta.input.peek(syn::Token![=]) {
                let _: syn::Expr = meta.value()?.parse()?;
            }
            Ok(())
        });
    }
    rename
}

fn read_source() -> RustSource {
    let mut structs = BTreeMap::new();
    let mut functions = BTreeMap::new();
    let src = root().join("src");
    let mut files: Vec<PathBuf> = std::fs::read_dir(&src)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "rs"))
        .collect();
    files.sort();
    for path in files {
        let module = path.file_stem().unwrap().to_string_lossy().into_owned();
        let file = syn::parse_file(&std::fs::read_to_string(&path).unwrap())
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        for item in file.items {
            match item {
                syn::Item::Struct(item) if matches!(item.vis, syn::Visibility::Public(_)) => {
                    let fields = match &item.fields {
                        syn::Fields::Named(named) => named
                            .named
                            .iter()
                            .filter(|f| matches!(f.vis, syn::Visibility::Public(_)))
                            .map(|f| FieldInfo {
                                name: f.ident.as_ref().unwrap().to_string(),
                                wire: serde_rename(&f.attrs),
                                doc: doc_of(&f.attrs).0,
                                ty: f.ty.clone(),
                            })
                            .collect(),
                        _ => Vec::new(),
                    };
                    structs.insert(
                        (module.clone(), item.ident.to_string()),
                        StructInfo { fields },
                    );
                }
                syn::Item::Fn(item) if matches!(item.vis, syn::Visibility::Public(_)) => {
                    let (doc, includes) = doc_of(&item.attrs);
                    let inputs = item
                        .sig
                        .inputs
                        .iter()
                        .filter_map(|input| match input {
                            syn::FnArg::Typed(typed) => match &*typed.pat {
                                syn::Pat::Ident(ident) => {
                                    Some((ident.ident.to_string(), (*typed.ty).clone()))
                                }
                                _ => None,
                            },
                            syn::FnArg::Receiver(_) => None,
                        })
                        .collect();
                    let output = match &item.sig.output {
                        syn::ReturnType::Type(_, ty) => Some((**ty).clone()),
                        syn::ReturnType::Default => None,
                    };
                    functions.insert(
                        (module.clone(), item.sig.ident.to_string()),
                        FnInfo {
                            doc,
                            includes,
                            inputs,
                            output,
                        },
                    );
                }
                _ => {}
            }
        }
    }
    RustSource { structs, functions }
}

/// The last path segment of a type, through references: `&readers::ReadOptions` -> `ReadOptions`.
fn last_segment(ty: &syn::Type) -> Option<&syn::PathSegment> {
    match ty {
        syn::Type::Path(path) => path.path.segments.last(),
        syn::Type::Reference(reference) => last_segment(&reference.elem),
        syn::Type::Group(group) => last_segment(&group.elem),
        syn::Type::Paren(paren) => last_segment(&paren.elem),
        _ => None,
    }
}

fn generic_argument(segment: &syn::PathSegment) -> Option<&syn::Type> {
    match &segment.arguments {
        syn::PathArguments::AngleBracketed(args) => args.args.iter().find_map(|arg| match arg {
            syn::GenericArgument::Type(ty) => Some(ty),
            _ => None,
        }),
        _ => None,
    }
}

/// The struct a type names, through `Option`, `Vec`, `Box` and references.
fn base_struct(ty: &syn::Type) -> Option<String> {
    let segment = last_segment(ty)?;
    let name = segment.ident.to_string();
    match name.as_str() {
        "Option" | "Vec" | "Box" => base_struct(generic_argument(segment)?),
        _ => Some(name),
    }
}

impl RustSource {
    fn find_struct(&self, module: &str, name: &str) -> Option<(&str, &StructInfo)> {
        self.structs
            .get_key_value(&(module.to_owned(), name.to_owned()))
            .map(|((m, _), info)| (m.as_str(), info))
            .or_else(|| {
                self.structs
                    .iter()
                    .find(|((_, n), _)| n == name)
                    .map(|((m, _), info)| (m.as_str(), info))
            })
    }

    /// A struct's fields, and those of every struct-typed field inside it, depth first; `Table`
    /// is a leaf.
    fn fields_deep(&self, module: &str, name: &str) -> Vec<FieldInfo> {
        let mut seen = BTreeSet::new();
        let mut out = Vec::new();
        self.collect_fields(module, name, &mut seen, &mut out);
        out
    }

    fn collect_fields(
        &self,
        module: &str,
        name: &str,
        seen: &mut BTreeSet<String>,
        out: &mut Vec<FieldInfo>,
    ) {
        if name == "Table" || !seen.insert(name.to_owned()) {
            return;
        }
        let Some((module, info)) = self.find_struct(module, name) else {
            return;
        };
        let module = module.to_owned();
        for f in &info.fields {
            out.push(f.clone());
        }
        for f in &info.fields {
            // Only an options or result struct's own sub-structs, not a list element type: a
            // Vec<Format> field describes entries, which the columns check handles.
            if last_segment(&f.ty).is_some_and(|s| s.ident == "Vec") {
                continue;
            }
            if let Some(inner) = base_struct(&f.ty) {
                self.collect_fields(&module, &inner, seen, out);
            }
        }
    }
}

/// `mzlib::readers::read_spectra_with` -> `("readers", "read_spectra_with")`; the crate root is
/// found through its re-exports, so `mzlib::bridge_version` resolves to `bridge`.
fn split_path(path: &str) -> (String, String) {
    let parts: Vec<&str> = path.split("::").collect();
    match parts.as_slice() {
        [_, name] => ("lib".to_owned(), (*name).to_owned()),
        [_, module, name] => ((*module).to_owned(), (*name).to_owned()),
        _ => (String::new(), path.to_owned()),
    }
}

fn lookup_fn<'a>(rust: &'a RustSource, path: &str) -> Option<(String, &'a FnInfo)> {
    let (module, name) = split_path(path);
    if module == "lib" {
        return rust
            .functions
            .iter()
            .find(|((_, n), _)| *n == name)
            .map(|((m, _), info)| (m.clone(), info));
    }
    rust.functions
        .get(&(module.clone(), name))
        .map(|info| (module, info))
}

fn binding_function<'a>(spec: &Spec, rust: &'a RustSource) -> Option<(String, &'a FnInfo)> {
    let name = deviation(spec.verb(), "name")
        .and_then(|d| d.rust.map(str::to_owned))
        .or_else(|| spec.rust("name"))?;
    lookup_fn(rust, &name)
}

fn bulk_function<'a>(spec: &Spec, rust: &'a RustSource) -> Option<(String, &'a FnInfo)> {
    let name = deviation(spec.verb(), "bulk")
        .and_then(|d| d.rust.map(str::to_owned))
        .or_else(|| spec.rust("bulk"))?;
    lookup_fn(rust, &name)
}

/// What a function returns inside `Result<..>`, and whether it is a list.
fn returned(info: &FnInfo) -> Option<(String, bool)> {
    let output = info.output.as_ref()?;
    let segment = last_segment(output)?;
    let inner = if segment.ident == "Result" {
        generic_argument(segment)?
    } else {
        output
    };
    let inner_segment = last_segment(inner)?;
    if inner_segment.ident == "Vec" {
        Some((base_struct(generic_argument(inner_segment)?)?, true))
    } else {
        Some((inner_segment.ident.to_string(), false))
    }
}

// =================================================================================================
// Units
// =================================================================================================

fn unit_aliases(head: &str) -> &'static [&'static str] {
    match head {
        "min" => &["minute"],
        "s" => &["second"],
        "ms" => &["millisecond"],
        "da" => &["dalton"],
        "v" => &["volt"],
        _ => &[],
    }
}

/// Whether `text` names `unit`: the unit itself, its singular, or a spelled-out alias. The same
/// rule as pyMzLib's `test_spec_docs.unit_mentioned`, so the two bindings agree on what counts.
fn unit_mentioned(unit: &str, text: &str) -> bool {
    let text = text.to_lowercase();
    let unit = unit.trim().to_lowercase();
    let head = unit.split(" (").next().unwrap_or("").trim().to_owned();
    let mut candidates = vec![unit.clone(), head.clone()];
    if head.ends_with('s') && head.len() > 3 {
        candidates.push(head[..head.len() - 1].to_owned());
    }
    candidates.extend(unit_aliases(&head).iter().map(|a| (*a).to_owned()));

    let bytes = text.as_bytes();
    let is_letter = |i: usize| bytes.get(i).is_some_and(u8::is_ascii_lowercase);
    candidates
        .iter()
        .filter(|c| !c.is_empty())
        .any(|candidate| {
            let short_or_symbolic =
                candidate.len() <= 2 || !candidate.chars().all(char::is_alphabetic);
            text.match_indices(candidate.as_str()).any(|(at, _)| {
                let before_ok = at == 0 || !is_letter(at - 1);
                let after_ok = !short_or_symbolic || !is_letter(at + candidate.len());
                before_ok && after_ok
            })
        })
}

// =================================================================================================
// The tests
// =================================================================================================

#[test]
fn specs_are_vendored_with_their_source() {
    assert!(
        !load_specs().is_empty(),
        "docs/specs/ holds no specs; run scripts/sync-specs.ps1"
    );
    assert!(
        specs_dir().join("SOURCE").exists(),
        "docs/specs/SOURCE is missing; re-run scripts/sync-specs.ps1"
    );
}

/// The keys a param or a result field may carry (bridge design/verbs/README.md). Anything else is
/// nearly always YAML's flow-mapping trap: an unquoted comma inside `{... doc: a, b}` ends the doc
/// at the comma and turns the rest into a stray key, silently truncating what every binding shows.
const FIELD_KEYS: &[&str] = &[
    "name",
    "wire",
    "type",
    "required",
    "unit",
    "range",
    "default",
    "doc",
    "nullable",
    "null_means",
    "present_when",
];

#[test]
fn no_spec_field_was_truncated_by_yaml() {
    let mut found = Vec::new();
    for spec in load_specs() {
        let mut items = spec.list("/params");
        items.extend(spec.list("/result/envelope_fields"));
        items.extend(spec.list("/result/columns"));
        items.extend(spec.list("/result/bulk/params"));
        items.extend(spec.list("/result/bulk/envelope_fields"));
        for item in items {
            if let Some(map) = item.as_object() {
                let stray: Vec<&String> = map
                    .keys()
                    .filter(|k| !FIELD_KEYS.contains(&k.as_str()))
                    .collect();
                if !stray.is_empty() {
                    found.push(format!(
                        "{}: {} has stray keys {stray:?}",
                        spec.file,
                        field(&item, "wire")
                            .or_else(|| field(&item, "name"))
                            .unwrap_or_default()
                    ));
                }
            }
        }
    }
    assert!(
        found.is_empty(),
        "An unquoted comma inside a flow mapping truncated a doc. Quote it in the bridge's spec and \
         re-sync; never edit docs/specs/ by hand:\n  {}",
        found.join("\n  ")
    );
}

#[test]
fn rendered_fragments_are_current() {
    let specs = load_specs();
    let rust = read_source();
    let want = expected_fragments(&specs, &rust);
    let dir = reference_dir();
    let have: BTreeSet<PathBuf> = std::fs::read_dir(&dir)
        .map(|entries| {
            entries
                .map(|e| e.unwrap().path())
                .filter(|p| p.extension().is_some_and(|ext| ext == "md"))
                .collect()
        })
        .unwrap_or_default();

    if std::env::var_os("MZLIB_RENDER_SPEC_DOCS").is_some() {
        std::fs::create_dir_all(&dir).unwrap();
        for stale in have.iter().filter(|p| !want.contains_key(*p)) {
            std::fs::remove_file(stale).unwrap();
        }
        for (path, text) in &want {
            let current = std::fs::read_to_string(path).unwrap_or_default();
            if current.replace("\r\n", "\n") != *text {
                std::fs::write(path, text).unwrap();
            }
        }
        return;
    }

    let mut problems = Vec::new();
    for stale in have.iter().filter(|p| !want.contains_key(*p)) {
        problems.push(format!("{}: no spec renders it any more", stale.display()));
    }
    for (path, text) in &want {
        match std::fs::read_to_string(path) {
            Err(_) => problems.push(format!("{}: missing", path.display())),
            Ok(current) if current.replace("\r\n", "\n") != *text => {
                problems.push(format!("{}: differs from its spec", path.display()));
            }
            Ok(_) => {}
        }
    }
    assert!(
        problems.is_empty(),
        "The reference facts are stale. Regenerate them with\n    \
         MZLIB_RENDER_SPEC_DOCS=1 cargo test --test spec_docs\n  {}",
        problems.join("\n  ")
    );
}

#[test]
fn every_deviation_names_a_real_param_or_field() {
    let specs = load_specs();
    let by_verb: HashMap<&str, &Spec> = specs.iter().map(|s| (s.verb(), s)).collect();
    let mut problems = Vec::new();
    for d in DEVIATIONS.iter().flat_map(|module| module.iter()) {
        if d.why.trim().is_empty() {
            problems.push(format!("{} {}: no reason given", d.verb, d.key));
        }
        if d.verb == "*" {
            continue;
        }
        let Some(spec) = by_verb.get(d.verb) else {
            problems.push(format!("'{}' has no vendored spec", d.verb));
            continue;
        };
        if !spec_keys(spec).contains(d.key) {
            problems.push(format!("{} names nothing in {}", d.key, spec.file));
        }
    }
    for (verb, key) in PENDING {
        match by_verb.get(verb) {
            None => problems.push(format!("PENDING names '{verb}', which has no spec")),
            Some(spec) if !spec_keys(spec).contains(*key) => {
                problems.push(format!("PENDING {key} names nothing in {}", spec.file));
            }
            Some(_) => {}
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

fn spec_keys(spec: &Spec) -> BTreeSet<String> {
    let mut keys = BTreeSet::from(["name".to_owned(), "bulk".to_owned()]);
    for p in spec.list("/params") {
        keys.insert(format!("param.{}", field(&p, "name").unwrap_or_default()));
    }
    for f in spec.list("/result/envelope_fields") {
        keys.insert(format!("field.{}", field(&f, "wire").unwrap_or_default()));
    }
    for f in spec.list("/result/columns") {
        keys.insert(format!("field.{}", field(&f, "wire").unwrap_or_default()));
    }
    for p in spec.list("/result/bulk/params") {
        keys.insert(format!(
            "bulk.param.{}",
            field(&p, "name").unwrap_or_default()
        ));
    }
    for f in spec.list("/result/bulk/envelope_fields") {
        keys.insert(format!(
            "bulk.field.{}",
            field(&f, "wire").unwrap_or_default()
        ));
    }
    if let Some(Value::Object(tables)) = spec.body.pointer("/result/tables") {
        for (name, fields) in tables {
            for f in fields.as_array().cloned().unwrap_or_default() {
                keys.insert(format!(
                    "table.{name}.{}",
                    field(&f, "wire").unwrap_or_default()
                ));
            }
        }
    }
    keys
}

/// What the lint found for one spec, before [`PENDING`] is applied.
struct Finding {
    verb: String,
    key: String,
    problem: String,
}

#[test]
fn every_param_and_field_is_documented_with_its_unit() {
    let specs = load_specs();
    let rust = read_source();
    let mut findings: Vec<Finding> = Vec::new();
    let mut checked: BTreeSet<(String, String)> = BTreeSet::new();

    for spec in &specs {
        let verb = spec.verb().to_owned();
        let Some(rust_name) = spec.rust("name") else {
            continue;
        };
        let Some((module, function)) = binding_function(spec, &rust) else {
            if spec.since("mzlibrust").is_some() {
                findings.push(Finding {
                    verb: verb.clone(),
                    key: "name".to_owned(),
                    problem: format!(
                        "{} says mzLibRust ships {rust_name} (since {}), but it does not exist",
                        spec.file,
                        spec.since("mzlibrust").unwrap()
                    ),
                });
            }
            continue;
        };

        // ---- parameters: the options struct, or the function's own arguments ------------------
        let options = spec.rust("options");
        let option_fields = options
            .as_deref()
            .map(|o| rust.fields_deep(&module, o))
            .unwrap_or_default();
        for p in spec.list("/params") {
            let key = format!("param.{}", field(&p, "name").unwrap_or_default());
            let owner = options.as_deref().unwrap_or("the options struct");
            check_param(
                &verb,
                key,
                &p,
                &option_fields,
                owner,
                function,
                &mut findings,
                &mut checked,
            );
        }

        // ---- the bulk form ---------------------------------------------------------------------
        let bulk_params = spec.list("/result/bulk/params");
        let bulk_fields = spec.list("/result/bulk/envelope_fields");
        if spec.rust("bulk").is_some() {
            match bulk_function(spec, &rust) {
                None => {
                    for p in &bulk_params {
                        let key = format!("bulk.param.{}", field(p, "name").unwrap_or_default());
                        findings.push(Finding {
                            verb: verb.clone(),
                            key,
                            problem: format!("{} does not exist", spec.rust("bulk").unwrap()),
                        });
                    }
                    for f in &bulk_fields {
                        let key = format!("bulk.field.{}", field(f, "wire").unwrap_or_default());
                        findings.push(Finding {
                            verb: verb.clone(),
                            key,
                            problem: format!("{} does not exist", spec.rust("bulk").unwrap()),
                        });
                    }
                }
                Some((bulk_module, bulk)) => {
                    // The options struct is the argument whose type is a struct named *Options.
                    let bulk_options = bulk
                        .inputs
                        .iter()
                        .filter_map(|(_, ty)| base_struct(ty))
                        .find(|name| name.ends_with("Options"));
                    let fields = bulk_options
                        .as_deref()
                        .map(|o| rust.fields_deep(&bulk_module, o))
                        .unwrap_or_default();
                    for p in &bulk_params {
                        let key = format!("bulk.param.{}", field(p, "name").unwrap_or_default());
                        let owner = bulk_options.as_deref().unwrap_or("the bulk options struct");
                        check_param(
                            &verb,
                            key,
                            p,
                            &fields,
                            owner,
                            bulk,
                            &mut findings,
                            &mut checked,
                        );
                    }
                    if let Some((returned, _)) = returned(bulk) {
                        let fields = rust.fields_deep(&bulk_module, &returned);
                        for f in &bulk_fields {
                            let key =
                                format!("bulk.field.{}", field(f, "wire").unwrap_or_default());
                            check_field(
                                &verb,
                                key,
                                f,
                                &fields,
                                &returned,
                                &mut findings,
                                &mut checked,
                            );
                        }
                    }
                }
            }
        }

        // ---- the result --------------------------------------------------------------------------
        let Some((result, is_list)) = returned(function) else {
            continue;
        };
        let envelope = spec.list("/result/envelope_fields");
        let columns = spec.list("/result/columns");
        if is_list {
            // A function that returns a list returns the entries; the envelope has nowhere to be
            // documented, so each of its fields must be a declared deviation.
            for f in &envelope {
                let key = format!("field.{}", field(f, "wire").unwrap_or_default());
                if deviation(&verb, &key).is_none() {
                    findings.push(Finding {
                        verb: verb.clone(),
                        key,
                        problem: format!(
                            "{rust_name} returns a list, so this envelope field needs a DEVIATIONS \
                             entry"
                        ),
                    });
                }
            }
            let fields = rust.fields_deep(&module, &result);
            for c in &columns {
                let key = format!("field.{}", field(c, "wire").unwrap_or_default());
                check_field(&verb, key, c, &fields, &result, &mut findings, &mut checked);
            }
        } else {
            let fields = rust.fields_deep(&module, &result);
            for f in &envelope {
                let key = format!("field.{}", field(f, "wire").unwrap_or_default());
                check_field(&verb, key, f, &fields, &result, &mut findings, &mut checked);
            }
            // Record-list results (`results`, `proteins`): the columns describe the entries of
            // the envelope field whose spec doc says its fields are "under result.columns". A
            // result that holds a columnar Table carries its columns as data, not as fields.
            if !columns.is_empty() {
                let container = envelope.iter().find(|f| {
                    field(f, "doc").is_some_and(|doc| doc.contains("under result.columns"))
                });
                let entry = container
                    .and_then(|f| {
                        let wire = field(f, "wire").unwrap_or_default();
                        self::rust_name(&verb, &format!("field.{wire}"), &wire)
                    })
                    .and_then(|name| {
                        rust.find_struct(&module, &result).and_then(|(_, info)| {
                            info.fields
                                .iter()
                                .find(|f| f.name == name)
                                .and_then(|f| base_struct(&f.ty))
                        })
                    });
                if let Some(entry) = entry {
                    let fields = rust.fields_deep(&module, &entry);
                    for c in &columns {
                        let key = format!("field.{}", field(c, "wire").unwrap_or_default());
                        check_field(&verb, key, c, &fields, &entry, &mut findings, &mut checked);
                    }
                }
            }
            // Typed sub-tables: `tables.<name>` describes the entries of the field `<name>`.
            for (table_name, table_fields) in spec_tables(spec) {
                let element = rust.find_struct(&module, &result).and_then(|(_, info)| {
                    info.fields
                        .iter()
                        .find(|f| f.name == table_name)
                        .and_then(|f| base_struct(&f.ty))
                });
                let Some(element) = element else { continue };
                let fields = rust.fields_deep(&module, &element);
                if fields.is_empty() || holds_table(&fields) {
                    continue; // a columnar Table: its columns are data, not struct fields
                }
                for f in &table_fields {
                    let key = format!(
                        "table.{table_name}.{}",
                        field(f, "wire").unwrap_or_default()
                    );
                    check_field(
                        &verb,
                        key,
                        f,
                        &fields,
                        &element,
                        &mut findings,
                        &mut checked,
                    );
                }
            }
        }
    }

    // ---- apply PENDING, and make sure it can only shrink --------------------------------------
    let pending: BTreeSet<(String, String)> = PENDING
        .iter()
        .map(|(v, k)| ((*v).to_owned(), (*k).to_owned()))
        .collect();
    let mut problems: Vec<String> = findings
        .iter()
        .filter(|f| !pending.contains(&(f.verb.clone(), f.key.clone())))
        .map(|f| format!("{} {}: {}", f.verb, f.key, f.problem))
        .collect();
    let failing: BTreeSet<(String, String)> = findings
        .iter()
        .map(|f| (f.verb.clone(), f.key.clone()))
        .collect();
    for entry in &pending {
        if !failing.contains(entry) && checked.contains(entry) {
            problems.push(format!(
                "{} {}: projected and documented now; delete it from PENDING",
                entry.0, entry.1
            ));
        }
    }
    problems.sort();
    assert!(
        problems.is_empty(),
        "{} problem(s) against the specs in docs/specs/:\n  {}",
        problems.len(),
        problems.join("\n  ")
    );
}

/// Whether a struct's fields include a columnar [`Table`]-typed one (flattened or not).
fn holds_table(fields: &[FieldInfo]) -> bool {
    fields
        .iter()
        .any(|f| base_struct(&f.ty).as_deref() == Some("Table"))
}

fn spec_tables(spec: &Spec) -> Vec<(String, Vec<Value>)> {
    match spec.body.pointer("/result/tables") {
        Some(Value::Object(tables)) => tables
            .iter()
            .map(|(name, fields)| (name.clone(), fields.as_array().cloned().unwrap_or_default()))
            .collect(),
        _ => Vec::new(),
    }
}

#[allow(clippy::too_many_arguments)]
fn check_param(
    verb: &str,
    key: String,
    param: &Value,
    fields: &[FieldInfo],
    owner: &str,
    function: &FnInfo,
    findings: &mut Vec<Finding>,
    checked: &mut BTreeSet<(String, String)>,
) {
    let wire = field(param, "name").unwrap_or_default();
    let Some(name) = rust_name(verb, &key, &wire) else {
        return;
    };
    checked.insert((verb.to_owned(), key.clone()));
    let unit = field(param, "unit");
    let mut report = |problem: String| {
        findings.push(Finding {
            verb: verb.to_owned(),
            key: key.clone(),
            problem,
        })
    };
    if let Some(found) = fields.iter().find(|x| x.name == name) {
        match unit.as_deref() {
            Some(unit) if !unit_mentioned(unit, &found.doc) => report(format!(
                "{owner}::{name} is documented without its unit '{unit}': {:?}",
                found.doc.trim()
            )),
            None if found.doc.trim().is_empty() => report(format!("{owner}::{name} has no doc")),
            _ => {}
        }
    } else if function.inputs.iter().any(|(n, _)| *n == name) {
        if let Some(unit) = unit.as_deref() {
            if !unit_mentioned(unit, &function.doc) {
                report(format!(
                    "argument `{name}` is documented without its unit '{unit}'"
                ));
            }
        }
    } else {
        report(format!(
            "`{name}` (wire --{wire}) is neither a field of {owner} nor an argument"
        ));
    }
}

fn check_field(
    verb: &str,
    key: String,
    spec_field: &Value,
    fields: &[FieldInfo],
    owner: &str,
    findings: &mut Vec<Finding>,
    checked: &mut BTreeSet<(String, String)>,
) {
    let wire = field(spec_field, "wire").unwrap_or_default();
    let Some(name) = rust_name(verb, &key, &wire) else {
        return;
    };
    checked.insert((verb.to_owned(), key.clone()));
    let found = fields.iter().find(|f| {
        f.name == name
            || (f.wire.as_deref() == Some(wire.as_str()) && deviation(verb, &key).is_none())
    });
    let Some(found) = found else {
        findings.push(Finding {
            verb: verb.to_owned(),
            key,
            problem: format!("`{name}` is not a field of {owner}"),
        });
        return;
    };
    if found.doc.trim().is_empty() {
        findings.push(Finding {
            verb: verb.to_owned(),
            key,
            problem: format!("{owner}::{} has no doc", found.name),
        });
        return;
    }
    if let Some(unit) = field(spec_field, "unit") {
        if !unit_mentioned(&unit, &found.doc) {
            findings.push(Finding {
                verb: verb.to_owned(),
                key,
                problem: format!(
                    "{owner}::{} is documented without its unit '{unit}': {:?}",
                    found.name,
                    found.doc.trim()
                ),
            });
        }
    }
}

#[test]
fn every_shipped_verb_includes_its_reference_facts() {
    let specs = load_specs();
    let rust = read_source();
    let mut problems = Vec::new();
    for spec in &specs {
        let Some((_, function)) = binding_function(spec, &rust) else {
            continue;
        };
        for wanted in [
            format!("{}.md", spec.stem()),
            format!("{}.see-also.md", spec.stem()),
        ] {
            if !function.includes.iter().any(|i| i.ends_with(&wanted)) {
                problems.push(format!(
                    "{}: add #[doc = include_str!(\"../docs/reference/{wanted}\")]",
                    spec.rust("name").unwrap_or_default()
                ));
            }
        }
        if spec.body.pointer("/result/bulk").is_some() {
            if let Some((_, bulk)) = bulk_function(spec, &rust) {
                let wanted = format!("{}.bulk.md", spec.stem());
                if !bulk.includes.iter().any(|i| i.ends_with(&wanted)) {
                    problems.push(format!(
                        "{}: add #[doc = include_str!(\"../docs/reference/{wanted}\")]",
                        spec.rust("bulk").unwrap_or_default()
                    ));
                }
            }
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

#[test]
fn the_unit_rule_matches_pymzlib() {
    // The same cases pyMzLib's rule is written for: a unit, its singular, an alias, and a short
    // symbol that must not match inside a word.
    assert!(unit_mentioned("scans", "Skip this many scans."));
    assert!(unit_mentioned("records", "one record per line"));
    assert!(unit_mentioned("min", "Retention time in minutes."));
    assert!(unit_mentioned("m/z", "The m/z of the ion."));
    assert!(!unit_mentioned("ms", "the msalign format"));
    assert!(unit_mentioned("ms", "in ms, as the file writes it"));
    assert!(unit_mentioned(
        "intensity (instrument units)",
        "Apex intensity."
    ));
    assert!(!unit_mentioned("scans", "Skip this many."));
}
