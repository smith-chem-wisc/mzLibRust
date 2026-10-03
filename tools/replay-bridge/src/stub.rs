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
//! * for `quant flashlfq` (`STDIN_ECHO`), the runs sent on stdin must be the recording's runs, and
//!   for `peptidoform convert` the sequences sent must be the recording's `input` column, in order.
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
    /// The recording's `input` column, for `peptidoform convert` (`STDIN_ECHO`).
    inputs: &'static [&'static str],
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
const STDIN_ECHO: &[&str] = &["quant flashlfq", "peptidoform convert"];

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
        } else if recorded != *value {
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
    if verb == "peptidoform convert" {
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
        let why = mismatch(recording, &options).or_else(|| {
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
