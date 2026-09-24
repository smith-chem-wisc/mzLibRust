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
//!   no top-level `path`, per BULK.md), and a one-path call only a one-document recording.
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

fn option<'a>(options: &'a [(String, Option<String>)], name: &str) -> Option<&'a Option<String>> {
    options.iter().find(|(n, _)| n == name).map(|(_, v)| v)
}

/// Why this recording cannot be the answer to these options, or `None` if it can.
fn mismatch(recording: &Recording, options: &[(String, Option<String>)]) -> Option<String> {
    for (name, value) in options {
        let key = name.replace('-', "_");
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

fn answer(argv: &[String]) -> (bool, String) {
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
        match mismatch(recording, &options) {
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
    let (ok, envelope) = answer(&argv);
    let mut stdout = std::io::stdout();
    let _ = stdout.write_all(envelope.as_bytes());
    let _ = stdout.flush();
    std::process::exit(if ok { 0 } else { 2 });
}
