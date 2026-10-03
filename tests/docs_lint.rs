//! Hold the guides to the promises the docs make about themselves. pyMzLib's
//! `pkg/python/tests/test_docs_lint.py`, for this crate.
//!
//! A module's `//!` documentation is its guide: it is what docs.rs shows when a reader opens the
//! module, and it is where the worked examples live.
//!
//! * **Every guide opens with a question -> function -> mzLib type table**, so a reader can find
//!   the call for their question before reading the prose.
//! * **Every guide says what to cite**, by including the `cite.<module>.md` fragment that
//!   `tests/spec_docs.rs` renders from the specs' DOIs.
//! * **Every example in a guide runs, or says why it does not.** Doctests run in CI against the
//!   replay bridge. A `no_run` or `ignore` block is allowed only when the prose just before it says
//!   `Not run:` and why, so a reader can tell an example that was checked from one that was not.
//! * **No counts or mzLib versions in prose.** "36 formats" was once written in seven places and
//!   went stale in all of them with one mzLib release. A count belongs in an executed example or a
//!   generated table; a version belongs in the changelog.
//!
//! std only, so it runs in every test job.

use std::path::{Path, PathBuf};

/// The modules whose `//!` documentation is a guide. `tests/spec_docs.rs::guide_module` maps each
/// spec module onto one of these.
const GUIDES: &[&str] = &[
    "flashlfq",
    "isobaric",
    "peptidoform",
    "pride",
    "proteins",
    "readers",
    "sdrf",
    "stats",
];

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The guides that exist: a module is held to the rules from the commit that adds it.
fn guides() -> Vec<(&'static str, PathBuf)> {
    GUIDES
        .iter()
        .map(|name| (*name, root().join("src").join(format!("{name}.rs"))))
        .filter(|(_, path)| path.exists())
        .collect()
}

/// The `//!` lines of a source file, with the marker and one following space removed, numbered
/// by their line in the file.
fn module_doc(path: &Path) -> Vec<(usize, String)> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .enumerate()
        .filter_map(|(i, line)| {
            let rest = line.trim_start().strip_prefix("//!")?;
            Some((i + 1, rest.strip_prefix(' ').unwrap_or(rest).to_owned()))
        })
        .collect()
}

fn rel(path: &Path) -> String {
    path.strip_prefix(root())
        .unwrap_or(path)
        .display()
        .to_string()
        .replace('\\', "/")
}

/// Each fenced block's opening line, as (line number, info string), and every prose line outside a
/// block, as (line number, text).
#[allow(clippy::type_complexity)]
fn split_blocks(lines: &[(usize, String)]) -> (Vec<(usize, String, usize)>, Vec<(usize, String)>) {
    let mut fences = Vec::new();
    let mut prose = Vec::new();
    let mut inside: Option<usize> = None;
    for (index, (number, line)) in lines.iter().enumerate() {
        let trimmed = line.trim_start();
        let ticks = trimmed.chars().take_while(|c| *c == '`').count();
        if ticks >= 3 {
            let info = trimmed[ticks..].trim();
            match inside {
                None => {
                    inside = Some(ticks);
                    fences.push((*number, info.to_owned(), index));
                }
                Some(open) if ticks == open && info.is_empty() => inside = None,
                Some(_) => {}
            }
            continue;
        }
        if inside.is_none() {
            prose.push((*number, line.clone()));
        }
    }
    (fences, prose)
}

#[test]
fn every_guide_opens_with_a_question_table() {
    let mut problems = Vec::new();
    for (name, path) in guides() {
        let doc = module_doc(&path);
        let opens = doc
            .iter()
            .take(40)
            .any(|(_, line)| line.trim_start().starts_with("| You want"));
        if !opens {
            problems.push(format!(
                "{}: the {name} guide must open with a '| You want to ... | Call | mzLib type |' \
                 table in its first 40 lines of //! documentation",
                rel(&path)
            ));
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

#[test]
fn every_guide_says_what_to_cite() {
    let mut problems = Vec::new();
    for (name, path) in guides() {
        let text = std::fs::read_to_string(&path).unwrap();
        let include = format!("#![doc = include_str!(\"../docs/reference/cite.{name}.md\")]");
        if !text.contains(&include) || !text.contains("//! ## Cite") {
            problems.push(format!(
                "{}: end the guide with a '## Cite' heading followed by {include}",
                rel(&path)
            ));
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

#[test]
fn every_guide_example_runs_or_says_why_not() {
    let mut problems = Vec::new();
    let mut pages: Vec<PathBuf> = guides().into_iter().map(|(_, path)| path).collect();
    pages.push(root().join("src").join("lib.rs"));
    for path in pages {
        let doc = module_doc(&path);
        let (fences, _) = split_blocks(&doc);
        for (number, info, index) in fences {
            let words: Vec<&str> = info.split([',', ' ']).map(str::trim).collect();
            if !words.iter().any(|w| *w == "no_run" || *w == "ignore") {
                continue;
            }
            let says_why = doc[index.saturating_sub(4)..index]
                .iter()
                .any(|(_, line)| line.contains("Not run:"));
            if !says_why {
                problems.push(format!(
                    "{}:{number}: a `{info}` example must run, or the prose just before it must \
                     say 'Not run:' and why",
                    rel(&path)
                ));
            }
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

const NUMBER_WORDS: &[&str] = &[
    "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve",
    "thirteen", "nineteen",
];
const COUNTED: &[&str] = &[
    "file types",
    "file formats",
    "formats",
    "wire verbs",
    "verbs",
    "readers",
];

/// pyMzLib's `COUNT` pattern: a number, then a noun that counts formats or verbs.
fn counts_in(line: &str) -> Vec<String> {
    let lower = line.to_ascii_lowercase();
    let words: Vec<&str> = lower
        .split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|w| !w.is_empty())
        .collect();
    let mut found = Vec::new();
    for (i, word) in words.iter().enumerate() {
        let is_number = word.chars().all(|c| c.is_ascii_digit()) || NUMBER_WORDS.contains(word);
        if !is_number {
            continue;
        }
        let rest = words[i + 1..].join(" ");
        if let Some(noun) = COUNTED
            .iter()
            .find(|noun| rest == **noun || rest.starts_with(&format!("{noun} ")))
        {
            found.push(format!("{word} {noun}"));
        }
    }
    found
}

/// pyMzLib's `VERSION` pattern: `mzLib 1.0.<n>`.
fn versions_in(line: &str) -> Vec<String> {
    let mut found = Vec::new();
    for (start, _) in line.match_indices("mzLib") {
        let rest = line[start + 5..].trim_start();
        let rest = rest.strip_prefix('v').unwrap_or(rest);
        if let Some(tail) = rest.strip_prefix("1.0.") {
            if tail.starts_with(|c: char| c.is_ascii_digit()) {
                let digits: String = tail.chars().take_while(char::is_ascii_digit).collect();
                found.push(format!("mzLib 1.0.{digits}"));
            }
        }
    }
    found
}

#[test]
fn the_count_and_version_rules_match_pymzlib() {
    assert_eq!(counts_in("all 36 file types mzLib knows"), ["36 file types"]);
    assert_eq!(counts_in("Nineteen formats offer no view"), ["nineteen formats"]);
    assert!(counts_in("19 of the 38 have it").is_empty());
    assert!(counts_in("36 columns").is_empty());
    assert_eq!(versions_in("since mzLib 1.0.593 (#1388)"), ["mzLib 1.0.593"]);
    assert_eq!(versions_in("mzLib v1.0.592"), ["mzLib 1.0.592"]);
    assert!(versions_in("mzLib #1388").is_empty());
}

#[test]
fn no_counts_or_mzlib_versions_in_prose() {
    let mut problems = Vec::new();
    let mut pages: Vec<(PathBuf, Vec<(usize, String)>)> = guides()
        .into_iter()
        .map(|(_, path)| {
            let doc = module_doc(&path);
            (path, doc)
        })
        .collect();
    let lib = root().join("src").join("lib.rs");
    let lib_doc = module_doc(&lib);
    pages.push((lib, lib_doc));
    let readme = root().join("README.md");
    let readme_lines = std::fs::read_to_string(&readme)
        .unwrap()
        .lines()
        .enumerate()
        .map(|(i, line)| (i + 1, line.to_owned()))
        .collect();
    pages.push((readme, readme_lines));

    for (path, lines) in pages {
        let (_, prose) = split_blocks(&lines);
        for (number, line) in prose {
            if line.trim_start().starts_with("<!--") {
                continue;
            }
            for found in counts_in(&line) {
                problems.push(format!(
                    "{}:{number}: a count of formats or verbs in prose: '{found}'",
                    rel(&path)
                ));
            }
            for found in versions_in(&line) {
                problems.push(format!(
                    "{}:{number}: an mzLib version in prose: '{found}'",
                    rel(&path)
                ));
            }
        }
    }
    assert!(
        problems.is_empty(),
        "Counts and versions go stale with the next mzLib release. Show the count in an executed \
         example or a generated table; put a version in CHANGELOG.md.\n{}",
        problems.join("\n")
    );
}
