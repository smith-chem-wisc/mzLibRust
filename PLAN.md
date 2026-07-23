# mzLibRust — plan

**mzLib callable from Rust.** A separate repo and crate, the Rust sibling of pyMzLib, built on the
**same language-neutral bridge** pyMzLib already uses. Written 2026-07-23 from everything the pyMzLib
build taught us — so this is not a green-field design, it's a port of a proven one.

> Pick-up for a fresh chat: read this, then read the three pyMzLib reference docs it points to
> (`code/pyMzLib/docs/contributing/conventions.md`, `design/expansion-plans/flashlfq-tranche.md`,
> `design/bakeoff-flashlfq/DESIGN.md` under `E:\CodeReview\pyMzLib`). The single most important idea
> is in §1.

---

## 1. The whole thing is cheap because the bridge is already language-neutral (D6)

pyMzLib's architecture was deliberately designed, from day one, so a Rust binding would be cheap. The
`.NET` **bridge** is a self-contained executable that speaks a **versioned JSON envelope over
stdin/stdout** and assumes nothing about its caller. It already does all the genuinely hard parts:

- the mzLib interop and the composition of mzLib's own methods (digestion, fragmentation, PRIDE,
  `MakeIdentifications` + `FlashLfqEngine`);
- the availability-vs-correctness **error classification** (timeouts/5xx → `ServiceUnavailable`, the
  rest by exception type);
- keeping engine chatter off stdout so the envelope is clean;
- carrying its own .NET runtime so no .NET install is needed.

**mzLibRust does none of that again.** It is a thin, idiomatic Rust crate that spawns the bridge,
writes stdin, reads one JSON line, and deserializes it. The port is: transport module + typed structs
+ an error enum + tests + docs. That's it.

The wire contract to target is exactly what pyMzLib consumes — see `_bridge.py` and the three verb
handlers (`Program.cs`, `Peptidoform.cs`, `Quantification.cs`) in the pyMzLib repo. Every field name is
the snake_case of an mzLib/FlashLFQ name; do not rename on the Rust side either.

## 2. Architecture map (pyMzLib → mzLibRust)

| pyMzLib | mzLibRust | notes |
|---|---|---|
| `_bridge.py` (only transport-aware module) | a `bridge` module | `std::process::Command`, write stdin, read stdout, `serde_json` the envelope. The ONLY place that knows the bridge exists. |
| typed dataclasses (`PrideFile`, `Digest`, `FlashLfqResults`, …) | `#[derive(Deserialize)]` structs (serde) | Mirror the wire, mzLib-named. `rename_all = "snake_case"` is unnecessary — the wire is already snake_case. |
| exception hierarchy (`PyMzLibError` → `UsageError`/`BridgeError`/`ServiceUnavailableError`/…) | one `MzLibError` enum (thiserror) | Map `envelope.error.type`: `"usage"`→`Usage`, `"ServiceUnavailable"`→`ServiceUnavailable`, else `Bridge{ error_type, message }`. Plus `Timeout`, `BridgeNotFound`, `Io`. |
| `pymzlib.pride.list_files(...)` | `mzlib::pride::list_files(...) -> Result<Vec<PrideFile>, MzLibError>` | Functions return `Result`, not exceptions — the natural Rust shape. |
| `PYMZLIB_BRIDGE` env override | `MZLIB_BRIDGE` env override | Same escape hatch for dev/offline. |

### The None/0 lesson maps perfectly onto Rust's type system

The hardest-won FlashLFQ lesson today: a **peptide** intensity is `0.0` when missing, a **protein**
intensity is `None` (NaN, "could not be resolved"). In Rust this is not a footnote, it's the types:

```rust
struct Peptide { intensities: HashMap<String, f64>, /* 0.0 = missing */ }
struct ProteinGroup { intensities: HashMap<String, Option<f64>>, /* None = unresolvable */ }
```

`Option<f64>` makes the distinction unignorable at the call site — a strictly better outcome than the
Python docs having to warn about it. Bake this in from the first line; it's a selling point.

## 3. The one genuinely new decision: distribution (grill this first)

Python solved the ~130 MB payload with a wheel that carries the binary. Rust has no wheel, and
crates.io rejects large binaries (~10 MB soft limit). This is the **only** part with no proven answer,
so it's the thing to `/grill-me` on before writing code. Options:

- **A. Download-at-build (`build.rs`)** — the build script fetches the platform bridge binary from a
  GitHub Release, checksum-verified, cached in `OUT_DIR`/a user cache. Crate stays tiny and
  publishable; keeps the zero-.NET-install promise. Cost: builds need network once (mitigate with the
  env override + a vendored/offline mode). This is what `ort`, `onnxruntime`-style crates do.
  **Recommended default.**
- **B. Locate-at-runtime + `MZLIB_BRIDGE`** — ship no binary; the user installs the bridge separately.
  Simplest crate, most user friction (contradicts the zero-friction goal). Good as the dev/offline
  fallback, not as the primary story.
- **C. Embed via `include_bytes!` / a `-sys` crate** — not viable at 130 MB.

Recommendation: **A with B as the escape hatch.** And note the two payload-shrink levers we already
found: mzLib issue **#1103** (TorchSharp/libtorch is ~238 MB of the payload and is dragged in
transitively — making it optional shrinks *both* bindings) and the **mzML-only native-reader prune**
(~20 MB of Thermo/Bruker readers). Both matter more for a download-at-build model.

### Where does the bridge binary come from?

Sub-decision under A: the bridge lives in `smith-chem-wisc/pyMzLib` today, and its CI already builds
it for all four platforms. Options, cleanest last:

1. mzLibRust's CI **builds the bridge itself** from the pyMzLib bridge source (clone + `publish-bridge.ps1`) — self-contained but duplicates build logic.
2. pyMzLib CI **also uploads the raw self-contained bridge binaries** as release assets (small
   addition to `wheels.yml`); mzLibRust's `build.rs` downloads those. Low-friction, shared source of truth.
3. **Extract the bridge into its own repo** `smith-chem-wisc/mzlib-bridge` that both pyMzLib and
   mzLibRust consume. The clean long-term end state (D6's logical conclusion), but a refactor of
   pyMzLib — defer until both bindings are real.

Recommend **(2) for the first cut, (3) as the eventual home.**

## 4. What we can leapfrog (because pyMzLib already paid the discovery cost)

pyMzLib discovered its bugs and doc gaps by shipping and running biologist bake-offs. mzLibRust starts
knowing all of them, so it can implement all three tranches *correctly from the start*:

- **FlashLFQ: expose `peaks` from day one** as the match-between-runs surface — the peptide roll-up
  drops most MBR transfers (140 peaks vs 52 at peptide level; a whole run's transfers invisible).
- **`Option<f64>` for protein intensity, `f64` for peptide** (§2).
- **Document MBR's requirements up front**: a complete/balanced design (same fractions across every
  condition & replicate), and `mbr_q_value_threshold` as the FDR control that makes transfers
  trustworthy (a pyteomics bake-off arm proved you get ~80% false transfers without it).
- **Follow mzLib/FlashLFQ names**; disclose the traps (census, isoform cap, fixed charges).
- Full `detection_type` vocabulary (`MSMS`, `MBR`, `MSMSIdentifiedButNotQuantified`,
  `MSMSAmbiguousPeakfinding`, `NotDetected`).

## 5. Conventions (from pyMzLib's `conventions.md` — most carry verbatim)

Carry as-is: one JSON envelope per call; error classification lives in the bridge (free for Rust);
large/variadic input on **stdin** not argv; **follow mzLib naming**; **disclose the traps**;
language-neutral wire; the bridge is a translation layer, not a second implementation; offline
fixtures by default + live canaries that skip (not fail) on outage.

Changes for Rust idiom:
- "Zero third-party runtime deps" (a Python-specific goal) → **minimal, well-established deps**:
  `serde`, `serde_json`, `thiserror`; nothing heavy or exotic. State a low MSRV.
- Errors are `Result<T, MzLibError>`, not exceptions.
- Prefer `Option<T>` where the wire can be null; `f64` where it is a real 0.

## 6. Roadmap (mirrors pyMzLib's milestones, but faster)

- **M0 — bridge round-trips.** `mzlib::bridge_version()` spawns the bridge, parses the envelope,
  checks the protocol version. Proves the whole transport + distribution story end to end (incl. a
  "no .NET installed" container test, like pyMzLib's).
- **M1 — PRIDE tranche.** The narrow vignette (`list_files`, `download`, `total_size_bytes`,
  `PrideFile`) — same "prove the packaging story" role pyMzLib's PRIDE M0 played.
- **M2 — peptidoform tranche.** `mzlib::peptidoform::fragments(accession)` → typed `Digest`.
- **M3 — FlashLFQ tranche.** `mzlib::flashlfq::quantify(...)` → typed `FlashLfqResults` with `peaks`,
  `Option<f64>` proteins, `mbr_rescued_peptide_count`. Correct from the start (§4).
- **Distribution/CI:** `build.rs` download-at-build; cargo build+test matrix (win/linux/macos);
  `docs.rs`; publish to **crates.io**. A biologist/rustacean bake-off to validate + find doc gaps,
  same methodology as `design/bakeoff-flashlfq/DESIGN.md`.

## 7. Rust-ecosystem prior-art check (do this before M2/M3 — the oracle move, Rust side)

Rust proteomics is nascent but not empty. Before building, check what exists so mzLibRust brings what's
missing rather than duplicating:
- **`mzdata`** — mzML/mzMLb/raw reading in Rust (the read layer a hand-rolled quant would use).
- **`rustyms`** — peptide fragmentation + ProForma in Rust. **Overlaps the peptidoform tranche** —
  check carefully; the value is mzLib's UniProt annotation census + digestion, not re-fragmenting.
- **`sage`** — a Rust search engine (identifications), not quant.
There is **no** mzLib-equivalent breadth (FlashLFQ, PRIDE client, deconvolution, the tested chemistry).
So the value prop — and the bake-off contrast — is even sharper than pyteomics-vs-pymzlib was.

## 7a. Project organization (SETTLED 2026-07-23)

- Each binding is its **own sibling folder + its own public repo** — `E:\CodeReview\mzLibRust`,
  later `E:\CodeReview\mzLibR` — **peers** of pyMzLib, never nested under it. Different toolchains
  (PyPI / crates.io / CRAN), CI, and release cadences; nesting would tangle git + CI and confuse
  `/project` (cwd = identity, so one folder must be one project).
- A **multi-language paper** ("mzLib for Python, Rust, and R — a language-neutral bridge") is its
  **own `/project` folder + repo** (e.g. `mzlib-bindings-paper`), referencing the binding repos as
  sources — the same way the lab already keeps code and papers in separate repos. Prefer **one**
  such paper over three per-language notes; pyMzLib's tentative Technical Note folds into it.
- The real shared substrate is the **bridge** → eventually extract `smith-chem-wisc/mzlib-bridge`
  (see §3), consumed by all bindings. Not a shared folder; a shared repo.
- R is the same pattern and about as cheap: call the bridge via `processx` + `jsonlite`.

## 8. Decisions to make in the new chat (the grill list)

1. **Distribution** (§3): download-at-build vs env-locate — and where the bridge binary comes from.
2. **Crate name** on crates.io (`mzlib`? check availability; the repo is `mzLibRust`). Repo:
   `smith-chem-wisc/mzLibRust` public, workspace private like pyMzLib.
3. **Sync vs async** API — the bridge is a subprocess; sync is fine and simplest; an optional async
   (tokio) layer can come later.
4. **MSRV** and the dependency floor.
5. Whether to **version the crate with the bridge** or independently (pyMzLib versions the wheel
   independently of the bridge; mirror that).

## 9. First steps for the new chat

1. `/project init` this folder (`E:\CodeReview\mzLibRust`); create the public product repo
   `smith-chem-wisc/mzLibRust` + the private workspace, mirroring pyMzLib's two-repo split.
2. `/grill-me` on the **distribution decision** (§3) — it shapes everything.
3. Prior-art check (§7): `mzdata`, `rustyms`, `sage`, and the mzLib seam.
4. M0: bridge round-trip + the no-.NET proof, proving distribution before any tranche.
5. Then M1 PRIDE, reusing the exact wire pyMzLib already consumes.

The user writes no Rust either (as with Python) — Claude owns all Rust authorship, idiom, packaging,
and CI. Explain tradeoffs in terms the user knows (C#/mzLib), the way `PYTHON_PRIMER.md` did for Python.
