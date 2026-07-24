<!-- Rendered by /project. Do not hand-edit — edits are overwritten on the next status/advance/next. -->

# RESUME — mzLibRust

**Phase:** BUILD · created 2026-07-23 · adopted 2026-07-24 · last rendered 2026-07-24

**Goal:** mzLib callable from Rust — a separate crate, the Rust sibling of pyMzLib, on the SAME
language-neutral bridge (D6). A port of a proven design. This workspace also hosts the emerging
**mzLibR** (R) sibling.

**Part of the `bridge` umbrella** (`E:/CodeReview/bridge`) — shared contract, principle, and the
upstream-mzLib fix queue live there; this project owns the Rust product.

## Pick up at

Read `PLAN.md` (§1 = the core idea: it's cheap because the bridge is already language-neutral). Then
the **§8 grill list**, led by the one genuinely new decision vs pyMzLib:

- **`RG-dist`** — distribution. A crate can't carry the ~130 MB .NET payload, so the bridge binary is
  **fetched at build time** (PLAN §3). Grill + implement.
- **`RG-grill`** — the rest of the §8 open decisions, before M2/M3.
- **`RG-mzlibr`** — mzLibR (R sibling) is being edited from this workspace; decide if it splits into
  its own project.

## What this session contributed upstream

Two mzLib fixes surfaced from mzLibRust's parity/bake-off work and are tracked in `bridge/UPSTREAM.md`:
- [#1114](https://github.com/smith-chem-wisc/mzLib/pull/1114) — ETD/ECD product types (a Rust parity
  test caught `y` ions that shouldn't exist). Reviewed + fixed.
- [#1119](https://github.com/smith-chem-wisc/mzLib/pull/1119) / issue #1113 — digestion cleavage-blocking
  modifications (a biologist bake-off arm caught an impossible peptide). New issue + PR.

## Locked decisions

- **RUST-D6** — reuse pyMzLib's language-neutral wire contract unchanged.
- **RUST-DIST** — download-at-build distribution (crate can't bundle ~130 MB).
- **RUST-TYPES** — `Option<f64>` vs `f64` is the Rust projection of the bridge availability principle.

## Phase gate (BUILD) — advisory

- Crate scaffolded (Cargo.toml/lock present) ✓ · repo live (`smith-chem-wisc/mzLibRust`) ✓ · PLAN written ✓
- Not yet: distribution decision settled + implemented; §8 grill worked; M2/M3.
