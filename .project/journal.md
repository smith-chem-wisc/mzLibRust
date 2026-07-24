# Journal — mzLibRust

Append-only transition log: date · phase change · decision · skip.

- **2026-07-24** · `adopt` (non-destructive) · phase inferred BUILD. Retrofitted the existing Rust
  crate + repo as a `/project`. Seeded state from `PLAN.md` (RUST-D6 reuse the language-neutral
  contract; RUST-DIST download-at-build; RUST-TYPES Option<f64> projection). Linked under the `bridge`
  umbrella. Recorded that mzLibRust's parity/bake-off work surfaced mzLib #1114 and #1113/#1119
  (tracked in bridge/UPSTREAM.md). Nothing existing moved or clobbered; PLAN.md left at root.
