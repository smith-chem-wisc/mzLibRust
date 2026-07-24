# Project: mzLibRust

This folder is a `/project`-managed project (adopted 2026-07-24, non-destructive). **You are de
facto working on it** when cwd is here.

- **Phase:** BUILD
- **What it is:** mzLib callable from **Rust** — a separate crate, the Rust sibling of pyMzLib, built
  on the SAME language-neutral bridge (D6). A port of a proven design. This workspace also hosts the
  emerging **mzLibR** (R) sibling.
- **Part of the `bridge` umbrella** (`E:/CodeReview/bridge`): the shared wire contract, the projection
  principle, and the cross-cutting upstream-mzLib fix queue live there. This project owns the
  Rust-specific product (crate + distribution).
- **Pick up at:** read `PLAN.md` (§1 is the key idea). Then the §8 grill list — chiefly the one new
  decision vs pyMzLib: **distribution** (a crate can't carry the ~130 MB payload → download-at-build,
  §3). Gaps: `RG-dist`, `RG-grill`, `RG-mzlibr`.

Memory key `E--CodeReview`. Live state: `.project/state.yaml`; human-readable `RESUME.md` (rendered by
`/project`, do not hand-edit). Canonical design doc is `PLAN.md` at root. The Rust type projection of
the bridge principle: `Option<f64>` where a format may not provide a value, `f64` where it always does.
