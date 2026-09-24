# AGENTS.md

`pmsim` is a Rust FEM library for solids, structures, and porous media. It is a thin layer over the sibling
`russell_*` crates (linear algebra, nonlinear solver, PDE, sparse) and `gemlab` (meshes). Much behavior lives in
those crates; check them (usually `../russell/russell_*`) when something isn't in this repo.

## Setup

- Native libraries are required to link: `liblapacke-dev`, `libopenblas-dev`, `libsuitesparse-dev` (Ubuntu names
  from `.github/workflows/ubuntu.yml`).
- The default linear solver is UMFPACK from SuiteSparse (`Config::lin_sol_genie` defaults to `Genie::Umfpack`).
- `Cargo.toml` features swap the solver/BLAS backend: `intel_mkl`, `local_sparse`, `cudss`. Do **not** run
  `cargo test --all-features` unless Intel MKL / CUDA / locally built sparse libs are installed; use plain
  `cargo test` (what CI uses). `all.bash` assumes the full-feature environment.
- To build against local `russell` checkout, add the `[patch.crates-io]` block from README "Development".

## Commands

- Format: `cargo fmt` (rustfmt `max_width = 120`, see `rustfmt.toml`). No clippy/lint job exists in CI.
- All tests: `cargo test` (add `-- --nocapture` to see solver output).
- Unit tests only: `cargo test --lib`.
- One integration group: `cargo test --test <group>`, where `<group>` is a directory name under `tests/`:
  `boundary_conditions`, `elasticity`, `extrapolation`, `frames`, `heat`, `material`, `plasticity`, `seepage`, `spo`.
- One test: `cargo test --test spo spo_755_tensile -- --nocapture`. Function names mostly match the filename
  (a few are prefixed `test_`, e.g. `test_von_mises_single_element_2d`).
- **`yscripts/*.bash` are partly stale**: they reference non-existent test targets like
  `--test test_spo_751_pres_cylin`. Use the `cargo test --test <group> <fn>` form instead.
- Examples: `cargo run --release --example spo_751_pres_cylin -- -g mumps --arclength`. `-g` selects the solver
  (`mumps`/`umfpack`/`cudss`); other flags include `--lmm`, `--arclength`, `--bordering`, `--residual`.

## Architecture

- `src/lib.rs` exposes modules `analytical`, `base`, `fem`, `material`, `util`, plus `prelude` (the public entry:
  `use pmsim::prelude::*`).
- Simulation flow:
  - Nonlinear/general: `Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)` then
    `sim.steady(&mut data, IniDir::Pos, stop, dll)?`.
  - Linear-only: `SimulatorLin::new(...)` then `sim.steady(&mut data, post_compute_second_values)?`.
  - Build order: `Schema` (element type + DOF numbering + per-cell attributes) -> `Param*` material structs ->
    `Config` -> essential/natural BCs (`BcEssential`, `BcNatural`).
- Functional areas: `src/fem/` (elements, state, assembly, simulator, Paraview/PostProc output),
  `src/material/` (linear elastic, von Mises, hardening-softening, elastoplastic explicit/implicit),
  `src/base/` (config, schema, BCs, parameters, sample meshes), `src/util/` (reference data for tests).
- Errors are `Result<_, pmsim::StrError>` (a `&'static str`).

## Tests and data

- Integration tests live at `tests/<group>/main.rs` (a `mod` list) plus one file per test. Each test reads meshes
  from `data/` relative to the crate root (e.g. `data/spo/spo_755_tensile.msh`) and writes results to `/tmp/pmsim/...`
  (the `DIR` constant). No external services; re-running is safe.
- Numerical checks use `russell_lab::approx_eq` against expected values in README and/or `src/util/reference_data*.rs`.
- `src/bin/pmsim2pv2d.rs` / `pmsim2pv3d.rs` convert simulation output into VTU/PVD files for ParaView.
