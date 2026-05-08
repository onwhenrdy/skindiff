# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`skindiff` — a C++14 numerical simulator for diffusion of active pharmaceutical ingredients (APIs) through skin. Models a stack of compartments (vehicle / donor → stratum corneum → deeper skin layers → sink) and solves the resulting tri-diagonal diffusion system with Crank-Nicolson time stepping. All code lives in `namespace sc` (Scientific Consilience GmbH).

There are two front-ends to the same simulation core:
- **`dskincmd`** — a standalone CLI executable built from this directory (`main.cpp` → `Session` → `SystemCmd`).
- **`dskin`** — an R/Rcpp package in `scripts/dskin/` that wraps the same core (`SystemR` in `scripts/dskin/src/dskin_rbinding.cpp`).

Both front-ends extend the `System` base class and only override `progressCallback` / `initRun` / `tearDownRun` / `testForStop`.

## Build & run

### CLI (`dskincmd`)
```
cmake -DCMAKE_BUILD_TYPE=Release .
make
./dskincmd --template          # write dskin_config.json template
./dskincmd dskin_config.json   # run from JSON config
./dskincmd --version
```

CMake quirks to be aware of:
- The Windows compiler-flag block is gated on `if (WINDOWS)`, which is **not** a standard CMake variable (the conventional one is `WIN32`); on Windows builds the inner flags are skipped unless you pass `-DWINDOWS=1`.
- Build number is read via `hg id -i` (Mercurial); if `hg` isn't on PATH the build still succeeds with `BUILD_NR=unknown`.
- `include_directories(SYSTEM includes)` references a directory that doesn't exist in the tree — harmless but misleading.
- Linking requires zlib (`-lz`) for the gzipped log output (`zstr.h`).

### CLI positional-argument form
Besides `--template` / `--version` / config-file, `CmdLineParser` accepts 19, 20, 21, or 23 positional parameters (last is always the file tag). The full list is printed by running `dskincmd` with no args; key fields are `C_0`, `D_Donor`, `D_SC`, `D_DSL`, `K_SC/Don`, `K_DSL/Don`, app area, lipid/DSL cross-sections, layer heights, sim time, resolution, scaling (mg/ug/ng), discretization (`EQUIDIST`/`BK`), matrix-builder method (`DSkin_1_3`/`DSkin_1_4`), finite-dose (yes/no), optional remove-at, optional replicate-after, optional `Vd` + `t_half` (enables PK sink).

### Tests (Catch2)
Tests are a separate qmake project at `tests/Tests.pro` and only exercise `tdmatrix.cpp` and the linear-solver routines in `algorithms.h` — they are **not** wired into the CMake build.
```
cd tests
qmake Tests.pro && make
./Tests                       # all tests
./Tests "[ThomasAlgorithm]"   # single test by tag
./Tests "[GaussReUsePivot]"
```
Catch2 is vendored as `tests/catch.hpp`; `tests/main.cpp` is just `CATCH_CONFIG_MAIN`.

### R package (`dskin`)
Built from `scripts/create_package.R`, which:
1. Copies `*.cpp`/`*.h` from the repo root into `scripts/dskin/src/`, **excluding** `main.cpp`, `session.{cpp,h}`, `systemcmd.{cpp,h}`, `consoleprogressbar.{cpp,h}` (these are CLI-only).
2. Runs `pkgbuild::compile_dll()` and `devtools::check()`.
3. Calls `devtools::build()`, then `clean.up.src()` removes the copied sources from `src/` (keeping only `RcppExports.cpp`, `dskin_rbinding.cpp`, `Makevars`).

Run from the `scripts/` directory:
```r
source("create_package.R")
```
A usage example (parameter-pack construction, `dskin.simulate`, plot/movie helpers) lives in `scripts/example.R`.

If you edit core C++ files, the changes are only picked up by the R package on the next `create_package.R` run — there is no symlink. The CLI (CMake) build is unaffected by the R-package layout.

## Architecture

### Simulation core (compiled into both CLI and R package)
- **`Parameter`** (`parameter.{h,cpp}`) — aggregate of `SystemParameter`, `LogParameter`, `PKParameter`, `SinkParameter`, `VehicleParameter`, and `std::vector<LayerParameter>`. Two parsers fill it: `CmdLineParser` (positional args) and `JsonParser` (JSON config; uses vendored `json.h` from nlohmann).
- **`Geometry`** (`geometry.{h,cpp}`) — discretizes the compartment stack into space steps. Two methods: `EQUI_DIST` (uniform) and `B_AND_K` (refined; Babucke & Kloker, 2009). Stores per-cell `m_space_steps`.
- **`Compartment`** / **`Sink`** — one compartment per layer plus optional sink (`Perfect_Sink` or `PK_Compartment` with `Vd` and `t_half`).
- **`MatrixBuilder`** (`matrixbuilder.{h,cpp}`) — assembles RHS and LHS tri-diagonal matrices for Crank-Nicolson. Three methods exposed via `Method` enum:
  - `DSkin_1_3` — central element concentrations with back-flux correction.
  - `DSkin_1_4` — element-edge concentrations (Crank "method of differences").
  - `DSkin_1_5` — fast variant of 1_4.
  Inline helpers `backFluxCorrection`, `areaCorrection`, `harmMeanFromIdx` implement the partition/area handling at layer boundaries.
- **`TDMatrix`** (`tdmatrix.{h,cpp}`) — tri-diagonal matrix with optional super-upper diagonal and pivot index for the LU-reuse path. Provides `inlineMultiply` to avoid allocation in the hot loop.
- **`algorithms.h`** — header-only inline solvers used inside the time loop:
  - `thomasIP` / `thomasReUseIP` — Thomas algorithm; the "ReUse" variant caches the LU factorization via `TDMatrix::isPrepared()`.
  - `gaussPivotIP` / `gaussReUsePivotIP` — Gauss with partial pivoting (used when the matrix isn't diagonally dominant).
- **`System`** (`system.{h,cpp}`) — owns compartments, geometry, matrix builder, and loggers. Constructor builds the geometry and matrices from the `Parameter`. `run()` executes the time loop: for each minute it does `n_ts` sub-steps of `rhs_matrix.inlineMultiply(c) → thomasReUseIP(lhs_matrix, c)`, then optionally replaces the donor compartment (`replace_after`) or removes it (`remove_at`, which forces a matrix rebuild and shifted geometry indices).
- **Loggers** — `CompartmentLog2D` writes mass-vs-time per compartment (and the sink); `CompartmentLog3D` writes concentration-depth profiles. Both can gzip output via `zstr.h` (`strict_fstream.h`).

### Front-end specifics (CLI only)
`session.cpp`, `systemcmd.cpp`, `consoleprogressbar.cpp` provide the command-line entry, the progress bar, and the after-run log-file flush. These are excluded from the R package on purpose (R provides its own progress reporting and interrupt handling in `SystemR`).

### Front-end specifics (R only)
`scripts/dskin/src/dskin_rbinding.cpp` defines `SystemR : public System`, exposes `.dskin.binding.simulate(json, write_to_R, write_to_files)` and `.dskin.binding.geometry(json)` via Rcpp, and converts the C++ loggers into R `data.frame`-like lists. The R-side glue (`scripts/dskin/R/dskin_bindings.R`) parallelizes a list of parameter packs via `foreach %dopar%`.

### Where parameters flow
```
CLI args / JSON   ─► CmdLineParser / JsonParser ─► Parameter
                                                       │
                                          ┌────────────┴────────────┐
                                          │                          │
                                       SystemCmd                  SystemR
                                          │                          │
                                          └─────► System ◄───────────┘
                                                    │
                              Geometry, Compartments, Sink, MatrixBuilder
                                                    │
                                                  run()  ── time loop ──► loggers
```
