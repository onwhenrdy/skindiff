# skindiff

Numerical simulation of drug diffusion through skin, exposed as an R
package. A C++17 engine solves the one-dimensional diffusion equation
through a stack of compartments
(vehicle / stratum corneum / deeper skin layers / sink) and returns
time-resolved mass and concentration-depth profiles. The R surface
covers parameter construction, result accessors, ggplot2 plotting,
and parameter fitting from permeation / penetration data.

Strict-units throughout: every unit-bearing argument is a
[`units`](https://r-quantities.github.io/units/) object built from a
package helper (`um()`, `hours()`, `mg_per_ml()`, ...), and every
numeric output column carries a unit. Bare numerics are rejected at
unit-bearing interfaces, by design.

Current version: **0.10.0**.

## Install

Hard dependencies: `Rcpp (>= 1.0.0)`, `cli (>= 3.6.0)`, `units (>= 0.8)`.
Suggests (only needed for plotting and tests): `ggplot2`, `testthat`.

```r
install.packages(c("Rcpp", "cli", "units"))
remotes::install_local(".")
# or
devtools::install()
```

You need a C++17-capable toolchain (Rtools 4.x on Windows, the matching
toolchain on macOS/Linux).

## Quick start

```r
library(skindiff)

# A vehicle, two skin layers, a perfect Franz-cell receptor.
veh <- vehicle(c_init = mg_per_ml(1.0), height = um(50),  D = um2_per_min(1000))
sc  <- layer("Stratum corneum", height = um(20),  D = um2_per_min(2),
             K = 50,  cross_section = 1.0, log_cdp = TRUE)
der <- layer("Dermis",          height = um(120), D = um2_per_min(200),
             K = 1.0, cross_section = 1.0, log_cdp = TRUE)

p   <- skin_params(area = cm2(1), vehicle = veh, layers = list(sc, der),
                   sink = perfect_sink("Receptor"),
                   duration = hours(24L), resolution = 4L, scaling = "ug")
res <- skin_simulate(p)

# Derived quantities
metrics(res)               # one-row data.frame: J_ss, t_lag, K_p, AUC, ...
permeated(res)             # data.frame(time, Q) at the receptor
profile_at(res, hours(8))  # depth profile slice at t = 8 h

# Plots (requires ggplot2 in Suggests)
library(ggplot2)
autoplot(res, what = "permeated")
autoplot(res, what = "profile", n_times = 6)

# Fitting layer D and K from observed data
perm <- permeation_obs(data.frame(
  time       = hours(c(1, 2, 4, 8, 24)),
  q_per_area = ug_per_cm2(c(0.05, 0.18, 0.42, 0.81, 1.6))
))
fit <- skin_fit(template     = p,
                observations = list(permeation = perm),
                fit_pars     = list("Stratum corneum" = c("D", "K")))
coef(fit)
autoplot(fit)
```

## Algorithm

The engine is a cell-centred **finite-volume** discretisation of the
1-D diffusion equation, time-stepped with **Crank-Nicolson**, solved
with a **Thomas tri-diagonal solver**. Three design choices are worth
calling out:

- **State variable is activity, not concentration.** Each cell stores
  `u = c/K`. At a partition discontinuity the equilibrium condition
  `c_L/K_L = c_R/K_R` says `u` is continuous even though `c` jumps —
  so the K-jump moves from the variable into the FVM coefficient
  `κ = K·A·D`, where the standard harmonic-mean face treatment
  handles it cleanly. Two-layer steady-state K-jump profiles come out
  exact to machine precision.

- **Per-compartment graded mesh.** Each compartment uses uniform cells
  of size `dx_i = (1/resolution) * sqrt(D_i / D_min)`. The smallest-D
  compartment uses cells of size `1/resolution` µm; higher-D
  compartments get proportionally coarser cells, since their gradient
  scale in `u` is correspondingly wider. One knob, no transition zones.

- **Boundary conditions are cell-edge Dirichlet.** Both top (donor
  surface, infinite-dose mode) and bottom (membrane↔sink) use the
  `κ_outside → ∞` limit of the harmonic mean, giving a true Dirichlet
  BC at the interface independent of the donor / sink properties.

The hot path in `crankNicolsonStepIP` fuses the matrix-vector product
with the forward Thomas sweep in a single pass, and stores the
prepared LHS diagonal as its reciprocal so the sweep is multiply-only.

The scheme is validated against closed-form solutions from Crank
(*Mathematics of Diffusion*) and Kasting 2001 — see
`tests/testthat/test-analytical.R` and the helpers in
`tests/testthat/helper-analytical.R`. Convergence is clean
second-order on a uniform mesh.

## Result structure

`skin_simulate()` returns a `skin_result` list. **All numeric fields
carry units.**

- `mass` — data.frame. `time` is `[min]`; per-compartment columns carry
  the chosen scaling unit (e.g. `[ng]` for `scaling = "ng"`).
- `concentration` — derived data.frame. `time` is `[min]`;
  per-compartment columns are `[scaling/ml]`. The sink column is `NA`
  for `perfect_sink()` (mass / fictitious Vd would be misleading).
- `cdp` — named list, one entry per compartment with `log_cdp = TRUE`.
  Each entry has `time` (`[min]`), `depth` (`[um]`), and a units-bearing
  matrix `conc[depth, time]` carrying `[scaling/ml]`.
- `geometry` — `min_step` (`[um]`), `max_step` (`[um]`), `n_cells`.
- `params`, `runtime` (`[s]`), `status`, `scaling`.

Operations that don't compose with units (`rowSums` over units columns
segfaults R; `lm`/`approx` need bare numerics) require an explicit
`as.numeric()` strip at the boundary.

## API surface

Twelve exported functions plus unit helpers and S3 methods:

| Group | Functions |
|---|---|
| Unit helpers | `um`, `mm`, `cm`, `cm2`, `mm2`, `ml`, `mg_per_ml`, `ug_per_ml`, `ng_per_ml`, `mg_per_cm2`, `ug_per_cm2`, `ng_per_cm2`, `um2_per_min`, `cm2_per_s`, `seconds`, `minutes`, `hours`, `days` |
| Compartment builders | `vehicle()`, `layer()`, `perfect_sink()`, `finite_sink()` |
| Composer + runner | `skin_params()`, `skin_simulate()` |
| Result accessors | `permeated()`, `flux()`, `permeated_at()`, `profile_at()`, `metrics()` |
| Observations + fit | `permeation_obs()`, `penetration_obs()`, `skin_fit()`, `skin_params_from_fit()` |
| S3 methods | `print` / `summary` for the classed objects; `coef` / `residuals` / `fitted` for `skin_fit`; `autoplot` for `skin_result` and `skin_fit` (via `ggplot2::autoplot` in Suggests) |

## Tests

```r
devtools::test()
# or
testthat::test_local()
```

Both the C++ engine and the R wrapper are tested. The C++ tests use
testthat's Catch wrapper and run via `tests/testthat/test-cpp.R`. The
analytical battery validates against Crank single-slab `Q(t)` and
transient CDP, two-layer steady-state K-jump CDP, semi-infinite erfc,
and Kasting finite-dose `M_t / M_∞` and membrane CDP.

## License

GPL (>= 3). See `LICENSE`.
