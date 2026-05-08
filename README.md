# skindiff

Numerical simulation of drug diffusion through skin, exposed as an R
package. The C++ engine builds a one-dimensional Crank-Nicolson
finite-difference model of a stack of compartments
(vehicle / stratum corneum / deeper skin layers / sink) and returns
time-resolved mass and concentration-depth profiles.

## Install

The package builds like a normal Rcpp package:

```r
# from the repository root
remotes::install_local(".")
# or
devtools::install()
```

You need a C++17-capable toolchain (Rtools 4.x on Windows, the matching
toolchain on macOS/Linux).

## Quick start

```r
library(skindiff)

p <- skin_params(
  vehicle = list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0,
                 log_cdp = TRUE),
  layers = data.frame(
    name          = c("Stratum corneum", "Deeper skin layers"),
    height        = c(190L, 200L),
    D             = c(28.25, 5767.78),
    K             = c(421.5, 0.0472),
    cross_section = c(0.001, 0.3),
    log_cdp       = c(TRUE, TRUE)
  ),
  sink     = list(Vd = 1.875e6),         # PK distribution volume, ml
  pk       = list(enabled = TRUE, thalf = 24),
  sim_time = 24 * 60,                    # one day, in minutes
  scaling  = "ng"
)

res <- skin_simulate(p)

res$mass            # time x compartment + sink, in ng
res$concentration   # mass / compartment volume
res$cdp$`Stratum corneum`$conc  # [depth_um, time] matrix, ng/ml
```

## Result structure

`skin_simulate()` returns a `skin_result` object with:

- `mass` -- a data.frame with columns `time` (min) and one column per
  logged compartment, holding the integrated mass in the chosen scaling
  unit (`mg`, `ug`, or `ng`).
- `concentration` -- the same shape as `mass`, but each value is divided
  by the compartment volume. The vehicle column is the donor
  concentration over time.
- `cdp` -- per-compartment concentration-depth profiles. Each entry has
  a `time` vector, a `depth_um` vector, and a numeric matrix `conc`
  indexed `[depth, time]` in `<scaling>/ml`.
- `geometry` -- mesh information (`min_step_um`, `max_step_um`, `eta`,
  `n_cells`).
- `params`, `runtime_s`, `status`.

## Tests

```r
devtools::test()
# or
testthat::test_local()
```

Both the C++ engine and the R wrapper are tested. The C++ tests use
testthat's Catch wrapper and run via `tests/testthat/test-cpp.R`.

## License

GPL (>= 3). See `LICENSE`.
