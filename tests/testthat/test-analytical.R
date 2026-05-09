# Analytical-solution validation against Crank, "The Mathematics of Diffusion"
# (2nd ed.). The engine should reproduce these closed forms accurately on the
# meshes we ship.
#
# Baseline numbers (graded mesh, Crank-Nicolson, Thomas-reuse, resolution = 4,
# true Dirichlet donor BC: all donor cells clamped + cell-edge alpha):
#
#   slab-lag-time L_inf  1.0e-6
#   slab slope_rel       8.2e-5
#   slab lag_rel         7.0e-4
#   ss-2L L1 L_inf       8.3e-13   <- machine precision
#   ss-2L L2 L_inf       8.4e-13   <- machine precision
#   ss-2L slope_rel      3.5e-12
#   erfc L_inf           1.8e-2
#   slab convergence     rate=2.00 <- clean 2nd-order, abs_err halves at
#                                     each dx -> dx/2 refinement
#
# The steady-state two-layer test hits machine precision because the
# symmetric tri-diagonal in u recovers piecewise-linear profiles exactly.
# The lag-time and convergence numbers cleanly demonstrate second-order
# spatial accuracy thanks to the true Dirichlet donor BC.
#
# Each block prints the achieved error so the run output stays informative.

# ---- helpers ---------------------------------------------------------------

# Drive the simulator with a thin, "infinite-dose" donor (finite_dose = FALSE)
# so the donor-skin interface is a true Dirichlet boundary at C0. The donor's
# D is set fast relative to the membrane only as a sanity precaution -- the
# Dirichlet BC clamps every donor cell, so the value is largely irrelevant.
run_dirichlet_skin <- function(layers_df, C0_mg_per_ml,
                               sim_time_min, scaling = "ng",
                               resolution = 4L,
                               donor_D = 1e4,
                               sink_Vd_ml = 1.0) {
  skin_simulate(skin_params(
    vehicle = list(
      c_init = C0_mg_per_ml, height = 3L, D = donor_D, app_area = 1.0,
      finite_dose = FALSE, log_mass = TRUE, log_cdp = FALSE
    ),
    layers = layers_df,
    sink = list(name = "Sink", Vd = sink_Vd_ml, log_mass = TRUE),
    pk = list(enabled = FALSE),
    sim_time = sim_time_min,
    resolution = resolution,
    scaling = scaling,
    max_module = 50
  ))
}

# Maximum and L2 relative error of `actual` against `expected`, ignoring
# samples where |expected| is below `floor` (avoids blow-ups near zero).
rel_errors <- function(actual, expected, floor = 1e-30) {
  keep <- abs(expected) > floor
  if (!any(keep)) return(c(linf = 0, l2 = 0))
  diffs <- abs(actual[keep] - expected[keep]) / abs(expected[keep])
  c(linf = max(diffs), l2 = sqrt(mean(diffs^2)))
}

report_errors <- function(label, errs, ...) {
  message(sprintf("[analytical:%s] L_inf=%.3e  L2=%.3e   %s",
                  label, errs["linf"], errs["l2"],
                  paste(..., collapse = "  ")))
}


# ---- single-slab lag-time --------------------------------------------------
# Crank eq. 4.24a. Membrane held at C0 on the donor side, perfect sink at the
# bottom; cumulative mass crossed approaches a straight line of slope
# C0*D/l with x-intercept l^2/(6*D) (the lag time).
test_that("single-slab lag-time matches Crank 4.24a", {
  C0 <- 1.0
  l  <- 100L
  D  <- 100.0
  app_area_cm2 <- 1.0

  layers <- data.frame(name = "Membrane", height = l,
                       D = D, K = 1.0, cross_section = 1.0,
                       log_mass = TRUE, log_cdp = FALSE)
  res <- run_dirichlet_skin(layers, C0, sim_time = 200L, scaling = "ng",
                            resolution = 4L)

  C0_internal <- C0 * 1e-12
  Q_per_area <- crank_single_slab_Q(res$mass$time, l = l, D = D,
                                    C_donor = C0_internal)
  expected_ng <- mg_per_um2_to_total(Q_per_area, app_area_cm2, "ng")

  t_lag <- crank_single_slab_lag_time(l, D)
  late  <- res$mass$time >= 4 * t_lag
  errs_late <- rel_errors(res$mass$Sink[late], expected_ng[late])

  fit <- stats::lm(res$mass$Sink[late] ~ res$mass$time[late])
  slope_sim <- unname(stats::coef(fit)[2])
  slope_an  <- crank_single_slab_slope(l, D, C0_internal) *
               (app_area_cm2 * 1e8) * 1e6
  slope_rel <- abs(slope_sim - slope_an) / abs(slope_an)
  intercept_sim <- -unname(stats::coef(fit)[1]) / slope_sim
  lag_rel <- abs(intercept_sim - t_lag) / t_lag

  report_errors("slab-lag-time", errs_late,
                sprintf("slope_rel=%.3e", slope_rel),
                sprintf("lag_rel=%.3e", lag_rel),
                sprintf("(t_lag=%.3f, sim_lag=%.3f)", t_lag, intercept_sim))

  expect_lt(errs_late["linf"], 1e-5)
  expect_lt(slope_rel,         1e-3)
  expect_lt(lag_rel,           1e-2)
})


# ---- two-layer steady-state K-jump (the headline test) ---------------------
# Crank section 12 / Ash & Barrer: piecewise-linear concentration profile
# with a partition jump at the L1/L2 interface. The activity-FVM scheme is
# exact for piecewise-linear u, so we expect machine-precision agreement.
test_that("two-layer steady-state K-jump", {
  C0 <- 1.0
  l1 <- 50L; l2 <- 50L
  D1 <- 100.0; D2 <- 1000.0
  K1 <- 5.0;   K2 <- 0.5
  app_area_cm2 <- 1.0
  sim_time <- 800L

  layers <- data.frame(
    name = c("L1", "L2"),
    height = c(l1, l2), D = c(D1, D2), K = c(K1, K2),
    cross_section = c(1.0, 1.0),
    log_mass = c(TRUE, TRUE), log_cdp = c(TRUE, TRUE)
  )
  res <- run_dirichlet_skin(layers, C0, sim_time, scaling = "ng",
                            resolution = 4L)

  C0_int <- C0 * 1e-12
  ss <- crank_two_layer_steady(l1, D1, K1, l2, D2, K2, C0_int)

  cdp1 <- res$cdp$L1; cdp2 <- res$cdp$L2
  conc1_final <- cdp1$conc[, ncol(cdp1$conc)]
  conc2_final <- cdp2$conc[, ncol(cdp2$conc)]
  exp1 <- ss$profile(cdp1$depth_um)            * 1e18
  exp2 <- ss$profile(cdp2$depth_um + l1)       * 1e18

  l2_floor <- 0.05 * exp2[1]
  keep2 <- exp2 > l2_floor

  errs1 <- rel_errors(conc1_final, exp1)
  errs2 <- rel_errors(conc2_final[keep2], exp2[keep2])
  report_errors("ss-two-layer-L1", errs1)
  report_errors("ss-two-layer-L2", errs2,
                sprintf("(n_pts=%d / %d)", sum(keep2), length(keep2)))

  late <- res$mass$time >= 0.5 * sim_time
  fit <- stats::lm(res$mass$Sink[late] ~ res$mass$time[late])
  slope_sim <- unname(stats::coef(fit)[2])
  slope_an  <- ss$flux * (app_area_cm2 * 1e8) * 1e6
  slope_rel <- abs(slope_sim - slope_an) / abs(slope_an)
  message(sprintf("[analytical:ss-two-layer] slope_rel=%.3e", slope_rel))

  # Steady-state piecewise-linear profile in u is recovered exactly by the
  # symmetric tridiagonal -- machine precision is the right floor here.
  expect_lt(errs1["linf"], 1e-10)
  expect_lt(errs2["linf"], 1e-10)
  expect_lt(slope_rel,     1e-10)
})


# ---- semi-infinite erfc profile -------------------------------------------
# Crank eq. 3.13. Surface clamped at C0, initial concentration zero, sample
# the profile while the diffusion front has not yet reached the bottom.
test_that("semi-infinite erfc profile matches Crank 3.13", {
  C0 <- 1.0
  L  <- 200L
  D  <- 100.0
  t_check <- 20

  layers <- data.frame(name = "Slab", height = L,
                       D = D, K = 1.0, cross_section = 1.0,
                       log_mass = FALSE, log_cdp = TRUE)
  res <- run_dirichlet_skin(layers, C0, sim_time = 30L, scaling = "ng",
                            resolution = 4L)

  cdp <- res$cdp$Slab
  ti <- which.min(abs(cdp$time - t_check))
  conc_sim <- cdp$conc[, ti]
  exp_ng_per_ml <- crank_semi_infinite(cdp$depth_um, cdp$time[ti], D,
                                       C0 * 1e-12) * 1e18
  # Compare only where the analytical solution is comfortably above the
  # noise floor (the slab isn't truly semi-infinite, so the tail accumulates
  # finite-domain error).
  meaningful <- exp_ng_per_ml > 1e-2 * C0 * 1e6
  errs <- rel_errors(conc_sim[meaningful], exp_ng_per_ml[meaningful])
  report_errors("erfc-short-time", errs,
                sprintf("(t=%g, n=%d / %d)",
                        cdp$time[ti], sum(meaningful), length(meaningful)))

  expect_lt(errs["linf"], 0.02)
})


# ---- spatial convergence rate ----------------------------------------------
# Run the lag-time problem at a sequence of mesh refinements (uniform mesh
# for clean rates) and fit a power law to the absolute error. With a true
# Dirichlet donor BC, the FVM scheme should give clean second-order
# convergence: the absolute error halves at every dx -> dx/2 refinement.
test_that("spatial convergence rate is roughly second order", {
  C0 <- 1.0
  l  <- 100L
  D  <- 100.0
  app_area_cm2 <- 1.0
  sim_time <- 200L

  layers <- data.frame(name = "Membrane", height = l,
                       D = D, K = 1.0, cross_section = 1.0,
                       log_mass = TRUE, log_cdp = FALSE)

  resolutions <- c(1L, 2L, 4L, 8L)
  C0_internal <- C0 * 1e-12
  errs <- numeric(length(resolutions))
  for (k in seq_along(resolutions)) {
    res <- run_dirichlet_skin(layers, C0, sim_time, scaling = "ng",
                              resolution = resolutions[k])
    Q_per_area <- crank_single_slab_Q(res$mass$time, l = l, D = D,
                                      C_donor = C0_internal)
    expected_ng <- mg_per_um2_to_total(Q_per_area, app_area_cm2, "ng")
    late <- res$mass$time >= 4 * crank_single_slab_lag_time(l, D)
    errs[k] <- max(abs(res$mass$Sink[late] - expected_ng[late]))
  }
  dx <- 1 / resolutions
  fit <- stats::lm(log(errs) ~ log(dx))
  rate <- unname(stats::coef(fit)[2])
  Q_an_max <- max(crank_single_slab_Q(c(sim_time), l = l, D = D,
                                      C_donor = C0_internal) *
                  app_area_cm2 * 1e8 * 1e6)
  rel_err_max <- max(errs) / Q_an_max
  message(sprintf("[analytical:convergence] rate=%.3f  abs_err=[%s]  rel_err_max=%.3e",
                  rate, paste(sprintf("%.3e", errs), collapse = ", "), rel_err_max))

  expect_gt(rate, 1.9)
  expect_lt(rate, 2.1)
  expect_lt(rel_err_max, 1e-5)
})
