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
#
# `layers_df` is a data.frame with one row per layer (columns: name, height,
# D, K, cross_section, optional log_mass, log_cdp); we translate it to a list
# of layer() objects here.
run_dirichlet_skin <- function(layers_df, C0_mg_per_ml,
                               sim_time_min, scaling = "ng",
                               resolution = 4L,
                               donor_D = 1e4,
                               sink_Vd_ml = NULL) {
  layers <- lapply(seq_len(nrow(layers_df)), function(i) {
    args <- list(
      name          = as.character(layers_df$name[[i]]),
      height        = um(as.integer(layers_df$height[[i]])),
      D             = um2_per_min(as.numeric(layers_df$D[[i]])),
      K             = as.numeric(layers_df$K[[i]]),
      cross_section = as.numeric(layers_df$cross_section[[i]])
    )
    if (!is.null(layers_df$log_mass)) args$log_mass <- as.logical(layers_df$log_mass[[i]])
    if (!is.null(layers_df$log_cdp))  args$log_cdp  <- as.logical(layers_df$log_cdp[[i]])
    do.call(layer, args)
  })

  snk <- if (is.null(sink_Vd_ml)) perfect_sink("Sink")
         else                     finite_sink("Sink", Vd = ml(sink_Vd_ml))

  skin_simulate(skin_params(
    area = cm2(1.0),
    vehicle = vehicle(
      c_init = mg_per_ml(C0_mg_per_ml), height = um(3L),
      D = um2_per_min(donor_D), finite_dose = FALSE,
      log_mass = TRUE, log_cdp = FALSE
    ),
    layers = layers,
    sink = snk,
    duration = minutes(sim_time_min),
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
  t_min <- as.numeric(res$mass$time)
  m_sink <- as.numeric(res$mass$Sink)
  Q_per_area <- crank_single_slab_Q(t_min, l = l, D = D,
                                    C_donor = C0_internal)
  expected_ng <- mg_per_um2_to_total(Q_per_area, app_area_cm2, "ng")

  t_lag <- crank_single_slab_lag_time(l, D)
  late  <- t_min >= 4 * t_lag
  errs_late <- rel_errors(m_sink[late], expected_ng[late])

  fit <- stats::lm(m_sink[late] ~ t_min[late])
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
  conc1_final <- as.numeric(cdp1$conc[, ncol(cdp1$conc)])
  conc2_final <- as.numeric(cdp2$conc[, ncol(cdp2$conc)])
  depth1 <- as.numeric(cdp1$depth)
  depth2 <- as.numeric(cdp2$depth)
  exp1 <- ss$profile(depth1)            * 1e18
  exp2 <- ss$profile(depth2 + l1)       * 1e18

  l2_floor <- 0.05 * exp2[1]
  keep2 <- exp2 > l2_floor

  errs1 <- rel_errors(conc1_final, exp1)
  errs2 <- rel_errors(conc2_final[keep2], exp2[keep2])
  report_errors("ss-two-layer-L1", errs1)
  report_errors("ss-two-layer-L2", errs2,
                sprintf("(n_pts=%d / %d)", sum(keep2), length(keep2)))

  t_min <- as.numeric(res$mass$time)
  m_sink <- as.numeric(res$mass$Sink)
  late <- t_min >= 0.5 * sim_time
  fit <- stats::lm(m_sink[late] ~ t_min[late])
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


# ---- semi-infinite erfc profile, sweep across times ------------------------
# Crank eq. 3.13. Surface clamped at C0, initial concentration zero, sample
# the profile while the diffusion front has not yet reached the bottom.
# Sweep multiple times so the test exercises early-, mid-, and late-stage
# transient profiles -- not just one snapshot.
test_that("semi-infinite erfc profile matches Crank 3.13 across multiple times", {
  C0 <- 1.0
  L  <- 200L
  D  <- 100.0
  test_times <- c(5, 10, 15, 20)   # min; stop before the back-boundary
                                   # correction kills the semi-infinite
                                   # approximation (sqrt(D*t) approaches L)

  layers <- data.frame(name = "Slab", height = L,
                       D = D, K = 1.0, cross_section = 1.0,
                       log_mass = FALSE, log_cdp = TRUE)
  res <- run_dirichlet_skin(layers, C0, sim_time = 30L, scaling = "ng",
                            resolution = 4L)

  cdp <- res$cdp$Slab
  t_grid   <- as.numeric(cdp$time)
  depth_um <- as.numeric(cdp$depth)

  worst_linf <- 0
  for (t_check in test_times) {
    ti <- which.min(abs(t_grid - t_check))
    conc_sim <- as.numeric(cdp$conc[, ti])
    exp_ng_per_ml <- crank_semi_infinite(depth_um, t_grid[ti], D,
                                         C0 * 1e-12) * 1e18
    # Restrict to where the semi-infinite assumption holds: analytical
    # value above 1% of C0 (filters out the tail where finite-domain
    # error dominates).
    meaningful <- exp_ng_per_ml > 1e-2 * C0 * 1e6
    errs <- rel_errors(conc_sim[meaningful], exp_ng_per_ml[meaningful])
    report_errors(sprintf("erfc t=%g", t_grid[ti]), errs,
                  sprintf("(n=%d / %d)", sum(meaningful), length(meaningful)))
    worst_linf <- max(worst_linf, errs["linf"])
  }
  expect_lt(worst_linf, 0.02)
})


# ---- transient single-slab CDP, full closed-form ---------------------------
# Crank eq. 4.16. Same geometry as the lag-time test (Dirichlet donor,
# perfect sink, finite slab) but compares the *concentration profile*
# c(x, t) at multiple times, including post-breakthrough where the erfc
# approximation fails.
test_that("transient single-slab CDP matches Crank 4.16 across multiple times", {
  C0 <- 1.0
  l  <- 100L
  D  <- 100.0
  test_times <- c(2, 10, 30, 100, 200)   # span the transient

  layers <- data.frame(name = "Membrane", height = l,
                       D = D, K = 1.0, cross_section = 1.0,
                       log_mass = FALSE, log_cdp = TRUE)
  res <- run_dirichlet_skin(layers, C0, sim_time = 200L, scaling = "ng",
                            resolution = 4L)

  cdp <- res$cdp$Membrane
  t_grid   <- as.numeric(cdp$time)
  depth_um <- as.numeric(cdp$depth)
  C0_internal <- C0 * 1e-12

  worst_linf <- 0
  for (t_check in test_times) {
    ti <- which.min(abs(t_grid - t_check))
    conc_sim <- as.numeric(cdp$conc[, ti])
    exp_internal <- crank_single_slab_C(depth_um, t_grid[ti], l, D, C0_internal)
    exp_ng_per_ml <- exp_internal * 1e18
    keep <- exp_ng_per_ml > 1e-3 * C0 * 1e6
    if (!any(keep)) next
    errs <- rel_errors(conc_sim[keep], exp_ng_per_ml[keep])
    report_errors(sprintf("crank-cdp t=%g", t_grid[ti]), errs,
                  sprintf("(n=%d / %d)", sum(keep), length(keep)))
    worst_linf <- max(worst_linf, errs["linf"])
  }
  expect_lt(worst_linf, 0.02)
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
    t_min <- as.numeric(res$mass$time)
    m_sink <- as.numeric(res$mass$Sink)
    Q_per_area <- crank_single_slab_Q(t_min, l = l, D = D,
                                      C_donor = C0_internal)
    expected_ng <- mg_per_um2_to_total(Q_per_area, app_area_cm2, "ng")
    late <- t_min >= 4 * crank_single_slab_lag_time(l, D)
    errs[k] <- max(abs(m_sink[late] - expected_ng[late]))
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


# ---- Kasting (2001) finite-dose / well-mixed-vehicle ----------------------
# Kasting eq. 10. A vehicle of thickness h_v with concentration C_v0 sits on
# top of a homogeneous membrane (thickness h, diffusivity D, partition
# K_mv); the membrane releases into a perfect sink. The vehicle is assumed
# well-mixed (D_v -> infty), which we approximate by giving the donor a very
# large diffusion coefficient relative to the membrane.
#
# The analytical solution is in M_t / M_∞ (cumulative permeated fraction).
# β = K_mv · h / h_v is the dimensionless dose parameter; for β = 1 the
# system is in the strong finite-dose regime (clearly distinct from the
# infinite-dose Crank solution).
test_that("finite-dose cumulative absorption matches Kasting 2001 eq. 10", {
  C_v0  <- 1.0      # mg/ml
  h_v   <- 100L     # um
  h     <- 100L     # um
  D     <- 100.0    # um^2/min
  D_v   <- 1.0e6    # um^2/min: "well-mixed" approximation (h_v^2/D_v = 0.01 min)
  K_mv  <- 1.0
  area_cm2_val <- 1.0
  sim_time <- 600L  # min: D*t/h^2 = 6, well past full absorption
  beta <- K_mv * h / h_v

  res <- skin_simulate(skin_params(
    area = cm2(area_cm2_val),
    vehicle = vehicle(
      c_init = mg_per_ml(C_v0), height = um(h_v),
      D = um2_per_min(D_v),
      finite_dose = TRUE, log_mass = TRUE
    ),
    layers = list(
      layer("Membrane", height = um(h), D = um2_per_min(D),
            K = K_mv, cross_section = 1.0,
            log_mass = TRUE, log_cdp = FALSE)
    ),
    sink = perfect_sink("Sink"),
    duration = minutes(sim_time),
    resolution = 4L,
    scaling = "ng",
    max_module = 50
  ))

  # Total applied mass M_∞: vehicle volume * c_init.
  # h_v[um] * area[cm^2] * 1e-4 = volume[ml]; * c_init[mg/ml] = mass[mg];
  # then to ng = *1e6.
  M_inf_ng <- h_v * area_cm2_val * 1e-4 * C_v0 * 1e6

  t_min  <- as.numeric(res$mass$time)
  m_sink <- as.numeric(res$mass$Sink)
  M_t_over_Minf_sim <- m_sink / M_inf_ng

  M_t_over_Minf_an <- kasting_M_over_Minf(t_min, beta, D, h)

  # Compare from t > 0 (absolute mass error is tiny everywhere; relative
  # error blows up near t=0 because both numerator and denominator vanish.
  # Use a small floor on the analytical value to filter out the noise band).
  keep <- M_t_over_Minf_an >= 0.01
  errs <- rel_errors(M_t_over_Minf_sim[keep], M_t_over_Minf_an[keep])
  report_errors("kasting-finite-dose", errs,
                sprintf("(beta=%.2f, n=%d/%d)", beta, sum(keep), length(t_min)))

  expect_lt(errs["linf"], 5e-3)
  expect_lt(errs["l2"],   2e-3)
})


# ---- Kasting (2001) finite-dose membrane CDP ------------------------------
# Same BVP as the M_t test above. The c(x, t) form derived from the
# eigen-expansion is in helper-analytical.R::kasting_membrane_C. Here we
# log the membrane CDP and compare snapshot profiles at several times.
test_that("Kasting finite-dose membrane CDP matches eigen-expansion", {
  C_v0  <- 1.0
  h_v   <- 100L
  h     <- 100L
  D     <- 100.0
  D_v   <- 1.0e6
  K_mv  <- 1.0
  area_cm2_val <- 1.0
  sim_time <- 600L
  beta <- K_mv * h / h_v
  test_times <- c(10, 30, 100, 200, 400)   # spans peak-flux through tail

  res <- skin_simulate(skin_params(
    area = cm2(area_cm2_val),
    vehicle = vehicle(
      c_init = mg_per_ml(C_v0), height = um(h_v),
      D = um2_per_min(D_v),
      finite_dose = TRUE, log_mass = TRUE
    ),
    layers = list(
      layer("Membrane", height = um(h), D = um2_per_min(D),
            K = K_mv, cross_section = 1.0,
            log_mass = TRUE, log_cdp = TRUE)
    ),
    sink = perfect_sink("Sink"),
    duration = minutes(sim_time),
    resolution = 4L,
    scaling = "ng",
    max_module = 50
  ))

  cdp <- res$cdp$Membrane
  t_grid   <- as.numeric(cdp$time)
  depth_um <- as.numeric(cdp$depth)
  C_v0_int <- C_v0 * 1e-12

  worst_linf <- 0
  for (t_check in test_times) {
    ti <- which.min(abs(t_grid - t_check))
    conc_sim <- as.numeric(cdp$conc[, ti])
    exp_internal <- kasting_membrane_C(depth_um, t_grid[ti],
                                       beta, D, h, K_mv, C_v0_int)
    exp_ng_per_ml <- exp_internal * 1e18
    # Restrict to where the analytical value is comfortably above zero
    # (1e-3 of the K_mv·C_v0 reference scale).
    keep <- exp_ng_per_ml > 1e-3 * K_mv * C_v0 * 1e6
    if (!any(keep)) next
    errs <- rel_errors(conc_sim[keep], exp_ng_per_ml[keep])
    report_errors(sprintf("kasting-cdp t=%g", t_grid[ti]), errs,
                  sprintf("(n=%d / %d)", sum(keep), length(keep)))
    worst_linf <- max(worst_linf, errs["linf"])
  }
  expect_lt(worst_linf, 0.02)
})
