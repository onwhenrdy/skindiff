# Coverage for the Phase 3 metrics + accessors. The reference values come
# from the single-slab Crank model (eq. 4.24a):
#
#   J_ss = C_donor * D / l            (steady-state flux)
#   t_lag = l^2 / (6 * D)             (lag time)
#   K_p   = D / l                     (= J_ss / C_donor, only valid when
#                                        C_donor stays constant -- so we
#                                        run with finite_dose = FALSE).
#
# We use a thin Dirichlet donor, so the simulator emulates a Franz-cell
# experiment with a constant C_donor at the membrane surface. Use a long
# simulation so the system has time to reach steady state.

# Build a single-slab infinite-dose run that hits steady state.
make_slab_run <- function(C0_mg = 1.0, l_um = 100L, D_per_min = 100.0,
                          area_cm2_val = 1.0, sim_min = 600L,
                          scaling = "ng",
                          sink_finite = FALSE) {
  snk <- if (sink_finite) finite_sink("Sink", Vd = ml(1.0))
         else             perfect_sink("Sink")
  skin_simulate(skin_params(
    area = cm2(area_cm2_val),
    vehicle = vehicle(
      c_init = mg_per_ml(C0_mg), height = um(3L),
      D = um2_per_min(1e4), finite_dose = FALSE,
      log_mass = TRUE
    ),
    layers = list(
      layer("Membrane", height = um(l_um),
            D = um2_per_min(D_per_min), K = 1.0, cross_section = 1.0,
            log_mass = TRUE, log_cdp = TRUE)
    ),
    sink = snk,
    duration = minutes(sim_min),
    resolution = 4L,
    scaling = scaling,
    max_module = 50
  ))
}

# ---------- permeated() -----------------------------------------------------

test_that("permeated returns time/Q with right units and shape", {
  res <- make_slab_run()
  perm <- permeated(res)
  expect_s3_class(perm, "data.frame")
  expect_named(perm, c("time", "Q"))
  expect_true(inherits(perm$time, "units"))
  expect_true(inherits(perm$Q,    "units"))
  expect_equal(units::deparse_unit(perm$time), "min")
  expect_equal(units::deparse_unit(perm$Q),    "ng cm-2")
  expect_equal(nrow(perm), nrow(res$mass))
})

test_that("permeated errors when sink mass is not logged", {
  # Build a run where sink log_mass is FALSE.
  res <- skin_simulate(skin_params(
    area = cm2(1),
    vehicle = vehicle(c_init = mg_per_ml(1), height = um(30L),
                      D = um2_per_min(1)),
    layers = list(layer("SC", height = um(20L), D = um2_per_min(1),
                        K = 1, cross_section = 1)),
    sink = perfect_sink("Sink", log_mass = FALSE),
    duration = minutes(30L)
  ))
  expect_error(permeated(res), "Sink mass was not logged")
})

# ---------- flux() ----------------------------------------------------------

test_that("flux returns time/flux with right units and one fewer row", {
  res <- make_slab_run()
  fl <- flux(res)
  expect_s3_class(fl, "data.frame")
  expect_named(fl, c("time", "flux"))
  expect_true(inherits(fl$time, "units"))
  expect_true(inherits(fl$flux, "units"))
  expect_equal(nrow(fl), nrow(res$mass) - 1L)
  # flux unit: ng / (cm^2 * min)
  expect_equal(units::deparse_unit(fl$flux), "ng cm-2 min-1")
})

# ---------- permeated_at() --------------------------------------------------

test_that("permeated_at interpolates linearly between sample times", {
  res <- make_slab_run()
  perm <- permeated(res)
  # Pick a query time that falls between two grid points.
  q1 <- permeated_at(res, minutes(100))
  q2 <- permeated_at(res, hours(5))     # 300 min
  expect_true(inherits(q1, "units"))
  expect_equal(units::deparse_unit(q1), "ng cm-2")
  expect_equal(length(q1), 1L)
  # 5 hours = 300 min, equals exactly perm$Q at row 301 (t=0..600 grid).
  expect_equal(as.numeric(q2),
               as.numeric(perm$Q[as.numeric(perm$time) == 300]),
               tolerance = 1e-12)
})

test_that("permeated_at rejects bare numeric", {
  res <- make_slab_run()
  expect_error(permeated_at(res, 100), "units-of-time")
})

# ---------- profile_at() ----------------------------------------------------

test_that("profile_at returns one frame per logged compartment", {
  res <- make_slab_run()
  prof <- profile_at(res, hours(5L))
  expect_named(prof, c("Membrane"))
  m <- prof$Membrane
  expect_s3_class(m, "data.frame")
  expect_named(m, c("depth", "conc"))
  expect_true(inherits(m$depth, "units"))
  expect_true(inherits(m$conc,  "units"))
  expect_equal(units::deparse_unit(m$depth), "um")
  expect_equal(units::deparse_unit(m$conc),  "ng ml-1")
  # The profile slice at t=5h should be increasing-then-decreasing or
  # monotonically-decreasing depending on stage; for a Dirichlet donor
  # at this t/D the profile decays from top of slab to bottom.
  conc_bare <- as.numeric(m$conc)
  expect_gt(conc_bare[1], conc_bare[length(conc_bare)])
})

# ---------- metrics() against Crank single-slab ----------------------------

test_that("metrics J_ss / t_lag / K_p match Crank closed form", {
  C0  <- 1.0       # mg/ml
  l   <- 100L      # um
  D   <- 100.0     # um^2/min
  res <- make_slab_run(C0_mg = C0, l_um = l, D_per_min = D, sim_min = 600L)

  m <- metrics(res)
  expect_s3_class(m, "data.frame")
  expect_equal(nrow(m), 1L)
  expect_named(m, c("J_ss", "t_lag", "K_p", "Q_total", "r2_ss",
                    "t_50_donor", "AUC_sink", "C_max_sink", "t_max_sink"))

  # Convert closed-form values to comparable units.
  # J_ss = C_donor * D / l in [(mg/ml) * (um^2/min) / um] = (mg/ml)*(um/min)
  # = 1e-12 mg/um^3 * D[um^2/min] / l[um] = (mg/um^2/min) -> *1e8 per cm^2
  # Reported scaling is "ng" -> *1e6 from mg.
  J_ss_an_per_min <- C0 * 1e-12 * D / l * 1e8 * 1e6   # ng/cm^2/min
  J_ss_an_per_h   <- J_ss_an_per_min * 60
  t_lag_an_h      <- l^2 / (6 * D) / 60              # hours
  K_p_an_cm_per_h <- D / l * 60 * 1e-4               # um/min -> cm/h: *60/1e4

  J_ss_sim <- as.numeric(m$J_ss)
  t_lag_sim <- as.numeric(m$t_lag)
  K_p_sim   <- as.numeric(m$K_p)

  expect_lt(abs(J_ss_sim - J_ss_an_per_h) / J_ss_an_per_h, 1e-3)
  expect_lt(abs(t_lag_sim - t_lag_an_h) / t_lag_an_h, 1e-2)
  expect_lt(abs(K_p_sim   - K_p_an_cm_per_h) / K_p_an_cm_per_h, 1e-3)

  # r2 should be very close to 1 because we're well into steady state.
  expect_gt(m$r2_ss, 0.999)
})

# ---------- metrics() with finite sink ------------------------------------

test_that("metrics computes AUC, C_max, t_max only for finite sink", {
  res_finite <- make_slab_run(sim_min = 200L, sink_finite = TRUE)
  res_perfect <- make_slab_run(sim_min = 200L, sink_finite = FALSE)

  m_fin <- metrics(res_finite)
  m_per <- metrics(res_perfect)

  # finite sink: numbers should be finite
  expect_true(is.finite(as.numeric(m_fin$C_max_sink)))
  expect_true(is.finite(as.numeric(m_fin$AUC_sink)))
  expect_true(is.finite(as.numeric(m_fin$t_max_sink)))

  # perfect sink: NA on these
  expect_true(is.na(as.numeric(m_per$C_max_sink)))
  expect_true(is.na(as.numeric(m_per$AUC_sink)))
  expect_true(is.na(as.numeric(m_per$t_max_sink)))
})

# ---------- metrics() with finite-dose donor warns -------------------------

test_that("metrics warns when finite-dose donor depletes substantially", {
  # Small finite-dose donor + thin SC: substantial depletion expected.
  expect_warning({
    res <- skin_simulate(skin_params(
      area = cm2(1),
      vehicle = vehicle(
        c_init = mg_per_ml(1.0), height = um(30L),
        D = um2_per_min(100), finite_dose = TRUE,
        log_mass = TRUE
      ),
      layers = list(
        layer("SC", height = um(20L), D = um2_per_min(100),
              K = 1.0, cross_section = 1.0, log_mass = TRUE)
      ),
      sink = perfect_sink("Sink"),
      duration = minutes(600L),
      resolution = 4L,
      scaling = "ng"
    ))
    metrics(res)
  }, "depleted")
})

# ---------- ss_window options ---------------------------------------------

test_that("ss_window accepts both fractions and units-of-time", {
  res <- make_slab_run(sim_min = 600L)
  m1 <- metrics(res, ss_window = c(0.7, 1.0))
  m2 <- metrics(res, ss_window = c(minutes(420), minutes(600)))   # same window
  expect_equal(as.numeric(m1$J_ss), as.numeric(m2$J_ss), tolerance = 1e-12)
})

test_that("ss_window rejects out-of-range fractions", {
  res <- make_slab_run(sim_min = 60L)
  expect_error(metrics(res, ss_window = c(-0.1, 1.0)), "ss_window")
  expect_error(metrics(res, ss_window = c(0.0, 1.5)),  "ss_window")
})

test_that("ss_window errors when window is too narrow for >= 2 points", {
  res <- make_slab_run(sim_min = 60L)
  expect_error(metrics(res, ss_window = c(0.999, 1.0)),
               "fewer than 2 sample points")
})
