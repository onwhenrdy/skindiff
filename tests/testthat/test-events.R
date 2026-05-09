# Coverage for the donor-replacement and donor-removal events. Both events
# come from real-world use cases (the original example.R simulates 18 days
# with donor replaced every 24h and removed at day 16). They run for both
# matrix methods.
#
# Behaviour confirmed:
#  - replace_after resets every donor cell to c_init at t = k * replace_after.
#  - remove_at deletes the donor compartment from the simulation; subsequent
#    timesteps run on (skin layers + sink). The donor's recorded mass series
#    is not surfaced in the result after removal -- this is a known
#    limitation of the current logger indexing.

methods <- c("DSkin_1_4", "Activity_FVM")

base_params <- function(matrix_method, ...) {
  defaults <- list(
    vehicle = list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0,
                   finite_dose = TRUE, log_mass = TRUE),
    layers = data.frame(name = "SC", height = 20L, D = 1.0, K = 1.0,
                        cross_section = 1.0, log_mass = TRUE),
    sink = list(name = "Sink", Vd = 1.0, log_mass = TRUE),
    sim_time = 60L,
    matrix_method = matrix_method
  )
  args <- modifyList(defaults, list(...))
  do.call(skin_params, args)
}

test_that("replace_after resets the donor at every multiple of the period", {
  for (mm in methods) {
    p <- base_params(
      mm,
      vehicle = list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0,
                     finite_dose = TRUE, log_mass = TRUE,
                     replace_after = 30L)
    )
    res <- skin_simulate(p)
    initial <- res$mass$Vehicle[1]

    # Just before replace at t = 30: depleted.
    expect_lt(res$mass$Vehicle[30], initial)

    # Right at t = 30 (after the replace event): back to initial.
    expect_equal(res$mass$Vehicle[31], initial, tolerance = 1e-12)

    # And again at t = 60.
    expect_equal(res$mass$Vehicle[61], initial, tolerance = 1e-12)
  }
})

test_that("replace_after: mass is conserved between events, jumps at events", {
  # The replace event is an *unbounded source*: it tops the donor back up to
  # c_init, injecting mass (= initial - depleted) into the system. Therefore:
  #   - between events the system is closed and total mass is constant,
  #   - at each replace event, total mass jumps up by the refill amount.
  for (mm in methods) {
    p <- base_params(
      mm,
      vehicle = list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0,
                     finite_dose = TRUE, log_mass = TRUE,
                     replace_after = 30L)
    )
    res <- skin_simulate(p)
    totals <- rowSums(res$mass[, c("Vehicle", "SC", "Sink")])

    # Recorded rows 1..30 cover sim times 0..29 (no replace event yet).
    rel_drift_pre <- max(abs(totals[1:30] - totals[1])) / totals[1]
    expect_lt(rel_drift_pre, 1e-10)

    # Recorded rows 31..60 cover sim times 30..59, after the first replace
    # but before the second.
    rel_drift_mid <- max(abs(totals[31:60] - totals[31])) / totals[31]
    expect_lt(rel_drift_mid, 1e-10)

    # The replace events at t = 30 and t = 60 inject positive mass.
    expect_gt(totals[31] - totals[30], 0)
    expect_gt(totals[61] - totals[60], 0)

    # The size of the jump at the first replace equals the donor refill --
    # i.e. how much donor mass was missing just before the replace. Donor
    # mass at row 30 (sim t = 29) is ALMOST the pre-replace value: there's
    # a one-minute diffusion step between row 30 and the replace event at
    # row 31, but the donor lost only a tiny fraction in that minute. So
    # the jump equals (initial - depleted_at_t=30), and depleted_at_t=30
    # is close to (but slightly less than) vehicle[30].
    initial <- res$mass$Vehicle[1]
    refill  <- totals[31] - totals[30]
    naive   <- initial - res$mass$Vehicle[30]   # ignores 29 -> 30 diffusion
    expect_lt(abs(refill - naive) / naive, 0.05)  # within ~few % of naive
  }
})

test_that("remove_at deletes the donor and the sim continues", {
  for (mm in methods) {
    p <- base_params(
      mm,
      vehicle = list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0,
                     finite_dose = TRUE, log_mass = TRUE,
                     remove_at = 30L)
    )
    res <- skin_simulate(p)

    # Donor's mass series isn't in the post-removal output (logger
    # indexing limitation -- see test header). Skin and sink continue.
    expect_false("Vehicle" %in% names(res$mass))
    expect_true("SC" %in% names(res$mass))
    expect_true("Sink" %in% names(res$mass))
    expect_equal(nrow(res$mass), 61L)

    # Sink continues to accumulate from the SC reservoir even after the
    # donor is gone.
    sink_at_30 <- res$mass$Sink[31]
    sink_at_60 <- res$mass$Sink[61]
    expect_gt(sink_at_60, sink_at_30)

    # SC keeps releasing mass into the sink (monotonic decrease, modulo
    # reflection at the surface boundary that's now in place).
    sc_at_30 <- res$mass$SC[31]
    sc_at_60 <- res$mass$SC[61]
    expect_lt(sc_at_60, sc_at_30)
  }
})

test_that("combined replace + remove reproduces example.R-style setup", {
  # The original example.R: replace every 30 min, remove at 90 min, sim to
  # 120 min. This exercises multiple replace events plus a removal and
  # confirms both methods complete cleanly.
  for (mm in methods) {
    p <- base_params(
      mm,
      sim_time = 120L,
      vehicle = list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0,
                     finite_dose = TRUE, log_mass = TRUE,
                     replace_after = 30L,
                     remove_at = 90L)
    )
    res <- skin_simulate(p)

    expect_equal(res$status, "executed")
    # Donor was replaced at t = 30, 60 (before removal at 90); column gets
    # dropped on removal.
    expect_false("Vehicle" %in% names(res$mass))
    expect_equal(nrow(res$mass), 121L)

    # Sink should hold meaningful (positive) mass at the end.
    expect_gt(tail(res$mass$Sink, 1), 0)
  }
})
