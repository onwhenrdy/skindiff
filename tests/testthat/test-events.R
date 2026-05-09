# Coverage for the donor-replacement and donor-removal events. Both events
# come from real-world use cases (the original example.R simulates 18 days
# with donor replaced every 24h and removed at day 16).
#
# Behaviour confirmed:
#  - replace_after resets every donor cell to c_init at t = k * replace_after.
#  - remove_at deletes the donor compartment from the simulation; subsequent
#    timesteps run on (skin layers + sink). The donor's recorded mass series
#    is not surfaced in the result after removal -- a known limitation of
#    the current logger indexing.

base_params <- function(...) {
  defaults <- list(
    vehicle = list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0,
                   finite_dose = TRUE, log_mass = TRUE),
    layers = data.frame(name = "SC", height = 20L, D = 1.0, K = 1.0,
                        cross_section = 1.0, log_mass = TRUE),
    sink = list(name = "Sink", Vd = 1.0, log_mass = TRUE),
    sim_time = 60L
  )
  args <- modifyList(defaults, list(...))
  do.call(skin_params, args)
}

test_that("replace_after resets the donor at every multiple of the period", {
  p <- base_params(
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
})

test_that("replace_after: mass is conserved between events, jumps at events", {
  # The replace event is an *unbounded source*: it tops the donor back up to
  # c_init, injecting mass (= initial - depleted) into the system. Therefore:
  #   - between events the system is closed and total mass is constant,
  #   - at each replace event, total mass jumps up by the refill amount.
  p <- base_params(
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
  # how much donor mass was missing just before the replace. Donor mass at
  # row 30 (sim t = 29) is almost the pre-replace value (the donor only
  # loses a tiny amount during the 29 -> 30 diffusion sub-step), so the
  # jump matches `initial - vehicle[30]` to within a few percent.
  initial <- res$mass$Vehicle[1]
  refill  <- totals[31] - totals[30]
  naive   <- initial - res$mass$Vehicle[30]
  expect_lt(abs(refill - naive) / naive, 0.05)
})

test_that("remove_at deletes the donor and the sim continues", {
  p <- base_params(
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
  expect_gt(res$mass$Sink[61], res$mass$Sink[31])
  # SC keeps releasing mass into the sink.
  expect_lt(res$mass$SC[61], res$mass$SC[31])
})

test_that("combined replace + remove reproduces example.R-style setup", {
  # Replace every 30 min, remove at 90 min, sim to 120 min. Exercises
  # multiple replace events plus a removal in one run.
  p <- base_params(
    sim_time = 120L,
    vehicle = list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0,
                   finite_dose = TRUE, log_mass = TRUE,
                   replace_after = 30L,
                   remove_at = 90L)
  )
  res <- skin_simulate(p)

  expect_equal(res$status, "executed")
  expect_false("Vehicle" %in% names(res$mass))
  expect_equal(nrow(res$mass), 121L)
  expect_gt(tail(res$mass$Sink, 1), 0)
})
