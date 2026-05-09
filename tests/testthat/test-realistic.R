# End-to-end test on the realistic SC / DSL setup that originally
# motivated the project (cf. the old scripts/example.R). The K-jump at the
# SC/DSL interface is steep (K = 421 -> K = 0.047, ratio ~9000), which used
# to be a notorious source of discretization error in the old c-formulation
# scheme.
#
# The activity-FVM scheme keeps `u = c/K` continuous at the K-jump by
# construction. The test runs the realistic stack and checks:
#   - the run completes and produces finite, sensibly-ordered output,
#   - the activity-jump diagnostic at the SC/DSL interface stays small.

run_realistic <- function() {
  skin_simulate(skin_params(
    area = cm2(15.0),
    vehicle = vehicle(
      c_init = mg_per_ml(127.2727), height = um(110L),
      D = um2_per_min(9.266667), finite_dose = TRUE,
      log_mass = TRUE, log_cdp = TRUE
    ),
    layers = list(
      layer("Stratum corneum", height = um(190L), D = um2_per_min(28.2539),
            K = 421.543, cross_section = 0.001,
            log_mass = TRUE, log_cdp = TRUE),
      layer("Deeper skin layers", height = um(200L), D = um2_per_min(5767.783),
            K = 0.04719648, cross_section = 0.3,
            log_mass = TRUE, log_cdp = TRUE)
    ),
    sink = finite_sink("Blood", Vd = ml(25 * 1000 * 75), log_mass = TRUE),
    duration   = hours(24L),
    resolution = 4L,
    scaling    = "ng",
    max_module = 200
  ))
}

activity_jump <- function(res) {
  K_sc  <- 421.543
  K_dsl <- 0.04719648
  # Strip units for activity-jump diagnostics; we just want bare ratios.
  sc  <- as.matrix(unclass(res$cdp$`Stratum corneum`$conc))
  dsl <- as.matrix(unclass(res$cdp$`Deeper skin layers`$conc))
  last <- ncol(sc)
  c_sc_bot  <- sc[nrow(sc),  last]
  c_dsl_top <- dsl[1,        last]
  u_sc_bot  <- c_sc_bot  / K_sc
  u_dsl_top <- c_dsl_top / K_dsl
  list(
    u_sc_bot  = u_sc_bot,
    u_dsl_top = u_dsl_top,
    rel_jump  = abs(u_sc_bot - u_dsl_top) / max(u_sc_bot, 1e-30)
  )
}

test_that("realistic SC/DSL example: outputs sane and activity continuous", {
  res <- run_realistic()

  # Status / finiteness / order.
  expect_equal(res$status, "executed")
  expect_true(all(is.finite(res$mass$Blood)))
  expect_true(all(is.finite(res$mass$Vehicle)))
  blood_bare <- as.numeric(res$mass$Blood)
  vehicle_bare <- as.numeric(res$mass$Vehicle)
  expect_true(all(blood_bare   >= -1e-10))
  expect_true(all(vehicle_bare >= -1e-10))
  # No replace event here, so vehicle mass is monotonically non-increasing.
  expect_true(all(diff(vehicle_bare) <= 1e-6 * max(vehicle_bare)))

  # Activity continuity at the SC/DSL K-jump. Baseline ~2e-2 on this stack
  # and resolution; the jump is largely a cell-center sampling artefact
  # (cells closest to the interface sit at finite distance from it).
  jmp <- activity_jump(res)
  message(sprintf("[realistic:activity] u_sc_bot=%.4e  u_dsl_top=%.4e  rel_jump=%.3e",
                  jmp$u_sc_bot, jmp$u_dsl_top, jmp$rel_jump))
  expect_lt(jmp$rel_jump, 0.05)
})
