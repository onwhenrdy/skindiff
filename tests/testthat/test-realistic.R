# End-to-end test on the realistic SC / DSL / PK setup that originally
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
    vehicle = list(
      c_init   = 127.2727, height = 110L, D = 9.266667, app_area = 15.0,
      finite_dose = TRUE, log_mass = TRUE, log_cdp = TRUE
    ),
    layers = data.frame(
      name          = c("Stratum corneum", "Deeper skin layers"),
      height        = c(190L, 200L),
      D             = c(28.2539, 5767.783),
      K             = c(421.543, 0.04719648),
      cross_section = c(0.001, 0.3),
      log_mass = c(TRUE, TRUE),
      log_cdp  = c(TRUE, TRUE)
    ),
    sink = list(name = "Blood", Vd = 25 * 1000 * 75, log_mass = TRUE),
    pk   = list(enabled = TRUE, thalf = 24),
    sim_time      = 24 * 60L,
    resolution    = 4L,
    scaling       = "ng",
    max_module    = 200
  ))
}

activity_jump <- function(res) {
  K_sc  <- 421.543
  K_dsl <- 0.04719648
  sc  <- res$cdp$`Stratum corneum`$conc
  dsl <- res$cdp$`Deeper skin layers`$conc
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

test_that("realistic SC/DSL/PK example: outputs sane and activity continuous", {
  res <- run_realistic()

  # Status / finiteness / order.
  expect_equal(res$status, "executed")
  expect_true(all(is.finite(res$mass$Blood)))
  expect_true(all(is.finite(res$mass$Vehicle)))
  expect_true(all(res$mass$Blood   >= -1e-10))
  expect_true(all(res$mass$Vehicle >= -1e-10))
  # No replace event here, so vehicle mass is monotonically non-increasing.
  expect_true(all(diff(res$mass$Vehicle) <= 1e-6 * max(res$mass$Vehicle)))

  # Activity continuity at the SC/DSL K-jump. Baseline ~2e-2 on this stack
  # and resolution; the jump is largely a cell-center sampling artefact
  # (cells closest to the interface sit at finite distance from it).
  jmp <- activity_jump(res)
  message(sprintf("[realistic:activity] u_sc_bot=%.4e  u_dsl_top=%.4e  rel_jump=%.3e",
                  jmp$u_sc_bot, jmp$u_dsl_top, jmp$rel_jump))
  expect_lt(jmp$rel_jump, 0.05)
})
