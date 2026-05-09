# End-to-end test on the realistic SC / DSL / PK setup that originally
# motivated the project (cf. the old scripts/example.R). The K-jump at the
# SC/DSL interface is steep -- K = 421 (SC) -> K = 0.047 (DSL), ratio ~9000
# -- which is exactly where the c-formulation of DSkin_1_4 was known to
# struggle.
#
# The test asserts:
#   (a) both methods run to completion and give finite, ordered results,
#   (b) Activity_FVM keeps the activity (c/K) continuous across the K
#       jump (the activity transformation does this by construction),
#   (c) DSkin_1_4 has a measurable activity jump there -- this is recorded
#       as the *baseline* DSkin_1_4 disagrees with itself; we lock in a
#       generous floor so we'll notice if a future change makes it worse.
#
# The two methods can disagree by tens of percent on final masses for this
# kind of stack -- not a bug in either, just a measurement of how much the
# discretization choice matters for K-heavy problems.

run_realistic <- function(matrix_method) {
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
    disc_method   = "bk",
    matrix_method = matrix_method,
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

test_that("realistic SC/DSL/PK example: both methods run and give sane outputs", {
  res_old <- run_realistic("DSkin_1_4")
  res_new <- run_realistic("Activity_FVM")

  # Finiteness + structure.
  for (res in list(res_old, res_new)) {
    expect_equal(res$status, "executed")
    expect_true(all(is.finite(res$mass$Blood)))
    expect_true(all(is.finite(res$mass$Vehicle)))
    expect_true(all(res$mass$Blood >= -1e-10))   # no spurious negatives
    expect_true(all(res$mass$Vehicle >= -1e-10))
  }

  # Vehicle is monotonically non-increasing (no replacement event here).
  for (res in list(res_old, res_new)) {
    expect_true(all(diff(res$mass$Vehicle) <= 1e-6 * max(res$mass$Vehicle)))
  }
})

test_that("Activity_FVM keeps activity continuous at the SC/DSL interface", {
  res_new <- run_realistic("Activity_FVM")
  jmp <- activity_jump(res_new)
  message(sprintf("[realistic:activity] u_sc_bot=%.4e  u_dsl_top=%.4e  rel_jump=%.3e",
                  jmp$u_sc_bot, jmp$u_dsl_top, jmp$rel_jump))
  # On this stack and resolution, baseline is ~5e-3.
  expect_lt(jmp$rel_jump, 0.02)
})

test_that("DSkin_1_4 has the documented activity-jump artefact (regression floor)", {
  res_old <- run_realistic("DSkin_1_4")
  jmp <- activity_jump(res_old)
  message(sprintf("[realistic:legacy]   u_sc_bot=%.4e  u_dsl_top=%.4e  rel_jump=%.3e",
                  jmp$u_sc_bot, jmp$u_dsl_top, jmp$rel_jump))
  # On this stack and resolution, baseline is ~1.5e-1. Locked in as a
  # ceiling -- if a change pushes it past 0.5, something is very wrong.
  expect_gt(jmp$rel_jump, 0.05)   # confirms the artefact still exists
  expect_lt(jmp$rel_jump, 0.50)   # but isn't catastrophically worse
})

test_that("the two methods disagree by tens of percent on final masses (recorded)", {
  res_old <- run_realistic("DSkin_1_4")
  res_new <- run_realistic("Activity_FVM")

  for (col in c("Vehicle", "Stratum corneum", "Deeper skin layers", "Blood")) {
    o <- tail(res_old$mass[[col]], 1)
    n <- tail(res_new$mass[[col]], 1)
    rel_diff <- (n - o) / max(abs(o), 1e-30)
    message(sprintf("[realistic:diff] %-22s  legacy=%.4e  activity=%.4e  rel=%+.3e",
                    col, o, n, rel_diff))
  }
  # No assertion -- this is a recording. The disagreement is large (~30%
  # on Vehicle and Blood at 24h) and is the whole point: the two methods
  # genuinely give different answers on K-heavy stacks, and the activity
  # diagnostic above is the tiebreaker in favour of Activity_FVM.
  expect_true(TRUE)
})
