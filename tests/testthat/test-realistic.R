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

run_realistic <- function(disc_method = "bk") {
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
    disc_method   = disc_method,
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

test_that("realistic SC/DSL/PK example runs and gives sane outputs", {
  res <- run_realistic()
  expect_equal(res$status, "executed")
  expect_true(all(is.finite(res$mass$Blood)))
  expect_true(all(is.finite(res$mass$Vehicle)))
  expect_true(all(res$mass$Blood   >= -1e-10))
  expect_true(all(res$mass$Vehicle >= -1e-10))
  # No replace event here, so vehicle mass is monotonically non-increasing.
  expect_true(all(diff(res$mass$Vehicle) <= 1e-6 * max(res$mass$Vehicle)))
})

test_that("activity is continuous at the SC/DSL K-jump interface", {
  res <- run_realistic("bk")
  jmp <- activity_jump(res)
  message(sprintf("[realistic:activity] u_sc_bot=%.4e  u_dsl_top=%.4e  rel_jump=%.3e",
                  jmp$u_sc_bot, jmp$u_dsl_top, jmp$rel_jump))
  # On this stack and resolution, baseline is ~5e-3.
  expect_lt(jmp$rel_jump, 0.02)
})

test_that("graded mesh agrees with bk on the realistic SC/DSL/PK example", {
  # Graded uses far fewer cells in the high-D DSL and runs ~100x faster
  # than bk while giving essentially the same physics. The activity-jump
  # diagnostic at the SC/DSL interface is naturally larger with graded
  # (the cells closest to the interface are coarser so cell-center u
  # samples sit farther from the interface), but the cumulative blood
  # mass at the end of the run agrees with bk to ~0.01%.
  res_bk     <- run_realistic("bk")
  res_graded <- run_realistic("graded")

  blood_bk     <- tail(res_bk$mass$Blood, 1)
  blood_graded <- tail(res_graded$mass$Blood, 1)
  rel <- abs(blood_graded - blood_bk) / blood_bk
  message(sprintf("[realistic:graded] blood_bk=%.4e  blood_graded=%.4e  rel=%.3e",
                  blood_bk, blood_graded, rel))
  expect_lt(rel, 1e-3)

  # Activity-jump diagnostic on graded is larger (coarser cells at the
  # interface) but still well below the c-formulation's ~15% artefact.
  jmp_graded <- activity_jump(res_graded)
  message(sprintf("[realistic:graded-jump] rel_jump=%.3e", jmp_graded$rel_jump))
  expect_lt(jmp_graded$rel_jump, 0.05)
})
