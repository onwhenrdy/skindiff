make_minimal <- function(vehicle = NULL, layers = NULL, sink = NULL,
                         pk = NULL, sim_time = NULL, ...) {
  v <- vehicle %||% list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0,
                         log_mass = TRUE, log_cdp = TRUE)
  l <- if (is.null(layers))
         data.frame(name = "SC", height = 20L, D = 1.0, K = 1.0,
                    cross_section = 1.0,
                    log_mass = TRUE, log_cdp = TRUE)
       else layers
  s <- sink %||% list(Vd = 1.0, log_mass = TRUE)
  pk_arg <- pk %||% list(enabled = FALSE, thalf = 1)
  st <- sim_time %||% 30L
  skin_params(vehicle = v, layers = l, sink = s, pk = pk_arg,
              sim_time = st, ...)
}

run_minimal <- function(...) skin_simulate(make_minimal(...))

`%||%` <- function(a, b) if (is.null(a)) b else a

test_that("skin_simulate returns the documented structure", {
  res <- run_minimal()
  expect_s3_class(res, "skin_result")
  expect_equal(res$status, "executed")
  expect_named(res, c("status", "scaling", "mass", "concentration",
                      "cdp", "geometry", "params", "runtime_s"))
  expect_s3_class(res$mass, "data.frame")
  expect_s3_class(res$concentration, "data.frame")
  expect_true(is.list(res$cdp))
  expect_true(is.list(res$geometry))
})

test_that("mass data.frame has the expected columns", {
  res <- run_minimal()
  expect_setequal(colnames(res$mass), c("time", "Vehicle", "SC", "Sink"))
  expect_setequal(colnames(res$concentration),
                  c("time", "Vehicle", "SC", "Sink"))
  expect_equal(nrow(res$mass), 31L)  # 0..30 inclusive
})

test_that("mass is conserved for a perfect-sink run (relative err < 1e-10)", {
  res <- run_minimal()
  totals <- rowSums(res$mass[, c("Vehicle", "SC", "Sink")])
  rel <- max(abs(totals - totals[1])) / totals[1]
  expect_lt(rel, 1e-10)
})

test_that("vehicle mass is non-increasing (K = 1, perfect sink)", {
  res <- run_minimal()
  diffs <- diff(res$mass$Vehicle)
  expect_true(all(diffs <= 1e-12))
})

test_that("sink mass is non-decreasing for a perfect sink", {
  res <- run_minimal()
  diffs <- diff(res$mass$Sink)
  expect_true(all(diffs >= -1e-12))
})

test_that("vehicle initial concentration matches the input c_init", {
  res <- run_minimal()
  # c_init = 1 mg/ml -> in default scaling 'mg', concentration[1] should be 1.
  expect_equal(res$concentration$Vehicle[1], 1.0, tolerance = 1e-12)
})

test_that("CDP profiles have the right shape and decay over depth", {
  res <- run_minimal()
  expect_named(res$cdp, c("Vehicle", "SC"))
  for (nm in names(res$cdp)) {
    s <- res$cdp[[nm]]
    expect_equal(length(s$depth_um), nrow(s$conc))
    expect_equal(length(s$time),    ncol(s$conc))
  }
  # At t = 0 the SC profile should be all zeros (c_init = 0 there).
  expect_true(all(res$cdp$SC$conc[, 1] == 0))
  # By the end, top of the SC has more mass than the bottom.
  end <- res$cdp$SC$conc[, ncol(res$cdp$SC$conc)]
  expect_gt(end[1], end[length(end)])
})

test_that("PK sink with finite Vd accumulates mass smoothly", {
  res <- skin_simulate(skin_params(
    vehicle = list(c_init = 1, height = 30L, D = 1, app_area = 1),
    layers  = data.frame(name = "SC", height = 20L, D = 1, K = 1, cross_section = 1),
    sink    = list(Vd = 10.0),
    pk      = list(enabled = TRUE, thalf = 1),  # 1-hour elimination
    sim_time = 60L
  ))
  expect_equal(res$status, "executed")
  expect_true(any(res$mass$Sink > 0))
})

test_that("infinite dose delivers more mass to the sink than finite dose", {
  finite <- skin_simulate(skin_params(
    vehicle = list(c_init = 1, height = 30L, D = 1, app_area = 1,
                   finite_dose = TRUE),
    layers  = data.frame(name = "SC", height = 20L, D = 1, K = 1, cross_section = 1),
    sink    = list(Vd = 1.0),
    sim_time = 60L
  ))
  infinite <- skin_simulate(skin_params(
    vehicle = list(c_init = 1, height = 30L, D = 1, app_area = 1,
                   finite_dose = FALSE),
    layers  = data.frame(name = "SC", height = 20L, D = 1, K = 1, cross_section = 1),
    sink    = list(Vd = 1.0),
    sim_time = 60L
  ))
  # Infinite dose holds the donor surface at c_init, so sink mass grows
  # at least as fast as in the finite-dose case.
  expect_gte(tail(infinite$mass$Sink, 1), tail(finite$mass$Sink, 1))
})

test_that("infinite-dose donor keeps the top donor cell at c_init", {
  res <- skin_simulate(skin_params(
    vehicle = list(c_init = 1, height = 30L, D = 1, app_area = 1,
                   finite_dose = FALSE, log_cdp = TRUE),
    layers  = data.frame(name = "SC", height = 20L, D = 1, K = 1,
                         cross_section = 1, log_cdp = FALSE),
    sink    = list(Vd = 1.0),
    sim_time = 30L
  ))
  top_cell_over_time <- res$cdp$Vehicle$conc[1, ]
  expect_true(all(abs(top_cell_over_time - 1.0) < 1e-9))
})

test_that("finer mesh converges to the coarser solution at the sink", {
  base <- run_minimal(resolution = 1L)
  fine <- run_minimal(resolution = 4L)
  rel <- abs(tail(fine$mass$Sink, 1) - tail(base$mass$Sink, 1)) /
         max(tail(fine$mass$Sink, 1), 1e-30)
  expect_lt(rel, 0.05)
})

test_that("print and summary methods work on a real result", {
  res <- run_minimal()
  expect_output(print(res), "skin_result")
  expect_output(summary(res), "summary")
})
