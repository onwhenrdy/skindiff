make_minimal <- function(vehicle = NULL, layers = NULL, sink = NULL,
                         duration = NULL, area = NULL, ...) {
  v <- vehicle  %||% vehicle_default()
  l <- layers   %||% list(layer_default())
  s <- sink     %||% perfect_sink()
  d <- duration %||% minutes(30L)
  a <- area     %||% cm2(1.0)
  skin_params(area = a, vehicle = v, layers = l, sink = s,
              duration = d, ...)
}

vehicle_default <- function(...) {
  defaults <- list(c_init = mg_per_ml(1.0), height = um(30L),
                   D = um2_per_min(1.0), log_mass = TRUE, log_cdp = TRUE)
  args <- modifyList(defaults, list(...))
  do.call(vehicle, args)
}

layer_default <- function(...) {
  defaults <- list(name = "SC", height = um(20L), D = um2_per_min(1.0),
                   K = 1.0, cross_section = 1.0,
                   log_mass = TRUE, log_cdp = TRUE)
  args <- modifyList(defaults, list(...))
  do.call(layer, args)
}

run_minimal <- function(...) skin_simulate(make_minimal(...))

`%||%` <- function(a, b) if (is.null(a)) b else a

test_that("skin_simulate returns the documented structure", {
  res <- run_minimal()
  expect_s3_class(res, "skin_result")
  expect_equal(res$status, "executed")
  expect_named(res, c("status", "scaling", "mass", "concentration",
                      "cdp", "geometry", "params", "runtime"))
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

test_that("perfect-sink concentration column is NA", {
  res <- run_minimal()
  expect_true(all(is.na(res$concentration$Sink)))
})

test_that("finite-sink concentration column is finite and non-decreasing", {
  res <- skin_simulate(make_minimal(sink = finite_sink("Sink", Vd = ml(1.0))))
  expect_true(all(is.finite(res$concentration$Sink)))
  expect_true(all(as.numeric(diff(res$concentration$Sink)) >= -1e-12))
})

test_that("mass is conserved for a perfect-sink run (relative err < 1e-10)", {
  res <- run_minimal()
  # rowSums on units-bearing columns is unsafe; sum directly instead.
  totals <- res$mass$Vehicle + res$mass$SC + res$mass$Sink
  rel <- as.numeric(max(abs(totals - totals[1])) / totals[1])
  expect_lt(rel, 1e-10)
})

test_that("vehicle mass is non-increasing (K = 1, perfect sink)", {
  res <- run_minimal()
  diffs <- as.numeric(diff(res$mass$Vehicle))
  expect_true(all(diffs <= 1e-12))
})

test_that("sink mass is non-decreasing for a perfect sink", {
  res <- run_minimal()
  diffs <- as.numeric(diff(res$mass$Sink))
  expect_true(all(diffs >= -1e-12))
})

test_that("vehicle initial concentration matches the input c_init", {
  res <- run_minimal()
  # c_init = mg_per_ml(1) -> in default scaling 'mg', concentration[1] should be 1.
  expect_equal(as.numeric(res$concentration$Vehicle[1]), 1.0,
               tolerance = 1e-12)
})

test_that("CDP profiles have the right shape and decay over depth", {
  res <- run_minimal()
  expect_named(res$cdp, c("Vehicle", "SC"))
  for (nm in names(res$cdp)) {
    s <- res$cdp[[nm]]
    expect_equal(length(s$depth), nrow(s$conc))
    expect_equal(length(s$time),  ncol(s$conc))
  }
  # At t = 0 the SC profile should be all zeros (c_init = 0 there).
  expect_true(all(as.numeric(res$cdp$SC$conc[, 1]) == 0))
  # By the end, top of the SC has more mass than the bottom.
  end <- as.numeric(res$cdp$SC$conc[, ncol(res$cdp$SC$conc)])
  expect_gt(end[1], end[length(end)])
})

test_that("infinite dose delivers more mass to the sink than finite dose", {
  finite <- skin_simulate(make_minimal(
    duration = minutes(60L),
    vehicle = vehicle_default(finite_dose = TRUE)
  ))
  infinite <- skin_simulate(make_minimal(
    duration = minutes(60L),
    vehicle = vehicle_default(finite_dose = FALSE)
  ))
  # Infinite dose holds the donor surface at c_init, so sink mass grows
  # at least as fast as in the finite-dose case.
  expect_gte(tail(infinite$mass$Sink, 1), tail(finite$mass$Sink, 1))
})

test_that("infinite-dose donor keeps the top donor cell at c_init", {
  res <- run_minimal(
    vehicle = vehicle_default(finite_dose = FALSE, log_cdp = TRUE),
    layers  = list(layer_default(log_cdp = FALSE))
  )
  top_cell_over_time <- as.numeric(res$cdp$Vehicle$conc[1, ])
  expect_true(all(abs(top_cell_over_time - 1.0) < 1e-9))
})

test_that("finer mesh converges to the coarser solution at the sink", {
  base <- run_minimal(resolution = 1L)
  fine <- run_minimal(resolution = 4L)
  fine_final <- as.numeric(tail(fine$mass$Sink, 1))
  base_final <- as.numeric(tail(base$mass$Sink, 1))
  rel <- abs(fine_final - base_final) / max(fine_final, 1e-30)
  expect_lt(rel, 0.05)
})

test_that("print and summary methods work on a real result", {
  res <- run_minimal()
  expect_output(print(res), "skin_result")
  expect_output(summary(res), "summary")
})
