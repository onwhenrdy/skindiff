# Coverage for the ggplot2-based autoplot methods. ggplot2 is in Suggests,
# so each test guards on availability.

make_plot_run <- function(log_cdp_layer = TRUE, sink_finite = FALSE) {
  snk <- if (sink_finite) finite_sink("Sink", Vd = ml(1.0))
         else             perfect_sink("Sink")
  skin_simulate(skin_params(
    area = cm2(1.0),
    vehicle = vehicle(
      c_init = mg_per_ml(1.0), height = um(30L),
      D = um2_per_min(1.0),
      log_mass = TRUE, log_cdp = TRUE
    ),
    layers = list(
      layer("SC", height = um(20L), D = um2_per_min(1.0),
            K = 1.0, cross_section = 1.0,
            log_mass = TRUE, log_cdp = log_cdp_layer)
    ),
    sink = snk,
    duration = minutes(60L),
    resolution = 4L,
    scaling = "ng"
  ))
}

# ---------- dispatch and basic structure ------------------------------------

test_that("autoplot returns a ggplot for each `what` value", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  for (w in c("mass", "concentration", "permeated", "flux", "profile")) {
    p <- ggplot2::autoplot(res, what = w)
    expect_s3_class(p, "ggplot")
  }
})

test_that("default `what` is mass", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p_default <- ggplot2::autoplot(res)
  p_mass    <- ggplot2::autoplot(res, what = "mass")
  expect_equal(class(p_default), class(p_mass))
})

# ---------- per-plot structure ----------------------------------------------

test_that("mass plot has one line geom and labels in scaling units", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "mass")
  expect_length(p$layers, 1L)
  expect_s3_class(p$layers[[1]]$geom, "GeomLine")
  expect_match(p$labels$y, "ng")           # scaling = "ng" set above
  expect_match(p$labels$x, "min")          # canonical time unit
})

test_that("concentration plot drops perfect-sink (all-NA) column", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run(sink_finite = FALSE)
  p <- ggplot2::autoplot(res, what = "concentration")
  comps <- levels(p$data$compartment)
  # Sink column is all-NA for perfect_sink; should be dropped.
  expect_false("Sink" %in% comps)
  expect_true(all(c("Vehicle", "SC") %in% comps))
})

test_that("concentration plot keeps finite sink", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run(sink_finite = TRUE)
  p <- ggplot2::autoplot(res, what = "concentration")
  comps <- levels(p$data$compartment)
  expect_true("Sink" %in% comps)
})

test_that("permeated plot has Q-vs-time labels", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "permeated")
  expect_match(p$labels$y, "Q")
  expect_match(p$labels$y, "ng")
  expect_match(p$labels$y, "cm-2")         # ng/cm^2 in canonical udunits
})

test_that("flux plot has flux-units label", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "flux")
  expect_match(p$labels$y, "flux")
  expect_match(p$labels$y, "ng cm-2 min-1")
})

test_that("profile errors when no cdp is logged", {
  skip_if_not_installed("ggplot2")
  res <- skin_simulate(skin_params(
    area = cm2(1.0),
    vehicle = vehicle(c_init = mg_per_ml(1.0), height = um(30L),
                      D = um2_per_min(1.0),
                      log_mass = TRUE, log_cdp = FALSE),
    layers = list(layer("SC", height = um(20L), D = um2_per_min(1.0),
                        K = 1.0, cross_section = 1.0,
                        log_mass = TRUE, log_cdp = FALSE)),
    sink = perfect_sink("Sink"),
    duration = minutes(30L)
  ))
  expect_error(ggplot2::autoplot(res, what = "profile"), "log_cdp")
})

# ---------- profile (line plot) ---------------------------------------------

test_that("profile plot is a single-axis line plot with one geom_line", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "profile")
  # Geoms: dashed compartment-boundary vlines + one line layer
  geom_classes <- vapply(p$layers, function(l) class(p$layers[[1]]$geom)[1L],
                         character(1L))
  expect_true(any(vapply(p$layers,
                         function(l) inherits(l$geom, "GeomLine"), logical(1L))))
  expect_match(p$labels$x, "depth")
  expect_match(p$labels$x, "skin surface = 0")
})

test_that("profile default uses 6 equidistant times", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "profile")
  # 6 distinct time levels in the underlying data
  expect_equal(length(levels(p$data$time_label)), 6L)
})

test_that("profile respects explicit `times`", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "profile",
                        times = minutes(c(10, 30, 50)))
  expect_equal(length(levels(p$data$time_label)), 3L)
})

test_that("profile respects `n_times`", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "profile", n_times = 4)
  expect_equal(length(levels(p$data$time_label)), 4L)
})

test_that("profile excludes vehicle when include_vehicle = FALSE", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  vehicle_name <- res$params$vehicle$name
  p_with    <- ggplot2::autoplot(res, what = "profile", include_vehicle = TRUE)
  p_without <- ggplot2::autoplot(res, what = "profile", include_vehicle = FALSE)
  expect_true (vehicle_name %in% levels(p_with$data$compartment))
  expect_false(vehicle_name %in% levels(p_without$data$compartment))
})

test_that("profile uses skin-surface = 0 convention (vehicle at negative depth)", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "profile", include_vehicle = TRUE)
  vehicle_data <- p$data[p$data$compartment == res$params$vehicle$name, ]
  expect_true(all(vehicle_data$depth_global < 0))

  sc_data <- p$data[p$data$compartment == "SC", ]
  expect_true(all(sc_data$depth_global >= 0))
})

test_that("profile accepts an explicit subset of compartments", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "profile", compartments = "SC")
  expect_equal(levels(p$data$compartment), "SC")
})

test_that("profile errors on unknown compartment name", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  expect_error(
    ggplot2::autoplot(res, what = "profile", compartments = "garbage"),
    "Not logged or unknown"
  )
})

test_that("profile rejects bare-numeric `times`", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  expect_error(
    ggplot2::autoplot(res, what = "profile", times = c(10, 30)),
    "units-of-time"
  )
})

test_that("profile rejects `times` outside the simulation range", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()           # 60-min sim
  expect_error(
    ggplot2::autoplot(res, what = "profile", times = hours(2)),
    "outside the simulation"
  )
})

# ---------- composability ---------------------------------------------------

test_that("returned ggplot can be extended with + layers", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  p <- ggplot2::autoplot(res, what = "permeated") +
       ggplot2::ggtitle("custom title")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "custom title")
})

# ---------- bad input handling ----------------------------------------------

test_that("autoplot rejects unknown `what`", {
  skip_if_not_installed("ggplot2")
  res <- make_plot_run()
  expect_error(ggplot2::autoplot(res, what = "garbage"), "should be one of")
})
