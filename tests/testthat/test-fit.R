# Coverage for skin_fit(): synthetic round-trip recovery (simulate with
# known D/K, fit, check we recover them), multi-subject, joint
# permeation+penetration, partial fitting (D fixed, K fitted), and the
# main sanity-error paths.

# Build a single-layer template with known D and K, run it, sample a few
# permeation points, and use those to fit D and K back. The simulator
# is deterministic so the round-trip should converge exactly to the
# starting values (within optimizer tolerance).

make_one_layer_template <- function(D = 100, K = 1, h = 100,
                                    c_init = 1, h_v = 50,
                                    duration_min = 480L,
                                    log_cdp = FALSE,
                                    finite_dose = TRUE) {
  skin_params(
    area = cm2(1.0),
    vehicle = vehicle(
      c_init = mg_per_ml(c_init), height = um(h_v),
      D = um2_per_min(1000),
      finite_dose = finite_dose,
      log_mass = TRUE, log_cdp = FALSE
    ),
    layers = list(
      layer("Skin", height = um(h), D = um2_per_min(D),
            K = K, cross_section = 1.0,
            log_mass = TRUE, log_cdp = log_cdp)
    ),
    sink = perfect_sink("Receptor"),
    duration = minutes(duration_min),
    resolution = 4L,
    scaling = "ng"
  )
}

make_two_layer_template <- function(D_sc = 1, K_sc = 50, h_sc = 20,
                                    D_dsl = 200, K_dsl = 1, h_dsl = 100,
                                    c_init = 1, h_v = 50,
                                    duration_min = 480L,
                                    log_cdp = FALSE) {
  skin_params(
    area = cm2(1.0),
    vehicle = vehicle(
      c_init = mg_per_ml(c_init), height = um(h_v),
      D = um2_per_min(1000),
      finite_dose = TRUE, log_mass = TRUE
    ),
    layers = list(
      layer("Stratum corneum", height = um(h_sc),
            D = um2_per_min(D_sc), K = K_sc, cross_section = 1.0,
            log_mass = TRUE, log_cdp = log_cdp),
      layer("Dermis", height = um(h_dsl),
            D = um2_per_min(D_dsl), K = K_dsl, cross_section = 1.0,
            log_mass = TRUE, log_cdp = log_cdp)
    ),
    sink = perfect_sink("Receptor"),
    duration = minutes(duration_min),
    resolution = 4L,
    scaling = "ng"
  )
}

# Simulate the truth, sample times, build a permeation_obs.
sample_permeation <- function(truth_template,
                              sample_times_min = c(30, 60, 120, 180, 240, 360, 480)) {
  res <- skin_simulate(truth_template)
  q_at <- permeated_at(res, minutes(sample_times_min))
  permeation_obs(data.frame(
    time       = minutes(sample_times_min),
    q_per_area = units::set_units(as.numeric(q_at), "ng/cm^2", mode = "standard")
  ))
}

sample_penetration <- function(truth_template,
                               t_min = 240, depths = c(0, 5, 10, 20, 40, 80)) {
  res <- skin_simulate(truth_template)
  prof <- profile_at(res, minutes(t_min))
  rows <- list()
  for (nm in names(prof)) {
    df <- prof[[nm]]
    for (i in seq_along(depths[-length(depths)])) {
      d_top <- depths[i]
      d_bot <- depths[i + 1L]
      # Need to know the layer offset to keep within-layer
      # depths in skin frame. But profile_at returns layer-local depth.
      # Easier: use the cdp data directly and pick midpoints in skin frame.
    }
  }
  # Simpler: build penetration data from the cdp directly in skin-frame.
  cdp_t_idx <- which.min(abs(as.numeric(res$cdp[[1L]]$time) - t_min))
  band_rows <- list()
  cum_top <- 0
  for (l in res$params$layers) {
    nm <- l$name
    s <- res$cdp[[nm]]
    ds <- as.numeric(s$depth)
    cs <- as.numeric(s$conc[, cdp_t_idx])
    # Build a strip per cell in skin frame
    if (length(ds) > 1L) {
      dx <- ds[2L] - ds[1L]
    } else {
      dx <- l$height
    }
    skin_mid <- cum_top + ds
    band_rows[[length(band_rows) + 1L]] <- data.frame(
      time = rep(minutes(t_min), length(skin_mid)),
      depth_top = um(skin_mid - dx / 2),
      depth_bottom = um(skin_mid + dx / 2),
      concentration = units::set_units(cs, "ng/ml", mode = "standard")
    )
    cum_top <- cum_top + l$height
  }
  d <- do.call(rbind, band_rows)
  penetration_obs(d)
}


# ---------- single-layer permeation round-trip ------------------------------

test_that("permeation-only fit recovers known D and K (single layer)", {
  truth <- make_one_layer_template(D = 100, K = 1)
  obs   <- sample_permeation(truth)

  # Start from a deliberately wrong guess, ensure we converge close to truth.
  template <- make_one_layer_template(D = 10, K = 5)
  fit <- skin_fit(
    template = template,
    observations = list(permeation = obs),
    fit_pars = list("Skin" = c("D", "K"))
  )
  expect_s3_class(fit, "skin_fit")
  expect_equal(fit$convergence, 0L)

  est <- coef(fit)
  expect_equal(est[["D[Skin]"]], 100, tolerance = 0.05)   # 5% relative
  expect_equal(est[["K[Skin]"]], 1,   tolerance = 0.05)
})


# ---------- partial fitting (only D) ----------------------------------------

test_that("permeation fit with only D (K fixed) recovers D", {
  truth <- make_one_layer_template(D = 50, K = 2)

  # Build obs from truth.
  obs <- sample_permeation(truth)

  # Template has correct K but wrong D.
  template <- make_one_layer_template(D = 5, K = 2)
  fit <- skin_fit(
    template = template,
    observations = list(permeation = obs),
    fit_pars = list("Skin" = "D")
  )
  est <- coef(fit)
  expect_named(est, "D[Skin]")
  expect_equal(est[["D[Skin]"]], 50, tolerance = 0.05)
})


# ---------- two-layer fit, fix DSL via QSAR --------------------------------

test_that("fit only SC's D/K with DSL fixed at QSAR values", {
  truth <- make_two_layer_template(D_sc = 1, K_sc = 50,
                                    D_dsl = 200, K_dsl = 1)
  obs <- sample_permeation(truth, sample_times_min = c(30, 60, 120, 240, 360, 480))

  template <- make_two_layer_template(D_sc = 0.5, K_sc = 20,
                                       D_dsl = 200, K_dsl = 1)
  fit <- skin_fit(
    template = template,
    observations = list(permeation = obs),
    fit_pars = list("Stratum corneum" = c("D", "K"))
  )
  est <- coef(fit)
  expect_equal(est[["D[Stratum corneum]"]], 1,  tolerance = 0.10)
  expect_equal(est[["K[Stratum corneum]"]], 50, tolerance = 0.10)
})


# ---------- penetration-only fit -------------------------------------------

test_that("penetration-only fit recovers D for a single layer", {
  truth <- make_one_layer_template(D = 50, K = 1, h = 100,
                                   log_cdp = TRUE)
  obs <- sample_penetration(truth, t_min = 60)

  template <- make_one_layer_template(D = 5, K = 1, h = 100, log_cdp = TRUE)
  fit <- skin_fit(
    template = template,
    observations = list(penetration = obs),
    fit_pars = list("Skin" = "D")
  )
  est <- coef(fit)
  expect_equal(est[["D[Skin]"]], 50, tolerance = 0.10)
})


# ---------- multi-subject fit ----------------------------------------------

test_that("multi-subject fit shares D/K across subjects with different geometry", {
  truth1 <- make_one_layer_template(D = 100, K = 1, h = 100)
  truth2 <- make_one_layer_template(D = 100, K = 1, h = 150)   # thicker SC

  obs1 <- sample_permeation(truth1)
  obs2 <- sample_permeation(truth2)

  perm <- permeation_obs(data.frame(
    subject    = rep(c("s1", "s2"), each = obs1$n),
    time       = minutes(c(obs1$time_min, obs2$time_min)),
    q_per_area = units::set_units(c(obs1$q_per_area_ng_cm2,
                                     obs2$q_per_area_ng_cm2),
                                   "ng/cm^2", mode = "standard")
  ))

  template <- list(
    s1 = make_one_layer_template(D = 10, K = 5, h = 100),
    s2 = make_one_layer_template(D = 10, K = 5, h = 150)
  )
  fit <- skin_fit(
    template = template,
    observations = list(permeation = perm),
    fit_pars = list("Skin" = c("D", "K"))
  )
  est <- coef(fit)
  expect_equal(est[["D[Skin]"]], 100, tolerance = 0.10)
  expect_equal(est[["K[Skin]"]], 1,   tolerance = 0.10)
})


# ---------- input validation ----------------------------------------------

test_that("skin_fit errors on unknown layer name in fit_pars", {
  template <- make_one_layer_template()
  obs <- sample_permeation(template)
  expect_error(
    skin_fit(template, list(permeation = obs),
             fit_pars = list("Unknown" = "D")),
    "not present in the template"
  )
})

test_that("skin_fit errors on non-D/K parameter name", {
  template <- make_one_layer_template()
  obs <- sample_permeation(template)
  expect_error(
    skin_fit(template, list(permeation = obs),
             fit_pars = list("Skin" = "h")),
    "Only.*D.*K.*fittable"
  )
})

test_that("skin_fit errors when sink mass is not logged but permeation given", {
  tpl <- skin_params(
    area = cm2(1.0),
    vehicle = vehicle(c_init = mg_per_ml(1), height = um(50L),
                      D = um2_per_min(1000),
                      log_mass = TRUE),
    layers = list(layer("Skin", height = um(100L), D = um2_per_min(100),
                        K = 1, cross_section = 1, log_mass = TRUE)),
    sink = perfect_sink("Sink", log_mass = FALSE),    # not logged
    duration = minutes(60L), scaling = "ng"
  )
  obs <- permeation_obs(data.frame(
    time = minutes(c(30, 60)),
    q_per_area = ng_per_cm2(c(1, 2))
  ))
  expect_error(
    skin_fit(tpl, list(permeation = obs),
             fit_pars = list("Skin" = "D")),
    "sink mass is not logged"
  )
})

test_that("skin_fit errors when penetration overlaps non-cdp layer", {
  tpl <- make_one_layer_template(log_cdp = FALSE)
  pen <- penetration_obs(data.frame(
    time = minutes(60), depth_top = um(0), depth_bottom = um(50),
    concentration = ng_per_ml(1000)
  ))
  expect_error(
    suppressWarnings(skin_fit(tpl, list(penetration = pen),
                               fit_pars = list("Skin" = "D"))),
    "not CDP-logged"
  )
})

test_that("skin_fit warns when bounds exclude the template value", {
  template <- make_one_layer_template(D = 100, K = 1)
  obs <- sample_permeation(template)
  expect_warning(
    fit <- skin_fit(
      template, list(permeation = obs),
      fit_pars = list("Skin" = "D"),
      bounds   = list("Skin" = list(D = c(um2_per_min(0.001),
                                          um2_per_min(0.01))))
    ),
    "exclude the template"
  )
})

test_that("skin_fit warns when bounds are supplied for a non-fitted parameter", {
  template <- make_one_layer_template()
  obs <- sample_permeation(template)
  expect_warning(
    skin_fit(
      template, list(permeation = obs),
      fit_pars = list("Skin" = "D"),
      bounds   = list("Skin" = list(K = c(1, 100)))   # K not in fit_pars
    ),
    "ignored"
  )
})

test_that("skin_fit errors when subject in obs is missing from template", {
  template <- list(s1 = make_one_layer_template())
  obs <- permeation_obs(data.frame(
    subject = c("s1", "s99"),
    time    = minutes(c(60, 60)),
    q_per_area = ng_per_cm2(c(1, 2))
  ))
  expect_error(
    skin_fit(template, list(permeation = obs),
             fit_pars = list("Skin" = "D")),
    "no template"
  )
})

test_that("skin_fit S3 methods produce valid output", {
  truth <- make_one_layer_template(D = 100, K = 1)
  obs <- sample_permeation(truth)
  template <- make_one_layer_template(D = 50, K = 2)
  fit <- skin_fit(template, list(permeation = obs),
                  fit_pars = list("Skin" = c("D", "K")))

  expect_output(print(fit), "skin_fit")
  expect_output(summary(fit), "RMSE")

  cf <- coef(fit)
  expect_named(cf, c("D[Skin]", "K[Skin]"))

  res <- residuals(fit)
  expect_true(is.list(res))
  expect_true("permeation" %in% names(res))

  fp <- fitted(fit)
  expect_true(is.list(fp))

  best_params <- skin_params_from_fit(fit)
  expect_s3_class(best_params, "skin_params")
})
