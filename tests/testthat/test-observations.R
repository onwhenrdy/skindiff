# Coverage for the observation constructors permeation_obs() and
# penetration_obs(). Strict-units validation, column requirements,
# subject handling.

# ---------- permeation_obs ----------------------------------------------------

test_that("permeation_obs accepts a tidy data.frame with required columns", {
  d <- data.frame(
    time       = hours(c(1, 2, 4)),
    q_per_area = ng_per_cm2(c(10, 30, 80))
  )
  o <- permeation_obs(d)
  expect_s3_class(o, "permeation_obs")
  expect_equal(o$n, 3L)
  expect_equal(o$time_min, c(60, 120, 240))           # canonical: min
  expect_equal(o$q_per_area_ng_cm2, c(10, 30, 80))    # canonical: ng/cm^2
  expect_true(all(is.na(o$sd_ng_cm2)))
  expect_equal(unique(o$subject), "default")
})

test_that("permeation_obs converts compatible units", {
  d <- data.frame(
    time       = minutes(c(60, 120)),
    q_per_area = ug_per_cm2(c(0.1, 0.2))
  )
  o <- permeation_obs(d)
  expect_equal(o$time_min, c(60, 120))
  expect_equal(o$q_per_area_ng_cm2, c(100, 200))      # ug -> ng = *1000
})

test_that("permeation_obs preserves sd and subject when supplied", {
  d <- data.frame(
    subject    = rep(c("s1", "s2"), each = 2),
    time       = rep(hours(c(1, 4)), times = 2),
    q_per_area = ng_per_cm2(c(10, 80, 12, 90)),
    sd         = ng_per_cm2(c(1, 5, 1.2, 6))
  )
  o <- permeation_obs(d)
  expect_equal(o$n, 4L)
  expect_equal(o$subject, c("s1", "s1", "s2", "s2"))
  expect_equal(o$sd_ng_cm2, c(1, 5, 1.2, 6))
})

test_that("permeation_obs rejects missing required columns", {
  expect_error(
    permeation_obs(data.frame(time = hours(1))),
    "missing required column"
  )
})

test_that("permeation_obs rejects bare-numeric time", {
  d <- data.frame(time = c(1, 2), q_per_area = ng_per_cm2(c(1, 2)))
  expect_error(permeation_obs(d), "units-aware")
})

test_that("permeation_obs rejects incompatible units on q_per_area", {
  d <- data.frame(time = hours(c(1, 2)),
                  q_per_area = ng_per_ml(c(1, 2)))   # wrong: should be ng/cm^2
  expect_error(permeation_obs(d), "incompatible unit")
})

test_that("permeation_obs rejects negative q_per_area", {
  d <- data.frame(time = hours(c(1, 2)),
                  q_per_area = ng_per_cm2(c(-1, 2)))
  expect_error(permeation_obs(d), "non-negative")
})

test_that("permeation_obs rejects empty data.frame", {
  expect_error(
    permeation_obs(data.frame(time = hours(numeric(0)),
                              q_per_area = ng_per_cm2(numeric(0)))),
    "at least one row"
  )
})


# ---------- penetration_obs --------------------------------------------------

test_that("penetration_obs accepts a tidy data.frame with required columns", {
  d <- data.frame(
    time          = rep(hours(24), 3),
    depth_top     = um(c(0,  2,  4)),
    depth_bottom  = um(c(2,  4,  6)),
    concentration = ug_per_ml(c(2.0, 1.0, 0.5))
  )
  o <- penetration_obs(d)
  expect_s3_class(o, "penetration_obs")
  expect_equal(o$n, 3L)
  expect_equal(o$time_min, rep(24 * 60, 3))
  expect_equal(o$depth_top_um, c(0, 2, 4))
  expect_equal(o$depth_bottom_um, c(2, 4, 6))
  expect_equal(o$conc_ng_ml, c(2000, 1000, 500))      # ug -> ng = *1000
})

test_that("penetration_obs accepts point measurements (top == bottom)", {
  d <- data.frame(
    time          = hours(c(24, 24)),
    depth_top     = um(c(5, 10)),
    depth_bottom  = um(c(5, 10)),
    concentration = ng_per_ml(c(800, 200))
  )
  o <- penetration_obs(d)
  expect_equal(o$n, 2L)
  expect_equal(o$depth_top_um, o$depth_bottom_um)
})

test_that("penetration_obs rejects depth_top > depth_bottom", {
  d <- data.frame(
    time          = hours(c(24, 24)),
    depth_top     = um(c(0, 4)),
    depth_bottom  = um(c(2, 2)),                      # second row inverted
    concentration = ng_per_ml(c(1000, 500))
  )
  expect_error(penetration_obs(d), "depth_top.*<=.*depth_bottom")
})

test_that("penetration_obs rejects negative depths", {
  d <- data.frame(
    time          = hours(24),
    depth_top     = um(-2),
    depth_bottom  = um(0),
    concentration = ng_per_ml(1000)
  )
  expect_error(penetration_obs(d), "skin surface = 0")
})

test_that("penetration_obs preserves sd and subject", {
  d <- data.frame(
    subject       = c("s1", "s1", "s2"),
    time          = hours(rep(24, 3)),
    depth_top     = um(c(0, 2, 0)),
    depth_bottom  = um(c(2, 4, 2)),
    concentration = ng_per_ml(c(2000, 1000, 1500)),
    sd            = ng_per_ml(c(200, 100, 200))
  )
  o <- penetration_obs(d)
  expect_equal(o$subject, c("s1", "s1", "s2"))
  expect_equal(o$sd_ng_ml, c(200, 100, 200))
})

test_that("penetration_obs rejects bare-numeric depth", {
  d <- data.frame(
    time          = rep(hours(24), 2),
    depth_top     = c(0, 2),
    depth_bottom  = c(2, 4),
    concentration = ng_per_ml(c(1, 1))
  )
  expect_error(penetration_obs(d), "units-aware")
})


# ---------- print methods ---------------------------------------------------

test_that("print methods run without error", {
  o1 <- permeation_obs(data.frame(time = hours(1), q_per_area = ng_per_cm2(1)))
  o2 <- penetration_obs(data.frame(
    time = hours(24), depth_top = um(0), depth_bottom = um(2),
    concentration = ng_per_ml(1000)
  ))
  expect_output(print(o1), "permeation_obs")
  expect_output(print(o2), "penetration_obs")
})
