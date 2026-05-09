# Build a fresh, valid set of parameters using the new constructor API.
# Every unit-bearing argument carries an explicit units object.

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
                   D = um2_per_min(1.0))
  args <- modifyList(defaults, list(...))
  do.call(vehicle, args)
}

layer_default <- function(...) {
  defaults <- list(name = "SC", height = um(20L), D = um2_per_min(1.0),
                   K = 1.0, cross_section = 1.0)
  args <- modifyList(defaults, list(...))
  do.call(layer, args)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------- constructor smoke tests ----------

test_that("vehicle() returns a skin_vehicle object", {
  v <- vehicle(c_init = mg_per_ml(1), height = um(30), D = um2_per_min(1))
  expect_s3_class(v, "skin_vehicle")
  expect_equal(v$c_init_mg_per_ml, 1)
  expect_equal(v$height_um, 30L)
  expect_true(v$finite_dose)
})

test_that("layer() returns a skin_layer object", {
  l <- layer("SC", height = um(20), D = um2_per_min(1), K = 1, cross_section = 1)
  expect_s3_class(l, "skin_layer")
  expect_equal(l$name, "SC")
  expect_equal(l$K, 1)
})

test_that("perfect_sink() returns a skin_sink with type 'perfect'", {
  s <- perfect_sink("Receptor")
  expect_s3_class(s, "skin_sink")
  expect_equal(s$type, "perfect")
  expect_equal(s$name, "Receptor")
})

test_that("finite_sink() returns a skin_sink with type 'finite'", {
  s <- finite_sink("Blood", Vd = ml(100))
  expect_s3_class(s, "skin_sink")
  expect_equal(s$type, "finite")
  expect_equal(s$Vd_ml, 100)
})

# ---------- unit conversion through the constructors ----------

test_that("vehicle accepts compatible units in non-canonical scale", {
  # mm, um -> equivalent for height
  v1 <- vehicle(c_init = mg_per_ml(1), height = um(110), D = um2_per_min(1))
  v2 <- vehicle(c_init = mg_per_ml(1), height = mm(0.110), D = um2_per_min(1))
  expect_equal(v1$height_um, 110L)
  expect_equal(v2$height_um, 110L)
})

test_that("layer height accepts cm and converts to integer micrometres", {
  # mm(0.2055) = 205.5 um -> rounds with an informational message
  expect_message(
    l <- layer("L", height = mm(0.2055), D = um2_per_min(1),
               K = 1, cross_section = 1),
    "rounding"
  )
  expect_equal(l$height_um, 206L)

  # cm(0.0205) = 205 um exactly -> no rounding, no message
  expect_silent(
    l2 <- layer("L", height = cm(0.0205), D = um2_per_min(1),
                K = 1, cross_section = 1)
  )
  expect_equal(l2$height_um, 205L)
})

test_that("finite_sink Vd accepts millilitre alternative units", {
  s1 <- finite_sink("S", Vd = ml(100))
  s2 <- finite_sink("S", Vd = units::set_units(0.1, "L"))
  expect_equal(s1$Vd_ml, 100)
  expect_equal(s2$Vd_ml, 100)
})

test_that("duration accepts hours and days", {
  p1 <- make_minimal(duration = minutes(60L))
  p2 <- make_minimal(duration = hours(1L))
  expect_equal(p1$sys$simulation_time, 60L)
  expect_equal(p2$sys$simulation_time, 60L)
})

# ---------- strict-units rejection ----------

test_that("vehicle rejects bare numeric c_init", {
  expect_error(
    vehicle(c_init = 1, height = um(30), D = um2_per_min(1)),
    "units-aware"
  )
})

test_that("vehicle rejects incompatible unit on height", {
  expect_error(
    vehicle(c_init = mg_per_ml(1), height = mg_per_ml(30), D = um2_per_min(1)),
    "incompatible unit"
  )
})

test_that("vehicle rejects bare numeric D", {
  expect_error(
    vehicle(c_init = mg_per_ml(1), height = um(30), D = 1),
    "units-aware"
  )
})

test_that("layer rejects units value for K (dimensionless)", {
  expect_error(
    layer("L", height = um(20), D = um2_per_min(1),
          K = um(1), cross_section = 1),
    "dimensionless"
  )
})

# ---------- skin_params composition ----------

test_that("skin_params builds a list with class skin_params", {
  p <- make_minimal()
  expect_s3_class(p, "skin_params")
  expect_true(is.list(p))
  expect_named(p, c("sys", "log", "sink", "vehicle", "layers", ".meta"))
})

test_that("skin_params accepts multiple layers", {
  p <- make_minimal(layers = list(
    layer("SC",  height = um(20L), D = um2_per_min(1), K = 0.5, cross_section = 1.0),
    layer("DSL", height = um(30L), D = um2_per_min(2), K = 0.1, cross_section = 0.5)
  ))
  expect_length(p$layers, 2L)
  expect_equal(p$layers[[1]]$name, "SC")
  expect_equal(p$layers[[2]]$cross_section, 0.5)
})

test_that("skin_params requires a vehicle() object", {
  expect_error(
    skin_params(area = cm2(1), vehicle = list(),
                layers = list(layer_default()), sink = perfect_sink(),
                duration = minutes(30L)),
    "skin_vehicle"
  )
})

test_that("skin_params requires a sink object", {
  expect_error(
    skin_params(area = cm2(1), vehicle = vehicle_default(),
                layers = list(layer_default()), sink = list(),
                duration = minutes(30L)),
    "skin_sink"
  )
})

test_that("skin_params requires layers to be non-empty", {
  expect_error(
    skin_params(area = cm2(1), vehicle = vehicle_default(),
                layers = list(), sink = perfect_sink(),
                duration = minutes(30L)),
    "non-empty"
  )
})

test_that("skin_params requires every layers element to be a layer()", {
  expect_error(
    skin_params(area = cm2(1), vehicle = vehicle_default(),
                layers = list(list(name = "SC")), sink = perfect_sink(),
                duration = minutes(30L)),
    "skin_layer"
  )
})

# ---------- range checks ----------

test_that("layer cross_section out of range names the constraint", {
  expect_error(
    layer("x", height = um(10), D = um2_per_min(1), K = 1, cross_section = 1.5),
    "cross_section"
  )
})

test_that("vehicle$D negative is rejected", {
  expect_error(
    vehicle(c_init = mg_per_ml(1), height = um(30L),
            D = um2_per_min(-1)),
    "out of range"
  )
})

test_that("invalid duration is rejected", {
  expect_error(make_minimal(duration = minutes(0L)), "duration")
})

test_that("invalid area is rejected", {
  expect_error(make_minimal(area = cm2(-1)), "area")
})

test_that("vehicle$height below 3 um is rejected", {
  expect_error(
    vehicle(c_init = mg_per_ml(1), height = um(2L), D = um2_per_min(1)),
    "out of range"
  )
})

# ---------- scaling normalisation ----------

test_that("scaling values are normalised", {
  for (s in c("mg", "ug", "ng")) {
    p <- make_minimal(scaling = s)
    expect_equal(p$log$scaling, s)
  }
})

# ---------- print methods ----------

test_that("print methods run without error and show units", {
  expect_output(print(vehicle_default()), "skindiff vehicle")
  expect_output(print(vehicle_default()), "\\[mg/ml\\]")
  expect_output(print(vehicle_default()), "\\[um\\]")
  expect_output(print(layer_default()), "skindiff layer")
  expect_output(print(perfect_sink()), "perfect")
  expect_output(print(finite_sink("S", Vd = ml(1))), "finite")
  expect_output(print(make_minimal()), "skin_params")
})
