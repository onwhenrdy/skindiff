# Build a fresh, valid set of parameter args, override individual fields by
# replacement (no merging of data.frames).
make_minimal <- function(vehicle = NULL, layers = NULL, sink = NULL,
                         pk = NULL, sim_time = NULL, ...) {
  v <- vehicle %||% list(c_init = 1.0, height = 30L, D = 1.0, app_area = 1.0)
  l <- if (is.null(layers))
         data.frame(name = "SC", height = 20L, D = 1.0, K = 1.0,
                    cross_section = 1.0)
       else layers
  s <- sink %||% list(Vd = 1.0)
  pk_arg <- pk %||% list(enabled = FALSE, thalf = 1)
  st <- sim_time %||% 30L
  skin_params(vehicle = v, layers = l, sink = s, pk = pk_arg,
              sim_time = st, ...)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

test_that("skin_params builds a list with class skin_params", {
  p <- make_minimal()
  expect_s3_class(p, "skin_params")
  expect_true(is.list(p))
  expect_named(p, c("sys", "log", "pk", "sink", "vehicle", "layers"))
})

test_that("skin_params accepts a data.frame for layers", {
  p <- make_minimal(layers = data.frame(
    name = c("SC", "DSL"),
    height = c(20L, 30L),
    D = c(1, 2),
    K = c(0.5, 0.1),
    cross_section = c(1.0, 0.5)
  ))
  expect_length(p$layers, 2L)
  expect_equal(p$layers[[1]]$name, "SC")
  expect_equal(p$layers[[2]]$cross_section, 0.5)
})

test_that("skin_params accepts a list-of-lists for layers", {
  p <- make_minimal(layers = list(
    list(name = "A", height = 10L, D = 1, K = 1, cross_section = 1)
  ))
  expect_length(p$layers, 1L)
  expect_equal(p$layers[[1]]$name, "A")
})

test_that("invalid layer cross_section is rejected", {
  expect_error(
    make_minimal(layers = data.frame(
      name = "x", height = 10L, D = 1, K = 1, cross_section = 1.5
    )),
    "cross_section"
  )
})

test_that("vehicle without app_area is rejected", {
  expect_error(make_minimal(vehicle = list(c_init = 1, height = 30L,
                                           D = 1, app_area = -1)),
               "app_area")
})

test_that("invalid sim_time is rejected", {
  expect_error(make_minimal(sim_time = 0L), "sim_time")
})

test_that("PK with non-positive thalf is rejected", {
  expect_error(
    make_minimal(pk = list(enabled = TRUE, thalf = 0)),
    "thalf"
  )
})

test_that("removing vehicle without layers is rejected", {
  expect_error(
    skin_params(
      vehicle = list(c_init = 1, height = 30L, D = 1, app_area = 1, remove_at = 10L),
      layers = list(),
      sink = list(Vd = 1.0),
      sim_time = 30L
    ),
    "layer"
  )
})

test_that("scaling values are normalized", {
  for (s in c("mg", "ug", "ng")) {
    p <- make_minimal(scaling = s)
    expect_equal(p$log$scaling, s)
  }
})

test_that("print method runs without error", {
  expect_output(print(make_minimal()), "skin_params")
})
