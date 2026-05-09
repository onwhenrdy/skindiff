#' Unit helpers for skindiff
#'
#' Thin wrappers around [units::set_units()] that build values with the
#' right unit attached. Use these (or call `units::set_units()` directly)
#' anywhere a `skindiff` constructor expects a unit-bearing argument.
#'
#' Every constructor in `skindiff` rejects bare numerics for unit-bearing
#' arguments -- that is the whole point of unit safety. If you have a value
#' in a non-canonical unit (e.g. millimetres for a layer thickness), pass
#' `mm(0.11)` and the package converts internally to micrometres.
#'
#' Note that `minutes()` is preferred over `min()` to avoid masking
#' [base::min()].
#'
#' @param x A numeric vector.
#' @return A `units` object carrying the indicated unit.
#'
#' @examples
#' um(110)               # 110 [um]
#' mm(0.11)              # 0.11 [mm]   -- equivalent to um(110)
#' mg_per_ml(127.27)
#' hours(24)
#'
#' @name skin_units
NULL

#' @rdname skin_units
#' @export
um <- function(x) units::set_units(x, "um", mode = "standard")

#' @rdname skin_units
#' @export
mm <- function(x) units::set_units(x, "mm", mode = "standard")

#' @rdname skin_units
#' @export
cm <- function(x) units::set_units(x, "cm", mode = "standard")

#' @rdname skin_units
#' @export
cm2 <- function(x) units::set_units(x, "cm^2", mode = "standard")

#' @rdname skin_units
#' @export
mm2 <- function(x) units::set_units(x, "mm^2", mode = "standard")

#' @rdname skin_units
#' @export
ml <- function(x) units::set_units(x, "ml", mode = "standard")

#' @rdname skin_units
#' @export
mg_per_ml <- function(x) units::set_units(x, "mg/ml", mode = "standard")

#' @rdname skin_units
#' @export
ug_per_ml <- function(x) units::set_units(x, "ug/ml", mode = "standard")

#' @rdname skin_units
#' @export
ng_per_ml <- function(x) units::set_units(x, "ng/ml", mode = "standard")

#' @rdname skin_units
#' @export
um2_per_min <- function(x) units::set_units(x, "um^2/min", mode = "standard")

#' @rdname skin_units
#' @export
cm2_per_s <- function(x) units::set_units(x, "cm^2/s", mode = "standard")

#' @rdname skin_units
#' @export
seconds <- function(x) units::set_units(x, "s", mode = "standard")

#' @rdname skin_units
#' @export
minutes <- function(x) units::set_units(x, "min", mode = "standard")

#' @rdname skin_units
#' @export
hours <- function(x) units::set_units(x, "h", mode = "standard")

#' @rdname skin_units
#' @export
days <- function(x) units::set_units(x, "d", mode = "standard")

# Per-area mass helpers (for permeation observations: cumulative mass / area).
#' @rdname skin_units
#' @export
mg_per_cm2 <- function(x) units::set_units(x, "mg/cm^2", mode = "standard")

#' @rdname skin_units
#' @export
ug_per_cm2 <- function(x) units::set_units(x, "ug/cm^2", mode = "standard")

#' @rdname skin_units
#' @export
ng_per_cm2 <- function(x) units::set_units(x, "ng/cm^2", mode = "standard")


# ---------- internal validators ----------------------------------------------

# Convert a units-bearing argument to a bare numeric in `target_unit`. Errors
# if x is not a units object, or if x's dimension is incompatible with
# `target_unit`. The C++ engine receives bare numerics in canonical units, so
# we strip units at the boundary.
#
# `arg`  : human-readable argument name for error messages
# `call` : calling environment, used by cli to report the right frame
.ensure_units <- function(x, target_unit, arg, call = parent.frame()) {
  if (missing(x) || is.null(x)) {
    cli::cli_abort("{.arg {arg}} is required.", call = call)
  }
  if (!inherits(x, "units")) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a units-aware quantity.",
      "x" = "Got {.obj_type_friendly {x}}.",
      "i" = "Use a {.pkg skindiff} unit helper (e.g. {.code um(110)}, {.code mg_per_ml(1.2)}, {.code hours(24)}) or {.fn units::set_units} directly."
    ), call = call)
  }
  if (length(x) != 1L) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a single value.",
      "x" = "Got a length-{length(x)} units vector."
    ), call = call)
  }
  if (!is.finite(as.numeric(x))) {
    cli::cli_abort("{.arg {arg}} must be finite.", call = call)
  }
  converted <- tryCatch(
    units::set_units(x, target_unit, mode = "standard"),
    error = function(e) {
      cli::cli_abort(c(
        "{.arg {arg}} has an incompatible unit.",
        "x" = "Got {.val {format(x)}}.",
        "i" = "Expected something convertible to {.val {target_unit}}."
      ), call = call, parent = e)
    }
  )
  as.numeric(converted)
}

# As .ensure_units, plus optional min/max bounds (in the canonical unit).
.ensure_units_range <- function(x, target_unit, arg,
                                min = NULL, max = NULL,
                                exclusive_min = FALSE, exclusive_max = FALSE,
                                call = parent.frame()) {
  val <- .ensure_units(x, target_unit, arg, call = call)
  if (!is.null(min)) {
    bad <- if (exclusive_min) val <= min else val < min
    if (bad) {
      cli::cli_abort(c(
        "{.arg {arg}} is out of range.",
        "x" = "Got {.val {val}} [{target_unit}], must be {if (exclusive_min) '>' else '>='} {.val {min}} [{target_unit}]."
      ), call = call)
    }
  }
  if (!is.null(max)) {
    bad <- if (exclusive_max) val >= max else val > max
    if (bad) {
      cli::cli_abort(c(
        "{.arg {arg}} is out of range.",
        "x" = "Got {.val {val}} [{target_unit}], must be {if (exclusive_max) '<' else '<='} {.val {max}} [{target_unit}]."
      ), call = call)
    }
  }
  val
}

# Same as .ensure_units, then round to integer in the canonical unit. Used
# for fields like `height` where the engine wants integer micrometres.
.ensure_units_int <- function(x, target_unit, arg, min = NULL,
                              call = parent.frame()) {
  val <- .ensure_units(x, target_unit, arg, call = call)
  rounded <- as.integer(round(val))
  if (abs(val - rounded) > 1e-9) {
    cli::cli_inform(c(
      "i" = "{.arg {arg}}: rounding {.val {val}} [{target_unit}] to {.val {rounded}} [{target_unit}]."
    ))
  }
  if (!is.null(min) && rounded < min) {
    cli::cli_abort(c(
      "{.arg {arg}} is out of range.",
      "x" = "Got {.val {format(x)}} (= {.val {rounded}} [{target_unit}] after conversion), must be at least {.val {min}} [{target_unit}]."
    ), call = call)
  }
  rounded
}

# A duration, in minutes, where NULL means "disabled" (returns 0L).
.ensure_duration_or_null <- function(x, arg, call = parent.frame()) {
  if (is.null(x)) return(0L)
  .ensure_units_int(x, "min", arg, min = 1L, call = call)
}

# A bare-numeric guard with cli errors.
.ensure_dimensionless <- function(x, arg, min = NULL, max = NULL,
                                  exclusive_min = FALSE, exclusive_max = FALSE,
                                  call = parent.frame()) {
  if (missing(x) || length(x) != 1L || !is.numeric(x) || is.na(x)) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a single numeric value.",
      "x" = "Got {.obj_type_friendly {x}}."
    ), call = call)
  }
  if (inherits(x, "units")) {
    cli::cli_abort(c(
      "{.arg {arg}} is dimensionless and must be a bare numeric.",
      "x" = "Got {.val {format(x)}} (units vector)."
    ), call = call)
  }
  val <- as.numeric(x)
  if (!is.null(min)) {
    bad <- if (exclusive_min) val <= min else val < min
    if (bad) cli::cli_abort(c(
      "{.arg {arg}} is out of range.",
      "x" = "Got {.val {val}}, must be {if (exclusive_min) '>' else '>='} {.val {min}}."
    ), call = call)
  }
  if (!is.null(max)) {
    bad <- if (exclusive_max) val >= max else val > max
    if (bad) cli::cli_abort(c(
      "{.arg {arg}} is out of range.",
      "x" = "Got {.val {val}}, must be {if (exclusive_max) '<' else '<='} {.val {max}}."
    ), call = call)
  }
  val
}

.ensure_int <- function(x, arg, min = NULL, max = NULL, call = parent.frame()) {
  v <- .ensure_dimensionless(x, arg, call = call)
  if (abs(v - round(v)) > 1e-9) {
    cli::cli_abort(c(
      "{.arg {arg}} must be an integer.",
      "x" = "Got {.val {v}}."
    ), call = call)
  }
  vi <- as.integer(round(v))
  if (!is.null(min) && vi < min) {
    cli::cli_abort(c(
      "{.arg {arg}} is out of range.",
      "x" = "Got {.val {vi}}, must be at least {.val {min}}."
    ), call = call)
  }
  if (!is.null(max) && vi > max) {
    cli::cli_abort(c(
      "{.arg {arg}} is out of range.",
      "x" = "Got {.val {vi}}, must be at most {.val {max}}."
    ), call = call)
  }
  vi
}

.ensure_lgl <- function(x, arg, call = parent.frame()) {
  if (missing(x) || length(x) != 1L || !is.logical(x) || is.na(x)) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a single TRUE/FALSE.",
      "x" = "Got {.obj_type_friendly {x}}."
    ), call = call)
  }
  as.logical(x)
}

.ensure_chr <- function(x, arg, call = parent.frame()) {
  if (missing(x) || length(x) != 1L || !is.character(x) || is.na(x) || !nzchar(x)) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a single non-empty string.",
      "x" = "Got {.obj_type_friendly {x}}."
    ), call = call)
  }
  as.character(x)
}
