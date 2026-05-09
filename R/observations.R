#' Permeation observations (cumulative mass per area at the receptor)
#'
#' Wraps a tidy data.frame of permeation measurements (Franz cell, in-vitro
#' release, etc.) for use with [skin_fit()]. Permeation = mass that has
#' crossed the membrane and accumulated at the receptor side, normalised
#' by the application area.
#'
#' @param data A data.frame with required columns:
#'   * `time` (units of time)
#'   * `q_per_area` (units of mass/area, e.g. `ng_per_cm2(...)`)
#'
#'   Optional columns:
#'   * `sd` (units of mass/area; per-point standard deviation for
#'     variance weighting)
#'   * `subject` (character; grouping for multi-subject fits)
#'
#' @return A `permeation_obs` object (classed list) consumed by
#'   [skin_fit()]. Internal storage is in canonical units (`min`,
#'   `ng/cm^2`); the user-facing data.frame's units are converted on
#'   construction.
#' @export
permeation_obs <- function(data) {
  if (!is.data.frame(data)) {
    cli::cli_abort(c(
      "{.arg data} must be a {.cls data.frame}.",
      "x" = "Got {.obj_type_friendly {data}}."
    ))
  }
  required <- c("time", "q_per_area")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(c(
      "Permeation data is missing required column(s).",
      "x" = "Missing: {.val {missing_cols}}",
      "i" = "Required: {.val {required}}; optional: {.val sd}, {.val subject}."
    ))
  }
  if (nrow(data) == 0L) {
    cli::cli_abort("{.arg data} must have at least one row.")
  }

  time_min  <- .ensure_obs_units_col(data$time,       "min",      "time")
  q_canon   <- .ensure_obs_units_col(data$q_per_area, "ng/cm^2",  "q_per_area")
  if (any(q_canon < 0)) {
    cli::cli_abort(c(
      "{.field q_per_area} must be non-negative.",
      "x" = "Got negative value(s)."
    ))
  }

  sd_canon <- if ("sd" %in% names(data)) {
    sd_v <- .ensure_obs_units_col(data$sd, "ng/cm^2", "sd")
    if (any(sd_v < 0, na.rm = TRUE)) {
      cli::cli_abort("{.field sd} must be non-negative.")
    }
    sd_v
  } else {
    rep(NA_real_, length(time_min))
  }

  subject <- if ("subject" %in% names(data)) {
    s <- as.character(data$subject)
    if (any(is.na(s) | !nzchar(s))) {
      cli::cli_abort("{.field subject} values must be non-empty strings.")
    }
    s
  } else {
    rep("default", length(time_min))
  }

  out <- list(
    time_min          = time_min,
    q_per_area_ng_cm2 = q_canon,
    sd_ng_cm2         = sd_canon,
    subject           = subject,
    n                 = length(time_min)
  )
  class(out) <- c("permeation_obs", "list")
  out
}

#' Penetration observations (concentration averaged over a depth band)
#'
#' Wraps a tidy data.frame of penetration measurements (tape strip,
#' cryocut, etc.) for use with [skin_fit()]. Penetration = drug
#' concentration averaged over a defined depth band of the skin, in
#' the skin-surface = 0 convention (the top of the stratum corneum is
#' depth 0; the vehicle is not measured here).
#'
#' For point measurements (single depth, no band), set `depth_top`
#' equal to `depth_bottom`.
#'
#' @param data A data.frame with required columns:
#'   * `time` (units of time)
#'   * `depth_top`, `depth_bottom` (units of length; `depth_top <=
#'     depth_bottom`; both measured from the skin surface)
#'   * `concentration` (units of concentration; the strip's
#'     volume-averaged drug concentration)
#'
#'   Optional columns:
#'   * `sd` (units of concentration; per-point standard deviation)
#'   * `subject` (character; grouping for multi-subject fits)
#'
#' @return A `penetration_obs` object (classed list) consumed by
#'   [skin_fit()]. Internal storage is in canonical units (`min`,
#'   `um`, `ng/ml`).
#' @export
penetration_obs <- function(data) {
  if (!is.data.frame(data)) {
    cli::cli_abort(c(
      "{.arg data} must be a {.cls data.frame}.",
      "x" = "Got {.obj_type_friendly {data}}."
    ))
  }
  required <- c("time", "depth_top", "depth_bottom", "concentration")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(c(
      "Penetration data is missing required column(s).",
      "x" = "Missing: {.val {missing_cols}}",
      "i" = "Required: {.val {required}}; optional: {.val sd}, {.val subject}."
    ))
  }
  if (nrow(data) == 0L) {
    cli::cli_abort("{.arg data} must have at least one row.")
  }

  time_min     <- .ensure_obs_units_col(data$time,         "min",   "time")
  depth_top    <- .ensure_obs_units_col(data$depth_top,    "um",    "depth_top")
  depth_bottom <- .ensure_obs_units_col(data$depth_bottom, "um",    "depth_bottom")
  conc         <- .ensure_obs_units_col(data$concentration,"ng/ml", "concentration")

  if (any(depth_top > depth_bottom + 1e-9)) {
    cli::cli_abort(c(
      "{.field depth_top} must be <= {.field depth_bottom} for every row.",
      "i" = "For point measurements use equal values."
    ))
  }
  if (any(depth_top < 0)) {
    cli::cli_abort(c(
      "{.field depth_top} must be >= 0 (skin surface = 0).",
      "i" = "Penetration data is by definition skin-only; vehicle data is not yet supported."
    ))
  }
  if (any(conc < 0)) {
    cli::cli_abort("{.field concentration} must be non-negative.")
  }

  sd_canon <- if ("sd" %in% names(data)) {
    sd_v <- .ensure_obs_units_col(data$sd, "ng/ml", "sd")
    if (any(sd_v < 0, na.rm = TRUE)) {
      cli::cli_abort("{.field sd} must be non-negative.")
    }
    sd_v
  } else {
    rep(NA_real_, length(time_min))
  }

  subject <- if ("subject" %in% names(data)) {
    s <- as.character(data$subject)
    if (any(is.na(s) | !nzchar(s))) {
      cli::cli_abort("{.field subject} values must be non-empty strings.")
    }
    s
  } else {
    rep("default", length(time_min))
  }

  out <- list(
    time_min        = time_min,
    depth_top_um    = depth_top,
    depth_bottom_um = depth_bottom,
    conc_ng_ml      = conc,
    sd_ng_ml        = sd_canon,
    subject         = subject,
    n               = length(time_min)
  )
  class(out) <- c("penetration_obs", "list")
  out
}

# ---------- print methods ---------------------------------------------------

#' @export
print.permeation_obs <- function(x, ...) {
  n_subj <- length(unique(x$subject))
  cat(sprintf("<permeation_obs>  %d row(s), %d subject(s)\n", x$n, n_subj))
  cat(sprintf("  time range  : %g .. %g min\n",
              min(x$time_min), max(x$time_min)))
  cat(sprintf("  q range     : %g .. %g ng/cm^2\n",
              min(x$q_per_area_ng_cm2), max(x$q_per_area_ng_cm2)))
  cat(sprintf("  sd present  : %s\n", !all(is.na(x$sd_ng_cm2))))
  invisible(x)
}

#' @export
print.penetration_obs <- function(x, ...) {
  n_subj <- length(unique(x$subject))
  cat(sprintf("<penetration_obs> %d row(s), %d subject(s)\n", x$n, n_subj))
  cat(sprintf("  time range  : %g .. %g min\n",
              min(x$time_min), max(x$time_min)))
  cat(sprintf("  depth range : %g .. %g um\n",
              min(x$depth_top_um), max(x$depth_bottom_um)))
  cat(sprintf("  conc range  : %g .. %g ng/ml\n",
              min(x$conc_ng_ml), max(x$conc_ng_ml)))
  cat(sprintf("  sd present  : %s\n", !all(is.na(x$sd_ng_ml))))
  invisible(x)
}

# ---------- internal validators ---------------------------------------------

# Validate a single data.frame column: must be a units object that's
# convertible to `target_unit`. Returns the bare numeric vector in
# canonical units.
.ensure_obs_units_col <- function(x, target_unit, col_name,
                                  call = parent.frame()) {
  if (!inherits(x, "units")) {
    cli::cli_abort(c(
      "Column {.field {col_name}} must be units-aware.",
      "x" = "Got {.obj_type_friendly {x}}.",
      "i" = "Use a {.pkg skindiff} unit helper (e.g. {.code hours()}, {.code um()}, {.code ng_per_cm2()}, {.code ng_per_ml()}) or {.fn units::set_units}."
    ), call = call)
  }
  if (any(is.na(as.numeric(x)))) {
    cli::cli_abort("Column {.field {col_name}} contains NA values.",
                   call = call)
  }
  converted <- tryCatch(
    units::set_units(x, target_unit, mode = "standard"),
    error = function(e) {
      cli::cli_abort(c(
        "Column {.field {col_name}} has incompatible unit.",
        "x" = "Got unit {.val {units::deparse_unit(x)}}.",
        "i" = "Expected something convertible to {.val {target_unit}}."
      ), call = call, parent = e)
    }
  )
  as.numeric(converted)
}
