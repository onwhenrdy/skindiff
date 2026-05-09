#' Run a skindiff simulation
#'
#' @param params A `skin_params` object built with [skin_params()].
#' @param show_progress If `TRUE`, prints a textual progress indicator while
#'   the simulation runs. Defaults to `FALSE`.
#'
#' @return An object of class `"skin_result"` -- a list containing:
#'
#'   * `status`:    `"executed"`, `"stopped"`, or `"failed"`.
#'   * `scaling`:   character; the mass unit reported in the result
#'                  (`"mg"`, `"ug"`, or `"ng"`).
#'   * `mass`:      data.frame with columns `time` (units of time) and one
#'                  column per logged compartment (vehicle, layers, sink),
#'                  each carrying the integrated mass with its scaling unit.
#'   * `concentration`: data.frame with the same columns as `mass`, holding
#'                  the average concentration (mass per ml) with units. The
#'                  vehicle column gives donor concentration vs time; the
#'                  sink column is `NA` for perfect sinks (use the cumulative
#'                  mass column instead).
#'   * `cdp`:       named list, one entry per compartment with `log_cdp =
#'                  TRUE`. Each entry has `time` (units of time), `depth`
#'                  (units of length), and a numeric matrix `conc` indexed
#'                  `[depth, time]` carrying its scaling/ml unit.
#'   * `geometry`:  list with `min_step` (units of length), `max_step`,
#'                  and `n_cells` (bare integer).
#'   * `params`:    the input parameters (unchanged).
#'   * `runtime`:   wall-clock runtime, units of time.
#'
#' @export
skin_simulate <- function(params, show_progress = FALSE) {
  if (!inherits(params, "skin_params")) {
    cli::cli_abort(c(
      "{.arg params} must be a {.cls skin_params} object.",
      "i" = "Build it with {.fn skin_params}."
    ))
  }
  show_progress <- .ensure_lgl(show_progress, "show_progress")

  t0 <- Sys.time()
  raw <- .cpp_simulate(unclass(params), show_progress = show_progress)
  runtime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  scaling_unit <- raw$scaling                 # "mg" / "ug" / "ng"
  conc_unit    <- paste0(scaling_unit, "/ml")

  mass_df <- .mass_to_df(raw$mass, scaling_unit)
  conc_df <- .mass_to_concentration(mass_df, params, conc_unit)
  cdp     <- .cdp_with_units(raw$cdp, conc_unit)

  geometry <- list(
    min_step = units::set_units(raw$geometry$min_step_um, "um"),
    max_step = units::set_units(raw$geometry$max_step_um, "um"),
    n_cells  = raw$geometry$n_cells
  )

  out <- list(
    status        = raw$status,
    scaling       = scaling_unit,
    mass          = mass_df,
    concentration = conc_df,
    cdp           = cdp,
    geometry      = geometry,
    params        = params,
    runtime       = units::set_units(runtime_s, "s")
  )
  class(out) <- "skin_result"
  out
}

#' @export
print.skin_result <- function(x, ...) {
  cat("<skin_result>\n")
  cat(sprintf("  status      : %s\n", x$status))
  cat(sprintf("  scaling     : %s\n", x$scaling))
  cat(sprintf("  runtime     : %s\n", format(x$runtime)))
  cat(sprintf("  time points : %d (mass)\n", nrow(x$mass)))
  cdp_names <- names(x$cdp)
  if (length(cdp_names) > 0) {
    cat(sprintf("  cdp recorded: %s\n", paste(cdp_names, collapse = ", ")))
  } else {
    cat("  cdp recorded: <none>\n")
  }
  cat(sprintf("  geometry    : %d cells, min step %s\n",
              x$geometry$n_cells, format(x$geometry$min_step)))
  invisible(x)
}

#' @export
summary.skin_result <- function(object, ...) {
  cat(sprintf("skindiff simulation summary (%s)\n", object$status))
  cat(sprintf("  scaling: %s\n", object$scaling))
  cat(sprintf("  end time: %s\n",
              format(units::set_units(max(as.numeric(object$mass$time)), "min"))))
  cat("  final masses:\n")
  for (nm in setdiff(names(object$mass), "time")) {
    final <- object$mass[[nm]][nrow(object$mass)]
    cat(sprintf("    %-20s %s\n", nm, format(final)))
  }
  invisible(object)
}

# ---------- internal: result reshaping ----------

.mass_to_df <- function(mass_list, scaling_unit) {
  if (length(mass_list) == 0L) {
    return(data.frame(time = units::set_units(numeric(0), "min")))
  }
  # Different series may have different lengths if a compartment was
  # removed mid-run (the donor's series ends at remove_at, the layers and
  # sink continue). Use the longest series as the canonical time grid and
  # pad shorter ones with NA.
  lengths <- vapply(mass_list, function(s) length(s$time), integer(1))
  canonical <- mass_list[[which.max(lengths)]]$time
  cols <- list(time = units::set_units(canonical, "min"))
  for (nm in names(mass_list)) {
    s <- mass_list[[nm]]
    matched <- match(canonical, s$time)
    cols[[nm]] <- units::set_units(s$value[matched], scaling_unit,
                                   mode = "standard")
  }
  data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
}

.mass_to_concentration <- function(mass_df, params, conc_unit) {
  if (nrow(mass_df) == 0L) return(mass_df)
  out <- data.frame(time = mass_df$time)
  vehicle_name    <- params$vehicle$name
  area_cm2        <- params$.meta$area_cm2
  vehicle_vol_ml  <- area_cm2 * params$vehicle$height * 1e-4   # cm^2 * um * 1e-4 = ml
  sink_is_perfect <- isTRUE(params$.meta$sink_is_perfect)

  na_col <- units::set_units(rep(NA_real_, nrow(mass_df)), conc_unit,
                             mode = "standard")

  for (nm in setdiff(names(mass_df), "time")) {
    out[[nm]] <- if (nm == vehicle_name) {
      mass_df[[nm]] / units::set_units(vehicle_vol_ml, "ml")
    } else if (nm == params$sink$name) {
      # Perfect sink: mass/Vd ~ 0 would be misleading. Flag with NA.
      if (sink_is_perfect) na_col
      else mass_df[[nm]] / units::set_units(params$sink$Vd, "ml")
    } else {
      layer <- .find_layer(params$layers, nm)
      if (is.null(layer)) {
        na_col
      } else {
        layer_vol_ml <- area_cm2 * layer$cross_section * layer$height * 1e-4
        mass_df[[nm]] / units::set_units(layer_vol_ml, "ml")
      }
    }
  }
  out
}

.cdp_with_units <- function(cdp, conc_unit) {
  for (nm in names(cdp)) {
    s <- cdp[[nm]]
    cdp[[nm]] <- list(
      time  = units::set_units(s$time, "min"),
      depth = units::set_units(s$depth_um, "um"),
      conc  = units::set_units(s$conc, conc_unit, mode = "standard")
    )
  }
  cdp
}

.find_layer <- function(layers, name) {
  for (l in layers) {
    if (identical(l$name, name)) return(l)
  }
  NULL
}
