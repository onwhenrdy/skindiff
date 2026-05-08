#' Run a skindiff simulation
#'
#' @param params A `skin_params` object built with [skin_params()].
#' @param show_progress If `TRUE`, prints a textual progress indicator while
#'   the simulation runs. Defaults to `FALSE`.
#'
#' @return An object of class `"skin_result"` — a list containing:
#'
#'   * `status`:    `"executed"`, `"stopped"`, or `"failed"`.
#'   * `scaling`:   character; the mass unit reported in the result
#'                  (`"mg"`, `"ug"`, or `"ng"`).
#'   * `mass`:      data.frame with columns `time` (min) and one column per
#'                  logged compartment (vehicle, layers, sink), holding the
#'                  integrated mass in the chosen scaling unit.
#'   * `concentration`: data.frame with the same columns as `mass` but
#'                  holding the average concentration (mass per ml). The
#'                  vehicle column gives donor concentration vs time.
#'   * `cdp`:       named list, one entry per compartment with `log_cdp =
#'                  TRUE`. Each entry has `time`, `depth_um`, and a numeric
#'                  matrix `conc` indexed `[depth, time]`, in
#'                  `<scaling>/ml`.
#'   * `geometry`:  list with `min_step_um`, `max_step_um`, `eta` (NA for
#'                  the equidistant mesh), and `n_cells`.
#'   * `params`:    the input parameters (unchanged).
#'   * `runtime_s`: wall-clock runtime in seconds.
#'
#' @export
skin_simulate <- function(params, show_progress = FALSE) {
  if (!inherits(params, "skin_params")) {
    stop("`params` must be a skin_params object built with skin_params()",
         call. = FALSE)
  }
  show_progress <- .as_lgl(show_progress, "show_progress")

  t0 <- Sys.time()
  raw <- .cpp_simulate(unclass(params), show_progress = show_progress)
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  mass_df <- .mass_to_df(raw$mass)
  conc_df <- .mass_to_concentration(mass_df, params)

  out <- list(
    status        = raw$status,
    scaling       = raw$scaling,
    mass          = mass_df,
    concentration = conc_df,
    cdp           = raw$cdp,
    geometry      = raw$geometry,
    params        = params,
    runtime_s     = runtime
  )
  class(out) <- "skin_result"
  out
}

#' @export
print.skin_result <- function(x, ...) {
  cat("<skin_result>\n")
  cat(sprintf("  status      : %s\n", x$status))
  cat(sprintf("  scaling     : %s\n", x$scaling))
  cat(sprintf("  runtime     : %.3f s\n", x$runtime_s))
  cat(sprintf("  time points : %d (mass)\n", nrow(x$mass)))
  cdp_names <- names(x$cdp)
  if (length(cdp_names) > 0) {
    cat(sprintf("  cdp recorded: %s\n", paste(cdp_names, collapse = ", ")))
  } else {
    cat("  cdp recorded: <none>\n")
  }
  cat(sprintf("  geometry    : %d cells, min step %.4g um\n",
              x$geometry$n_cells, x$geometry$min_step_um))
  invisible(x)
}

#' @export
summary.skin_result <- function(object, ...) {
  unit <- object$scaling
  cat(sprintf("skindiff simulation summary (%s)\n", object$status))
  cat(sprintf("  scaling: %s\n", unit))
  cat(sprintf("  end time: %g min\n", max(object$mass$time)))
  cat("  final masses:\n")
  finals <- vapply(setdiff(names(object$mass), "time"),
                   function(nm) object$mass[[nm]][nrow(object$mass)], numeric(1))
  for (nm in names(finals)) {
    cat(sprintf("    %-20s %.6g %s\n", nm, finals[[nm]], unit))
  }
  invisible(object)
}

# ---------- internal: result reshaping ----------

.mass_to_df <- function(mass_list) {
  if (length(mass_list) == 0L) {
    return(data.frame(time = numeric(0)))
  }
  times <- mass_list[[1]]$time
  cols <- list(time = times)
  for (nm in names(mass_list)) {
    s <- mass_list[[nm]]
    if (length(s$time) != length(times) || !isTRUE(all.equal(s$time, times))) {
      stop(sprintf("Mass time vectors disagree for compartment '%s' — log intervals must match.",
                   nm), call. = FALSE)
    }
    cols[[nm]] <- s$value
  }
  as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
}

.mass_to_concentration <- function(mass_df, params) {
  if (nrow(mass_df) == 0L) return(mass_df)
  out <- data.frame(time = mass_df$time)
  vehicle_name <- params$vehicle$name
  vehicle_vol_ml <- params$vehicle$app_area * params$vehicle$height * 1e-4
  for (nm in setdiff(names(mass_df), "time")) {
    if (nm == vehicle_name) {
      out[[nm]] <- mass_df[[nm]] / vehicle_vol_ml
      next
    }
    if (nm == params$sink$name) {
      out[[nm]] <- mass_df[[nm]] / params$sink$Vd
      next
    }
    layer <- .find_layer(params$layers, nm)
    if (is.null(layer)) {
      out[[nm]] <- NA_real_
      next
    }
    layer_vol_ml <- params$vehicle$app_area * layer$cross_section *
                    layer$height * 1e-4
    out[[nm]] <- mass_df[[nm]] / layer_vol_ml
  }
  out
}

.find_layer <- function(layers, name) {
  for (l in layers) {
    if (identical(l$name, name)) return(l)
  }
  NULL
}
