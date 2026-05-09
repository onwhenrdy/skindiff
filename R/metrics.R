#' Cumulative permeated mass per unit area at the receptor
#'
#' Returns the standard Franz-cell output curve: cumulative drug that has
#' crossed the membrane per unit application area, vs time.
#'
#' @param res A `skin_result` object.
#' @return A data.frame with columns `time` (units of time) and `Q`
#'   (units of mass per area).
#' @export
permeated <- function(res) {
  if (!inherits(res, "skin_result")) {
    cli::cli_abort("{.arg res} must be a {.cls skin_result} object.")
  }
  sink_name <- res$params$sink$name
  if (!sink_name %in% names(res$mass)) {
    cli::cli_abort(c(
      "Sink mass was not logged.",
      "i" = "Set {.code log_mass = TRUE} on the sink in your call."
    ))
  }
  area <- units::set_units(res$params$.meta$area_cm2, "cm^2")
  data.frame(
    time = res$mass$time,
    Q    = res$mass[[sink_name]] / area
  )
}

#' Permeated mass per unit area at specific time(s)
#'
#' Linearly interpolates the cumulative permeated curve [permeated()] to
#' the requested time(s).
#'
#' @param res A `skin_result` object.
#' @param t Time(s) at which to evaluate. Must have units of time.
#' @return A units vector of permeated amounts in mass-per-area.
#' @export
permeated_at <- function(res, t) {
  perm <- permeated(res)
  t_min <- .ensure_units_vec_min(t, "t")
  t_grid <- as.numeric(perm$time)
  Q_unit <- units::deparse_unit(perm$Q)
  q_at   <- stats::approx(t_grid, as.numeric(perm$Q),
                          xout = t_min, rule = 2)$y
  units::set_units(q_at, Q_unit, mode = "standard")
}

#' Flux at the receptor
#'
#' Numerical derivative of the cumulative permeated curve. Returns flux at
#' the midpoints of consecutive sample times.
#'
#' @param res A `skin_result` object.
#' @return A data.frame with columns `time` (midpoint times, units of
#'   time) and `flux` (units of mass per area per time).
#' @export
flux <- function(res) {
  perm <- permeated(res)
  if (nrow(perm) < 2) {
    cli::cli_abort("Need at least two sample points to compute flux.")
  }
  t_mid <- (utils::head(perm$time, -1) + utils::tail(perm$time, -1)) / 2
  flux_vals <- diff(perm$Q) / diff(perm$time)
  data.frame(time = t_mid, flux = flux_vals)
}

#' Concentration-depth profile slice at a specific time
#'
#' For each compartment that was logged with `log_cdp = TRUE`, return a
#' depth-vs-concentration data.frame interpolated to the requested time.
#'
#' @param res A `skin_result` object.
#' @param t Time at which to evaluate (single value, units of time).
#' @return A list with one data.frame (`depth`, `conc`) per logged
#'   compartment.
#' @export
profile_at <- function(res, t) {
  if (!inherits(res, "skin_result")) {
    cli::cli_abort("{.arg res} must be a {.cls skin_result} object.")
  }
  if (length(res$cdp) == 0L) {
    cli::cli_abort(c(
      "No concentration-depth profiles were logged.",
      "i" = "Set {.code log_cdp = TRUE} on the compartments you want to inspect."
    ))
  }
  t_min <- .ensure_units_scalar_min(t, "t")
  out <- list()
  for (nm in names(res$cdp)) {
    s <- res$cdp[[nm]]
    t_grid <- as.numeric(s$time)
    conc_at <- apply(unclass(s$conc), 1, function(row) {
      stats::approx(t_grid, row, xout = t_min, rule = 2)$y
    })
    conc_unit <- units::deparse_unit(s$conc)
    out[[nm]] <- data.frame(
      depth = s$depth,
      conc  = units::set_units(conc_at, conc_unit, mode = "standard")
    )
  }
  out
}

#' Standard skin-permeation metrics for a run
#'
#' @description
#' Computes a one-row summary of standard skin-permeation outputs from a
#' completed simulation: steady-state flux, lag time, permeability
#' coefficient, cumulative permeated amount, donor half-life, and (for
#' finite-Vd sinks) receptor AUC, C_max, and t_max.
#'
#' @details
#' Reported columns:
#' \itemize{
#'   \item `J_ss`       : steady-state flux (mass / area / hour)
#'   \item `t_lag`      : lag time, x-intercept of the J_ss tangent (hours)
#'   \item `K_p`        : permeability coefficient `J_ss / c_donor` (cm/h)
#'   \item `Q_total`    : cumulative permeated at end of run (mass / area)
#'   \item `r2_ss`      : R-squared of the steady-state line fit
#'   \item `t_50_donor` : time for vehicle mass to halve, NA if never reached
#'   \item `AUC_sink`   : integral of receptor concentration, NA for perfect sinks
#'   \item `C_max_sink` : peak receptor concentration, NA for perfect sinks
#'   \item `t_max_sink` : time of `C_max_sink`, NA for perfect sinks
#' }
#'
#' For a finite-dose run, `K_p` is reported but only physically meaningful
#' if the donor concentration stayed roughly constant. If the donor
#' depleted by more than `depletion_warn` (default 10\%), a warning is
#' emitted.
#'
#' @param res A `skin_result` object.
#' @param ss_window Time window for the steady-state line fit. Either a
#'   length-2 numeric vector of fractions in `[0, 1]` of the total
#'   duration (default `c(0.7, 1.0)` for the last 30\%), or a length-2
#'   units-of-time vector (e.g. `c(hours(12), hours(24))`).
#' @param depletion_warn Threshold above which a finite-dose donor's
#'   end-of-run depletion triggers a `K_p` validity warning. Default 0.1
#'   (10\%).
#' @return A one-row data.frame with units-bearing columns.
#' @export
metrics <- function(res, ss_window = c(0.7, 1.0), depletion_warn = 0.1) {
  if (!inherits(res, "skin_result")) {
    cli::cli_abort("{.arg res} must be a {.cls skin_result} object.")
  }

  perm <- permeated(res)
  t_grid_min <- as.numeric(perm$time)             # in minutes
  Q_vals     <- as.numeric(perm$Q)
  Q_unit     <- units::deparse_unit(perm$Q)        # e.g. "ng/cm^2"

  # ---- steady-state window ----
  win_min <- .ss_window_to_minutes(ss_window, t_grid_min)
  in_window <- t_grid_min >= win_min[1] & t_grid_min <= win_min[2]
  if (sum(in_window) < 2L) {
    cli::cli_abort(c(
      "Steady-state window contains fewer than 2 sample points.",
      "i" = "Got {.val {sum(in_window)}} point(s) in [{.val {win_min[1]}}, {.val {win_min[2]}}] min.",
      "i" = "Run for longer or widen {.arg ss_window}."
    ))
  }

  # ---- linear fit Q ~ t in window ----
  fit   <- stats::lm(Q_vals[in_window] ~ t_grid_min[in_window])
  slope <- unname(stats::coef(fit)[2L])    # Q-units per minute
  inter <- unname(stats::coef(fit)[1L])
  r2    <- summary(fit)$r.squared

  # J_ss: flux per hour (slope is per minute -> *60)
  J_ss <- units::set_units(slope * 60, paste0(Q_unit, "/h"),
                           mode = "standard")
  # t_lag: -intercept/slope, then convert to hours
  t_lag_min <- if (slope > 0) -inter / slope else NA_real_
  t_lag <- units::set_units(t_lag_min / 60, "h", mode = "standard")

  # ---- K_p = J_ss / c_donor ----
  c_donor_mg_per_ml <- res$params$vehicle$c_init     # canonical
  scaling_unit      <- res$params$log$scaling
  if (c_donor_mg_per_ml > 0) {
    c_donor <- units::set_units(
      units::set_units(c_donor_mg_per_ml, "mg/ml"),
      paste0(scaling_unit, "/ml"),
      mode = "standard"
    )
    K_p <- units::set_units(J_ss / c_donor, "cm/h", mode = "standard")
  } else {
    K_p <- units::set_units(NA_real_, "cm/h", mode = "standard")
  }

  # K_p validity check for finite dose
  vehicle_name <- res$params$vehicle$name
  vehicle_mass <- res$mass[[vehicle_name]]
  if (!is.null(vehicle_mass) && length(vehicle_mass) > 0L) {
    v       <- as.numeric(vehicle_mass)
    initial <- v[1]
    finite_idx <- which(!is.na(v))
    if (length(finite_idx) > 0L && initial > 0) {
      final <- v[utils::tail(finite_idx, 1L)]
      depletion <- (initial - final) / initial
      if (isTRUE(res$params$vehicle$finite_dose) && depletion > depletion_warn) {
        cli::cli_warn(c(
          "{.val K_p} may be unreliable: finite-dose donor depleted by {.val {round(depletion * 100, 1)}}% by end-of-run.",
          "i" = "K_p assumes a roughly constant donor concentration."
        ))
      }
    }
  }

  # ---- Q_total ----
  Q_total <- utils::tail(perm$Q, 1L)

  # ---- t_50_donor ----
  t_50_donor <- units::set_units(NA_real_, "h", mode = "standard")
  if (!is.null(vehicle_mass) && length(vehicle_mass) > 0L) {
    v       <- as.numeric(vehicle_mass)
    initial <- v[1]
    if (initial > 0) {
      threshold <- initial / 2
      below_idx <- which(!is.na(v) & v <= threshold)
      if (length(below_idx) > 0L) {
        first <- below_idx[1]
        t_50_min <- if (first == 1L) {
          0
        } else {
          t1 <- t_grid_min[first - 1L]; t2 <- t_grid_min[first]
          v1 <- v[first - 1L];          v2 <- v[first]
          t1 + (threshold - v1) / (v2 - v1) * (t2 - t1)
        }
        t_50_donor <- units::set_units(t_50_min / 60, "h", mode = "standard")
      }
    }
  }

  # ---- AUC, C_max, t_max for finite sink ----
  sink_is_perfect <- isTRUE(res$params$.meta$sink_is_perfect)
  AUC_sink   <- units::set_units(NA_real_, paste0(scaling_unit, "*h/ml"),
                                 mode = "standard")
  C_max_sink <- units::set_units(NA_real_, paste0(scaling_unit, "/ml"),
                                 mode = "standard")
  t_max_sink <- units::set_units(NA_real_, "h", mode = "standard")

  if (!sink_is_perfect) {
    sink_name <- res$params$sink$name
    c_sink    <- res$concentration[[sink_name]]
    if (!is.null(c_sink) && nrow(perm) >= 2L) {
      cs <- as.numeric(c_sink)
      C_unit <- units::deparse_unit(c_sink)
      auc_min <- sum((cs[-1L] + cs[-length(cs)]) / 2 * diff(t_grid_min))
      AUC_sink <- units::set_units(
        units::set_units(auc_min, paste0(C_unit, "*min"), mode = "standard"),
        paste0(C_unit, "*h"),
        mode = "standard"
      )
      i_max <- which.max(cs)
      C_max_sink <- units::set_units(cs[i_max], C_unit, mode = "standard")
      t_max_sink <- units::set_units(t_grid_min[i_max] / 60, "h",
                                     mode = "standard")
    }
  }

  data.frame(
    J_ss        = J_ss,
    t_lag       = t_lag,
    K_p         = K_p,
    Q_total     = Q_total,
    r2_ss       = r2,
    t_50_donor  = t_50_donor,
    AUC_sink    = AUC_sink,
    C_max_sink  = C_max_sink,
    t_max_sink  = t_max_sink
  )
}

# ---------- internal helpers ----------

.ss_window_to_minutes <- function(ss_window, t_grid_min,
                                  call = parent.frame()) {
  if (length(ss_window) != 2L) {
    cli::cli_abort("{.arg ss_window} must have exactly 2 elements.",
                   call = call)
  }
  if (inherits(ss_window, "units")) {
    converted <- tryCatch(
      units::set_units(ss_window, "min", mode = "standard"),
      error = function(e) {
        cli::cli_abort(c(
          "{.arg ss_window} has incompatible unit.",
          "x" = "Got {.val {format(ss_window)}}.",
          "i" = "Expected something convertible to {.val min}."
        ), call = call, parent = e)
      }
    )
    win <- as.numeric(converted)
  } else {
    if (!is.numeric(ss_window)) {
      cli::cli_abort(
        "{.arg ss_window} must be numeric or a units-of-time vector.",
        call = call
      )
    }
    if (any(ss_window < 0) || any(ss_window > 1)) {
      cli::cli_abort(c(
        "{.arg ss_window} fractions must be in {.val [0, 1]}.",
        "x" = "Got {.val {ss_window}}."
      ), call = call)
    }
    duration_min <- max(t_grid_min)
    win <- ss_window * duration_min
  }
  if (win[2L] <= win[1L]) {
    cli::cli_abort(c(
      "{.arg ss_window} must be strictly increasing.",
      "x" = "Got [{.val {win[1L]}}, {.val {win[2L]}}] min."
    ), call = call)
  }
  win
}

.ensure_units_vec_min <- function(x, arg, call = parent.frame()) {
  if (!inherits(x, "units")) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a units-of-time vector.",
      "i" = "Use a {.pkg skindiff} helper like {.code minutes(...)} or {.code hours(...)}."
    ), call = call)
  }
  converted <- tryCatch(
    units::set_units(x, "min", mode = "standard"),
    error = function(e) {
      cli::cli_abort(c(
        "{.arg {arg}} has incompatible unit.",
        "x" = "Got {.val {format(x)}}.",
        "i" = "Expected something convertible to {.val min}."
      ), call = call, parent = e)
    }
  )
  as.numeric(converted)
}

.ensure_units_scalar_min <- function(x, arg, call = parent.frame()) {
  v <- .ensure_units_vec_min(x, arg, call = call)
  if (length(v) != 1L) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a single value.",
      "x" = "Got length {.val {length(v)}}."
    ), call = call)
  }
  v
}
