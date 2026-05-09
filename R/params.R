#' Build a vehicle (donor) compartment
#'
#' Describes the donor formulation that sits on top of the skin. Every
#' unit-bearing argument requires a `units` object built with
#' [units::set_units()] or one of the `skindiff` helpers (e.g. [um()],
#' [mg_per_ml()]). Bare numerics are rejected.
#'
#' @param c_init Initial donor concentration. Must have units of
#'   concentration (e.g. `mg_per_ml(127.27)`, `ng_per_ml(50)`).
#' @param height Donor thickness. Must have units of length (e.g.
#'   `um(110)`, `mm(0.11)`). Rounded to integer micrometres internally.
#' @param D Diffusion coefficient inside the donor. Must have units of
#'   area-per-time (e.g. `um2_per_min(9.27)`, `cm2_per_s(1e-7)`).
#' @param finite_dose If `FALSE`, the donor concentration is clamped at
#'   `c_init` for the entire run (true Dirichlet boundary at the
#'   donor/skin interface).
#' @param replace_after Donor refresh period (units of time, e.g.
#'   `hours(24)`), or `NULL` to disable. When set, every donor cell is
#'   reset to `c_init` at every multiple of this duration.
#' @param remove_at Donor removal time (units of time), or `NULL` to
#'   disable. When set, the donor is deleted at this time and the
#'   simulation continues on (skin layers + sink).
#' @param name Compartment label used in result output.
#' @param log_mass,log_cdp Whether to record the mass time-series and/or
#'   the concentration-depth profile for this compartment.
#'
#' @return A `skin_vehicle` object (a classed list) ready for [skin_params()].
#' @export
vehicle <- function(c_init,
                    height,
                    D,
                    finite_dose   = TRUE,
                    replace_after = NULL,
                    remove_at     = NULL,
                    name          = "Vehicle",
                    log_mass      = TRUE,
                    log_cdp       = FALSE) {
  out <- list(
    name              = .ensure_chr(name, "name"),
    c_init_mg_per_ml  = .ensure_units_range(c_init, "mg/ml", "c_init",
                                            min = 0),
    height_um         = .ensure_units_int(height, "um", "height", min = 3L),
    D_um2_per_min     = .ensure_units_range(D, "um^2/min", "D", min = 0),
    finite_dose       = .ensure_lgl(finite_dose, "finite_dose"),
    replace_after_min = .ensure_duration_or_null(replace_after, "replace_after"),
    remove_at_min     = .ensure_duration_or_null(remove_at, "remove_at"),
    log_mass          = .ensure_lgl(log_mass, "log_mass"),
    log_cdp           = .ensure_lgl(log_cdp,  "log_cdp")
  )
  class(out) <- c("skin_vehicle", "list")
  out
}

#' Build a skin-layer compartment
#'
#' Layers are passed top-to-bottom (closest to vehicle first). The layer's
#' physical area is `area * cross_section`, where `area` is set on
#' [skin_params()].
#'
#' @param name Layer label (e.g. `"Stratum corneum"`).
#' @param height Layer thickness. Must have units of length.
#' @param D Diffusion coefficient. Must have units of area-per-time.
#' @param K Partition coefficient relative to the vehicle (dimensionless,
#'   strictly positive). At equilibrium, `c_layer = K * c_vehicle`.
#' @param cross_section Fraction of the application area occupied by this
#'   layer (dimensionless, in `(0, 1]`).
#' @param c_init Initial concentration in this layer. Must have units of
#'   concentration. Defaults to `mg_per_ml(0)`. Almost always 0.
#' @param log_mass,log_cdp Whether to record the mass time-series and/or
#'   the concentration-depth profile for this compartment.
#'
#' @return A `skin_layer` object (a classed list) ready for [skin_params()].
#' @export
layer <- function(name,
                  height,
                  D,
                  K,
                  cross_section,
                  c_init   = mg_per_ml(0),
                  log_mass = TRUE,
                  log_cdp  = FALSE) {
  out <- list(
    name              = .ensure_chr(name, "name"),
    height_um         = .ensure_units_int(height, "um", "height", min = 3L),
    D_um2_per_min     = .ensure_units_range(D, "um^2/min", "D", min = 0),
    K                 = .ensure_dimensionless(K, "K",
                                              min = 0, exclusive_min = TRUE),
    cross_section     = .ensure_dimensionless(cross_section, "cross_section",
                                              min = 0, max = 1,
                                              exclusive_min = TRUE),
    c_init_mg_per_ml  = .ensure_units_range(c_init, "mg/ml", "c_init", min = 0),
    log_mass          = .ensure_lgl(log_mass, "log_mass"),
    log_cdp           = .ensure_lgl(log_cdp,  "log_cdp")
  )
  class(out) <- c("skin_layer", "list")
  out
}

#' Build a perfect-sink receptor
#'
#' A perfect sink has effectively infinite volume, so the receptor
#' concentration stays at zero and incoming mass simply accumulates. This
#' is the right model for typical Franz-cell experiments where the
#' receptor is large and well-stirred. The reported sink concentration
#' (`mass / volume`) is `NA` for a perfect sink because the value is not
#' physically meaningful; use the cumulative mass column instead.
#'
#' @param name Compartment label used in result output.
#' @param log_mass Whether to record the mass time-series.
#'
#' @return A `skin_sink` object (a classed list) ready for [skin_params()].
#' @export
perfect_sink <- function(name = "Sink", log_mass = TRUE) {
  out <- list(
    name             = .ensure_chr(name, "name"),
    type             = "perfect",
    Vd_ml            = 1.0e9,   # large enough that c_sink stays near zero
    c_init_mg_per_ml = 0.0,
    log_mass         = .ensure_lgl(log_mass, "log_mass")
  )
  class(out) <- c("skin_sink", "list")
  out
}

#' Build a finite-volume receptor
#'
#' A finite sink has a fixed volume of distribution `Vd`. The receptor
#' concentration `mass / Vd` is meaningful and reported. The sink does
#' *not* back-diffuse into the membrane: the membrane sees a perfect
#' Dirichlet boundary regardless of `Vd`, so the only physical effect of
#' `Vd` is on the *reported* sink concentration.
#'
#' @param name Compartment label used in result output.
#' @param Vd Volume of distribution (units of volume, strictly positive).
#' @param c_init Initial receptor concentration (units of concentration).
#'   Defaults to `mg_per_ml(0)`. Almost always 0.
#' @param log_mass Whether to record the mass time-series.
#'
#' @return A `skin_sink` object (a classed list) ready for [skin_params()].
#' @export
finite_sink <- function(name, Vd, c_init = mg_per_ml(0), log_mass = TRUE) {
  out <- list(
    name             = .ensure_chr(name, "name"),
    type             = "finite",
    Vd_ml            = .ensure_units_range(Vd, "ml", "Vd",
                                           min = 0, exclusive_min = TRUE),
    c_init_mg_per_ml = .ensure_units_range(c_init, "mg/ml", "c_init", min = 0),
    log_mass         = .ensure_lgl(log_mass, "log_mass")
  )
  class(out) <- c("skin_sink", "list")
  out
}

#' Build a validated set of skindiff simulation parameters
#'
#' Composes the experimental geometry (`area`), the donor (`vehicle`), a
#' top-to-bottom stack of skin layers, and a receptor sink into a single
#' validated parameter object that can be passed to [skin_simulate()].
#'
#' All compartment-level quantities live on the constructor objects
#' ([vehicle()], [layer()], [perfect_sink()], [finite_sink()]). All
#' system-level quantities (application area, simulation duration, log
#' intervals, mesh resolution, output scaling) live here.
#'
#' Unit-bearing arguments require `units` objects (use [um()],
#' [mg_per_ml()], [hours()], etc., or call [units::set_units()]
#' directly).
#'
#' @param area Application area (skin patch area). Must have units of area.
#' @param vehicle A [vehicle()] object.
#' @param layers A list of [layer()] objects, ordered top-to-bottom.
#' @param sink A [perfect_sink()] or [finite_sink()] object.
#' @param duration Total simulation time (units of time, integer minutes
#'   internally).
#' @param resolution Mesh refinement: integer cells per micrometre in the
#'   smallest-D compartment (dimensionless integer >= 1). Higher-D
#'   compartments get proportionally coarser cells.
#' @param max_module Stability target for the implicit sub-step count
#'   (dimensionless, > 0).
#' @param scaling Output mass scaling: `"mg"`, `"ug"`, or `"ng"`.
#' @param mass_log_interval Sample interval for mass time-series (units
#'   of time, integer minutes internally).
#' @param cdp_log_interval Sample interval for concentration-depth
#'   profiles (units of time, integer minutes internally).
#'
#' @return A `skin_params` object ready for [skin_simulate()].
#' @export
skin_params <- function(area,
                        vehicle,
                        layers,
                        sink,
                        duration,
                        resolution        = 1L,
                        max_module        = 50,
                        scaling           = c("mg", "ug", "ng"),
                        mass_log_interval = minutes(1L),
                        cdp_log_interval  = minutes(1L)) {
  if (missing(vehicle) || !inherits(vehicle, "skin_vehicle")) {
    cli::cli_abort(c(
      "{.arg vehicle} must be a {.cls skin_vehicle} object.",
      "i" = "Build it with {.fn vehicle}."
    ))
  }
  if (missing(sink) || !inherits(sink, "skin_sink")) {
    cli::cli_abort(c(
      "{.arg sink} must be a {.cls skin_sink} object.",
      "i" = "Build it with {.fn perfect_sink} or {.fn finite_sink}."
    ))
  }
  if (missing(layers) || !is.list(layers) || length(layers) == 0L) {
    cli::cli_abort(c(
      "{.arg layers} must be a non-empty list of {.cls skin_layer} objects.",
      "i" = "Build each one with {.fn layer}."
    ))
  }
  for (i in seq_along(layers)) {
    if (!inherits(layers[[i]], "skin_layer")) {
      cli::cli_abort(c(
        "{.arg layers[[{i}]]} must be a {.cls skin_layer} object.",
        "i" = "Build it with {.fn layer}."
      ))
    }
  }
  scaling <- match.arg(scaling)

  area_cm2_val   <- .ensure_units_range(area, "cm^2", "area",
                                        min = 0, exclusive_min = TRUE)
  duration_min_  <- .ensure_units_int(duration, "min", "duration", min = 1L)
  mass_log_min   <- .ensure_units_int(mass_log_interval, "min",
                                      "mass_log_interval", min = 1L)
  cdp_log_min    <- .ensure_units_int(cdp_log_interval, "min",
                                      "cdp_log_interval", min = 1L)
  resolution_int <- .ensure_int(resolution, "resolution", min = 1L)
  max_module_val <- .ensure_dimensionless(max_module, "max_module",
                                          min = 0, exclusive_min = TRUE)

  # Internal nested-list shape consumed by the C++ binding. Field names on
  # the C++ side stay short (c_init, D, height, Vd, ...) so the binding
  # layer doesn't need to change.
  params <- list(
    sys = list(
      resolution      = resolution_int,
      max_module      = max_module_val,
      simulation_time = duration_min_
    ),
    log = list(
      scaling           = scaling,
      mass_log_interval = mass_log_min,
      cdp_log_interval  = cdp_log_min
    ),
    sink     = .sink_to_internal(sink),
    vehicle  = .vehicle_to_internal(vehicle, area_cm2_val),
    layers   = lapply(layers, .layer_to_internal),
    .meta    = list(area_cm2 = area_cm2_val,
                    sink_is_perfect = identical(sink$type, "perfect"))
  )

  res <- .cpp_validate(params)
  if (!isTRUE(res$ok)) {
    cli::cli_abort(c(
      "Invalid parameters.",
      "x" = "{res$error}"
    ))
  }

  class(params) <- "skin_params"
  params
}

# ---------- print methods ----------------------------------------------------

#' @export
print.skin_vehicle <- function(x, ...) {
  cat(sprintf("<skindiff vehicle> %s\n", x$name))
  cat(sprintf("  c_init        : %s\n", format(mg_per_ml(x$c_init_mg_per_ml))))
  cat(sprintf("  height        : %s\n", format(um(x$height_um))))
  cat(sprintf("  D             : %s\n", format(um2_per_min(x$D_um2_per_min))))
  cat(sprintf("  finite_dose   : %s\n", x$finite_dose))
  if (x$replace_after_min > 0) {
    cat(sprintf("  replace_after : %s\n", format(minutes(x$replace_after_min))))
  }
  if (x$remove_at_min > 0) {
    cat(sprintf("  remove_at     : %s\n", format(minutes(x$remove_at_min))))
  }
  invisible(x)
}

#' @export
print.skin_layer <- function(x, ...) {
  cat(sprintf("<skindiff layer> %s\n", x$name))
  cat(sprintf("  height        : %s\n", format(um(x$height_um))))
  cat(sprintf("  D             : %s\n", format(um2_per_min(x$D_um2_per_min))))
  cat(sprintf("  K             : %g\n", x$K))
  cat(sprintf("  cross_section : %g\n", x$cross_section))
  if (x$c_init_mg_per_ml != 0) {
    cat(sprintf("  c_init        : %s\n", format(mg_per_ml(x$c_init_mg_per_ml))))
  }
  invisible(x)
}

#' @export
print.skin_sink <- function(x, ...) {
  cat(sprintf("<skindiff sink: %s> %s\n", x$type, x$name))
  if (identical(x$type, "finite")) {
    cat(sprintf("  Vd            : %s\n", format(ml(x$Vd_ml))))
    if (x$c_init_mg_per_ml != 0) {
      cat(sprintf("  c_init        : %s\n", format(mg_per_ml(x$c_init_mg_per_ml))))
    }
  }
  invisible(x)
}

#' @export
print.skin_params <- function(x, ...) {
  cat("<skin_params>\n")
  cat(sprintf("  area               : %s\n",     format(cm2(x$.meta$area_cm2))))
  cat(sprintf("  duration           : %s\n",     format(minutes(x$sys$simulation_time))))
  cat(sprintf("  resolution         : %d cells/um\n", x$sys$resolution))
  cat(sprintf("  scaling            : %s\n",     x$log$scaling))
  cat("\n")
  cat(sprintf("  vehicle            : %s (h=%s, c0=%s, D=%s)\n",
              x$vehicle$name,
              format(um(x$vehicle$height)),
              format(mg_per_ml(x$vehicle$c_init)),
              format(um2_per_min(x$vehicle$D))))
  cat(sprintf("  layers (%d):\n", length(x$layers)))
  for (i in seq_along(x$layers)) {
    l <- x$layers[[i]]
    cat(sprintf("    [%d] %s (h=%s, D=%s, K=%g, cs=%g)\n",
                i, l$name,
                format(um(l$height)),
                format(um2_per_min(l$D)),
                l$K, l$cross_section))
  }
  if (isTRUE(x$.meta$sink_is_perfect)) {
    cat(sprintf("  sink               : %s (perfect)\n", x$sink$name))
  } else {
    cat(sprintf("  sink               : %s (finite, Vd=%s)\n",
                x$sink$name, format(ml(x$sink$Vd))))
  }
  invisible(x)
}

# ---------- internal: constructor objects -> serialised list shape ----------

# The C++ binding layer reads short field names (c_init, D, height, Vd, ...).
# These helpers translate the constructor objects' long names to the canonical
# short names. All values are bare numerics in canonical units at this point.

.vehicle_to_internal <- function(v, area_cm2) {
  list(
    name          = v$name,
    c_init        = v$c_init_mg_per_ml,
    app_area      = area_cm2,
    D             = v$D_um2_per_min,
    height        = v$height_um,
    replace_after = v$replace_after_min,
    remove_at     = v$remove_at_min,
    finite_dose   = v$finite_dose,
    log_mass      = v$log_mass,
    log_cdp       = v$log_cdp
  )
}

.layer_to_internal <- function(l) {
  list(
    name          = l$name,
    c_init        = l$c_init_mg_per_ml,
    D             = l$D_um2_per_min,
    K             = l$K,
    cross_section = l$cross_section,
    height        = l$height_um,
    log_mass      = l$log_mass,
    log_cdp       = l$log_cdp
  )
}

.sink_to_internal <- function(s) {
  list(
    name     = s$name,
    c_init   = s$c_init_mg_per_ml,
    Vd       = s$Vd_ml,
    log_mass = s$log_mass
  )
}
