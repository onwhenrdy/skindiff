#' Build a validated set of skindiff simulation parameters
#'
#' Constructs the parameter object consumed by [skin_simulate()]. All units
#' are stated explicitly. The returned list is validated both at the R level
#' and again inside the C++ engine.
#'
#' @param vehicle Named list describing the vehicle / donor compartment.
#'   Fields:
#'   * `c_init` (numeric, mg/ml): initial donor concentration.
#'   * `app_area` (numeric, cm^2): application area.
#'   * `D` (numeric, um^2/min): diffusion coefficient inside the donor.
#'   * `height` (integer, um): donor thickness.
#'   * `finite_dose` (logical): if `FALSE`, the donor concentration stays
#'     clamped at `c_init` for the entire run.
#'   * `replace_after` (integer, min): if > 0, donor is reset to `c_init`
#'     every `replace_after` minutes.
#'   * `remove_at` (integer, min): if > 0, donor is removed at this time.
#'   * `name` (character): label used in result output.
#'   * `log_mass`, `log_cdp` (logical): whether to record this compartment.
#' @param layers `data.frame` with one row per skin layer (top to bottom).
#'   Required columns: `name`, `height` (um), `D` (um^2/min), `K`
#'   (partition coefficient relative to vehicle), `cross_section` (in
#'   `(0, 1]`). Optional: `c_init`, `log_mass`, `log_cdp`.
#' @param sink Named list with fields `Vd` (ml), `c_init` (mg/ml), `name`,
#'   `log_mass`. Only `Vd` is consulted when PK is disabled.
#' @param pk Named list with fields `enabled` (logical) and `thalf`
#'   (numeric, hours). When `enabled = FALSE` the sink behaves as a perfect
#'   sink.
#' @param sim_time Total simulation time, integer, minutes.
#' @param resolution Mesh refinement: number of cells per micrometre.
#' @param disc_method One of `"bk"` (refined Babucke-Kloker mesh) or
#'   `"equidist"`.
#' @param matrix_method Discretization scheme used by the engine. One of:
#'   `"Activity_FVM"` (default) -- a cell-centred finite-volume scheme in
#'   the activity variable `u = c/K`. Symmetric, second-order in the
#'   bulk, and handles partition jumps cleanly because `u` is continuous
#'   at interfaces by construction. Or `"DSkin_1_4"` -- the legacy Crank-
#'   style finite-difference scheme in the c (concentration) variable,
#'   with `K_l/K_c` corrections in the stencil. The legacy scheme has a
#'   measurable activity-continuity artefact at K-jumps (~15% on a
#'   realistic SC/DSL stack vs ~0.5% for `"Activity_FVM"`). See
#'   `tests/testthat/test-analytical.R` and `test-realistic.R` for the
#'   side-by-side comparison.
#' @param eta BK transition scaling factor, in `(0, 1]`. Ignored for
#'   `disc_method = "equidist"`.
#' @param max_module Stability target for the implicit sub-step count.
#' @param scaling Output mass scaling: `"mg"`, `"ug"`, or `"ng"`.
#' @param mass_log_interval Sample interval for mass time-series, minutes.
#' @param cdp_log_interval Sample interval for concentration-depth profiles,
#'   minutes.
#'
#' @return A list with class `"skin_params"` ready for [skin_simulate()].
#'
#' @export
skin_params <- function(
  vehicle,
  layers,
  sink,
  pk = list(enabled = FALSE, thalf = 1),
  sim_time = 600L,
  resolution = 1L,
  disc_method = c("bk", "equidist"),
  matrix_method = c("Activity_FVM", "DSkin_1_4"),
  eta = 0.6,
  max_module = 50,
  scaling = c("mg", "ug", "ng"),
  mass_log_interval = 1L,
  cdp_log_interval = 1L
) {
  disc_method <- match.arg(disc_method)
  matrix_method <- match.arg(matrix_method)
  scaling <- match.arg(scaling)

  vehicle_l <- .normalize_vehicle(vehicle)
  layers_l <- .normalize_layers(layers)
  sink_l <- .normalize_sink(sink)
  pk_l <- .normalize_pk(pk)

  params <- list(
    sys = list(
      disc_method     = disc_method,
      matrix_method   = matrix_method,
      resolution      = .as_int(resolution, "resolution", min_val = 1L),
      max_module      = .as_dbl(max_module, "max_module", min_val = 0, exclusive = TRUE),
      eta             = .as_dbl(eta, "eta", min_val = 0, max_val = 1, exclusive_min = TRUE),
      simulation_time = .as_int(sim_time, "sim_time", min_val = 1L)
    ),
    log = list(
      scaling           = scaling,
      mass_log_interval = .as_int(mass_log_interval, "mass_log_interval", min_val = 1L),
      cdp_log_interval  = .as_int(cdp_log_interval,  "cdp_log_interval",  min_val = 1L)
    ),
    pk = pk_l,
    sink = sink_l,
    vehicle = vehicle_l,
    layers = layers_l
  )

  # Defer to C++ for the canonical validation pass.
  res <- .cpp_validate(params)
  if (!isTRUE(res$ok)) {
    stop(sprintf("Invalid parameters: %s", res$error), call. = FALSE)
  }

  class(params) <- "skin_params"
  params
}

#' @export
print.skin_params <- function(x, ...) {
  cat("<skin_params>\n")
  cat(sprintf("  vehicle       : %s, %s ml, h=%d um, D=%g\n",
              x$vehicle$name,
              format(x$vehicle$app_area * x$vehicle$height / 1e4),
              x$vehicle$height, x$vehicle$D))
  cat(sprintf("  layers        : %d\n", length(x$layers)))
  for (i in seq_along(x$layers)) {
    l <- x$layers[[i]]
    cat(sprintf("    [%d] %s: h=%d um, D=%g, K=%g, cs=%g\n",
                i, l$name, l$height, l$D, l$K, l$cross_section))
  }
  cat(sprintf("  sink          : %s, Vd=%g ml%s\n",
              x$sink$name, x$sink$Vd,
              if (isTRUE(x$pk$enabled)) sprintf(" (PK, t1/2=%g h)", x$pk$thalf)
              else " (perfect)"))
  cat(sprintf("  sim_time      : %d min\n", x$sys$simulation_time))
  cat(sprintf("  disc          : %s, resolution=%d, max_module=%g\n",
              x$sys$disc_method, x$sys$resolution, x$sys$max_module))
  cat(sprintf("  scaling       : %s\n", x$log$scaling))
  invisible(x)
}

# ---------- internal normalization helpers ----------

.normalize_vehicle <- function(v) {
  if (!is.list(v)) stop("`vehicle` must be a named list", call. = FALSE)
  list(
    name          = .as_chr(v$name %||% "Vehicle", "vehicle$name"),
    c_init        = .as_dbl(v$c_init   %||% 1.0, "vehicle$c_init",   min_val = 0),
    app_area      = .as_dbl(v$app_area %||% 1.0, "vehicle$app_area", min_val = 0, exclusive = TRUE),
    D             = .as_dbl(v$D        %||% 1.0, "vehicle$D",        min_val = 0),
    height        = .as_int(v$height   %||% 10L, "vehicle$height",   min_val = 3L),
    replace_after = .as_int(v$replace_after %||% 0L, "vehicle$replace_after", min_val = 0L),
    remove_at     = .as_int(v$remove_at     %||% 0L, "vehicle$remove_at",     min_val = 0L),
    finite_dose   = .as_lgl(v$finite_dose %||% TRUE,  "vehicle$finite_dose"),
    log_mass      = .as_lgl(v$log_mass    %||% TRUE,  "vehicle$log_mass"),
    log_cdp       = .as_lgl(v$log_cdp     %||% FALSE, "vehicle$log_cdp")
  )
}

.normalize_layers <- function(layers) {
  if (is.data.frame(layers)) {
    required <- c("name", "height", "D", "K", "cross_section")
    missing <- setdiff(required, names(layers))
    if (length(missing) > 0) {
      stop(sprintf("`layers` data.frame missing columns: %s",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
    n <- nrow(layers)
    if (n == 0L) stop("`layers` must have at least one row", call. = FALSE)
    out <- vector("list", n)
    for (i in seq_len(n)) {
      out[[i]] <- list(
        name          = .as_chr(layers$name[[i]],   sprintf("layers$name[%d]", i)),
        c_init        = .as_dbl(if (is.null(layers$c_init)) 0.0 else layers$c_init[[i]],
                                sprintf("layers$c_init[%d]", i), min_val = 0),
        D             = .as_dbl(layers$D[[i]],     sprintf("layers$D[%d]", i),     min_val = 0),
        K             = .as_dbl(layers$K[[i]],     sprintf("layers$K[%d]", i),     min_val = 0, exclusive = TRUE),
        cross_section = .as_dbl(layers$cross_section[[i]],
                                sprintf("layers$cross_section[%d]", i),
                                min_val = 0, max_val = 1, exclusive_min = TRUE),
        height        = .as_int(layers$height[[i]], sprintf("layers$height[%d]", i), min_val = 3L),
        log_mass      = .as_lgl(if (is.null(layers$log_mass)) TRUE  else layers$log_mass[[i]],
                                sprintf("layers$log_mass[%d]", i)),
        log_cdp       = .as_lgl(if (is.null(layers$log_cdp))  FALSE else layers$log_cdp[[i]],
                                sprintf("layers$log_cdp[%d]", i))
      )
    }
    return(out)
  }
  if (is.list(layers) && length(layers) > 0) {
    return(lapply(seq_along(layers), function(i) {
      l <- layers[[i]]
      if (!is.list(l)) stop(sprintf("layers[[%d]] must be a list", i), call. = FALSE)
      list(
        name          = .as_chr(l$name,  sprintf("layers[[%d]]$name", i)),
        c_init        = .as_dbl(l$c_init %||% 0.0, sprintf("layers[[%d]]$c_init", i), min_val = 0),
        D             = .as_dbl(l$D,     sprintf("layers[[%d]]$D", i),     min_val = 0),
        K             = .as_dbl(l$K,     sprintf("layers[[%d]]$K", i),     min_val = 0, exclusive = TRUE),
        cross_section = .as_dbl(l$cross_section, sprintf("layers[[%d]]$cross_section", i),
                                min_val = 0, max_val = 1, exclusive_min = TRUE),
        height        = .as_int(l$height, sprintf("layers[[%d]]$height", i), min_val = 3L),
        log_mass      = .as_lgl(l$log_mass %||% TRUE,  sprintf("layers[[%d]]$log_mass", i)),
        log_cdp       = .as_lgl(l$log_cdp  %||% FALSE, sprintf("layers[[%d]]$log_cdp", i))
      )
    }))
  }
  stop("`layers` must be a data.frame or non-empty list of lists", call. = FALSE)
}

.normalize_sink <- function(s) {
  if (!is.list(s)) stop("`sink` must be a named list", call. = FALSE)
  list(
    name     = .as_chr(s$name %||% "Sink",    "sink$name"),
    c_init   = .as_dbl(s$c_init %||% 0.0,     "sink$c_init", min_val = 0),
    Vd       = .as_dbl(s$Vd     %||% 1.0,     "sink$Vd",     min_val = 0, exclusive = TRUE),
    log_mass = .as_lgl(s$log_mass %||% TRUE,  "sink$log_mass")
  )
}

.normalize_pk <- function(p) {
  if (!is.list(p)) stop("`pk` must be a named list", call. = FALSE)
  list(
    enabled = .as_lgl(p$enabled %||% FALSE, "pk$enabled"),
    thalf   = .as_dbl(p$thalf   %||% 1.0,   "pk$thalf", min_val = 0)
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.as_int <- function(x, name, min_val = NULL, max_val = NULL) {
  if (length(x) != 1 || !is.numeric(x) || is.na(x)) {
    stop(sprintf("%s must be a single integer value", name), call. = FALSE)
  }
  if (abs(x - round(x)) > 1e-9) {
    stop(sprintf("%s must be an integer (got %g)", name, x), call. = FALSE)
  }
  xi <- as.integer(round(x))
  if (!is.null(min_val) && xi < min_val) {
    stop(sprintf("%s must be >= %d (got %d)", name, as.integer(min_val), xi), call. = FALSE)
  }
  if (!is.null(max_val) && xi > max_val) {
    stop(sprintf("%s must be <= %d (got %d)", name, as.integer(max_val), xi), call. = FALSE)
  }
  xi
}

.as_dbl <- function(x, name, min_val = NULL, max_val = NULL,
                    exclusive = FALSE, exclusive_min = FALSE, exclusive_max = FALSE) {
  if (length(x) != 1 || !is.numeric(x) || is.na(x)) {
    stop(sprintf("%s must be a single numeric value", name), call. = FALSE)
  }
  xd <- as.numeric(x)
  excl_lo <- exclusive || exclusive_min
  excl_hi <- exclusive || exclusive_max
  if (!is.null(min_val)) {
    bad <- if (excl_lo) xd <= min_val else xd < min_val
    if (bad) stop(sprintf("%s must be %s %g (got %g)",
                          name, if (excl_lo) ">" else ">=", min_val, xd),
                  call. = FALSE)
  }
  if (!is.null(max_val)) {
    bad <- if (excl_hi) xd >= max_val else xd > max_val
    if (bad) stop(sprintf("%s must be %s %g (got %g)",
                          name, if (excl_hi) "<" else "<=", max_val, xd),
                  call. = FALSE)
  }
  xd
}

.as_lgl <- function(x, name) {
  if (length(x) != 1 || !is.logical(x) || is.na(x)) {
    stop(sprintf("%s must be a single TRUE/FALSE", name), call. = FALSE)
  }
  as.logical(x)
}

.as_chr <- function(x, name) {
  if (length(x) != 1 || !is.character(x) || is.na(x) || !nzchar(x)) {
    stop(sprintf("%s must be a single non-empty string", name), call. = FALSE)
  }
  as.character(x)
}
