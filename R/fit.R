#' Fit D / K parameters of a skin model to observed data
#'
#' Optimises selected layer-level diffusion coefficients (`D`) and
#' partition coefficients (`K`) so the simulator's predictions match
#' permeation and/or penetration observations. Other compartment
#' parameters (heights, cross-sections, application area, durations,
#' vehicle/sink properties) stay at the values supplied in `template`.
#'
#' For multi-subject fits, supply `template` as a named list of
#' `skin_params` (one per subject) -- typically subjects differ in skin
#' geometry. Fitted parameters (`fit_pars`) are *shared* across subjects;
#' subject-level random effects are not yet supported.
#'
#' @param template Either a single [skin_params()] object (single-subject
#'   fit) or a named list of `skin_params` keyed by subject id
#'   (multi-subject fit). Layer name sets must match across templates.
#' @param observations A list with one or both entries:
#'   * `permeation`  : a [permeation_obs()] object
#'   * `penetration` : a [penetration_obs()] object
#'
#'   Subjects in each obs object's `subject` column must appear as keys
#'   of `template` (or, for single-template fits, all observations must
#'   share the default subject).
#' @param fit_pars Named list mapping layer names to character vectors
#'   of the parameters to fit (`"D"`, `"K"`, or both). Layers not
#'   listed -- and parameters not listed within a layer -- stay at
#'   their template values.
#' @param bounds Optional named list with the same outer structure as
#'   `fit_pars`. Each entry is a list with `D` and/or `K` mapping to a
#'   length-2 vector. `D` bounds must carry units (e.g.
#'   `c(um2_per_min(0.01), um2_per_min(100))`); `K` bounds are
#'   dimensionless. Defaults: `D` in `[1e-3, 1e6] um^2/min`, `K` in
#'   `[1e-3, 1e6]`.
#' @param weights Per-block weighting strategy:
#'   * `"auto"` (default): use `1/sd^2` per point if `sd` is supplied,
#'     otherwise scale each block by `1/range(observed)^2` so blocks
#'     contribute comparably.
#'   * `"uniform"`: ignore SDs and block scaling.
#'   * `list(permeation = w_p, penetration = w_n)`: explicit per-block
#'     scalar multipliers.
#' @param transform Per-modality residual transform:
#'   `list(permeation = "linear", penetration = "log")` by default.
#'   `"log"` uses `log(predicted) - log(observed)`; `"linear"` uses
#'   `predicted - observed`.
#' @param optimizer Optimisation method.
#' @param n_starts Number of random-start optimisations (best is
#'   reported). Default 1.
#' @param n_boot Number of bootstrap resamples for confidence
#'   intervals. Default 0 (no bootstrap; only Hessian-based asymptotic
#'   SEs).
#' @param control Optional list passed to [stats::optim()] as `control`.
#'
#' @return A `skin_fit` object.
#' @export
skin_fit <- function(template,
                     observations,
                     fit_pars,
                     bounds    = list(),
                     weights   = "auto",
                     transform = list(permeation = "linear",
                                      penetration = "log"),
                     optimizer = c("L-BFGS-B", "nlminb"),
                     n_starts  = 1L,
                     n_boot    = 0L,
                     control   = list(maxit = 200, trace = 0)) {

  the_call <- match.call()
  optimizer <- match.arg(optimizer)

  # ---- Validate and normalise template ----
  template <- .normalise_template(template)
  subjects <- names(template)

  # ---- Validate observations ----
  obs <- .normalise_observations(observations, subjects)

  # ---- Validate fit_pars / bounds / sanity warnings ----
  spec <- .validate_fit_spec(template, fit_pars, bounds, obs)

  # ---- Build internal parameter index and bounds in log-space ----
  par_idx <- spec$par_idx
  log_lo  <- log(spec$bounds_lo)
  log_hi  <- log(spec$bounds_hi)
  theta0  <- log(spec$start)

  # ---- Compose loss closure ----
  loss_fn <- .make_loss(template, obs, par_idx, weights, transform)

  # ---- Optimise (optionally with multi-start) ----
  best <- .run_optim(theta0, log_lo, log_hi, loss_fn,
                     optimizer = optimizer,
                     n_starts  = n_starts,
                     control   = control)

  # ---- Hessian-based asymptotic SE (delta method, log -> linear) ----
  ses <- .compute_se(best$hessian, best$par, n_obs = obs$n_total,
                     n_par = length(theta0), residual_var = best$value)

  # ---- Build the skin_fit object ----
  estimate <- exp(best$par)
  par_table <- spec$par_meta
  par_table$estimate <- estimate
  par_table$se       <- ses$linear
  par_table$lower    <- spec$bounds_lo
  par_table$upper    <- spec$bounds_hi

  # ---- Final predictions and residuals at the optimum ----
  pred <- .predict_all(template, obs, par_idx, best$par)

  # ---- Optional bootstrap CIs ----
  boot <- if (n_boot > 0L) {
    .bootstrap_ci(template, obs, par_idx, weights, transform,
                  log_lo, log_hi, best$par, n_boot,
                  optimizer = optimizer, control = control)
  } else NULL

  out <- list(
    call         = the_call,
    template     = template,
    observations = obs,
    fit_pars     = fit_pars,
    par_table    = par_table,
    theta_hat    = best$par,
    loss         = best$value,
    iterations   = best$counts,
    convergence  = best$convergence,
    message      = best$message,
    hessian      = best$hessian,
    predictions  = pred,
    bootstrap    = boot
  )
  class(out) <- c("skin_fit", "list")
  out
}


# ============================================================================
#  Internal: template / observations normalisation
# ============================================================================

.normalise_template <- function(template) {
  if (inherits(template, "skin_params")) {
    list(default = template)
  } else if (is.list(template) && length(template) > 0L &&
             all(vapply(template, inherits, logical(1L), "skin_params"))) {
    if (is.null(names(template)) || any(!nzchar(names(template)))) {
      cli::cli_abort(c(
        "{.arg template} must be a named list of {.cls skin_params}.",
        "i" = "Use {.code list(s1 = template_s1, s2 = template_s2, ...)}."
      ))
    }
    # Verify all templates have the same layer name set
    layer_sets <- lapply(template, function(t) vapply(t$layers, function(l) l$name, character(1L)))
    if (length(unique(lapply(layer_sets, sort))) > 1L) {
      cli::cli_abort(c(
        "Multi-subject templates must share the same layer name set.",
        "i" = "Found different sets across subjects: {.val {unique(layer_sets)}}"
      ))
    }
    template
  } else {
    cli::cli_abort(c(
      "{.arg template} must be a {.cls skin_params} or a named list of them.",
      "x" = "Got {.obj_type_friendly {template}}."
    ))
  }
}

.normalise_observations <- function(observations, subjects) {
  if (!is.list(observations)) {
    cli::cli_abort("{.arg observations} must be a list.")
  }
  known <- c("permeation", "penetration")
  unknown <- setdiff(names(observations), known)
  if (length(unknown) > 0L) {
    cli::cli_abort(c(
      "{.arg observations} has unknown entries: {.val {unknown}}",
      "i" = "Allowed entries: {.val {known}}."
    ))
  }
  if (length(observations) == 0L) {
    cli::cli_abort("{.arg observations} must contain at least one of {.val permeation} or {.val penetration}.")
  }
  perm <- observations$permeation
  pen  <- observations$penetration
  if (!is.null(perm) && !inherits(perm, "permeation_obs")) {
    cli::cli_abort("{.field observations$permeation} must be a {.cls permeation_obs}.")
  }
  if (!is.null(pen) && !inherits(pen, "penetration_obs")) {
    cli::cli_abort("{.field observations$penetration} must be a {.cls penetration_obs}.")
  }

  obs_subjects <- unique(c(if (!is.null(perm)) perm$subject,
                            if (!is.null(pen)) pen$subject))
  if (length(subjects) == 1L && identical(subjects, "default")) {
    if (any(obs_subjects != "default")) {
      cli::cli_abort(c(
        "Multi-subject observations require a multi-subject {.arg template}.",
        "x" = "Observation subjects: {.val {obs_subjects}}",
        "i" = "Pass {.arg template} as a named list keyed by subject id."
      ))
    }
  } else {
    missing_subj <- setdiff(obs_subjects, subjects)
    if (length(missing_subj) > 0L) {
      cli::cli_abort(c(
        "Some observation subjects have no template.",
        "x" = "Missing template for subject(s): {.val {missing_subj}}",
        "i" = "Templates provided for: {.val {subjects}}"
      ))
    }
  }

  n_total <- (if (!is.null(perm)) perm$n else 0L) +
             (if (!is.null(pen)) pen$n else 0L)

  list(permeation  = perm,
       penetration = pen,
       subjects    = obs_subjects,
       n_total     = n_total)
}


# ============================================================================
#  Internal: fit_spec validation, default bounds, sanity warnings
# ============================================================================

.validate_fit_spec <- function(template, fit_pars, bounds, obs) {
  if (!is.list(fit_pars) || length(fit_pars) == 0L) {
    cli::cli_abort(c(
      "{.arg fit_pars} must be a non-empty named list.",
      "i" = "Example: {.code list(\"Stratum corneum\" = c(\"D\", \"K\"))}"
    ))
  }
  if (is.null(names(fit_pars)) || any(!nzchar(names(fit_pars)))) {
    cli::cli_abort("{.arg fit_pars} must be a named list (layer name -> parameters).")
  }

  # All templates share the same layer name set; pick the first.
  tpl <- template[[1L]]
  layer_names <- vapply(tpl$layers, function(l) l$name, character(1L))

  # Validate each entry of fit_pars
  par_meta_rows <- list()
  start_vals    <- numeric()
  bounds_lo     <- numeric()
  bounds_hi     <- numeric()
  par_names     <- character()
  par_idx_list  <- list()  # per layer: list of par-name -> index

  default_D_lo <- 1e-3                                     # um^2/min
  default_D_hi <- 1e6
  default_K_lo <- 1e-3
  default_K_hi <- 1e6

  idx <- 0L
  for (lname in names(fit_pars)) {
    if (!lname %in% layer_names) {
      cli::cli_abort(c(
        "Layer in {.arg fit_pars} is not present in the template.",
        "x" = "Unknown layer: {.val {lname}}",
        "i" = "Available layers: {.val {layer_names}}"
      ))
    }
    pars <- fit_pars[[lname]]
    if (!is.character(pars) || length(pars) == 0L) {
      cli::cli_abort("{.field fit_pars[[\"{lname}\"]]} must be a non-empty character vector.")
    }
    bad_pars <- setdiff(pars, c("D", "K"))
    if (length(bad_pars) > 0L) {
      cli::cli_abort(c(
        "Only {.val D} and {.val K} are fittable in v1.",
        "x" = "Got: {.val {bad_pars}} for layer {.val {lname}}"
      ))
    }
    # Get the layer's current values
    layer <- tpl$layers[[which(layer_names == lname)]]
    layer_idx <- which(layer_names == lname)

    layer_par_idx <- list()
    for (par in pars) {
      idx <- idx + 1L
      start <- if (par == "D") layer$D else layer$K
      lo_default <- if (par == "D") default_D_lo else default_K_lo
      hi_default <- if (par == "D") default_D_hi else default_K_hi

      # User-supplied bounds?
      lo <- lo_default
      hi <- hi_default
      if (!is.null(bounds[[lname]]) && !is.null(bounds[[lname]][[par]])) {
        b <- bounds[[lname]][[par]]
        if (length(b) != 2L) {
          cli::cli_abort("Bounds for {.val {lname}}/{.val {par}} must be length 2.")
        }
        # D bounds must carry units (length); K bounds are bare numeric
        if (par == "D") {
          b_canon <- .ensure_obs_units_col(b, "um^2/min",
                                           sprintf("bounds[[\"%s\"]]$D", lname))
          lo <- b_canon[1L]; hi <- b_canon[2L]
        } else {
          if (inherits(b, "units")) {
            cli::cli_abort("Bounds for {.val K} are dimensionless; do not pass a units object.")
          }
          if (!is.numeric(b)) cli::cli_abort("K bounds must be numeric.")
          lo <- b[1L]; hi <- b[2L]
        }
      }
      if (!(lo > 0 && hi > lo)) {
        cli::cli_abort(c(
          "Bounds for {.val {lname}}/{.val {par}} must be positive and ordered.",
          "x" = "Got [{.val {lo}}, {.val {hi}}]"
        ))
      }

      # Warn: bounds exclude the start value
      if (start < lo || start > hi) {
        cli::cli_warn(c(
          "Bounds for {.val {lname}}/{.val {par}} exclude the template value.",
          "x" = "Template = {.val {start}}; bounds = [{.val {lo}}, {.val {hi}}]",
          "i" = "Optimisation will start at the clamped boundary."
        ))
        start <- max(lo, min(hi, start))
      }

      layer_par_idx[[par]] <- list(idx = idx, layer_idx = layer_idx)

      par_meta_rows[[idx]] <- data.frame(
        layer = lname, par = par, stringsAsFactors = FALSE
      )
      start_vals[idx] <- start
      bounds_lo[idx]  <- lo
      bounds_hi[idx]  <- hi
      par_names[idx]  <- sprintf("%s[%s]", par, lname)
    }
    par_idx_list[[lname]] <- layer_par_idx
  }

  # Warn: bounds supplied for things not in fit_pars
  for (lname in names(bounds)) {
    if (!lname %in% names(fit_pars)) {
      cli::cli_warn("Bounds supplied for layer {.val {lname}} not in {.arg fit_pars}; ignored.")
      next
    }
    extra <- setdiff(names(bounds[[lname]]), fit_pars[[lname]])
    if (length(extra) > 0L) {
      cli::cli_warn("Bounds supplied for {.val {lname}}/{.val {extra}} not in {.arg fit_pars}; ignored.")
    }
  }

  # Warn: under-determined fit
  n_par   <- idx
  n_data  <- obs$n_total
  if (n_par >= n_data) {
    cli::cli_warn(c(
      "Fit may be under-determined.",
      "x" = "Number of fitted parameters ({.val {n_par}}) >= number of observations ({.val {n_data}}).",
      "i" = "Standard errors will likely be unreliable."
    ))
  }

  # Hard error: required logging missing.
  if (!is.null(obs$permeation)) {
    if (!isTRUE(tpl$sink$log_mass)) {
      cli::cli_abort(c(
        "Permeation data was supplied but the sink mass is not logged.",
        "i" = "Set {.code log_mass = TRUE} on the sink in the template."
      ))
    }
  }
  if (!is.null(obs$penetration)) {
    # Every layer that any strip overlaps must have log_cdp = TRUE.
    layer_tops <- c(0, cumsum(vapply(tpl$layers, function(l) l$height,
                                     numeric(1L))))
    n_layers <- length(tpl$layers)
    layer_top_um <- layer_tops[seq_len(n_layers)]
    layer_bot_um <- layer_tops[seq(2L, n_layers + 1L)]
    needed <- character()
    for (i in seq_along(obs$penetration$time_min)) {
      st <- obs$penetration$depth_top_um[i]
      sb <- obs$penetration$depth_bottom_um[i]
      for (j in seq_len(n_layers)) {
        if (max(layer_top_um[j], st) < min(layer_bot_um[j], sb)) {
          needed <- union(needed, tpl$layers[[j]]$name)
        }
      }
    }
    not_logged <- needed[!vapply(needed, function(nm) {
      tpl$layers[[which(vapply(tpl$layers, function(l) l$name == nm,
                               logical(1L)))]]$log_cdp
    }, logical(1L))]
    if (length(not_logged) > 0L) {
      cli::cli_abort(c(
        "Penetration data overlaps layer(s) that are not CDP-logged.",
        "x" = "Missing {.code log_cdp = TRUE} for: {.val {not_logged}}",
        "i" = "The fit cannot predict the strip concentration without the layer's depth profile."
      ))
    }
  }

  par_meta <- do.call(rbind, par_meta_rows)
  par_meta$par_name <- par_names

  list(
    par_idx   = par_idx_list,
    par_meta  = par_meta,
    start     = start_vals,
    bounds_lo = bounds_lo,
    bounds_hi = bounds_hi
  )
}


# ============================================================================
#  Internal: apply theta to a template (per subject) and simulate
# ============================================================================

.apply_theta <- function(tpl, par_idx, theta_log) {
  layer_names <- vapply(tpl$layers, function(l) l$name, character(1L))
  for (lname in names(par_idx)) {
    li <- which(layer_names == lname)
    for (par in names(par_idx[[lname]])) {
      val <- exp(theta_log[par_idx[[lname]][[par]]$idx])
      tpl$layers[[li]][[par]] <- val
    }
  }
  # Force scaling to "ng" so prediction units are always canonical.
  tpl$log$scaling <- "ng"
  tpl
}

# Simulate one subject and return canonical-unit predictions.
.simulate_subject <- function(tpl) {
  raw <- .cpp_simulate(unclass(tpl), show_progress = FALSE)
  raw  # raw cpp result; keep structure flat for fast access
}

# ============================================================================
#  Internal: predict permeation / penetration at canonical units (ng/cm^2,
#  ng/ml). The simulator runs with scaling = "ng" so values are already in
#  ng on the mass side.
# ============================================================================

.predict_permeation_subject <- function(raw, sink_name, area_cm2, times_min) {
  sink_series <- raw$mass[[sink_name]]
  q_grid <- sink_series$value / area_cm2     # ng / cm^2
  t_grid <- sink_series$time
  stats::approx(t_grid, q_grid, xout = times_min, rule = 2)$y
}

# Compute predicted strip concentrations for a set of penetration rows.
# `cdp_raw` is the cpp-returned cdp list: each entry has $time, $depth_um,
# $conc (matrix, [depth, time]).
.predict_penetration_subject <- function(raw, layer_meta,
                                         times_min, depth_top_um, depth_bottom_um) {
  cdp <- raw$cdp
  n <- length(times_min)
  out <- numeric(n)

  # layer_meta: data.frame with columns name, top_um, bottom_um, dx_um (skin frame)
  # for each skin layer (vehicle excluded).

  for (k in seq_len(n)) {
    tk    <- times_min[k]
    s_top <- depth_top_um[k]
    s_bot <- depth_bottom_um[k]
    strip_thickness <- max(s_bot - s_top, .Machine$double.eps)

    int_mass <- 0  # ng/ml * um (per unit area-relative)
    for (i in seq_len(nrow(layer_meta))) {
      lm <- layer_meta[i, ]
      lt <- lm$top_um; lb <- lm$bottom_um
      ovl_top <- max(lt, s_top)
      ovl_bot <- min(lb, s_bot)
      if (ovl_top >= ovl_bot) next  # no overlap
      if (!lm$name %in% names(cdp)) next  # not logged
      s <- cdp[[lm$name]]
      s_t <- s$time
      mid_local <- s$depth_um            # cell midpoints, layer-local
      conc_mat  <- s$conc                # [depth, time]
      # Interpolate each row across time
      conc_at_tk <- vapply(seq_along(mid_local), function(d) {
        stats::approx(s_t, conc_mat[d, ], xout = tk, rule = 2)$y
      }, numeric(1L))
      # Map to skin-frame depth: skin_depth = lt + mid_local
      skin_mid <- lt + mid_local
      half_dx  <- lm$dx_um / 2
      cell_top <- skin_mid - half_dx
      cell_bot <- skin_mid + half_dx
      cell_overlap <- pmax(0, pmin(cell_bot, ovl_bot) -
                                pmax(cell_top, ovl_top))
      int_mass <- int_mass + sum(conc_at_tk * cell_overlap)
    }
    out[k] <- int_mass / strip_thickness   # ng/ml
  }
  out
}

# Build per-subject layer metadata used by the penetration predictor.
.layer_meta <- function(tpl, raw_cdp) {
  cum <- 0
  rows <- list()
  for (l in tpl$layers) {
    h <- l$height
    # Cell width: derive from cdp midpoints if logged, else fall back to layer height.
    dx <- if (l$name %in% names(raw_cdp) && length(raw_cdp[[l$name]]$depth_um) > 1L) {
      raw_cdp[[l$name]]$depth_um[2L] - raw_cdp[[l$name]]$depth_um[1L]
    } else {
      h
    }
    rows[[length(rows) + 1L]] <- data.frame(
      name = l$name, top_um = cum, bottom_um = cum + h, dx_um = dx,
      stringsAsFactors = FALSE
    )
    cum <- cum + h
  }
  do.call(rbind, rows)
}


# ============================================================================
#  Internal: loss function closure
# ============================================================================

.make_loss <- function(template, obs, par_idx, weights, transform) {

  perm_obs <- obs$permeation
  pen_obs  <- obs$penetration

  # Pre-compute per-subject indexers
  perm_subj_idx <- if (!is.null(perm_obs)) split(seq_along(perm_obs$subject), perm_obs$subject)
                   else NULL
  pen_subj_idx  <- if (!is.null(pen_obs))  split(seq_along(pen_obs$subject),  pen_obs$subject)
                   else NULL

  perm_transform <- if (is.null(transform$permeation)) "linear" else transform$permeation
  pen_transform  <- if (is.null(transform$penetration)) "log"   else transform$penetration

  # Pre-compute block-level weights. For "auto" without SDs, use range-based
  # scaling on linear residuals; log residuals are already unitless and need
  # no block scaling.
  block_weight_for <- function(values, transform_kind, modality) {
    if (is.list(weights) && !is.null(weights[[modality]])) return(weights[[modality]])
    if (identical(weights, "uniform")) return(1)
    if (identical(transform_kind, "log"))  return(1)
    1 / max(diff(range(values)), .Machine$double.eps)^2
  }
  perm_w_block <- if (!is.null(perm_obs)) {
    block_weight_for(perm_obs$q_per_area_ng_cm2, perm_transform, "permeation")
  } else NA
  pen_w_block <- if (!is.null(pen_obs)) {
    block_weight_for(pen_obs$conc_ng_ml, pen_transform, "penetration")
  } else NA

  perm_use_var <- !is.null(perm_obs) && identical(weights, "auto") &&
                  !all(is.na(perm_obs$sd_ng_cm2))
  pen_use_var  <- !is.null(pen_obs)  && identical(weights, "auto") &&
                  !all(is.na(pen_obs$sd_ng_ml))

  function(theta_log) {
    total <- 0
    # For each subject, simulate once and use the result for both modalities
    for (subj in names(template)) {
      tpl_s <- .apply_theta(template[[subj]], par_idx, theta_log)
      raw   <- .simulate_subject(tpl_s)
      area_cm2 <- tpl_s$.meta$area_cm2

      # Permeation contribution
      if (!is.null(perm_obs) && subj %in% names(perm_subj_idx)) {
        rows <- perm_subj_idx[[subj]]
        pred <- .predict_permeation_subject(raw, tpl_s$sink$name,
                                            area_cm2, perm_obs$time_min[rows])
        obs_v <- perm_obs$q_per_area_ng_cm2[rows]
        sd_v  <- perm_obs$sd_ng_cm2[rows]
        if (perm_transform == "log") {
          eps <- 1e-30
          r <- log(pmax(pred, eps)) - log(pmax(obs_v, eps))
        } else {
          r <- pred - obs_v
        }
        if (perm_use_var) {
          # Use 1/sd^2 weighting per-point
          w <- 1 / pmax(sd_v, .Machine$double.eps)^2
          total <- total + sum(w * r^2)
        } else {
          total <- total + perm_w_block * sum(r^2)
        }
      }

      # Penetration contribution
      if (!is.null(pen_obs) && subj %in% names(pen_subj_idx)) {
        rows <- pen_subj_idx[[subj]]
        lm <- .layer_meta(tpl_s, raw$cdp)
        pred <- .predict_penetration_subject(
          raw, lm, pen_obs$time_min[rows],
          pen_obs$depth_top_um[rows], pen_obs$depth_bottom_um[rows]
        )
        obs_v <- pen_obs$conc_ng_ml[rows]
        sd_v  <- pen_obs$sd_ng_ml[rows]
        if (pen_transform == "log") {
          eps <- 1e-30
          r <- log(pmax(pred, eps)) - log(pmax(obs_v, eps))
        } else {
          r <- pred - obs_v
        }
        if (pen_use_var) {
          w <- 1 / pmax(sd_v, .Machine$double.eps)^2
          total <- total + sum(w * r^2)
        } else {
          total <- total + pen_w_block * sum(r^2)
        }
      }
    }
    total
  }
}


# ============================================================================
#  Internal: optimiser wrapper
# ============================================================================

.run_optim <- function(theta0, lower, upper, loss_fn,
                       optimizer = "L-BFGS-B",
                       n_starts  = 1L,
                       control   = list()) {

  one_start <- function(start) {
    res <- tryCatch(
      stats::optim(par = start, fn = loss_fn, method = optimizer,
                   lower = lower, upper = upper,
                   control = control, hessian = TRUE),
      error = function(e) NULL
    )
    res
  }

  starts <- list(theta0)
  if (n_starts > 1L) {
    extras <- replicate(n_starts - 1L,
                        stats::runif(length(theta0), min = lower, max = upper),
                        simplify = FALSE)
    starts <- c(starts, extras)
  }

  results <- lapply(starts, one_start)
  results <- results[!vapply(results, is.null, logical(1L))]
  if (length(results) == 0L) {
    cli::cli_abort("Optimiser failed for all starts.")
  }
  values <- vapply(results, function(r) r$value, numeric(1L))
  best   <- results[[which.min(values)]]
  best
}


# ============================================================================
#  Internal: Hessian -> SE
# ============================================================================

.compute_se <- function(hessian, par_log, n_obs, n_par, residual_var) {
  cov_log <- tryCatch(solve(0.5 * hessian), error = function(e) NULL)
  if (is.null(cov_log)) {
    se_log    <- rep(NA_real_, length(par_log))
    se_linear <- rep(NA_real_, length(par_log))
  } else {
    se_log    <- sqrt(pmax(diag(cov_log), 0))
    # Delta method: SE(exp(theta)) ~ exp(theta) * SE(theta)
    se_linear <- exp(par_log) * se_log
  }
  list(log = se_log, linear = se_linear)
}


# ============================================================================
#  Internal: predictions at the optimum, for residual / fitted reporting
# ============================================================================

.predict_all <- function(template, obs, par_idx, theta_log) {
  perm_pred <- NULL
  pen_pred  <- NULL
  for (subj in names(template)) {
    tpl_s <- .apply_theta(template[[subj]], par_idx, theta_log)
    raw   <- .simulate_subject(tpl_s)
    area_cm2 <- tpl_s$.meta$area_cm2

    if (!is.null(obs$permeation)) {
      rows <- which(obs$permeation$subject == subj)
      if (length(rows) > 0L) {
        pr <- .predict_permeation_subject(raw, tpl_s$sink$name,
                                          area_cm2, obs$permeation$time_min[rows])
        perm_pred <- rbind(perm_pred, data.frame(
          subject  = subj,
          time_min = obs$permeation$time_min[rows],
          observed = obs$permeation$q_per_area_ng_cm2[rows],
          predicted = pr,
          residual  = pr - obs$permeation$q_per_area_ng_cm2[rows],
          stringsAsFactors = FALSE
        ))
      }
    }
    if (!is.null(obs$penetration)) {
      rows <- which(obs$penetration$subject == subj)
      if (length(rows) > 0L) {
        lm <- .layer_meta(tpl_s, raw$cdp)
        pr <- .predict_penetration_subject(
          raw, lm, obs$penetration$time_min[rows],
          obs$penetration$depth_top_um[rows], obs$penetration$depth_bottom_um[rows]
        )
        pen_pred <- rbind(pen_pred, data.frame(
          subject     = subj,
          time_min    = obs$penetration$time_min[rows],
          depth_top   = obs$penetration$depth_top_um[rows],
          depth_bot   = obs$penetration$depth_bottom_um[rows],
          observed    = obs$penetration$conc_ng_ml[rows],
          predicted   = pr,
          residual    = pr - obs$penetration$conc_ng_ml[rows],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  list(permeation = perm_pred, penetration = pen_pred)
}


# ============================================================================
#  Internal: simple bootstrap CIs (block-resampling per modality)
# ============================================================================

.bootstrap_ci <- function(template, obs, par_idx, weights, transform,
                          log_lo, log_hi, theta_hat, n_boot,
                          optimizer, control) {
  results <- matrix(NA_real_, nrow = n_boot, ncol = length(theta_hat))
  for (b in seq_len(n_boot)) {
    obs_b <- obs
    if (!is.null(obs$permeation)) {
      idx <- sample.int(obs$permeation$n, replace = TRUE)
      obs_b$permeation <- .subset_obs(obs$permeation, idx, "permeation_obs")
    }
    if (!is.null(obs$penetration)) {
      idx <- sample.int(obs$penetration$n, replace = TRUE)
      obs_b$penetration <- .subset_obs(obs$penetration, idx, "penetration_obs")
    }
    loss_b <- .make_loss(template, obs_b, par_idx, weights, transform)
    fit_b <- tryCatch(
      .run_optim(theta_hat, log_lo, log_hi, loss_b,
                 optimizer = optimizer, n_starts = 1L,
                 control = control),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) results[b, ] <- exp(fit_b$par)
  }
  ci <- t(apply(results, 2, stats::quantile,
                probs = c(0.025, 0.975), na.rm = TRUE))
  colnames(ci) <- c("ci_lo_2.5", "ci_hi_97.5")
  list(samples = results, ci = ci)
}

.subset_obs <- function(o, idx, cls) {
  out <- o
  for (nm in setdiff(names(o), c("n"))) {
    if (length(o[[nm]]) == o$n) out[[nm]] <- o[[nm]][idx]
  }
  out$n <- length(idx)
  class(out) <- c(cls, "list")
  out
}
