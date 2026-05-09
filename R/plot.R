#' Plot a skindiff simulation result
#'
#' ggplot2 visualisations of a `skin_result`. The package keeps `ggplot2` as
#' a `Suggests` dependency, so this method requires `ggplot2` to be
#' installed but does not pull it in by default.
#'
#' @param object A `skin_result` from [skin_simulate()].
#' @param what Which view to render. One of:
#'   \describe{
#'     \item{`"mass"`}{Mass time-series, one line per logged compartment.}
#'     \item{`"concentration"`}{Concentration time-series, one line per
#'       logged compartment. Compartments whose concentration is `NA`
#'       throughout (e.g. a `perfect_sink()`) are dropped.}
#'     \item{`"permeated"`}{Cumulative permeated mass per unit area at the
#'       receptor, `Q(t)`.}
#'     \item{`"flux"`}{Instantaneous flux at the receptor, `dQ/dt`.}
#'     \item{`"profile"`}{Concentration-depth profile *line plot*: a
#'       single x-axis stacks all selected compartments along physical
#'       depth (skin surface = 0; vehicle, when included, sits at
#'       negative depth), with one coloured line per requested time
#'       slice. K-jumps at compartment boundaries appear as
#'       discontinuities, with dashed vertical lines marking the
#'       interfaces. Linear y by default; override with
#'       `+ scale_y_log10()`.}
#'   }
#' @param times A units-of-time vector of explicit time slices to render
#'   for `what = "profile"`. If `NULL`, `n_times` equally-spaced slices
#'   from 0 to the simulation end are used.
#' @param n_times Number of equally-spaced time slices for `what =
#'   "profile"` when `times` is `NULL`. Default 6.
#' @param compartments Character vector of compartment names to include
#'   in `what = "profile"`. If `NULL` (default), all compartments with
#'   `log_cdp = TRUE` are used (subject to `include_vehicle`). Listing
#'   compartments explicitly bypasses `include_vehicle`.
#' @param include_vehicle If `TRUE` (default), the vehicle is included in
#'   `what = "profile"` when it has `log_cdp = TRUE`. Set `FALSE` to
#'   plot only the skin layers. Ignored if `compartments` is non-`NULL`.
#' @param ... Reserved for future arguments.
#'
#' @return A `ggplot` object. Compose further with `+` as usual.
#' @exportS3Method ggplot2::autoplot
autoplot.skin_result <- function(object,
                                 what = c("mass", "concentration",
                                          "permeated", "flux", "profile"),
                                 times           = NULL,
                                 n_times         = 6L,
                                 compartments    = NULL,
                                 include_vehicle = TRUE,
                                 ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg ggplot2} is required for plotting.",
      "i" = "Install with {.code install.packages('ggplot2')}."
    ))
  }
  what <- match.arg(what)
  switch(what,
    mass          = .plot_mass(object),
    concentration = .plot_concentration(object),
    permeated     = .plot_permeated(object),
    flux          = .plot_flux(object),
    profile       = .plot_profile_lines(object,
                                        times = times,
                                        n_times = n_times,
                                        compartments = compartments,
                                        include_vehicle = include_vehicle)
  )
}

# ---------- per-view implementations ----------------------------------------

.plot_mass <- function(res) {
  df <- res$mass
  long_df <- .pivot_longer_units(df, "time")
  unit_x <- .unit_label(df$time)
  unit_y <- .unit_label(df[[setdiff(names(df), "time")[1]]])

  ggplot2::ggplot(long_df,
                  ggplot2::aes(x = time, y = value, color = compartment)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(
      x     = sprintf("time [%s]", unit_x),
      y     = sprintf("mass [%s]", unit_y),
      color = "compartment"
    )
}

.plot_concentration <- function(res) {
  df <- res$concentration
  # Drop columns that are all NA (e.g. perfect_sink concentration).
  value_cols <- setdiff(names(df), "time")
  has_data <- vapply(value_cols,
                     function(nm) any(!is.na(as.numeric(df[[nm]]))),
                     logical(1))
  df <- df[, c("time", value_cols[has_data]), drop = FALSE]

  long_df <- .pivot_longer_units(df, "time")
  unit_x <- .unit_label(df$time)
  unit_y <- .unit_label(df[[setdiff(names(df), "time")[1]]])

  ggplot2::ggplot(long_df,
                  ggplot2::aes(x = time, y = value, color = compartment)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(
      x     = sprintf("time [%s]", unit_x),
      y     = sprintf("concentration [%s]", unit_y),
      color = "compartment"
    )
}

.plot_permeated <- function(res) {
  perm <- permeated(res)
  ggplot2::ggplot(
    data.frame(time = as.numeric(perm$time),
               Q    = as.numeric(perm$Q)),
    ggplot2::aes(x = time, y = Q)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(
      x = sprintf("time [%s]", .unit_label(perm$time)),
      y = sprintf("Q [%s]",    .unit_label(perm$Q))
    )
}

.plot_flux <- function(res) {
  fl <- flux(res)
  ggplot2::ggplot(
    data.frame(time = as.numeric(fl$time),
               flux = as.numeric(fl$flux)),
    ggplot2::aes(x = time, y = flux)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(
      x = sprintf("time [%s]", .unit_label(fl$time)),
      y = sprintf("flux [%s]", .unit_label(fl$flux))
    )
}

.plot_profile_lines <- function(res, times, n_times,
                                compartments, include_vehicle) {
  if (length(res$cdp) == 0L) {
    cli::cli_abort(c(
      "No concentration-depth profiles were logged.",
      "i" = "Set {.code log_cdp = TRUE} on the compartments you want to inspect."
    ))
  }

  available    <- names(res$cdp)
  vehicle_name <- res$params$vehicle$name

  # Resolve compartment selection.
  if (is.null(compartments)) {
    selected <- available
    if (!isTRUE(include_vehicle)) {
      selected <- setdiff(selected, vehicle_name)
    }
  } else {
    not_avail <- setdiff(compartments, available)
    if (length(not_avail) > 0L) {
      cli::cli_abort(c(
        "Some compartments are not available for the profile plot.",
        "x" = "Not logged or unknown: {.val {not_avail}}",
        "i" = "Compartments with {.code log_cdp = TRUE}: {.val {available}}"
      ))
    }
    # Preserve physical order from `available`, not user-supplied order.
    selected <- intersect(available, compartments)
  }
  if (length(selected) == 0L) {
    cli::cli_abort("No compartments selected for the profile plot.")
  }

  # Compartment heights from params, indexed by name.
  heights_um <- list()
  heights_um[[vehicle_name]] <- res$params$vehicle$height
  for (l in res$params$layers) heights_um[[l$name]] <- l$height

  # Cumulative offsets (top of stack = 0).
  offsets_top <- list()
  cum_top <- 0
  for (nm in available) {
    offsets_top[[nm]] <- cum_top
    cum_top <- cum_top + heights_um[[nm]]
  }
  # Shift so the skin surface (top of first non-vehicle compartment) sits at 0.
  if (vehicle_name %in% available) {
    h_v <- heights_um[[vehicle_name]]
    for (nm in names(offsets_top)) {
      offsets_top[[nm]] <- offsets_top[[nm]] - h_v
    }
  }

  # Time slices: explicit (units required) or n equidistant from 0..end.
  t_grid_min <- as.numeric(res$cdp[[available[1L]]]$time)
  end_min    <- max(t_grid_min)
  if (is.null(times)) {
    n_times <- as.integer(n_times)
    if (length(n_times) != 1L || is.na(n_times) || n_times < 2L) {
      cli::cli_abort("{.arg n_times} must be a single integer >= 2.")
    }
    times_min <- seq(0, end_min, length.out = n_times)
  } else {
    if (!inherits(times, "units")) {
      cli::cli_abort(c(
        "{.arg times} must be a units-of-time vector.",
        "i" = "Use a {.pkg skindiff} unit helper, e.g. {.code hours(c(1, 2, 4))}."
      ))
    }
    times_min <- as.numeric(units::set_units(times, "min", mode = "standard"))
    if (any(times_min < 0) || any(times_min > end_min + 1e-9)) {
      cli::cli_abort(c(
        "{.arg times} contains values outside the simulation range.",
        "x" = "Got values up to {.val {max(times_min)}} min; sim ends at {.val {end_min}} min."
      ))
    }
  }

  # Build long data: one row per (compartment, depth, time slice).
  parts <- list()
  for (nm in selected) {
    s        <- res$cdp[[nm]]
    s_t      <- as.numeric(s$time)
    s_d      <- as.numeric(s$depth)        # midpoints, relative to compartment top
    conc_mat <- unclass(s$conc)             # [depth, time]
    offset   <- offsets_top[[nm]]
    for (tt in times_min) {
      conc_at_t <- apply(conc_mat, 1L, function(row) {
        stats::approx(s_t, row, xout = tt, rule = 2)$y
      })
      parts[[length(parts) + 1L]] <- data.frame(
        depth_global = offset + s_d,
        conc         = conc_at_t,
        time_min     = tt,
        compartment  = nm,
        stringsAsFactors = FALSE
      )
    }
  }
  long_df <- do.call(rbind, parts)
  long_df$compartment <- factor(long_df$compartment, levels = selected)

  # Time labels: hours if duration spans an hour or more, minutes otherwise.
  use_hours <- end_min >= 60
  unit_lbl  <- if (use_hours) "h" else "min"
  vals      <- if (use_hours) times_min / 60 else times_min
  time_levels <- as.character(times_min)
  long_df$time_label <- factor(as.character(long_df$time_min),
                               levels = time_levels,
                               labels = sprintf("%g %s", vals, unit_lbl))

  # Compartment boundary positions for dashed vlines (interior boundaries
  # only -- skip the outermost edges of the plot).
  edges <- c(
    vapply(selected, function(nm) offsets_top[[nm]], numeric(1L)),
    offsets_top[[utils::tail(selected, 1L)]] +
      heights_um[[utils::tail(selected, 1L)]]
  )
  edges <- sort(unique(edges))
  interior_edges <- edges[-c(1L, length(edges))]

  unit_depth <- .unit_label(res$cdp[[1L]]$depth)
  unit_conc  <- .unit_label(res$cdp[[1L]]$conc)

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = depth_global, y = conc,
                                    color = time_label,
                                    group = interaction(compartment, time_label)))

  if (length(interior_edges) > 0L) {
    p <- p + ggplot2::geom_vline(xintercept = interior_edges,
                                 linetype = "dashed", color = "grey60")
  }

  p +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_color_viridis_d(option = "magma", end = 0.85) +
    ggplot2::labs(
      x     = sprintf("depth [%s] (skin surface = 0)", unit_depth),
      y     = sprintf("concentration [%s]",            unit_conc),
      color = "time"
    )
}

# ---------- internal helpers ------------------------------------------------

# Wide -> long with units stripped; the wide df has an id column
# (e.g. "time") plus one column per compartment.
.pivot_longer_units <- function(df, id_col) {
  value_cols <- setdiff(names(df), id_col)
  parts <- vector("list", length(value_cols))
  id_bare <- as.numeric(df[[id_col]])
  for (i in seq_along(value_cols)) {
    nm <- value_cols[i]
    parts[[i]] <- data.frame(
      time        = id_bare,
      compartment = nm,
      value       = as.numeric(df[[nm]]),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, parts)
  out$compartment <- factor(out$compartment, levels = value_cols)
  out
}

.unit_label <- function(x) {
  if (inherits(x, "units")) units::deparse_unit(x) else "?"
}

# Suppress R CMD check NOTEs about ggplot2 NSE column references in aes().
utils::globalVariables(c("time", "value", "compartment", "Q", "flux",
                         "conc", "depth_global", "time_label",
                         "observed", "predicted", "subject", "depth_mid",
                         "modality", "time_min", "x"))


#' Plot a `skin_fit`: data + model overlay
#'
#' For each fitted modality (`permeation` and/or `penetration`), draws
#' the observations as points (with SD bars if available) and the
#' best-fit model prediction. Faceted by subject.
#'
#' @param object A `skin_fit` from [skin_fit()].
#' @param what Optional character: one of `"permeation"`,
#'   `"penetration"`, or `"observed_vs_predicted"`. If `NULL` (default),
#'   all available modalities are shown stacked.
#' @param ... Reserved.
#'
#' @return A `ggplot` object.
#' @exportS3Method ggplot2::autoplot
autoplot.skin_fit <- function(object,
                              what = NULL,
                              ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg ggplot2} is required for plotting.",
      "i" = "Install with {.code install.packages('ggplot2')}."
    ))
  }
  available <- character()
  if (!is.null(object$predictions$permeation))  available <- c(available, "permeation")
  if (!is.null(object$predictions$penetration)) available <- c(available, "penetration")

  what <- if (is.null(what)) "all"
          else match.arg(what, c("permeation", "penetration",
                                 "observed_vs_predicted", "all"))

  if (identical(what, "observed_vs_predicted")) return(.plot_fit_obs_vs_pred(object))

  views <- if (identical(what, "all")) available else what
  views <- intersect(views, available)
  if (length(views) == 0L) {
    cli::cli_abort("No predictions available to plot.")
  }
  # Build per-modality plots and combine vertically via patchwork-like row.
  if (length(views) == 1L) {
    return(switch(views,
      permeation  = .plot_fit_permeation(object),
      penetration = .plot_fit_penetration(object)
    ))
  }
  # Two modalities: combine via faceting on a 'modality' column.
  parts <- lapply(views, function(v) {
    df <- object$predictions[[v]]
    df$modality <- v
    df
  })
  # Render via two stacked plots returned as a list-of-ggplots? We commit
  # to a single ggplot output: facet by modality with free scales.
  perm_df <- object$predictions$permeation
  pen_df  <- object$predictions$penetration

  # Build a long observed/predicted pair per facet by harmonising columns.
  perm_long <- if (!is.null(perm_df)) {
    base <- data.frame(
      modality = "permeation",
      subject  = perm_df$subject,
      x        = perm_df$time_min,
      observed = perm_df$observed,
      predicted = perm_df$predicted,
      stringsAsFactors = FALSE
    )
    base
  } else NULL
  pen_long <- if (!is.null(pen_df)) {
    data.frame(
      modality = "penetration",
      subject  = pen_df$subject,
      x        = (pen_df$depth_top + pen_df$depth_bot) / 2,   # depth midpoint, um
      observed = pen_df$observed,
      predicted = pen_df$predicted,
      stringsAsFactors = FALSE
    )
  } else NULL
  long <- rbind(perm_long, pen_long)

  ggplot2::ggplot(long, ggplot2::aes(x = x)) +
    ggplot2::geom_point(ggplot2::aes(y = observed, color = subject)) +
    ggplot2::geom_line(ggplot2::aes(y = predicted, color = subject,
                                    group = subject), linewidth = 0.7) +
    ggplot2::facet_wrap(~ modality, scales = "free", ncol = 1L) +
    ggplot2::labs(x = "time [min] / depth midpoint [um]",
                  y = "value (ng/cm^2 for permeation, ng/ml for penetration)",
                  color = "subject")
}

.plot_fit_permeation <- function(fit) {
  p <- fit$predictions$permeation
  ggplot2::ggplot(p, ggplot2::aes(x = time_min)) +
    ggplot2::geom_point(ggplot2::aes(y = observed, color = subject)) +
    ggplot2::geom_line(ggplot2::aes(y = predicted, color = subject,
                                    group = subject), linewidth = 0.7) +
    ggplot2::labs(x = "time [min]",
                  y = "Q [ng/cm^2]",
                  color = "subject",
                  title = "permeation: data + best-fit model")
}

.plot_fit_penetration <- function(fit) {
  p <- fit$predictions$penetration
  p$depth_mid <- (p$depth_top + p$depth_bot) / 2
  ggplot2::ggplot(p, ggplot2::aes(x = depth_mid)) +
    ggplot2::geom_point(ggplot2::aes(y = observed, color = subject)) +
    ggplot2::geom_line(ggplot2::aes(y = predicted, color = subject,
                                    group = interaction(subject, time_min)),
                       linewidth = 0.7) +
    ggplot2::facet_wrap(~ time_min, scales = "free_y",
                        labeller = ggplot2::label_both) +
    ggplot2::labs(x = "depth midpoint [um]",
                  y = "concentration [ng/ml]",
                  color = "subject",
                  title = "penetration: data + best-fit model")
}

.plot_fit_obs_vs_pred <- function(fit) {
  parts <- list()
  if (!is.null(fit$predictions$permeation)) {
    df <- fit$predictions$permeation
    df$modality <- "permeation"
    parts[[length(parts) + 1L]] <- df[, c("subject", "modality", "observed", "predicted")]
  }
  if (!is.null(fit$predictions$penetration)) {
    df <- fit$predictions$penetration
    df$modality <- "penetration"
    parts[[length(parts) + 1L]] <- df[, c("subject", "modality", "observed", "predicted")]
  }
  long <- do.call(rbind, parts)
  ggplot2::ggplot(long, ggplot2::aes(x = observed, y = predicted,
                                     color = subject, shape = modality)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          color = "grey50") +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "observed", y = "predicted",
                  color = "subject", shape = "modality",
                  title = "fit quality: predicted vs observed")
}
