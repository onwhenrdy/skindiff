#' @export
print.skin_fit <- function(x, ...) {
  cat("<skin_fit>\n")
  cat(sprintf("  loss          : %.6g\n", x$loss))
  cat(sprintf("  convergence   : %d (%s)\n",
              x$convergence,
              if (isTRUE(x$convergence == 0)) "converged" else "see $message"))
  cat(sprintf("  fitted pars   : %d\n", nrow(x$par_table)))
  for (i in seq_len(nrow(x$par_table))) {
    row <- x$par_table[i, ]
    se_str <- if (is.na(row$se)) "(SE unavailable)"
              else sprintf("+/- %.4g", row$se)
    cat(sprintf("    %-20s = %.6g %s\n", row$par_name, row$estimate, se_str))
  }
  if (!is.null(x$bootstrap)) {
    ci <- x$bootstrap$ci
    cat(sprintf("  bootstrap CIs (%d resamples):\n", nrow(x$bootstrap$samples)))
    for (i in seq_len(nrow(ci))) {
      cat(sprintf("    %-20s [%.4g, %.4g]\n",
                  x$par_table$par_name[i],
                  ci[i, 1L], ci[i, 2L]))
    }
  }
  invisible(x)
}

#' @export
summary.skin_fit <- function(object, ...) {
  print(object)
  cat("\n  observations  : ", object$observations$n_total, " row(s)\n", sep = "")
  if (!is.null(object$predictions$permeation)) {
    p <- object$predictions$permeation
    cat(sprintf("  permeation    : RMSE = %.4g ng/cm^2 over %d points\n",
                sqrt(mean(p$residual^2)), nrow(p)))
  }
  if (!is.null(object$predictions$penetration)) {
    p <- object$predictions$penetration
    cat(sprintf("  penetration   : RMSE = %.4g ng/ml over %d points\n",
                sqrt(mean(p$residual^2)), nrow(p)))
  }
  invisible(object)
}

#' @export
coef.skin_fit <- function(object, ...) {
  out <- as.list(object$par_table$estimate)
  names(out) <- object$par_table$par_name
  out
}

#' @export
residuals.skin_fit <- function(object, ...) {
  parts <- list()
  if (!is.null(object$predictions$permeation)) {
    p <- object$predictions$permeation
    parts$permeation  <- data.frame(
      modality = "permeation", subject = p$subject,
      time_min = p$time_min,
      observed = p$observed, predicted = p$predicted,
      residual = p$residual,
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(object$predictions$penetration)) {
    p <- object$predictions$penetration
    parts$penetration <- data.frame(
      modality = "penetration", subject = p$subject,
      time_min = p$time_min,
      depth_top = p$depth_top, depth_bot = p$depth_bot,
      observed = p$observed, predicted = p$predicted,
      residual = p$residual,
      stringsAsFactors = FALSE
    )
  }
  parts
}

#' @export
fitted.skin_fit <- function(object, ...) {
  res <- stats::residuals(object)
  for (m in names(res)) res[[m]]$residual <- NULL
  res
}

#' Extract the best-fit `skin_params` from a `skin_fit`
#'
#' Reconstructs a `skin_params` (or list of them, for multi-subject fits)
#' with the fitted `D` and `K` values substituted into the original
#' template(s). The resulting object is ready for [skin_simulate()] or
#' further analysis with [metrics()] / [autoplot()].
#'
#' @param fit A `skin_fit` object.
#' @param subject For multi-subject fits, the subject id whose template
#'   to extract. If `NULL` (default), returns a named list of all
#'   subjects' templates.
#'
#' @return A `skin_params` object or a named list of them.
#' @export
skin_params_from_fit <- function(fit, subject = NULL) {
  if (!inherits(fit, "skin_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls skin_fit} object.")
  }
  par_idx <- list()
  layer_names_global <- vapply(fit$template[[1L]]$layers,
                               function(l) l$name, character(1L))
  # Reconstruct par_idx structure from par_table
  pt <- fit$par_table
  for (i in seq_len(nrow(pt))) {
    lname <- pt$layer[i]; par <- pt$par[i]
    li <- which(layer_names_global == lname)
    par_idx[[lname]] <- par_idx[[lname]] %||% list()
    par_idx[[lname]][[par]] <- list(idx = i, layer_idx = li)
  }
  out_list <- lapply(names(fit$template), function(s) {
    tpl <- fit$template[[s]]
    tpl <- .apply_theta(tpl, par_idx, fit$theta_hat)
    # Restore original scaling from the user's template
    tpl$log$scaling <- fit$template[[s]]$log$scaling
    class(tpl) <- "skin_params"
    tpl
  })
  names(out_list) <- names(fit$template)

  if (length(out_list) == 1L && identical(names(out_list), "default")) {
    out_list[[1L]]
  } else if (!is.null(subject)) {
    if (!subject %in% names(out_list)) {
      cli::cli_abort(c(
        "Unknown subject {.val {subject}}.",
        "i" = "Available: {.val {names(out_list)}}"
      ))
    }
    out_list[[subject]]
  } else {
    out_list
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a
