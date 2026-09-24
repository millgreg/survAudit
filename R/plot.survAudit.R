# Diagnostic plotting for survAudit

#' Plot Cox Proportional Hazards Model Diagnostics
#'
#' Produces ggplot2-based diagnostic displays for a \code{survAudit} object.
#' Five diagnostic panels are available:
#' \itemize{
#'   \item \strong{Proportional Hazards (\code{which = "ph"}):} Plots scaled
#'   Schoenfeld residuals against transformed time for each covariate, with a
#'   LOESS smoother and horizontal reference lines at the estimated coefficients.
#'   \item \strong{Functional Form (\code{which = "functional"}):} Plots
#'   martingale residuals from a null (covariate-free) Cox model against
#'   each continuous covariate, with a LOESS smoother and a linear fit
#'   reference line to screen for non-linear effects.
#'   \item \strong{Influence Diagnostics (\code{which = "influence"}):} Plots
#'   standardized DFBETA statistics for each covariate across observation
#'   indices, with diagnostic threshold reference lines at
#'   \eqn{\pm 2 / \sqrt{n}}.
#'   \item \strong{Outlier Assessment (\code{which = "outliers"}):} A four-panel
#'   display showing martingale, deviance, log-odds, and normal deviate residuals
#'   against the linear predictor, with reference threshold lines.
#'   \item \strong{Global Goodness-of-Fit (\code{which = "gof"}):} Plots the cumulative
#'   hazard of the Cox-Snell residuals against the residuals themselves, with a
#'   dashed 45-degree reference line.
#' }
#'
#' \strong{Large-Sample Subsampling:}
#' For large datasets (e.g., \eqn{n > 5000}), \code{plot.survAudit} subsamples
#' data to \code{max_points} observations per panel (default 5000) to prevent
#' overplotting and rendering delays:
#' \itemize{
#'   \item \strong{\code{"ph"} and \code{"functional"}:} Randomly subsamples points
#'   per covariate; smoothers are computed on this subsample for fast display.
#'   \item \strong{\code{"influence"} and \code{"outliers"}:} Extreme-preserving
#'   subsampling. All observations exceeding diagnostic thresholds
#'   (\eqn{\pm 2/\sqrt{n}} for DFBETAs, \eqn{\pm 1.96} for deviance and normal deviate,
#'   \eqn{\pm 3.66} for log-odds, and \eqn{|r| > 2} for martingale) are strictly
#'   retained; only non-flagged points are subsampled.
#'   \item \strong{\code{"gof"}:} Systematic quantile subsampling across the
#'   Cox-Snell cumulative hazard curve.
#' }
#' A plot caption indicates whenever subsampling is active. To display all points
#' without subsampling, set \code{max_points = NULL}.
#'
#' @param x An object of class \code{survAudit}.
#' @param which Character vector specifying which plots to produce.
#'   Subset of \code{c("ph", "functional", "influence", "outliers", "gof")}.
#'   Defaults to all available diagnostics.
#' @param vars Optional character or numeric vector specifying covariates to
#'   include in covariate-specific plots (\code{"ph"}, \code{"functional"}, and
#'   \code{"influence"}). If \code{NULL} (the default), all available covariates
#'   are included.
#' @param max_points Integer specifying the maximum number of observations to
#'   display per panel (and to use for smoothing). Defaults to 5000. If \code{NULL}
#'   or non-positive, all observations are used without subsampling.
#' @param ask Logical. If \code{TRUE} (the default in interactive
#'   sessions), the user is prompted between plots.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot} object (if a single plot was requested) or a
#'   named list of \code{ggplot} objects (if multiple plots were produced),
#'   invisibly.
#'
#' @export
#'
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ trt + celltype + karno + age,
#'              data = veteran)
#' audit <- survAudit(fit, data = veteran)
#' plot(audit, which = "ph")
#' plot(audit, which = "functional", vars = "age")
#' plot(audit, which = "influence", max_points = 500)
plot.survAudit <- function(x,
                           which = c("ph", "functional",
                                     "influence", "outliers", "gof"),
                           vars = NULL,
                           max_points = 5000L,
                           ask = interactive(),
                           ...) {

  which <- match.arg(which, choices = c("ph", "functional",
                                        "influence", "outliers", "gof"),
                     several.ok = TRUE)

  # Colour palette ──────────────────────────────────────────────
  col_point  <- "#999999"
  col_smooth <- "#2166AC"
  col_ref    <- "#B2182B"
  col_thresh <- "#D6604D"
  alpha_pt   <- 0.4

  plots <- list()

  # ── PH ──────────────────────────────────────────────────────
  if ("ph" %in% which) {
    if (is.null(x$ph) || is.null(x$ph$zph)) {
      message("PH diagnostics not available; skipping 'ph' plot.")
    } else {
      p_ph <- .plot_ph(x$ph, col_point, col_smooth,
                       col_ref, alpha_pt, vars = vars, max_points = max_points)
      if (!is.null(p_ph)) plots[["ph"]] <- p_ph
    }
  }

  # ── Functional form ────────────────────────────────────────
  if ("functional" %in% which) {
    if (is.null(x$functional_form) ||
        length(x$functional_form$results) == 0L) {
      message("Functional form diagnostics not available; ",
              "skipping 'functional' plot.")
    } else {
      p_ff <- .plot_functional(
        x$functional_form, col_point, col_smooth, col_ref, alpha_pt,
        vars = vars, max_points = max_points
      )
      if (!is.null(p_ff)) plots[["functional"]] <- p_ff
    }
  }

  # ── Influence ──────────────────────────────────────────────
  if ("influence" %in% which) {
    if (is.null(x$influence)) {
      message("Influence diagnostics not available; ",
              "skipping 'influence' plot.")
    } else {
      p_inf <- .plot_influence(
        x$influence, col_point, col_smooth, col_ref, col_thresh, alpha_pt,
        vars = vars, max_points = max_points
      )
      if (!is.null(p_inf)) plots[["influence"]] <- p_inf
    }
  }

  # ── Outliers ───────────────────────────────────────────────
  if ("outliers" %in% which) {
    if (is.null(x$outliers)) {
      message("Outlier diagnostics not available; ",
              "skipping 'outliers' plot.")
    } else {
      plots[["outliers"]] <- .plot_outliers(
        x$outliers, col_point, col_ref, alpha_pt, max_points = max_points
      )
    }
  }

  # ── Global Goodness-of-Fit ───────────────────────────────────
  if ("gof" %in% which) {
    if (is.null(x$gof) || length(x$gof$plot_x) == 0L) {
      message("Goodness-of-fit diagnostics not available; ",
              "skipping 'gof' plot.")
    } else {
      plots[["gof"]] <- .plot_gof(
        x$gof, col_point, col_ref, alpha_pt, max_points = max_points
      )
    }
  }

  # ── Display ────────────────────────────────────────────────
  if (length(plots) == 0L) {
    message("No diagnostic plots available.")
    return(invisible(NULL))
  }

  if (isTRUE(ask) && length(plots) > 1L) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }

  for (i in seq_along(plots)) {
    print(plots[[i]])
  }

  if (length(plots) == 1L) {
    return(invisible(plots[[1L]]))
  }
  invisible(plots)
}


# ═══════════════════════════════════════════════════════════════════
# Internal plotting helpers
# ═══════════════════════════════════════════════════════════════════

#' Safe sampling helper avoiding the base R length-1 integer trap
#'
#' @param idx Integer vector of candidate indices.
#' @param size Integer number of samples to draw.
#' @return Sampled vector of indices.
#' @keywords internal
#' @noRd
.sample_indices <- function(idx, size) {
  if (length(idx) <= 1L) {
    idx[seq_len(min(length(idx), size))]
  } else {
    sample(idx, size)
  }
}


#' Plot PH diagnostics (Schoenfeld residuals)
#'
#' @param ph PH diagnostics list from the survAudit object.
#' @param col_point Point colour.
#' @param col_smooth Smooth line colour.
#' @param col_ref Reference line colour.
#' @param alpha_pt Point alpha.
#' @param vars Optional character or integer vector of covariates to plot.
#' @param max_points Optional maximum number of points to plot/smooth per covariate.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_ph <- function(ph, col_point, col_smooth, col_ref, alpha_pt, vars = NULL,
                     max_points = 5000L) {

  zph <- ph$zph

  # Build long-format data.frame using un-aggregated residuals for full transparency
  if (!is.null(ph$raw_residuals)) {
    y_mat <- ph$raw_residuals
    time_vals <- ph$raw_time
  } else {
    # Fallback for older survAudit objects that don't have raw_residuals saved
    y_mat <- zph$y
    time_vals <- zph$x
  }
  
  if (is.null(dim(y_mat))) {
    # Single covariate
    y_mat <- matrix(y_mat, ncol = 1L)
    colnames(y_mat) <- names(stats::coef(zph))
    if (is.null(colnames(y_mat))) colnames(y_mat) <- "covariate"
  }
  
  var_names <- colnames(y_mat)

  if (!is.null(vars)) {
    if (is.numeric(vars)) {
      vars <- var_names[vars[!is.na(vars) & vars >= 1 & vars <= length(var_names)]]
    }
    matched_vars <- intersect(var_names, vars)
    if (length(matched_vars) == 0L) {
      warning("None of the requested covariates in 'vars' were found for the PH plot.", call. = FALSE)
      return(NULL)
    }
    var_names <- matched_vars
    y_mat <- y_mat[, var_names, drop = FALSE]
  }

  n_t <- length(time_vals)
  n_vars <- ncol(y_mat)

  df_full <- data.frame(
    time     = rep(time_vals, times = n_vars),
    residual = as.vector(y_mat),
    variable = rep(var_names, each = n_t),
    stringsAsFactors = FALSE
  )

  # Subsample points for plotting and smoothing in large datasets to prevent overplotting
  is_subsampled <- FALSE
  df_plot <- df_full
  if (!is.null(max_points) && is.finite(max_points) && max_points > 0 && n_t > max_points) {
    set.seed(123)
    df_plot <- do.call(rbind, lapply(split(df_full, df_full$variable), function(sub_df) {
      if (nrow(sub_df) > max_points) {
        sub_df[.sample_indices(seq_len(nrow(sub_df)), max_points), ]
      } else {
        sub_df
      }
    }))
    is_subsampled <- TRUE
  }

  # Calculate reference horizontal lines (cox.zph$y contains scaled Schoenfeld + beta_hat, so colMeans recovers beta_hat)
  beta_hat <- colMeans(y_mat, na.rm = TRUE)
  df_ref <- data.frame(
    variable = var_names,
    beta = beta_hat,
    stringsAsFactors = FALSE
  )

  tname <- switch(ph$transform,
                  "km" = "Kaplan-Meier",
                  "log" = "Log",
                  "rank" = "Rank",
                  "identity" = "Identity",
                  ph$transform)

  plot_caption <- if (is_subsampled) {
    paste0("Points subsampled to ", format(max_points, big.mark = ","),
           " observations per covariate.")
  } else {
    NULL
  }

  ggplot(df_plot, aes(x = .data$time, y = .data$residual)) +
    geom_point(colour = col_point, alpha = alpha_pt, size = 1) +
    geom_smooth(method = "loess", formula = y ~ x,
                se = TRUE, colour = col_smooth, linewidth = 0.8,
                fill = col_smooth, alpha = 0.15) +
    geom_hline(data = df_ref, aes(yintercept = .data$beta), linetype = "dashed", colour = col_ref,
               linewidth = 0.5) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(
      title   = "Proportional Hazards Diagnostics: Scaled Schoenfeld Residuals",
      x       = paste0("Transformed Time (", tname, ")"),
      y       = "Scaled Schoenfeld Residual",
      caption = plot_caption
    ) +
    theme_minimal() +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 9, colour = "#666666", face = "italic"),
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 10),
      strip.text   = element_text(size = 12, face = "bold")
    )
}


#' Plot functional form diagnostics
#'
#' @param ff Functional form diagnostics list.
#' @param col_point Point colour.
#' @param col_smooth Smooth line colour.
#' @param col_ref Reference line colour.
#' @param alpha_pt Point alpha.
#' @param vars Optional character or integer vector of covariates to plot.
#' @param max_points Optional maximum number of points to plot/smooth per covariate.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_functional <- function(ff, col_point, col_smooth,
                             col_ref, alpha_pt, vars = NULL,
                             max_points = 5000L) {

  target_vars <- names(ff$results)
  if (!is.null(vars)) {
    if (is.numeric(vars)) {
      vars <- target_vars[vars[!is.na(vars) & vars >= 1 & vars <= length(target_vars)]]
    }
    matched_vars <- intersect(target_vars, vars)
    if (length(matched_vars) == 0L) {
      warning("None of the requested covariates in 'vars' were found for the functional form plot.", call. = FALSE)
      return(NULL)
    }
    target_vars <- matched_vars
  }

  dfs <- lapply(target_vars, function(vname) {
    res <- ff$results[[vname]]
    data.frame(
      covariate_value = res$covariate_values,
      residual        = res$residuals,
      variable        = vname,
      stringsAsFactors = FALSE
    )
  })

  df_full <- do.call(rbind, dfs)

  # Subsample points for plotting and smoothing in large datasets to prevent overplotting
  is_subsampled <- FALSE
  df_plot <- df_full
  if (nrow(df_full) > 0 && !is.null(max_points) && is.finite(max_points) && max_points > 0) {
    counts <- table(df_full$variable)
    if (any(counts > max_points)) {
      set.seed(123)
      df_plot <- do.call(rbind, lapply(split(df_full, df_full$variable), function(sub_df) {
        if (nrow(sub_df) > max_points) {
          sub_df[.sample_indices(seq_len(nrow(sub_df)), max_points), ]
        } else {
          sub_df
        }
      }))
      is_subsampled <- TRUE
    }
  }

  plot_caption <- if (is_subsampled) {
    paste0("Points subsampled to ", format(max_points, big.mark = ","),
           " observations per covariate.")
  } else {
    NULL
  }

  ggplot(df_plot, aes(x = .data$covariate_value, y = .data$residual)) +
    geom_point(colour = col_point, alpha = alpha_pt, size = 1) +
    geom_smooth(method = "loess", formula = y ~ x,
                se = TRUE, colour = col_smooth, linewidth = 0.8,
                fill = col_smooth, alpha = 0.15) +
    geom_smooth(method = "lm", formula = y ~ x,
                se = FALSE, colour = col_ref,
                linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~ variable, scales = "free") +
    labs(
      title   = "Functional Form Assessment: Martingale Residuals",
      x       = "Covariate Value",
      y       = "Martingale Residual",
      caption = plot_caption
    ) +
    theme_minimal() +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 9, colour = "#666666", face = "italic"),
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 10),
      strip.text   = element_text(size = 12, face = "bold")
    )
}


#' Plot influence diagnostics
#'
#' @param inf Influence diagnostics list.
#' @param col_point Point colour.
#' @param col_smooth Smooth line colour (unused here).
#' @param col_ref Reference line colour.
#' @param col_thresh Threshold line colour.
#' @param alpha_pt Point alpha.
#' @param vars Optional character or integer vector of covariates to plot.
#' @param max_points Optional maximum number of points to plot per covariate.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_influence <- function(inf, col_point, col_smooth,
                            col_ref, col_thresh, alpha_pt, vars = NULL,
                            max_points = 5000L) {

  # DFBETAs in long format
  dfb <- inf$dfbetas
  n <- nrow(dfb)
  var_names <- colnames(dfb)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(ncol(dfb)))

  if (!is.null(vars)) {
    if (is.numeric(vars)) {
      vars <- var_names[vars[!is.na(vars) & vars >= 1 & vars <= length(var_names)]]
    }
    matched_vars <- intersect(var_names, vars)
    if (length(matched_vars) == 0L) {
      warning("None of the requested covariates in 'vars' were found for the influence plot.", call. = FALSE)
      return(NULL)
    }
    var_names <- matched_vars
    dfb <- dfb[, var_names, drop = FALSE]
  }

  p <- ncol(dfb)

  df_full <- data.frame(
    obs      = rep(seq_len(n), times = p),
    value    = as.vector(dfb),
    variable = rep(var_names, each = n),
    stringsAsFactors = FALSE
  )

  df_plot <- df_full
  threshold <- inf$threshold

  # Subsample non-flagged points for plotting in large datasets while strictly retaining all flagged points
  is_subsampled <- FALSE
  if (!is.null(max_points) && is.finite(max_points) && max_points > 0 && n > max_points) {
    set.seed(123)
    df_plot <- do.call(rbind, lapply(split(df_full, df_full$variable), function(sub_df) {
      if (nrow(sub_df) <= max_points) return(sub_df)
      is_extreme <- abs(sub_df$value) > threshold
      n_extreme <- sum(is_extreme)
      if (n_extreme >= max_points) {
        # Strict preservation: retain all flagged cases; add zero normal cases
        sub_df[which(is_extreme), ]
      } else {
        n_sample <- max_points - n_extreme
        idx_normal <- which(!is_extreme)
        sampled_normal <- .sample_indices(idx_normal, n_sample)
        sub_df[sort(c(which(is_extreme), sampled_normal)), ]
      }
    }))
    is_subsampled <- TRUE
  }

  # Threshold data for DFBETAs panels only
  thresh_df <- data.frame(
    variable  = rep(var_names, each = 2L),
    yintercept = rep(c(threshold, -threshold), times = p),
    stringsAsFactors = FALSE
  )

  plot_caption <- if (is_subsampled) {
    paste0("Points subsampled to ", format(max_points, big.mark = ","),
           " observations per covariate (all flagged cases retained).")
  } else {
    NULL
  }

  ggplot(df_plot, aes(x = .data$obs, y = .data$value)) +
    geom_point(colour = col_point, alpha = alpha_pt, size = 1) +
    geom_hline(
      data        = thresh_df,
      aes(yintercept = .data$yintercept),
      linetype    = "dotted",
      colour      = col_thresh,
      linewidth   = 0.5
    ) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(
      title   = "Influence Diagnostics: DFBETAs",
      x       = "Observation Index",
      y       = "Standardized DFBETA",
      caption = plot_caption
    ) +
    theme_minimal() +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 9, colour = "#666666", face = "italic"),
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 10),
      strip.text   = element_text(size = 12, face = "bold")
    )
}


#' Plot outlier diagnostics
#'
#' @param ol Outlier diagnostics list.
#' @param col_point Point colour.
#' @param col_ref Reference line colour.
#' @param alpha_pt Point alpha.
#' @param max_points Optional maximum number of points to plot per residual type.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_outliers <- function(ol, col_point, col_ref, alpha_pt, max_points = 5000L) {

  lp <- ol$linear_predictor
  n <- length(lp)

  df_full <- data.frame(
    linear_predictor = rep(lp, 4L),
    residual         = c(ol$martingale, ol$deviance,
                         ol$log_odds, ol$normal_deviate),
    type             = rep(c("Martingale", "Deviance",
                             "Log-Odds", "Normal Deviate"),
                            each = n),
    stringsAsFactors = FALSE
  )

  df_plot <- df_full

  # Subsample non-flagged points for plotting in large datasets while strictly retaining all flagged points
  is_subsampled <- FALSE
  if (!is.null(max_points) && is.finite(max_points) && max_points > 0 && n > max_points) {
    set.seed(123)
    df_plot <- do.call(rbind, lapply(split(df_full, df_full$type), function(sub_df) {
      if (nrow(sub_df) <= max_points) return(sub_df)
      t_name <- unique(sub_df$type)
      thresh_val <- switch(t_name,
        "Deviance"       = 1.96,
        "Normal Deviate" = 1.96,
        "Log-Odds"       = 3.66,
        "Martingale"     = 2.0,
        2.0
      )
      is_extreme <- abs(sub_df$residual) > thresh_val
      n_extreme <- sum(is_extreme)
      if (n_extreme >= max_points) {
        # Strict preservation: retain all flagged cases; add zero normal cases
        sub_df[which(is_extreme), ]
      } else {
        n_sample <- max_points - n_extreme
        idx_normal <- which(!is_extreme)
        sampled_normal <- .sample_indices(idx_normal, n_sample)
        sub_df[sort(c(which(is_extreme), sampled_normal)), ]
      }
    }))
    is_subsampled <- TRUE
  }

  # Reference lines: +/-1.96 for deviance & normal deviate, +/-3.66 for log-odds
  ref_lines <- data.frame(
    type       = c("Deviance", "Deviance",
                   "Normal Deviate", "Normal Deviate",
                   "Log-Odds", "Log-Odds"),
    yintercept = c(1.96, -1.96, 1.96, -1.96, 3.66, -3.66),
    stringsAsFactors = FALSE
  )

  plot_caption <- if (is_subsampled) {
    paste0("Points subsampled to ", format(max_points, big.mark = ","),
           " observations per residual type (all flagged cases retained).")
  } else {
    NULL
  }

  ggplot(df_plot, aes(x = .data$linear_predictor, y = .data$residual)) +
    geom_point(colour = col_point, alpha = alpha_pt, size = 1) +
    geom_hline(
      data      = ref_lines,
      aes(yintercept = .data$yintercept),
      linetype  = "dashed",
      colour    = col_ref,
      linewidth = 0.5
    ) +
    facet_wrap(~ type, scales = "free_y") +
    labs(
      title   = "Outlier Assessment: Residual Diagnostics",
      x       = "Linear Predictor",
      y       = "Residual",
      caption = plot_caption
    ) +
    theme_minimal() +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 9, colour = "#666666", face = "italic"),
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 10),
      strip.text   = element_text(size = 12, face = "bold")
    )
}

#' Plot goodness-of-fit diagnostics
#'
#' @param gof Goodness-of-fit diagnostics list.
#' @param col_point Point colour.
#' @param col_ref Reference line colour.
#' @param alpha_pt Point alpha.
#' @param max_points Optional maximum number of points to plot.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_gof <- function(gof, col_point, col_ref, alpha_pt, max_points = 5000L) {

  df_full <- data.frame(
    x = gof$plot_x,
    y = gof$plot_y,
    stringsAsFactors = FALSE
  )

  # Ensure sorted by Cox-Snell residual so systematic index sampling corresponds to quantiles
  if (is.unsorted(df_full$x)) {
    df_full <- df_full[order(df_full$x), ]
  }

  df_plot <- df_full

  # Subsample points systematically across quantiles in large datasets
  is_subsampled <- FALSE
  if (!is.null(max_points) && is.finite(max_points) && max_points > 0 && nrow(df_full) > max_points) {
    idx <- round(seq(1L, nrow(df_full), length.out = max_points))
    df_plot <- df_full[idx, ]
    is_subsampled <- TRUE
  }

  plot_caption <- if (is_subsampled) {
    paste0("Points systematically subsampled to ", format(max_points, big.mark = ","),
           " observations.")
  } else {
    NULL
  }

  ggplot(df_plot, aes(x = .data$x, y = .data$y)) +
    geom_point(colour = col_point, alpha = alpha_pt, size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                colour = col_ref, linewidth = 0.8) +
    labs(
      title   = "Global Goodness-of-Fit: Cox-Snell Residuals",
      x       = "Cox-Snell Residual",
      y       = "Cumulative Hazard",
      caption = plot_caption
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.caption  = element_text(size = 9, colour = "#666666", face = "italic"),
      axis.title    = element_text(size = 12),
      axis.text     = element_text(size = 10)
    )
}
