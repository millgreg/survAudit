#' Plot survAudit Diagnostics
#'
#' Produces \code{ggplot2} diagnostic panels for a Cox PH model audit.
#' Four plot types are available: proportional hazards (Schoenfeld
#' residuals), functional form (martingale residuals), influence
#' diagnostics (DFBETAs and likelihood displacement), and outlier
#' assessment (multiple residual types).
#'
#' @param x An object of class \code{survAudit}.
#' @param which Character vector specifying which plots to produce.
#'   Valid values are \code{"ph"}, \code{"functional"},
#'   \code{"influence"}, \code{"outliers"}, and \code{"calibration"}.
#'   Defaults to all five.
#' @param ask Logical. If \code{TRUE} (the default in interactive
#'   sessions), the user is prompted between plots.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot} object (invisibly). If multiple plots are
#'   requested, a named list of \code{ggplot} objects is returned
#'   invisibly.
#'
#' @export
#'
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ trt + celltype + karno + age,
#'              data = veteran)
#' audit <- survAudit(fit)
#' plot(audit, which = "ph")
plot.survAudit <- function(x,
                           which = c("ph", "functional",
                                     "influence", "outliers", "calibration"),
                           ask = interactive(),
                           ...) {

  which <- match.arg(which, choices = c("ph", "functional",
                                        "influence", "outliers", "calibration"),
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
      plots[["ph"]] <- .plot_ph(x$ph, col_point, col_smooth,
                                col_ref, alpha_pt)
    }
  }

  # ── Functional form ────────────────────────────────────────
  if ("functional" %in% which) {
    if (is.null(x$functional_form) ||
        length(x$functional_form$results) == 0L) {
      message("Functional form diagnostics not available; ",
              "skipping 'functional' plot.")
    } else {
      plots[["functional"]] <- .plot_functional(
        x$functional_form, col_point, col_smooth, col_ref, alpha_pt
      )
    }
  }

  # ── Influence ──────────────────────────────────────────────
  if ("influence" %in% which) {
    if (is.null(x$influence)) {
      message("Influence diagnostics not available; ",
              "skipping 'influence' plot.")
    } else {
      plots[["influence"]] <- .plot_influence(
        x$influence, col_point, col_smooth, col_ref, col_thresh, alpha_pt
      )
    }
  }

  # ── Outliers ───────────────────────────────────────────────
  if ("outliers" %in% which) {
    if (is.null(x$outliers)) {
      message("Outlier diagnostics not available; ",
              "skipping 'outliers' plot.")
    } else {
      plots[["outliers"]] <- .plot_outliers(
        x$outliers, col_point, col_ref, alpha_pt
      )
    }
  }

  # ── Calibration ──────────────────────────────────────────────
  if ("calibration" %in% which) {
    if (is.null(x$calibration) || length(x$calibration$plot_x) == 0L) {
      message("Calibration diagnostics not available; ",
              "skipping 'calibration' plot.")
    } else {
      plots[["calibration"]] <- .plot_calibration(
        x$calibration, col_point, col_ref, alpha_pt
      )
    }
  }

  # ── Display ────────────────────────────────────────────────
  if (length(plots) == 0L) {
    message("No diagnostic plots available.")
    return(invisible(NULL))
  }

  if (isTRUE(ask)) {
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

#' Plot PH diagnostics (Schoenfeld residuals)
#'
#' @param ph PH diagnostics list from the survAudit object.
#' @param col_point Point colour.
#' @param col_smooth Smooth line colour.
#' @param col_ref Reference line colour.
#' @param alpha_pt Point alpha.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_ph <- function(ph, col_point, col_smooth, col_ref, alpha_pt) {

  zph <- ph$zph

  # Build long-format data.frame
  # zph$y is a matrix (time-points x covariates), zph$x is the time axis
  y_mat <- zph$y
  if (is.null(dim(y_mat))) {
    # Single covariate
    y_mat <- matrix(y_mat, ncol = 1L)
    colnames(y_mat) <- names(coef(zph))
    if (is.null(colnames(y_mat))) colnames(y_mat) <- "covariate"
  }

  time_vals <- zph$x
  n_t <- length(time_vals)
  n_vars <- ncol(y_mat)
  var_names <- colnames(y_mat)

  df <- data.frame(
    time     = rep(time_vals, times = n_vars),
    residual = as.vector(y_mat),
    variable = rep(var_names, each = n_t),
    stringsAsFactors = FALSE
  )

  # Subsample points for plotting in large datasets to prevent overplotting
  df_points <- df
  if (n_t > 2000) {
    set.seed(123)
    df_points <- do.call(rbind, lapply(split(df, df$variable), function(sub_df) {
      if (nrow(sub_df) > 2000) {
        sub_df[sample(seq_len(nrow(sub_df)), 2000), ]
      } else {
        sub_df
      }
    }))
  }

  ggplot(df, aes(x = .data$time, y = .data$residual)) +
    geom_point(data = df_points, colour = col_point, alpha = alpha_pt, size = 1) +
    geom_smooth(method = "loess", formula = y ~ x,
                se = TRUE, colour = col_smooth, linewidth = 0.8,
                fill = col_smooth, alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = col_ref,
               linewidth = 0.5) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(
      title = "Proportional Hazards Diagnostics: Scaled Schoenfeld Residuals",
      x     = paste0("Transformed Time (", ph$transform, ")"),
      y     = "Scaled Schoenfeld Residual"
    ) +
    theme_minimal() +
    theme(
      plot.title  = element_text(size = 12, face = "bold"),
      strip.text  = element_text(face = "bold")
    )
}


#' Plot functional form diagnostics
#'
#' @param ff Functional form diagnostics list.
#' @param col_point Point colour.
#' @param col_smooth Smooth line colour.
#' @param col_ref Reference line colour.
#' @param alpha_pt Point alpha.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_functional <- function(ff, col_point, col_smooth,
                             col_ref, alpha_pt) {

  dfs <- lapply(names(ff$results), function(vname) {
    res <- ff$results[[vname]]
    data.frame(
      covariate_value = res$covariate_values,
      residual        = res$residuals,
      variable        = vname,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, dfs)

  # Subsample points for plotting in large datasets to prevent overplotting
  df_points <- df
  # Check count per covariate
  if (nrow(df) > 0) {
    n_vars <- length(unique(df$variable))
    n_obs_per_var <- nrow(df) / n_vars
    if (n_obs_per_var > 2000) {
      set.seed(123)
      df_points <- do.call(rbind, lapply(split(df, df$variable), function(sub_df) {
        if (nrow(sub_df) > 2000) {
          sub_df[sample(seq_len(nrow(sub_df)), 2000), ]
        } else {
          sub_df
        }
      }))
    }
  }

  ggplot(df, aes(x = .data$covariate_value, y = .data$residual)) +
    geom_point(data = df_points, colour = col_point, alpha = alpha_pt, size = 1) +
    geom_smooth(method = "loess", formula = y ~ x,
                se = TRUE, colour = col_smooth, linewidth = 0.8,
                fill = col_smooth, alpha = 0.15) +
    geom_smooth(method = "lm", formula = y ~ x,
                se = FALSE, colour = col_ref,
                linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~ variable, scales = "free") +
    labs(
      title = "Functional Form Assessment: Martingale Residuals",
      x     = "Covariate Value",
      y     = "Martingale Residual"
    ) +
    theme_minimal() +
    theme(
      plot.title  = element_text(size = 12, face = "bold"),
      strip.text  = element_text(face = "bold")
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
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_influence <- function(inf, col_point, col_smooth,
                            col_ref, col_thresh, alpha_pt) {

  # DFBETAs in long format
  dfb <- inf$dfbetas
  n <- nrow(dfb)
  p <- ncol(dfb)
  var_names <- colnames(dfb)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(p))

  df_dfb <- data.frame(
    obs      = rep(seq_len(n), times = p),
    value    = as.vector(dfb),
    variable = rep(var_names, each = n),
    stringsAsFactors = FALSE
  )

  # Add likelihood displacement as another panel
  df_ld <- data.frame(
    obs      = seq_len(n),
    value    = inf$likelihood_displacement,
    variable = "Likelihood Displacement",
    stringsAsFactors = FALSE
  )

  df <- rbind(df_dfb, df_ld)

  threshold <- inf$threshold

  # Threshold data for DFBETAs panels only
  thresh_df <- data.frame(
    variable  = rep(var_names, each = 2L),
    yintercept = rep(c(threshold, -threshold), times = p),
    stringsAsFactors = FALSE
  )

  ggplot(df, aes(x = .data$obs, y = .data$value)) +
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
      title = "Influence Diagnostics",
      x     = "Observation Index",
      y     = "Value"
    ) +
    theme_minimal() +
    theme(
      plot.title  = element_text(size = 12, face = "bold"),
      strip.text  = element_text(face = "bold")
    )
}


#' Plot outlier diagnostics
#'
#' @param ol Outlier diagnostics list.
#' @param col_point Point colour.
#' @param col_ref Reference line colour.
#' @param alpha_pt Point alpha.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_outliers <- function(ol, col_point, col_ref, alpha_pt) {

  lp <- ol$linear_predictor
  n <- length(lp)

  df <- data.frame(
    linear_predictor = rep(lp, 4L),
    residual         = c(ol$martingale, ol$deviance,
                         ol$log_odds, ol$normal_deviate),
    type             = rep(c("Martingale", "Deviance",
                             "Log-Odds", "Normal Deviate"),
                           each = n),
    stringsAsFactors = FALSE
  )

  # Reference lines: ±1.96 for deviance & normal deviate, ±3.66 for log-odds
  ref_lines <- data.frame(
    type       = c("Deviance", "Deviance",
                   "Normal Deviate", "Normal Deviate",
                   "Log-Odds", "Log-Odds"),
    yintercept = c(1.96, -1.96, 1.96, -1.96, 3.66, -3.66),
    stringsAsFactors = FALSE
  )

  ggplot(df, aes(x = .data$linear_predictor, y = .data$residual)) +
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
      title = "Outlier Assessment: Residual Diagnostics",
      x     = "Linear Predictor",
      y     = "Residual"
    ) +
    theme_minimal() +
    theme(
      plot.title  = element_text(size = 12, face = "bold"),
      strip.text  = element_text(face = "bold")
    )
}

#' Plot calibration diagnostics
#'
#' @param cal Calibration diagnostics list.
#' @param col_point Point colour.
#' @param col_ref Reference line colour.
#' @param alpha_pt Point alpha.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_calibration <- function(cal, col_point, col_ref, alpha_pt) {

  df <- data.frame(
    x = cal$plot_x,
    y = cal$plot_y,
    stringsAsFactors = FALSE
  )

  ggplot(df, aes(x = .data$x, y = .data$y)) +
    geom_point(colour = col_point, alpha = alpha_pt, size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                colour = col_ref, linewidth = 0.8) +
    labs(
      title = "Overall Model Calibration",
      subtitle = "Nelson-Aalen Cumulative Hazard of Cox-Snell Residuals",
      x     = "Cox-Snell Residual",
      y     = "Cumulative Hazard"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, face = "italic")
    )
}
