#' Summarise a survAudit Object
#'
#' Produces a detailed summary of all diagnostic results for a Cox PH
#' model audit, including the full proportional hazards table,
#' per-covariate functional form assessment, top influential
#' observations, top outliers, data context, and the complete
#' assumption ontology.
#'
#' @param object An object of class \code{survAudit}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class \code{summary.survAudit}, printed via its
#'   own \code{\link{print.summary.survAudit}} method.
#'
#' @export
#'
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ trt + celltype + karno + age,
#'              data = veteran)
#' audit <- survAudit(fit)
#' summary(audit)
summary.survAudit <- function(object, ...) {

  out <- list(
    model_info      = object$model_info,
    data_context    = object$data_context,
    ph              = object$ph,
    functional_form = object$functional_form,
    influence       = object$influence,
    outliers        = object$outliers,
    epv             = object$epv,
    vif             = object$vif,
    gof             = object$gof,
    assumptions     = object$assumptions,
    alpha           = object$alpha,
    audit_time      = object$audit_time
  )

  # Pre-compute top influential observations (by max |dfbetas|)
  if (!is.null(object$influence)) {
    dfb <- object$influence$dfbetas
    if (!is.null(dfb) && nrow(dfb) > 0) {
      max_abs_dfb <- apply(abs(dfb), 1, max)
      n_show <- min(5L, length(max_abs_dfb))
      top_idx <- order(max_abs_dfb, decreasing = TRUE)[seq_len(n_show)]
      out$top_influence <- data.frame(
        obs            = top_idx,
        max_abs_dfbetas = max_abs_dfb[top_idx],
        stringsAsFactors = FALSE
      )
    }
  }

  # Pre-compute top deviance residuals
  if (!is.null(object$outliers)) {
    dev_r <- object$outliers$deviance
    n_show <- min(5L, length(dev_r))
    top_idx <- order(abs(dev_r), decreasing = TRUE)[seq_len(n_show)]
    out$top_deviance <- data.frame(
      obs              = top_idx,
      deviance_residual = dev_r[top_idx],
      stringsAsFactors  = FALSE
    )
  }

  class(out) <- "summary.survAudit"
  out
}


#' Print a summary.survAudit Object
#'
#' Renders the detailed diagnostic summary to the console.
#'
#' @param x An object of class \code{summary.survAudit}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return \code{x}, invisibly.
#'
#' @export
#'
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ trt + celltype + karno + age,
#'              data = veteran)
#' audit <- survAudit(fit)
#' s <- summary(audit)
#' print(s)
print.summary.survAudit <- function(x, ...) {

  alpha <- x$alpha

  # ═══════════════════════════════════════════════════════════════
  # Title
  # ═══════════════════════════════════════════════════════════════
  cat(.section_header("survAudit: Detailed Diagnostic Report"), "\n\n")
  cat("  Audit time: ", format(x$audit_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

  # ═══════════════════════════════════════════════════════════════
  # Data Context
  # ═══════════════════════════════════════════════════════════════
  cat(.section_header("Data Context"), "\n\n")

  if (!is.null(x$data_context)) {
    dc <- x$data_context
    cat("  Sample size:           ", dc$n, "\n")
    cat("  Events:                ", dc$n_events, "\n")
    cat("  Censored:              ", dc$n_censored,
        " (", sprintf("%.1f%%", dc$censoring_rate * 100), ")\n", sep = "")
    cat("  Time range:            ",
        sprintf("[%.2f, %.2f]", dc$time_range[1], dc$time_range[2]), "\n")
    cat("  Median time:           ",
        sprintf("%.2f", dc$time_median), "\n")
    cat("  Time IQR:              ",
        sprintf("[%.2f, %.2f]", dc$time_iqr[1], dc$time_iqr[2]), "\n")
    cat("  Tied event times:      ", dc$n_ties,
        " (", sprintf("%.1f%%", dc$tie_fraction * 100), ")\n", sep = "")
    if (!is.null(dc$counting_process)) {
      cat("  Counting process:      ",
          if (isTRUE(dc$counting_process)) "yes" else "no", "\n")
    }
    if (!is.null(dc$missing_data)) {
      if (is.character(dc$missing_data)) {
        cat("  Missing data:          ", dc$missing_data, "\n")
      } else if (is.data.frame(dc$missing_data)) {
        total_missing <- sum(dc$missing_data$n_missing)
        if (total_missing == 0) {
          cat("  Missing data:          none (complete case data)\n")
        } else {
          cat("  Missing data:\n")
          for (i in seq_len(nrow(dc$missing_data))) {
            row <- dc$missing_data[i, ]
            if (row$n_missing > 0) {
              cat(sprintf("    %s: %d missing (%.1f%%)\n", row$variable, row$n_missing, row$pct_missing))
            }
          }
        }
      }
    }
  } else {
    cat("  Data context not available.\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Model Info
  # ═══════════════════════════════════════════════════════════════
  cat("\n")
  cat(.section_header("Model Information"), "\n\n")

  mi <- x$model_info
  formula_str <- if (!is.null(mi$formula)) {
    deparse(mi$formula, width.cutoff = 500L)
  } else {
    "unknown"
  }
  cat("  Formula:      ", formula_str, "\n")
  cat("  Coefficients: ", mi$n_coef, "\n")
  cat("  Tie method:   ", if (!is.null(mi$method)) mi$method else "unknown",
      "\n")
  cat("  Converged:    ",
      if (isTRUE(mi$converged)) "yes" else "no / unknown", "\n")

  # ═══════════════════════════════════════════════════════════════
  # Global Goodness-of-Fit (Cox-Snell)
  # ═══════════════════════════════════════════════════════════════
  cat("\n")
  cat(.section_header("Global Goodness-of-Fit"), "\n\n")

  if (!is.null(x$gof)) {
    cat("  Cox-Snell residuals computed.\n")
    cat("  (See plot(audit, which = \"gof\") for the Nelson-Aalen cumulative hazard plot.)\n")
  } else {
    cat("  Goodness-of-fit diagnostics not available.\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Proportional Hazards
  # ═══════════════════════════════════════════════════════════════
  cat("\n")
  cat(.section_header("Proportional Hazards (cox.zph)"), "\n\n")

  if (!is.null(x$ph)) {
    cat("  Transform: ", x$ph$transform, "\n\n")

    tbl <- x$ph$table
    # Sort tbl by p-value ascending, keeping "GLOBAL" at the bottom
    if (!is.null(tbl) && nrow(tbl) > 0) {
      has_global <- "GLOBAL" %in% rownames(tbl)
      if (has_global) {
        covar_rows <- tbl[rownames(tbl) != "GLOBAL", , drop = FALSE]
        global_row <- tbl["GLOBAL", , drop = FALSE]
      } else {
        covar_rows <- tbl
        global_row <- NULL
      }
      
      if (nrow(covar_rows) > 0) {
        sort_ord <- order(covar_rows[, "p"], na.last = TRUE)
        covar_rows <- covar_rows[sort_ord, , drop = FALSE]
      }
      
      if (!is.null(global_row)) {
        tbl_sorted <- rbind(covar_rows, global_row)
      } else {
        tbl_sorted <- covar_rows
      }
    } else {
      tbl_sorted <- tbl
    }

    # Print full table with aligned columns
    cat("  ", sprintf("%-25s %8s %6s %10s", "Variable", "chisq",
                      "df", "p"), "\n")
    cat("  ", .rule_line(52), "\n")
    for (rn in rownames(tbl_sorted)) {
      cat("  ", sprintf("%-25s %8.3f %6d %10s",
                        rn, tbl_sorted[rn, "chisq"], as.integer(tbl_sorted[rn, "df"]),
                        .format_p(tbl_sorted[rn, "p"])), "\n")
    }

    # Trends
    if (!is.null(x$ph$trends) && nrow(x$ph$trends) > 0L) {
      trends <- x$ph$trends[x$ph$trends$direction != "none", ]
      if (nrow(trends) > 0L) {
        cat("\n  Detected trends:\n")
        for (i in seq_len(nrow(trends))) {
          cat("    ", trends$variable[i], ": ",
              trends$direction[i], " effect over time\n", sep = "")
        }
      }
    }
  } else {
    cat("  Proportional hazards diagnostics not available.\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Functional Form
  # ═══════════════════════════════════════════════════════════════
  cat("\n")
  cat(.section_header("Functional Form Assessment"), "\n\n")

  if (!is.null(x$functional_form) &&
      length(x$functional_form$results) > 0L) {
    cat("  Continuous variables assessed: ",
        paste(x$functional_form$continuous_vars, collapse = ", "), "\n\n")
    for (vname in names(x$functional_form$results)) {
      res <- x$functional_form$results[[vname]]
      status <- if (isTRUE(res$departure_detected)) {
        "DEPARTURE detected"
      } else {
        "no departure detected"
      }
      cat("  ", vname, ": ", status, "\n", sep = "")
    }
  } else {
    cat("  Functional form diagnostics not available.\n")
    cat("  (Requires continuous covariates and model data.)\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Influence Diagnostics
  # ═══════════════════════════════════════════════════════════════
  cat("\n")
  cat(.section_header("Influence Diagnostics"), "\n\n")

  if (!is.null(x$influence)) {
    inf <- x$influence
    cat("  Threshold (2/sqrt(n)): ",
        sprintf("%.4f", inf$threshold), "\n")
        
    if (!is.null(x$data_context) && x$data_context$n > 10000) {
      cat("  *Note: In large datasets (N > 10,000), this threshold may be overly sensitive.\n")
      cat("         Evaluate the absolute magnitude of coefficient change.*\n")
    }
    
    cat("  Max |dfbetas|:         ",
        sprintf("%.4f", inf$max_dfbetas$value),
        " (obs #", inf$max_dfbetas$obs,
        ", ", inf$max_dfbetas$variable, ")\n", sep = "")
    cat("  Flagged observations:  ",
        length(inf$flagged_obs), "\n")

    # Top influential obs table
    if (!is.null(x$top_influence)) {
      cat("\n  Top 5 most influential observations (by max |dfbetas|):\n\n")
      cat("  ", sprintf("%-8s %14s",
                        "Obs", "Max |dfbetas|"), "\n")
      cat("  ", .rule_line(24), "\n")
      ti <- x$top_influence
      for (i in seq_len(nrow(ti))) {
        cat("  ", sprintf("%-8d %14.4f",
                          ti$obs[i], ti$max_abs_dfbetas[i]), "\n")
      }
    }
  } else {
    cat("  Influence diagnostics not available.\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Outlier Assessment
  # ═══════════════════════════════════════════════════════════════
  cat("\n")
  cat(.section_header("Outlier Assessment"), "\n\n")

  if (!is.null(x$outliers)) {
    ol <- x$outliers
    
    cat("  *Note: Deviance residuals are the most robust symmetric measure for outliers,\n")
    cat("         especially in heavily censored datasets.*\n\n")
    
    cat("  Flagged deviance residuals (|resid| > 1.96):       ",
        length(ol$flagged$deviance), "\n")
    cat("  Flagged log-odds residuals (|resid| > 3.66):       ",
        length(ol$flagged$log_odds), "\n")
    cat("  Flagged normal deviate residuals (|resid| > 1.96): ",
        length(ol$flagged$normal_deviate), "\n")

    # Top deviance residuals
    if (!is.null(x$top_deviance)) {
      cat("\n  Top 5 largest deviance residuals:\n\n")
      cat("  ", sprintf("%-8s %14s", "Obs", "Deviance Resid"), "\n")
      cat("  ", .rule_line(24), "\n")
      td <- x$top_deviance
      for (i in seq_len(nrow(td))) {
        cat("  ", sprintf("%-8d %14.4f",
                          td$obs[i], td$deviance_residual[i]), "\n")
      }
    }
  } else {
    cat("  Outlier diagnostics not available.\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Event Sufficiency (Events Per Variable)
  # ═══════════════════════════════════════════════════════════════
  cat("\n")
  cat(.section_header("Event Sufficiency (Events Per Variable - EPV)"), "\n\n")

  if (!is.null(x$epv)) {
    cat("  Events:                       ", x$epv$n_events, "\n")
    cat("  Estimated parameters (coefs): ", x$epv$n_parameters, "\n")
    cat("  Events Per Variable (EPV):    ", sprintf("%.1f", x$epv$ratio), "\n")
    cat("  Assessment:                   ", x$epv$classification, "\n")
  } else {
    cat("  EPV diagnostics not available.\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Collinearity Diagnostics (VIF)
  # ═══════════════════════════════════════════════════════════════
  cat("\n")
  cat(.section_header("Collinearity Diagnostics (VIF)"), "\n\n")
  cat("  Threshold: VIF > 5 (or scaled GVIF^(1/(2*Df)) > 2.236 for multi-Df terms)\n\n")

  if (!is.null(x$vif) && !is.null(x$vif$vif)) {
    vif_mat <- x$vif$vif
    all_df1 <- all(vif_mat[, "Df"] == 1)

    if (all_df1) {
      cat("  ", sprintf("%-25s %8s", "Term", "VIF"), "\n")
      cat("  ", .rule_line(35), "\n")
      for (rn in rownames(vif_mat)) {
        cat("  ", sprintf("%-25s %8.3f", rn, vif_mat[rn, "GVIF"]), "\n")
      }
    } else {
      cat("  ", sprintf("%-25s %8s %6s %18s", "Term", "GVIF", "Df", "GVIF^(1/(2*Df))"), "\n")
      cat("  ", .rule_line(60), "\n")
      for (rn in rownames(vif_mat)) {
        cat("  ", sprintf("%-25s %8.3f %6d %18.3f",
                          rn, vif_mat[rn, "GVIF"], as.integer(vif_mat[rn, "Df"]),
                          vif_mat[rn, "GVIF^(1/(2*Df))"]), "\n")
      }
    }

    if (length(x$vif$flagged) > 0) {
      cat("\n  Flagged terms (VIF > 5 or GVIF^(1/(2*Df)) > 2.236):\n")
      cat("    ", paste(x$vif$flagged, collapse = ", "), "\n")
    } else {
      cat("\n  No collinearity issues detected (all thresholds satisfied).\n")
    }
  } else if (!is.null(x$vif) && !is.null(x$vif$message)) {
    cat("  VIF computation skipped: ", x$vif$message, "\n")
  } else {
    cat("  Collinearity diagnostics not available.\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Full Assumption Classification
  # ═══════════════════════════════════════════════════════════════
  cat("\n")
  cat(.section_header("Assumption Classification"), "\n\n")

  if (!is.null(x$assumptions)) {
    .print_assumption_group("Non-Identifiable",
                            x$assumptions$non_identifiable)
    .print_assumption_group("Partially Assessable",
                            x$assumptions$partially_assessable)
    .print_assumption_group("Statistically Assessable",
                            x$assumptions$assessable)
  } else {
    cat("  Assumption classification not available.\n")
  }

  cat("\n")
  invisible(x)
}


#' Print a group of assumptions
#'
#' Helper to print a named group of assumptions from the ontology.
#'
#' @param group_name Character string naming the group.
#' @param items List of assumption items.
#' @keywords internal
.print_assumption_group <- function(group_name, items) {
  cat("  ", group_name, ":\n", sep = "")

  if (is.null(items) || length(items) == 0L) {
    cat("    (none)\n\n")
    return(invisible(NULL))
  }

  for (item in items) {
    label <- if (!is.null(item$label)) item$label else item$name
    label <- if (is.null(label)) "unnamed" else label
    cat("    - ", label, sep = "")

    # Evidence summary if available
    if (!is.null(item$evidence_summary)) {
      cat(": ", item$evidence_summary, sep = "")
    }

    # Justification for non-identifiable
    if (!is.null(item$justification)) {
      cat(" [justified: ", item$justification, "]", sep = "")
    }

    cat("\n")
  }
  cat("\n")
}
