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
  cat("=================================================================\n")
  cat("             survAudit: Detailed Diagnostic Report               \n")
  cat("=================================================================\n\n")

  # ═══════════════════════════════════════════════════════════════
  # Data Context & Model Info
  # ═══════════════════════════════════════════════════════════════
  cat(.section_header("Data Context & Model Info"), "\n\n")
  
  mi <- x$model_info
  formula_str <- if (!is.null(mi$formula)) {
    deparse(mi$formula, width.cutoff = 500L)
  } else {
    "unknown"
  }
  cat("  Formula:      ", formula_str, "\n")
  
  if (!is.null(x$data_context)) {
    dc <- x$data_context
    cat("  Sample size:  ", dc$n, "  |  Events: ", dc$n_events, 
        "  |  Censored: ", dc$n_censored,
        " (", sprintf("%.1f%%", dc$censoring_rate * 100), ")\n", sep = "")
    
    missing_str <- "none"
    if (!is.null(dc$missing_data) && is.data.frame(dc$missing_data)) {
      missing_rows <- dc$missing_data[dc$missing_data$n_missing > 0, ]
      if (nrow(missing_rows) > 0) {
        missing_parts <- lapply(seq_len(nrow(missing_rows)), function(i) {
          sprintf("%s: %d missing (%.1f%%)", missing_rows$variable[i], missing_rows$n_missing[i], missing_rows$pct_missing[i])
        })
        missing_str <- paste(unlist(missing_parts), collapse = ", ")
      }
    } else if (is.character(dc$missing_data)) {
      missing_str <- dc$missing_data
    }
    
    cat("  Missing data: ", missing_str, 
        "  |  Tied events: ", dc$n_ties,
        " (", sprintf("%.1f%%", dc$tie_fraction * 100), ")\n", sep = "")
  } else {
    cat("  Data context not available.\n")
  }
  cat("\n")

  # ═══════════════════════════════════════════════════════════════
  # 1. Non-Identifiable Assumptions
  # ═══════════════════════════════════════════════════════════════
  cat(.section_header("1. Non-Identifiable Assumptions"), "\n\n")
  cat("  These cannot be tested statistically and require domain knowledge:\n")
  if (!is.null(x$assumptions) && !is.null(x$assumptions$non_identifiable)) {
    for (item in x$assumptions$non_identifiable) {
      label <- if (!is.null(item$label)) item$label else item$name
      cat("    - ", label, "\n", sep = "")
    }
  } else {
    cat("    - Independent censoring\n")
    cat("    - Absence of unmeasured confounding\n")
  }
  cat("\n")

  # ═══════════════════════════════════════════════════════════════
  # 2. Partially Assessable Assumptions
  # ═══════════════════════════════════════════════════════════════
  cat(.section_header("2. Partially Assessable Assumptions"), "\n\n")

  # --- Goodness of Fit ---
  cat("  [ Global Goodness-of-Fit ]\n")
  if (!is.null(x$gof)) {
    cat("  Cox-Snell residuals computed.\n")
    cat("  (Use plot(audit, which = \"gof\") to visually assess macro-calibration.)\n")
  } else {
    cat("  Goodness-of-fit diagnostics not available.\n")
  }
  cat("\n")

  # --- Functional Form ---
  cat("  [ Functional Form Assessment ]\n")
  if (!is.null(x$functional_form) && length(x$functional_form$results) > 0L) {
    cat("  Continuous variables assessed: ", paste(x$functional_form$continuous_vars, collapse = ", "), "\n")
    cat("  (Use plot(audit, which = \"functional\") to visually inspect non-linearity.)\n")
  } else {
    cat("  Functional form diagnostics not available.\n")
  }
  cat("\n")

  # --- Influence Diagnostics ---
  cat("  [ Influence Diagnostics ]\n")
  if (!is.null(x$influence)) {
    inf <- x$influence
    cat("  Threshold (2/sqrt(n)): ", sprintf("%.4f", inf$threshold), "\n")
    cat("  Flagged observations (|dfbetas| > threshold): ", length(inf$flagged_obs), "\n")
    
    if (!is.null(x$top_influence)) {
      cat("\n  Top 5 most influential observations (by max |dfbetas|):\n")
      cat("  ", sprintf("%-8s %14s", "Obs", "Max |dfbetas|"), "\n")
      cat("  ", .rule_line(24), "\n")
      ti <- x$top_influence
      for (i in seq_len(nrow(ti))) {
        cat("  ", sprintf("%-8d %14.4f", ti$obs[i], ti$max_abs_dfbetas[i]), "\n")
      }
    }
    cat("\n  (Use plot(audit, which = \"influence\") to visually inspect highly influential cases.)\n")
  } else {
    cat("  Influence diagnostics not available.\n")
  }
  cat("\n")

  # --- Outliers ---
  cat("  [ Outlier Assessment ]\n")
  if (!is.null(x$outliers)) {
    ol <- x$outliers
    cat("  Flagged deviance residuals (|resid| > 1.96):       ", length(ol$flagged$deviance), "\n")
    cat("  Flagged log-odds residuals (|resid| > 3.66):       ", length(ol$flagged$log_odds), "\n")
    if (!is.null(x$top_deviance)) {
      cat("\n  Top 5 largest deviance residuals:\n")
      cat("  ", sprintf("%-8s %14s", "Obs", "Deviance Resid"), "\n")
      cat("  ", .rule_line(24), "\n")
      td <- x$top_deviance
      for (i in seq_len(nrow(td))) {
        cat("  ", sprintf("%-8d %14.4f", td$obs[i], td$deviance_residual[i]), "\n")
      }
    }
    cat("\n  (Use plot(audit, which = \"outliers\") to visually inspect distributions.)\n")
  } else {
    cat("  Outlier diagnostics not available.\n")
  }
  cat("\n")

  # ═══════════════════════════════════════════════════════════════
  # 3. Statistically Assessable Assumptions
  # ═══════════════════════════════════════════════════════════════
  cat(.section_header("3. Statistically Assessable Assumptions"), "\n\n")

  # --- Proportional Hazards ---
  cat("  [ Proportional Hazards ]\n")
  if (!is.null(x$ph)) {
    tbl <- x$ph$table
    
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
    cat("  ", sprintf("%-25s %8s %6s %10s", "Variable", "chisq", "df", "p"), "\n")
    cat("  ", .rule_line(52), "\n")
    for (rn in rownames(tbl_sorted)) {
      cat("  ", sprintf("%-25s %8.3f %6d %10s",
                        rn, tbl_sorted[rn, "chisq"], as.integer(tbl_sorted[rn, "df"]),
                        .format_p(tbl_sorted[rn, "p"])), "\n")
    }

    # Summary text
    if (!is.null(global_row)) {
      global_p_str <- .format_p(global_row[1, "p"])
      is_sig <- if (global_row[1, "p"] < alpha) "significant" else "not significant"
      cat(sprintf("\n  Global test: %s (p = %s)\n", is_sig, global_p_str))
    }
    
    if (!is.null(nrow(covar_rows)) && nrow(covar_rows) > 0) {
      flagged_vars <- rownames(covar_rows)[covar_rows[, "p"] < alpha]
      if (length(flagged_vars) > 0) {
        cat(sprintf("  Flagged covariates (p < %s): %s\n", alpha, paste(flagged_vars, collapse = ", ")))
      } else {
        cat(sprintf("  Flagged covariates (p < %s): none\n", alpha))
      }
    } else if (!is.null(tbl)) {
      # Single covariate model fallback
      flagged_vars <- rownames(tbl)[tbl[, "p"] < alpha]
      if (length(flagged_vars) > 0) {
        cat(sprintf("  Flagged covariates (p < %s): %s\n", alpha, paste(flagged_vars, collapse = ", ")))
      } else {
        cat(sprintf("  Flagged covariates (p < %s): none\n", alpha))
      }
    }
    
    cat("  (Use plot(audit, which = \"ph\") to visually assess time-varying effects.)\n")
  } else {
    cat("  Proportional hazards diagnostics not available.\n")
  }
  cat("\n")

  # --- EPV ---
  cat("  [ Event Sufficiency (EPV) ]\n")
  if (!is.null(x$epv)) {
    cat("  Events Per Variable (EPV): ", sprintf("%.1f", x$epv$ratio), "\n")
    cat("  Assessment:                ", x$epv$classification, " (target: >= 20)\n")
  } else {
    cat("  EPV diagnostics not available.\n")
  }
  cat("\n")

  # --- Collinearity ---
  cat("  [ Collinearity Diagnostics (VIF) ]\n")
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
      cat("\n  No collinearity issues detected.\n")
    }
  } else if (!is.null(x$vif) && !is.null(x$vif$message)) {
    cat("  VIF computation skipped: ", x$vif$message, "\n")
  } else {
    cat("  Collinearity diagnostics not available.\n")
  }

  cat("\n")
  invisible(x)
}
