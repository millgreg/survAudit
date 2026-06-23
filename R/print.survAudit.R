#' Print a survAudit Object
#'
#' Displays a compact, human-readable summary of a Cox PH model
#' diagnostic audit, including model overview, assessable assumptions
#' with flagged covariates, and non-identifiable assumptions that
#' require analyst justification.
#'
#' @param x An object of class \code{survAudit}.
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
#' print(audit)
print.survAudit <- function(x, ...) {

  alpha <- x$alpha

  # ── Header ─────────────────────────────────────────────────────
  cat(.section_header("survAudit: Cox PH Model Diagnostic Audit"), "\n\n")

  # Model formula
  formula_str <- if (!is.null(x$model_info$formula)) {
    deparse(x$model_info$formula, width.cutoff = 500L)
  } else {
    "unknown"
  }
  cat("  Model:      ", formula_str, "\n")

  # Events / sample size
  if (!is.null(x$data_context)) {
    dc <- x$data_context
    cens_pct <- sprintf("%.1f", dc$censoring_rate * 100)
    cat("  Events:     ", dc$n_events, " / ", dc$n,
        " (", cens_pct, "% censored)\n", sep = "")
  } else {
    cat("  Events:      not available\n")
  }

  # EPV
  if (!is.null(x$epv)) {
    cat("  EPV:        ", sprintf("%.1f", x$epv$ratio),
        " (", x$epv$classification, ")\n", sep = "")
  } else {
    cat("  EPV:         not available\n")
  }

  # ── Non-Identifiable Assumptions ───────────────────────────────
  cat("\n")
  cat(.section_header(
    "Non-Identifiable Assumptions (analyst justification required)"
  ), "\n\n")

  ni_items <- .get_non_identifiable_items(x)
  for (item in ni_items) {
    mark <- if (is.null(item$justification)) "[ ]" else "[x]"
    cat("  ", mark, " ", item$label, "\n", sep = "")
  }

  # ── Partially Assessable Assumptions ───────────────────────────
  cat("\n")
  cat(.section_header("Partially Assessable Assumptions"), "\n\n")

  # -- Outlier Impact --
  cat("  Outlier Impact\n")
  if (!is.null(x$outliers)) {
    cat("    ", .summarize_outliers(x$outliers), "\n", sep = "")
  } else {
    cat("    not available\n")
  }

  # -- Missing Data --
  cat("  Missing Data\n")
  if (!is.null(x$data_context)) {
    cat("    ", .summarize_missing(x$data_context), "\n", sep = "")
  } else {
    cat("    not available\n")
  }

  # -- Functional Form --
  cat("  Functional Form (visual assessment)\n")
  if (!is.null(x$functional_form) &&
      length(x$functional_form$results) > 0L) {
    for (vname in names(x$functional_form$results)) {
      res <- x$functional_form$results[[vname]]
      if (isTRUE(res$departure_detected)) {
        cat("    ", vname, ": possible departure from linearity\n",
            sep = "")
      } else {
        cat("    ", vname, ": no apparent departure from linearity\n",
            sep = "")
      }
    }
  } else {
    cat("    not available\n")
  }

  # ── Statistically Assessable Assumptions ───────────────────────
  cat("\n")
  cat(.section_header("Statistically Assessable Assumptions"), "\n\n")

  # -- Proportional Hazards --
  cat("  Proportional Hazards\n")
  if (!is.null(x$ph)) {
    cat("    Global test: p = ", .format_p(x$ph$global_p), "\n", sep = "")

    # Per-covariate results where p < alpha
    tbl <- x$ph$table
    # Exclude GLOBAL row
    covar_rows <- rownames(tbl)
    covar_rows <- covar_rows[covar_rows != "GLOBAL"]

    found_any <- FALSE
    for (cv in covar_rows) {
      p_val <- tbl[cv, "p"]
      if (!is.na(p_val) && p_val < alpha) {
        cat("    ", cv, ": p = ", .format_p(p_val), "\n", sep = "")
        found_any <- TRUE
      }
    }
    if (!found_any) {
      cat("    No covariate-level violations detected (p >= ", alpha, ")\n", sep = "")
    }
  } else {
    cat("    not available\n")
  }

  cat("\n")

  # -- Influence Stability --
  cat("  Influence Stability\n")
  if (!is.null(x$influence)) {
    inf <- x$influence
    cat("    Max |dfbetas|: ",
        sprintf("%.2f", inf$max_dfbetas$value),
        " (obs #", inf$max_dfbetas$obs,
        ", covariate: ", inf$max_dfbetas$variable, ")\n", sep = "")
  } else {
    cat("    not available\n")
  }

  cat("\n")

  # -- Event Sufficiency --
  cat("  Event Sufficiency\n")
  if (!is.null(x$epv)) {
    cat("    EPV = ", sprintf("%.1f", x$epv$ratio),
        " (", x$epv$classification, ")\n", sep = "")
  } else {
    cat("    not available\n")
  }

  cat("\n")

  # -- Collinearity --
  cat("  Collinearity\n")
  if (!is.null(x$vif) && !is.null(x$vif$vif)) {
    cat("    ", .summarize_vif(x$vif), "\n", sep = "")
  } else if (!is.null(x$vif) && !is.null(x$vif$message)) {
    cat("    ", x$vif$message, "\n", sep = "")
  } else {
    cat("    not available\n")
  }

  cat("\n  Use summary() for detailed diagnostics.",
      " Use plot() for visual assessment.\n", sep = "")

  invisible(x)
}


#' Retrieve non-identifiable assumption items
#'
#' Returns the list of non-identifiable assumptions with their labels
#' and current justification status, drawn from the assumption classification
#' or sensible defaults.
#'
#' @param x A \code{survAudit} object.
#' @return A list of lists, each with elements \code{label} and
#'   \code{justification}.
#' @keywords internal
.get_non_identifiable_items <- function(x) {

  # Default non-identifiable assumptions
  defaults <- list(
    list(label = "Independent censoring",            justification = NULL),
    list(label = "Absence of unmeasured confounding", justification = NULL)
  )

  # If the assumption classification has non_identifiable entries, use them
 if (!is.null(x$assumptions) &&
      !is.null(x$assumptions$non_identifiable) &&
      length(x$assumptions$non_identifiable) > 0L) {
    ni <- x$assumptions$non_identifiable
    items <- lapply(ni, function(a) {
      list(
        label         = if (!is.null(a$label)) a$label else a$name,
        justification = a$justification
      )
    })
    return(items)
  }

  defaults
}
