# Assumption framework builder for survAudit

#' Build the structured assumption classification
#'
#' Classifies Cox model assumptions into three categories:
#' \strong{assessable} (testable from data), \strong{partially assessable}
#' (partially informed by data), and \strong{non-identifiable} (cannot be
#' verified statistically). Each assessable assumption includes a concise
#' evidence summary derived from the diagnostic results.
#'
#' @param ph List of proportional hazards diagnostics from
#'   \code{.compute_ph_diagnostics}.
#' @param functional_form List of functional form diagnostics from
#'   \code{.compute_functional_form}.
#' @param influence List of influence diagnostics from
#'   \code{.compute_influence}.
#' @param outliers List of outlier diagnostics from
#'   \code{.compute_outliers}.
#' @param epv List of events-per-variable results from
#'   \code{.compute_epv}.
#' @param data_context List of data context from
#'   \code{.compute_data_context}.
#' @param alpha Numeric significance level.
#'
#' @return A list with three sublists: \code{assessable},
#'   \code{partially_assessable}, and \code{non_identifiable}.
#'
#' @keywords internal
#' @noRd
.build_assumptions <- function(ph, functional_form, influence, outliers,
                               epv, data_context, alpha, vif = NULL) {

  # --- Assessable assumptions ---

  # Proportional hazards
  ph_evidence <- .summarize_ph(ph, alpha)
  proportional_hazards <- list(
    name = "proportional_hazards",
    type = "assessable",
    evidence_summary = ph_evidence,
    diagnostic_key = "ph"
  )

  # Functional form
  ff_evidence <- .summarize_functional_form(functional_form)
  functional_form_assumption <- list(
    name = "functional_form",
    type = "partially_assessable",
    evidence_summary = ff_evidence,
    diagnostic_key = "functional_form"
  )

  # Influence stability
  inf_evidence <- .summarize_influence(influence)
  influence_stability <- list(
    name = "influence_stability",
    type = "assessable",
    evidence_summary = inf_evidence,
    diagnostic_key = "influence"
  )

  # Event sufficiency
  epv_evidence <- .summarize_epv(epv)
  event_sufficiency <- list(
    name = "event_sufficiency",
    type = "assessable",
    evidence_summary = epv_evidence,
    diagnostic_key = "epv"
  )

  # Collinearity
  vif_evidence <- .summarize_vif(vif)
  collinearity <- list(
    name = "collinearity",
    label = "Collinearity",
    type = "assessable",
    evidence_summary = vif_evidence,
    diagnostic_key = "vif"
  )

  assessable <- list(
    proportional_hazards = proportional_hazards,
    influence_stability = influence_stability,
    collinearity = collinearity,
    event_sufficiency = event_sufficiency
  )

  # --- Partially assessable assumptions ---

  # Outlier impact
  outlier_evidence <- .summarize_outliers(outliers)
  outlier_impact <- list(
    name = "outlier_impact",
    type = "partially_assessable",
    evidence_summary = outlier_evidence,
    diagnostic_key = "outliers"
  )

  # Missing data
  missing_evidence <- .summarize_missing(data_context)
  missing_data <- list(
    name = "missing_data",
    type = "partially_assessable",
    evidence_summary = missing_evidence,
    diagnostic_key = "data_context"
  )

  partially_assessable <- list(
    functional_form = functional_form_assumption,
    outlier_impact = outlier_impact,
    missing_data = missing_data
  )

  # --- Non-identifiable assumptions ---

  independent_censoring <- list(
    name = "independent_censoring",
    label = "Independent censoring",
    type = "non_identifiable",
    description = paste0(
      "Censoring mechanism is independent of the event process ",
      "conditional on covariates. Cannot be tested from observed data; ",
      "requires study design justification (Tsiatis, 1975)."
    ),
    justification = NULL
  )

  unmeasured_confounding <- list(
    name = "unmeasured_confounding",
    label = "Absence of unmeasured confounding",
    type = "non_identifiable",
    description = paste0(
      "No unmeasured variables confound the covariate-outcome ",
      "relationship. Cannot be verified statistically; requires domain ",
      "knowledge and study design arguments."
    ),
    justification = NULL
  )

  non_identifiable <- list(
    independent_censoring = independent_censoring,
    unmeasured_confounding = unmeasured_confounding
  )

  list(
    non_identifiable = non_identifiable,
    partially_assessable = partially_assessable,
    assessable = assessable
  )
}


# --- Helper functions for evidence summaries ---

#' Summarize proportional hazards evidence
#' @param ph PH diagnostics list.
#' @param alpha Significance level.
#' @return Character string.
#' @keywords internal
#' @noRd
.summarize_ph <- function(ph, alpha) {
  global_p <- ph$global_p
  p_str <- .format_p(global_p)

  # Count per-covariate violations
  table_mat <- ph$table
  # Exclude the GLOBAL row (last row) for per-covariate checks
  if (nrow(table_mat) > 1) {
    covariate_rows <- table_mat[-nrow(table_mat), , drop = FALSE]
    n_flagged <- sum(covariate_rows[, "p"] < alpha, na.rm = TRUE)
    n_total <- nrow(covariate_rows)
  } else {
    # Single covariate model: the table has only one row which is both
    # the covariate and the global test
    n_flagged <- sum(table_mat[, "p"] < alpha, na.rm = TRUE)
    n_total <- nrow(table_mat)
  }

  global_str <- if (global_p < alpha) {
    paste0("Global cox.zph test significant (p = ", p_str, ")")
  } else {
    paste0("Global cox.zph test not significant (p = ", p_str, ")")
  }

  covar_str <- if (n_flagged > 0) {
    paste0(n_flagged, " of ", n_total, " covariate(s) flagged at alpha = ", alpha, ".")
  } else {
    paste0("no covariate-level violations detected at alpha = ", alpha, ".")
  }

  paste0(global_str, "; ", covar_str)
}

#' Summarize functional form evidence
#' @param functional_form Functional form diagnostics list.
#' @return Character string.
#' @keywords internal
#' @noRd
.summarize_functional_form <- function(functional_form) {
  results <- functional_form$results
  vars <- functional_form$continuous_vars

  if (length(vars) == 0 || length(results) == 0) {
    return("No continuous covariates to assess for functional form.")
  }

  paste0(
    "Assessment of ", length(results), " continuous covariate(s) is inherently visual. ",
    "Use plot(audit, which = 'functional') to inspect martingale residual LOESS smooths for non-linearity."
  )
}

#' Summarize influence evidence
#' @param influence Influence diagnostics list.
#' @return Character string.
#' @keywords internal
#' @noRd
.summarize_influence <- function(influence) {
  max_dfb <- influence$max_dfbetas
  n_flagged <- length(influence$flagged_obs)

  paste0(
    "Max |dfbetas| = ", formatC(max_dfb$value, format = "f", digits = 4),
    " (obs ", max_dfb$obs, ", variable '", max_dfb$variable, "'); ",
    n_flagged, " observation(s) exceed dfbetas threshold."
  )
}

#' Summarize EPV evidence
#' @param epv EPV diagnostics list.
#' @return Character string.
#' @keywords internal
#' @noRd
.summarize_epv <- function(epv) {
  paste0(
    "Events-per-variable ratio = ",
    formatC(epv$ratio, format = "f", digits = 1),
    " (", epv$n_events, " events / ", epv$n_parameters,
    " parameters); classified as '", epv$classification, "'."
  )
}

#' Summarize outlier evidence
#' @param outliers Outlier diagnostics list.
#' @return Character string.
#' @keywords internal
#' @noRd
.summarize_outliers <- function(outliers) {
  n_dev <- length(outliers$flagged$deviance)
  n_lo <- length(outliers$flagged$log_odds)
  n_nd <- length(outliers$flagged$normal_deviate)

  paste0(
    "Flagged observations: ", n_dev, " by deviance residuals, ",
    n_lo, " by log-odds, ", n_nd, " by normal deviate."
  )
}

#' Summarize missing data evidence
#' @param data_context Data context list.
#' @return Character string.
#' @keywords internal
#' @noRd
.summarize_missing <- function(data_context) {
  md <- data_context$missing_data

  if (is.character(md) && md == "not available") {
    return("Missing data status not available (original data not provided).")
  }

  if (is.data.frame(md) && nrow(md) == 0) {
    return("No covariates identified to check for missing data.")
  }

  total_missing <- sum(md$n_missing)
  n_vars_missing <- sum(md$n_missing > 0)

  if (total_missing == 0) {
    paste0("Complete data: no missing values across ", nrow(md), " covariate(s).")
  } else {
    paste0(
      n_vars_missing, " of ", nrow(md),
      " covariate(s) have missing values (",
      total_missing, " total missing entries)."
    )
  }
}

#' Summarize collinearity evidence
#' @param vif VIF diagnostics list.
#' @return Character string.
#' @keywords internal
#' @noRd
.summarize_vif <- function(vif) {
  if (is.null(vif) || is.null(vif$vif)) {
    return("Collinearity diagnostics not available.")
  }

  vifs <- vif$vif
  if (is.vector(vifs)) {
    max_vif <- max(vifs, na.rm = TRUE)
  } else if (is.matrix(vifs)) {
    # If Df > 1, use GVIF^(1/(2*Df)), otherwise use GVIF
    vals <- ifelse(vifs[, 2] > 1, vifs[, 3], vifs[, 1])
    max_vif <- max(vals, na.rm = TRUE)
  } else {
    return("Collinearity diagnostics not available.")
  }

  n_flagged <- length(vif$flagged)
  if (n_flagged == 0) {
    paste0("No collinearity detected; max VIF/scaled GVIF = ",
           formatC(max_vif, format = "f", digits = 2), ".")
  } else {
    paste0("Collinearity detected for ", n_flagged, " term(s): ",
           paste(vif$flagged, collapse = ", "),
           " (max VIF/scaled GVIF = ",
           formatC(max_vif, format = "f", digits = 2), ").")
  }
}
