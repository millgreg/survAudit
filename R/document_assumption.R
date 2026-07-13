#' Document Qualitative Justifications for Assumptions
#'
#' A helper function to record qualitative, text-based justifications for
#' non-identifiable model assumptions directly onto the \code{survAudit} object.
#' This allows the clinical reasoning behind untestable assumptions (e.g., 
#' independent censoring) to be saved alongside the statistical diagnostics.
#'
#' @param audit An object of class \code{survAudit}.
#' @param assumption A character string matching the name of a non-identifiable 
#'   assumption. Valid options are typically \code{"independent_censoring"}, 
#'   \code{"unmeasured_confounding"}, and \code{"causal_interpretability"}.
#' @param justification A character string containing the justification text.
#'
#' @return The updated \code{survAudit} object.
#' @export
#'
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ age, data = veteran)
#' audit <- survAudit(fit, data = veteran)
#' 
#' audit <- document_assumption(
#'   audit, 
#'   assumption = "independent_censoring", 
#'   justification = "Censoring is administrative at the end of the study."
#' )
#' 
document_assumption <- function(audit, assumption, justification) {
  if (!inherits(audit, "survAudit")) {
    stop("The 'audit' argument must be an object of class 'survAudit'.", call. = FALSE)
  }
  
  if (!is.character(assumption) || length(assumption) != 1L) {
    stop("'assumption' must be a single character string.", call. = FALSE)
  }
  
  if (!is.character(justification) || length(justification) != 1L) {
    stop("'justification' must be a single character string.", call. = FALSE)
  }
  
  non_id_list <- audit$assumptions$non_identifiable
  
  if (is.null(non_id_list) || !(assumption %in% names(non_id_list))) {
    valid_opts <- paste(names(non_id_list), collapse = ", ")
    stop(sprintf(
      "Assumption '%s' not found in the non-identifiable list. Valid options are: %s",
      assumption, valid_opts
    ), call. = FALSE)
  }
  
  # Update the justification
  audit$assumptions$non_identifiable[[assumption]]$justification <- justification
  
  # Return the updated object
  audit
}
