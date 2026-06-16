# S3 class definition for survAudit objects

#' Low-level constructor for survAudit objects
#'
#' Creates a new \code{survAudit} object from pre-computed diagnostic
#' components. This is the internal constructor; users should call
#' \code{\link{survAudit}} instead.
#'
#' @param model_info List of model metadata from \code{.extract_model_info}.
#' @param data_context List of data context from \code{.compute_data_context}.
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
#' @param assumptions List of structured assumptions from
#'   \code{.build_assumptions}.
#' @param alpha Numeric significance level used for testing.
#' @param audit_time POSIXct timestamp of when the audit was performed.
#'
#' @return A list with class \code{"survAudit"}.
#'
#' @keywords internal
#' @noRd
.new_survAudit <- function(model_info,
                           data_context,
                           ph,
                           functional_form,
                           influence,
                           outliers,
                           epv,
                           vif,
                           gof,
                           assumptions,
                           alpha,
                           audit_time) {
  structure(
    list(
      model_info = model_info,
      data_context = data_context,
      ph = ph,
      functional_form = functional_form,
      influence = influence,
      outliers = outliers,
      epv = epv,
      vif = vif,
      gof = gof,
      assumptions = assumptions,
      alpha = alpha,
      audit_time = audit_time
    ),
    class = "survAudit"
  )
}

#' Validate a survAudit object
#'
#' Checks that all required components are present in a \code{survAudit}
#' object and that they have the correct types.
#'
#' @param x An object to validate as a \code{survAudit}.
#'
#' @return Invisible \code{TRUE} if validation passes. Stops with an error
#'   if any required component is missing or has an incorrect type.
#'
#' @keywords internal
#' @noRd
.validate_survAudit <- function(x) {
  if (!inherits(x, "survAudit")) {
    stop("Object does not have class 'survAudit'.", call. = FALSE)
  }

  required_components <- c(
    "model_info", "data_context", "ph", "functional_form",
    "influence", "outliers", "epv", "vif", "gof", "assumptions", "alpha", "audit_time"
  )

  missing <- setdiff(required_components, names(x))
  if (length(missing) > 0) {
    stop(
      "survAudit object is missing required components: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }

  # Type checks — required components must be lists
  required_lists <- c("model_info", "assumptions")
  for (comp in required_lists) {
    if (!is.list(x[[comp]])) {
      stop(
        "Component '", comp, "' must be a list, but is ",
        class(x[[comp]])[1], ".",
        call. = FALSE
      )
    }
  }

  # Optional diagnostic components: must be list or NULL
  optional_lists <- c("data_context", "ph", "functional_form",
                       "influence", "outliers", "epv", "vif", "gof")
  for (comp in optional_lists) {
    if (!is.null(x[[comp]]) && !is.list(x[[comp]])) {
      stop(
        "Component '", comp, "' must be a list or NULL, but is ",
        class(x[[comp]])[1], ".",
        call. = FALSE
      )
    }
  }

  if (!is.numeric(x$alpha) || length(x$alpha) != 1) {
    stop("Component 'alpha' must be a single numeric value.", call. = FALSE)
  }

  if (!inherits(x$audit_time, "POSIXct")) {
    stop("Component 'audit_time' must be a POSIXct object.", call. = FALSE)
  }

  invisible(TRUE)
}
