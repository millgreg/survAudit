# Internal utility functions for survAudit
# These are not exported and are used by other package functions.

#' Validate that an object is a coxph fit
#'
#' Checks that \code{fit} inherits from class \code{"coxph"} and stops with
#' an informative error message if it does not.
#'
#' @param fit An object to validate.
#'
#' @return Invisible \code{TRUE} if validation passes.
#'
#' @keywords internal
#' @noRd
.validate_coxph <- function(fit) {
  if (!inherits(fit, "coxph")) {
    stop(
      "Expected a 'coxph' object (from survival::coxph), but received an ",
      "object of class: ", paste(class(fit), collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Get continuous covariate names from a Cox model
#'
#' Examines the model frame to identify numeric (non-factor, non-character)
#' variables, excluding the response (Surv) variable.
#'
#' @param fit A fitted \code{coxph} object.
#' @param data Optional data frame. If \code{NULL}, the model frame is used.
#'
#' @return A character vector of continuous covariate names.
#'
#' @keywords internal
#' @noRd
.get_continuous_vars <- function(fit, data = NULL) {
  mf <- stats::model.frame(fit)

  # Identify the response column (Surv object) — always the first column

  response_idx <- attr(stats::terms(mf), "response")
  if (response_idx > 0) {
    mf <- mf[, -response_idx, drop = FALSE]
  }

  is_continuous <- vapply(mf, function(x) {
    is.numeric(x) && !is.factor(x) && length(unique(x)) > 5
  }, logical(1))

  names(which(is_continuous))
}

#' Get categorical covariate names from a Cox model
#'
#' Examines the model frame to identify factor or character variables,
#' excluding the response (Surv) variable.
#'
#' @param fit A fitted \code{coxph} object.
#' @param data Optional data frame. If \code{NULL}, the model frame is used.
#'
#' @return A character vector of categorical covariate names.
#'
#' @keywords internal
#' @noRd
.get_categorical_vars <- function(fit, data = NULL) {
  mf <- stats::model.frame(fit)

  # Remove the response column
  response_idx <- attr(stats::terms(mf), "response")
  if (response_idx > 0) {
    mf <- mf[, -response_idx, drop = FALSE]
  }

  is_categorical <- vapply(mf, function(x) {
    is.factor(x) || is.character(x)
  }, logical(1))

  names(which(is_categorical))
}

#' Format a p-value for display
#'
#' Returns \code{"<0.001"} if the p-value is very small, otherwise rounds
#' to 3 decimal places.
#'
#' @param p A numeric p-value.
#'
#' @return A character string representing the formatted p-value.
#'
#' @keywords internal
#' @noRd
.format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  formatC(round(p, 3), format = "f", digits = 3)
}

#' Create a formatted section header
#'
#' Produces a section header in the form \code{"-- text --------"} padded
#' to the specified width with box-drawing characters.
#'
#' @param text Character string for the header label.
#' @param width Integer total width of the header line.
#'
#' @return A character string.
#'
#' @keywords internal
#' @noRd
.section_header <- function(text, width = 65) {
  prefix <- paste0("\u2500\u2500 ", text, " ")
  remaining <- max(0, width - nchar(prefix))
  tail <- paste(rep("\u2500", remaining), collapse = "")
  paste0(prefix, tail)
}

#' Create a horizontal rule line
#'
#' Produces a horizontal rule using box-drawing characters.
#'
#' @param width Integer width of the rule line.
#'
#' @return A character string.
#'
#' @keywords internal
#' @noRd
.rule_line <- function(width = 65) {
  paste(rep("\u2500", width), collapse = "")
}

#' Extract model metadata from a coxph fit
#'
#' Extracts key information about the fitted Cox model into a named list.
#'
#' @param fit A fitted \code{coxph} object.
#'
#' @return A list with components:
#'   \describe{
#'     \item{call}{The original model call.}
#'     \item{formula}{The model formula.}
#'     \item{coefficients}{Named numeric vector of coefficients.}
#'     \item{n_coef}{Integer count of coefficients.}
#'     \item{method}{Character string: \code{"efron"}, \code{"breslow"},
#'       or \code{"exact"}.}
#'     \item{converged}{Logical indicating whether the model converged.}
#'   }
#'
#' @keywords internal
#' @noRd
.extract_model_info <- function(fit) {
  .validate_coxph(fit)

  coefficients <- stats::coef(fit)


  # Extract the tie-handling method

  method <- fit$method
  if (is.null(method)) {
    method <- "efron"
  }

  # Check convergence
  converged <- TRUE
  if (!is.null(fit$info) && !is.null(fit$info$convergence)) {
    converged <- fit$info$convergence == 0
  }
  # More robust: check if iter reached maxiter without convergence flag

  if (!is.null(fit$fail) && fit$fail) {
    converged <- FALSE
  }

  list(
    call = fit$call,
    formula = fit$formula,
    coefficients = coefficients,
    n_coef = length(coefficients),
    method = method,
    converged = converged
  )
}
