# Variance Inflation Factors (VIF) diagnostic for survAudit

#' Compute Generalized Variance Inflation Factors (GVIF)
#'
#' Computes the GVIF for a fitted Cox model to assess multicollinearity
#' among the predictors. This implementation uses the GVIF method from
#' the \code{car} package, which correctly groups dummy variables of
#' factor terms to prevent false collinearity flags.
#'
#' @param fit A fitted \code{coxph} object.
#'
#' @return A list with components:
#'   \describe{
#'     \item{vif}{A matrix with columns \code{"GVIF"}, \code{"Df"}, and
#'       \code{"GVIF^(1/(2*Df))"}, or \code{NULL} if not computable.}
#'     \item{flagged}{Character vector of flagged terms.}
#'     \item{message}{Character error message or \code{NULL}.}
#'   }
#'
#' @keywords internal
#' @noRd
.compute_vif <- function(fit) {
  .validate_coxph(fit)

  coefs <- stats::coef(fit)
  if (any(is.na(coefs))) {
    return(list(
      vif = NULL,
      flagged = character(0),
      message = "Aliased coefficients detected in the model."
    ))
  }

  v <- stats::vcov(fit)
  has_intercept <- names(coefs)[1] == "(Intercept)"
  if (has_intercept) {
    v <- v[-1, -1, drop = FALSE]
  }

  x_mat <- stats::model.matrix(fit)
  assign <- attr(x_mat, "assign")
  if (is.null(assign)) {
    return(list(
      vif = NULL,
      flagged = character(0),
      message = "Could not retrieve term assignments from model matrix."
    ))
  }

  if (has_intercept) {
    assign <- assign[-1]
  }

  tms <- labels(stats::terms(fit))
  n_terms <- length(tms)

  if (n_terms < 2) {
    return(list(
      vif = NULL,
      flagged = character(0),
      message = "Model contains fewer than 2 terms."
    ))
  }

  R <- tryCatch({
    stats::cov2cor(v)
  }, error = function(e) {
    NULL
  })

  if (is.null(R)) {
    return(list(
      vif = NULL,
      flagged = character(0),
      message = "Could not compute correlation matrix of coefficients."
    ))
  }

  detR <- tryCatch(det(R), error = function(e) NA)
  if (is.na(detR) || detR <= 0) {
    return(list(
      vif = NULL,
      flagged = character(0),
      message = "Correlation matrix of coefficients is singular or negative definite."
    ))
  }

  result <- matrix(0, n_terms, 3)
  rownames(result) <- tms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")

  for (term_idx in 1:n_terms) {
    subs <- which(assign == term_idx)
    if (length(subs) == 0) next

    det_subs <- tryCatch(det(as.matrix(R[subs, subs, drop = FALSE])), error = function(e) NA)
    det_nsubs <- tryCatch(det(as.matrix(R[-subs, -subs, drop = FALSE])), error = function(e) NA)

    if (is.na(det_subs) || is.na(det_nsubs)) {
      result[term_idx, 1] <- NA
    } else {
      result[term_idx, 1] <- det_subs * det_nsubs / detR
    }
    result[term_idx, 2] <- length(subs)
  }

  for (term_idx in 1:n_terms) {
    gvif <- result[term_idx, 1]
    df <- result[term_idx, 2]
    if (df > 0 && !is.na(gvif)) {
      result[term_idx, 3] <- gvif^(1 / (2 * df))
    } else {
      result[term_idx, 3] <- NA
    }
  }

  flagged <- c()
  for (term_idx in 1:n_terms) {
    df <- result[term_idx, 2]
    gvif <- result[term_idx, 1]
    gvif_scaled <- result[term_idx, 3]

    if (df == 1) {
      if (!is.na(gvif) && gvif > 5) {
        flagged <- c(flagged, tms[term_idx])
      }
    } else {
      if (!is.na(gvif_scaled) && gvif_scaled > 2.236) { # sqrt(5) approx 2.236
        flagged <- c(flagged, tms[term_idx])
      }
    }
  }

  list(
    vif = result,
    flagged = flagged,
    message = NULL
  )
}
