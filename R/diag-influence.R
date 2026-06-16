# Influence diagnostics for survAudit

#' Compute influence diagnostics
#'
#' Calculates dfbeta and dfbetas (standardized) residuals for each
#' observation. Identifies influential observations using a
#' \code{2 / sqrt(n)} threshold on dfbetas.
#'
#' @param fit A fitted \code{coxph} object.
#'
#' @return A list with components:
#'   \describe{
#'     \item{dfbeta}{Matrix of dfbeta residuals (n x p).}
#'     \item{dfbetas}{Matrix of standardized dfbetas residuals (n x p).}
#'     \item{max_dfbetas}{A list with \code{value} (the maximum
#'       \code{|dfbetas|} value), \code{obs} (the observation index),
#'       and \code{variable} (the variable name).}
#'     \item{threshold}{Numeric. The \code{2 / sqrt(n)} threshold used
#'       for flagging.}
#'     \item{flagged_obs}{Integer vector of observation indices where any
#'       \code{|dfbetas|} exceeds the threshold.}
#'   }
#'
#' @keywords internal
#' @noRd
.compute_influence <- function(fit) {
  .validate_coxph(fit)

  # Extract residual matrices
  dfbeta <- stats::residuals(fit, type = "dfbeta")
  dfbetas <- stats::residuals(fit, type = "dfbetas")

  # Ensure matrices and assign column names if NULL
  if (!is.matrix(dfbeta)) {
    dfbeta <- matrix(dfbeta, ncol = 1)
  }
  if (is.null(colnames(dfbeta))) {
    colnames(dfbeta) <- names(stats::coef(fit))
  }

  if (!is.matrix(dfbetas)) {
    dfbetas <- matrix(dfbetas, ncol = 1)
  }
  if (is.null(colnames(dfbetas))) {
    colnames(dfbetas) <- names(stats::coef(fit))
  }

  # Number of observations
  nobs <- nrow(dfbetas)

  # Threshold for dfbetas flagging
  threshold <- 2 / sqrt(nobs)

  # Find the maximum |dfbetas|
  abs_dfbetas <- abs(dfbetas)
  max_idx <- which.max(abs_dfbetas)  # linear index
  max_row <- ((max_idx - 1) %% nrow(abs_dfbetas)) + 1
  max_col <- ((max_idx - 1) %/% nrow(abs_dfbetas)) + 1

  max_dfbetas <- list(
    value = abs_dfbetas[max_row, max_col],
    obs = as.integer(max_row),
    variable = colnames(dfbetas)[max_col]
  )

  # Flagged observations: any |dfbetas| > threshold
  flagged_rows <- which(apply(abs_dfbetas, 1, function(row) {
    any(row > threshold)
  }))
  flagged_obs <- as.integer(flagged_rows)

  list(
    dfbeta = dfbeta,
    dfbetas = dfbetas,
    max_dfbetas = max_dfbetas,
    threshold = threshold,
    flagged_obs = flagged_obs
  )
}
