# Influence diagnostics for survAudit

#' Compute influence diagnostics
#'
#' Calculates dfbeta, dfbetas (standardized), score residuals, and
#' likelihood displacement for each observation. Identifies influential
#' observations using a \code{2 / sqrt(n)} threshold on dfbetas.
#'
#' @param fit A fitted \code{coxph} object.
#'
#' @return A list with components:
#'   \describe{
#'     \item{dfbeta}{Matrix of dfbeta residuals (n x p).}
#'     \item{dfbetas}{Matrix of standardized dfbetas residuals (n x p).}
#'     \item{likelihood_displacement}{Numeric vector of length n. The
#'       likelihood displacement for each observation, computed as
#'       \code{diag(score \%*\% V \%*\% t(score))} where V is the
#'       variance-covariance matrix.}
#'     \item{max_dfbetas}{A list with \code{value} (the maximum
#'       \code{|dfbetas|} value), \code{obs} (the observation index),
#'       and \code{variable} (the variable name).}
#'     \item{max_ld}{A list with \code{value} (the maximum likelihood
#'       displacement) and \code{obs} (the observation index).}
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
  score_residuals <- stats::residuals(fit, type = "score")

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

  if (!is.matrix(score_residuals)) {
    score_residuals <- matrix(score_residuals, ncol = 1)
  }
  if (is.null(colnames(score_residuals))) {
    colnames(score_residuals) <- names(stats::coef(fit))
  }

  # Likelihood displacement: diag(score %*% vcov(fit) %*% t(score))
  # For efficiency, compute row-wise: ld_i = score_i %*% V %*% score_i'
  V <- stats::vcov(fit)
  # Compute score %*% V first (n x p), then element-wise multiply and rowSums
  score_V <- score_residuals %*% V
  likelihood_displacement <- rowSums(score_V * score_residuals)

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

  # Find the maximum likelihood displacement
  max_ld_idx <- which.max(likelihood_displacement)
  max_ld <- list(
    value = likelihood_displacement[max_ld_idx],
    obs = as.integer(max_ld_idx)
  )

  # Flagged observations: any |dfbetas| > threshold
  flagged_rows <- which(apply(abs_dfbetas, 1, function(row) {
    any(row > threshold)
  }))
  flagged_obs <- as.integer(flagged_rows)

  list(
    dfbeta = dfbeta,
    dfbetas = dfbetas,
    likelihood_displacement = likelihood_displacement,
    max_dfbetas = max_dfbetas,
    max_ld = max_ld,
    threshold = threshold,
    flagged_obs = flagged_obs
  )
}
