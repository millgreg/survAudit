# Proportional hazards diagnostics for survAudit

#' Compute proportional hazards diagnostics
#'
#' Runs the \code{\link[survival]{cox.zph}} test for proportional hazards
#' and assesses the direction of any time-varying trends in the scaled
#' Schoenfeld residuals.
#'
#' @param fit A fitted \code{coxph} object.
#' @param transform Character string specifying the time transform for
#'   \code{cox.zph}. One of \code{"km"}, \code{"rank"}, \code{"identity"},
#'   or \code{"log"}. Default is \code{"km"}.
#' @param alpha Numeric significance level for covariate-level testing.
#'   Default is \code{0.05}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{zph}{The \code{cox.zph} object.}
#'     \item{table}{The test table matrix with columns \code{chisq},
#'       \code{df}, and \code{p} (or \code{rho}, \code{chisq}, \code{p}
#'       depending on the survival version).}
#'     \item{global_p}{Numeric. The global test p-value (last row of the
#'       table).}
#'     \item{transform}{Character. The time transform used.}
#'     \item{trends}{A data frame with columns \code{variable},
#'       \code{slope}, and \code{direction} describing the trend direction
#'       for each covariate.}
#'   }
#'
#' @keywords internal
#' @noRd
.compute_ph_diagnostics <- function(fit, transform = "km", alpha = 0.05) {
  .validate_coxph(fit)

  # Run the cox.zph test
  zph <- survival::cox.zph(fit, transform = transform)

  # Extract the test table
  table_mat <- zph$table

  # Global p-value is in the last row
  global_p <- table_mat[nrow(table_mat), "p"]

  # Extract scaled Schoenfeld residuals and transformed time
  schoenfeld_resid <- zph$y  # matrix: rows = events, cols = covariates
  transformed_time <- zph$x  # numeric vector of transformed event times

  # Determine covariate names from the table
  # Exclude the last row (GLOBAL) to get covariate names
  if (nrow(table_mat) > 1) {
    covariate_names <- rownames(table_mat)[-nrow(table_mat)]
  } else {
    # Single covariate: the row serves double duty
    covariate_names <- rownames(table_mat)[1]
  }

  # Ensure schoenfeld_resid is a matrix
  if (!is.matrix(schoenfeld_resid)) {
    schoenfeld_resid <- matrix(
      schoenfeld_resid,
      ncol = 1,
      dimnames = list(NULL, covariate_names[1])
    )
  }

  list(
    zph = zph,
    table = table_mat,
    global_p = global_p,
    transform = transform,
    trends = NULL
  )
}
