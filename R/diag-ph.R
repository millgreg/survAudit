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
#'   }
#'
#' @keywords internal
#' @noRd
.compute_ph_diagnostics <- function(fit, transform = "km") {
  .validate_coxph(fit)

  # Run the cox.zph test
  zph <- survival::cox.zph(fit, transform = transform)

  # Extract the test table
  table_mat <- zph$table

  # Global p-value is in the last row
  global_p <- table_mat[nrow(table_mat), "p"]

  # Extract transformed time
  transformed_time <- zph$x  # numeric vector of transformed event times

  # Extract the un-aggregated scaled Schoenfeld residuals for each dummy coefficient
  # We do this instead of zph$y because newer versions of survival aggregate multi-df terms,
  # which masks opposing time-varying effects within a categorical factor.
  raw_resid <- residuals(fit, type = "scaledsch")
  if (is.null(dim(raw_resid))) {
    raw_resid <- matrix(raw_resid, ncol = 1L)
    colnames(raw_resid) <- names(stats::coef(fit))
  }

  list(
    zph = zph,
    table = table_mat,
    global_p = global_p,
    transform = transform,
    raw_residuals = raw_resid,
    raw_time = transformed_time
  )
}
