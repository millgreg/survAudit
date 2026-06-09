# Outlier detection diagnostics for survAudit

#' Compute outlier diagnostics
#'
#' Detects potential outliers using multiple residual types: martingale,
#' deviance, log-odds, and normal deviate residuals. Observations are
#' flagged when residuals exceed standard thresholds.
#'
#' @param fit A fitted \code{coxph} object.
#'
#' @return A list with components:
#'   \describe{
#'     \item{martingale}{Numeric vector of martingale residuals.}
#'     \item{deviance}{Numeric vector of deviance residuals.}
#'     \item{log_odds}{Numeric vector of log-odds residuals
#'       (\code{log(S / (1 - S))} where S is predicted survival).}
#'     \item{normal_deviate}{Numeric vector of normal deviate residuals
#'       (\code{qnorm(S)}).}
#'     \item{linear_predictor}{Numeric vector of linear predictor values
#'       (for use as the x-axis in diagnostic plots).}
#'     \item{flagged}{A list with integer vectors of observation indices:
#'       \describe{
#'         \item{deviance}{Indices where \code{|deviance| > 1.96}.}
#'         \item{log_odds}{Indices where \code{|log_odds| > 3.66}.}
#'         \item{normal_deviate}{Indices where \code{|normal_deviate| > 1.96}.}
#'       }
#'     }
#'   }
#'
#' @keywords internal
#' @noRd
.compute_outliers <- function(fit) {
  .validate_coxph(fit)

  # Extract residuals
  martingale <- stats::residuals(fit, type = "martingale")
  deviance <- stats::residuals(fit, type = "deviance")

  # Predicted survival: S = exp(-H) where H = expected cumulative hazard
  expected_cumhaz <- stats::predict(fit, type = "expected")
  predicted_survival <- exp(-expected_cumhaz)

  # Linear predictor for x-axis in plots
  linear_predictor <- stats::predict(fit, type = "lp")

  # Clamp predicted survival away from 0 and 1 to avoid Inf/-Inf
  eps <- .Machine$double.eps^0.5
  S_clamped <- pmin(pmax(predicted_survival, eps), 1 - eps)

  # Log-odds: log(S / (1 - S))
  log_odds <- log(S_clamped / (1 - S_clamped))

  # Normal deviate: qnorm(S)
  normal_deviate <- stats::qnorm(S_clamped)

  # Flag observations exceeding thresholds
  flagged_deviance <- as.integer(which(abs(deviance) > 1.96))
  flagged_log_odds <- as.integer(which(abs(log_odds) > 3.66))
  flagged_normal_deviate <- as.integer(which(abs(normal_deviate) > 1.96))

  list(
    martingale = martingale,
    deviance = deviance,
    log_odds = log_odds,
    normal_deviate = normal_deviate,
    linear_predictor = linear_predictor,
    flagged = list(
      deviance = flagged_deviance,
      log_odds = flagged_log_odds,
      normal_deviate = flagged_normal_deviate
    )
  )
}
