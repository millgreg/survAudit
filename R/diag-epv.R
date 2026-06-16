# Events-per-variable diagnostic for survAudit

#' Compute events-per-variable (EPV) ratio
#'
#' Calculates the events-per-variable ratio for a fitted Cox model, a key
#' metric for assessing whether the sample size is sufficient relative to
#' the model complexity.
#'
#' @param fit A fitted \code{coxph} object.
#'
#' @return A list with components:
#'   \describe{
#'     \item{n_events}{Integer. Number of events in the model.}
#'     \item{n_parameters}{Integer. Number of estimated parameters
#'       (coefficients).}
#'     \item{ratio}{Numeric. The events-per-variable ratio.}
#'     \item{classification}{Character. One of \code{"adequate"} (ratio
#'       >= 20), \code{"marginal"} (10 <= ratio < 20), \code{"low"}
#'       (5 <= ratio < 10), or \code{"unreliable"} (ratio < 5).}
#'   }
#'
#' @keywords internal
#' @noRd
.compute_epv <- function(fit) {
  .validate_coxph(fit)

  n_events <- fit$nevent
  n_parameters <- length(stats::coef(fit))

  ratio <- n_events / n_parameters

  classification <- if (ratio >= 20) {
    "adequate"
  } else if (ratio >= 10) {
    "marginal"
  } else if (ratio >= 5) {
    "low"
  } else {
    "unreliable"
  }

  list(
    n_events = as.integer(n_events),
    n_parameters = as.integer(n_parameters),
    ratio = ratio,
    classification = classification
  )
}
