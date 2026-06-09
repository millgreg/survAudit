# Calibration diagnostics for survAudit

#' Compute overall model calibration using Cox-Snell residuals
#'
#' Calculates Cox-Snell residuals as (Event Status - Martingale Residual).
#' Estimates the Nelson-Aalen cumulative hazard of these residuals to
#' support the calibration diagnostic plot.
#'
#' @param fit A fitted \code{coxph} object.
#'
#' @return A list with components:
#'   \describe{
#'     \item{coxsnell}{Numeric vector of Cox-Snell residuals.}
#'     \item{plot_x}{Numeric vector of unique residual values (for plotting).}
#'     \item{plot_y}{Numeric vector of the Nelson-Aalen cumulative hazard (for plotting).}
#'   }
#'
#' @keywords internal
#' @noRd
.compute_calibration <- function(fit) {
  .validate_coxph(fit)

  mf <- stats::model.frame(fit)
  surv_obj <- stats::model.response(mf)
  
  # Extract status (1 = event, 0 = censored). The status is the last column of Surv object
  status <- surv_obj[, ncol(surv_obj)] 
  
  mart_resid <- stats::residuals(fit, type = "martingale")
  
  coxsnell <- status - mart_resid
  
  # Compute empirical cumulative hazard of Cox-Snell residuals
  plot_data <- tryCatch({
    cs_surv <- survival::Surv(time = coxsnell, event = status)
    cs_fit <- survival::survfit(cs_surv ~ 1)
    
    # Return time (x-axis) and cumhaz (y-axis)
    list(
      x = cs_fit$time,
      y = cs_fit$cumhaz
    )
  }, error = function(e) {
    list(x = numeric(0), y = numeric(0))
  })

  list(
    coxsnell = coxsnell,
    plot_x = plot_data$x,
    plot_y = plot_data$y
  )
}
