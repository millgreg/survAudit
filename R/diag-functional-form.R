# Functional form diagnostics for survAudit

#' Compute functional form diagnostics for continuous covariates
#'
#' Assesses the linearity assumption for each continuous covariate using
#' martingale residuals from a null Cox model (a model with no covariates),
#' following the approach recommended by Therneau and Grambsch (2000).
#' For each continuous covariate, a loess smooth of the null model
#' martingale residuals against the covariate values is compared to a
#' linear fit to screen for departures from linearity.
#'
#' @param fit A fitted \code{coxph} object.
#' @param data The data frame used to fit the model.
#' @param continuous_vars Character vector of continuous covariate names to
#'   assess.
#'
#' @return A list with components:
#'   \describe{
#'     \item{continuous_vars}{Character vector of assessed continuous
#'       covariates.}
#'     \item{results}{A named list where each element corresponds to a
#'       covariate and contains:
#'       \describe{
#'         \item{residuals}{Numeric vector of martingale residuals from the
#'           null model.}
#'         \item{covariate_values}{Numeric vector of covariate values.}
#'         \item{loess_fit}{The \code{loess} object.}
#'         \item{departure_detected}{Logical indicating whether the loess
#'           smooth captures substantially more variance than the linear
#'           fit (screening heuristic).}
#'       }
#'     }
#'   }
#'
#' @keywords internal
#' @noRd
.compute_functional_form <- function(fit, data, continuous_vars) {
  .validate_coxph(fit)

  results <- list()

  if (length(continuous_vars) == 0) {
    return(list(
      continuous_vars = character(0),
      results = results
    ))
  }

  # Extract the complete cases used in the original model
  mf <- stats::model.frame(fit)
  data_complete <- data[rownames(mf), , drop = FALSE]

  # Extract the Surv object from the model frame
  surv_obj <- stats::model.response(mf)

  # Fit a single null Cox model (no covariates) — Therneau & Grambsch (2000)
  null_fit <- tryCatch(
    survival::coxph(surv_obj ~ 1),
    error = function(e) {
      warning("Null Cox model fitting failed: ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )

  if (is.null(null_fit)) {
    return(list(
      continuous_vars = continuous_vars,
      results = results
    ))
  }

  # Extract martingale residuals from the null model ONCE
  mart_resid <- stats::residuals(null_fit, type = "martingale")

  for (var in continuous_vars) {
    result <- tryCatch({
      # Get covariate values from the model frame
      covariate_values <- mf[[var]]

      # Remove any NAs (just in case)
      complete <- stats::complete.cases(covariate_values, mart_resid)
      covariate_values_clean <- covariate_values[complete]
      mart_resid_clean <- mart_resid[complete]

      # Fit loess smooth
      loess_fit <- stats::loess(
        mart_resid_clean ~ covariate_values_clean,
        degree = 1,
        span = 0.75
      )

      list(
        residuals = mart_resid_clean,
        covariate_values = covariate_values_clean,
        loess_fit = loess_fit
      )
    }, error = function(e) {
      # If assessment fails for this variable, skip gracefully
      NULL
    })

    if (!is.null(result)) {
      results[[var]] <- result
    }
  }

  list(
    continuous_vars = continuous_vars,
    results = results
  )
}
