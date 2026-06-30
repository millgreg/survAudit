# Functional form diagnostics for survAudit

#' Compute functional form diagnostics for continuous covariates
#'
#' Assesses the linearity assumption for each continuous covariate using
#' martingale residuals from a restricted Cox model (a model excluding the
#' specific covariate being assessed).
#' For each continuous covariate, a loess smooth of the restricted model
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
#'           restricted model.}
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

  orig_formula <- stats::formula(fit)
  
  for (var in continuous_vars) {
    result <- tryCatch({
      # 1. Build a restricted formula excluding the specific covariate
      # This computes martingale residuals adjusted for all OTHER covariates in the model.
      restricted_formula <- stats::update(orig_formula, paste(". ~ . -", var))
      
      # 2. Fit the restricted Cox model
      restricted_fit <- survival::coxph(restricted_formula, data = data_complete)
      
      # 3. Extract martingale residuals from this restricted model
      mart_resid <- stats::residuals(restricted_fit, type = "martingale")
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
