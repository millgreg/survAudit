# Functional form diagnostics for survAudit

#' Compute functional form diagnostics for continuous covariates
#'
#' For each continuous covariate, fits a reduced Cox model excluding that
#' covariate, extracts martingale residuals, and compares a loess smooth
#' to a linear fit to detect departures from linearity.
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
#'           reduced model.}
#'         \item{covariate_values}{Numeric vector of covariate values.}
#'         \item{loess_fit}{The \code{loess} object.}
#'         \item{departure_detected}{Logical indicating whether the loess
#'           smooth captures substantially more variance than the linear
#'           fit.}
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

  original_formula <- stats::formula(fit)

  # Extract the complete cases used in the original model
  mf <- stats::model.frame(fit)
  data_complete <- data[rownames(mf), , drop = FALSE]

  for (var in continuous_vars) {
    result <- tryCatch({
      # Build a reduced formula excluding this covariate
      reduced_formula <- .remove_term(original_formula, var)

      # Fit the reduced model (could be a null model: Surv(time, status) ~ 1)
      null_fit <- survival::coxph(reduced_formula, data = data_complete)

      # Extract martingale residuals from the null model
      mart_resid <- stats::residuals(null_fit, type = "martingale")

      # Get covariate values from the complete cases data
      covariate_values <- data_complete[[var]]

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

      # Fit linear model
      lm_fit <- stats::lm(mart_resid_clean ~ covariate_values_clean)

      # Compare: residual variance ratio
      loess_resid_var <- mean((mart_resid_clean - stats::fitted(loess_fit))^2)
      lm_resid_var <- mean((mart_resid_clean - stats::fitted(lm_fit))^2)

      # If loess captures substantially more variance, flag departure
      # ratio > 1.1 means linear residual variance is >10% larger than
      # loess residual variance
      if (loess_resid_var > 0) {
        variance_ratio <- lm_resid_var / loess_resid_var
      } else {
        variance_ratio <- 1.0
      }
      departure_detected <- variance_ratio > 1.1

      list(
        residuals = mart_resid_clean,
        covariate_values = covariate_values_clean,
        loess_fit = loess_fit,
        departure_detected = departure_detected
      )
    }, error = function(e) {
      # If the null model can't be fitted, skip gracefully
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


#' Remove a term from a formula
#'
#' Removes a single term from a formula by name, returning the modified
#' formula.
#'
#' @param formula A formula object.
#' @param term Character string naming the term to remove.
#'
#' @return A modified formula with the term removed.
#'
#' @keywords internal
#' @noRd
.remove_term <- function(formula, term) {
  # Use update to remove the term
  drop_formula <- stats::as.formula(paste0("~ . - ", term))
  stats::update(formula, drop_formula)
}
