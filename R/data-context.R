# Data context computation for survAudit

#' Compute comprehensive data context from a fitted coxph model
#'
#' Extracts sample sizes, event counts, censoring information, time
#' distributions, tie statistics, and missing data summaries from a
#' fitted Cox proportional hazards model.
#'
#' @param fit A fitted \code{coxph} object.
#' @param data Optional data frame used to fit the model. If provided, used
#'   to check for missing values in the original data. If \code{NULL}, missing
#'   data information will be reported as \code{"not available"}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{n}{Integer. Total sample size.}
#'     \item{n_events}{Integer. Number of events.}
#'     \item{n_censored}{Integer. Number of censored observations.}
#'     \item{censoring_rate}{Numeric. Proportion of censored observations.}
#'     \item{time_range}{Numeric vector of length 2: \code{c(min, max)}.}
#'     \item{time_median}{Numeric. Median of event/censoring times.}
#'     \item{time_iqr}{Numeric vector of length 2: \code{c(Q1, Q3)}.}
#'     \item{n_ties}{Integer. Number of tied event times.}
#'     \item{tie_fraction}{Numeric. Proportion of events with tied times.}
#'     \item{missing_data}{A data frame or character string. If \code{data} is
#'       provided, a data frame with columns \code{variable}, \code{n_missing},
#'       and \code{pct_missing}. Otherwise the string \code{"not available"}.}
#'     \item{counting_process}{Logical. \code{TRUE} if the model uses
#'       start-stop (counting process) format.}
#'   }
#'
#' @keywords internal
#' @noRd
.compute_data_context <- function(fit, data = NULL) {
  .validate_coxph(fit)

  # Extract the Surv response object
  surv_obj <- fit$y
  if (is.null(surv_obj)) {
    mf <- stats::model.frame(fit)
    surv_obj <- stats::model.response(mf)
  }

  # Detect counting process format (start-stop has 3 columns)
  counting_process <- ncol(surv_obj) == 3

  # Extract times and status
  if (counting_process) {
    # For counting process: columns are (start, stop, status)
    times <- surv_obj[, 2]
    status <- surv_obj[, 3]
  } else {
    # Standard format: columns are (time, status)
    times <- surv_obj[, 1]
    status <- surv_obj[, 2]
  }

  n <- length(status)
  n_events <- sum(status == 1)
  n_censored <- n - n_events

  # Censoring rate
  censoring_rate <- n_censored / n

  # Time distribution
  time_range <- c(min(times), max(times))
  if (counting_process) {
    # Median and IQR are meaningless for fragmented counting process rows
    time_median <- NA_real_
    time_iqr <- c(NA_real_, NA_real_)
  } else {
    time_median <- stats::median(times)
    time_iqr <- as.numeric(stats::quantile(times, probs = c(0.25, 0.75)))
  }

  # Tied event times
  event_times <- times[status == 1]
  if (length(event_times) > 0) {
    time_table <- table(event_times)
    n_tied_groups <- sum(time_table > 1)
    n_tied_events <- sum(time_table[time_table > 1])
    n_ties <- n_tied_events
    tie_fraction <- n_ties / n_events
  } else {
    n_ties <- 0L
    tie_fraction <- 0
  }

  # Missing data
  if (!is.null(data)) {
    # Get covariate names from the model terms
    model_terms <- stats::terms(fit)
    term_labels <- attr(model_terms, "term.labels")

    # Identify all base variable names in the model (both LHS and RHS)
    # This ensures we flag if patients were dropped due to missing survival time or status
    all_vars <- all.vars(stats::formula(fit))

    # Only check variables that exist in the provided data
    available_vars <- intersect(all_vars, names(data))

    if (length(available_vars) > 0) {
      n_missing <- vapply(available_vars, function(v) {
        sum(is.na(data[[v]]))
      }, integer(1))

      pct_missing <- round(n_missing / nrow(data) * 100, 2)

      missing_data <- data.frame(
        variable = available_vars,
        n_missing = n_missing,
        pct_missing = pct_missing,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    } else {
      missing_data <- data.frame(
        variable = character(0),
        n_missing = integer(0),
        pct_missing = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  } else {
    missing_data <- "not available"
  }

  list(
    n = n,
    n_events = n_events,
    n_censored = n_censored,
    censoring_rate = censoring_rate,
    time_range = time_range,
    time_median = time_median,
    time_iqr = time_iqr,
    n_ties = as.integer(n_ties),
    tie_fraction = tie_fraction,
    missing_data = missing_data,
    counting_process = counting_process
  )
}
