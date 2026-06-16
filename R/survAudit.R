#' Audit a Cox Proportional Hazards Model
#'
#' Performs a comprehensive diagnostic audit of a fitted Cox proportional
#' hazards model, including proportional hazards testing, functional form
#' assessment, influence diagnostics, outlier detection, and event
#' sufficiency checks.
#'
#' @param fit A fitted \code{\link[survival]{coxph}} object from the
#'   \code{survival} package.
#' @param data Optional data frame used to fit the model. If \code{NULL}
#'   (the default), the function attempts to extract it from
#'   \code{fit$model} or via \code{\link[stats]{model.frame}}.
#' @param alpha Significance level for flagging diagnostics (default
#'   \code{0.05}). Used for display and interpretation only; does not
#'   affect the underlying tests.
#' @param ph_transform Transform for the \code{\link[survival]{cox.zph}}
#'   test. One of \code{"km"} (default), \code{"rank"},
#'   \code{"identity"}, or \code{"log"}.
#'
#' @return An object of class \code{survAudit} containing:
#'   \describe{
#'     \item{model_info}{List of model metadata (call, formula,
#'       coefficients, etc.).}
#'     \item{data_context}{List summarising the data (sample size, events,
#'       censoring, time distribution, ties, missing data).}
#'     \item{ph}{Proportional hazards diagnostic results or \code{NULL}.}
#'     \item{functional_form}{Functional form assessment results or
#'       \code{NULL}.}
#'     \item{influence}{Influence diagnostic results or \code{NULL}.}
#'     \item{outliers}{Outlier assessment results or \code{NULL}.}
#'     \item{epv}{Events-per-variable ratio and classification or
#'       \code{NULL}.}
#'     \item{gof}{Overall model goodness-of-fit results or \code{NULL}.}
#'     \item{assumptions}{Structured assumption classification.}
#'     \item{alpha}{The significance level used.}
#'     \item{audit_time}{Timestamp of the audit.}
#'   }
#'
#' @export
#'
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ trt + celltype + karno + age,
#'              data = veteran)
#' audit <- survAudit(fit)
#' print(audit)
survAudit <- function(fit, data = NULL, alpha = 0.05, ph_transform = "km") {

 # ── 1. Validate input ─────────────────────────────────────────────
 .validate_coxph(fit)

 ph_transform <- match.arg(ph_transform,
                            choices = c("km", "rank", "identity", "log"))

 if (!is.numeric(alpha) || length(alpha) != 1L ||
     alpha <= 0 || alpha >= 1) {
   stop("`alpha` must be a single numeric value in (0, 1).", call. = FALSE)
 }

 audit_time <- Sys.time()

 # ── 2. Extract data from model if not provided ────────────────────
 data_available <- TRUE
 if (is.null(data)) {
   data <- tryCatch(
     {
       d <- fit$model
       if (is.null(d)) d <- model.frame(fit)
       d
     },
     error = function(e) NULL
   )
   if (is.null(data)) {
     warning(
       "Could not extract data from model object. ",
       "Some diagnostics (functional form) may be unavailable. ",
       "Provide the data explicitly via the `data` argument.",
       call. = FALSE
     )
     data_available <- FALSE
   }
 }

 # ── 3. Extract model info ─────────────────────────────────────────
 model_info <- .extract_model_info(fit)

 # ── 4. Compute data context ───────────────────────────────────────
 data_context <- tryCatch(
   .compute_data_context(fit, data),
   error = function(e) {
     warning("Data context computation failed: ", conditionMessage(e),
             call. = FALSE)
     NULL
   }
 )

 # ── 5. Identify continuous variables ──────────────────────────────
 continuous_vars <- if (data_available) {
   tryCatch(
     .get_continuous_vars(fit, data),
     error = function(e) {
       warning("Could not identify continuous variables: ",
               conditionMessage(e), call. = FALSE)
       character(0)
     }
   )
 } else {
   character(0)
 }

 # ── 6. Run all diagnostic engines ─────────────────────────────────
 ph <- tryCatch(
   .compute_ph_diagnostics(fit, transform = ph_transform, alpha = alpha),
   error = function(e) {
     warning("PH diagnostics failed: ", conditionMessage(e),
             call. = FALSE)
     NULL
   }
 )

 functional_form <- tryCatch(
   {
     if (!data_available || length(continuous_vars) == 0L) {
       NULL
     } else {
       .compute_functional_form(fit, data, continuous_vars)
     }
   },
   error = function(e) {
     warning("Functional form diagnostics failed: ",
             conditionMessage(e), call. = FALSE)
     NULL
   }
 )

 influence <- tryCatch(
   .compute_influence(fit),
   error = function(e) {
     warning("Influence diagnostics failed: ", conditionMessage(e),
             call. = FALSE)
     NULL
   }
 )

 outliers <- tryCatch(
   .compute_outliers(fit),
   error = function(e) {
     warning("Outlier diagnostics failed: ", conditionMessage(e),
             call. = FALSE)
     NULL
   }
 )

 epv <- tryCatch(
   .compute_epv(fit),
   error = function(e) {
     warning("EPV computation failed: ", conditionMessage(e),
             call. = FALSE)
     NULL
   }
 )

 vif <- tryCatch(
   .compute_vif(fit),
   error = function(e) {
     warning("VIF computation failed: ", conditionMessage(e),
             call. = FALSE)
     NULL
   }
 )

 gof <- tryCatch(
   .compute_gof(fit),
   error = function(e) {
     warning("Goodness-of-Fit computation failed: ", conditionMessage(e),
             call. = FALSE)
     NULL
   }
 )

 # ── 7. Build assumption classification ────────────────────────────
 assumptions <- tryCatch(
   .build_assumptions(ph, functional_form, influence, outliers, epv,
                      data_context, alpha, vif),
   error = function(e) {
     warning("Assumption classification construction failed: ",
             conditionMessage(e), call. = FALSE)
     list(assessable = list(), partially_assessable = list(),
          non_identifiable = list())
   }
 )

 # ── 8. Construct and return survAudit object ──────────────────────
 obj <- .new_survAudit(
   model_info    = model_info,
   data_context  = data_context,
   ph            = ph,
   functional_form = functional_form,
   influence     = influence,
   outliers      = outliers,
   epv           = epv,
   vif           = vif,
   gof           = gof,
   assumptions   = assumptions,
   alpha         = alpha,
   audit_time    = audit_time
 )

 .validate_survAudit(obj)
 obj
}
