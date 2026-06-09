#' survAudit: Auditing and Transparency Framework for Cox PH Models
#'
#' Provides a unified diagnostic auditing framework for Cox proportional
#' hazards models. Automates proportional hazards testing, functional form
#' assessment, influence diagnostics, outlier detection, and event sufficiency
#' checks.
#'
#' The package implements a structured assumption ontology that explicitly
#' classifies model assumptions as statistically assessable, partially
#' assessable, or non-identifiable, encouraging transparent reporting of
#' unverifiable assumptions alongside data-driven diagnostics.
#'
#' @section Main function:
#' \code{\link{survAudit}} takes a fitted \code{coxph} object and returns
#' a structured audit object with \code{print}, \code{summary}, and
#' \code{plot} methods.
#'
#' @import ggplot2
#' @import survival
#' @importFrom stats model.frame coef residuals vcov median quantile formula terms qnorm lm predict
#'
#' @docType package
#' @name survAudit-package
#' @aliases survAudit-package
"_PACKAGE"
