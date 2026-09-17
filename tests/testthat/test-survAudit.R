# ─────────────────────────────────────────────────────────────────
# Tests for the survAudit package core functionality
# Uses the veteran dataset from the survival package
# ─────────────────────────────────────────────────────────────────

library(survival)

veteran <- survival::veteran
fit <- coxph(Surv(time, status) ~ trt + celltype + karno + age,
             data = veteran)

# ── Test 1: survAudit() returns correct class ────────────────────
test_that("survAudit() returns an object of class 'survAudit'", {
  audit <- survAudit(fit, data = veteran)
  expect_s3_class(audit, "survAudit")
})

# ── Test 2: print() runs without error ───────────────────────────
test_that("print.survAudit() runs without error", {
  audit <- survAudit(fit, data = veteran)
  expect_output(print(audit))
})

# ── Test 3: summary() runs without error ─────────────────────────
test_that("summary.survAudit() runs without error", {
  audit <- survAudit(fit, data = veteran)
  s <- summary(audit)
  expect_s3_class(s, "summary.survAudit")
  expect_output(print(s))
})

# ── Test 4: plot() runs without error for each 'which' value ─────
test_that("plot.survAudit() runs without error for 'ph'", {
  audit <- survAudit(fit, data = veteran)
  p <- plot(audit, which = "ph", ask = FALSE)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plot.survAudit() runs without error for 'influence'", {
  audit <- survAudit(fit, data = veteran)
  p <- plot(audit, which = "influence", ask = FALSE)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plot.survAudit() runs without error for 'outliers'", {
  audit <- survAudit(fit, data = veteran)
  p <- plot(audit, which = "outliers", ask = FALSE)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plot.survAudit() runs without error for 'functional'", {
  audit <- survAudit(fit, data = veteran)
  # functional may be NULL if no continuous vars detected;
  # should not error either way
  expect_no_error(plot(audit, which = "functional", ask = FALSE))
})

test_that("plot.survAudit() runs without error for 'gof'", {
  audit <- survAudit(fit, data = veteran)
  p <- plot(audit, which = "gof", ask = FALSE)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

# ── Test 5: EPV is computed correctly ────────────────────────────
test_that("EPV is computed correctly against manual calculation", {
  audit <- survAudit(fit, data = veteran)
  skip_if(is.null(audit$epv), "EPV diagnostics not available")

  # Manual calculation
  n_events <- sum(veteran$status == 1)
  # celltype is a factor with 4 levels => 3 dummy variables
  # trt, karno, age => 1 parameter each
  # Total: 3 + 1 + 1 + 1 = 6 parameters
  n_params <- length(coef(fit))
  expected_epv <- n_events / n_params

  expect_equal(audit$epv$n_events, n_events)
  expect_equal(audit$epv$n_parameters, n_params)
  expect_equal(audit$epv$ratio, expected_epv, tolerance = 0.01)
})

# ── Test 6: All components of the audit object are present ───────
test_that("survAudit object has all expected components", {
  audit <- survAudit(fit, data = veteran)
  expected_names <- c("model_info", "data_context", "ph",
                      "functional_form", "influence", "outliers",
                      "epv", "vif", "gof", "assumptions", "alpha", "audit_time")
  for (nm in expected_names) {
    expect_true(nm %in% names(audit),
                info = paste("Missing component:", nm))
  }
})

# ── Test 7: Input validation — non-coxph object produces error ───
test_that("survAudit() errors on non-coxph input", {
  expect_error(survAudit(lm(mpg ~ wt, data = mtcars)))
  expect_error(survAudit("not a model"))
  expect_error(survAudit(42))
})

# ── Test 8: survAudit works with a single-covariate model ────────
test_that("survAudit works with a single-covariate model", {
  fit_single <- coxph(Surv(time, status) ~ karno, data = veteran)
  audit <- survAudit(fit_single, data = veteran)
  expect_s3_class(audit, "survAudit")
  expect_output(print(audit))
})

# ── Test 9: survAudit() errors when data argument is omitted ─────
test_that("survAudit() errors when data argument is omitted", {
  expect_error(
    survAudit(fit),
    "The `data` argument must be provided to survAudit()"
  )
})

# ── Test 10: Non-identifiable assumptions have NULL justification ─
test_that("Non-identifiable assumptions have NULL justification by default", {
  audit <- survAudit(fit, data = veteran)
  skip_if(is.null(audit$assumptions),
          "Assumption ontology not available")

  ni <- audit$assumptions$non_identifiable
  skip_if(is.null(ni) || length(ni) == 0L,
          "No non-identifiable assumptions listed")

  for (item in ni) {
    expect_null(item$justification,
                info = paste("Assumption:", item$label %||% item$name))
  }
})

# ── Test 11: alpha parameter is respected ────────────────────────
test_that("alpha parameter is stored correctly", {
  audit_default <- survAudit(fit, data = veteran)
  expect_equal(audit_default$alpha, 0.05)

  audit_strict <- survAudit(fit, data = veteran, alpha = 0.01)
  expect_equal(audit_strict$alpha, 0.01)
})

# ── Test 12: Invalid alpha produces error ────────────────────────
test_that("Invalid alpha values produce errors", {
  expect_error(survAudit(fit, data = veteran, alpha = 0))
  expect_error(survAudit(fit, data = veteran, alpha = 1))
  expect_error(survAudit(fit, data = veteran, alpha = -0.5))
  expect_error(survAudit(fit, data = veteran, alpha = "abc"))
})

# ── Test 13: VIF is calculated correctly ─────────────────────────
test_that("VIF is calculated correctly and handles factor terms", {
  audit <- survAudit(fit, data = veteran)
  expect_false(is.null(audit$vif))
  expect_true(is.matrix(audit$vif$vif))
  expect_equal(colnames(audit$vif$vif), c("GVIF", "Df", "GVIF^(1/(2*Df))"))
  # trt, celltype, karno, age are 4 terms
  expect_equal(nrow(audit$vif$vif), 4)
  expect_true(all(audit$vif$vif[, "Df"] >= 1))
})

# ── Test 14: survAudit works with models containing interaction terms ───
test_that("survAudit works with interaction terms", {
  fit_int <- coxph(Surv(time, status) ~ trt * age + karno, data = veteran)
  audit_int <- survAudit(fit_int, data = veteran)
  expect_s3_class(audit_int, "survAudit")
  # VIF should gracefully handle or flag interaction terms
  expect_output(print(audit_int))
})

# ── Test 15: survAudit handles missing data (survival::lung dataset) ────
test_that("survAudit handles missing data properly", {
  # lung dataset has missing values in covariates like meal.cal, wt.loss
  lung_data <- survival::lung
  fit_miss <- coxph(Surv(time, status) ~ age + sex + ph.ecog + meal.cal + wt.loss, 
                    data = lung_data)
  
  audit_miss <- survAudit(fit_miss, data = lung_data)
  expect_s3_class(audit_miss, "survAudit")
  
  # Data context should document missingness
  expect_true(is.data.frame(audit_miss$data_context$missing_data))
  expect_true(sum(audit_miss$data_context$missing_data$n_missing) > 0)
  
  # Ensure plots run without error despite missingness
  expect_no_error(plot(audit_miss, which = "ph", ask = FALSE))
  expect_no_error(plot(audit_miss, which = "functional", ask = FALSE))
})

# ── Test 16: Graceful degradation on engine failures ──────────────
test_that("survAudit gracefully handles engine failures via tryCatch", {
  fit_broken <- fit
  fit_broken$var <- matrix("string", nrow = 1, ncol = 1)
  
  suppressWarnings(
    expect_warning(audit_broken <- survAudit(fit_broken, data = veteran), "VIF computation failed|Influence diagnostics failed|PH diagnostics failed")
  )
  
  expect_s3_class(audit_broken, "survAudit")
  expect_null(audit_broken$influence)
  expect_null(audit_broken$vif)
})

# ── Test 17: Zero covariate (null) models ─────────────────────────
test_that("survAudit handles zero covariate models gracefully", {
  fit_null <- coxph(Surv(time, status) ~ 1, data = veteran)
  audit_null <- survAudit(fit_null, data = veteran)
  
  expect_s3_class(audit_null, "survAudit")
  expect_null(audit_null$ph)
  expect_null(audit_null$functional_form)
  expect_null(audit_null$influence)
  expect_null(audit_null$vif)
  
  expect_false(is.null(audit_null$outliers))
  expect_false(is.null(audit_null$gof))
  
  expect_output(print(audit_null))
  expect_output(print(summary(audit_null)))
  expect_no_error(plot(audit_null, which = "outliers", ask = FALSE))
})

# ── Test 18: Models with only categorical covariates ──────────────
test_that("survAudit handles models with no continuous covariates", {
  fit_cat <- coxph(Surv(time, status) ~ celltype + factor(trt), data = veteran)
  audit_cat <- survAudit(fit_cat, data = veteran)
  
  expect_s3_class(audit_cat, "survAudit")
  expect_null(audit_cat$functional_form)
  expect_no_error(plot(audit_cat, which = "functional", ask = FALSE))
})

# ── Test 19: Display methods degrade gracefully ───────────────────
test_that("Display methods degrade gracefully when components are missing", {
  audit_partial <- survAudit(fit, data = veteran)
  
  audit_partial$vif <- NULL
  audit_partial$ph <- NULL
  audit_partial$influence <- NULL
  
  expect_output(print(audit_partial))
  expect_output(print(summary(audit_partial)))
  
  expect_message(plot(audit_partial, which = "ph", ask = FALSE), "PH diagnostics not available")
  expect_message(plot(audit_partial, which = "influence", ask = FALSE), "Influence diagnostics not available")
})

# ── Test 20: plot() with vars argument ─────────────────────────────
test_that("plot() supports vars argument for targeted covariate plotting", {
  audit_test <- survAudit(fit, data = veteran)

  # Functional form with single var
  p_ff_single <- plot(audit_test, which = "functional", vars = "age", ask = FALSE)
  expect_s3_class(p_ff_single, "ggplot")

  # PH with subset of vars
  p_ph_sub <- plot(audit_test, which = "ph", vars = c("age", "karno"), ask = FALSE)
  expect_s3_class(p_ph_sub, "ggplot")

  # Influence with subset of vars
  p_inf_sub <- plot(audit_test, which = "influence", vars = "age", ask = FALSE)
  expect_s3_class(p_inf_sub, "ggplot")

  # Non-existent variable warns and handles gracefully
  expect_warning(
    plot(audit_test, which = "functional", vars = "nonexistent", ask = FALSE),
    "None of the requested covariates"
  )
})

