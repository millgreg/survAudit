# ─────────────────────────────────────────────────────────────────
# Tests for survAudit edge cases and graceful degradation
# ─────────────────────────────────────────────────────────────────

library(survival)

veteran <- survival::veteran
fit <- coxph(Surv(time, status) ~ trt + celltype + karno + age, data = veteran)

# ── Test 1: Graceful degradation on engine failures ──────────────
test_that("survAudit gracefully handles engine failures via tryCatch", {
  
  # We can simulate engine failures by passing a severely corrupted model
  # that passes the basic class check but fails during computation.
  fit_broken <- fit
  # Break the variance-covariance matrix (causes influence and VIF to fail)
  fit_broken$var <- matrix("string", nrow = 1, ncol = 1)
  
  # When the engine fails, survAudit should issue a warning and return NULL for that component
  suppressWarnings(
    expect_warning(audit_broken <- survAudit(fit_broken, data = veteran), "VIF computation failed|Influence diagnostics failed|PH diagnostics failed")
  )
  
  expect_s3_class(audit_broken, "survAudit")
  # Influence should be NULL because vcov(fit) throws an error
  expect_null(audit_broken$influence)
  # VIF should be NULL
  expect_null(audit_broken$vif)
})

# ── Test 2: Zero covariate (null) models ─────────────────────────
test_that("survAudit handles zero covariate models gracefully", {
  fit_null <- coxph(Surv(time, status) ~ 1, data = veteran)
  
  # Should not error, and should return NULL for covariate-specific diagnostics
  audit_null <- survAudit(fit_null, data = veteran)
  
  expect_s3_class(audit_null, "survAudit")
  expect_null(audit_null$ph)
  expect_null(audit_null$functional_form)
  expect_null(audit_null$influence)
  expect_null(audit_null$vif)
  
  # But global diagnostics like GOF and Outliers should still exist
  expect_false(is.null(audit_null$outliers))
  expect_false(is.null(audit_null$gof))
  
  # Display methods should not crash
  expect_output(print(audit_null))
  expect_output(print(summary(audit_null)))
  expect_no_error(plot(audit_null, which = "outliers", ask = FALSE))
})

# ── Test 3: Models with only categorical covariates ──────────────
test_that("survAudit handles models with no continuous covariates", {
  fit_cat <- coxph(Surv(time, status) ~ celltype + factor(trt), data = veteran)
  audit_cat <- survAudit(fit_cat, data = veteran)
  
  expect_s3_class(audit_cat, "survAudit")
  # Functional form should be NULL
  expect_null(audit_cat$functional_form)
  
  # Ensure plots don't crash if we ask for functional form when it's NULL
  expect_no_error(plot(audit_cat, which = "functional", ask = FALSE))
})

# ── Test 4: Display methods degrade gracefully ───────────────────
test_that("Display methods degrade gracefully when components are missing", {
  audit_partial <- survAudit(fit, data = veteran)
  
  # Manually strip components to simulate engine failures
  audit_partial$vif <- NULL
  audit_partial$ph <- NULL
  audit_partial$influence <- NULL
  
  expect_output(print(audit_partial))
  
  s <- summary(audit_partial)
  expect_output(print(s))
  
  # Attempting to plot a missing component should throw a message, not an error
  expect_message(plot(audit_partial, which = "ph", ask = FALSE), "PH diagnostics not available")
  expect_message(plot(audit_partial, which = "influence", ask = FALSE), "Influence diagnostics not available")
})
