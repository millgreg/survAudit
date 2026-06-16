# survAudit

[![R-CMD-check](https://github.com/millgreg/survAudit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/millgreg/survAudit/actions/workflows/R-CMD-check.yaml)

> An Auditing and Transparency Framework for Cox Proportional Hazards Models

`survAudit` provides a unified diagnostic auditing framework for Cox proportional hazards (PH) models. Rather than relying on fragmented diagnostic plots and tests, `survAudit` aggregates critical statistical diagnostics into a structured **Assumption Classification** and provides a mechanism to document qualitative justifications for non-testable assumptions directly on the R model audit object.

Full documentation and a reproducible case study are available at: **[https://millgreg.github.io/survAudit](https://millgreg.github.io/survAudit)**

## Features

1. **Unified Assumption Classification**: Explicitly organizes model assumptions into:
   * **Non-Identifiable**: Assumptions that cannot be verified statistically (e.g., Independent Censoring, Absence of Unmeasured Confounding) and require qualitative justification.
   * **Partially Assessable**: Assumptions informed by metrics but requiring clinical/domain judgment (e.g., Outlier Impact, Missing Data).
   * **Statistically Assessable**: Assumptions rigorously testable from data (e.g., Proportional Hazards, Functional Form, Influence Stability, Event Sufficiency).
2. **Advanced Diagnostics**:
   * **Collinearity (GVIF)**: Implements the Generalized Variance Inflation Factor (GVIF) algorithm (matching `car` package behavior) to correctly group dummy variables of factor predictors and prevent false collinearity flags.
   * **Functional Form (Linearity)**: Martingale residuals from multivariate reduced models are compared against continuous predictors using LOESS smooths with 95% confidence bands to identify departures from linearity visually.
   * **Visual Diagnostics**: Overall goodness-of-fit via Cox-Snell residuals, Outlier panels (Deviance, Log-Odds), and Influence diagnostics (DFBETAs and Likelihood Displacement).
3. **Auditable R Objects**: Allows documenting qualitative justifications directly on the R object. Saving the object via `saveRDS()` preserves the complete, clinical-grade audit trail alongside the model.

---

## Installation

You can install the development version of `survAudit` directly from GitHub using the `remotes` package:

```R
# Install from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("millgreg/survAudit")
```

Alternatively, if you have the source archive downloaded, you can install it manually:

```R
# Install the built source archive locally
install.packages("survAudit_1.0.0.tar.gz", repos = NULL, type = "source")
```

---

## Quick Start

```R
library(survival)
library(survAudit)

# 1. Fit a naive Cox Proportional Hazards Model
data(pbc, package = "survival")
pbc$status <- ifelse(pbc$status == 2, 1, 0)
fit <- coxph(
  Surv(time, status) ~ age + bili + albumin + protime + factor(edema),
  data = pbc
)

# 2. Run the Diagnostic Audit
audit <- survAudit(fit, data = pbc)

# 3. View Compact Console Output
print(audit)

# 4. View Detailed Diagnostic Report
summary(audit)

# 5. Plot Diagnostics (ggplot2 panels)
plot(audit, which = "gof")         # Cox-Snell Overall Calibration
plot(audit, which = "ph")          # Proportional Hazards (Schoenfeld)
plot(audit, which = "functional")  # Linearity (Martingale)
plot(audit, which = "influence")   # DFBETAs & Likelihood Displacement
```

For a comprehensive walkthrough of interpreting these outputs, please read the [Introduction Vignette / Case Study](https://millgreg.github.io/survAudit/articles/survAudit-introduction.html).

### Documenting Qualitative Justifications

Check off outstanding non-identifiable assumptions by writing qualitative justifications directly into the audit object:

```R
# Document independent censoring
audit$assumptions$non_identifiable$independent_censoring$justification <- 
  "Censoring is administrative (end of study period) and patient drop-out is unrelated to disease severity."

# Document unmeasured confounding
audit$assumptions$non_identifiable$unmeasured_confounding$justification <- 
  "Baseline clinical confounders (performance score, age, celltype) were controlled."

# Print audit to see the checked results [x]
print(audit)

# Save the audit trail
saveRDS(audit, "cox_model_audit.rds")
```

---

## License

This package is licensed under the **GPL-3** License.
