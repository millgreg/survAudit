# survAudit

[![R-CMD-check](https://github.com/millgreg/survAudit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/millgreg/survAudit/actions/workflows/R-CMD-check.yaml)

> Comprehensive Diagnostics for Cox Proportional Hazards Models

`survAudit` provides a unified diagnostic auditing framework for Cox proportional hazards (PH) models. It aggregates statistical diagnostics into a structured **Assumption Classification** and provides a mechanism to document qualitative justifications for non-identifiable assumptions directly on the R model audit object.

Full documentation and a reproducible case study are available at: **[https://millgreg.github.io/survAudit](https://millgreg.github.io/survAudit)**

## Features

1. **Unified Assumption Classification**: Explicitly organizes model assumptions into:
   * **Non-Identifiable**: Assumptions that cannot be verified statistically (e.g., Independent Censoring, Absence of Unmeasured Confounding) and require qualitative justification.
   * **Partially Assessable**: Assumptions informed by metrics but requiring clinical/domain judgment (e.g., Outlier Impact, Missing Data).
   * **Statistically Assessable**: Assumptions rigorously testable from data (e.g., Proportional Hazards, Functional Form, Influence Stability, Event Sufficiency).
2. **Aggregated Diagnostics**:
   * **Collinearity (GVIF)**: Implements the Generalized Variance Inflation Factor (GVIF) algorithm (Fox & Monette, 1992) to group dummy variables of factor predictors.
   * **Functional Form (Linearity)**: Martingale residuals from multivariate reduced models are compared against continuous predictors using LOESS smooths.
   * **Visual Diagnostics**: Overall goodness-of-fit via Cox-Snell residuals, outlier panels, and influence diagnostics.
3. **Auditable R Objects**: Allows documenting qualitative justifications directly on the R object. Saving the object via `saveRDS()` preserves the audit trail alongside the model.

---

## Methodological Notes

The diagnostics in `survAudit` are statistical tests and heuristics. They have inherent limitations that require critical scientific judgment:

* **Proportional Hazards**: The global and variable-specific PH tests (via `cox.zph()`) only check for *linear* deviations from proportional hazards. Furthermore, these tests can be overly sensitive in large sample sizes, flagging minor, clinically irrelevant violations. Visual inspection of the scaled Schoenfeld residuals is strictly recommended.
* **Functional Form**: LOESS smoothing of martingale residuals is highly sensitive to the distribution of covariates. Confidence bands widen significantly in the tails where data is sparse, so perceived non-linearities at the extremes should be interpreted with caution.
* **Non-testability**: Statistical models cannot verify their own causal structure. The assumptions of independent censoring and unmeasured confounding *cannot* be verified from the observed data and rely entirely on study design and domain knowledge.

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
install.packages("survAudit_1.1.0.tar.gz", repos = NULL, type = "source")
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
plot(audit, which = "influence")   # DFBETAs
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

## Citation

If you use `survAudit` in your research, please cite it:

```bibtex
@Manual{,
  title = {survAudit: Comprehensive Diagnostics for Cox Proportional Hazards Models},
  author = {Gregor Miller},
  year = {2026},
  note = {R package version 1.1.0},
  url = {https://millgreg.github.io/survAudit},
}
```

## Dependencies

`survAudit` builds heavily on the foundations provided by the [`survival`](https://cran.r-project.org/package=survival) package (Therneau, 2024) for model fitting and residual calculations. Visual diagnostics are implemented using `ggplot2` and `patchwork`.

---

## License

This package is licensed under the **GPL-3** License.
