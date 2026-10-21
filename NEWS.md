# survAudit 1.3.1

* **Ontology Alignment**: Reclassified influence stability from statistically assessable to partially assessable in `audit$assumptions` and console `print()`, aligning with the framework distinction that DFBETAs threshold screening requires qualitative scrutiny rather than formal hypothesis testing.

# survAudit 1.3.0

* **Smoothing**: Added `smooth_method` parameter to `plot.survAudit()` (default `"auto"`). Smoothers now automatically switch from LOESS to GAM for panels with more than 1,000 observations, matching `ggplot2` conventions. Smoothers are always fitted on the full dataset, independent of `max_points` dot subsampling.
* **Performance**: Removed redundant base LOESS fitting in `survAudit()` functional form diagnostics, eliminating unnecessary compute overhead and memory bloat.
* **Internal Cleanup**: Added `@noRd` to internal plotting and printing helpers to remove orphan `.Rd` files, removed unused internal functions and arguments, and aligned documentation with mandatory `data` argument requirements.
* **Documentation**: Updated vignette, help pages, and internal summaries to use method-agnostic "smoothed curve" terminology.

# survAudit 1.2.0

* **Subsampling**: Added `max_points` (default 5000) to `plot.survAudit()` to optimize plotting in large datasets while retaining all flagged outliers and influential cases.
* **Console Output**: Wrapped outlier summary in `print.survAudit()` across two lines to fit standard console widths.
* **Documentation**: Refined vignette and case study terminology, added section links, and updated assumption justifications.

# survAudit 1.1.0

* **Bug Fix**: Fixed indexing logic for all influence and outlier diagnostic outputs to use the original row names instead of complete-case row numbers, preventing incorrect matching when models drop missing data.
* **Bug Fix**: Functional form logic now requires the `data` argument and successfully refits restricted models even when the dataset contains missing rows.
* **Bug Fix**: Proportional hazards visualization now successfully plots un-aggregated scaled Schoenfeld residuals for all dummy variables of a categorical covariate.
* **Documentation**: The primary vignette and case study were rewritten to focus strictly on practical diagnostic interpretation, completely removing promotional language. Added explicit methodological caveats regarding the sample-size sensitivity and linear restrictions of the `cox.zph` test.

# survAudit 0.1.0

* Initial release of `survAudit`, providing a unified diagnostic auditing framework for Cox Proportional Hazards models.
* Implements a structured assumption ontology categorizing assumptions into Statistically Assessable, Partially Assessable, and Non-Identifiable.
* Features automatic collinearity diagnostics (VIF/GVIF) with clean tabular outputs and appropriate scaling for multi-degree-of-freedom categorical variables.
* Contains advanced adjusted functional form checks using martingale residuals from reduced models compared with LOESS smooths.
* Includes influence diagnostics (DFBETAs) and multi-residual outlier detection (Martingale, Deviance, Log-Odds, Normal Deviate).
* Provides comprehensive S3 methods for `print()`, `summary()`, and `plot()` (generating faceted `ggplot2` diagnostic panels).
* Native support for plot pacing using standard R graphical device paging (`grDevices::devAskNewPage(TRUE)`).
