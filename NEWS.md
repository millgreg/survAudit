# survAudit 0.1.0

* Initial release of `survAudit`, providing a unified diagnostic auditing framework for Cox Proportional Hazards models.
* Implements a structured assumption ontology categorizing assumptions into Statistically Assessable, Partially Assessable, and Non-Identifiable.
* Features automatic collinearity diagnostics (VIF/GVIF) with clean tabular outputs and appropriate scaling for multi-degree-of-freedom categorical variables.
* Contains advanced adjusted functional form checks using martingale residuals from reduced models compared with LOESS smooths.
* Includes influence diagnostics (DFBETAs and Likelihood Displacement) and multi-residual outlier detection (Martingale, Deviance, Log-Odds, Normal Deviate).
* Provides comprehensive S3 methods for `print()`, `summary()`, and `plot()` (generating faceted `ggplot2` diagnostic panels).
* Native support for plot pacing using standard R graphical device paging (`grDevices::devAskNewPage(TRUE)`).
