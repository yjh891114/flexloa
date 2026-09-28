# flexloa 0.1.0

* First release.
* `loa_fit()`, `loa_compare()`: conventional Bland-Altman limits on identity,
  square-root, cube-root, general n-th-root and log scales, with the classical
  diagnostics (Shapiro-Wilk, heteroscedasticity correlation, proportional-bias
  regression, empirical coverage overall and by thirds of the pair means)
  collected in a decision table ranked by residual heteroscedasticity.
* `loa_band()`: exact back-transformation of the limits by numerical inversion
  of the defining equation (`method = "exact"`, default) or the first-order
  formulas of Yoon, Yoon & Hahn (2019) (`method = "mvt"`).
* `loa_pairs()`: intra- and inter-reader difference pairs for replicated
  multi-reader designs (`replicates = "matched"`, `"mean"` or `"first"`).
* `loa_components()`: variance components (REML via lme4, or expected mean
  squares) with intra-reader, inter-reader and session-averaged inter-reader
  limits, optional residual trimming and residual diagnostics.
* `simulate_agreement()`: simulated paired data with error homoscedastic on a
  chosen transformed scale.
