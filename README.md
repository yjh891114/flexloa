# flexloa

**Conventional Bland–Altman limits of agreement (LOA), extended to error that grows with the magnitude of the measurement.**

`flexloa` operationalises the exploratory strategy used in Yoon, Yoon & Hahn
(2019, *BMC Medical Research Methodology* 19:90): instead of assuming constant
Bland–Altman limits, it computes the conventional limits under several candidate
variance-stabilising transformations (identity, √, ∛, general n-th root, log),
back-transforms each set of limits into a **curved LOA band on the
original scale**, and reports the **goodness-of-fit evidence needed to choose
between them**:

| Criterion | What it tells you |
|---|---|
| Shapiro–Wilk p (transformed differences) | Normality of errors on the chosen scale |
| Heteroscedasticity ρ, p (Spearman, \|d − bias\| vs mean; Bland & Altman 1999) | Whether the transformation removed the error's size-dependence |
| Proportional-bias slope, p (regression of d on mean) | Whether bias changes with the measurement scale |
| Empirical coverage, overall and by thirds of the pair means | % of observed differences inside the back-transformed band (target ≈ level). Overall coverage is ≈ level for any transformation fitted to the same data; the coverage in the lower / middle / upper third of the measurement range shows whether the band is too wide at one end and too narrow at the other |

## Worked example

The package ships a synthetic multi-reader CT dataset, `recist_readings`, and a
vignette that walks through it from the conventional (constant) limits to the
magnitude-dependent limits, and shows where the two lead to different clinical
decisions:

```r
vignette("recist-example", package = "flexloa")
```

## Installation

```r
# from GitHub
# install.packages("remotes")
remotes::install_github("yjh891114/flexloa")

# or from a source tarball
install.packages("flexloa_0.1.0.tar.gz", repos = NULL, type = "source")
```

Only base R is required. `lme4` (optional) enables REML variance components
for replicated multi-reader designs.

## Quick start

```r
library(flexloa)

d   <- simulate_agreement(nsubj = 250, transform = "sqrt", seed = 2026)
cmp <- loa_compare(d$m1, d$m2)   # identity, sqrt, cbrt, log by default

cmp$gof        # decision table, sorted by |hetero_rho|
plot(cmp)      # Bland-Altman scatter + curved LOA per transformation
best <- cmp$fits[[cmp$best]]
plot(best); summary(best)
loa_band(best, mean = c(20, 60, 120))   # limits at chosen sizes
loa_band(best, mean = c(20, 60, 120), method = "mvt")  # first-order formulas
```

## Replicated multi-reader designs

```r
# exploratory: reduce the reader study to difference pairs
p_intra <- loa_pairs(value, id, reader, session, type = "intra")
p_inter <- loa_pairs(value, id, reader, session, type = "inter")
loa_compare(p_intra); loa_compare(p_inter)

# definitive: variance components on the chosen scale
vc <- loa_components(value, id = lesion, reader = radiologist,
                     transform = "sqrt", trim = 0)
vc                               # components, residual diagnostics (REML if lme4 present)
loa_band(vc$intra, mean = 50)    # repeat-reading (intra-reader) limits
loa_band(vc$inter, mean = 50)    # between-reader limits, single readings
loa_band(vc$inter_mean, mean = 50)  # between-reader limits, session-averaged readings
```

## Notes

* Man pages are written as roxygen2 comments; run `devtools::document()` to
  regenerate `man/` before `R CMD check`.
* Transformations are ranked by the classical heteroscedasticity criterion
  (Spearman correlation between absolute centred differences and means):
  the preferred scale is the one on which the error no longer depends on
  the size of the measurement.
* `loa_band()` returns `NA` where a limit is not attainable within the domain
  of the transformation (e.g. very small means under a root transform).

## References

* Yoon JH, Yoon SH, Hahn S. BMC Med Res Methodol. 2019;19:90.
* Bland JM, Altman DG. Lancet. 1986;327:307–10.
* Euser AM, Dekker FW, le Cessie S. J Clin Epidemiol. 2008;61:978–82.
