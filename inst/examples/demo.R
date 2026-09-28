## flexloa demonstration -----------------------------------------------
library(flexloa)

## 1. Simulate paired readings whose error grows with lesion size
##    (homoscedastic on the square-root scale)
d <- simulate_agreement(nsubj = 250, transform = "sqrt", seed = 2026)

## 2. Compare candidate transformations with a formal decision table
cmp <- loa_compare(d$m1, d$m2,
                   transforms = c("identity", "sqrt", "cbrt", "log"))
print(cmp)          # normality, heteroscedasticity, proportional bias, coverage
plot(cmp)           # curved LOA overlays (cf. Yoon et al. 2019, Fig. 1)

## 3. Inspect the selected model
best <- cmp$fits[[cmp$best]]
print(best)
plot(best)
summary(best)                       # band at pretty grid of means
loa_band(best, mean = c(20, 60, 120))

## 4. Replicated multi-reader design: intra vs inter-reader LOA
set.seed(7)
truth <- rlnorm(60, 3.5, 0.5)
dat <- expand.grid(id = 1:60, reader = 1:4, rep = 1:2)
reader_eff <- rnorm(4, 0, 0.15)
dat$value <- (sqrt(truth[dat$id]) + reader_eff[dat$reader] +
              rnorm(nrow(dat), 0, 0.3))^2
vc <- loa_components(dat$value, dat$id, dat$reader, transform = "sqrt")
print(vc)
loa_band(vc$inter, mean = c(20, 60, 120))
