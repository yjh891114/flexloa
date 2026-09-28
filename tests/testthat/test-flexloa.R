test_that("identity band reproduces constant Bland-Altman limits", {
  d <- simulate_agreement(120, transform = "identity",
                          error_sd = 3, meanlog = 4, seed = 11)
  f <- loa_fit(d$m1, d$m2, transform = "identity")
  b <- loa_band(f, mean = c(20, 50, 100))
  expect_equal(b$lower, rep(f$loa_t[["lower"]], 3), tolerance = 1e-6)
  expect_equal(b$upper, rep(f$loa_t[["upper"]], 3), tolerance = 1e-6)
})

test_that("band solver satisfies the defining equation", {
  f <- list(spec = loa_transform("sqrt"),
            loa_t = c(lower = -0.6, upper = 0.6), bias = 0)
  b <- loa_band(f, mean = c(25, 64, 100))
  for (i in seq_len(nrow(b))) {
    m <- b$mean[i]; D <- b$upper[i]
    expect_equal(sqrt(m + D / 2) - sqrt(m - D / 2), 0.6, tolerance = 1e-6)
  }
})

test_that("log band matches its closed form", {
  spec <- loa_transform("log")
  L <- 0.4; m <- 50
  D_closed <- 2 * m * (exp(L) - 1) / (exp(L) + 1)
  b <- loa_band(list(spec = spec, loa_t = c(lower = -L, upper = L)), m)
  expect_equal(b$upper, D_closed, tolerance = 1e-6)
  expect_equal(b$lower, -D_closed, tolerance = 1e-6)
})

test_that("heteroscedasticity criterion selects the generating transformation (usually)", {
  d <- simulate_agreement(400, transform = "sqrt", seed = 99)
  cmp <- loa_compare(d$m1, d$m2)
  expect_true(cmp$best %in% c("sqrt", "cbrt"))
  expect_true(all(c("shapiro_p", "hetero_rho", "hetero_p",
                    "prop_bias_slope", "prop_bias_p", "coverage",
                    "cov_low", "cov_mid", "cov_high") %in% names(cmp$gof)))
})

test_that("empirical coverage is near nominal", {
  d <- simulate_agreement(500, transform = "sqrt", seed = 3)
  f <- loa_fit(d$m1, d$m2, transform = "sqrt")
  expect_gt(f$gof$coverage, 0.90)
  expect_lt(f$gof$coverage, 0.99)
  ## coverage by thirds: computed from the same band, so the overall
  ## value is the (weighted) mean of the three
  m <- f$data$mean
  th <- cut(m, stats::quantile(m, c(0, 1/3, 2/3, 1)), include.lowest = TRUE)
  w <- as.numeric(table(th)) / length(m)
  expect_equal(sum(w * c(f$gof$cov_low, f$gof$cov_mid, f$gof$cov_high)),
               f$gof$coverage, tolerance = 1e-8)
  ## constant limits on size-dependent errors: too wide for small,
  ## too narrow for large values
  fi <- loa_fit(d$m1, d$m2, transform = "identity")
  expect_gt(fi$gof$cov_low, fi$gof$cov_high)
})

test_that("loa_components returns sensible variance components", {
  set.seed(7)
  truth <- rlnorm(40, 3.5, 0.5)
  dat <- expand.grid(id = 1:40, reader = 1:3, rep = 1:2)
  reader_eff <- rnorm(3, 0, 0.15)
  dat$value <- (sqrt(truth[dat$id]) + reader_eff[dat$reader] +
                rnorm(nrow(dat), 0, 0.3))^2
  vc <- loa_components(dat$value, dat$id, dat$reader, transform = "sqrt")
  expect_true(all(vc$varcomp >= 0))
  expect_gt(vc$inter$loa_t[["upper"]], vc$intra$loa_t[["upper"]])
})

test_that("loa_pairs forms the expected numbers of pairs", {
  set.seed(3)
  truth <- rlnorm(20, 3.5, 0.5)
  dat <- expand.grid(id = 1:20, reader = 1:4, session = 1:2)
  dat$value <- (sqrt(truth[dat$id]) + rnorm(nrow(dat), 0, 0.3))^2
  p_intra <- loa_pairs(dat$value, dat$id, dat$reader, dat$session, type = "intra")
  p_inter <- loa_pairs(dat$value, dat$id, dat$reader, dat$session, type = "inter")
  p_mean  <- loa_pairs(dat$value, dat$id, dat$reader, dat$session,
                       type = "inter", replicates = "mean")
  p_first <- loa_pairs(dat$value, dat$id, dat$reader, dat$session,
                       type = "inter", replicates = "first")
  expect_equal(nrow(p_intra), 20 * 4)          # one pair per unit x reader
  expect_equal(nrow(p_inter), 20 * 6 * 2)      # 6 reader pairs, 2 sessions
  expect_equal(nrow(p_mean), 20 * 6)
  expect_equal(nrow(p_first), 20 * 6)
  expect_s3_class(p_intra, "loa_pairs")
  cmp <- loa_compare(p_intra)
  expect_s3_class(cmp, "loa_compare")
  expect_equal(loa_fit(p_intra)$gof$n_pairs, 80)
})

test_that("method = 'mvt' reproduces the first-order root formulas", {
  f <- list(spec = loa_transform("sqrt"),
            loa_t = c(lower = -0.6, upper = 0.6), bias = 0.1)
  m <- c(25, 64, 100)
  b <- loa_band(f, m, method = "mvt")
  expect_equal(b$upper, 2 * sqrt(m) * 0.6, tolerance = 1e-10)
  expect_equal(b$bias,  2 * sqrt(m) * 0.1, tolerance = 1e-10)
  ## agrees with the exact solution to first order
  e <- loa_band(f, m)
  expect_true(all(abs(b$upper / e$upper - 1) < 0.02))
  ## identity: exact and mvt coincide
  g <- list(spec = loa_transform("identity"), loa_t = c(lower = -3, upper = 3))
  expect_equal(loa_band(g, m, method = "mvt")$upper, rep(3, 3))
})

test_that("loa_components: trim, inter_mean and residual diagnostics", {
  set.seed(21)
  truth <- rlnorm(30, 3.5, 0.5)
  dat <- expand.grid(id = 1:30, reader = 1:3, session = 1:2)
  dat$value <- (sqrt(truth[dat$id]) + rnorm(3, 0, 0.2)[dat$reader] +
                rnorm(nrow(dat), 0, 0.3))^2
  vc <- loa_components(dat$value, dat$id, dat$reader, transform = "sqrt")
  expect_equal(vc$k_mean, 2)
  expect_equal(vc$n_trimmed, 0L)
  expect_true(all(c("shapiro_p", "hetero_rho", "hetero_p") %in% names(vc$diag)))
  ## session-averaged inter-reader limits are narrower than single-read
  expect_lt(vc$inter_mean$loa_t[["upper"]], vc$inter$loa_t[["upper"]])
  expect_gt(vc$inter_mean$loa_t[["upper"]], 0)
  vt <- loa_components(dat$value, dat$id, dat$reader, transform = "sqrt",
                       trim = 0.05)
  expect_equal(vt$n_trimmed, as.integer(floor(0.05 * nrow(dat))))
  expect_equal(vt$n_used, nrow(dat) - vt$n_trimmed)
  expect_error(loa_components(dat$value, dat$id, dat$reader, trim = 0.7))
})

test_that("loa_components handles a single-session design", {
  set.seed(5)
  truth <- rlnorm(40, 3.5, 0.5)
  dat <- expand.grid(id = 1:40, reader = 1:3)
  dat$value <- (sqrt(truth[dat$id]) + rnorm(3, 0, 0.2)[dat$reader] +
                rnorm(nrow(dat), 0, 0.3))^2
  vc <- loa_components(dat$value, dat$id, dat$reader, transform = "sqrt")
  expect_equal(unname(vc$varcomp[["id_reader"]]), 0)
  expect_equal(vc$k_mean, 1)
  expect_gt(vc$inter$loa_t[["upper"]], 0)
})
