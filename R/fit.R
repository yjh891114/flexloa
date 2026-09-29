#' Fit a transformation-based limits-of-agreement model to paired data
#'
#' Fits a measurement-error model on a chosen transformed scale to paired
#' measurements (two readings of the same quantity: two readers, two
#' occasions, or two methods), derives limits of agreement (LOA) on the
#' transformed scale, and back-transforms them into a (generally curved)
#' LOA band on the original scale. Alongside the limits, the function
#' reports the classical Bland-Altman diagnostics needed to decide
#' between candidate transformations:
#' \itemize{
#'   \item \strong{Shapiro-Wilk test} of the transformed-scale differences
#'     (normality of the error distribution on the working scale).
#'   \item \strong{Heteroscedasticity (size-dependency) test}: the
#'     Spearman rank correlation between the absolute centred differences
#'     and the pair means, as recommended by Bland and Altman (1999); a
#'     well-chosen variance-stabilising transformation should render this
#'     correlation negligible and non-significant.
#'   \item \strong{Proportional-bias test}: the ordinary least-squares
#'     regression of the transformed differences on the pair means
#'     (slope estimate and its p-value).
#'   \item \strong{Empirical coverage}: the proportion of observed
#'     original-scale differences falling inside the back-transformed
#'     band (target = \code{level}), reported overall
#'     (\code{coverage}) and within the lower, middle and upper thirds
#'     of the pair means (\code{cov_low}, \code{cov_mid},
#'     \code{cov_high}). The overall coverage of limits estimated from
#'     the same data is close to \code{level} for any transformation,
#'     so it serves only as a check; the coverage by thirds shows
#'     whether the band is too wide at one end of the measurement range
#'     and too narrow at the other, which is what a wrongly chosen scale
#'     produces.
#' }
#'
#' @param x,y Numeric vectors: the first and second measurement of each
#'   subject/lesion (original scale). Alternatively \code{x} may be a
#'   \code{\link{loa_pairs}} object, in which case \code{y} is omitted.
#' @param transform Transformation passed to \code{\link{loa_transform}}.
#' @param level Nominal LOA coverage level (default 0.95).
#' @param n Passed to \code{\link{loa_transform}} for \code{"root"}.
#'
#' @return An object of class \code{"loa_fit"}: a list with elements
#'   \code{spec}, \code{bias}, \code{sd}, \code{loa_t} (transformed
#'   scale), \code{level}, \code{data} (means and differences), and
#'   \code{gof} (a one-row data frame of the criteria above).
#'
#' @references
#' Bland JM, Altman DG. Statistical methods for assessing agreement
#' between two methods of clinical measurement. Lancet. 1986;327:307-10.
#'
#' Yoon JH, Yoon SH, Hahn S. Development of an algorithm for evaluating
#' the impact of measurement variability on response categorization in
#' oncology trials. BMC Med Res Methodol. 2019;19:90.
#'
#' @examples
#' d <- simulate_agreement(nsubj = 150, transform = "sqrt", seed = 1)
#' f <- loa_fit(d$m1, d$m2, transform = "sqrt")
#' f
#' head(loa_band(f, mean = c(20, 50, 100)))
#' @export
loa_fit <- function(x, y = NULL, transform = "sqrt", level = 0.95,
                    n = NULL) {
  xy <- .resolve_xy(x, y); x <- xy$x; y <- xy$y
  spec <- loa_transform(transform, n = n)
  if (length(x) != length(y))
    stop("`x` and `y` must have the same length.")
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  N <- length(x)
  if (N < 3L) stop("At least 3 complete pairs are required.")
  if (spec$positive_only && any(c(x, y) <= 0))
    stop(sprintf("Transformation '%s' requires strictly positive data.",
                 spec$name))

  tx <- spec$f(x); ty <- spec$f(y)
  d  <- ty - tx
  bias <- mean(d)
  s    <- stats::sd(d)
  z    <- stats::qnorm(1 - (1 - level) / 2)
  loa_t <- c(lower = bias - z * s, upper = bias + z * s)

  ## --- classical Bland-Altman diagnostics ------------------------------
  means <- (x + y) / 2

  ## (i) normality of the transformed differences
  ## shapiro.test() accepts at most 5000 values; use a fixed random
  ## subsample beyond that so that the test remains reproducible
  sw_p <- if (N >= 3) {
    dd <- if (N > 5000) {
      idx <- withr_free_sample(N, 5000)
      d[idx]
    } else d
    stats::shapiro.test(dd)$p.value
  } else NA_real_

  ## (ii) heteroscedasticity: Spearman correlation between |d - bias|
  ##      and the pair means (Bland & Altman 1999)
  ct <- suppressWarnings(
    stats::cor.test(abs(d - bias), means, method = "spearman"))

  ## (iii) proportional bias: regression of transformed differences on
  ##       the pair means
  pb <- stats::lm(d ~ means)
  pb_coef <- summary(pb)$coefficients
  pb_slope <- pb_coef["means", "Estimate"]
  pb_p     <- pb_coef["means", "Pr(>|t|)"]

  ## (iv) empirical coverage of the back-transformed band
  band <- loa_band(list(spec = spec, loa_t = loa_t, bias = bias), means)
  diffs <- y - x
  inside <- diffs >= band$lower & diffs <= band$upper
  coverage <- mean(inside, na.rm = TRUE)
  ## coverage within thirds of the pair means (lower / middle / upper)
  third <- .thirds(means)
  cov3 <- vapply(1:3, function(g) mean(inside[third == g], na.rm = TRUE),
                 numeric(1))

  gof <- data.frame(
    transform     = spec$name,
    n_pairs       = N,
    shapiro_p     = sw_p,
    hetero_rho    = unname(ct$estimate),
    hetero_p      = ct$p.value,
    prop_bias_slope = pb_slope,
    prop_bias_p   = pb_p,
    coverage      = coverage,
    cov_low       = cov3[1L],
    cov_mid       = cov3[2L],
    cov_high      = cov3[3L],
    stringsAsFactors = FALSE)

  out <- list(spec = spec, bias = bias, sd = s, loa_t = loa_t,
              level = level,
              data = data.frame(mean = means, diff = diffs, tdiff = d),
              gof = gof)
  class(out) <- "loa_fit"
  out
}

#' Compare candidate transformations for limits of agreement
#'
#' Fits \code{\link{loa_fit}} under each candidate transformation and
#' assembles a decision table of the classical Bland-Altman diagnostics
#' (normality of differences, heteroscedasticity test, proportional-bias
#' test, empirical coverage overall and by thirds of the pair means),
#' ordered by the absolute heteroscedasticity
#' correlation \code{|hetero_rho|} (smaller is better: the preferred
#' transformation is the one that removes the size-dependence of the
#' error). Transformations that cannot be applied to the data (e.g.
#' \code{"log"} with zeros) are skipped with a warning.
#'
#' @param x,y As in \code{\link{loa_fit}}; \code{x} may be a
#'   \code{\link{loa_pairs}} object.
#' @param transforms Character vector of candidate transformations.
#' @param level Nominal coverage level.
#' @param ... Passed to \code{\link{loa_fit}} (e.g. \code{n}).
#'
#' @return An object of class \code{"loa_compare"}: a list with
#'   \code{fits} (named list of \code{loa_fit} objects), \code{gof}
#'   (the combined decision table, sorted by \code{|hetero_rho|}), and
#'   \code{best} (name of the transformation with the smallest residual
#'   heteroscedasticity).
#'
#' @examples
#' d <- simulate_agreement(nsubj = 200, transform = "sqrt", seed = 42)
#' cmp <- loa_compare(d$m1, d$m2)
#' cmp$gof
#' plot(cmp)
#' @export
loa_compare <- function(x, y = NULL,
                        transforms = c("identity", "sqrt", "cbrt", "log"),
                        level = 0.95, ...) {
  xy <- .resolve_xy(x, y); x <- xy$x; y <- xy$y
  fits <- list()
  for (tr in transforms) {
    f <- tryCatch(loa_fit(x, y, transform = tr, level = level, ...),
                  error = function(e) {
                    warning(sprintf("Transformation '%s' skipped: %s",
                                    tr, conditionMessage(e)), call. = FALSE)
                    NULL
                  })
    if (!is.null(f)) fits[[f$gof$transform]] <- f
  }
  if (!length(fits)) stop("No transformation could be fitted.")
  gof <- do.call(rbind, lapply(fits, `[[`, "gof"))
  rownames(gof) <- NULL
  gof <- gof[order(abs(gof$hetero_rho)), ]
  out <- list(fits = fits, gof = gof, best = gof$transform[1L],
              level = level)
  class(out) <- "loa_compare"
  out
}

## group the pair means into lower / middle / upper thirds (1, 2, 3)
.thirds <- function(means) {
  br <- stats::quantile(means, c(1/3, 2/3), names = FALSE, type = 7)
  g <- 1L + (means > br[1L]) + (means > br[2L])
  as.integer(g)
}

## deterministic subsample of `k` indices out of `n` (does not touch the
## user's RNG stream)
withr_free_sample <- function(n, k) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(20190502)
  sample.int(n, k)
}
