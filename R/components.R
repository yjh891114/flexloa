#' Intra- and inter-reader limits of agreement from replicated designs
#'
#' Estimates variance components on a transformed scale from a replicated
#' multi-reader design (each unit measured repeatedly by several readers)
#' and returns intra-reader and inter-reader limits of agreement, each
#' back-transformable to a curved band with \code{\link{loa_band}}.
#' Following the hierarchical model of Yoon et al. (2019), the model is
#' \deqn{f(Y_{ijk}) = \mu + a_i + b_j + g_{ij} + e_{ijk},}
#' with unit (\code{id}), reader, and unit-by-reader interaction random
#' effects. The intra-reader error variance of a difference between two
#' replicate readings by the same reader is \eqn{2\sigma^2_e}; the
#' inter-reader error variance of a difference between single readings by
#' two readers is \eqn{2(\sigma^2_b + \sigma^2_g + \sigma^2_e)}; and the
#' inter-reader error variance of a difference between two readers'
#' session-averaged readings is
#' \eqn{2(\sigma^2_b + \sigma^2_g + \sigma^2_e / \bar{k})}, where
#' \eqn{\bar{k}} is the mean number of replicate readings per
#' unit-by-reader cell (\code{inter_mean}).
#'
#' If \pkg{lme4} is installed the components are estimated by REML;
#' otherwise a method-of-moments (expected mean squares) estimator is
#' used, which is exact for balanced designs and approximate otherwise.
#' When no unit is read more than once by the same reader (a single
#' session), the unit-by-reader interaction cannot be separated from the
#' residual; it is then omitted from the model and reported as zero, and
#' the intra-reader limits are not meaningful.
#'
#' Because the definitive limits derive from the fitted model, two
#' residual-based diagnostics of that model are reported alongside the
#' components: the Shapiro-Wilk p-value of the residuals and the Spearman
#' rank correlation between absolute residuals and fitted values (a check
#' that the chosen scale has removed the size-dependence of the error at
#' the level of the hierarchical model, not only of the exploratory
#' pairs).
#'
#' @param value Numeric vector of measurements (original scale).
#' @param id Factor-like: measured unit (lesion, subject, sample).
#' @param reader Factor-like: reader/observer/device.
#' @param transform,n Passed to \code{\link{loa_transform}}.
#' @param level Nominal coverage level.
#' @param trim Fraction (0 to 0.5) of observations with the largest
#'   absolute standardized residuals from an initial fit to remove before
#'   the final fit; \code{0} (default) fits all observations once.
#'
#' @return An object of class \code{"loa_components"}: a list with
#'   \code{spec}, \code{varcomp} (named numeric: \code{id},
#'   \code{reader}, \code{id_reader}, \code{residual}), \code{method}
#'   ("REML (lme4)" or "Method of moments (EMS)"), \code{k_mean} (mean
#'   replicates per unit-by-reader cell), \code{n_used},
#'   \code{n_trimmed}, \code{diag} (a one-row data frame with
#'   \code{shapiro_p}, \code{hetero_rho}, \code{hetero_p} computed from
#'   the residuals of the final fit), \code{fit} (the fitted
#'   \code{lmer} or \code{lm} object of the final fit), and three sub-objects
#'   \code{intra}, \code{inter} and \code{inter_mean}, each a list with
#'   \code{spec}, \code{bias = 0} and \code{loa_t}, suitable for
#'   \code{\link{loa_band}}.
#'
#' @references
#' Yoon JH, Yoon SH, Hahn S. Development of an algorithm for evaluating
#' the impact of measurement variability on response categorization in
#' oncology trials. BMC Med Res Methodol. 2019;19:90.
#'
#' @examples
#' set.seed(7)
#' truth <- rlnorm(40, 3.5, 0.5)
#' dat <- expand.grid(id = 1:40, reader = 1:3, rep = 1:2)
#' dat$value <- (sqrt(truth[dat$id]) + rnorm(3, 0, 0.15)[dat$reader] +
#'               rnorm(nrow(dat), 0, 0.3))^2
#' vc <- loa_components(dat$value, dat$id, dat$reader, transform = "sqrt")
#' vc
#' loa_band(vc$inter, mean = c(20, 50, 100))
#' loa_band(vc$inter_mean, mean = c(20, 50, 100))
#' @export
loa_components <- function(value, id, reader, transform = "sqrt",
                           level = 0.95, n = NULL, trim = 0) {
  spec <- loa_transform(transform, n = n)
  if (!is.numeric(trim) || length(trim) != 1L || trim < 0 || trim >= 0.5)
    stop("`trim` must be a single number in [0, 0.5).")
  ok <- is.finite(value) & !is.na(id) & !is.na(reader)
  value <- value[ok]
  id <- factor(id[ok]); reader <- factor(reader[ok])
  if (spec$positive_only && any(value <= 0))
    stop(sprintf("Transformation '%s' requires strictly positive data.",
                 spec$name))
  ty <- spec$f(value)
  use_lme4 <- requireNamespace("lme4", quietly = TRUE)

  fit_once <- function(ty, id, reader) {
    id <- droplevels(id); reader <- droplevels(reader)
    ## without replicate readings the unit-by-reader interaction is
    ## confounded with the residual: drop it and report it as zero
    replicated <- any(table(id, reader) > 1)
    if (use_lme4) {
      dat <- data.frame(ty = ty, id = id, reader = reader)
      form <- if (replicated)
        ty ~ 1 + (1 | id) + (1 | reader) + (1 | id:reader)
      else ty ~ 1 + (1 | id) + (1 | reader)
      m <- lme4::lmer(form, data = dat, REML = TRUE)
      vc <- as.data.frame(lme4::VarCorr(m))
      getv <- function(g) {
        i <- match(g, vc$grp)
        if (is.na(i)) 0 else vc$vcov[i]
      }
      list(v_id = getv("id"), v_rd = getv("reader"),
           v_int = getv("id:reader"), v_e = getv("Residual"),
           resid = stats::residuals(m), fitted = stats::fitted(m),
           fit = m, method = "REML (lme4)")
    } else {
      ## Method of moments via expected mean squares (balanced two-way
      ## random model with replication).
      fit <- if (replicated) stats::lm(ty ~ id * reader) else stats::lm(ty ~ id + reader)
      an  <- stats::anova(fit)
      ms  <- an[["Mean Sq"]]; names(ms) <- rownames(an)
      s_lev <- nlevels(id); r_lev <- nlevels(reader)
      k <- length(ty) / (s_lev * r_lev)   # average replicates per cell
      v_e   <- ms[["Residuals"]]
      if (replicated) {
        v_int <- max(0, (ms[["id:reader"]] - v_e) / k)
        v_rd  <- max(0, (ms[["reader"]] - ms[["id:reader"]]) / (s_lev * k))
        v_id  <- max(0, (ms[["id"]] - ms[["id:reader"]]) / (r_lev * k))
      } else {
        v_int <- 0
        v_rd  <- max(0, (ms[["reader"]] - v_e) / s_lev)
        v_id  <- max(0, (ms[["id"]] - v_e) / r_lev)
      }
      list(v_id = v_id, v_rd = v_rd, v_int = v_int, v_e = v_e,
           resid = stats::residuals(fit), fitted = stats::fitted(fit),
           fit = fit, method = "Method of moments (EMS)")
    }
  }

  ## optional trimming on standardized residuals of an initial fit
  n_trimmed <- 0L
  if (trim > 0) {
    f0 <- fit_once(ty, id, reader)
    z  <- abs(f0$resid) / stats::sd(f0$resid)
    n_trimmed <- as.integer(floor(trim * length(ty)))
    if (n_trimmed > 0L) {
      drop <- order(z, decreasing = TRUE)[seq_len(n_trimmed)]
      keep <- setdiff(seq_along(ty), drop)
      ty <- ty[keep]; id <- id[keep]; reader <- reader[keep]
    }
  }
  f <- fit_once(ty, id, reader)
  k_mean <- length(ty) / (nlevels(droplevels(id)) *
                          nlevels(droplevels(reader)))

  ## residual diagnostics of the final model
  r <- f$resid; fv <- f$fitted
  sw_p <- if (length(r) >= 3 && length(r) <= 5000)
    stats::shapiro.test(r)$p.value else NA_real_
  ct <- suppressWarnings(stats::cor.test(abs(r), fv, method = "spearman"))
  diag <- data.frame(shapiro_p = sw_p, hetero_rho = unname(ct$estimate),
                     hetero_p = ct$p.value)

  z <- stats::qnorm(1 - (1 - level) / 2)
  mk <- function(v) list(spec = spec, bias = 0,
                         loa_t = c(lower = -z * sqrt(v), upper = z * sqrt(v)))
  intra      <- mk(2 * f$v_e)
  inter      <- mk(2 * (f$v_rd + f$v_int + f$v_e))
  inter_mean <- mk(2 * (f$v_rd + f$v_int + f$v_e / k_mean))

  out <- list(spec = spec,
              varcomp = c(id = f$v_id, reader = f$v_rd,
                          id_reader = f$v_int, residual = f$v_e),
              method = f$method, level = level, k_mean = k_mean,
              n_used = length(ty), n_trimmed = n_trimmed, diag = diag,
              fit = f$fit,
              intra = intra, inter = inter, inter_mean = inter_mean)
  class(out) <- "loa_components"
  out
}
