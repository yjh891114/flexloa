# Internal: solve, for a given mean m on the ORIGINAL scale, the difference
# D such that f(m + D/2) - f(m - D/2) = L, where L is a limit of agreement
# (or the bias) on the TRANSFORMED scale. Signed: a negative L yields a
# negative D by symmetry of the solver.
.solve_band <- function(spec, L, m) {
  if (!is.finite(m) || !is.finite(L)) return(NA_real_)
  if (spec$positive_only && m <= 0) return(NA_real_)
  if (L == 0) return(0)
  aL <- abs(L)
  g <- function(D) spec$f(m + D / 2) - spec$f(m - D / 2) - aL
  if (spec$positive_only) {
    upperD <- 2 * m * (1 - 1e-9)
    g_up <- suppressWarnings(g(upperD))
    if (!is.finite(g_up) || g_up < 0) {
      # For log-type transforms g -> +Inf as D -> 2m, so a non-finite value
      # means the limit is attainable arbitrarily close to the boundary;
      # a finite negative value means the limit is not attainable.
      if (is.finite(g_up)) return(NA_real_)
      # shrink until finite
      shrink <- 1 - 1e-6
      while (!is.finite(g(upperD)) && shrink > 0) {
        upperD <- 2 * m * shrink
        shrink <- shrink - 0.05
      }
      if (!is.finite(g(upperD)) || g(upperD) < 0) return(NA_real_)
    }
  } else {
    upperD <- 2 * aL
    while (g(upperD) < 0) upperD <- upperD * 2
  }
  root <- tryCatch(stats::uniroot(g, lower = 0, upper = upperD,
                                  tol = .Machine$double.eps^0.5)$root,
                   error = function(e) NA_real_)
  sign(L) * root
}

#' Curved limits of agreement on the original measurement scale
#'
#' Given a fitted \code{\link{loa_fit}} object (or any list carrying a
#' transformation \code{spec}, transformed-scale limits \code{loa_t} and
#' \code{bias}), computes the back-transformed bias curve and LOA band as
#' a function of the mean measurement on the original scale. For the
#' identity transformation the band reduces to the classical constant
#' Bland-Altman limits; for root or log transformations it is curved,
#' widening with the measurement scale.
#'
#' @param object A \code{loa_fit} object, one element of a
#'   \code{loa_compare} object, or a list with components \code{spec}
#'   (see \code{\link{loa_transform}}), \code{loa_t} (length-2 named
#'   numeric, \code{lower}/\code{upper}) and optionally \code{bias}.
#' @param mean Numeric vector of means (original scale) at which to
#'   evaluate the band.
#' @param method \code{"exact"} (default) solves the defining equation
#'   \eqn{f(m + D/2) - f(m - D/2) = L} numerically; \code{"mvt"} uses
#'   the first-order (mean value theorem) approximation
#'   \eqn{D = L / f'(m)}, which reproduces the closed-form expressions
#'   for root transformations given in Yoon et al. (2019).
#'
#' @references
#' Yoon JH, Yoon SH, Hahn S. Development of an algorithm for evaluating
#' the impact of measurement variability on response categorization in
#' oncology trials. BMC Med Res Methodol. 2019;19:90.
#'
#' @return A \code{data.frame} with columns \code{mean}, \code{bias},
#'   \code{lower}, \code{upper} (all on the original scale). \code{NA}
#'   is returned where a limit is not attainable within the domain of
#'   the transformation.
#' @export
loa_band <- function(object, mean, method = c("exact", "mvt")) {
  method <- match.arg(method)
  spec <- object$spec
  loa_t <- object$loa_t
  bias <- if (is.null(object$bias)) 0 else object$bias
  if (method == "mvt") {
    ## first-order: D = L / f'(m), with f'(m) = exp(dlog(m))
    fprime <- function(m) {
      out <- rep(NA_real_, length(m))
      ok <- is.finite(m) & (!spec$positive_only | m > 0)
      out[ok] <- exp(spec$dlog(m[ok]))
      out
    }
    fp <- fprime(mean)
    return(data.frame(mean = mean, bias = bias / fp,
                      lower = loa_t[["lower"]] / fp,
                      upper = loa_t[["upper"]] / fp))
  }
  data.frame(
    mean  = mean,
    bias  = vapply(mean, function(m) .solve_band(spec, bias, m), numeric(1)),
    lower = vapply(mean, function(m) .solve_band(spec, loa_t[["lower"]], m),
                   numeric(1)),
    upper = vapply(mean, function(m) .solve_band(spec, loa_t[["upper"]], m),
                   numeric(1)))
}
