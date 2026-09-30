#' @export
print.loa_fit <- function(x, digits = 4, ...) {
  cat("Transformation-based limits of agreement (flexloa)\n")
  cat(sprintf("  Transformation : %s\n", x$spec$name))
  cat(sprintf("  Pairs          : %d\n", x$gof$n_pairs))
  cat(sprintf("  Bias (transf.) : %s\n", format(x$bias, digits = digits)))
  cat(sprintf("  SD   (transf.) : %s\n", format(x$sd, digits = digits)))
  cat(sprintf("  %.0f%% LOA (transf.): [%s, %s]\n", 100 * x$level,
              format(x$loa_t[["lower"]], digits = digits),
              format(x$loa_t[["upper"]], digits = digits)))
  if (!is.null(x$ci_t)) {
    cat(sprintf("  %.0f%% CI: bias [%s, %s]; lower limit [%s, %s]; upper limit [%s, %s]\n",
                100 * x$conf,
                format(x$ci_t["bias", "lo"], digits = digits),
                format(x$ci_t["bias", "hi"], digits = digits),
                format(x$ci_t["lower", "lo"], digits = digits),
                format(x$ci_t["lower", "hi"], digits = digits),
                format(x$ci_t["upper", "lo"], digits = digits),
                format(x$ci_t["upper", "hi"], digits = digits)))
    cat(sprintf("  CI method: %s\n", x$ci_method))
  }
  cat("\nClassical Bland-Altman diagnostics:\n")
  print(x$gof, row.names = FALSE, digits = digits)
  invisible(x)
}

#' @export
summary.loa_fit <- function(object, at = NULL, ...) {
  if (is.null(at)) {
    rng <- range(object$data$mean)
    at <- pretty(rng, n = 5)
    at <- at[at >= rng[1] & at <= rng[2]]
    if (!length(at)) at <- mean(rng)
  }
  band <- loa_band(object, at, ci = !is.null(object$ci_t))
  cat("Back-transformed limits of agreement on the original scale:\n\n")
  print(band, row.names = FALSE, digits = 4)
  invisible(band)
}

#' @export
print.loa_compare <- function(x, digits = 4, ...) {
  cat("Comparison of candidate transformations (flexloa)\n")
  cat(sprintf("  Preferred transformation (smallest |hetero_rho|): %s\n\n", x$best))
  print(x$gof, row.names = FALSE, digits = digits)
  cat("\nGuidance: prefer non-significant heteroscedasticity (|rho| near 0),\n")
  cat("non-significant proportional bias, Shapiro p not small, and coverage\n")
  cat("near the nominal level in each third of the measurement range\n")
  cat("(cov_low, cov_mid, cov_high), not only overall.\n")
  invisible(x)
}

#' @export
print.loa_components <- function(x, digits = 4, ...) {
  cat("Variance components on the transformed scale (flexloa)\n")
  cat(sprintf("  Transformation : %s\n", x$spec$name))
  cat(sprintf("  Estimation     : %s\n", x$method))
  cat(sprintf("  Observations   : %d used", x$n_used))
  if (!is.null(x$n_trimmed) && x$n_trimmed > 0)
    cat(sprintf(", %d trimmed", x$n_trimmed))
  cat(sprintf("; %.2f replicates per unit-by-reader cell\n", x$k_mean))
  if (is.null(x$varcomp_ci)) {
    print(round(x$varcomp, digits))
  } else {
    tab <- cbind(estimate = x$varcomp, x$varcomp_ci)
    print(round(tab, digits))
    cat(sprintf("  (%.0f%% percentile intervals from %d cluster-bootstrap resamples of units)\n",
                100 * x$conf, x$nboot))
  }
  fmt <- function(o) {
    s <- sprintf("[%s, %s]", format(o$loa_t[["lower"]], digits = digits),
                 format(o$loa_t[["upper"]], digits = digits))
    if (!is.null(o$ci_t))
      s <- sprintf("%s  (upper limit %.0f%% CI %s to %s)", s, 100 * x$conf,
                   format(o$ci_t["upper", "lo"], digits = digits),
                   format(o$ci_t["upper", "hi"], digits = digits))
    s
  }
  lv <- 100 * x$level
  cat(sprintf("\n  Intra-reader %.0f%% LOA (transf.)               : %s\n",
              lv, fmt(x$intra)))
  cat(sprintf("  Inter-reader %.0f%% LOA (transf.), single read: %s\n",
              lv, fmt(x$inter)))
  cat(sprintf("  Inter-reader %.0f%% LOA (transf.), session mean: %s\n",
              lv, fmt(x$inter_mean)))
  if (!is.null(x$diag)) {
    cat("\nResidual diagnostics of the fitted model:\n")
    print(x$diag, row.names = FALSE, digits = digits)
  }
  cat("\nUse loa_band(x$intra, mean), loa_band(x$inter, mean) or\n")
  cat("loa_band(x$inter_mean, mean) for curved limits on the original scale.\n")
  invisible(x)
}

#' Bland-Altman plot with model-based curved limits of agreement
#'
#' @param x A \code{loa_fit} object.
#' @param grid_n Number of grid points for the band.
#' @param ci Logical; draw the pointwise confidence bands of the limits
#'   (and of the bias) as shaded ribbons (see \code{\link{loa_band}}).
#' @param ci_col Fill color of the ribbons.
#' @param xlab,ylab Axis labels.
#' @param ... Passed to \code{plot}.
#' @export
plot.loa_fit <- function(x, grid_n = 200, ci = FALSE,
                         ci_col = grDevices::adjustcolor("grey30", 0.2),
                         xlab = "Average of the two measurements",
                         ylab = "Difference", ...) {
  m <- x$data$mean; d <- x$data$diff
  graphics::plot(m, d, xlab = xlab, ylab = ylab, pch = 1, ...)
  graphics::abline(h = 0, lty = 3, col = "grey60")
  g <- seq(max(min(m), 1e-9), max(m), length.out = grid_n)
  b <- loa_band(x, g, ci = ci)
  if (ci) {
    .ribbon(b$mean, b$lower_lo, b$lower_hi, ci_col)
    .ribbon(b$mean, b$upper_lo, b$upper_hi, ci_col)
    .ribbon(b$mean, b$bias_lo, b$bias_hi, ci_col)
  }
  graphics::lines(b$mean, b$bias, lwd = 1.5)
  graphics::lines(b$mean, b$lower, lty = 2, lwd = 1.5)
  graphics::lines(b$mean, b$upper, lty = 2, lwd = 1.5)
  leg <- c(sprintf("bias (%s)", x$spec$name), sprintf("%.0f%% LOA", 100 * x$level))
  if (ci) leg[2] <- sprintf("%s with %.0f%% CI", leg[2], 100 * x$conf)
  graphics::legend("topleft", bty = "n", legend = leg, lty = c(1, 2), lwd = 1.5)
  invisible(x)
}

## shaded ribbon between two curves (NA segments are skipped)
.ribbon <- function(x, lo, hi, col) {
  ok <- is.finite(lo) & is.finite(hi)
  if (!any(ok)) return(invisible())
  r <- rle(ok); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  for (j in which(r$values)) {
    i <- starts[j]:ends[j]
    graphics::polygon(c(x[i], rev(x[i])), c(lo[i], rev(hi[i])),
                      col = col, border = NA)
  }
}

#' Overlay curved LOA bands from several transformations
#'
#' Reproduces the exploratory display of Yoon et al. (2019, Fig. 1):
#' one Bland-Altman scatter with the back-transformed LOA of every
#' candidate transformation superimposed.
#'
#' @param x A \code{loa_compare} object.
#' @param grid_n Number of grid points for the bands.
#' @param cols,ltys Optional colors / line types per transformation.
#' @param ci Logical; shade the pointwise confidence bands of the limits
#'   of the preferred transformation (\code{x$best}).
#' @param xlab,ylab Axis labels.
#' @param ... Passed to \code{plot}.
#' @export
plot.loa_compare <- function(x, grid_n = 200, cols = NULL, ltys = NULL,
                             ci = FALSE,
                             xlab = "Average of the two measurements",
                             ylab = "Difference", ...) {
  f1 <- x$fits[[1L]]
  m <- f1$data$mean; d <- f1$data$diff
  nt <- length(x$fits)
  if (is.null(cols)) cols <- seq_len(nt)
  if (is.null(ltys)) ltys <- seq_len(nt)
  graphics::plot(m, d, xlab = xlab, ylab = ylab, pch = 1, col = "grey40", ...)
  graphics::abline(h = 0, lty = 3, col = "grey70")
  g <- seq(max(min(m), 1e-9), max(m), length.out = grid_n)
  if (ci) {
    ib <- match(x$best, names(x$fits))
    bb <- loa_band(x$fits[[ib]], g, ci = TRUE)
    cc <- grDevices::adjustcolor(cols[ib], 0.2)
    .ribbon(bb$mean, bb$lower_lo, bb$lower_hi, cc)
    .ribbon(bb$mean, bb$upper_lo, bb$upper_hi, cc)
  }
  for (i in seq_len(nt)) {
    b <- loa_band(x$fits[[i]], g)
    graphics::lines(b$mean, b$lower, col = cols[i], lty = ltys[i], lwd = 1.6)
    graphics::lines(b$mean, b$upper, col = cols[i], lty = ltys[i], lwd = 1.6)
  }
  graphics::legend("topleft", bty = "n", legend = names(x$fits),
                   col = cols, lty = ltys, lwd = 1.6,
                   title = sprintf("%.0f%% LOA", 100 * x$level))
  invisible(x)
}
