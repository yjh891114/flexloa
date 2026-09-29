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
  band <- loa_band(object, at)
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
  print(round(x$varcomp, digits))
  fmt <- function(o) sprintf("[%s, %s]",
                             format(o$loa_t[["lower"]], digits = digits),
                             format(o$loa_t[["upper"]], digits = digits))
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
#' @param xlab,ylab Axis labels.
#' @param ... Passed to \code{plot}.
#' @export
plot.loa_fit <- function(x, grid_n = 200,
                         xlab = "Average of the two measurements",
                         ylab = "Difference", ...) {
  m <- x$data$mean; d <- x$data$diff
  graphics::plot(m, d, xlab = xlab, ylab = ylab, pch = 1, ...)
  graphics::abline(h = 0, lty = 3, col = "grey60")
  g <- seq(max(min(m), 1e-9), max(m), length.out = grid_n)
  b <- loa_band(x, g)
  graphics::lines(b$mean, b$bias, lwd = 1.5)
  graphics::lines(b$mean, b$lower, lty = 2, lwd = 1.5)
  graphics::lines(b$mean, b$upper, lty = 2, lwd = 1.5)
  graphics::legend("topleft", bty = "n",
                   legend = c(sprintf("bias (%s)", x$spec$name),
                              sprintf("%.0f%% LOA", 100 * x$level)),
                   lty = c(1, 2), lwd = 1.5)
  invisible(x)
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
#' @param xlab,ylab Axis labels.
#' @param ... Passed to \code{plot}.
#' @export
plot.loa_compare <- function(x, grid_n = 200, cols = NULL, ltys = NULL,
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
