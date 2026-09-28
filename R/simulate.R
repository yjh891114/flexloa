#' Simulate paired agreement data with size-dependent measurement error
#'
#' Generates paired measurements whose error is homoscedastic on a chosen
#' transformed scale (and therefore size-dependent on the original scale),
#' mimicking, e.g., tumour-diameter readings on CT.
#'
#' @param nsubj Number of subjects/units.
#' @param meanlog,sdlog Log-normal parameters for the true sizes.
#' @param error_sd Error SD on the transformed scale.
#' @param transform Scale on which the error is additive.
#' @param bias Systematic bias of the second reading (transformed scale).
#' @param seed Optional RNG seed.
#'
#' @return A data frame with columns \code{id}, \code{m1}, \code{m2}.
#' @examples
#' d <- simulate_agreement(100, transform = "sqrt", seed = 1)
#' @export
simulate_agreement <- function(nsubj = 100, meanlog = 3.5, sdlog = 0.6,
                               error_sd = 0.35, transform = "sqrt",
                               bias = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  spec <- loa_transform(transform)
  truth <- stats::rlnorm(nsubj, meanlog, sdlog)
  t1 <- spec$f(truth) + stats::rnorm(nsubj, 0, error_sd)
  t2 <- spec$f(truth) + bias + stats::rnorm(nsubj, 0, error_sd)
  eps <- 1e-6
  if (spec$positive_only) {
    t1 <- pmax(t1, spec$f(eps))
    t2 <- pmax(t2, spec$f(eps))
  }
  data.frame(id = seq_len(nsubj), m1 = spec$finv(t1), m2 = spec$finv(t2))
}
