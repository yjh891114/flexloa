#' Transformation specifications for limits-of-agreement modeling
#'
#' Builds the transformation object used throughout \pkg{flexloa}. Each
#' specification carries the forward transformation \code{f}, its inverse
#' \code{finv}, and \code{dlog}, the log of the absolute derivative of
#' \code{f} (retained for completeness and for user-supplied extensions).
#'
#' @param transform Character, one of \code{"identity"}, \code{"sqrt"},
#'   \code{"cbrt"}, \code{"log"}, \code{"root"}; or a user-supplied list
#'   with elements \code{name}, \code{f}, \code{finv}, \code{dlog}.
#' @param n Root order, required when \code{transform = "root"} (e.g.
#'   \code{n = 2} reproduces \code{"sqrt"}).
#'
#' @return A list of class \code{"loa_transform"} with elements
#'   \code{name}, \code{f}, \code{finv}, \code{dlog}, and
#'   \code{positive_only} (logical; whether the transform requires
#'   strictly positive data).
#'
#' @references
#' Euser AM, Dekker FW, le Cessie S. A practical approach to Bland-Altman
#' plots and variation coefficients for log transformed variables.
#' J Clin Epidemiol. 2008;61(10):978-82.
#'
#' @examples
#' sp <- loa_transform("sqrt")
#' sp$f(9)      # 3
#' sp$finv(3)   # 9
#' @export
loa_transform <- function(transform = c("identity", "sqrt", "cbrt",
                                        "log", "root"),
                          n = NULL) {
  if (is.list(transform)) {
    stopifnot(all(c("name", "f", "finv", "dlog") %in% names(transform)))
    if (is.null(transform$positive_only)) transform$positive_only <- TRUE
    class(transform) <- "loa_transform"
    return(transform)
  }
  transform <- match.arg(transform)
  spec <- switch(transform,
    identity = list(
      name = "identity",
      f    = function(y) y,
      finv = function(z) z,
      dlog = function(y) rep(0, length(y)),
      positive_only = FALSE),
    sqrt = list(
      name = "sqrt",
      f    = function(y) sqrt(y),
      finv = function(z) z^2,
      dlog = function(y) log(0.5) - 0.5 * log(y),
      positive_only = TRUE),
    cbrt = list(
      name = "cbrt",
      f    = function(y) y^(1/3),
      finv = function(z) z^3,
      dlog = function(y) log(1/3) - (2/3) * log(y),
      positive_only = TRUE),
    log = list(
      name = "log",
      f    = function(y) log(y),
      finv = function(z) exp(z),
      dlog = function(y) -log(y),
      positive_only = TRUE),
    root = {
      if (is.null(n) || n <= 0)
        stop("`n` (a positive root order) must be supplied for transform = 'root'.")
      force(n)
      list(
        name = paste0("root", n),
        f    = function(y) y^(1/n),
        finv = function(z) z^n,
        dlog = function(y) log(1/n) + (1/n - 1) * log(y),
        positive_only = TRUE)
    })
  class(spec) <- "loa_transform"
  spec
}
