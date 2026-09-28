#' Convert a general reader study into difference pairs
#'
#' Reduces long-format data from a general reader study -- any number of
#' units, readers and replicate sessions -- to the paired form used by
#' \code{\link{loa_fit}} and \code{\link{loa_compare}}. Two kinds of
#' pairs are formed:
#' \itemize{
#'   \item \code{type = "intra"}: within each reader, all pairs of that
#'     reader's replicate readings of the same unit (repeat-reading
#'     variability).
#'   \item \code{type = "inter"}: within each unit, all pairs of readers
#'     (between-reader variability). With replicate sessions, the
#'     \code{replicates} argument decides how the readers' replicates
#'     enter: session-matched pairs (\code{"matched"}, the default:
#'     reader A session s versus reader B session s), the per-reader
#'     average over sessions (\code{"mean"}), or the first session only
#'     (\code{"first"}).
#' }
#' Pairs sharing a unit are not mutually independent, so the diagnostics
#' computed from them by \code{\link{loa_compare}} are exploratory --
#' their role is to choose the transformation -- while definitive
#' variance estimation should use \code{\link{loa_components}}.
#'
#' @param value Numeric vector of measurements (original scale).
#' @param id Factor-like: measured unit (lesion, subject, sample).
#' @param reader Factor-like: reader/observer/device.
#' @param session Factor-like: replicate session/occasion. May be
#'   \code{NULL} when every reader measured each unit once (then only
#'   \code{type = "inter"} is possible).
#' @param type \code{"intra"} or \code{"inter"}; see Details.
#' @param replicates How replicate sessions enter inter-reader pairs:
#'   \code{"matched"}, \code{"mean"} or \code{"first"}.
#'
#' @return A \code{data.frame} of class \code{"loa_pairs"} with columns
#'   \code{x}, \code{y} (the two measurements of a pair, original
#'   scale), \code{id}, and either \code{reader} plus \code{session_x},
#'   \code{session_y} (intra) or \code{reader_x}, \code{reader_y} plus
#'   \code{session} (inter). It can be passed directly as the first
#'   argument of \code{\link{loa_fit}} or \code{\link{loa_compare}}.
#'
#' @examples
#' set.seed(1)
#' truth <- rlnorm(30, 3.5, 0.5)
#' dat <- expand.grid(id = 1:30, reader = 1:3, session = 1:2)
#' dat$value <- (sqrt(truth[dat$id]) + rnorm(nrow(dat), 0, 0.3))^2
#' p_intra <- loa_pairs(dat$value, dat$id, dat$reader, dat$session,
#'                      type = "intra")   # 30 x 3 = 90 pairs
#' p_inter <- loa_pairs(dat$value, dat$id, dat$reader, dat$session,
#'                      type = "inter")   # 30 x 3 x 2 = 180 pairs
#' loa_compare(p_intra)
#' @export
loa_pairs <- function(value, id, reader, session = NULL,
                      type = c("intra", "inter"),
                      replicates = c("matched", "mean", "first")) {
  type <- match.arg(type)
  replicates <- match.arg(replicates)
  n <- length(value)
  if (length(id) != n || length(reader) != n)
    stop("`value`, `id` and `reader` must have the same length.")
  if (is.null(session)) {
    if (type == "intra")
      stop("`session` is required for type = 'intra'.")
    session <- rep(1L, n)
  }
  if (length(session) != n)
    stop("`session` must have the same length as `value`.")
  ok <- is.finite(value) & !is.na(id) & !is.na(reader) & !is.na(session)
  d <- data.frame(value = value[ok], id = id[ok], reader = reader[ok],
                  session = session[ok], stringsAsFactors = FALSE)
  d <- d[order(d$id, d$reader, d$session), ]

  pairs_of <- function(k) {
    if (k < 2L) return(matrix(integer(0), nrow = 2))
    utils::combn(k, 2)
  }

  if (type == "intra") {
    out <- list()
    for (g in split(d, list(d$id, d$reader), drop = TRUE)) {
      cb <- pairs_of(nrow(g))
      if (!ncol(cb)) next
      out[[length(out) + 1L]] <- data.frame(
        x = g$value[cb[1, ]], y = g$value[cb[2, ]],
        id = g$id[1], reader = g$reader[1],
        session_x = g$session[cb[1, ]], session_y = g$session[cb[2, ]],
        stringsAsFactors = FALSE)
    }
  } else {
    if (replicates == "mean") {
      d <- stats::aggregate(value ~ id + reader, data = d, FUN = mean)
      d$session <- "mean"
    } else if (replicates == "first") {
      d <- d[!duplicated(d[c("id", "reader")]), ]
    }
    out <- list()
    for (g in split(d, list(d$id, d$session), drop = TRUE)) {
      cb <- pairs_of(nrow(g))
      if (!ncol(cb)) next
      out[[length(out) + 1L]] <- data.frame(
        x = g$value[cb[1, ]], y = g$value[cb[2, ]],
        id = g$id[1], reader_x = g$reader[cb[1, ]],
        reader_y = g$reader[cb[2, ]], session = g$session[1],
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) stop("No pairs could be formed from the data.")
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "type") <- type
  attr(res, "replicates") <- if (type == "inter") replicates else NA
  class(res) <- c("loa_pairs", "data.frame")
  res
}

#' @export
print.loa_pairs <- function(x, ...) {
  cat(sprintf("%s-reader difference pairs (flexloa): %d pairs, %d units\n",
              if (attr(x, "type") == "intra") "Intra" else "Inter",
              nrow(x), length(unique(x$id))))
  if (attr(x, "type") == "inter")
    cat(sprintf("  replicates: %s\n", attr(x, "replicates")))
  print(utils::head(as.data.frame(x), ...), row.names = FALSE)
  if (nrow(x) > 6L) cat("  ...\n")
  invisible(x)
}

# Internal: accept either (x, y) vectors or a loa_pairs object as `x`.
.resolve_xy <- function(x, y) {
  if (inherits(x, "loa_pairs")) {
    if (!missing(y) && !is.null(y))
      warning("`y` ignored because `x` is a loa_pairs object.")
    return(list(x = x$x, y = x$y))
  }
  if (missing(y) || is.null(y))
    stop("`y` must be supplied unless `x` is a loa_pairs object.")
  list(x = x, y = y)
}
