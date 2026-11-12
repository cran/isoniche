#' Niche breadth indices from pairwise isodar parameters
#'
#' Summarises the intercepts and slopes of the pairwise isodars directly
#' (Appendix S1 of Granot, Dubiner & Belmaker), as an alternative to the
#' isodar-adjusted inequality ([iai()]), and checks that the isodars are
#' consistent with each other. A perfect generalist has intercepts of 0 and
#' slopes of 1 between all habitat pairs.
#'
#' @details
#' **Orientation invariance.** An isodar between habitats \eqn{i} and
#' \eqn{j} can be written with either habitat on the x-axis:
#' \eqn{N_j = a + b N_i} or, equivalently, \eqn{N_i = -a/b + (1/b) N_j}. With
#' `scale = "log"` (default) every index is the same whichever way each isodar
#' is written:
#' * slope deviations are \eqn{(\log b)^2}{(log b)^2}, since
#'   \eqn{\log(1/b) = -\log b}{log(1/b) = -log b};
#' * intercept deviations use the perpendicular distance of the isodar from
#'   the origin, \eqn{d = |a| / \sqrt{1 + b^2}}{d = |a| / sqrt(1 + b^2)},
#'   which is the same for both ways of writing the line.
#'
#' The variance-based indices are
#' \deqn{1 / \sqrt{\sum w (\log b)^2 / \sum w}}{1 / sqrt(sum(w * log(b)^2) / sum(w))}
#' for slopes and
#' \deqn{1 / \sqrt{\sum w d^2 / \sum w}}{1 / sqrt(sum(w * d^2) / sum(w))}
#' for intercepts, where \eqn{w} are the isodar weights. They measure
#' deviations from the values that indicate no habitat preference
#' (\eqn{b = 1}, \eqn{d = 0}). Larger values indicate broader niches; extreme
#' specialists approach 0, and a perfect generalist gives `Inf`. The mean
#' indices are the weighted means of \eqn{|\log b|}{|log b|} and of \eqn{d}.
#'
#' With `weights = "1/var"` and `scale = "log"`, each isodar is weighted by
#' \eqn{1/sd^2}, the inverse variance of its perpendicular residuals, which is
#' also orientation-invariant.
#'
#' **Non-positive slopes** (for example, a habitat that loses individuals as
#' another gains them) have no logarithm. Such isodars are left out of the
#' slope indices, with a warning naming them, but still used in the intercept
#' indices.
#'
#' **Linear scale.** `scale = "linear"` gives the earlier indices: each isodar
#' is first expressed so that its intercept is non-negative (habitats swapped
#' if needed), then the mean intercept and mean slope, and
#' \eqn{1 / \sqrt{\sum w (b - 1)^2 / \sum w}}{1 / sqrt(sum(w * (b - 1)^2) / sum(w))}
#' and \eqn{1 / \sqrt{\sum w a^2 / \sum w}}{1 / sqrt(sum(w * a^2) / sum(w))},
#' with `"1/var"` weights \eqn{1/(sd^2(1+b^2))}. These depend on how each
#' isodar is oriented, and slopes of 2 and 0.5 (preferences of equal strength
#' in opposite directions) count as different deviations.
#'
#' **Consistency check.** If the isodars describe one consistent pattern of
#' habitat selection, the slopes around any cycle of three habitats multiply
#' to 1: \eqn{b_{AB} b_{BC} b_{CA} = 1}, where \eqn{b_{AB}} is the slope of
#' habitat \eqn{B} on habitat \eqn{A}. Each cycle's product is returned in
#' `cycles`, and the result is flagged as inconsistent when
#' \eqn{|\log(product)|}{|log(product)|} exceeds `cycle_tolerance` for any
#' cycle (or the product is not positive).
#'
#' **Units.** The intercept indices are in units of abundance per sampling
#' unit (see `reference_effort` in [fit_isodars()]), so they depend on survey
#' effort and are comparable among species only when surveys are
#' standardised. The slope indices are unitless.
#'
#' @param isodars Pairwise isodars from `fit_isodars(dim = "pairwise")`, or a
#'   data frame with columns `habitat_x`, `habitat_y`, `intercept`, `slope`
#'   (and `sd`, the SD of perpendicular residuals, for `weights = "1/var"`).
#' @param weights `NULL` (equal weights), `"1/var"` (inverse variance of each
#'   isodar; see Details), or a numeric vector with one non-negative weight per
#'   isodar.
#' @param scale `"log"` (default; orientation-invariant indices) or
#'   `"linear"` (the earlier indices); see Details.
#' @param cycle_tolerance Largest acceptable \eqn{|\log(product)|}{|log(product)|}
#'   of slopes around a cycle of three habitats (default `log(1.5)`).
#' @param log_slope Deprecated: use `scale`. `TRUE` is `scale = "log"` and
#'   `FALSE` is `scale = "linear"`.
#'
#' @return An object of class `pairwise_niche`, a list with
#'   * `indices`: a one-row tibble. For `scale = "log"`: `mean_abs_log_slope`,
#'     `mean_distance`, `slope_niche_breadth`, `intercept_niche_breadth`,
#'     `n_isodars`, `n_slope` (isodars used for the slope indices) and `note`.
#'     For `scale = "linear"`: `mean_intercept`, `mean_slope`,
#'     `slope_niche_breadth`, `intercept_niche_breadth`, `n_isodars` and
#'     `note`.
#'   * `cycles`: a tibble with one row per cycle of three habitats with
#'     isodars between all three pairs: `habitat_1`, `habitat_2`, `habitat_3`,
#'     `product` (\eqn{b_{12} b_{23} b_{31}}), `log_product` and
#'     `inconsistent`.
#'   * `inconsistent`: `TRUE` if any cycle is inconsistent.
#'   * `isodars`: the isodars used, with their `weight`, and for
#'     `scale = "log"` the `distance` \eqn{d} and `log_slope`; for
#'     `scale = "linear"` they are oriented to non-negative intercepts.
#'   * `scale`, `cycle_tolerance`.
#'
#' @seealso [fit_isodars()], [iai()]
#' @examples
#' set.seed(1)
#' sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
#' pw <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
#' pn <- pairwise_niche(pw, weights = "1/var")
#' pn
#' pn$cycles
#' pairwise_niche(pw, scale = "linear")$indices
#' @export
pairwise_niche <- function(isodars, weights = NULL, scale = c("log", "linear"),
                           cycle_tolerance = log(1.5), log_slope = NULL) {
  if (!is.data.frame(isodars) ||
      !all(c("habitat_x", "habitat_y", "intercept", "slope") %in% names(isodars))) {
    stop("'isodars' must be pairwise isodars from fit_isodars().", call. = FALSE)
  }
  if (inherits(isodars, "isodar_line")) {
    stop("pairwise_niche() needs pairwise isodars (dim = \"pairwise\").",
         call. = FALSE)
  }
  if (nrow(isodars) == 0) {
    stop("'isodars' has no rows (fewer than two habitats used?).", call. = FALSE)
  }
  if (!is.null(log_slope)) {
    if (!missing(scale)) {
      stop("Give 'scale' or the deprecated 'log_slope', not both.", call. = FALSE)
    }
    warning("'log_slope' is deprecated; use scale = \"log\" or \"linear\".",
            call. = FALSE)
    scale <- if (isTRUE(as.logical(log_slope))) "log" else "linear"
  }
  scale <- match.arg(scale)
  if (!is.numeric(cycle_tolerance) || length(cycle_tolerance) != 1 ||
      !is.finite(cycle_tolerance) || cycle_tolerance < 0) {
    stop("'cycle_tolerance' must be a single non-negative number.", call. = FALSE)
  }
  iso <- if (scale == "linear") orient_isodars(isodars) else
    tibble::as_tibble(as.data.frame(isodars, stringsAsFactors = FALSE))
  a <- iso$intercept
  b <- iso$slope
  K <- nrow(iso)
  w <- if (is.null(weights)) {
    rep(1, K)
  } else if (identical(weights, "1/var")) {
    if (!"sd" %in% names(iso)) {
      stop("weights = \"1/var\" needs an 'sd' column in 'isodars'.", call. = FALSE)
    }
    if (scale == "log") 1 / pmax(iso$sd^2, 1e-12) else
      1 / pmax(iso$sd^2 * (1 + b^2), 1e-12)
  } else if (is.numeric(weights) && length(weights) == K &&
             all(is.finite(weights)) && all(weights >= 0) && sum(weights) > 0) {
    weights
  } else {
    stop("'weights' must be NULL, \"1/var\" or one non-negative number per isodar.",
         call. = FALSE)
  }
  iso$weight <- w
  wmean <- function(x, use = rep(TRUE, K)) sum(w[use] * x[use]) / sum(w[use])

  if (scale == "log") {
    d <- abs(a) / sqrt(1 + b^2)
    pos <- is.finite(b) & b > 0
    note <- ""
    if (any(!pos)) {
      pairs <- paste0(iso$habitat_x[!pos], "-", iso$habitat_y[!pos])
      note <- paste0("non-positive slope(s) left out of the slope indices: ",
                     paste(pairs, collapse = ", "))
      warning("Isodar(s) with non-positive slopes have no logarithm and are left ",
              "out of the slope indices: ", paste(pairs, collapse = ", "), ".",
              call. = FALSE)
    }
    lb <- ifelse(pos, log(pmax(b, 1e-300)), NA_real_)
    indices <- tibble::tibble(
      mean_abs_log_slope = if (any(pos)) wmean(abs(lb), pos) else NA_real_,
      mean_distance = wmean(d),
      slope_niche_breadth = if (any(pos)) 1 / sqrt(wmean(lb^2, pos)) else NA_real_,
      intercept_niche_breadth = 1 / sqrt(wmean(d^2)),
      n_isodars = K,
      n_slope = sum(pos),
      note = note
    )
    iso$distance <- d
    iso$log_slope <- lb
  } else {
    indices <- tibble::tibble(
      mean_intercept = wmean(a),
      mean_slope = wmean(b),
      slope_niche_breadth = 1 / sqrt(wmean((b - 1)^2)),
      intercept_niche_breadth = 1 / sqrt(wmean(a^2)),
      n_isodars = K,
      note = ""
    )
  }
  cycles <- isodar_cycles(isodars, cycle_tolerance)
  out <- list(indices = indices, cycles = cycles,
              inconsistent = any(cycles$inconsistent), isodars = iso,
              scale = scale, cycle_tolerance = cycle_tolerance)
  class(out) <- "pairwise_niche"
  out
}

# Product of the slopes around every cycle of three habitats
isodar_cycles <- function(isodars, tolerance) {
  none <- tibble::tibble(habitat_1 = character(0), habitat_2 = character(0),
                         habitat_3 = character(0), product = numeric(0),
                         log_product = numeric(0), inconsistent = logical(0))
  # sorted, so that cycles do not depend on how the isodars are written
  habs <- sort(unique(c(isodars$habitat_x, isodars$habitat_y)))
  H <- length(habs)
  if (H < 3) return(none)
  # S[i, j]: slope of habitat j on habitat i (TLS lines are symmetric, so the
  # reverse orientation has slope 1 / b)
  S <- matrix(NA_real_, H, H, dimnames = list(habs, habs))
  for (k in seq_len(nrow(isodars))) {
    hx <- isodars$habitat_x[k]
    hy <- isodars$habitat_y[k]
    S[hx, hy] <- isodars$slope[k]
    S[hy, hx] <- if (isodars$slope[k] != 0) 1 / isodars$slope[k] else NA_real_
  }
  out <- list()
  fitted <- function(p, q) !is.na(S[p, q]) || !is.na(S[q, p])
  for (i in 1:(H - 2)) for (j in (i + 1):(H - 1)) for (k in (j + 1):H) {
    if (!fitted(i, j) || !fitted(j, k) || !fitted(k, i)) next
    # NA when a zero slope makes the reverse orientation undefined
    p <- S[i, j] * S[j, k] * S[k, i]
    lp <- if (is.finite(p) && p > 0) log(p) else NA_real_
    out[[length(out) + 1]] <- data.frame(
      habitat_1 = habs[i], habitat_2 = habs[j], habitat_3 = habs[k],
      product = p, log_product = lp,
      inconsistent = is.na(lp) || abs(lp) > tolerance,
      stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(none)
  tibble::as_tibble(do.call(rbind, out))
}

#' @export
print.pairwise_niche <- function(x, ...) {
  cat("Niche breadth from pairwise isodars (", x$scale, " scale)\n", sep = "")
  print(x$indices, ...)
  if (nrow(x$cycles) > 0) {
    n_bad <- sum(x$cycles$inconsistent)
    cat(nrow(x$cycles), " cycle(s) of three habitats; ", n_bad,
        " inconsistent (|log(product of slopes)| > ",
        format(signif(x$cycle_tolerance, 3)), "); see $cycles\n", sep = "")
  }
  invisible(x)
}

# Express each isodar with a non-negative intercept (swap habitats if needed)
orient_isodars <- function(isodars) {
  iso <- as.data.frame(isodars, stringsAsFactors = FALSE)
  if (is.null(iso$flipped)) iso$flipped <- FALSE
  neg <- is.finite(iso$intercept) & iso$intercept < 0
  if (any(neg & !(iso$slope > 0))) {
    stop("An isodar with a negative intercept and a non-positive slope cannot ",
         "be expressed with a non-negative intercept.", call. = FALSE)
  }
  if (any(neg)) {
    hx <- iso$habitat_x[neg]
    iso$habitat_x[neg] <- iso$habitat_y[neg]
    iso$habitat_y[neg] <- hx
    iso$intercept[neg] <- -iso$intercept[neg] / iso$slope[neg]
    iso$slope[neg] <- 1 / iso$slope[neg]
    iso$flipped[neg] <- !iso$flipped[neg]
  }
  tibble::as_tibble(iso)
}
