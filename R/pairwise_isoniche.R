#' Pairwise isodar-based niche breadth indices
#'
#' Computes niche breadth indices from a set of pairwise isodars (slopes and
#' intercepts) in the format of the output table of \code{fit_isodar()}.
#' Two methods are available: a weighted mean of isodar components and a
#' weighted inverse variance of the isodars' divergence from neutrality.
#'
#' @param data A data frame returned by \code{fit_isodar()}.
#' @param method Character. One of \code{"mean"} or \code{"inverse_variance"}.
#' @param weights Weighting scheme for pairwise isodars. One of:
#' \itemize{
#'   \item \code{NULL}: uniform weights
#'   \item numeric vector of length \code{nrow(data)}
#'   \item \code{"1/var"}: weights proportional to 1/sd^2 (sd from \code{fit_isodar})
#'   \item \code{"sig"}: weights based on significance of intercepts, using \code{sig_weights}
#' }
#' @param alpha Numeric in (0, 1). Significance threshold used when \code{weights = "sig"}.
#' @param sig_weights Numeric vector of length 2 giving weights for significant and non-significant
#'   isodars when \code{weights = "sig"}. Can be named \code{c(sig = ..., nonsig = ...)} or unnamed
#'   \code{c(sig, nonsig)}. Default is \code{c(sig = 1, nonsig = 0)}.
#'
#' @return A one-row tibble with columns \code{intercept}, \code{slope}, and \code{method}.
#'
#' @examples
#' set.seed(1)
#' isod <- simulate_isodars(1, 2, 5, 1, noise = 2, n = 10)
#' INB <- fit_isodar(isod, n_habitats = 3)
#' pairwise_isoniche(INB, method = "inverse_variance", weights = "1/var")
#' pairwise_isoniche(INB, method = "inverse_variance", weights = "sig",
#'              sig_weights = c(sig = 1, nonsig = 0.25), alpha = 0.1)
#'
#' @export
pairwise_isoniche <- function(
    data,
    method = c("inverse_variance", "mean"),
    weights = NULL,
    alpha = 0.1,
    sig_weights = c(sig = 1, nonsig = 0)
) {
  method <- match.arg(method)
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) >= 2)

  required <- c("slope", "intercept")
  if (!all(required %in% names(data))) {
    stop("Pairwise methods require columns: slope, intercept")
  }

  stopifnot(is.numeric(alpha), length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)

  # Validate sig_weights
  stopifnot(is.numeric(sig_weights), length(sig_weights) == 2L, all(is.finite(sig_weights)))
  if (is.null(names(sig_weights))) {
    w_sig <- sig_weights[1L]
    w_nonsig <- sig_weights[2L]
  } else {
    if (!all(c("sig", "nonsig") %in% names(sig_weights))) {
      stop("sig_weights must be length-2 with names c('sig','nonsig') or unnamed c(sig, nonsig).")
    }
    w_sig <- unname(sig_weights["sig"])
    w_nonsig <- unname(sig_weights["nonsig"])
  }
  if (w_sig < 0 || w_nonsig < 0) stop("sig_weights must be non-negative.")

  n <- nrow(data)

  # construct raw weights
  if (is.null(weights)) {
    w <- rep(1, n)

  } else if (is.character(weights)) {
    if (length(weights) != 1L) stop("weights must be a single string when character.")

    if (identical(weights, "1/var")) {
      if (!"sd" %in% names(data)) stop("weights = '1/var' requires column 'sd'.")
      eps <- 1e-12
      w <- 1 / pmax(data$sd, eps)^2

    } else if (identical(weights, "sig")) {
      if (!"p_intercept" %in% names(data)) {
        stop("weights = 'sig' requires column 'p_intercept'.")
      }
      is_sig <- is.finite(data$p_intercept) & (data$p_intercept < alpha)
      w <- ifelse(is_sig, w_sig, w_nonsig)

      if (all(w <= 0, na.rm = TRUE)) {
        stop("All weights are zero under weights = 'sig' (check sig_weights and alpha).")
      }

    } else {
      stop("Unknown weights specification: ", weights)
    }

  } else {
    stopifnot(is.numeric(weights), length(weights) == n)
    w <- weights
  }

  # filter invalid rows / zero weights
  ok <- is.finite(w) & (w > 0) &
    is.finite(data$slope) & is.finite(data$intercept)

  if (!any(ok)) {
    stop("No usable isodars after filtering missing values / non-positive weights.")
  }

  w <- w[ok]
  slope <- data$slope[ok]
  intercept <- data$intercept[ok]

  s <- sum(w)
  if (!is.finite(s) || s <= 0) stop("Weights sum to 0 after filtering.")
  w <- w / s

  # compute indices
  if (method == "mean") {
    INB_slope     <- sum(w * slope)
    INB_intercept <- sum(w * intercept)
  } else { # inverse_variance
    INB_slope     <- 1 / sqrt(sum(w * (slope - 1)^2))
    INB_intercept <- 1 / sqrt(sum(w * (intercept)^2))
  }

  tibble::tibble(
    intercept = INB_intercept,
    slope = INB_slope,
    method = method
  )
}
