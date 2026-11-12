#' Simulate abundances in habitats linked by isodars
#'
#' Simulates abundances at `n` sites along a chain of isodars: habitat 1 is
#' drawn uniformly between 0 and `hab1_max`, and each following habitat
#' \eqn{k + 1} = `intercepts[k]` + `slopes[k]` x habitat \eqn{k} + noise.
#' Values are rounded to integers and negative values set to 0. The default
#' is three habitats, set by `slope1`, `slope2`, `int1` and `int2`; the
#' sequence of random numbers is the same as in version 1.0.0.
#'
#' @param slope1,int1 Isodar of habitat 2 on habitat 1.
#' @param slope2,int2 Isodar of habitat 3 on habitat 2.
#' @param noise Standard deviation of the Gaussian noise.
#' @param n Number of sites.
#' @param hab1_max Upper bound of abundance in habitat 1. The default,
#'   `100 / slopes[1]^1.2 - intercepts[1]`, keeps total abundances in a
#'   comparable range across parameter values.
#' @param slopes,intercepts Any number of habitats: vectors of length
#'   \eqn{H - 1} giving the isodar of each habitat on the previous one
#'   (default: `c(slope1, slope2)` and `c(int1, int2)`).
#' @param p_missing Probability that a habitat is missing (not sampled) at a
#'   site; at least `min_habitats` habitats are kept at every site.
#' @param min_habitats Minimum number of habitats recorded per site.
#' @param format `"wide"` (one column per habitat, `NA` where missing) or
#'   `"long"` (columns `site`, `habitat`, `abundance`; missing habitats have
#'   no row).
#' @return Wide: a tibble with columns `hab1`, ..., `habH` and
#'   `total_abundance`. Long: a tibble with `site`, `habitat` and `abundance`.
#' @examples
#' set.seed(1)
#' simulate_isodars(1, 1.1, 1, 0, noise = 2, n = 10)
#' simulate_isodars(slopes = c(1.5, 1, 0.8), intercepts = c(5, 0, 2), n = 5,
#'                  p_missing = 0.3)
#' @export
simulate_isodars <- function(slope1 = 1, slope2 = 1, int1 = 0, int2 = 0,
                             noise = 1, n = 30,
                             hab1_max = 100 / (slopes[1]^1.2) - intercepts[1],
                             slopes = c(slope1, slope2),
                             intercepts = c(int1, int2), p_missing = 0,
                             min_habitats = 2, format = c("wide", "long")) {
  format <- match.arg(format)
  stopifnot(is.numeric(noise), noise >= 0, n > 0)
  if (!is.numeric(slopes) || !is.numeric(intercepts) || length(slopes) < 1 ||
      length(slopes) != length(intercepts)) {
    stop("'slopes' and 'intercepts' must be numeric vectors of the same length.",
         call. = FALSE)
  }
  H <- length(slopes) + 1
  if (!is.numeric(p_missing) || length(p_missing) != 1 || p_missing < 0 ||
      p_missing >= 1) {
    stop("'p_missing' must be a number in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(min_habitats) || min_habitats < 1 || min_habitats > H) {
    stop("'min_habitats' must be between 1 and the number of habitats.", call. = FALSE)
  }
  if (!is.finite(hab1_max) || hab1_max <= 0) {
    stop("'hab1_max' must be positive; try different parameter values.", call. = FALSE)
  }
  X <- matrix(0, n, H, dimnames = list(NULL, paste0("hab", seq_len(H))))
  X[, 1] <- stats::runif(n, 0, hab1_max)
  for (k in seq_len(H - 1)) {
    X[, k + 1] <- slopes[k] * X[, k] + intercepts[k] + stats::rnorm(n, 0, noise)
  }
  X <- pmax(round(X), 0)
  if (p_missing > 0) {
    miss <- matrix(stats::runif(n * H) < p_missing, n, H)
    for (i in seq_len(n)) {
      while (sum(!miss[i, ]) < min_habitats) {
        miss[i, which(miss[i, ])[sample.int(sum(miss[i, ]), 1)]] <- FALSE
      }
    }
    X[miss] <- NA
  }
  if (format == "long") {
    ok <- !is.na(X)
    return(tibble::tibble(
      site = paste0("site", row(X)[ok]), habitat = colnames(X)[col(X)[ok]],
      abundance = X[ok]
    )[order(row(X)[ok], col(X)[ok]), ])
  }
  out <- tibble::as_tibble(as.data.frame(X))
  out$total_abundance <- rowSums(X, na.rm = TRUE)
  out
}
