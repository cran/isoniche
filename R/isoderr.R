#' Bootstrap uncertainty of the isodar-adjusted inequality
#'
#' Resamples sites with replacement (keeping all habitats of a site together),
#' refits the isodars and recomputes the IAI for each resample. Absent
#' habitats are identified once, on the full data, and kept fixed in every
#' resample.
#'
#' @inheritParams fit_isodars
#' @inheritParams iai
#' @param n_boot Number of bootstrap replicates.
#' @param output `"SD"` or `"SE"` (both the bootstrap standard deviation, which
#'   is the standard error of the estimate) or `"CI"` (percentile interval).
#' @param conf_level Confidence level for `output = "CI"`.
#' @param seed Optional random-number seed.
#' @return A tibble with `total_abundance`, `niche_breadth` (estimate from the
#'   full data), `uncertainty` (SD/SE) or `lower` and `upper` (CI), and
#'   `n_valid` (bootstrap replicates that produced a value).
#' @examples
#' set.seed(1)
#' sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
#' isoderr(sim, hab_cols = c("hab1", "hab2", "hab3"), abundances = 100,
#'         n_boot = 50, seed = 1)
#' @export
isoderr <- function(data, hab_cols = NULL, n_habitats = NULL, dim = "pairwise",
                    abundances = NULL, n_boot = 1000, output = c("SD", "SE", "CI"),
                    conf_level = 0.95, weights = NULL, flip_intercept = TRUE,
                    min_occupancy = 1, max_search = 10000, seed = NULL,
                    ndim_method = "auto", absent = "auto", absent_max_count = 0,
                    absent_min_sites = 5) {
  output <- match.arg(output)
  if (!is.null(seed)) {
    old_rng <- save_rng()
    on.exit(restore_rng(old_rng), add = TRUE)
    set.seed(seed)
  }
  X <- as.data.frame(habitat_matrix(data, hab_cols, n_habitats))
  habs <- names(X)
  one <- function(d, absent) {
    f <- fit_isodars(d, hab_cols = habs, dim = dim, flip_intercept = flip_intercept,
                     ndim_method = ndim_method, absent = absent,
                     absent_max_count = absent_max_count,
                     absent_min_sites = absent_min_sites)
    r <- iai(f, abundances = abundances, weights = if (dim == "pairwise") weights else NULL,
             min_occupancy = min_occupancy, max_search = max_search)
    attr(r, "absent") <- attr(f, "absent")
    r
  }
  est <- one(X, absent)
  # absent habitats are decided once, on the full data
  absent_set <- attr(est, "absent")
  n_ab <- nrow(est)
  boot <- matrix(NA_real_, n_boot, n_ab)
  for (b in seq_len(n_boot)) {
    d <- X[sample.int(nrow(X), replace = TRUE), , drop = FALSE]
    r <- tryCatch(suppressWarnings(one(d, absent_set)), error = function(e) NULL)
    if (!is.null(r) && nrow(r) == n_ab) boot[b, ] <- r$niche_breadth
  }
  res <- tibble::tibble(total_abundance = est$total_abundance,
                        niche_breadth = est$niche_breadth)
  if (output %in% c("SD", "SE")) {
    res$uncertainty <- apply(boot, 2, stats::sd, na.rm = TRUE)
  } else {
    a <- (1 - conf_level) / 2
    res$lower <- apply(boot, 2, stats::quantile, probs = a, na.rm = TRUE, names = FALSE)
    res$upper <- apply(boot, 2, stats::quantile, probs = 1 - a, na.rm = TRUE, names = FALSE)
  }
  res$n_valid <- colSums(is.finite(boot))
  res
}
