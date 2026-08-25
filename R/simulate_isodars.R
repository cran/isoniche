#' Simulate abundances across three habitats
#'
#' Simulates abundance data for three habitats across \code{n} sites using
#' sequential linear relationships (isodars) with Gaussian noise:
#' \itemize{
#' \item hab2 is generated from hab1 using \code{slope1} and \code{int1}
#' \item hab3 is generated from hab2 using \code{slope2} and \code{int2}
#' }
#'
#' The output includes simple per-site niche metrics (CV and Gini) computed
#' from the three habitat abundances.
#'
#' @param slope1 Numeric. Isodar slope for \code{hab2 ~ hab1}.
#' @param slope2 Numeric. Isodar slope for \code{hab3 ~ hab2}.
#' @param int1 Numeric. Isodar intercept for \code{hab2 ~ hab1}.
#' @param int2 Numeric. Isodar intercept for \code{hab3 ~ hab2}.
#' @param noise Numeric (>= 0). Standard deviation of Gaussian noise.
#' @param hab1_max Integer. maximum abundance value in habitat 1 (other maxima calculated from hab1_max).
#' Default is derived to generally suit the isodar parameters, but is arbitrary.
#' @param n Integer (> 0). Number of simulated sites.
#'
#' @return A tibble with one row per simulated site containing:
#' \itemize{
#' \item hab1, hab2, hab3: simulated abundances (nonnegative integers)
#' \item total_abundance: hab1 + hab2 + hab3
#' \item mean_abundance: total_abundance divided  by 3
#' \item sd: scaled standard deviation of habitat abundances
#' \item cv: coefficient of variation (NA if mean_abundance == 0)
#' \item gini: Gini coefficient of habitat abundances
#' }
#'
#' @examples
#' set.seed(1)
#' simulate_isodars(1, 1.1, 1, 0, noise = 2, n = 10)
#' simulate_isodars(1.5, 2, 5, 2, noise = 2, n = 10)
#'
#' @export
simulate_isodars <- function(slope1, slope2, int1, int2, noise = 1, n = 30,
  hab1_max = 100 / (slope1^1.2) - int1) {
  stopifnot(
 is.numeric(slope1),length(slope1) == 1,
 is.numeric(slope2),length(slope2) == 1,
 is.numeric(int1),  length(int1) == 1,
 is.numeric(int2),  length(int2) == 1,
 is.numeric(noise), length(noise) == 1, noise >= 0,
 is.numeric(n),  length(n) == 1, n > 0,
 is.numeric(hab1_max), length(hab1_max) == 1, hab1_max > 0
  )

  simulation <- tibble::tibble(
 hab1 = rep(NA_real_, n),
 hab2 = rep(NA_real_, n),
 hab3 = rep(NA_real_, n)
  )

  if (!is.finite(hab1_max) || hab1_max <= 0) {
 stop("Upper bound for hab1 is non-positive or non-finite; try different values.")
  }

  simulation$hab1 <- stats::runif(n, min = 0, max = hab1_max)
  simulation$hab2 <- slope1 * simulation$hab1 + int1 + stats::rnorm(n, mean = 0, sd = noise)
  simulation$hab3 <- slope2 * simulation$hab2 + int2 + stats::rnorm(n, mean = 0, sd = noise)

  simulation <- round(simulation, 0)
  simulation[simulation < 0] <- 0

  h1 <- simulation$hab1
  h2 <- simulation$hab2
  h3 <- simulation$hab3

  total_abundance <- h1 + h2 + h3
  mean_abundance  <- total_abundance / 3

  sd_raw <- vapply(seq_len(n), function(i) stats::sd(c(h1[i], h2[i], h3[i])), numeric(1))
  sd_scaled <- sqrt(2 / 3) * sd_raw

  cv <- ifelse(mean_abundance > 0, sd_scaled / mean_abundance, NA_real_)
  gini <- function(x) {
 x <- as.numeric(x)
 x <- x[is.finite(x)]
 if (!length(x)) return(NA_real_)
 x <- sort(x)
 n <- length(x)
 (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
  }
  gini <- vapply(seq_len(n), function(i) gini(c(h1[i], h2[i], h3[i])), numeric(1))

  dplyr::bind_cols(
 simulation,
 tibble::tibble(
total_abundance = total_abundance,
mean_abundance  = mean_abundance,
sd  = sd_scaled,
cv  = cv,
gini= gini
 )
  )
}
