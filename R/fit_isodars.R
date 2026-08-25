#' Fit pairwise or n-dimensional isodars using total least squares
#'
#' Fits isodar relationships using total least squares (TLS), which minimizes
#' orthogonal distances between observations and the fitted relationship.
#'
#' With \code{dim = "pairwise"} (default), TLS relationships are fitted between
#' every pair of habitats.
#'
#' With \code{dim = "ndim"}, a single multidimensional TLS hyperplane is fitted
#' using all habitats simultaneously.
#'
#' @param data A data frame containing abundance columns for the habitats.
#' @param hab_cols Character vector giving the habitat abundance columns.
#' If \code{NULL}, the first \code{n_habitats} columns are used.
#' @param n_habitats Number of habitats to use when \code{hab_cols = NULL}.
#' @param flip_intercept Logical. If TRUE, pairwise relationships with a
#' negative intercept are reversed so that the reported intercept is positive.
#' @param dim Either \code{"pairwise"} (default) or \code{"ndim"}.
#'
#' @return For \code{dim = "pairwise"}, a tibble with one row for every
#' habitat pair, containing the coefficients of the TLS regression.
#' For \code{dim = "ndim"}, a tibble containing the coefficients
#' of the n-dimensional TLS hyperplane.
#' @examples
#' set.seed(1)
#' populations1 <- simulate_isodars(1, 1.1, 1, 0, noise = 2, n = 10)
#'   fit_isodars(populations1, hab_cols=c("hab1","hab2","hab3"))
#' populations2 <- simulate_isodars(1.5, 2, 5, 2, noise = 2, n = 10)
#'   fit_isodars(populations2, hab_cols=c("hab1","hab2","hab3"))
#'
#' @export
fit_isodars <- function(
 data,
 hab_cols = NULL,
 n_habitats = ncol(data),
 flip_intercept = TRUE,
 dim = "pairwise"
) {

  if (!is.data.frame(data)) {
 stop("'data' must be a data frame.")
  }
  dim <- match.arg(dim, c("pairwise", "ndim"))
 if (is.null(hab_cols)) {
 if (n_habitats < 2 ||
  n_habitats > ncol(data)) {
stop("'n_habitats' must be between 2 and the number of columns in 'data'.")
 }
 habs <- names(data)[seq_len(n_habitats)]
} else {
 if (!is.character(hab_cols) || length(hab_cols) < 2) {
stop("'hab_cols' must contain at least two column names.")
 }
 if (!all(hab_cols %in% names(data))) {
stop("Some values in 'hab_cols' are not columns in 'data'.")
 }
 habs <- hab_cols
 n_habitats <- length(habs)
  }
  if (any(is.na(habs)) || any(habs == "")) {
 stop("All habitat columns must have non-empty column names.")
  }

  # Pairwise TLS

  fit_pair <- function(x, y, name_x, name_y) {

 ok <- is.finite(x) & is.finite(y)
 x <- x[ok]
 y <- y[ok]
 n <- length(x)
 if (n < 3 ||
  sd(x) == 0 ||
  sd(y) == 0 ||
  all(x == 0) ||
  all(y == 0)) {

return(tibble::tibble(
  habitat_x = name_x,
  habitat_y = name_y,
  intercept = NA_real_,
  slope = NA_real_,
  sd = NA_real_,
  flipped = NA
))
 }
 # Center the data
 x_mean <- mean(x)
 y_mean <- mean(y)
 X <- cbind(
x - x_mean,
y - y_mean
 )

 # TLS is obtained from the eigenvector corresponding
 # to the smallest eigenvalue of the covariance matrix.
 S <- crossprod(X) / (n - 1)
 eig <- eigen(S, symmetric = TRUE)
 normal <- eig$vectors[, 2]
 # Equation of line:
 # normal[1] * (x - x_mean) +
 # normal[2] * (y - y_mean) = 0
 # Rearranged as:
 # y = intercept + slope * x
 if (abs(normal[2]) < 1e-12) {
# Vertical TLS relationship cannot be expressed as y = a + bx
return(tibble::tibble(
  habitat_x = name_x,
  habitat_y = name_y,
  intercept = NA_real_,
  slope = NA_real_,
  sd = NA_real_,
  flipped = NA
))
 }

 slope <- -normal[1] / normal[2]
 intercept <- y_mean - slope * x_mean

 # Orthogonal distances to TLS line
 distances <- (
normal[1] * (x - x_mean) +
  normal[2] * (y - y_mean)
 )
 sd_tls <- sqrt(
sum(distances^2) / (n - 2)
 )

 # Flip axes if intercept is negative

 if (isTRUE(flip_intercept) &&
  is.finite(intercept) &&
  intercept < 0) {

flipped_fit <- fit_pair(
  y,
  x,
  name_y,
  name_x
)

return(tibble::tibble(
  habitat_x = flipped_fit$habitat_x,
  habitat_y = flipped_fit$habitat_y,
  intercept = flipped_fit$intercept,
  slope = flipped_fit$slope,
  sd = flipped_fit$sd,
  flipped = TRUE
))
 }

 tibble::tibble(
habitat_x = name_x,
habitat_y = name_y,
intercept = intercept,
slope = slope,
sd = sd_tls,
flipped = FALSE
 )

  }

  # PAIRWISE

  if (dim == "pairwise") {

 results <- list()
 k <- 1
 for (i in seq_len(n_habitats - 1)) {
for (j in (i + 1):n_habitats) {
  results[[k]] <- fit_pair(
 data[[habs[i]]],
 data[[habs[j]]],
 habs[i],
 habs[j]
  )
  k <- k + 1
}
 }
 return(dplyr::bind_rows(results))

  }

  # N-DIMENSIONAL

  X <- as.matrix(data[, habs, drop = FALSE])

  # Keep only complete observations
  X <- X[stats::complete.cases(X), , drop = FALSE]
  n <- nrow(X)
  p <- ncol(X)

  if (n < p + 1) {
 stop(
"Not enough complete observations for multidimensional TLS."
 )
  }

  # Check that all columns are numeric
  if (!all(vapply(
 data[habs],
 is.numeric,
 logical(1)
  ))) {
 stop("All habitat columns must be numeric.")
  }

  # Check for zero variance
  if (any(apply(X, 2, sd) == 0)) {
 stop("At least one habitat has zero variance.")
  }

  # Center the data
  means <- colMeans(X)
  X_centered <- sweep(
 X,
 2,
 means,
 FUN = "-"
  )

  # Covariance matrix
 S <- crossprod(X_centered) / (n - 1)
 # Eigenvector corresponding to the smallest eigenvalue
 eig <- eigen(
 S,
 symmetric = TRUE
  )
 normal <- eig$vectors[, p]
 # Make sign reproducible
 if (normal[1] < 0) {
 normal <- -normal
  }

  # Hyperplane:
  # b1*H1 + b2*H2 + ... + bp*Hp + intercept = 0
  intercept <- -sum(normal * means)
  # Orthogonal distances
  distances <- as.vector(
  X_centered %*% normal
  )

  # Residual SD
 sd_tls <- sqrt(
 sum(distances^2) / (n - p)
  )

  tibble::tibble(
 habitat = habs,
 coefficient = as.numeric(normal),
 intercept = intercept,
 sd = sd_tls
  )
}
