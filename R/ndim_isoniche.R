#' Multidimensional isodar-based niche breadth indices
#'
#' Computes vector- and orthogonal isodar-based niche breadth components
#' from abundance data across all habitats simultaneously, generalizing the
#' slope (density-dependent) and intercept/baseline (density-independent)
#' components of isodars to \eqn{n} dimensions.
#'
#' Both indices return values between 0 and 1, where 1 indicates a perfect
#' generalist and 0 a specialist. Habitat axes can optionally be weighted
#' by inverse residual variance.
#'
#' @param data A data frame of abundance data (one column per habitat).
#' @param n_habitats Integer (>=2). Number of habitat columns to use, starting at column 1.
#' @param weights Logical. If \code{TRUE}, habitats are weighted by inverse residual variance.
#'   If \code{FALSE}, uniform weights are used.
#' @param scale Logical. If \code{TRUE}, results are raised to the power of
#'    \eqn{n\_habitats - 1} to compensate for reduced sensitivity due to added dimensions.
#'
#' @return A tibble with \code{vector}, \code{orthogonal}, and \code{n_habitats}.
#'
#' @examples
#' set.seed(1)
#' isod <- simulate_isodars(1, 2, 5, 1, noise = 2, n = 10)
#'
#' @export
ndim_isoniche <- function(data,
                          n_habitats = ncol(data),
                          scale = TRUE,
                          weights = TRUE) {

  if (scale==TRUE) {
    m = n_habitats - 1
  }
  else{
    m = 1
  }

  stopifnot(is.data.frame(data))
  stopifnot(is.numeric(n_habitats), length(n_habitats) == 1L,
            n_habitats >= 2L, n_habitats <= ncol(data))
  stopifnot(is.logical(weights), length(weights) == 1L)

  eps <- 1e-8

  X <- as.matrix(data[, seq_len(n_habitats), drop = FALSE])
  X <- X[stats::complete.cases(X), , drop = FALSE]

  if (nrow(X) < 2L) {
    stop("Not enough complete observations to calculate niche breadth.")
  }

  n <- ncol(X)

  # mean and centering

  mu <- colMeans(X)
  Xc <- sweep(X, 2L, mu, "-")

  # dominant density-dependent axis

  s <- svd(Xc)
  v <- s$v[, 1L]

  # residuals for vector component

  scores_v <- drop(Xc %*% v)
  residuals_v <- Xc - scores_v %*% t(v)

  # neutral axis

  v0 <- rep(1 / sqrt(n), n)

  # residuals for orthogonal component

  scores_v0 <- drop(Xc %*% v0)
  residuals_v0 <- Xc - scores_v0 %*% t(v0)

  # weight function

  make_weights <- function(residuals_mat, use_weights) {
    if (!use_weights) {
      return(rep(1 / n, n))
    }
    sigma2 <- apply(residuals_mat, 2L, stats::var)
    w_raw <- 1 / (sigma2 + eps)
    w_sum <- sum(w_raw)
    if (!is.finite(w_sum) || w_sum <= 0) {
      stop("Failed to compute weights.")
    }
    w_raw / w_sum
  }

  w_vec <- make_weights(residuals_v, weights)
  w_ort <- make_weights(residuals_v0, weights)

  # -------------------------

  # Density-dependent (vector)

  # -------------------------

  v_abs <- abs(v)
  v_neutral <- rep(1 / sqrt(n), n)

  delta_v <- v_abs / v_neutral - 1
  Q_vec <- sum(w_vec * delta_v^2)
  INB_vector <- 1 / (1 + Q_vec) ^ m

  # -------------------------

  # Density-independent (orthogonal, scale-invariant)

  # -------------------------

  mu_sum <- sum(mu)

  if (mu_sum <= .Machine$double.eps) {
    INB_orthog <- NA_real_
  } else {

    # compositional normalization
    p <- mu / mu_sum
    p0 <- rep(1 / n, n)

    r <- p - p0
    r_norm <- sqrt(sum(r^2))

    if (r_norm < .Machine$double.eps) {
      INB_orthog <- 1
    } else {
      r_abs <- abs(r)
      r_neutral <- rep(r_norm / sqrt(n), n)

      delta_r <- r_abs / r_neutral - 1
      Q_ort <- sum(w_ort * delta_r^2)
      INB_orthog <- 1 / (1 + Q_ort) ^ m
    }

  }

  tibble::tibble(
    vector = INB_vector,
    orthogonal = INB_orthog,
    n_habitats = n_habitats
  )
}
