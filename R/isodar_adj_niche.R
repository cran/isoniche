#' Compute the isodar-adjusted inequality index for one or more total abundances
#'
#' Reconstructs habitat abundances for each requested total abundance by solving a
#' weighted least-squares system defined by the pairwise isodars:
#' \deqn{habitat_y - slope \cdot habitat_x = intercept}
#'
#' The reconstruction is constrained to be nonnegative and to sum to the requested
#' total abundance (active-set heuristic). Niche breadth is then computed
#' from the reconstructed habitat vector; currently defined as:
#' \deqn{niche breadth = 1 - Gini(x)} of isodar-reconstructed abundances at a given N.
#'
#' @param isodars A data frame returned by \code{fit_isodar()}
#' @param abundances Numeric vector of total abundances. If NULL, the minimal total
#'   abundance yielding occupancy in all habitats ("baseline") is chosen when possible; otherwise
#'   a default exploratory sequence is used
#' @param weights Weighting scheme. One of:
#' \itemize{
#'   \item \code{NULL}: equal weights
#'   \item numeric vector of length \code{nrow(isodars)}
#'   \item \code{"1/var"}: weights proportional to 1/sd^2
#'   \item \code{"sig"}: weights based on significance of intercepts, using \code{sig_weights}
#' }
#' @param method Character. Currently only \code{"gini"}
#' @param plot Logical. If TRUE, plots adjusted niche breadth versus total abundance
#' @param max_search Integer. Maximum total abundance to search when \code{abundances = NULL},, in order to find baseline n
#' @param alpha Numeric in (0, 1). Significance threshold used when \code{weights = "sig"}
#' @param sig_weights Numeric vector of length 2 giving weights for significant and non-significant
#'   isodars when \code{weights = "sig"}. Can be named \code{c(sig = , nonsig = )} or unnamed
#'   \code{c(sig, nonsig)}. Default is \code{c(sig = 1, nonsig = 0)}.
#'
#' @return A tibble with one row per total abundance containing:
#' \itemize{
#'   \item \code{total_abundance}
#'   \item niche_breadth
#'   \item one column per habitat with reconstructed abundance
#' }
#'
#' @examples
#' set.seed(1)
#' isod <- simulate_isodars(1, 2, 5, 1, noise = 2, n = 10)
#' INB <- fit_isodar(isod, n_habitats = 3)
#' isodar_adj_niche(INB, max_search = 100) # automatically checks only baseline abundance
#' isodar_adj_niche(INB, abundances = c(15, 30, 45), max_search = 100) # set specific values
#'
#' @export
isodar_adj_niche <- function(
    isodars,
    abundances = NULL,
    weights = NULL,
    alpha = 0.1,
    sig_weights = c(sig = 1, nonsig = 0),
    method = c("gini"),
    plot = TRUE,
    max_search = 10000
) {
  method <- match.arg(method)

  stopifnot(is.data.frame(isodars))
  req <- c("habitat_x", "habitat_y", "slope", "intercept")
  missing <- setdiff(req, names(isodars))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

  habitats <- sort(unique(c(isodars$habitat_x, isodars$habitat_y)))
  H <- length(habitats)
  if (H < 2) stop("Need at least 2 habitats.")

  stopifnot(
    is.numeric(alpha), length(alpha) == 1L, is.finite(alpha),
    alpha > 0, alpha < 1
  )

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

  # Weight vector for abundance reconstruction
  w <- rep(1, nrow(isodars))

  if (is.null(weights)) {
    # equal weights (already set)
  } else if (is.character(weights)) {
    if (length(weights) != 1L) stop("weights character option must be length 1.")

    if (weights == "1/var") {
      if (!"sd" %in% names(isodars)) stop("weights = '1/var' requires column 'sd'.")
      w <- 1 / pmax(isodars$sd, 1e-12)^2

    } else if (weights == "sig") {
      if (!"p_intercept" %in% names(isodars)) {
        stop("weights = 'sig' requires column 'p_intercept'.")
      }

      is_sig <- is.finite(isodars$p_intercept) & (isodars$p_intercept < alpha)
      w <- ifelse(is_sig, w_sig, w_nonsig)

      if (all(w <= 0, na.rm = TRUE)) {
        stop("All weights are zero under weights = 'sig' (check sig_weights and alpha).")
      }

    } else {
      stop("Unknown weights specification: ", weights,
           ". Use NULL, numeric, '1/var', or 'sig'.")
    }

  } else if (is.numeric(weights)) {
    if (length(weights) != nrow(isodars)) {
      stop("Numeric weights must have length nrow(isodars).")
    }
    if (any(!is.finite(weights))) stop("Numeric weights must be finite.")
    if (any(weights < 0)) stop("Numeric weights must be non-negative.")
    w <- weights

  } else {
    stop("weights must be NULL, a numeric vector, or one of '1/var', 'sig'.")
  }

  ok_row <- is.finite(isodars$slope) & is.finite(isodars$intercept) & is.finite(w) & (w > 0)
  is2 <- isodars[ok_row, , drop = FALSE]
  w2  <- w[ok_row]

  if (nrow(is2) == 0) {
    warning("No usable isodars after filtering missing values / zero weights; returning NA tibble.")
    totals <- if (is.null(abundances)) NA_real_ else as.numeric(abundances)
    extra <- as.list(rep(NA_real_, H))
    names(extra) <- habitats
    return(dplyr::bind_cols(
      tibble::tibble(total_abundance = totals, niche_breadth = NA_real_),
      tibble::as_tibble(extra)
    ))
  }

  # De-duplicate directed pairs: keep most-weighted
  key <- paste(is2$habitat_x, is2$habitat_y, sep = "->")
  if (any(duplicated(key))) {
    o <- order(key, -w2)
    is2 <- is2[o, , drop = FALSE]
    w2  <- w2[o]
    keep <- !duplicated(paste(is2$habitat_x, is2$habitat_y, sep = "->"))
    is2 <- is2[keep, , drop = FALSE]
    w2  <- w2[keep]
  }

  # y - slope * x = intercept
  m <- nrow(is2)
  A <- matrix(0, nrow = m, ncol = H)
  b <- as.numeric(is2$intercept)

  for (i in seq_len(m)) {
    xi <- match(is2$habitat_x[i], habitats)
    yi <- match(is2$habitat_y[i], habitats)
    A[i, yi] <- 1
    A[i, xi] <- -is2$slope[i]
  }

  # Weighted normal equations (using diag(w^2))
  w2sq <- w2^2
  Q   <- t(A) %*% (A * w2sq)
  rhs <- t(A) %*% (b * w2sq)

  # constrained solver
  solve_total <- function(Total) {
    Total <- as.numeric(Total)
    if (!is.finite(Total) || Total < 0) return(rep(NA_real_, H))

    free <- rep(TRUE, H)
    x <- rep(0, H)

    for (iter in seq_len(H + 2L)) {
      idx <- which(free)
      if (length(idx) == 0) break

      Qf <- Q[idx, idx, drop = FALSE]
      rf <- rhs[idx]
      ones <- rep(1, length(idx))

      M <- rbind(
        cbind(Qf, ones),
        c(ones, 0)
      )
      bb <- c(as.numeric(rf), Total)

      sol <- tryCatch(solve(M, bb), error = function(e) NULL)
      if (is.null(sol)) return(rep(NA_real_, H))

      x[idx] <- sol[seq_along(idx)]

      neg <- which(x < -1e-10 & free)
      if (length(neg) == 0) break
      free[neg] <- FALSE
      x[neg] <- 0
    }

    x[x < 0] <- 0
    s <- sum(x)
    if (s > 0) x <- x / s * Total
    x
  }

  # choose baseline abundances (minimal full occupancy)
  if (is.null(abundances)) {
    candidates <- seq_len(max_search)
    occupied <- vapply(candidates, function(T) {
      x <- solve_total(T)
      if (!all(is.finite(x))) return(FALSE)
      if (abs(sum(x) - T) > 1e-6) return(FALSE)
      all(x > 1e-9)
    }, logical(1))

    if (!any(occupied)) {
      warning(
        "No total abundance up to max_search=", max_search,
        " yields occupancy in all habitats. Using a default sequence."
      )
      abundances <- seq(1, 1000, length.out = 30)
    } else {
      abundances <- min(candidates[occupied])
    }
  } else {
    abundances <- as.numeric(abundances)
    if (any(!is.finite(abundances)) || any(abundances < 0)) {
      stop("abundances must be finite and >= 0.")
    }
  }

  # local Gini
  gini_fun <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    if (all(x == 0)) return(NA_real_)
    x <- sort(x)
    n <- length(x)
    (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
  }

  # compute isodar-adjusted inequality
  out <- lapply(abundances, function(Total) {
    x <- solve_total(Total)
    names(x) <- habitats

    iai <- switch(
      method,
      gini = 1 - gini_fun(x)
    )

    tibble::as_tibble(c(
      list(total_abundance = Total, niche_breadth = iai),
      as.list(x)
    ))
  })

  res <- dplyr::bind_rows(out)

  if (isTRUE(plot)) {
    graphics::plot(
      res$total_abundance, res$niche_breadth,
      type = "l",
      xlab = "Total abundance",
      ylab = "Isodar-adjusted inequality",
      ylim = c(0, 1),
      lwd = 2,
      col = "red"
    )
  }

  res
}
