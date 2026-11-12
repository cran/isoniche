#' Isodar-adjusted inequality (IAI)
#'
#' Computes niche breadth from fitted isodars by allocating a given total
#' abundance among habitats according to the isodars and measuring the
#' inequality of that allocation (1 - Gini). Evaluating it at several totals
#' gives a rarefaction-like curve of niche breadth against abundance.
#'
#' @details
#' **Allocation.** For pairwise isodars, the predicted abundances are those that
#' minimise the (optionally weighted) sum of squared deviations from all
#' pairwise isodar equations, \eqn{y - b x = a}, subject to summing to the
#' total, with no habitat below zero. This constrained least-squares problem
#' is solved exactly: if some habitats would be negative, the solution keeps
#' the best set of habitats at zero (an active-set method, which can bring
#' back a habitat dropped in an earlier step). For
#' an n-dimensional isodar (a line), the allocation is the single point at
#' which the line reaches the total; negative habitats are set to zero and
#' the others rescaled in proportion to reach the total (so a clipped
#' allocation is no longer on the line).
#'
#' **Density-independent niche breadth** (`abundances = NULL`) is evaluated at
#' the smallest total (searched over integers) at which every used habitat
#' holds at least `min_occupancy` individuals. If no total up to `max_search`
#' achieves this, a single row with `NA` is returned.
#'
#' **Absent habitats** (see [fit_isodars()]) are fixed at 0 and count as 0 in
#' the inequality, so they lower niche breadth; with a single used habitat,
#' niche breadth is 0.
#'
#' **Density-dependent niche breadth** is the value at a high total abundance,
#' supplied through `abundances`.
#'
#' @param isodars Output of [fit_isodars()].
#' @param abundances Numeric vector of total abundances, or `NULL` for the
#'   density-independent value.
#' @param weights Pairwise only: `NULL` (equal weights), `"1/var"` (inverse
#'   variance of each isodar's vertical residuals, \eqn{1/(sd^2(1+b^2))}), or a
#'   numeric vector with one weight per isodar.
#' @param scaled If `TRUE`, 1 - Gini is rescaled to 0-1, so that a species
#'   found in a single habitat scores 0 and an even spread scores 1.
#' @param min_occupancy Individuals required in every habitat for the
#'   density-independent value.
#' @param max_search Largest total searched for the density-independent value.
#' @param dim Optional; inferred from `isodars`. Kept for compatibility.
#' @param plot If `TRUE` and several abundances are given, plots the curve.
#' @param method Inequality measure; only `"gini"` is currently implemented.
#'
#' @return A tibble with `total_abundance`, `niche_breadth`, one column per
#'   habitat (predicted abundance), and `note`.
#' @seealso [fit_isodars()], [isoderr()], [pairwise_niche()]
#' @examples
#' set.seed(1)
#' sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
#' pw <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
#' iai(pw)                          # density-independent
#' iai(pw, abundances = c(50, 500)) # density-dependent
#' ln <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"), dim = "ndim")
#' iai(ln, abundances = 500)
#' @export
iai <- function(isodars, abundances = NULL, weights = NULL, scaled = TRUE,
                min_occupancy = 1, max_search = 10000, dim = NULL,
                plot = FALSE, method = "gini") {
  method <- match.arg(method, "gini")
  scaled <- isTRUE(as.logical(scaled))
  if (!is.null(abundances) &&
      (!is.numeric(abundances) || !all(is.finite(abundances)) || any(abundances < 0))) {
    stop("'abundances' must be NULL or finite, non-negative numbers.", call. = FALSE)
  }
  kind <- if (inherits(isodars, "isodar_line")) "ndim" else
    if (inherits(isodars, "isodars_pairwise")) "pairwise" else NA
  if (is.na(kind)) {
    # plain data frames from older code
    kind <- if (all(c("habitat_x", "habitat_y") %in% names(isodars))) "pairwise" else
      stop("'isodars' must be the output of fit_isodars().", call. = FALSE)
  }
  if (!is.null(dim) && !identical(dim, kind)) {
    stop("'dim = \"", dim, "\"' does not match the fitted isodars (", kind, ").",
         call. = FALSE)
  }
  alloc <- if (kind == "pairwise") pairwise_allocator(isodars, weights) else
    line_allocator(isodars)
  habitats <- attr(alloc, "habitats")
  H <- length(habitats)
  absent <- attr(isodars, "absent")
  used <- !habitats %in% absent
  absent_note <- if (any(!used)) {
    paste0("absent habitat(s) fixed at 0: ", paste(habitats[!used], collapse = ", "))
  } else ""

  if (is.null(abundances)) {
    first_ok <- NA_real_
    for (Total in seq_len(max_search)) {
      x <- alloc(Total)
      if (all(is.finite(x)) && all(x[used] >= min_occupancy - 1e-9)) {
        first_ok <- Total
        break
      }
    }
    if (is.na(first_ok)) {
      note <- paste0("no total up to ", max_search, " occupies all ",
                     if (any(!used)) "used " else "", "habitats")
      out <- tibble::as_tibble(c(list(total_abundance = NA_real_, niche_breadth = NA_real_),
                                 stats::setNames(as.list(rep(NA_real_, H)), habitats),
                                 list(note = note)))
      return(out)
    }
    abundances <- first_ok
  }
  rows <- lapply(as.numeric(abundances), function(Total) {
    x <- alloc(Total)
    nb <- 1 - gini(x)
    if (scaled) nb <- min(max((nb - 1 / H) * H / (H - 1), 0), 1) # rounding
    note <- if (any(x[used] == 0)) "habitat(s) at zero" else ""
    note <- paste(c(note, absent_note)[c(note, absent_note) != ""], collapse = "; ")
    c(list(total_abundance = Total, niche_breadth = nb),
      stats::setNames(as.list(x), habitats),
      list(note = note))
  })
  out <- tibble::as_tibble(do.call(rbind.data.frame, rows))
  if (plot && nrow(out) > 1) {
    graphics::plot(out$total_abundance, out$niche_breadth, type = "l",
                   xlab = "Total abundance", ylab = "Niche breadth (IAI)")
  }
  out
}

gini <- function(x) {
  x <- sort(as.numeric(x))
  n <- length(x)
  if (sum(x) <= 0) return(NA_real_)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
}

# Allocation along an n-dimensional isodar line
line_allocator <- function(fit) {
  mu <- fit$mean; v <- fit$direction; habs <- fit$habitat
  f <- function(Total) {
    x <- mu + (Total - sum(mu)) / sum(v) * v
    if (any(x < 0)) {
      x <- pmax(x, 0)
      if (sum(x) > 0) x <- x * Total / sum(x)
    }
    stats::setNames(x, habs)
  }
  attr(f, "habitats") <- habs
  f
}

# Pairwise isodar equations y - b x = a as a weighted linear system A x = b
pairwise_system <- function(isodars, weights = NULL) {
  habs <- attr(isodars, "habitats")
  if (is.null(habs)) habs <- unique(c(isodars$habitat_x, isodars$habitat_y))
  H <- length(habs); K <- nrow(isodars)
  A <- matrix(0, K, H, dimnames = list(NULL, habs))
  for (k in seq_len(K)) {
    A[k, isodars$habitat_y[k]] <- 1
    A[k, isodars$habitat_x[k]] <- -isodars$slope[k]
  }
  w <- if (is.null(weights)) rep(1, K) else if (identical(weights, "1/var")) {
    1 / pmax(isodars$sd^2 * (1 + isodars$slope^2), 1e-12)
  } else if (is.numeric(weights) && length(weights) == K) weights else
    stop("'weights' must be NULL, \"1/var\" or one number per isodar.", call. = FALSE)
  # absent habitats are fixed at 0
  list(A = A, b = isodars$intercept, w = w, habs = habs,
       start = which(!habs %in% attr(isodars, "absent")))
}

# Weighted least-squares solution over the free habitats, summing to Total
solve_free <- function(sys, Total, free) {
  Af <- sys$A[, free, drop = FALSE]
  m <- length(free)
  Q <- crossprod(Af, Af * sys$w)
  r <- crossprod(Af, sys$b * sys$w)
  M <- rbind(cbind(2 * Q, 1), c(rep(1, m), 0))
  s <- tryCatch(solve(M, c(2 * r, Total)), error = function(e) {
    tryCatch(qr.solve(M, c(2 * r, Total)), error = function(e2) {
      stop("The pairwise isodars do not determine the allocation among ",
           "habitats (singular system). This usually means a habitat with ",
           "(almost) no records was kept: declare it absent with 'absent', ",
           "or use dim = \"ndim\".", call. = FALSE)
    })
  })
  out <- s[seq_len(m)]
  # Lagrange multiplier of the sum constraint: gradient = -lambda on free habitats
  attr(out, "lambda") <- s[m + 1]
  out
}

# Exact constrained least-squares allocation (x >= 0, sum = Total) by a
# primal active-set method, started from a feasible allocation x
refine_allocation <- function(sys, Total, x) {
  H <- length(x)
  Q <- crossprod(sys$A, sys$A * sys$w)
  r <- crossprod(sys$A, sys$b * sys$w)
  candidates <- sys$start
  free <- candidates[x[candidates] > 0]
  tol <- 1e-10 * max(1, Total)
  for (iter in seq_len(10 * H + 10)) {
    if (length(free) == 0) break
    xs <- numeric(H)
    xs[free] <- solve_free(sys, Total, free)
    if (all(xs[free] >= -tol)) {
      x <- pmax(xs, 0)
      lambda <- attr(solve_free(sys, Total, free), "lambda")
      g <- as.vector(2 * (Q %*% x - r))
      # a habitat held at 0 should be freed if that lowers the objective
      out <- setdiff(candidates, free)
      viol <- g[out] + lambda
      if (length(out) == 0 || all(viol >= -1e-9 * max(1, abs(lambda)))) break
      free <- sort(c(free, out[which.min(viol)]))
    } else {
      # step towards the free solution until a habitat reaches 0
      neg <- free[xs[free] < 0]
      alpha <- min(x[neg] / (x[neg] - xs[neg]))
      x <- x + alpha * (xs - x)
      x[abs(x) < tol] <- 0
      free <- free[x[free] > 0]
    }
  }
  pmax(x, 0)
}

# Constrained least-squares allocation from pairwise isodars
pairwise_allocator <- function(isodars, weights = NULL) {
  sys <- pairwise_system(isodars, weights)
  habs <- sys$habs
  H <- length(habs)
  f <- function(Total) {
    free <- sys$start
    x <- numeric(H)
    clipped <- FALSE
    for (iter in seq_len(H)) {
      xf <- solve_free(sys, Total, free)
      if (all(xf >= 0) || length(free) == 1) { x[free] <- pmax(xf, 0); break }
      clipped <- TRUE
      free <- free[xf >= 0]
      if (length(free) == 0) break
    }
    if (sum(x) > 0) x <- x * Total / sum(x)
    # dropping habitats can overshoot: move to the exact optimum
    if (clipped && Total > 0 && sum(x) > 0) x <- refine_allocation(sys, Total, x)
    stats::setNames(x, habs)
  }
  attr(f, "habitats") <- habs
  f
}
