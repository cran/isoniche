# Internal helpers shared by niche_breadth() and iai_curve().
#
# Before any habitat is clipped to zero, the allocation of a total T among
# habitats is affine in T: x(T) = c + T u, with sum(c) = 0 and sum(u) = 1.
# Absent habitats have c = u = 0.

allocator <- function(fit, weights = NULL) {
  if (inherits(fit, "isodar_line")) line_allocator(fit) else
    pairwise_allocator(fit, weights)
}

affine_parts <- function(fit, weights = NULL) {
  if (inherits(fit, "isodar_line")) {
    mu <- fit$mean
    v <- fit$direction
    u <- v / sum(v)
    return(list(c = mu - sum(mu) * u, u = u, habs = fit$habitat))
  }
  sys <- pairwise_system(fit, weights)
  x0 <- x1 <- numeric(length(sys$habs))
  x0[sys$start] <- solve_free(sys, 0, sys$start)
  x1[sys$start] <- solve_free(sys, 1, sys$start)
  list(c = x0, u = x1 - x0, habs = sys$habs)
}

# Smallest integer total (>= 1) at which every used habitat holds at least
# m individuals; the same rule as iai(abundances = NULL), solved directly.
di_total <- function(af, used, m = 1, max_search = 10000) {
  ok <- function(Total) all(af$c[used] + Total * af$u[used] >= m - 1e-9)
  lo <- 1
  hi <- max_search
  for (h in which(used)) {
    if (af$u[h] > 1e-12) {
      lo <- max(lo, (m - 1e-9 - af$c[h]) / af$u[h])
    } else if (af$u[h] < -1e-12) {
      hi <- min(hi, (m - 1e-9 - af$c[h]) / af$u[h])
    } else if (af$c[h] < m - 1e-9) {
      return(NA_real_)
    }
  }
  if (lo > hi + 1e-6 || lo > max_search) return(NA_real_)
  Total <- max(1, ceiling(lo - 1e-6))
  while (Total > 1 && ok(Total - 1)) Total <- Total - 1
  while (Total <= max_search && !ok(Total)) Total <- Total + 1
  if (Total > max_search) NA_real_ else Total
}

# Allocation share as total abundance tends to infinity
asymptote_allocation <- function(fit, weights = NULL) {
  if (inherits(fit, "isodar_line")) {
    u <- pmax(fit$direction, 0)
    return(stats::setNames(u / sum(u), fit$habitat))
  }
  fit$intercept <- 0
  pairwise_allocator(fit, weights)(1)
}

scaled_nb <- function(x, scaled = TRUE) {
  H <- length(x)
  nb <- 1 - gini(x)
  if (scaled) nb <- min(max((nb - 1 / H) * H / (H - 1), 0), 1)
  nb
}

# Totals (>= lo) at which two used habitats swap rank in predicted abundance
find_reversals <- function(fit, weights, lo, hi, scaled = TRUE) {
  none <- data.frame(total_abundance = numeric(0), habitat_1 = character(0),
                     habitat_2 = character(0), niche_breadth = numeric(0),
                     stringsAsFactors = FALSE)
  habs <- attr(fit, "habitats")
  used <- which(!habs %in% attr(fit, "absent"))
  if (length(used) < 2 || !is.finite(lo) || hi <= lo) return(none)
  f <- allocator(fit, weights)
  grid <- 10^seq(log10(lo), log10(hi), length.out = 400)
  X <- vapply(grid, f, numeric(length(habs)))
  out <- list()
  for (a in seq_along(used)[-length(used)]) {
    for (b in (a + 1):length(used)) {
      i <- used[a]; j <- used[b]
      d <- X[i, ] - X[j, ]
      # habitats tied along the curve (e.g. symmetric isodars) do not swap
      # rank: ignore sign changes that are only rounding error
      nz <- which(abs(d) > 1e-9 * pmax(abs(X[i, ]), abs(X[j, ]), 1))
      dn <- d[nz]
      for (m in which(dn[-length(dn)] * dn[-1] < 0)) {
        k <- nz[m]
        r <- stats::uniroot(function(t) {
          x <- f(t)
          x[i] - x[j]
        }, grid[c(k, nz[m + 1])], tol = 1e-8)$root
        first <- if (d[k] > 0) c(i, j) else c(j, i)
        out[[length(out) + 1]] <- data.frame(
          total_abundance = r, habitat_1 = habs[first[1]],
          habitat_2 = habs[first[2]], niche_breadth = scaled_nb(f(r), scaled),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(out) == 0) return(none)
  out <- do.call(rbind, out)
  out[order(out$total_abundance), , drop = FALSE]
}

# Refit to a bootstrap resample of sites, keeping the absent habitats and the
# line method of the original fit; NULL if the fit fails
boot_fit <- function(D, fit, args) {
  nm <- attr(fit, "ndim_method")
  nm <- if (!is.null(nm) && nm %in% c("pca", "latent")) nm else "auto"
  Db <- D[sample.int(nrow(D), replace = TRUE), , drop = FALSE]
  tryCatch(suppressWarnings(fit_isodars(
    as.data.frame(Db, optional = TRUE), hab_cols = colnames(D), dim = args$dim,
    flip_intercept = args$flip_intercept, ndim_method = nm,
    absent = attr(fit, "absent"))), error = function(e) NULL)
}

# Save and restore the user's random-number state around a seeded call
save_rng <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
}

restore_rng <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
}

# Total of each site (complete sites when there are any)
site_totals <- function(D) {
  complete <- stats::complete.cases(D)
  rowSums(D[if (any(complete)) complete else TRUE, , drop = FALSE], na.rm = TRUE)
}

# Niche breadth at the density-independent point and the density-dependent
# value for one fit
nb_core <- function(fit, weights, min_occupancy, max_search, scaled,
                    dependent_at, max_observed) {
  if (inherits(fit, "isodar_line")) weights <- NULL
  habs <- attr(fit, "habitats")
  used <- !habs %in% attr(fit, "absent")
  af <- affine_parts(fit, weights)
  T0 <- di_total(af, used, min_occupancy, max_search)
  nb_di <- if (is.na(T0)) NA_real_ else scaled_nb(af$c + T0 * af$u, scaled)
  x_dd <- if (identical(dependent_at, "asymptote")) {
    asymptote_allocation(fit, weights)
  } else {
    at <- if (identical(dependent_at, "max_observed")) max_observed else dependent_at
    allocator(fit, weights)(at)
  }
  list(nb_di = nb_di, spillover = T0, nb_dd = scaled_nb(x_dd, scaled),
       asymptote = x_dd, used = used, weights = weights)
}
