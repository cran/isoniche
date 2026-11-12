#' Niche breadth as a function of total abundance
#'
#' Computes the isodar-adjusted inequality along a range of total abundances
#' for each species, with optional bootstrap bands, and the points needed to
#' plot it: the density-independent value at each species' spillover
#' abundance, the asymptote, preference reversals and the totals at which the
#' curves of two species cross.
#'
#' @details
#' By default each species' curve runs from its spillover abundance
#' \eqn{T_0} to \eqn{1000 T_0} (see [niche_breadth()]), usually far beyond
#' the abundances observed; the plot draws the part beyond the largest
#' observed site total as a thinner line. `totals` sets a common range
#' instead; totals below a species' \eqn{T_0} are allowed (some used habitat
#' then holds less than one individual), so curves, reversals and crossings
#' can lie below \eqn{T_0}. Totals are in individuals per reference effort
#' (see [fit_isodars()]).
#'
#' The bands are pointwise percentile intervals at each fixed total. They
#' differ from the intervals of `nb_di` in [niche_breadth()], where the
#' spillover abundance is re-estimated in every bootstrap replicate; report
#' those for density-independent niche breadth.
#'
#' Where the curves of two species cross, their ranking by niche breadth
#' depends on the abundance at which they are compared. These rank reversals
#' are returned in attribute `"crossings"`, using the totals where both curves
#' are defined.
#'
#' @param fit Output of [niche_breadth()] (whose settings are reused), the
#'   long-format output of [fit_isodars()], or a single fit from
#'   [fit_isodars()] (no bootstrap).
#' @param totals Optional numeric vector of total abundances.
#' @param n_boot Number of bootstrap replicates for the bands (0 for none).
#' @param conf_level Confidence level of the bands.
#' @param weights,scaled,min_occupancy,max_search As in [niche_breadth()];
#'   ignored when `fit` is a `niche_breadth` result.
#' @param seed Optional random-number seed.
#'
#' @return A tibble of class `iai_curve` with `species` (if given),
#'   `total_abundance`, `relative_abundance` (total / \eqn{T_0}),
#'   `niche_breadth`, and `lower` and `upper` (bootstrap band; `NA` without
#'   bootstrap). Attributes: `"points"` (per species: `spillover_abundance`,
#'   `nb_di`, `nb_dd` at the asymptote, and `max_observed`, the largest
#'   site total in the data), `"reversals"` (preference reversals
#'   within the range) and `"crossings"` (`species_1`, `species_2`,
#'   `total_abundance`, `niche_breadth`).
#' @seealso [plot.iai_curve()], [niche_breadth()]
#' @examples
#' three <- moorea_fish[moorea_fish$species %in% c("Acanthurus nigrofuscus",
#'   "Scarus schlegeli", "Siganus argenteus"), ]
#' fits <- fit_isodars(three, species = "species", site = "site",
#'                     habitat = "habitat", abundance = "abundance", dim = "ndim")
#' curve <- iai_curve(fits, n_boot = 50, seed = 1)
#' plot(curve)
#' plot(curve, x_axis = "relative")
#' attr(curve, "crossings")
#' @export
iai_curve <- function(fit, totals = NULL, n_boot = 0, conf_level = 0.95,
                      weights = "1/var", scaled = TRUE, min_occupancy = 1,
                      seed = NULL, max_search = 10000) {
  if (inherits(fit, "niche_breadth")) {
    st <- attr(fit, "settings")
    fits <- attr(fit, "fits")
    weights <- st$weights
    scaled <- st$scaled
    min_occupancy <- st$min_occupancy
    max_search <- st$max_search
  } else if (inherits(fit, "isodar_fits")) {
    fits <- fit
  } else if (inherits(fit, c("isodars_pairwise", "isodar_line"))) {
    if (n_boot > 0) {
      stop("Bootstrap bands need the data: use the long-format output of ",
           "fit_isodars().", call. = FALSE)
    }
    fits <- tibble::tibble(species = "", fit = list(fit), counts = list(NULL),
                           effort = list(NULL))
    attr(fits, "args") <- list(dim = if (inherits(fit, "isodar_line")) "ndim" else
      "pairwise", flip_intercept = TRUE)
  } else {
    stop("'fit' must be the output of niche_breadth() or fit_isodars().",
         call. = FALSE)
  }
  if (!is.null(totals) && (!is.numeric(totals) || length(totals) == 0 ||
                           any(!is.finite(totals)) || any(totals <= 0))) {
    stop("'totals' must be NULL or positive numbers.", call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 0) {
    stop("'n_boot' must be a non-negative number.", call. = FALSE)
  }
  scaled <- isTRUE(as.logical(scaled))
  if (!is.null(seed)) {
    old_rng <- save_rng()
    on.exit(restore_rng(old_rng), add = TRUE)
    set.seed(seed)
  }
  args <- attr(fits, "args")
  a <- (1 - conf_level) / 2
  curves <- points <- revs <- list()
  for (i in seq_len(nrow(fits))) {
    f <- fits$fit[[i]]
    if (is.null(f)) next
    sp <- fits$species[i]
    w <- if (inherits(f, "isodar_line")) NULL else weights
    used <- !attr(f, "habitats") %in% attr(f, "absent")
    af <- affine_parts(f, w)
    T0 <- di_total(af, used, min_occupancy, max_search)
    grid <- if (!is.null(totals)) sort(unique(totals)) else
      (if (is.finite(T0)) T0 else 1) * 10^seq(0, 3, length.out = 101)
    nb_on <- function(fit_k) {
      alloc <- allocator(fit_k, if (inherits(fit_k, "isodar_line")) NULL else weights)
      vapply(grid, function(t) scaled_nb(alloc(t), scaled), numeric(1))
    }
    nb <- nb_on(f)
    lower <- upper <- rep(NA_real_, length(grid))
    D <- if (!is.null(fits$counts[[i]])) {
      density_matrix(fits$counts[[i]], fits$effort[[i]],
                     attr(fits, "reference_effort"))
    }
    max_obs <- if (is.null(D)) NA_real_ else max(site_totals(D))
    if (n_boot > 0 && !is.null(D)) {
      B <- matrix(NA_real_, n_boot, length(grid))
      for (b in seq_len(n_boot)) {
        fb <- boot_fit(D, f, args)
        if (!is.null(fb)) B[b, ] <- tryCatch(nb_on(fb), error = function(e) NA_real_)
      }
      if (any(is.finite(B))) {
        lower <- apply(B, 2, stats::quantile, probs = a, na.rm = TRUE, names = FALSE)
        upper <- apply(B, 2, stats::quantile, probs = 1 - a, na.rm = TRUE, names = FALSE)
      }
    }
    curves[[length(curves) + 1]] <- data.frame(
      species = sp, total_abundance = grid, relative_abundance = grid / T0,
      niche_breadth = nb, lower = lower, upper = upper, stringsAsFactors = FALSE)
    points[[length(points) + 1]] <- data.frame(
      species = sp, spillover_abundance = T0,
      nb_di = if (is.finite(T0)) scaled_nb(af$c + T0 * af$u, scaled) else NA_real_,
      nb_dd = scaled_nb(asymptote_allocation(f, w), scaled),
      max_observed = max_obs, stringsAsFactors = FALSE)
    rv <- find_reversals(f, w, min(grid), max(grid), scaled)
    if (nrow(rv) > 0) {
      revs[[length(revs) + 1]] <- data.frame(species = sp, rv,
                                             relative_abundance = rv$total_abundance / T0,
                                             stringsAsFactors = FALSE)
    }
  }
  if (length(curves) == 0) stop("No species could be fitted.", call. = FALSE)
  out <- tibble::as_tibble(do.call(rbind, curves))
  pts <- tibble::as_tibble(do.call(rbind, points))
  rv <- if (length(revs) > 0) tibble::as_tibble(do.call(rbind, revs)) else
    tibble::tibble(species = character(0), total_abundance = numeric(0),
                   habitat_1 = character(0), habitat_2 = character(0),
                   niche_breadth = numeric(0), relative_abundance = numeric(0))
  cr <- curve_crossings(out)
  if (!isTRUE(attr(fits, "species_given"))) {
    out$species <- NULL
    pts$species <- NULL
    rv$species <- NULL
  }
  attr(out, "points") <- pts
  attr(out, "reversals") <- rv
  attr(out, "crossings") <- cr
  attr(out, "scaled") <- scaled
  attr(out, "unit") <- if (!is.null(attr(fits, "reference_effort"))) {
    list(effort_given = isTRUE(attr(fits, "effort_given")),
         reference_effort = attr(fits, "reference_effort"))
  }
  class(out) <- c("iai_curve", class(out))
  out
}

# Totals at which the niche breadth curves of two species cross
curve_crossings <- function(curves) {
  none <- tibble::tibble(species_1 = character(0), species_2 = character(0),
                         total_abundance = numeric(0), niche_breadth = numeric(0))
  sp <- unique(curves$species)
  if (length(sp) < 2) return(none)
  lg <- log(curves$total_abundance)
  G <- seq(min(lg), max(lg), length.out = 400)
  Y <- vapply(sp, function(s) {
    d <- curves[curves$species == s, ]
    t <- log(d$total_abundance)
    y <- if (length(t) > 1) stats::approx(t, d$niche_breadth, G)$y else rep(NA_real_, length(G))
    y
  }, numeric(length(G)))
  out <- list()
  for (i in seq_along(sp)[-length(sp)]) {
    for (j in (i + 1):length(sp)) {
      d <- Y[, i] - Y[, j]
      k <- which(d[-length(d)] * d[-1] < 0)
      if (length(k) == 0) next
      t <- G[k] - d[k] * (G[k + 1] - G[k]) / (d[k + 1] - d[k])
      out[[length(out) + 1]] <- data.frame(
        species_1 = sp[i], species_2 = sp[j], total_abundance = exp(t),
        niche_breadth = stats::approx(G, Y[, i], t)$y, stringsAsFactors = FALSE)
    }
  }
  if (length(out) == 0) return(none)
  tibble::as_tibble(do.call(rbind, out))
}

#' Plot niche breadth curves
#'
#' Plots the curves from [iai_curve()]: one line per species with bootstrap
#' bands (dotted), a point at each species' density-independent value (at its
#' own spillover abundance), the asymptote (density-dependent value) as a
#' short horizontal dash at the right, and preference reversals (crosses).
#' Beyond the largest total observed at any site of a species, the curve is
#' drawn as a thinner line: there it extrapolates the fitted isodars.
#'
#' @param x Output of [iai_curve()].
#' @param species Optional species to plot (default: all).
#' @param x_axis `"log_abundance"` (total abundance, log scale) or `"relative"`
#'   (total abundance / spillover abundance, log scale; every species' DI
#'   value is then at 1).
#' @param abundances Optional total abundances to mark with vertical lines
#'   (`"log_abundance"` only).
#' @param crossings If `TRUE`, marks the totals where two species' curves
#'   cross (open circles; `"log_abundance"` only).
#' @param legend `TRUE` (default) adds a legend (up to 12 species) in a
#'   corner free of curves and bands, or above the plot if there is none; a
#'   position such as `"topleft"` or `"bottomright"` (see
#'   [graphics::legend()]), or `"outside"` (above the plot), places it there;
#'   `FALSE` omits it.
#' @param col Colours, one per species (default: the current palette).
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @seealso [iai_curve()]
#' @export
plot.iai_curve <- function(x, species = NULL, x_axis = c("log_abundance", "relative"),
                           abundances = NULL, crossings = FALSE, legend = TRUE,
                           col = NULL, ...) {
  x_axis <- match.arg(x_axis)
  sp_of <- function(d) if ("species" %in% names(d)) d[["species"]] else rep("", nrow(d))
  pts <- attr(x, "points")
  rv <- attr(x, "reversals")
  cr <- attr(x, "crossings")
  all_sp <- unique(sp_of(x))
  if (!is.null(species)) {
    bad <- setdiff(species, all_sp)
    if (length(bad) > 0) stop("Species not in 'x': ", paste(bad, collapse = ", "),
                              call. = FALSE)
    all_sp <- species
  }
  if (is.null(col)) col <- seq_along(all_sp)
  col <- rep_len(col, length(all_sp))
  xvar <- if (x_axis == "log_abundance") "total_abundance" else "relative_abundance"
  keep <- sp_of(x) %in% all_sp & is.finite(x[[xvar]])
  if (!any(keep)) stop("Nothing to plot (no spillover abundance for the relative axis?).",
                       call. = FALSE)
  xr <- range(x[[xvar]][keep])
  dash <- xr[2] * c(1.15, 1.6)
  ylim <- if (isTRUE(attr(x, "scaled"))) c(0, 1) else
    range(c(x$niche_breadth[keep], x$lower[keep], x$upper[keep]), na.rm = TRUE)
  graphics::plot(NA, xlim = c(xr[1], dash[2]), ylim = ylim, log = "x",
                 xaxt = "n",
                 xlab = if (x_axis == "log_abundance") curve_xlab(attr(x, "unit")) else
                   "Total abundance / spillover abundance",
                 ylab = "Niche breadth (IAI)", ...)
  decades <- 10^(floor(log10(xr[1])):ceiling(log10(dash[2])))
  at <- sort(c(decades, 2 * decades, 5 * decades))
  at <- at[at >= xr[1] / 1.01 & at <= xr[2] * 1.01]
  # over a wide range, label only powers of 10 so that no label is dropped
  lab <- if (log10(xr[2] / xr[1]) > 2.5) at %in% decades else rep(TRUE, length(at))
  graphics::axis(1, at = at, labels = FALSE)
  graphics::axis(1, at = at[lab], labels = format(at[lab], big.mark = ",",
                 scientific = FALSE, trim = TRUE, drop0trailing = TRUE), tick = FALSE)
  if (x_axis == "log_abundance" && !is.null(abundances)) {
    graphics::abline(v = abundances, lty = 2, col = "grey50")
  }
  for (k in seq_along(all_sp)) {
    s <- all_sp[k]
    d <- x[sp_of(x) == s, ]
    p <- pts[sp_of(pts) == s, ]
    # thinner beyond the largest observed site total (extrapolation)
    obs <- if ("max_observed" %in% names(p) && is.finite(p$max_observed)) {
      d$total_abundance <= p$max_observed
    } else rep(TRUE, nrow(d))
    graphics::lines(d[[xvar]], d$niche_breadth, col = col[k], lwd = 1)
    if (any(obs)) {
      last <- max(which(obs))
      seg <- seq_len(min(last + 1, nrow(d)))
      graphics::lines(d[[xvar]][seg], d$niche_breadth[seg], col = col[k], lwd = 2.5)
    }
    if (any(is.finite(d$lower))) {
      graphics::lines(d[[xvar]], d$lower, col = col[k], lty = 3)
      graphics::lines(d[[xvar]], d$upper, col = col[k], lty = 3)
    }
    px <- if (x_axis == "log_abundance") p$spillover_abundance else 1
    if (is.finite(px) && is.finite(p$nb_di)) {
      graphics::points(px, p$nb_di, pch = 19, col = col[k])
    }
    graphics::segments(dash[1], p$nb_dd, dash[2], p$nb_dd, col = col[k], lwd = 3)
    r <- rv[sp_of(rv) == s, ]
    if (nrow(r) > 0) graphics::points(r[[xvar]], r$niche_breadth, pch = 4, col = col[k])
  }
  if (crossings && x_axis == "log_abundance" && nrow(cr) > 0) {
    cr <- cr[cr$species_1 %in% all_sp & cr$species_2 %in% all_sp, ]
    graphics::points(cr$total_abundance, cr$niche_breadth, pch = 1, cex = 1.6)
  }
  pos <- if (is.character(legend)) legend else if (isTRUE(legend)) {
    lx <- log10(rep(x[[xvar]][keep], 3))
    ly <- c(x$niche_breadth[keep], x$lower[keep], x$upper[keep])
    emptiest_corner(lx, ly, log10(c(xr[1], dash[2])), ylim)
  } else NULL
  if (!is.null(pos) && length(all_sp) <= 12 && any(nzchar(all_sp))) {
    if (identical(pos, "outside")) {
      graphics::legend("bottom", legend = all_sp, col = col, lwd = 2, bty = "n",
                       cex = 0.8, ncol = min(3, length(all_sp)), inset = c(0, 1),
                       xpd = TRUE)
    } else {
      graphics::legend(pos, legend = all_sp, col = col, lwd = 2, bty = "n",
                       cex = 0.8)
    }
  }
  invisible(x)
}

#' Plot niche breadth against abundance
#'
#' Computes the niche breadth curves of a [niche_breadth()] result with
#' [iai_curve()] and plots them with [plot.iai_curve()].
#'
#' @param x Output of [niche_breadth()].
#' @param n_boot,seed Passed to [iai_curve()] (bootstrap bands; default none).
#' @param ... Passed to [plot.iai_curve()], e.g. `species` or `x_axis`.
#' @return The [iai_curve()] result, invisibly.
#' @seealso [iai_curve()], [plot.iai_curve()]
#' @examples
#' three <- moorea_fish[moorea_fish$species %in% c("Acanthurus nigrofuscus",
#'   "Scarus schlegeli", "Siganus argenteus"), ]
#' nb <- niche_breadth(three, species = "species", site = "site",
#'                     habitat = "habitat", abundance = "abundance",
#'                     n_boot = 0, diagnostics = list(n_null = 0))
#' plot(nb)
#' @export
plot.niche_breadth <- function(x, n_boot = 0, seed = NULL, ...) {
  curve <- iai_curve(x, n_boot = n_boot, seed = seed)
  plot(curve, ...)
  invisible(curve)
}

# x-axis label from attribute "unit" of an iai_curve (NULL: unknown unit)
curve_xlab <- function(unit) {
  if (is.null(unit)) return("Total abundance")
  if (isTRUE(unit$effort_given)) {
    paste0("Total abundance (individuals per reference effort of ",
           format(unit$reference_effort), ")")
  } else {
    "Total abundance (individuals per site)"
  }
}

# Corner of the plot without curve or band points (for the legend);
# "outside" (above the plot) when every corner has some
emptiest_corner <- function(x, y, xlim, ylim) {
  ok <- is.finite(x) & is.finite(y)
  u <- (x[ok] - xlim[1]) / diff(xlim)
  v <- (y[ok] - ylim[1]) / diff(ylim)
  corners <- c("topleft", "topright", "bottomleft", "bottomright")
  n <- c(sum(u < 0.45 & v > 0.6), sum(u > 0.55 & v > 0.6),
         sum(u < 0.45 & v < 0.4), sum(u > 0.55 & v < 0.4))
  if (min(n) > 0) "outside" else corners[which.min(n)]
}
