#' Density-independent and density-dependent niche breadth
#'
#' Fits isodars (unless given fitted ones) and returns, for each species, the
#' density-independent (DI) and density-dependent (DD) niche breadth with
#' bootstrap confidence intervals, the spillover abundance, preference
#' reversals and diagnostic flags. With a `species` column, everything runs
#' per species and the result has one row per species; otherwise the data are
#' one species.
#'
#' @details
#' **Niche breadth** is the isodar-adjusted inequality ([iai()]): total
#' abundance is allocated among habitats as the isodars predict and niche
#' breadth is 1 - Gini of that allocation (scaled to 0-1 by default). Absent
#' habitats count as 0 (see [fit_isodars()]).
#'
#' **Spillover abundance** \eqn{T_0} is the smallest total abundance (an
#' integer, in individuals per reference effort; see [fit_isodars()]) at which
#' every used habitat holds at least `min_occupancy` individuals: the point at
#' which the species starts to spill over into its least-favoured used
#' habitat. DI niche breadth is the value at \eqn{T_0}. It depends on the
#' sampling unit (`reference_effort`); DD niche breadth does not.
#'
#' **DD niche breadth.** With `dependent_at = "asymptote"` (default) it is the
#' limit as total abundance \eqn{T \to \infty}, where the allocation is
#' proportional to the direction of the isodar line (or follows the pairwise
#' slopes) and intercepts no longer matter; it is computed directly. A number
#' gives the value at that total, and `"max_observed"` the value at the
#' largest site total in the data.
#'
#' **The curve between them.** If no two habitats change rank and no habitat
#' is abandoned anywhere above \eqn{T_0} (no `reversal` or `abandoned` flag),
#' niche breadth follows
#' \deqn{NB(T) = NB_{DD} + (NB_{DI} - NB_{DD}) T_0 / T,}{NB(T) = NB_DD + (NB_DI - NB_DD) T0 / T,}
#' and is halfway between the two values at \eqn{T = 2 T_0}. Otherwise it
#' has this form only piecewise, with different constants between
#' successive reversals. Where the predicted abundances of two habitats
#' cross, the ranking of habitats changes and the curve can be
#' non-monotonic; these preference-reversal abundances are returned in
#' `reversals` (totals from \eqn{T_0} upward) and flagged. A habitat used at low density whose predicted abundance falls to
#' 0 at high totals is listed in `abandoned` (it still counts for DI niche
#' breadth). Use [iai_curve()] to compute and plot whole curves.
#'
#' **Bootstrap.** Whole sites are resampled with replacement and the isodars
#' refitted, keeping the absent habitats (and the line method) of the full
#' data. Intervals are percentile intervals.
#'
#' **Diagnostics.** Each `flag_*` column is `TRUE` when:
#' * `few_sites`: fewer than `min_sites` sites;
#' * `few_pair_sites`: some pair of used habitats shares fewer than
#'   `min_pair_sites` sites;
#' * `low_abundance`: the median site total is below `min_median_count`
#'   individuals (counts, not densities), or more than `max_zero_share` of the
#'   sampled cells of used habitats are zeros;
#' * `di_extrapolation`: \eqn{T_0} lies outside the observed range of site
#'   totals by more than a factor `extrapolation_ratio`;
#' * `dd_extrapolation`: for the asymptote, the largest site total is below
#'   `extrapolation_ratio` x \eqn{T_0} (with the default 2, the data do not
#'   reach the halfway point of the curve); for a numeric `dependent_at`, it
#'   exceeds the largest site total by more than that factor;
#' * `abandoned`: some used habitat is abandoned at high density
#'   (informational);
#' * `reversal`: habitats change rank along the curve (informational);
#' * `poor_fit`: the line explains less than `min_var_explained` of the
#'   variance (`ndim_method = "pca"`), or the pairwise isodars are
#'   inconsistent: for some three habitats \eqn{i, j, k} the slopes
#'   \eqn{b_{ij} b_{jk}} and \eqn{b_{ik}} differ in sign or by more than a
#'   factor `max_slope_ratio`;
#' * `disconnected`: habitats are not connected through shared sites;
#' * `intercept_noise`: no isodar intercept differs from those expected from
#'   sampling alone. Following the paper (section 1.2.3), `n_null` data sets
#'   are simulated without density-dependent habitat selection: each site's
#'   individuals are redistributed among its sampled habitats in proportion to
#'   the species' overall density in each habitat (multinomial). Two-sided
#'   p-values of the pairwise TLS intercepts are Holm-corrected within
#'   species; `intercept_p` is the smallest, and the flag is set when it
#'   exceeds `alpha`. The test needs integer counts (`NA` otherwise). Set
#'   `n_null = 0` to skip it;
#' * `unstable_bootstrap`: more than `max_failed_share` of bootstrap
#'   replicates failed, or a confidence interval is wider than `max_ci_width`;
#' * `undersampled`: a habitat had no records at too few sites to be declared
#'   absent (see [fit_isodars()]);
#' * `no_di`: no total up to `max_search` occupies every used habitat;
#' * `fit_failed`: the isodars could not be fitted (see `note`).
#'
#' Defaults: `min_sites = 10`, `min_pair_sites = 5`, `min_median_count = 5`,
#' `max_zero_share = 0.5`, `extrapolation_ratio = 2`,
#' `min_var_explained = 0.5`, `max_slope_ratio = 2`, `n_null = 199`,
#' `alpha = 0.05`, `max_failed_share = 0.1`, `max_ci_width = 0.5`. Change any
#' of them through `diagnostics`, e.g. `diagnostics = list(min_sites = 8)`.
#'
#' **Comparing species.** Absent habitats lower niche breadth, so species
#' compared with each other should be assessed over the same set of habitats:
#' the habitats available to them.
#'
#' @param fit_or_data Data in long or wide format (see [fit_isodars()]), or the
#'   long-format output of [fit_isodars()] (an `isodar_fits` table).
#' @param ... Passed to [fit_isodars()] when `fit_or_data` is data: the column
#'   names (`species`, `site`, `habitat`, `abundance`, `effort`) or `hab_cols`,
#'   and `dim`, `ndim_method`, `flip_intercept`, `absent_max_count`,
#'   `absent_min_sites`.
#' @param dependent_at `"asymptote"`, `"max_observed"` or a total abundance;
#'   see Details.
#' @param min_occupancy Individuals required in every used habitat at the
#'   spillover abundance.
#' @param reference_effort,absent Passed to [fit_isodars()] (data only).
#' @param n_boot Number of bootstrap replicates (0 for none).
#' @param conf_level Confidence level of the bootstrap intervals.
#' @param scaled If `TRUE`, niche breadth is scaled to 0-1 (see [iai()]).
#' @param weights Pairwise isodars only: `"1/var"` (default, as in the paper)
#'   or `NULL`; see [iai()].
#' @param diagnostics Named list of thresholds for the flags; see Details.
#' @param seed Optional random-number seed (bootstrap and null model). Each
#'   species then gets its own random streams, so its results do not depend on
#'   `n_boot` or on the other species; your random-number state is restored
#'   afterwards.
#' @param max_search Largest total searched for the spillover abundance.
#'
#' @return A tibble of class `niche_breadth`, one row per species, with
#'   `species` (if given), `method`, `n_sites`, `nb_di`, `nb_di_lower`,
#'   `nb_di_upper`, `spillover_abundance`, `spillover_lower`,
#'   `spillover_upper`, `nb_dd`, `nb_dd_lower`, `nb_dd_upper`, `dependent_at`
#'   (`Inf` for the asymptote), `var_explained` (line, `"pca"` only),
#'   `intercept_p`, `absent` and `abandoned` (habitats, comma-separated),
#'   `reversals` (list of data frames: `total_abundance`, `habitat_1` (more
#'   abundant before), `habitat_2`, `niche_breadth`), `n_boot_valid`, `flags`
#'   (names of the flags raised, comma-separated), one logical `flag_*` column
#'   per diagnostic, and `note`. Printing it adds one line giving the number of
#'   flagged species.
#' @seealso [iai_curve()], [fit_isodars()], [check_coverage()]
#' @examples
#' three <- moorea_fish[moorea_fish$species %in% c("Acanthurus nigrofuscus",
#'   "Scarus schlegeli", "Siganus argenteus"), ]
#' nb <- niche_breadth(three, species = "species", site = "site",
#'                     habitat = "habitat", abundance = "abundance",
#'                     n_boot = 99, seed = 1)
#' nb
#' nb[c("species", "nb_di", "spillover_abundance", "nb_dd", "flags")]
#' @export
niche_breadth <- function(fit_or_data, ..., dependent_at = "asymptote",
                          min_occupancy = 1, reference_effort = NULL,
                          absent = "auto", n_boot = 999, conf_level = 0.95,
                          scaled = TRUE, weights = "1/var", diagnostics = list(),
                          seed = NULL, max_search = 10000) {
  fits <- as_isodar_fits(fit_or_data, ..., reference_effort = reference_effort,
                         absent = absent, absent_given = !missing(absent))
  if (!(identical(dependent_at, "asymptote") || identical(dependent_at, "max_observed") ||
        (is.numeric(dependent_at) && length(dependent_at) == 1 &&
         is.finite(dependent_at) && dependent_at > 0))) {
    stop("'dependent_at' must be \"asymptote\", \"max_observed\" or a positive number.",
         call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 0) {
    stop("'n_boot' must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be between 0 and 1.", call. = FALSE)
  }
  if (!is.null(weights) && !identical(weights, "1/var")) {
    stop("'weights' must be \"1/var\" or NULL.", call. = FALSE)
  }
  settings <- list(dependent_at = dependent_at, min_occupancy = min_occupancy,
                   n_boot = as.integer(n_boot), conf_level = conf_level,
                   scaled = isTRUE(as.logical(scaled)), weights = weights,
                   diagnostics = nb_diagnostics(diagnostics),
                   max_search = max_search)
  ref <- attr(fits, "reference_effort")
  args <- attr(fits, "args")
  # with a seed, each species gets its own streams for the bootstrap and the
  # null model, so neither depends on n_boot or on the other species
  seeds <- NULL
  if (!is.null(seed)) {
    old_rng <- save_rng()
    on.exit(restore_rng(old_rng), add = TRUE)
    set.seed(seed)
    seeds <- sample.int(.Machine$integer.max, 2 * nrow(fits))
  }
  rows <- lapply(seq_len(nrow(fits)), function(i) {
    nb_species(fits$fit[[i]], fits$counts[[i]], fits$effort[[i]], ref, args,
               settings, fits$note[[i]], seeds = seeds[2 * i - c(1, 0)])
  })
  col <- function(name, type) vapply(rows, `[[`, type, name)
  out <- tibble::tibble(
    species = fits$species,
    method = col("method", character(1)),
    n_sites = fits$n_sites,
    nb_di = col("nb_di", numeric(1)),
    nb_di_lower = col("nb_di_lower", numeric(1)),
    nb_di_upper = col("nb_di_upper", numeric(1)),
    spillover_abundance = col("spillover", numeric(1)),
    spillover_lower = col("spillover_lower", numeric(1)),
    spillover_upper = col("spillover_upper", numeric(1)),
    nb_dd = col("nb_dd", numeric(1)),
    nb_dd_lower = col("nb_dd_lower", numeric(1)),
    nb_dd_upper = col("nb_dd_upper", numeric(1)),
    dependent_at = col("dependent_at", numeric(1)),
    var_explained = col("var_explained", numeric(1)),
    intercept_p = col("intercept_p", numeric(1)),
    absent = col("absent", character(1)),
    abandoned = col("abandoned", character(1)),
    reversals = lapply(rows, `[[`, "reversals"),
    n_boot_valid = col("n_boot_valid", integer(1)),
    flags = vapply(rows, function(r) paste(flag_names[r$flags], collapse = ", "),
                   character(1))
  )
  flag_mat <- do.call(rbind, lapply(rows, `[[`, "flags"))
  for (f in flag_names) out[[paste0("flag_", f)]] <- unname(flag_mat[, f])
  err <- vapply(rows, `[[`, character(1), "error")
  out$note <- ifelse(nzchar(err) & nzchar(fits$note), paste(fits$note, err, sep = "; "),
                     paste0(fits$note, err))
  if (!isTRUE(attr(fits, "species_given"))) out$species <- NULL
  attr(out, "fits") <- fits
  attr(out, "settings") <- settings
  class(out) <- c("niche_breadth", class(out))
  out
}

flag_names <- c("few_sites", "few_pair_sites", "low_abundance",
                "di_extrapolation", "dd_extrapolation", "abandoned", "reversal",
                "poor_fit", "disconnected", "intercept_noise",
                "unstable_bootstrap", "undersampled", "no_di", "fit_failed")

nb_diagnostics <- function(diagnostics = list()) {
  defaults <- list(min_sites = 10, min_pair_sites = 5, min_median_count = 5,
                   max_zero_share = 0.5, extrapolation_ratio = 2,
                   min_var_explained = 0.5, max_slope_ratio = 2, n_null = 199,
                   alpha = 0.05, max_failed_share = 0.1, max_ci_width = 0.5)
  if (!is.list(diagnostics) || (length(diagnostics) > 0 && is.null(names(diagnostics)))) {
    stop("'diagnostics' must be a named list.", call. = FALSE)
  }
  bad <- setdiff(names(diagnostics), names(defaults))
  if (length(bad) > 0) {
    stop("Unknown diagnostics: ", paste(bad, collapse = ", "), ". Known: ",
         paste(names(defaults), collapse = ", "), ".", call. = FALSE)
  }
  ok <- vapply(diagnostics, function(x) is.numeric(x) && length(x) == 1 && x >= 0,
               logical(1))
  if (!all(ok)) stop("Diagnostic thresholds must be single non-negative numbers.",
                     call. = FALSE)
  defaults[names(diagnostics)] <- diagnostics
  defaults
}

fit_args <- function(dim = "pairwise", flip_intercept = TRUE, ndim_method = "auto",
                     absent = "auto", absent_max_count = 0, absent_min_sites = 5) {
  list(dim = match.arg(dim, c("pairwise", "ndim")), flip_intercept = flip_intercept,
       ndim_method = match.arg(ndim_method, c("auto", "pca", "latent")),
       absent = absent, absent_max_count = absent_max_count,
       absent_min_sites = absent_min_sites)
}

# Data or fitted isodars -> isodar_fits
as_isodar_fits <- function(x, ..., reference_effort = NULL, absent = "auto",
                           absent_given = FALSE) {
  dots <- list(...)
  if (inherits(x, "isodar_fits")) {
    if (length(dots) > 0 || !is.null(reference_effort) || absent_given) {
      stop("Isodars are already fitted: set fitting options in fit_isodars().",
           call. = FALSE)
    }
    return(x)
  }
  if (inherits(x, c("isodars_pairwise", "isodar_line"))) {
    stop("Give the data, or the long-format output of fit_isodars(); a single ",
         "fit has no data for the bootstrap and diagnostics.", call. = FALSE)
  }
  if (!is.data.frame(x)) {
    stop("'fit_or_data' must be a data frame or the output of fit_isodars().",
         call. = FALSE)
  }
  if (!is.null(dots$habitat)) {
    warn_unused_effort(x, dots$effort)
    return(suppressWarnings(do.call(fit_isodars, c(list(x), dots,
      list(reference_effort = reference_effort, absent = absent)))))
  }
  if (!is.null(reference_effort)) {
    stop("'reference_effort' needs long-format data.", call. = FALSE)
  }
  wide <- c("hab_cols", "n_habitats")
  known <- c(wide, names(formals(fit_args)))
  bad <- setdiff(names(dots), known)
  if (length(bad) > 0) {
    stop("Unknown argument(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  args <- do.call(fit_args, c(dots[setdiff(names(dots), wide)], list(absent = absent)))
  wide_as_fits(x, dots$hab_cols, dots$n_habitats, args)
}

nb_species <- function(fit, C, E, ref, args, st, note, seeds = NULL) {
  dg <- st$diagnostics
  flags <- stats::setNames(rep(FALSE, length(flag_names)), flag_names)
  r <- list(method = NA_character_, nb_di = NA_real_, nb_di_lower = NA_real_,
            nb_di_upper = NA_real_, spillover = NA_real_,
            spillover_lower = NA_real_, spillover_upper = NA_real_,
            nb_dd = NA_real_, nb_dd_lower = NA_real_, nb_dd_upper = NA_real_,
            dependent_at = NA_real_, var_explained = NA_real_,
            intercept_p = NA_real_, absent = "", abandoned = "",
            reversals = find_reversals(NULL, NULL, NA, NA), n_boot_valid = 0L,
            error = "")
  habs <- colnames(C)
  flags["few_sites"] <- nrow(C) < dg$min_sites
  flags["undersampled"] <- grepl("undersampled habitat", note, fixed = TRUE)
  if (is.null(fit)) {
    flags[if (grepl("not connected", note)) "disconnected" else "fit_failed"] <- TRUE
    r$flags <- flags
    return(r)
  }
  D <- density_matrix(C, E, ref)
  absent <- attr(fit, "absent")
  used <- !habs %in% absent
  is_line <- inherits(fit, "isodar_line")
  nm <- attr(fit, "ndim_method")
  r$method <- if ((is_line && identical(nm, "single habitat")) ||
                  (!is_line && nrow(fit) == 0)) {
    "single habitat"
  } else if (is_line) paste0("line (", nm, ")") else "pairwise"
  totals <- site_totals(D)
  max_obs <- max(totals)
  min_obs <- min(totals[totals > 0])
  est <- tryCatch(nb_core(fit, st$weights, st$min_occupancy, st$max_search,
                          st$scaled, st$dependent_at, max_obs),
                  error = function(e) conditionMessage(e))
  if (is.character(est)) {
    flags["fit_failed"] <- TRUE
    r$error <- paste0("niche breadth failed: ", est)
    r$flags <- flags
    return(r)
  }
  r$nb_di <- est$nb_di
  r$spillover <- est$spillover
  r$nb_dd <- est$nb_dd
  r$dependent_at <- if (identical(st$dependent_at, "asymptote")) Inf else
    if (identical(st$dependent_at, "max_observed")) max_obs else st$dependent_at
  if (is_line) r$var_explained <- attr(fit, "var_explained")
  r$absent <- paste(absent, collapse = ", ")
  asym <- asymptote_allocation(fit, est$weights)
  abandoned <- habs[used & asym <= 1e-12]
  r$abandoned <- paste(abandoned, collapse = ", ")
  lo <- if (is.finite(est$spillover)) est$spillover else 1
  r$reversals <- find_reversals(fit, est$weights, lo, max(1e4 * lo, 100 * max_obs),
                                st$scaled)

  if (st$n_boot > 0) {
    if (!is.null(seeds)) set.seed(seeds[1])
    B <- matrix(NA_real_, st$n_boot, 3)
    for (b in seq_len(st$n_boot)) {
      fb <- boot_fit(D, fit, args)
      if (is.null(fb)) next
      rb <- tryCatch(nb_core(fb, st$weights, st$min_occupancy, st$max_search,
                             st$scaled, st$dependent_at, max_obs),
                     error = function(e) NULL)
      if (!is.null(rb)) B[b, ] <- c(rb$nb_di, rb$spillover, rb$nb_dd)
    }
    a <- (1 - st$conf_level) / 2
    q <- function(x) {
      if (!any(is.finite(x))) return(c(NA_real_, NA_real_))
      stats::quantile(x, c(a, 1 - a), na.rm = TRUE, names = FALSE)
    }
    r[c("nb_di_lower", "nb_di_upper")] <- as.list(q(B[, 1]))
    r[c("spillover_lower", "spillover_upper")] <- as.list(q(B[, 2]))
    r[c("nb_dd_lower", "nb_dd_upper")] <- as.list(q(B[, 3]))
    r$n_boot_valid <- sum(is.finite(B[, 3]))
    wide_ci <- c(r$nb_di_upper - r$nb_di_lower, r$nb_dd_upper - r$nb_dd_lower)
    flags["unstable_bootstrap"] <-
      1 - r$n_boot_valid / st$n_boot > dg$max_failed_share ||
      any(wide_ci > dg$max_ci_width, na.rm = TRUE)
  }

  # diagnostics
  uh <- which(used)
  if (length(uh) >= 2) {
    shared <- crossprod(!is.na(D[, uh, drop = FALSE]))
    flags["few_pair_sites"] <- min(shared[upper.tri(shared)]) < dg$min_pair_sites
  }
  cells <- C[, uh, drop = FALSE]
  flags["low_abundance"] <-
    stats::median(rowSums(C, na.rm = TRUE)) < dg$min_median_count ||
    mean(cells[!is.na(cells)] == 0) > dg$max_zero_share
  T0 <- est$spillover
  ratio <- dg$extrapolation_ratio
  if (is.finite(T0)) {
    flags["di_extrapolation"] <- T0 < min_obs / ratio || T0 > max_obs * ratio
  }
  flags["dd_extrapolation"] <- if (identical(st$dependent_at, "asymptote")) {
    is.finite(T0) && max_obs < ratio * T0
  } else {
    r$dependent_at > ratio * max_obs
  }
  flags["abandoned"] <- length(abandoned) > 0
  flags["reversal"] <- nrow(r$reversals) > 0
  pw <- if (!is_line) fit else tryCatch(suppressWarnings(fit_isodars(
    as.data.frame(D, optional = TRUE), hab_cols = habs, absent = absent)),
    error = function(e) NULL)
  flags["poor_fit"] <- if (is_line) {
    is.finite(r$var_explained) && r$var_explained < dg$min_var_explained
  } else {
    pairwise_inconsistent(fit, dg$max_slope_ratio)
  }
  if (!is.null(seeds)) set.seed(seeds[2])
  r$intercept_p <- intercept_null_p(C, E, ref, pw, dg$n_null)
  flags["intercept_noise"] <- is.finite(r$intercept_p) && r$intercept_p > dg$alpha
  flags["no_di"] <- is.na(T0)
  r$flags <- flags
  r
}

# Slopes of three habitats i, j, k should satisfy b_ij * b_jk = b_ik
pairwise_inconsistent <- function(pw, max_ratio = 2) {
  habs <- unique(c(pw$habitat_x, pw$habitat_y))
  H <- length(habs)
  if (H < 3) return(FALSE)
  S <- matrix(NA_real_, H, H, dimnames = list(habs, habs))
  for (k in seq_len(nrow(pw))) {
    S[pw$habitat_x[k], pw$habitat_y[k]] <- pw$slope[k]
    if (pw$slope[k] != 0) S[pw$habitat_y[k], pw$habitat_x[k]] <- 1 / pw$slope[k]
  }
  for (i in 1:(H - 2)) for (j in (i + 1):(H - 1)) for (k in (j + 1):H) {
    s1 <- S[i, j]; s2 <- S[j, k]; s3 <- S[i, k]
    if (anyNA(c(s1, s2, s3))) next
    p <- s1 * s2
    if (sign(p) != sign(s3)) return(TRUE)
    if (s3 != 0 && abs(log(p / s3)) > log(max_ratio)) return(TRUE)
  }
  FALSE
}

# Smallest Holm-corrected p-value of the pairwise intercepts against a
# multinomial null without density-dependent habitat selection
intercept_null_p <- function(C, E, ref, pw, n_null) {
  if (is.null(pw) || nrow(pw) == 0 || n_null < 1) return(NA_real_)
  if (any(abs(C - round(C)) > 1e-8, na.rm = TRUE)) return(NA_real_)
  S <- nrow(C)
  H <- ncol(C)
  q <- colSums(C, na.rm = TRUE) / colSums(E, na.rm = TRUE)
  n_s <- round(rowSums(C, na.rm = TRUE))
  sims <- array(NA_real_, c(S, H, n_null), dimnames = list(NULL, colnames(C), NULL))
  for (s in seq_len(S)) {
    ok <- !is.na(C[s, ])
    p <- q[ok] * E[s, ok]
    if (sum(p) > 0) sims[s, ok, ] <- stats::rmultinom(n_null, n_s[s], p)
  }
  p <- vapply(seq_len(nrow(pw)), function(k) {
    hx <- pw$habitat_x[k]
    hy <- pw$habitat_y[k]
    x <- matrix(sims[, hx, ], S, n_null)
    y <- matrix(sims[, hy, ], S, n_null)
    rows <- !is.na(x[, 1]) & !is.na(y[, 1])
    x <- x[rows, , drop = FALSE] * ref / E[rows, hx]
    y <- y[rows, , drop = FALSE] * ref / E[rows, hy]
    a <- tls_intercepts(x, y)
    a <- a[is.finite(a)]
    if (length(a) == 0) return(NA_real_)
    obs <- pw$intercept[k]
    min(1, 2 * (min(sum(a >= obs), sum(a <= obs)) + 1) / (length(a) + 1))
  }, numeric(1))
  if (all(is.na(p))) return(NA_real_)
  min(stats::p.adjust(p, "holm"), na.rm = TRUE)
}

# TLS intercepts for each column (replicate) of x and y
tls_intercepts <- function(x, y) {
  mx <- colMeans(x)
  my <- colMeans(y)
  xc <- sweep(x, 2, mx)
  yc <- sweep(y, 2, my)
  sxx <- colSums(xc^2)
  syy <- colSums(yc^2)
  sxy <- colSums(xc * yc)
  slope <- tan(0.5 * atan2(2 * sxy, sxx - syy))
  out <- my - slope * mx
  out[sxx == 0 | syy == 0] <- NA
  out
}

#' @export
print.niche_breadth <- function(x, ...) {
  NextMethod()
  # abandoned habitats and reversals are informational
  warn <- intersect(paste0("flag_", setdiff(flag_names, c("abandoned", "reversal"))),
                    names(x))
  if (length(warn) > 0) {
    fl <- as.matrix(as.data.frame(x)[warn])
    n <- sum(rowSums(fl, na.rm = TRUE) > 0)
    unit <- if ("species" %in% names(x)) " species" else ""
    cat(n, unit, " flagged (not counting abandoned habitats or reversals); ",
        "see summary()\n", sep = "")
  }
  fits <- attr(x, "fits")
  if (!is.null(fits)) cat(abundance_unit(fits), "\n", sep = "")
  invisible(x)
}

# The unit of abundance (spillover abundance, curve x-axis) in words
abundance_unit <- function(fits) {
  if (isTRUE(attr(fits, "effort_given"))) {
    paste0("Abundances are per reference effort of ",
           format(attr(fits, "reference_effort")), " in each habitat.")
  } else {
    "Abundances are per site (equal effort in every site and habitat)."
  }
}

# What each flag means and what to do about it (used by summary())
flag_info <- data.frame(
  flag = flag_names,
  type = c("data", "data", "data", "extrapolation", "extrapolation", "info",
           "info", "fit", "fit", "fit", "uncertainty", "data", "fit", "fit"),
  meaning = c(
    "fewer sites than min_sites",
    "a pair of used habitats shares fewer than min_pair_sites sites",
    "small counts per site, or mostly zeros",
    "spillover abundance outside the range of observed site totals",
    "the data do not reach the abundances that determine DD niche breadth",
    "a used habitat is abandoned at high density",
    "habitats change rank along the niche breadth curve",
    "a line does not describe the data well, or the pairwise isodars disagree",
    "habitats are not connected through shared sites",
    "isodar intercepts do not differ from sampling noise",
    "many bootstrap replicates failed, or a confidence interval is very wide",
    "zeros in a habitat sampled at too few sites to call it absent",
    "no total abundance occupies every used habitat",
    "isodars could not be fitted"
  ),
  advice = c(
    "group nearby sites (group_sites()) or treat as unreliable",
    "use dim = \"ndim\" with ndim_method = \"latent\", or group sites",
    "group sites, or report as unreliable; estimates rest on few individuals",
    "DI niche breadth is extrapolated; interpret with care",
    "DD niche breadth is extrapolated; interpret with care",
    "informational: the habitat counts for DI but not for DD niche breadth",
    "informational: compare species at a stated abundance (see iai_curve())",
    "check the isodars (fit_isodars()); try the other method",
    "add sites that sample habitats together, or drop a habitat",
    "no evidence that niche breadth changes with abundance; report DD niche breadth",
    "increase n_boot; if still wide, the data cannot resolve niche breadth",
    "sample more sites, or name absent habitats with 'absent'",
    "check the isodars; the least-favoured habitat may never fill",
    "see the 'note' column"
  ),
  stringsAsFactors = FALSE
)

#' Summarise niche breadth results and their diagnostics
#'
#' Summarises a [niche_breadth()] result: how many species have estimates,
#' the range of niche breadth values, and for each diagnostic flag how many
#' species raised it, what it means and what to do about it.
#'
#' @param object The output of [niche_breadth()].
#' @param x The output of `summary()`.
#' @param ... Ignored.
#' @return An object of class `summary.niche_breadth`: a list with
#'   `n_species`, `n_estimated` (species with DI and DD niche breadth),
#'   `n_clean` (estimated species with no warning flags), `clean` (their
#'   names, or row numbers without a `species` column), `estimates` (quartiles
#'   of `nb_di`, `spillover_abundance` and `nb_dd`), `flags` (a data frame
#'   with `flag`, `type`, `n`, `meaning` and `advice`, for the flags raised)
#'   and `settings`. Flags of type `"info"` (abandoned habitats, preference
#'   reversals) describe the species and are not warnings.
#' @seealso [niche_breadth()], `vignette("diagnostics")`
#' @examples
#' three <- moorea_fish[moorea_fish$species %in% c("Acanthurus nigrofuscus",
#'   "Scarus schlegeli", "Siganus argenteus"), ]
#' nb <- niche_breadth(three, species = "species", site = "site",
#'                     habitat = "habitat", abundance = "abundance",
#'                     n_boot = 0, seed = 1)
#' summary(nb)
#' @export
summary.niche_breadth <- function(object, ...) {
  x <- as.data.frame(object)
  fl <- intersect(paste0("flag_", flag_names), names(x))
  counts <- vapply(fl, function(f) sum(x[[f]], na.rm = TRUE), numeric(1))
  info <- flag_info[match(sub("^flag_", "", fl), flag_info$flag), ]
  info$n <- as.integer(counts)
  info <- info[info$n > 0, c("flag", "type", "n", "meaning", "advice")]
  info <- info[order(info$type == "info", -info$n), ]
  rownames(info) <- NULL
  warn <- fl[!sub("^flag_", "", fl) %in% flag_info$flag[flag_info$type == "info"]]
  estimated <- is.finite(x$nb_di) & is.finite(x$nb_dd)
  clean <- estimated & rowSums(as.matrix(x[warn]), na.rm = TRUE) == 0
  ids <- if ("species" %in% names(x)) x$species else seq_len(nrow(x))
  q <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0) return(rep(NA_real_, 5))
    stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  }
  est <- rbind(nb_di = q(x$nb_di), spillover_abundance = q(x$spillover_abundance),
               nb_dd = q(x$nb_dd))
  colnames(est) <- c("min", "q25", "median", "q75", "max")
  out <- list(n_species = nrow(x), n_estimated = sum(estimated),
              n_clean = sum(clean), clean = ids[clean], estimates = est,
              flags = info, settings = attr(object, "settings"),
              species_given = "species" %in% names(x))
  class(out) <- "summary.niche_breadth"
  out
}

#' @rdname summary.niche_breadth
#' @export
print.summary.niche_breadth <- function(x, ...) {
  unit <- if (x$species_given) " species" else " data set(s)"
  st <- x$settings
  dd <- if (is.null(st)) "" else if (identical(st$dependent_at, "asymptote")) {
    "asymptote"
  } else format(st$dependent_at)
  cat("Niche breadth for ", x$n_species, unit, "\n", sep = "")
  if (!is.null(st)) {
    cat("  DD niche breadth at: ", dd, "; bootstrap replicates: ", st$n_boot,
        "\n", sep = "")
  }
  cat("  with DI and DD estimates: ", x$n_estimated, "\n", sep = "")
  cat("  without warning flags:    ", x$n_clean, "\n\n", sep = "")
  if (x$n_estimated > 0) {
    cat("Estimates (quartiles):\n")
    print(round(x$estimates, 3))
    cat("\n")
  }
  if (nrow(x$flags) == 0) {
    cat("No flags raised.\n")
  } else {
    cat("Flags raised (species):\n")
    for (i in seq_len(nrow(x$flags))) {
      f <- x$flags[i, ]
      cat(sprintf("  %-19s %4d  %s\n", f$flag, f$n, f$meaning))
      cat(sprintf("  %-19s       -> %s\n", "", f$advice))
    }
    cat("\nSee vignette(\"diagnostics\") for details. Thresholds are set with",
        "niche_breadth(diagnostics = ).\n")
  }
  invisible(x)
}
