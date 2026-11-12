#' Fit isodars
#'
#' Fits isodars to abundances (or densities) recorded in two or more habitats
#' at a set of sites (or other replicates, such as years). Data are either in
#' wide format (one species; one row per site, one column per habitat) or in
#' long format (one or many species; one row per species x site x habitat).
#'
#' @details
#' **Input formats.** In wide format, give `hab_cols` (or `n_habitats`); the
#' values are used as they are, and the result is a single fit. In long
#' format, give the column names `site`, `habitat` and `abundance` (and
#' optionally `species` and `effort`); the result is an `isodar_fits` table
#' with one row per species, which [niche_breadth()] and [iai_curve()] take.
#' Replicates may be sites or time steps (for example, years at one site); in
#' long format they go in the `site` column. A site x habitat combination
#' that appears for any species is taken as sampled, and a species without a
#' row there has 0 individuals there; a habitat missing at a site simply has
#' no rows. Sites where a species has no individuals in any habitat carry no
#' information on its habitat selection and are dropped for that species.
#'
#' **Effort and densities.** Isodars relate densities. In long format,
#' density = abundance / effort x `reference_effort`. Effort may be an area, a
#' transect length (constant width) or a unitless relative effort, but must
#' be the same kind of measure in every row; it must be the same for every
#' row of a site x habitat combination. `reference_effort` (default: the
#' median effort of sampled site x habitat combinations) is the sampling unit
#' in which abundances are expressed: "one individual" in the
#' density-independent niche breadth, the spillover abundance and the x-axis
#' of [iai_curve()] all mean one individual per `reference_effort`.
#' Density-independent niche breadth and spillover abundance therefore depend
#' on this unit; slopes and density-dependent niche breadth do not. The
#' default depends on how samples happen to be pooled, so set
#' `reference_effort` to a natural unit (for example, the area of one standard
#' transect) and use the same value for every species you compare.
#'
#' In wide format, rows are used as given: unlike long format, a site with 0
#' in every habitat is kept (it pulls the isodars toward the origin), so
#' remove such rows if they mean the species was absent from the site.
#'
#' **Pairwise isodars** (`dim = "pairwise"`) are fitted for every pair of
#' habitats by total least squares (TLS), which is identical to major-axis
#' regression in two dimensions. Each pair uses every site where both habitats
#' were recorded (non-missing), so sites lacking some habitats still contribute.
#' When many sites are incomplete, both pairwise isodars and the latent line
#' (below) use all of them; fitting both is a useful check.
#'
#' **The n-dimensional isodar** (`dim = "ndim"`) is a line. Under the ideal
#' free distribution, equal fitness across \eqn{n} habitats imposes \eqn{n-1}
#' equalities, so abundances across habitats change along a single line as
#' total abundance changes, whatever the number of habitats. The line is
#' described by a point (`mean`) and a unit direction vector (`direction`),
#' with one component per habitat. Each component is the rate at which that
#' habitat gains individuals as total abundance rises (the analogue of the
#' isodar slope); where the line enters each habitat reflects its
#' density-independent preference (the analogue of the intercept). Every
#' pairwise isodar can be derived from it, and with two habitats it is the
#' classic isodar.
#'
#' With `ndim_method = "pca"` the line is the first principal axis of the
#' covariance matrix of complete sites (the TLS line). With
#' `ndim_method = "latent"` it is estimated from all sites with at least two
#' recorded habitats, as a one-factor model,
#' \eqn{N_{hs} = \mu_h + v_h t_s + \varepsilon_{hs}}, where \eqn{t_s} is the
#' position of site \eqn{s} along the isodar (requires \pkg{glmmTMB}).
#' `ndim_method = "auto"` uses `"pca"` when every site is complete and
#' `"latent"` otherwise.
#'
#' **Absent habitats.** A habitat that the species never uses carries no
#' information about habitat selection, and an isodar cannot be fitted to it.
#' Such a habitat is excluded from the fit and fixed at 0 abundance; [iai()]
#' still counts it, as 0, in the inequality, so it lowers niche breadth. With
#' `absent = "auto"` a habitat is declared absent if the species was recorded
#' in it no more than `absent_max_count` times in total (in the units of
#' wide `data`; as counts, before conversion to densities, in long format),
#' and the habitat was sampled at no fewer than `absent_min_sites`
#' sites where the species was present in another habitat. A habitat with
#' (almost) no records but fewer such sites is *undersampled*: it is kept in
#' the fit, with a warning, because its zeros may reflect too little sampling
#' rather than avoidance. A habitat used at low density and abandoned at high
#' density is not absent; it is fitted, and [iai()] clips it to 0 at totals
#' where the isodars predict negative abundance. With the default
#' `absent_max_count = 0`, a single individual makes a habitat used, which
#' can change the spillover abundance a lot; raise `absent_max_count` if a
#' few stray individuals should not count.
#'
#' Absent habitats lower niche breadth, so species should be compared over the
#' same set of habitats, namely the habitats available to them.
#'
#' @param data A data frame. Wide format: one row per site and one numeric
#'   column per habitat; `NA` marks a habitat that is absent from, or was not
#'   sampled at, a site. Long format: one row per species x site x habitat
#'   (see Details).
#' @param hab_cols Wide format: character vector of habitat column names.
#' @param n_habitats Alternatively, the number of habitat columns, taken from
#'   the start of `data`. One of `hab_cols` or `n_habitats` is required.
#' @param dim `"pairwise"` or `"ndim"`.
#' @param flip_intercept Pairwise only. If `TRUE`, an isodar with a negative
#'   intercept and positive slope is expressed in the reverse direction
#'   (habitats swapped), so that intercepts are non-negative.
#' @param ndim_method `"auto"`, `"pca"` or `"latent"`; see Details.
#' @param absent Habitats the species does not use: `"auto"` (detected from
#'   the data, see Details), `"none"`, or a character vector of habitat names
#'   (`character(0)` is the same as `"none"`).
#' @param absent_max_count,absent_min_sites Thresholds for `absent = "auto"`;
#'   see Details.
#' @param species,site,habitat,abundance,effort Long format: names of the
#'   columns holding species (optional; without it the data are one species),
#'   site (or other replicate), habitat, abundance (counts) and sampling
#'   effort (optional; without it every row has effort 1).
#' @param reference_effort Long format: the sampling unit in which densities
#'   are expressed (default: median effort); see Details.
#'
#' @return For `dim = "pairwise"`, a tibble of class `isodars_pairwise` with
#'   one row per pair of used habitats: `habitat_x`, `habitat_y`, `intercept`,
#'   `slope`, `sd` (SD of perpendicular residuals), `n_sites`, and `flipped`
#'   (whether the habitats were swapped). For `dim = "ndim"`, a tibble of class
#'   `isodar_line` with one row per habitat: `habitat`, `mean`, `direction`
#'   (both 0 for absent habitats); attributes `var_explained` (share of
#'   variance along the line, `"pca"` only), `eigenvalues`, `n_sites` and
#'   `ndim_method`. Both have attributes `habitats` (all habitats), `absent`
#'   and `undersampled`.
#'
#'   In long format, a tibble of class `isodar_fits` with one row per species:
#'   `species` (if given), `fit` (the fit as above, or `NULL` if it failed),
#'   `counts` and `effort` (site x habitat matrices), `n_sites` and `note`
#'   (warnings and errors from fitting). One summary warning is given if any
#'   species could not be fitted.
#'
#' @references Morris, D.W. (1987). Spatial scale and the cost of
#'   density-dependent habitat selection. *Evolutionary Ecology*, 1, 379–388.
#' @seealso [niche_breadth()], [iai()], [isoderr()]
#' @examples
#' set.seed(1)
#' sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
#' fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
#' fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"), dim = "ndim")
#'
#' # long format, several species
#' fish <- moorea_fish[moorea_fish$species %in% c("Acanthurus nigrofuscus",
#'   "Scarus schlegeli", "Zebrasoma scopas"), ]
#' fits <- fit_isodars(fish, species = "species", site = "site",
#'                     habitat = "habitat", abundance = "abundance")
#' fits
#' fits$fit[[which(fits$species == "Scarus schlegeli")]]
#' @export
fit_isodars <- function(data, hab_cols = NULL, n_habitats = NULL,
                        dim = c("pairwise", "ndim"), flip_intercept = TRUE,
                        ndim_method = c("auto", "pca", "latent"),
                        absent = "auto", absent_max_count = 0,
                        absent_min_sites = 5, species = NULL, site = NULL,
                        habitat = NULL, abundance = NULL, effort = NULL,
                        reference_effort = NULL) {
  dim <- match.arg(dim)
  ndim_method <- match.arg(ndim_method)
  if (!is.null(habitat)) {
    if (!is.null(hab_cols) || !is.null(n_habitats)) {
      stop("Give either 'hab_cols' (wide format) or 'habitat' (long format), ",
           "not both.", call. = FALSE)
    }
    args <- list(dim = dim, flip_intercept = flip_intercept,
                 ndim_method = ndim_method, absent = absent,
                 absent_max_count = absent_max_count,
                 absent_min_sites = absent_min_sites)
    return(fit_isodars_long(data, species, site, habitat, abundance, effort,
                            reference_effort, args))
  }
  if (!is.null(species) || !is.null(site) || !is.null(abundance)) {
    stop("Long format needs the 'habitat' column name.", call. = FALSE)
  }
  X <- habitat_matrix(data, hab_cols, n_habitats)
  habs <- colnames(X)
  ab <- find_absent(X, absent, absent_max_count, absent_min_sites)
  if (length(ab$undersampled) > 0) {
    warning("Habitat(s) with no records but too few informative sites to be ",
            "declared absent (kept in the fit): ",
            paste(ab$undersampled, collapse = ", "), ".", call. = FALSE)
  }
  used <- setdiff(habs, ab$absent)
  if (length(used) == 0) stop("Every habitat is absent.", call. = FALSE)
  Xu <- X[, used, drop = FALSE]
  res <- if (length(used) == 1) {
    single_habitat_fit(Xu, dim)
  } else if (dim == "pairwise") {
    fit_pairwise(Xu, flip_intercept)
  } else {
    fit_line(Xu, ndim_method)
  }
  if (dim == "ndim" && length(ab$absent) > 0) {
    # absent habitats sit at 0 and do not move along the line
    keep <- attributes(res)[c("var_explained", "eigenvalues", "n_sites",
                              "ndim_method")]
    full <- tibble::tibble(habitat = habs, mean = 0, direction = 0)
    full$mean[match(used, habs)] <- res$mean
    full$direction[match(used, habs)] <- res$direction
    res <- full
    attributes(res)[names(keep)] <- keep
    class(res) <- c("isodar_line", class(full))
  }
  attr(res, "habitats") <- habs
  attr(res, "absent") <- ab$absent
  attr(res, "undersampled") <- ab$undersampled
  res
}

# --- helpers ----------------------------------------------------------------

# Habitats never used by the species (fixed at 0), and habitats with no
# records but too few informative sites to tell
find_absent <- function(X, absent = "auto", max_count = 0, min_sites = 5) {
  habs <- colnames(X)
  if (is.character(absent) && length(absent) == 1 && absent %in% c("auto", "none")) {
    if (absent == "none") return(list(absent = character(0), undersampled = character(0)))
  } else if (is.character(absent)) {
    bad <- setdiff(absent, habs)
    if (length(bad) > 0) {
      stop("'absent' habitat(s) not among the habitats: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
    return(list(absent = habs[habs %in% absent], undersampled = character(0)))
  } else {
    stop("'absent' must be \"auto\", \"none\" or a character vector of habitats.",
         call. = FALSE)
  }
  if (!is.numeric(max_count) || length(max_count) != 1 || max_count < 0 ||
      !is.numeric(min_sites) || length(min_sites) != 1 || min_sites < 1) {
    stop("'absent_max_count' must be >= 0 and 'absent_min_sites' >= 1.",
         call. = FALSE)
  }
  absent <- undersampled <- character(0)
  for (h in seq_along(habs)) {
    total <- sum(X[, h], na.rm = TRUE)
    if (total > max_count) next
    others <- X[, -h, drop = FALSE]
    present_elsewhere <- rowSums(others > 0, na.rm = TRUE) > 0
    informative <- sum(!is.na(X[, h]) & present_elsewhere)
    if (informative >= min_sites) {
      absent <- c(absent, habs[h])
    } else if (informative > 0) {
      undersampled <- c(undersampled, habs[h])
    }
  }
  list(absent = absent, undersampled = undersampled)
}

# Only one used habitat: no isodar, every individual is in that habitat
single_habitat_fit <- function(X, dim) {
  h <- colnames(X)
  if (dim == "pairwise") {
    res <- tibble::tibble(habitat_x = character(0), habitat_y = character(0),
                          intercept = numeric(0), slope = numeric(0),
                          sd = numeric(0), n_sites = integer(0),
                          flipped = logical(0))
    class(res) <- c("isodars_pairwise", class(res))
  } else {
    res <- tibble::tibble(habitat = h, mean = mean(X[, 1], na.rm = TRUE),
                          direction = 1)
    attr(res, "var_explained") <- NA_real_
    attr(res, "eigenvalues") <- NA_real_
    attr(res, "n_sites") <- sum(!is.na(X[, 1]))
    attr(res, "ndim_method") <- "single habitat"
    class(res) <- c("isodar_line", class(res))
  }
  res
}

habitat_matrix <- function(data, hab_cols, n_habitats) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  if (is.null(hab_cols)) {
    if (is.null(n_habitats)) {
      stop("Specify the habitat columns with 'hab_cols' (or 'n_habitats').",
           call. = FALSE)
    }
    if (n_habitats > ncol(data)) stop("'n_habitats' exceeds the number of columns.", call. = FALSE)
    hab_cols <- names(data)[seq_len(n_habitats)]
  }
  missing_cols <- setdiff(hab_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Habitat column(s) not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (length(hab_cols) < 2) stop("At least two habitats are required.", call. = FALSE)
  X <- data[, hab_cols, drop = FALSE]
  if (!all(vapply(X, is.numeric, logical(1)))) {
    stop("Habitat columns must be numeric: ", paste(hab_cols, collapse = ", "), call. = FALSE)
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (any(X < 0, na.rm = TRUE)) stop("Abundances must be non-negative.", call. = FALSE)
  X
}

# TLS (major-axis) line for one pair
tls_pair <- function(x, y) {
  mx <- mean(x); my <- mean(y)
  sxx <- stats::var(x); syy <- stats::var(y); sxy <- stats::cov(x, y)
  theta <- 0.5 * atan2(2 * sxy, sxx - syy)
  slope <- tan(theta)
  intercept <- my - slope * mx
  resid <- ((y - my) * cos(theta) - (x - mx) * sin(theta))
  list(intercept = intercept, slope = slope,
       sd = sqrt(sum(resid^2) / max(length(x) - 2, 1)))
}

fit_pairwise <- function(X, flip_intercept = TRUE) {
  habs <- colnames(X)
  out <- list()
  for (i in seq_len(ncol(X) - 1)) {
    for (j in (i + 1):ncol(X)) {
      ok <- stats::complete.cases(X[, c(i, j)])
      x <- X[ok, i]; y <- X[ok, j]
      if (sum(ok) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
        warning("Isodar ", habs[i], "-", habs[j],
                " not fitted (fewer than 3 shared sites or no variation).", call. = FALSE)
        next
      }
      f <- tls_pair(x, y)
      hx <- habs[i]; hy <- habs[j]; flipped <- FALSE
      if (flip_intercept && is.finite(f$intercept) && f$intercept < 0 &&
          is.finite(f$slope) && f$slope > 0) {
        # TLS is symmetric in x and y: no refit needed
        f$intercept <- -f$intercept / f$slope
        f$slope <- 1 / f$slope
        hx <- habs[j]; hy <- habs[i]; flipped <- TRUE
      }
      out[[length(out) + 1]] <- data.frame(
        habitat_x = hx, habitat_y = hy, intercept = f$intercept, slope = f$slope,
        sd = f$sd, n_sites = sum(ok), flipped = flipped, stringsAsFactors = FALSE)
    }
  }
  if (length(out) == 0) stop("No pairwise isodar could be fitted.", call. = FALSE)
  res <- tibble::as_tibble(do.call(rbind, out))
  attr(res, "habitats") <- habs
  class(res) <- c("isodars_pairwise", class(res))
  res
}

fit_line <- function(X, ndim_method = "auto") {
  habs <- colnames(X)
  complete <- stats::complete.cases(X)
  if (ndim_method == "auto") ndim_method <- if (all(complete)) "pca" else "latent"
  if (ndim_method == "pca") {
    Xc <- X[complete, , drop = FALSE]
    if (nrow(Xc) < ncol(X) + 1) {
      stop("Not enough complete sites for ndim_method = 'pca'; try 'latent'.", call. = FALSE)
    }
    mu <- colMeans(Xc)
    e <- eigen(stats::cov(Xc), symmetric = TRUE)
    v <- e$vectors[, 1]
    ev <- e$values
    var_exp <- ev[1] / sum(ev)
    n_used <- nrow(Xc)
  } else {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("ndim_method = 'latent' requires the 'glmmTMB' package.", call. = FALSE)
    }
    long <- data.frame(site = factor(rep(seq_len(nrow(X)), ncol(X))),
                       habitat = factor(rep(habs, each = nrow(X)), levels = habs),
                       y = as.vector(X))
    long <- long[!is.na(long$y), ]
    keep <- names(which(table(long$site) >= 2))
    long <- droplevels(long[long$site %in% keep, ])
    check_connected(X[as.integer(keep), , drop = FALSE])
    m <- glmmTMB::glmmTMB(y ~ 0 + habitat + rr(0 + habitat | site, d = 1), data = long)
    L <- m$obj$env$report(m$fit$parfull)$fact_load[[1]][, 1]
    v <- L / sqrt(sum(L^2))
    mu <- as.numeric(glmmTMB::fixef(m)$cond)
    ev <- NA_real_; var_exp <- NA_real_
    n_used <- nlevels(long$site)
  }
  if (sum(v) < 0) v <- -v
  res <- tibble::tibble(habitat = habs, mean = as.numeric(mu), direction = as.numeric(v))
  attr(res, "var_explained") <- var_exp
  attr(res, "eigenvalues") <- ev
  attr(res, "n_sites") <- n_used
  attr(res, "ndim_method") <- ndim_method
  class(res) <- c("isodar_line", class(res))
  res
}

# habitats must be linked through sites where they co-occur
check_connected <- function(X) {
  H <- ncol(X)
  co <- crossprod(!is.na(X)) > 0
  reached <- 1
  repeat {
    new <- unique(c(reached, which(colSums(co[reached, , drop = FALSE]) > 0)))
    if (length(new) == length(reached)) break
    reached <- new
  }
  if (length(reached) < H) {
    stop("Habitats are not connected through shared sites: ",
         paste(colnames(X)[-reached], collapse = ", "),
         " never co-occur with the others.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @export
print.isodar_line <- function(x, ...) {
  cat("n-dimensional isodar (line), method:", attr(x, "ndim_method"),
      "| sites:", attr(x, "n_sites"), "\n")
  ve <- attr(x, "var_explained")
  if (!is.null(ve) && is.finite(ve)) cat("Variance along the line:", round(ve, 3), "\n")
  ab <- attr(x, "absent")
  if (length(ab) > 0) cat("Absent habitats (fixed at 0):", paste(ab, collapse = ", "), "\n")
  NextMethod()
}
