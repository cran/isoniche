#' Sampling coverage of each species and habitat
#'
#' Summarises, for each species and habitat, how many of the species' sites
#' sampled the habitat, how often the species was found there, and whether
#' [fit_isodars()] would treat the habitat as used, absent or undersampled.
#' Run it before fitting, to decide which species and habitats can be
#' analysed and whether sites should be grouped ([group_sites()]).
#'
#' @inheritParams fit_isodars
#' @param data Long-format data: one row per species x site x habitat.
#' @return A tibble with one row per species x habitat: `species` (if given),
#'   `habitat`, `n_sites` (sites where the species was recorded in any
#'   habitat), `n_complete` (of those, sites where every habitat was sampled),
#'   `n_sampled` (sites where this habitat was sampled), `n_occupied` (sites
#'   where the species was recorded in this habitat), `total_count` and
#'   `status` (`"used"`, `"absent"` or `"undersampled"`).
#' @seealso [fit_isodars()], [group_sites()]
#' @examples
#' cov <- check_coverage(moorea_fish, species = "species", site = "site",
#'                       habitat = "habitat", abundance = "abundance")
#' cov[cov$status != "used", ]
#' @export
check_coverage <- function(data, species = NULL, site, habitat, abundance,
                           effort = NULL, absent_max_count = 0,
                           absent_min_sites = 5) {
  prep <- prepare_long(data, species, site, habitat, abundance, effort)
  habs <- prep$habitats
  rows <- lapply(seq_along(prep$species), function(i) {
    s <- names(prep$species)[i]
    C <- prep$species[[i]]$counts
    ab <- find_absent(C, "auto", absent_max_count, absent_min_sites)
    status <- ifelse(habs %in% ab$absent, "absent",
                     ifelse(habs %in% ab$undersampled, "undersampled", "used"))
    data.frame(species = s, habitat = habs, n_sites = nrow(C),
               n_complete = sum(stats::complete.cases(C)),
               n_sampled = colSums(!is.na(C)),
               n_occupied = colSums(C > 0, na.rm = TRUE),
               total_count = colSums(C, na.rm = TRUE),
               status = status, stringsAsFactors = FALSE, row.names = NULL)
  })
  out <- tibble::as_tibble(do.call(rbind, rows))
  if (!prep$species_given) out$species <- NULL
  out
}

#' Group nearby sites into sampling units
#'
#' Merges sites within a given distance of each other into one sampling unit,
#' summing counts and efforts. Grouping helps when single sites lack some
#' habitats or hold too few individuals, at the cost of fewer units.
#'
#' @details
#' Sites are clustered by complete linkage on great-circle distance, so every
#' pair of sites in a group lies within `radius_km` of each other (single
#' linkage, `linkage = "single"`, can chain distant sites together). With
#' `radius_km = 0` only sites with identical coordinates are merged.
#'
#' Counts are summed over the sites of a group, and so are efforts: the
#' effort of a group in a habitat is the summed effort of the sites that
#' sampled that habitat. Without an `effort` column, every site x habitat has
#' effort 1, so the group's effort is the number of sites that sampled the
#' habitat. Groups are named after their sites (joined with `"+"`).
#'
#' The radius is not chosen automatically: choosing it to maximise
#' completeness or niche breadth would be circular. Use
#' [radius_sensitivity()] to see how the number of units, their completeness
#' and the niche breadth estimates change with the radius.
#'
#' @inheritParams fit_isodars
#' @param data Long-format data: one row per species x site x habitat.
#' @param lon,lat Names of the longitude and latitude columns (decimal
#'   degrees; one location per site).
#' @param radius_km Grouping radius in km: one number, or a named vector with
#'   one radius per species.
#' @param linkage `"complete"` or `"single"`.
#' @return A tibble in long format with the original column names: the
#'   species (if given), site (the group name), habitat, abundance and effort
#'   (named `"effort"` if there was no effort column), plus `lon` and `lat`
#'   (mean of the sites) and `n_sites` (number of sites merged). Every sampled
#'   group x habitat appears for every species, with 0 where the species was
#'   not recorded.
#' @seealso [radius_sensitivity()], [check_coverage()]
#' @examples
#' set.seed(1)
#' d <- simulate_isodars(1.5, 1.5, 5, 5, n = 12, format = "long")
#' sites <- unique(d$site)
#' xy <- data.frame(site = sites, lon = stats::runif(12, 0, 0.05),
#'                  lat = stats::runif(12, 0, 0.05))
#' d <- merge(d[c("site", "habitat", "abundance")], xy)
#' g <- group_sites(d, site = "site", lon = "lon", lat = "lat", radius_km = 2,
#'                  habitat = "habitat", abundance = "abundance")
#' head(g)
#' @export
group_sites <- function(data, site, lon, lat, radius_km, linkage = "complete",
                        species = NULL, habitat, abundance, effort = NULL) {
  linkage <- match.arg(linkage, c("complete", "single"))
  check_col_arg(lon, "lon", data)
  check_col_arg(lat, "lat", data)
  warn_unused_effort(data, effort)
  prep <- prepare_long(data, species, site, habitat, abundance, effort)
  habs <- prep$habitats
  # site coordinates
  st <- as.character(data[[site]])
  xy <- data.frame(lon = data[[lon]], lat = data[[lat]])
  if (!is.numeric(xy$lon) || !is.numeric(xy$lat) || anyNA(xy)) {
    stop("'lon' and 'lat' must be numeric without missing values.", call. = FALSE)
  }
  sites <- unique(st)
  first <- match(sites, st)
  coords <- xy[first, ]
  same <- abs(xy$lon - coords$lon[match(st, sites)]) < 1e-9 &
    abs(xy$lat - coords$lat[match(st, sites)]) < 1e-9
  if (!all(same)) stop("Each site must have a single location.", call. = FALSE)
  sp_names <- names(prep$species)
  if (!is.numeric(radius_km) || any(!is.finite(radius_km)) || any(radius_km < 0)) {
    stop("'radius_km' must be non-negative numbers.", call. = FALSE)
  }
  if (length(radius_km) == 1 && is.null(names(radius_km))) {
    radius_km <- rep(radius_km, length(sp_names))
  } else if (!prep$species_given || is.null(names(radius_km)) ||
             !all(sp_names %in% names(radius_km))) {
    stop("'radius_km' must be one number or a named vector with a radius for ",
         "every species.", call. = FALSE)
  } else {
    radius_km <- radius_km[sp_names]
  }
  # effort of every site x habitat, from all rows
  E <- matrix(NA_real_, length(sites), length(habs), dimnames = list(sites, habs))
  eff <- if (is.null(effort)) rep(1, nrow(data)) else data[[effort]]
  E[cbind(match(st, sites), match(as.character(data[[habitat]]), habs))] <- eff
  dist <- stats::as.dist(great_circle_km(coords$lon, coords$lat))
  tree <- if (length(sites) > 1) stats::hclust(dist, method = linkage) else NULL
  groups_at <- function(r) {
    if (is.null(tree)) return(1L)
    stats::cutree(tree, h = r + 1e-9)
  }
  cache <- list()
  rows <- lapply(seq_along(sp_names), function(i) {
    s <- sp_names[i]
    key <- format(radius_km[i])
    if (is.null(cache[[key]])) cache[[key]] <<- groups_at(radius_km[i])
    g <- cache[[key]]
    # counts of this species at every sampled site x habitat
    C <- E
    C[!is.na(C)] <- 0
    sc <- prep$species[[i]]$counts
    C[rownames(sc), ][!is.na(sc)] <- sc[!is.na(sc)]
    lab <- vapply(split(sites, g), paste, character(1), collapse = "+")
    out <- list()
    for (k in seq_along(lab)) {
      members <- which(g == as.integer(names(lab)[k]))
      for (h in seq_along(habs)) {
        e <- E[members, h]
        if (all(is.na(e))) next
        out[[length(out) + 1]] <- data.frame(
          species = s, site = lab[[k]], habitat = habs[h],
          abundance = sum(C[members, h], na.rm = TRUE), effort = sum(e, na.rm = TRUE),
          lon = mean(coords$lon[members]), lat = mean(coords$lat[members]),
          n_sites = length(members), stringsAsFactors = FALSE)
      }
    }
    do.call(rbind, out)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!prep$species_given) out$species <- NULL
  nm <- c(species = species, site = site, habitat = habitat, abundance = abundance,
          effort = if (is.null(effort)) "effort" else effort, lon = lon, lat = lat)
  names(out)[match(names(nm), names(out))] <- nm
  tibble::as_tibble(out)
}

great_circle_km <- function(lon, lat) {
  rad <- pi / 180
  p <- lat * rad
  l <- lon * rad
  dp <- outer(p, p, "-")
  dl <- outer(l, l, "-")
  a <- sin(dp / 2)^2 + outer(cos(p), cos(p)) * sin(dl / 2)^2
  a[a > 1] <- 1
  2 * 6371.0088 * asin(sqrt(a))
}

#' Sensitivity of niche breadth to the grouping radius
#'
#' Groups sites with [group_sites()] at each radius and estimates niche
#' breadth with [niche_breadth()], showing how the number of sampling units,
#' their completeness and the estimates change with the radius.
#'
#' @inheritParams group_sites
#' @param radii_km Radii (km) to compare.
#' @param ... Passed to [niche_breadth()] (for example `dim`, `n_boot`; the
#'   default here is `n_boot = 0`).
#' @return A tibble with one row per radius x species: `radius_km`, `species`
#'   (if given), `n_units` (sampling units with the species), `complete_share`
#'   (share of those units in which every habitat was sampled), `nb_di`,
#'   `spillover_abundance`, `nb_dd` and `flags`.
#' @seealso [group_sites()]
#' @examples
#' set.seed(1)
#' d <- simulate_isodars(1.5, 1.5, 5, 5, n = 20, p_missing = 0.3,
#'                       format = "long")
#' xy <- data.frame(site = unique(d$site), lon = stats::runif(20, 0, 0.1),
#'                  lat = stats::runif(20, 0, 0.1))
#' d <- merge(d, xy)
#' radius_sensitivity(d, site = "site", lon = "lon", lat = "lat",
#'                    habitat = "habitat", abundance = "abundance",
#'                    radii_km = c(0, 2, 5), diagnostics = list(n_null = 0))
#' @export
radius_sensitivity <- function(data, species = NULL, site, habitat, abundance,
                               effort = NULL, lon, lat,
                               radii_km = c(0, 0.5, 1, 2, 5),
                               linkage = "complete", ...) {
  dots <- list(...)
  if (is.null(dots$n_boot)) dots$n_boot <- 0
  out <- lapply(radii_km, function(r) {
    g <- group_sites(data, site = site, lon = lon, lat = lat, radius_km = r,
                     linkage = linkage, species = species, habitat = habitat,
                     abundance = abundance, effort = effort)
    nb <- do.call(niche_breadth, c(list(g, species = species, site = site,
                                        habitat = habitat, abundance = abundance,
                                        effort = if (is.null(effort)) "effort" else effort),
                                   dots))
    fits <- attr(nb, "fits")
    res <- tibble::tibble(
      radius_km = r, species = fits$species, n_units = fits$n_sites,
      complete_share = vapply(fits$counts, function(C) {
        if (nrow(C) == 0) NA_real_ else mean(stats::complete.cases(C))
      }, numeric(1)),
      nb_di = nb$nb_di, spillover_abundance = nb$spillover_abundance,
      nb_dd = nb$nb_dd, flags = nb$flags)
    if (is.null(species)) res$species <- NULL
    res
  })
  tibble::as_tibble(do.call(rbind, out))
}
