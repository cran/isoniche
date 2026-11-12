#' Convert a community matrix to long format
#'
#' Turns a community matrix (one row per sample, one column per species, as
#' used by vegan and most community-ecology software) and a table describing
#' the samples into the long format used by [check_coverage()],
#' [fit_isodars()] and [niche_breadth()]: one row per species x site x
#' habitat.
#'
#' @details
#' Each row of `comm` is one sample (a transect, quadrat, trap, count...) and
#' the matching row of `samples` says at which site and in which habitat it
#' was taken. Samples taken at the same site and habitat are pooled: counts
#' are summed, and so are efforts. Without an `effort` column, every sample
#' has effort 1, so the effort of a site x habitat is its number of samples;
#' sites and habitats sampled unevenly are then still compared as densities
#' per sample.
#'
#' Every sampled site x habitat appears for every species, with 0 where the
#' species was not recorded (which is what makes zeros informative). Species
#' never recorded are dropped.
#'
#' @param comm A matrix or data frame of counts: one row per sample, one
#'   column per species (column names are the species names). Counts must be
#'   non-negative and not missing.
#' @param samples A data frame with one row per row of `comm`, in the same
#'   order, giving the site, habitat and (optionally) effort of each sample.
#' @param site,habitat Names of the site and habitat columns of `samples`.
#' @param effort Optional name of the effort column of `samples` (area,
#'   transect length, trap nights...; positive).
#' @param keep Optional names of other site-level columns of `samples` to
#'   carry over, such as longitude and latitude for [group_sites()]. They must
#'   have one value per site.
#' @return A tibble with columns `species`, the site and habitat columns
#'   (named as in `samples`), `abundance`, `effort` (named as in `samples`, or
#'   `"effort"` = number of samples) and any `keep` columns, ready for
#'   `niche_breadth(x, species = "species", site = ..., habitat = ...,
#'   abundance = "abundance", effort = ...)`. Pass the effort column as
#'   `effort`: without it, every pooled site x habitat counts as equal effort
#'   (`fit_isodars()` and `niche_breadth()` warn when it is left out).
#' @seealso [check_coverage()], [niche_breadth()], `vignette("data-preparation")`
#' @examples
#' comm <- matrix(c(5, 0, 2,
#'                  3, 1, 0,
#'                  0, 4, 6,
#'                  1, 2, 2), nrow = 4, byrow = TRUE,
#'                dimnames = list(NULL, c("sp_a", "sp_b", "sp_c")))
#' samples <- data.frame(site = c("s1", "s1", "s1", "s2"),
#'                       habitat = c("reef", "reef", "lagoon", "reef"))
#' community_to_long(comm, samples, site = "site", habitat = "habitat")
#' @export
community_to_long <- function(comm, samples, site, habitat, effort = NULL,
                              keep = NULL) {
  if (!is.matrix(comm) && !is.data.frame(comm)) {
    stop("'comm' must be a matrix or data frame (samples x species).", call. = FALSE)
  }
  if (is.data.frame(comm)) {
    num <- vapply(comm, is.numeric, logical(1))
    if (!all(num)) {
      stop("Every column of 'comm' must be numeric (non-numeric: ",
           paste(names(comm)[!num], collapse = ", "),
           "); put site and habitat in 'samples'.", call. = FALSE)
    }
    comm <- as.matrix(comm)
  }
  if (!is.numeric(comm)) stop("'comm' must hold numbers.", call. = FALSE)
  sp <- colnames(comm)
  if (is.null(sp) || anyNA(sp) || any(!nzchar(sp)) || anyDuplicated(sp)) {
    stop("'comm' must have unique, non-empty column names (the species).",
         call. = FALSE)
  }
  if (anyNA(comm) || any(!is.finite(comm)) || any(comm < 0)) {
    stop("'comm' must hold finite, non-negative counts without missing values.",
         call. = FALSE)
  }
  if (!is.data.frame(samples)) stop("'samples' must be a data frame.", call. = FALSE)
  if (nrow(samples) != nrow(comm)) {
    stop("'samples' must have one row per row of 'comm' (", nrow(comm),
         " rows; 'samples' has ", nrow(samples), ").", call. = FALSE)
  }
  check_col_arg(site, "site", samples)
  check_col_arg(habitat, "habitat", samples)
  check_col_arg(effort, "effort", samples, required = FALSE)
  if (!is.null(keep)) {
    if (!is.character(keep)) stop("'keep' must be column names.", call. = FALSE)
    miss <- setdiff(keep, names(samples))
    if (length(miss) > 0) {
      stop("Column(s) not found in 'samples': ", paste(miss, collapse = ", "),
           call. = FALSE)
    }
    if (any(keep %in% c(site, habitat, effort, "species", "abundance"))) {
      stop("'keep' must not include the site, habitat, effort, species or ",
           "abundance columns.", call. = FALSE)
    }
  }
  st <- as.character(samples[[site]])
  hb <- samples[[habitat]]
  habs <- if (is.factor(hb)) levels(droplevels(hb)) else sort(unique(as.character(hb)))
  hb <- as.character(hb)
  if (anyNA(st) || anyNA(hb)) {
    stop("Missing values in the site or habitat column.", call. = FALSE)
  }
  eff <- if (is.null(effort)) rep(1, nrow(samples)) else samples[[effort]]
  if (!is.numeric(eff) || any(!is.finite(eff)) || any(eff <= 0)) {
    stop("Column '", effort, "' must hold finite, positive numbers.", call. = FALSE)
  }
  sites <- unique(st)
  cell <- paste(match(st, sites), match(hb, habs))
  cells <- unique(cell)
  first <- match(cells, cell)
  idx <- match(cell, cells)
  counts <- rowsum(comm, idx, reorder = TRUE)
  cell_effort <- as.vector(rowsum(eff, idx, reorder = TRUE))
  keep_vals <- lapply(keep, function(k) {
    v <- samples[[k]]
    ok <- vapply(split(v, st), function(x) length(unique(x)) == 1, logical(1))
    if (!all(ok)) {
      stop("Column '", k, "' has more than one value at site '",
           names(ok)[!ok][1], "'.", call. = FALSE)
    }
    v[match(st[first], st)]
  })
  names(keep_vals) <- keep
  used <- colSums(counts) > 0
  if (!any(used)) stop("No species was recorded in 'comm'.", call. = FALSE)
  counts <- counts[, used, drop = FALSE]
  n_cell <- length(cells)
  ord <- order(match(st[first], sites), match(hb[first], habs))
  out <- data.frame(
    species = rep(colnames(counts), each = n_cell),
    site = rep(st[first][ord], ncol(counts)),
    habitat = rep(hb[first][ord], ncol(counts)),
    abundance = as.vector(counts[ord, , drop = FALSE]),
    effort = rep(cell_effort[ord], ncol(counts)),
    stringsAsFactors = FALSE
  )
  for (k in keep) out[[k]] <- rep(keep_vals[[k]][ord], ncol(counts))
  if (is.factor(samples[[habitat]])) out$habitat <- factor(out$habitat, levels = habs)
  names(out)[2:3] <- c(site, habitat)
  names(out)[5] <- if (is.null(effort)) "effort" else effort
  out <- tibble::as_tibble(out)
  # lets fit_isodars() and niche_breadth() warn if the effort is not used
  attr(out, "effort_column") <- names(out)[5]
  out
}
