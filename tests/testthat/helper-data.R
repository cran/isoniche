moorea_habs <- c("barrier_reef", "fringing_reef", "outer_slope")

# One species of moorea_fish (the paper's Table S2) in wide format
moorea_species <- function(species) {
  d <- moorea_fish[moorea_fish$species == species, ]
  sites <- unique(d$site)
  out <- data.frame(Site = sites)
  for (h in moorea_habs) {
    dh <- d[d$habitat == h, ]
    out[[h]] <- dh$abundance[match(sites, dh$site)]
  }
  out
}

moorea_long <- function(species) moorea_fish[moorea_fish$species %in% species, ]

long_cols <- list(species = "species", site = "site", habitat = "habitat",
                  abundance = "abundance")

ref_species <- c("Acanthurus nigrofuscus", "Scarus schlegeli", "Siganus argenteus")

# Line through four habitats with ~45% of habitat-by-site cells missing and
# at least two habitats recorded at every site
sim_incomplete <- function(n = 300, p_miss = 0.57, noise = 0.5) {
  v <- c(1, 2, 3, 4) / sqrt(30)
  mu <- c(2, 5, 8, 11)
  t <- stats::runif(n, 0, 30)
  X <- outer(t, v) + matrix(mu, n, 4, byrow = TRUE) +
    stats::rnorm(n * 4, 0, noise)
  miss <- matrix(stats::runif(n * 4) < p_miss, n, 4)
  for (i in seq_len(n)) {
    while (sum(!miss[i, ]) < 2) miss[i, sample.int(4, 1)] <- FALSE
  }
  X[miss] <- NA
  colnames(X) <- paste0("h", 1:4)
  list(data = as.data.frame(X), direction = v)
}
