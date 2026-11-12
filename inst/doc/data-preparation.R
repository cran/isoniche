## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.5)
has_glmmTMB <- requireNamespace("glmmTMB", quietly = TRUE)

## ----setup--------------------------------------------------------------------
library(isoniche)

## ----community----------------------------------------------------------------
comm <- matrix(c(12, 0, 3,
                  8, 2, 0,
                  5, 1, 1,
                  0, 5, 9,
                  2, 7, 4),
               nrow = 5, byrow = TRUE,
               dimnames = list(NULL, c("sp_a", "sp_b", "sp_c")))
samples <- data.frame(site = c("north", "north", "north", "south", "south"),
                      habitat = c("reef", "reef", "lagoon", "reef", "lagoon"),
                      area_m2 = c(100, 50, 100, 100, 100),
                      lon = c(-149.80, -149.80, -149.80, -149.85, -149.85),
                      lat = c(-17.48, -17.48, -17.48, -17.55, -17.55))
long <- community_to_long(comm, samples, site = "site", habitat = "habitat",
                          effort = "area_m2", keep = c("lon", "lat"))
long

## ----effort-------------------------------------------------------------------
fish <- moorea_fish[moorea_fish$species == "Chaetodon ornatissimus", ]
fish$effort <- 1
unit1 <- niche_breadth(fish, site = "site", habitat = "habitat",
                       abundance = "abundance", effort = "effort",
                       n_boot = 0, diagnostics = list(n_null = 0))
unit2 <- niche_breadth(fish, site = "site", habitat = "habitat",
                       abundance = "abundance", effort = "effort",
                       reference_effort = 2, n_boot = 0,
                       diagnostics = list(n_null = 0))
data.frame(reference_effort = c(1, 2),
           rbind(unit1[c("nb_di", "spillover_abundance", "nb_dd")],
                 unit2[c("nb_di", "spillover_abundance", "nb_dd")]))

## ----absent-------------------------------------------------------------------
cov <- check_coverage(moorea_fish, species = "species", site = "site",
                      habitat = "habitat", abundance = "abundance")
cov[cov$species == "Chromis viridis", ]

## ----latent, eval = has_glmmTMB-----------------------------------------------
set.seed(1)
inc <- simulate_isodars(slopes = c(1.5, 1, 0.8), intercepts = c(5, 0, 2),
                        n = 60, p_missing = 0.3, format = "long")
head(inc)
niche_breadth(inc, site = "site", habitat = "habitat",
              abundance = "abundance", dim = "ndim", n_boot = 0,
              diagnostics = list(n_null = 0))[c("method", "n_sites", "nb_di",
                                                "nb_dd", "flags")]

## ----group--------------------------------------------------------------------
set.seed(2)
d <- simulate_isodars(1.5, 1.5, 5, 5, n = 30, p_missing = 0.3,
                      format = "long")
xy <- data.frame(site = unique(d$site), lon = stats::runif(30, 0, 0.1),
                 lat = stats::runif(30, 0, 0.1))
d <- merge(d, xy)
g <- group_sites(d, site = "site", lon = "lon", lat = "lat", radius_km = 2,
                 habitat = "habitat", abundance = "abundance")
head(g)

## ----radius-------------------------------------------------------------------
radius_sensitivity(d, site = "site", lon = "lon", lat = "lat",
                   habitat = "habitat", abundance = "abundance",
                   radii_km = c(0, 2, 5), diagnostics = list(n_null = 0))

## ----wide---------------------------------------------------------------------
set.seed(1)
sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
head(sim)
niche_breadth(sim, hab_cols = c("hab1", "hab2", "hab3"), n_boot = 0,
              diagnostics = list(n_null = 0))[c("nb_di", "nb_dd")]

