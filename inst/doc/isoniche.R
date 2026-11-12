## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.5)

## ----setup--------------------------------------------------------------------
library(isoniche)

## ----data---------------------------------------------------------------------
head(moorea_fish)

## ----community----------------------------------------------------------------
comm <- matrix(c(12, 0, 3,
                  8, 2, 0,
                  0, 5, 9),
               nrow = 3, byrow = TRUE,
               dimnames = list(NULL, c("sp_a", "sp_b", "sp_c")))
samples <- data.frame(site = c("north", "north", "south"),
                      habitat = c("reef", "lagoon", "reef"))
community_to_long(comm, samples, site = "site", habitat = "habitat")

## ----coverage-----------------------------------------------------------------
cov <- check_coverage(moorea_fish, species = "species", site = "site",
                      habitat = "habitat", abundance = "abundance")
cov[cov$species == "Chaetodon ornatissimus", ]
table(cov$status)

## ----nb-----------------------------------------------------------------------
species <- c("Balistapus undulatus", "Chaetodon ornatissimus",
             "Chaetodon ulietensis", "Ctenochaetus binotatus",
             "Dascyllus flavicaudus", "Forcipiger longirostris",
             "Halichoeres marginatus", "Heniochus chrysostomus",
             "Parupeneus multifasciatus", "Sargocentron microstoma",
             "Scarus schlegeli", "Zebrasoma scopas")
fish <- moorea_fish[moorea_fish$species %in% species, ]

nb <- niche_breadth(fish, species = "species", site = "site",
                    habitat = "habitat", abundance = "abundance",
                    n_boot = 99, seed = 1)
nb

## ----nb-table-----------------------------------------------------------------
res <- as.data.frame(nb[c("species", "nb_di", "nb_di_lower", "nb_di_upper",
                          "spillover_abundance", "nb_dd", "nb_dd_lower",
                          "nb_dd_upper")])
res[-1] <- round(res[-1], 2)
res[order(-res$nb_di), ]

## ----summary------------------------------------------------------------------
summary(nb)

## ----clean--------------------------------------------------------------------
s <- summary(nb)
clean <- nb[nb$species %in% s$clean, ]
nrow(clean)

## ----curve, fig.height = 5----------------------------------------------------
curves <- iai_curve(nb, n_boot = 49, seed = 1)
plot(curves, species = c("Chaetodon ornatissimus", "Forcipiger longirostris",
                         "Scarus schlegeli", "Zebrasoma scopas"))

## ----compare, fig.width = 5, fig.height = 5-----------------------------------
plot(nb$nb_di, nb$nb_dd, pch = 19, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "DI niche breadth", ylab = "DD niche breadth")
abline(0, 1, lty = 2)
cor(nb$nb_di, nb$nb_dd, method = "spearman")

## ----crossings----------------------------------------------------------------
head(attr(curves, "crossings"))

