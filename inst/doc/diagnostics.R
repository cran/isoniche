## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.5)

## ----setup--------------------------------------------------------------------
library(isoniche)

## ----all----------------------------------------------------------------------
nb <- niche_breadth(moorea_fish, species = "species", site = "site",
                    habitat = "habitat", abundance = "abundance",
                    n_boot = 0, seed = 1)
summary(nb)

## ----noise--------------------------------------------------------------------
noisy <- nb[nb$flag_intercept_noise & !nb$flag_few_sites, ]
noisy[c("species", "nb_di", "nb_dd", "intercept_p")]

## ----look---------------------------------------------------------------------
sp <- "Acanthurus nigrofuscus"
nb[nb$species == sp, c("nb_di", "nb_dd", "var_explained", "flags")]
fit <- fit_isodars(moorea_fish[moorea_fish$species == sp, ], site = "site",
                   habitat = "habitat", abundance = "abundance")
fit$fit[[1]]

## ----look-curve, fig.height = 5-----------------------------------------------
plot(iai_curve(fit))
nb$reversals[[which(nb$species == sp)]]

