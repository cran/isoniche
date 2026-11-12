## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.5)

## ----setup--------------------------------------------------------------------
library(isoniche)

## ----bias, fig.height = 3.5---------------------------------------------------
set.seed(2)
sim <- simulate_isodars(slope1 = 1, slope2 = 1, int1 = 8, int2 = 8, noise = 1,
                        n = 60)
classic <- apply(as.matrix(sim[1:3]), 1, function(x) {
  x <- sort(x)
  n <- length(x)
  g <- 2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
  (1 - g - 1 / n) * n / (n - 1)
})
plot(sim$total_abundance, classic, log = "x", pch = 19,
     xlab = "Total abundance at the site",
     ylab = "Classic niche breadth (1 - Gini)")

## ----pairwise-----------------------------------------------------------------
three <- c("Acanthurus nigrofuscus", "Scarus schlegeli", "Siganus argenteus")
fish3 <- moorea_fish[moorea_fish$species %in% three, ]
pw <- fit_isodars(fish3, species = "species", site = "site",
                  habitat = "habitat", abundance = "abundance")
pw$fit[[1]]

## ----line---------------------------------------------------------------------
ln <- fit_isodars(fish3, species = "species", site = "site",
                  habitat = "habitat", abundance = "abundance", dim = "ndim")
ln$fit[[3]]

## ----nb-----------------------------------------------------------------------
nb <- niche_breadth(fish3, species = "species", site = "site",
                    habitat = "habitat", abundance = "abundance",
                    n_boot = 199, seed = 1)
nb[c("species", "nb_di", "nb_di_lower", "nb_di_upper", "spillover_abundance",
     "nb_dd", "nb_dd_lower", "nb_dd_upper")]

## ----reversals----------------------------------------------------------------
nb$abandoned
nb$reversals[[1]]

## ----curve, fig.height = 5----------------------------------------------------
curve <- iai_curve(ln, n_boot = 99, seed = 1)
plot(curve, abundances = 1000, crossings = TRUE)

## ----curve-relative, fig.height = 5-------------------------------------------
plot(curve, x_axis = "relative")

## ----appendix-----------------------------------------------------------------
pn <- pairwise_niche(pw$fit[[1]], weights = "1/var")
pn

## ----cycles-------------------------------------------------------------------
pn$cycles

## ----simulate-----------------------------------------------------------------
set.seed(4)
sim <- simulate_isodars(2, 2, 0, 0, noise = 0.5, n = 60)
iai(fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"), dim = "ndim"),
    abundances = 100)

