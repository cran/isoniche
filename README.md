# isoniche

<!-- badges: start -->
[![R-CMD-check](https://github.com/itai811/isoniche/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/itai811/isoniche/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/isoniche)](https://CRAN.R-project.org/package=isoniche)
<!-- badges: end -->

**isoniche** estimates the habitat niche breadth of species from counts in
several habitats, corrected for abundance.

Classic indices of niche breadth (Levins' B, 1 - Gini or Shannon diversity
of habitat use) make rare species look like specialists: a rare species
keeps its few individuals in its preferred habitat, while the same species,
when abundant, spills into other habitats. isoniche uses *isodars*, the
relationships between a species' densities in pairs of habitats, to predict
where individuals settle at any total abundance. For each species it
reports:

* **density-independent niche breadth**: how widely the species uses
  habitats when rare;
* **density-dependent niche breadth**: how widely it uses them when
  abundant;
* the **niche breadth curve** between the two, with bootstrap confidence
  intervals and diagnostics that say when an estimate can be trusted.

The values can be compared among species that differ in abundance.

## Installation

```r
# install.packages("remotes")
remotes::install_github("itai811/isoniche", build_vignettes = TRUE)
```

## What you need

Counts of each species in two or more habitats at many replicates (sites,
or repeated surveys), ideally 10 or more, most of them sampling every
habitat. Data go in long format, one row per species x site x habitat; a
community matrix (samples x species) can be converted with
`community_to_long()`.

## Example

Reef fishes counted in three habitats at 13 sites around Moorea:

```r
library(isoniche)
head(moorea_fish)

# 1. how well was each species sampled in each habitat?
cov <- check_coverage(moorea_fish, species = "species", site = "site",
                      habitat = "habitat", abundance = "abundance")

# 2. niche breadth of every species, with bootstrap confidence intervals
fish <- moorea_fish[moorea_fish$species %in% c("Chaetodon ornatissimus",
  "Forcipiger longirostris", "Scarus schlegeli", "Zebrasoma scopas"), ]
nb <- niche_breadth(fish, species = "species", site = "site",
                    habitat = "habitat", abundance = "abundance",
                    n_boot = 199, seed = 1)
nb[c("species", "nb_di", "nb_di_lower", "nb_di_upper", "nb_dd")]

# 3. which estimates can be trusted, and what to do about the others
summary(nb)

# 4. niche breadth against abundance
plot(iai_curve(nb, n_boot = 99, seed = 1))
```

## Learn more

Documentation website: <https://itai811.github.io/isoniche/>


* `vignette("isoniche")`: get started, a complete analysis;
* `vignette("data-preparation")`: data formats, sampling effort, choosing
  habitats, absent habitats, grouping sites;
* `vignette("diagnostics")`: when to trust an estimate;
* `vignette("method")`: how the method works.

## Citation

Granot, I., Dubiner, S. & Belmaker, J. Why abundance-based indices of niche
breadth are biased, and what can be done to improve them. *Ecology Letters*
(in revision). Itai Granot and Shahar Dubiner contributed equally.

See `citation("isoniche")`.

## Version 1.1.0

Version 1.1.0 rewrites the core of the package, adds the functions above
and fixes several bugs of 1.0.0 (on CRAN); results can differ from 1.0.0.
See [NEWS.md](NEWS.md).

## License

GPL-3
