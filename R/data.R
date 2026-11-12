#' Reef fishes of Moorea in three habitats
#'
#' Counts of 325 reef fish species in three reef habitats at 13 sites around
#' Moorea, French Polynesia: the data of Granot, Dubiner & Belmaker (Table S2).
#' Fish were counted on belt transects by CRIOBE between 2004 and 2018; there
#' is one count per species, site and habitat. Only sites where
#' a species was recorded are included for that species; all three habitats
#' were surveyed at every site, and zeros are explicit.
#'
#' @format A tibble with one row per species x site x habitat:
#' \describe{
#'   \item{species}{Species name.}
#'   \item{site}{Site (13 sites).}
#'   \item{habitat}{`"barrier_reef"`, `"fringing_reef"` or `"outer_slope"`.}
#'   \item{abundance}{Number of individuals counted.}
#' }
#' @source "MCR LTER Coral Reef Long-term Population and Community Dynamics:
#'   Fishes" in BioTIME 2.0: Dornelas, M., Antão, L. H., Bates, A. E. et al.
#'   (2025). BioTIME 2.0: expanding and improving a database of biodiversity
#'   time series. *Global Ecology and Biogeography*, 34(5), e70003.
#'   \doi{10.1111/geb.70003}
#' @examples
#' head(moorea_fish)
#' length(unique(moorea_fish$species))
"moorea_fish"
