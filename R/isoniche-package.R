#' isoniche: abundance-independent habitat niche breadth
#'
#' Estimates the habitat niche breadth of species from counts in several
#' habitats, corrected for abundance with isodars (Granot, Dubiner & Belmaker
#' 2026). For each species it reports density-independent (DI) niche breadth,
#' how widely the species uses habitats when rare; density-dependent (DD)
#' niche breadth, how widely it uses them when abundant; and the curve of
#' niche breadth against abundance between the two.
#'
#' @section A typical analysis:
#' 1. Put the data in long format, one row per species x site x habitat
#'    ([community_to_long()] converts a community matrix).
#' 2. Check sampling coverage with [check_coverage()]; group nearby sites with
#'    [group_sites()] if sites are small or incomplete.
#' 3. Estimate niche breadth with [niche_breadth()] and check the diagnostics
#'    with `summary()`.
#' 4. Plot niche breadth against abundance with [iai_curve()] and
#'    [plot.iai_curve()].
#'
#' Start with `vignette("isoniche")`; see also `vignette("data-preparation")`,
#' `vignette("diagnostics")` and `vignette("method")`.
#'
#' @section Lower-level functions:
#' [fit_isodars()] fits pairwise or n-dimensional isodars; [iai()] computes
#' niche breadth from fitted isodars at given total abundances; [isoderr()]
#' bootstraps it for one species; [pairwise_niche()] gives indices from the
#' isodar parameters; [simulate_isodars()] simulates data from known isodars.
#'
#' @references Granot, I., Dubiner, S. & Belmaker, J. (2026). Why
#'   abundance-based indices of niche breadth are biased, and what can be done
#'   to improve them. *Ecology Letters* (in revision).
#' @keywords internal
"_PACKAGE"

#' @importFrom stats sd
NULL
