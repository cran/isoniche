# isoniche 1.1.0

New maintainer: Itai Granot. Author order is now Granot, Dubiner, Belmaker;
Itai Granot and Shahar Dubiner contributed equally.

## Bug fixes

These change results relative to 1.0.0.

* The n-dimensional isodar (`fit_isodars(dim = "ndim")`) is now a line (the
  first principal axis), not a hyperplane (the smallest eigenvector) combined
  with "closest to equal allocation". The old approach biased niche breadth
  toward generalism.
* Inverse-variance weights (`weights = "1/var"`) were applied twice
  (effectively 1 / sd^4). They are now 1 / (sd^2 (1 + slope^2)), applied once.
* `isoderr(output = "SE")` divided the bootstrap standard deviation by
  `sqrt(n_boot)`. The bootstrap SD is itself the standard error, so "SD" and
  "SE" now both return it.
* `fit_isodars()` no longer recurses forever when an isodar has a negative
  slope and a negative intercept.
* Density-independent niche breadth (`iai(abundances = NULL)`) is now evaluated
  at the first integer total with at least `min_occupancy` (default 1)
  individual in every habitat. Previously the threshold was > 1e-9
  individuals. If no total qualifies, one row of `NA` is returned.
* `hab_cols = NULL` silently used every column of `data`. Now `hab_cols` (or
  `n_habitats`) is required.

## New features

* A new user-facing layer works on one or many species at once:
  * `fit_isodars()` accepts long-format data (one row per species x site x
    habitat, columns named by `species`, `site`, `habitat`, `abundance` and
    optionally `effort`) and returns one fit per species (`isodar_fits`).
    Counts are converted to densities per `reference_effort` (default: the
    median effort). The wide-format interface is unchanged.
  * `niche_breadth()` returns, per species, density-independent and
    density-dependent niche breadth with site-bootstrap confidence intervals,
    the spillover abundance (the total at which the least-favoured used
    habitat reaches one individual), preference-reversal abundances,
    abandoned habitats and diagnostic flags as columns. Printing it gives
    one summary line and the unit of abundance. Density-dependent niche
    breadth is by default the limit at infinite abundance
    (`dependent_at = "asymptote"`), computed directly. The flags include a
    test of whether isodar intercepts differ from sampling noise
    (multinomial null model, as in the paper). With `seed`, each species
    gets its own random streams, and the user's random-number state is
    restored afterwards.
  * `summary()` of a `niche_breadth()` result counts the species with
    estimates and without warning flags, gives quartiles of the estimates
    and, for each flag raised, how many species raised it, what it means and
    what to do about it.
  * `iai_curve()` and its `plot()` method give niche breadth as a function of
    total abundance, with bootstrap bands, each species' density-independent
    value and asymptote, preference reversals, and the totals at which two
    species' curves cross. The curve is drawn thinner beyond the largest
    observed site total, the x-axis gives the unit of abundance, and the
    legend avoids the curves. `plot()` of a `niche_breadth()` result plots
    its curves.
  * `community_to_long()` converts a community matrix (samples x species) and
    a table of sample attributes into the long format, pooling samples from
    the same site and habitat and summing their effort.
  * `check_coverage()` summarises sampling coverage per species and habitat;
    `group_sites()` merges nearby sites (complete linkage on great-circle
    distance, summing counts and efforts); `radius_sensitivity()` shows how
    results change with the grouping radius.
  * `fit_isodars()` and `niche_breadth()` warn when the data have an effort
    column that is not used, and `reference_effort` without an `effort`
    column is an error.
* `simulate_isodars()` simulates any number of habitats (`slopes`,
  `intercepts`), optional missing habitats (`p_missing`) and long-format
  output; the three-habitat default gives the same numbers as before.
* New data set `moorea_fish`: reef fishes of Moorea in three habitats (the
  paper's Table S2).
* New vignettes: `vignette("isoniche")` (get started),
  `vignette("data-preparation")`, `vignette("diagnostics")` and
  `vignette("method")`; `citation("isoniche")`; a pkgdown website.
* Habitats a species never uses are handled. `fit_isodars()` (and
  `isoderr()`) gain `absent = "auto"`, which declares a habitat absent when
  the species has no more than `absent_max_count` (default 0) records in it
  across at least `absent_min_sites` (default 5) sites where it was present
  in another habitat. Absent habitats are excluded from the isodar fit, fixed
  at 0, and counted as 0 in the inequality, so they lower niche breadth.
  Previously such habitats broke the fit and `iai()` returned `NA`. Habitats
  with no records but fewer informative sites give a warning and are kept.
  Use `absent = "none"` for the old behaviour or name habitats directly.
  A species that uses a single habitat has niche breadth 0.
  `isoderr()` decides the absent set once, on the full data. Results for
  species that use every habitat are unchanged.
* Scaled niche breadth is kept within [0, 1] (rounding error could give
  values such as -1e-16).
* New `pairwise_niche()` summarises pairwise isodars into niche breadth
  indices (paper, revised Appendix S1). By default (`scale = "log"`) the
  indices do not depend on which habitat of each pair is on the x-axis:
  slope deviations are `log(b)^2` and intercept deviations use the
  perpendicular distance of the isodar from the origin,
  `abs(a) / sqrt(1 + b^2)`; it returns the inverse root-mean-square
  deviations and the weighted means of `abs(log(b))` and of that distance,
  with optional weights (`"1/var"`, the inverse perpendicular residual
  variance, or numeric). Isodars with non-positive slopes are left out of
  the slope indices with a warning. `$cycles` checks consistency: the
  product of the slopes around each cycle of three habitats, which should be
  1, flagged beyond `cycle_tolerance`. `scale = "linear"` gives the earlier
  indices (mean intercept and slope, deviations `(b - 1)^2` and `a^2`);
  `log_slope` is a deprecated alias.
* Pairwise isodars use every site where both habitats of a pair were recorded,
  so sites that lack some habitats still contribute.
* `fit_isodars(dim = "ndim", ndim_method = "latent")` estimates the line from
  incomplete sites (every site with at least two habitats) as a one-factor
  model fitted with glmmTMB, and checks that the habitats are connected through
  shared sites. `ndim_method = "auto"` picks `"pca"` when all sites are
  complete and `"latent"` otherwise.
* `fit_isodars()` returns classed tibbles (`isodars_pairwise`, `isodar_line`).
* `iai()` allocates abundance among habitats by weighted least squares on the
  pairwise isodar equations, constrained to the total and to non-negative
  abundances, solved exactly by an active-set method. (Setting negative
  habitats to 0 and re-solving, as first implemented, could miss the optimum
  when a habitat dropped early should have been kept. For the 325 species of
  `moorea_fish`, no density-independent value or spillover abundance
  changes, and no value in the paper's reference table; density-dependent
  niche breadth rises slightly for 3 species at the asymptote, by at most
  0.014, and for 1 species at a total of 1000, by 0.002.) When the isodars do not
  determine the allocation (for example, a habitat without records was
  kept), `iai()` says so instead of failing with "singular matrix".
* `isoderr(output = "CI")` returns percentile bootstrap intervals.
* `iai()` gives an informative error when `abundances` contains `NA`,
  non-finite or negative values.
* Unit tests cover the paper's reference values, a simulation with a known
  answer, each bug fix above, input validation and the latent method.

# isoniche 1.0.0

* First CRAN release (maintainer Shahar Dubiner).
