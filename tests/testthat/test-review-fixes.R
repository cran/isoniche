# Fixes after the pre-submission review

test_that("an unused effort column triggers a warning", {
  comm <- matrix(c(5, 0, 2, 3, 1, 0, 0, 4, 6), nrow = 3, byrow = TRUE,
                 dimnames = list(NULL, c("sp_a", "sp_b", "sp_c")))
  samples <- data.frame(site = c("s1", "s1", "s2"),
                        habitat = c("reef", "lagoon", "reef"),
                        area = c(10, 20, 10))
  long <- community_to_long(comm, samples, site = "site", habitat = "habitat",
                            effort = "area")
  expect_equal(attr(long, "effort_column"), "area")
  sub <- long[long$species != "sp_c", ]
  # the toy data cannot be fitted either; collect all warnings
  msgs <- character(0)
  withCallingHandlers(
    fit_isodars(sub, species = "species", site = "site", habitat = "habitat",
                abundance = "abundance"),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("effort column ('area')", msgs, fixed = TRUE)))
  d <- moorea_long(ref_species)
  d$effort <- 1
  expect_warning(
    niche_breadth(d, species = "species", site = "site", habitat = "habitat",
                  abundance = "abundance", n_boot = 0,
                  diagnostics = list(n_null = 0)),
    "effort = \"effort\"")
  expect_no_warning(
    niche_breadth(d, species = "species", site = "site", habitat = "habitat",
                  abundance = "abundance", effort = "effort", n_boot = 0,
                  diagnostics = list(n_null = 0)))
})

test_that("reference_effort without effort is an error", {
  expect_error(
    niche_breadth(moorea_long(ref_species), species = "species", site = "site",
                  habitat = "habitat", abundance = "abundance",
                  reference_effort = 100, n_boot = 0),
    "needs an 'effort' column")
})

test_that("printed results and curves give the unit of abundance", {
  d <- moorea_long(ref_species)
  nb <- niche_breadth(d, species = "species", site = "site",
                      habitat = "habitat", abundance = "abundance", n_boot = 0,
                      diagnostics = list(n_null = 0))
  expect_output(print(nb), "per site")
  d$area <- 50
  nb2 <- niche_breadth(d, species = "species", site = "site",
                       habitat = "habitat", abundance = "abundance",
                       effort = "area", reference_effort = 100, n_boot = 0,
                       diagnostics = list(n_null = 0))
  expect_output(print(nb2), "reference effort of 100")
  cv <- iai_curve(nb2)
  expect_equal(curve_xlab(attr(cv, "unit")),
               "Total abundance (individuals per reference effort of 100)")
  expect_true(all(is.finite(attr(cv, "points")$max_observed)))
  # densities are per 100 units, so totals double relative to counts per 50
  expect_equal(attr(cv, "points")$max_observed,
               2 * attr(iai_curve(nb), "points")$max_observed)
})

test_that("plot() of a niche_breadth result draws its curves", {
  nb <- niche_breadth(moorea_long(ref_species), species = "species",
                      site = "site", habitat = "habitat",
                      abundance = "abundance", n_boot = 0,
                      diagnostics = list(n_null = 0))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  cv <- plot(nb)
  expect_s3_class(cv, "iai_curve")
})

test_that("tied habitats do not give spurious reversals", {
  nb <- niche_breadth(moorea_long(c("Bothus mancus", "Platax orbicularis",
                                    "Acanthurus nigrofuscus")),
                      species = "species", site = "site", habitat = "habitat",
                      abundance = "abundance", n_boot = 0,
                      diagnostics = list(n_null = 0))
  n_rev <- vapply(nb$reversals, nrow, integer(1))
  names(n_rev) <- nb$species
  expect_equal(unname(n_rev[c("Bothus mancus", "Platax orbicularis")]), c(0L, 0L))
  expect_equal(unname(n_rev["Acanthurus nigrofuscus"]), 3L)
})

test_that("a seed gives per-species results and restores the RNG state", {
  d <- moorea_long(ref_species)
  set.seed(123)
  before <- .Random.seed
  a <- niche_breadth(d, species = "species", site = "site", habitat = "habitat",
                     abundance = "abundance", n_boot = 0, seed = 1)
  expect_identical(.Random.seed, before)
  b <- niche_breadth(d, species = "species", site = "site", habitat = "habitat",
                     abundance = "abundance", n_boot = 20, seed = 1)
  expect_equal(a$intercept_p, b$intercept_p)
  b2 <- niche_breadth(d, species = "species", site = "site", habitat = "habitat",
                      abundance = "abundance", n_boot = 20, seed = 1)
  expect_equal(b$nb_dd_lower, b2$nb_dd_lower)
  cv <- iai_curve(b, n_boot = 5, seed = 2)
  expect_identical(.Random.seed, before)
})

test_that("pairwise allocation is the exact constrained optimum", {
  d <- moorea_long("Chaetodon ulietensis")
  f <- fit_isodars(d, site = "site", habitat = "habitat",
                   abundance = "abundance")$fit[[1]]
  x <- iai(f, abundances = 5, weights = "1/var")
  # dropping both negative habitats at once gave 0 / 5 / 0
  expect_equal(c(x$barrier_reef, x$fringing_reef, x$outer_slope),
               c(1.916, 3.084, 0), tolerance = 1e-3)
  # no other non-negative allocation summing to 5 fits better
  sys <- pairwise_system(f, "1/var")
  obj <- function(v) sum(sys$w * (sys$A %*% v - sys$b)^2)
  best <- obj(unlist(x[c("barrier_reef", "fringing_reef", "outer_slope")]))
  grid <- expand.grid(p = seq(0, 1, 0.01), q = seq(0, 1, 0.01))
  grid <- grid[grid$p + grid$q <= 1, ]
  vals <- apply(grid, 1, function(g) obj(5 * c(g[1], g[2], 1 - g[1] - g[2])))
  expect_true(best <= min(vals) + 1e-8)
})

test_that("a missing column names the argument", {
  expect_error(
    check_coverage(moorea_fish, species = "species", site = "site",
                   habitat = "Habitat", abundance = "abundance"),
    "Column 'Habitat' (argument 'habitat') not found", fixed = TRUE)
})
