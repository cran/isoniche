nb_long <- function(d, ...) {
  do.call(niche_breadth, c(list(d), long_cols, list(...)))
}

test_that("niche_breadth() reproduces the reference values", {
  pw <- nb_long(moorea_long(ref_species), n_boot = 0)
  expect_s3_class(pw, "niche_breadth")
  expect_equal(pw$species, ref_species)
  expect_equal(pw$method, rep("pairwise", 3))
  expect_equal(round(pw$nb_di, 3), c(0.468, 0.143, 0.159))
  expect_equal(pw$spillover_abundance, c(133, 30, 42))
  dd <- nb_long(moorea_long(ref_species), n_boot = 0, dependent_at = 1000)
  expect_equal(round(dd$nb_dd, 3), c(0.520, 0.414, 0.011))
  expect_equal(dd$dependent_at, rep(1000, 3))
  ln <- nb_long(moorea_long(ref_species), n_boot = 0, dim = "ndim",
                dependent_at = 1000)
  expect_equal(ln$method, rep("line (pca)", 3))
  expect_equal(round(ln$nb_di, 3), c(0.387, 0.142, 0.150))
  expect_equal(ln$spillover_abundance, c(128, 31, 44))
  expect_equal(round(ln$nb_dd, 3), c(0.513, 0.412, 0.016))
  expect_equal(round(ln$var_explained, 3), c(0.666, 0.848, 0.996))
})

test_that("the spillover abundance equals the iai() search for every species", {
  nb <- nb_long(moorea_fish, n_boot = 0, diagnostics = list(n_null = 0))
  fits <- attr(nb, "fits")
  for (i in which(!is.na(nb$nb_di))) {
    di <- iai(fits$fit[[i]], weights = "1/var")
    expect_equal(nb$spillover_abundance[i], di$total_abundance)
    expect_equal(nb$nb_di[i], di$niche_breadth)
  }
  expect_true(all(nb$flag_fit_failed == is.na(nb$method) |
                    grepl("failed", nb$note)))
})

test_that("the asymptote is the limit of iai() at large totals", {
  nb <- nb_long(moorea_long(ref_species), n_boot = 0)
  fits <- attr(nb, "fits")
  for (i in 1:3) {
    far <- iai(fits$fit[[i]], abundances = 1e9, weights = "1/var")$niche_breadth
    expect_equal(nb$nb_dd[i], far, tolerance = 1e-6)
  }
  expect_equal(nb$dependent_at, rep(Inf, 3))
  mx <- nb_long(moorea_long(ref_species), n_boot = 0, dependent_at = "max_observed")
  expect_true(all(is.finite(mx$dependent_at)))
})

test_that("niche breadth follows the closed form between DI and DD", {
  set.seed(3)
  sim <- simulate_isodars(1.5, 1.5, 3, 3, noise = 1, n = 40)
  nb <- niche_breadth(sim, hab_cols = c("hab1", "hab2", "hab3"), dim = "ndim",
                      n_boot = 0)
  expect_equal(nrow(nb$reversals[[1]]), 0)
  expect_equal(nb$abandoned, "")
  fit <- attr(nb, "fits")$fit[[1]]
  T0 <- nb$spillover_abundance
  for (Tt in c(T0, 2 * T0, 10 * T0, 500)) {
    expect_equal(iai(fit, abundances = Tt)$niche_breadth,
                 nb$nb_dd + (nb$nb_di - nb$nb_dd) * T0 / Tt)
  }
})

test_that("preference reversals and abandoned habitats are reported", {
  nb <- nb_long(moorea_long(ref_species), n_boot = 0)
  expect_equal(nb$abandoned, c("", "outer_slope", "outer_slope"))
  rv <- nb$reversals[[1]]
  expect_gt(nrow(rv), 0)
  fit <- attr(nb, "fits")$fit[[1]]
  x <- iai(fit, abundances = rv$total_abundance[1], weights = "1/var")
  expect_equal(x[[rv$habitat_1[1]]], x[[rv$habitat_2[1]]], tolerance = 1e-6)
  expect_true(all(rv$total_abundance >= nb$spillover_abundance[1]))
  expect_true(nb$flag_reversal[1])
  expect_true(nb$flag_abandoned[2])
})

test_that("bootstrap intervals are reproducible and bracket the estimate", {
  d <- moorea_long(ref_species[2:3])
  a <- nb_long(d, n_boot = 60, seed = 5)
  b <- nb_long(d, n_boot = 60, seed = 5)
  expect_equal(a$nb_dd_lower, b$nb_dd_lower)
  expect_true(all(a$nb_di_lower <= a$nb_di & a$nb_di <= a$nb_di_upper))
  expect_true(all(a$n_boot_valid > 50))
  expect_true(all(a$spillover_lower <= a$spillover_upper))
})

test_that("diagnostic flags are columns and print one summary line", {
  nb <- nb_long(moorea_long(c(ref_species, "Acanthurus guttatus")), n_boot = 0)
  expect_true(all(paste0("flag_", isoniche:::flag_names) %in% names(nb)))
  k <- match("Acanthurus guttatus", nb$species)
  expect_true(nb$flag_fit_failed[k])
  expect_true(is.na(nb$nb_di[k]))
  expect_match(nb$flags[k], "fit_failed")
  expect_output(print(nb), "species flagged \\(not counting")
  strict <- nb_long(moorea_long(ref_species), n_boot = 0,
                    diagnostics = list(min_sites = 20, n_null = 0))
  expect_true(all(strict$flag_few_sites))
  expect_true(all(is.na(strict$intercept_p)))
  expect_error(nb_long(moorea_long(ref_species), diagnostics = list(min_site = 3)),
               "Unknown diagnostics")
})

test_that("the intercept test separates intercepts from sampling noise", {
  set.seed(3)
  n <- 15
  base <- stats::runif(n, 5, 60)
  # no habitat selection: individuals spread in fixed proportions
  noise <- data.frame(t(vapply(round(base * 3), function(N) {
    as.numeric(stats::rmultinom(1, N, c(0.2, 0.3, 0.5)))
  }, numeric(3))))
  names(noise) <- c("a", "b", "c")
  p0 <- niche_breadth(noise, hab_cols = c("a", "b", "c"), n_boot = 0, seed = 1)
  expect_true(p0$flag_intercept_noise) # about 5% of seeds would fail
  # strong intercepts
  sel <- data.frame(a = round(base), b = round(30 + base), c = round(60 + base))
  p1 <- niche_breadth(sel, hab_cols = c("a", "b", "c"), n_boot = 0, seed = 1)
  expect_false(p1$flag_intercept_noise)
  expect_lt(p1$intercept_p, 0.05)
})

test_that("wide data, a single species and argument checks", {
  set.seed(1)
  sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
  nb <- niche_breadth(sim, hab_cols = c("hab1", "hab2", "hab3"), n_boot = 10)
  expect_false("species" %in% names(nb))
  expect_equal(nb$nb_di, iai(fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3")),
                             weights = "1/var")$niche_breadth)
  one <- niche_breadth(moorea_long(ref_species[1])[-1], site = "site",
                       habitat = "habitat", abundance = "abundance", n_boot = 0)
  expect_false("species" %in% names(one))
  expect_equal(nrow(one), 1)
  expect_error(niche_breadth(sim, hab_cols = c("hab1", "hab2"), dependent_at = "max"),
               "dependent_at")
  expect_error(niche_breadth(sim, hab_cols = c("hab1", "hab2"), weights = "1/sd"),
               "weights")
  expect_error(niche_breadth(sim, hab_cols = c("hab1", "hab2"), reference_effort = 2),
               "long-format")
  fits <- do.call(fit_isodars, c(list(moorea_long(ref_species)), long_cols))
  expect_error(niche_breadth(fits, dim = "ndim"), "already fitted")
  expect_error(niche_breadth(fits$fit[[1]]), "Give the data")
  expect_equal(niche_breadth(fits, n_boot = 0)$nb_di,
               nb_long(moorea_long(ref_species), n_boot = 0)$nb_di)
})

test_that("summary() counts flags and explains them", {
  nb <- niche_breadth(moorea_long(ref_species), species = "species",
                      site = "site", habitat = "habitat",
                      abundance = "abundance", n_boot = 0, seed = 1)
  s <- summary(nb)
  expect_s3_class(s, "summary.niche_breadth")
  expect_equal(s$n_species, 3)
  expect_equal(s$n_estimated, 3)
  fl <- as.data.frame(nb)[grep("^flag_", names(nb))]
  for (i in seq_len(nrow(s$flags))) {
    expect_equal(s$flags$n[i], sum(fl[[paste0("flag_", s$flags$flag[i])]]))
  }
  expect_true(all(nzchar(s$flags$advice)))
  ok <- rowSums(fl[paste0("flag_", setdiff(flag_names, c("abandoned", "reversal")))]) == 0
  expect_setequal(s$clean, nb$species[ok])
  expect_equal(s$estimates["nb_di", "median"], stats::median(nb$nb_di))
  expect_output(print(s), "Flags raised")
  expect_output(print(nb), "see summary()", fixed = TRUE)
})

test_that("every flag has an explanation", {
  expect_setequal(flag_info$flag, flag_names)
  expect_true(all(nzchar(flag_info$meaning) & nzchar(flag_info$advice)))
})
