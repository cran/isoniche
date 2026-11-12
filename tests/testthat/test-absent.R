# Two used habitats and one never used (all zeros) at 12 sites
absent_data <- function(n = 12) {
  set.seed(11)
  a <- round(stats::runif(n, 2, 40))
  data.frame(a = a, b = pmax(round(3 + 1.5 * a + stats::rnorm(n, 0, 2)), 0),
             c = 0)
}
habs3 <- c("a", "b", "c")

test_that("a habitat never used is declared absent and fixed at 0", {
  d <- absent_data()
  pw <- expect_silent(fit_isodars(d, hab_cols = habs3))
  expect_equal(attr(pw, "absent"), "c")
  expect_equal(attr(pw, "habitats"), habs3)
  expect_equal(nrow(pw), 1)
  expect_equal(c(pw$habitat_x, pw$habitat_y), c("a", "b"))

  di <- iai(pw)
  expect_false(is.na(di$niche_breadth))
  expect_equal(di$c, 0)
  expect_gte(min(di$a, di$b), 1 - 1e-9)
  expect_match(di$note, "absent habitat\\(s\\) fixed at 0: c")

  # the zero counts in the Gini: same as the two-habitat allocation plus a 0
  dd <- iai(pw, abundances = 500)
  two <- fit_isodars(d, hab_cols = c("a", "b"))
  x <- unlist(iai(two, abundances = 500)[c("a", "b")])
  expect_equal(unname(unlist(dd[c("a", "b")])), unname(x))
  g <- 1 - isoniche:::gini(c(x, 0))
  expect_equal(dd$niche_breadth, (g - 1 / 3) * 3 / 2)
  expect_lt(dd$niche_breadth, iai(two, abundances = 500)$niche_breadth)
})

test_that("absent habitats work for the line", {
  d <- absent_data()
  ln <- fit_isodars(d, hab_cols = habs3, dim = "ndim")
  expect_s3_class(ln, "isodar_line")
  expect_equal(ln$habitat, habs3)
  expect_equal(ln$mean[3], 0)
  expect_equal(ln$direction[3], 0)
  expect_equal(sum(ln$direction^2), 1)
  expect_equal(attr(ln, "ndim_method"), "pca")
  expect_output(print(ln), "Absent habitats \\(fixed at 0\\): c")
  two <- fit_isodars(d, hab_cols = c("a", "b"), dim = "ndim")
  expect_equal(ln$direction[1:2], two$direction)
  expect_equal(attr(ln, "var_explained"), attr(two, "var_explained"))
  dd <- iai(ln, abundances = 500)
  expect_equal(dd$c, 0)
  expect_equal(dd$a + dd$b, 500)
})

test_that("zeros at too few sites are a warning, not absence", {
  d <- absent_data()[1:4, ]
  w <- capture_warnings(pw <- fit_isodars(d, hab_cols = habs3))
  expect_match(w, "too few informative", all = FALSE)
  expect_equal(attr(pw, "absent"), character(0))
  expect_equal(attr(pw, "undersampled"), "c")
  # the threshold is an argument
  pw3 <- expect_silent(fit_isodars(d, hab_cols = habs3, absent_min_sites = 4))
  expect_equal(attr(pw3, "absent"), "c")
  # only sites where the species is present elsewhere count
  d2 <- absent_data()
  d2[1:8, c("a", "b")] <- 0
  w2 <- capture_warnings(fit_isodars(d2, hab_cols = habs3))
  expect_match(w2, "too few informative", all = FALSE)
})

test_that("absent_max_count tolerates stray records", {
  d <- absent_data()
  d$c[3] <- 1
  pw <- fit_isodars(d, hab_cols = habs3)
  expect_equal(attr(pw, "absent"), character(0))
  pw1 <- fit_isodars(d, hab_cols = habs3, absent_max_count = 1)
  expect_equal(attr(pw1, "absent"), "c")
  expect_equal(iai(pw1, abundances = 100)$c, 0)
})

test_that("absent can be switched off or given by name", {
  d <- absent_data()
  w <- capture_warnings(pw <- fit_isodars(d, hab_cols = habs3, absent = "none"))
  expect_match(w, "not fitted")
  expect_equal(attr(pw, "absent"), character(0))
  expect_equal(
    fit_isodars(d, hab_cols = habs3, absent = character(0), dim = "ndim")$direction[3],
    0
  )
  d$c <- d$a %/% 3
  byname <- fit_isodars(d, hab_cols = habs3, absent = "c")
  expect_equal(attr(byname, "absent"), "c")
  expect_equal(nrow(byname), 1)
  expect_error(fit_isodars(d, hab_cols = habs3, absent = "z"), "not among")
  expect_error(fit_isodars(d, hab_cols = habs3, absent = TRUE), "'absent' must")
  expect_error(fit_isodars(d, hab_cols = habs3, absent_min_sites = 0),
               "absent_min_sites")
  expect_error(fit_isodars(d, hab_cols = habs3, absent = habs3), "Every habitat")
})

test_that("a single used habitat gives niche breadth 0", {
  d <- absent_data()
  d$b <- 0
  pw <- fit_isodars(d, hab_cols = habs3)
  expect_equal(attr(pw, "absent"), c("b", "c"))
  expect_equal(nrow(pw), 0)
  di <- iai(pw)
  expect_equal(di$total_abundance, 1)
  expect_equal(di$niche_breadth, 0)
  expect_equal(iai(pw, abundances = 100)$a, 100)
  expect_equal(iai(pw, abundances = 100, weights = "1/var")$niche_breadth, 0)
  ln <- fit_isodars(d, hab_cols = habs3, dim = "ndim")
  expect_equal(attr(ln, "ndim_method"), "single habitat")
  expect_equal(iai(ln, abundances = 100)$niche_breadth, 0)
  expect_error(pairwise_niche(pw), "fewer than two habitats")
})

test_that("habitats abandoned at high density are not absent", {
  pw <- fit_isodars(moorea_species("Scarus schlegeli"), hab_cols = moorea_habs)
  expect_equal(attr(pw, "absent"), character(0))
  expect_true(any(pw$slope < 0))
  for (sp in c("Acanthurus nigrofuscus", "Siganus argenteus")) {
    fit <- fit_isodars(moorea_species(sp), hab_cols = moorea_habs, dim = "ndim")
    expect_equal(attr(fit, "absent"), character(0))
  }
})

test_that("isoderr() keeps the absent set fixed across resamples", {
  d <- absent_data(20)
  out <- isoderr(d, hab_cols = habs3, abundances = 200, n_boot = 30, seed = 1)
  expect_equal(out$n_valid, 30)
  expect_equal(out$niche_breadth,
               iai(fit_isodars(d, hab_cols = habs3), abundances = 200)$niche_breadth)
  expect_true(is.finite(out$uncertainty))
})

test_that("scaled niche breadth stays within [0, 1]", {
  d <- absent_data()
  d$b <- 0
  pw <- fit_isodars(d, hab_cols = habs3)
  nb <- iai(pw, abundances = c(1, 7, 1000))$niche_breadth
  expect_true(all(nb >= 0 & nb <= 1))
})
