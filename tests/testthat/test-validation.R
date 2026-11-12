wide <- data.frame(a = c(1, 4, 6, 9, 12), b = c(2, 5, 9, 11, 15),
                   c = c(0, 2, 3, 6, 7), label = letters[1:5])

test_that("fit_isodars() checks its data", {
  expect_error(fit_isodars(as.matrix(wide[1:3]), hab_cols = c("a", "b")),
               "data frame")
  expect_error(fit_isodars(wide, hab_cols = c("a", "z")), "not found: z")
  expect_error(fit_isodars(wide, hab_cols = "a"), "two habitats")
  expect_error(fit_isodars(wide, hab_cols = c("a", "label")), "numeric")
  expect_error(fit_isodars(wide, n_habitats = 10), "exceeds")
  neg <- wide
  neg$a[2] <- -1
  expect_error(fit_isodars(neg, hab_cols = c("a", "b")), "non-negative")
  expect_error(fit_isodars(wide, hab_cols = c("a", "b"), dim = "plane"))
  expect_error(fit_isodars(wide, hab_cols = c("a", "b"), dim = "ndim",
                           ndim_method = "svd"))
})

test_that("pairs without enough shared sites are skipped with a warning", {
  d <- wide[1:3]
  d$c[1:3] <- NA
  expect_warning(
    expect_warning(pw <- fit_isodars(d, hab_cols = c("a", "b", "c")),
                   "a-c not fitted"),
    "b-c not fitted"
  )
  expect_equal(nrow(pw), 1)
  flat <- wide
  flat$b <- 3
  expect_error(suppressWarnings(fit_isodars(flat, hab_cols = c("a", "b"))),
               "No pairwise isodar")
})

test_that("ndim_method = 'pca' needs enough complete sites", {
  d <- wide[1:3]
  d$c <- c(1, NA, NA, NA, NA)
  expect_error(fit_isodars(d, hab_cols = c("a", "b", "c"), dim = "ndim",
                           ndim_method = "pca"), "complete sites")
})

test_that("the latent method requires connected habitats", {
  skip_if_not_installed("glmmTMB")
  d <- data.frame(a = c(1:6, NA, NA, NA, NA), b = c(2:7, NA, NA, NA, NA),
                  c = c(NA, NA, NA, NA, NA, NA, 1, 3, 5, 7),
                  e = c(NA, NA, NA, NA, NA, NA, 2, 5, 8, 12))
  expect_error(fit_isodars(d, hab_cols = c("a", "b", "c", "e"), dim = "ndim"),
               "not connected")
})

test_that("iai() checks its arguments", {
  pw <- fit_isodars(wide, hab_cols = c("a", "b", "c"))
  ln <- fit_isodars(wide, hab_cols = c("a", "b", "c"), dim = "ndim")
  expect_error(iai(data.frame(x = 1)), "output of fit_isodars")
  expect_error(iai(pw, abundances = 10, weights = c(1, 2)), "weights")
  expect_error(iai(pw, abundances = 10, weights = "1/sd"), "weights")
  expect_error(iai(pw, abundances = 10, dim = "ndim"), "does not match")
  expect_error(iai(ln, abundances = 10, dim = "pairwise"), "does not match")
  expect_error(iai(pw, abundances = 10, method = "shannon"))
  expect_error(iai(pw, abundances = NA), "finite, non-negative")
  expect_error(iai(pw, abundances = -5), "finite, non-negative")
  expect_error(iai(ln, abundances = "100"), "finite, non-negative")
})

test_that("iai() returns one row per total with habitat columns", {
  pw <- fit_isodars(wide, hab_cols = c("a", "b", "c"))
  out <- iai(pw, abundances = c(10, 50, 100))
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("total_abundance", "niche_breadth", "a", "b", "c", "note"))
  expect_equal(out$total_abundance, c(10, 50, 100))
  expect_true(all(out$niche_breadth >= 0 & out$niche_breadth <= 1))
  expect_equal(rowSums(out[c("a", "b", "c")]), c(10, 50, 100))
  # unscaled 1 - Gini ranges from 1/H to 1
  raw <- iai(pw, abundances = 100, scaled = FALSE)$niche_breadth
  expect_equal(iai(pw, abundances = 100)$niche_breadth, (raw - 1 / 3) * 3 / 2)
})

test_that("isoderr() returns intervals", {
  set.seed(1)
  sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
  ci <- isoderr(sim, hab_cols = c("hab1", "hab2", "hab3"),
                abundances = c(50, 500), n_boot = 30, output = "CI", seed = 1)
  expect_named(ci, c("total_abundance", "niche_breadth", "lower", "upper",
                     "n_valid"))
  expect_true(all(ci$lower <= ci$upper))
  expect_error(isoderr(sim, hab_cols = c("hab1", "hab2"), output = "IQR"))
})

test_that("simulate_isodars() checks its parameters", {
  expect_error(simulate_isodars(1, 1, 200, 0), "hab1_max")
  expect_error(simulate_isodars(1, 1, 0, 0, noise = -1))
  set.seed(1)
  s <- simulate_isodars(1, 1.1, 1, 0, noise = 2, n = 10)
  expect_named(s, c("hab1", "hab2", "hab3", "total_abundance"))
  expect_equal(nrow(s), 10)
  expect_true(all(as.matrix(s) >= 0))
  expect_equal(s$total_abundance, s$hab1 + s$hab2 + s$hab3)
})
