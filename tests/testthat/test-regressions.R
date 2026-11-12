# One test per bug fixed relative to CRAN 1.0.0 (see NEWS.md)

test_that("bug 1: the n-dimensional isodar is a line (first principal axis)", {
  set.seed(2)
  sim <- simulate_isodars(1.5, 3, 2, 0, noise = 1, n = 40)
  habs <- c("hab1", "hab2", "hab3")
  ln <- fit_isodars(sim, hab_cols = habs, dim = "ndim")
  X <- as.matrix(sim[habs])
  e <- eigen(stats::cov(X), symmetric = TRUE)
  expect_equal(sum(ln$direction^2), 1)
  # largest eigenvector, not the smallest (hyperplane normal)
  expect_equal(abs(sum(ln$direction * e$vectors[, 1])), 1)
  expect_equal(ln$mean, unname(colMeans(X)))

  # allocations lie on the line and sum to the total
  for (total in c(200, 1000)) {
    x <- unlist(iai(ln, abundances = total)[habs])
    expect_equal(sum(x), total)
    d <- unname((x - ln$mean) / ln$direction)
    expect_equal(d[2:3], rep(d[1], 2))
  }

  # far along the line, allocation is proportional to the direction, not
  # pulled toward equal allocation
  far <- unlist(iai(ln, abundances = 1e7)[habs])
  expect_equal(unname(far / sum(far)), ln$direction / sum(ln$direction),
               tolerance = 1e-4)
})

test_that("bug 2: inverse-variance weights are applied once", {
  pw <- fit_isodars(moorea_species("Acanthurus nigrofuscus"),
                    hab_cols = moorea_habs)
  w <- 1 / (pw$sd^2 * (1 + pw$slope^2))
  by_name <- iai(pw, abundances = 1000, weights = "1/var")
  expect_equal(by_name$niche_breadth,
               iai(pw, abundances = 1000, weights = w)$niche_breadth)
  # the 1.0.0 behaviour (weights squared) gives a different answer
  expect_false(isTRUE(all.equal(
    by_name$niche_breadth,
    iai(pw, abundances = 1000, weights = w^2)$niche_breadth
  )))
})

test_that("bug 3: isoderr() SD and SE are the bootstrap SD", {
  set.seed(1)
  sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
  habs <- c("hab1", "hab2", "hab3")
  n_boot <- 40
  sd_out <- isoderr(sim, hab_cols = habs, abundances = 100, n_boot = n_boot,
                    seed = 7)
  se_out <- isoderr(sim, hab_cols = habs, abundances = 100, n_boot = n_boot,
                    output = "SE", seed = 7)
  expect_equal(se_out$uncertainty, sd_out$uncertainty)

  # recompute the replicates by hand with the same random numbers
  set.seed(7)
  X <- as.data.frame(sim[habs])
  reps <- vapply(seq_len(n_boot), function(b) {
    d <- X[sample.int(nrow(X), replace = TRUE), , drop = FALSE]
    iai(fit_isodars(d, hab_cols = habs), abundances = 100)$niche_breadth
  }, numeric(1))
  expect_equal(sd_out$uncertainty, stats::sd(reps))
  expect_equal(sd_out$n_valid, n_boot)
})

test_that("bug 4: fit_isodars() terminates with negative slopes", {
  set.seed(3)
  x <- 0:20
  d <- data.frame(a = x, b = pmax(round(25 - x + stats::rnorm(21, 0, 1)), 0))
  expect_no_error(pw <- fit_isodars(d, hab_cols = c("a", "b")))
  expect_lt(pw$slope, 0)
  expect_false(pw$flipped)
  expect_true(is.finite(pw$intercept))
  expect_no_error(iai(pw, abundances = c(10, 100)))

  # a positive slope with a negative intercept is flipped once
  d2 <- data.frame(a = 3:20, b = 2 * (3:20) - 5 + c(0.3, -0.3))
  pw2 <- fit_isodars(d2, hab_cols = c("a", "b"))
  expect_true(pw2$flipped)
  expect_gte(pw2$intercept, 0)
  expect_equal(c(pw2$habitat_x, pw2$habitat_y), c("b", "a"))
})

test_that("bug 5: density-independent total puts >= min_occupancy everywhere", {
  pw <- fit_isodars(moorea_species("Acanthurus nigrofuscus"),
                    hab_cols = moorea_habs)
  for (m in c(1, 2, 5)) {
    di <- iai(pw, weights = "1/var", min_occupancy = m)
    expect_false(is.na(di$total_abundance))
    expect_gte(min(unlist(di[moorea_habs])), m - 1e-9)
    before <- iai(pw, abundances = di$total_abundance - 1, weights = "1/var")
    expect_lt(min(unlist(before[moorea_habs])), m)
  }
  # no qualifying total: one NA row. In Scarus schlegeli the outer slope
  # never holds 2 individuals (it is abandoned as total abundance rises).
  none <- iai(pw, weights = "1/var", max_search = 5)
  ss <- fit_isodars(moorea_species("Scarus schlegeli"), hab_cols = moorea_habs)
  expect_true(is.na(iai(ss, weights = "1/var", min_occupancy = 2)$niche_breadth))
  expect_equal(nrow(none), 1)
  expect_true(is.na(none$niche_breadth))
  expect_true(is.na(none$total_abundance))
  expect_match(none$note, "no total up to 5")
})

test_that("bug 6: habitat columns must be given", {
  d <- moorea_species("Scarus schlegeli")
  expect_error(fit_isodars(d), "hab_cols")
  expect_error(isoderr(d, n_boot = 2), "hab_cols")
  # n_habitats takes the first columns
  wide <- d[moorea_habs]
  expect_equal(fit_isodars(wide, n_habitats = 3),
               fit_isodars(wide, hab_cols = moorea_habs))
})
