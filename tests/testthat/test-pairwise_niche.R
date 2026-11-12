iso <- tibble::tibble(
  habitat_x = c("a", "a", "b"), habitat_y = c("b", "c", "c"),
  intercept = c(2, 0, 4), slope = c(2, 0.5, 1), sd = c(1, 2, 0.5)
)

# Write every isodar the other way round: x = -a/b + (1/b) y
swap_axes <- function(d) {
  out <- d
  out$habitat_x <- d$habitat_y
  out$habitat_y <- d$habitat_x
  out$intercept <- -d$intercept / d$slope
  out$slope <- 1 / d$slope
  out
}

test_that("log-scale indices match the revised Appendix S1 formulas", {
  out <- pairwise_niche(iso)
  expect_s3_class(out, "pairwise_niche")
  ind <- out$indices
  expect_named(ind, c("mean_abs_log_slope", "mean_distance",
                      "slope_niche_breadth", "intercept_niche_breadth",
                      "n_isodars", "n_slope", "note"))
  d <- abs(iso$intercept) / sqrt(1 + iso$slope^2)
  expect_equal(ind$mean_abs_log_slope, mean(abs(log(iso$slope))))
  expect_equal(ind$mean_distance, mean(d))
  expect_equal(ind$slope_niche_breadth, 1 / sqrt(mean(log(iso$slope)^2)))
  expect_equal(ind$intercept_niche_breadth, 1 / sqrt(mean(d^2)))
  expect_equal(ind$n_isodars, 3)
  expect_equal(ind$n_slope, 3)
  expect_equal(out$isodars$distance, d)
})

test_that("log-scale weights are applied, and 1/var uses perpendicular variance", {
  w <- c(1, 2, 3)
  out <- pairwise_niche(iso, weights = w)$indices
  d <- abs(iso$intercept) / sqrt(1 + iso$slope^2)
  expect_equal(out$slope_niche_breadth, 1 / sqrt(sum(w * log(iso$slope)^2) / 6))
  expect_equal(out$intercept_niche_breadth, 1 / sqrt(sum(w * d^2) / 6))
  expect_equal(out$mean_distance, sum(w * d) / 6)
  expect_equal(pairwise_niche(iso, weights = "1/var")$indices,
               pairwise_niche(iso, weights = 1 / iso$sd^2)$indices)
})

test_that("log-scale results do not depend on which habitat is on the x-axis", {
  sw <- swap_axes(iso)
  for (w in list(NULL, "1/var", c(1, 2, 3))) {
    expect_equal(pairwise_niche(sw, weights = w)$indices,
                 pairwise_niche(iso, weights = w)$indices)
    expect_equal(pairwise_niche(sw, weights = w)$cycles,
                 pairwise_niche(iso, weights = w)$cycles)
  }
  # swapping only some isodars
  part <- iso
  part[2, ] <- swap_axes(iso[2, ])
  expect_equal(pairwise_niche(part, weights = "1/var")$indices,
               pairwise_niche(iso, weights = "1/var")$indices)
  # on fitted isodars, whatever fit_isodars() chose
  pw <- fit_isodars(moorea_species("Acanthurus nigrofuscus"), hab_cols = moorea_habs)
  expect_equal(pairwise_niche(swap_axes(pw), weights = "1/var")$indices,
               pairwise_niche(pw, weights = "1/var")$indices)
})

test_that("reciprocal slopes count as equal deviations on the log scale", {
  sym <- tibble::tibble(habitat_x = c("a", "a"), habitat_y = c("b", "c"),
                        intercept = c(0, 0), slope = c(2, 0.5))
  expect_equal(pairwise_niche(sym)$indices$slope_niche_breadth, 1 / log(2))
  expect_equal(pairwise_niche(sym, scale = "linear")$indices$slope_niche_breadth,
               1 / sqrt((1 + 0.25) / 2))
})

test_that("cycles of three habitats report the product of slopes", {
  out <- pairwise_niche(iso)
  cy <- out$cycles
  expect_equal(nrow(cy), 1)
  expect_equal(c(cy$habitat_1, cy$habitat_2, cy$habitat_3), c("a", "b", "c"))
  # b_ab * b_bc * b_ca = 2 * 1 * (1 / 0.5)
  expect_equal(cy$product, 4)
  expect_equal(cy$log_product, log(4))
  expect_true(cy$inconsistent)
  expect_true(out$inconsistent)
  ok <- iso
  ok$slope <- c(2, 1, 0.5)
  out_ok <- pairwise_niche(ok)
  expect_equal(out_ok$cycles$product, 1)
  expect_false(out_ok$inconsistent)
  expect_true(pairwise_niche(ok, cycle_tolerance = 0)$cycles$log_product == 0)
  # a threshold just below or above |log(4)|
  expect_false(pairwise_niche(iso, cycle_tolerance = log(4.1))$inconsistent)
  expect_true(pairwise_niche(iso, cycle_tolerance = log(3.9))$inconsistent)
  # fewer than three habitats: no cycles
  expect_equal(nrow(pairwise_niche(iso[1, ])$cycles), 0)
  expect_output(print(out), "1 inconsistent")
})

test_that("non-positive slopes are left out of the slope index with a warning", {
  pw <- fit_isodars(moorea_species("Scarus schlegeli"), hab_cols = moorea_habs)
  expect_true(any(pw$slope < 0))
  neg <- pw$slope <= 0
  pairs <- paste0(pw$habitat_x[neg], "-", pw$habitat_y[neg])
  expect_warning(out <- pairwise_niche(pw, weights = "1/var"), pairs[1],
                 fixed = TRUE)
  ind <- out$indices
  expect_equal(ind$n_slope, sum(!neg))
  w <- 1 / pw$sd^2
  expect_equal(ind$slope_niche_breadth,
               1 / sqrt(sum(w[!neg] * log(pw$slope[!neg])^2) / sum(w[!neg])))
  d <- abs(pw$intercept) / sqrt(1 + pw$slope^2)
  expect_equal(ind$intercept_niche_breadth, 1 / sqrt(sum(w * d^2) / sum(w)))
  expect_match(ind$note, "non-positive slope")
  # the cycle still uses every isodar
  s <- function(from, to) {
    k <- which(pw$habitat_x == from & pw$habitat_y == to)
    if (length(k) == 1) pw$slope[k] else
      1 / pw$slope[pw$habitat_x == to & pw$habitat_y == from]
  }
  h <- sort(moorea_habs)
  expect_equal(out$cycles$product, s(h[1], h[2]) * s(h[2], h[3]) * s(h[3], h[1]))
})

test_that("a negative product around a cycle is inconsistent", {
  neg <- iso
  neg$slope[1] <- -2
  out <- suppressWarnings(pairwise_niche(neg))
  expect_true(out$cycles$product < 0)
  expect_true(is.na(out$cycles$log_product))
  expect_true(out$inconsistent)
})

test_that("the linear scale keeps the earlier indices", {
  out <- pairwise_niche(iso, scale = "linear")$indices
  expect_named(out, c("mean_intercept", "mean_slope", "slope_niche_breadth",
                      "intercept_niche_breadth", "n_isodars", "note"))
  expect_equal(out$mean_intercept, 2)
  expect_equal(out$mean_slope, 3.5 / 3)
  expect_equal(out$slope_niche_breadth, 1 / sqrt((1 + 0.25 + 0) / 3))
  expect_equal(out$intercept_niche_breadth, 1 / sqrt((4 + 0 + 16) / 3))
  w <- c(1, 2, 3)
  outw <- pairwise_niche(iso, weights = w, scale = "linear")$indices
  expect_equal(outw$mean_intercept, sum(w * iso$intercept) / 6)
  expect_equal(outw$slope_niche_breadth, 1 / sqrt(sum(w * (iso$slope - 1)^2) / 6))
  v <- 1 / (iso$sd^2 * (1 + iso$slope^2))
  expect_equal(pairwise_niche(iso, weights = "1/var", scale = "linear")$indices,
               pairwise_niche(iso, weights = v, scale = "linear")$indices)
  expect_equal(pairwise_niche(iso, weights = c(1, 1, 0),
                              scale = "linear")$indices$mean_slope, 1.25)
})

test_that("the linear scale orients isodars to non-negative intercepts", {
  neg <- iso
  neg$intercept[1] <- -2
  neg$slope[1] <- 2
  o <- pairwise_niche(neg, scale = "linear")$isodars
  expect_equal(c(o$habitat_x[1], o$habitat_y[1]), c("b", "a"))
  expect_equal(o$intercept[1], 1)
  expect_equal(o$slope[1], 0.5)
  expect_true(o$flipped[1])
  bad <- iso
  bad$intercept[1] <- -1
  bad$slope[1] <- -1
  expect_error(pairwise_niche(bad, scale = "linear"), "non-negative intercept")
})

test_that("log_slope is a deprecated alias of scale", {
  expect_warning(a <- pairwise_niche(iso, log_slope = TRUE), "deprecated")
  expect_equal(a$indices, pairwise_niche(iso)$indices)
  expect_warning(b <- pairwise_niche(iso, log_slope = FALSE), "deprecated")
  expect_equal(b$indices, pairwise_niche(iso, scale = "linear")$indices)
  expect_error(pairwise_niche(iso, scale = "log", log_slope = TRUE), "not both")
})

test_that("a perfect generalist gives infinite indices", {
  gen <- tibble::tibble(habitat_x = "a", habitat_y = "b", intercept = 0,
                        slope = 1)
  for (s in c("log", "linear")) {
    out <- pairwise_niche(gen, scale = s)$indices
    expect_equal(out$slope_niche_breadth, Inf)
    expect_equal(out$intercept_niche_breadth, Inf)
  }
})

test_that("pairwise_niche() checks its arguments", {
  ln <- fit_isodars(moorea_species("Scarus schlegeli"), hab_cols = moorea_habs,
                    dim = "ndim")
  expect_error(pairwise_niche(ln), "pairwise isodars")
  expect_error(pairwise_niche(data.frame(x = 1)), "pairwise isodars")
  expect_error(pairwise_niche(iso[0, ]), "no rows")
  expect_error(pairwise_niche(iso, weights = c(1, 2)), "weights")
  expect_error(pairwise_niche(iso, weights = c(1, -1, 1)), "weights")
  expect_error(pairwise_niche(iso[, 1:4], weights = "1/var"), "'sd' column")
  expect_error(pairwise_niche(iso, scale = "sqrt"), "should be one of")
  expect_error(pairwise_niche(iso, cycle_tolerance = -1), "cycle_tolerance")
})
