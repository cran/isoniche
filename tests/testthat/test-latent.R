test_that("the latent method recovers the line from incomplete sites", {
  skip_if_not_installed("glmmTMB")
  set.seed(2)
  sim <- sim_incomplete()
  d <- sim$data
  expect_gt(mean(is.na(d)), 0.4)
  expect_true(all(rowSums(!is.na(d)) >= 2))

  ln <- fit_isodars(d, hab_cols = paste0("h", 1:4), dim = "ndim")
  expect_equal(attr(ln, "ndim_method"), "latent")
  expect_equal(attr(ln, "n_sites"), nrow(d))
  expect_equal(sum(ln$direction^2), 1, tolerance = 1e-6)
  expect_lt(max(abs(ln$direction - sim$direction)), 0.01)

  # "pca" with only complete sites is available but uses far fewer sites
  pca <- fit_isodars(d, hab_cols = paste0("h", 1:4), dim = "ndim",
                     ndim_method = "pca")
  expect_lt(attr(pca, "n_sites"), attr(ln, "n_sites") / 5)
  expect_no_error(iai(ln, abundances = 100))
})
