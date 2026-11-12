fit_long <- function(d, ...) {
  do.call(fit_isodars, c(list(d), long_cols, list(...)))
}

test_that("long format gives the same fits as wide format", {
  fits <- fit_long(moorea_long(ref_species))
  expect_s3_class(fits, "isodar_fits")
  expect_equal(fits$species, ref_species)
  expect_equal(fits$n_sites, c(13L, 13L, 13L))
  expect_equal(attr(fits, "reference_effort"), 1)
  for (i in seq_along(ref_species)) {
    wide <- fit_isodars(moorea_species(ref_species[i]), hab_cols = moorea_habs)
    expect_equal(as.data.frame(fits$fit[[i]]), as.data.frame(wide))
  }
  ln <- fit_long(moorea_long(ref_species[1]), dim = "ndim")
  expect_equal(attr(ln$fit[[1]], "var_explained"),
               attr(fit_isodars(moorea_species(ref_species[1]), hab_cols = moorea_habs,
                                dim = "ndim"), "var_explained"))
  expect_output(print(fits), "Isodars \\(pairwise\\) for 3 species")
})

test_that("missing rows are zeros where the habitat was sampled", {
  d <- moorea_long("Scarus schlegeli")
  sparse <- d[d$abundance > 0, ]
  # sampled cells come from another species' rows
  other <- moorea_long("Acanthurus nigrofuscus")
  both <- rbind(sparse, other)
  fits <- do.call(fit_isodars, c(list(both), long_cols))
  full <- do.call(fit_isodars, c(list(rbind(d, other)), long_cols))
  expect_equal(fits$counts, full$counts)
  expect_equal(fits$fit, full$fit)
})

test_that("habitats missing at a site are NA, and empty sites are dropped", {
  d <- data.frame(site = rep(c("s1", "s2", "s3"), each = 2),
                  habitat = c("a", "b", "a", "c", "b", "c"),
                  abundance = c(3, 0, 0, 0, 4, 5))
  prep <- isoniche:::prepare_long(d, NULL, "site", "habitat", "abundance")
  C <- prep$species[[1]]$counts
  expect_equal(rownames(C), c("s1", "s3"))
  expect_equal(unname(C["s1", ]), c(3, 0, NA))
  expect_equal(unname(C["s3", ]), c(NA, 4, 5))
})

test_that("effort converts counts to densities per reference effort", {
  d <- moorea_long(ref_species)
  d$effort <- ifelse(d$habitat == "outer_slope", 2, 1)
  fits <- do.call(fit_isodars, c(list(d), long_cols,
                                 list(effort = "effort", reference_effort = 1)))
  D <- fits$counts[[1]] / fits$effort[[1]]
  expect_equal(D[, "outer_slope"], fits$counts[[1]][, "outer_slope"] / 2)
  # same effort everywhere: densities scale, slopes do not change
  d$effort <- 4
  f4 <- do.call(fit_isodars, c(list(d), long_cols,
                               list(effort = "effort", reference_effort = 1)))
  f1 <- do.call(fit_isodars, c(list(d), long_cols, list(effort = "effort")))
  expect_equal(attr(f1, "reference_effort"), 4)
  expect_equal(f4$fit[[1]]$slope, f1$fit[[1]]$slope)
  expect_equal(f4$fit[[1]]$intercept, f1$fit[[1]]$intercept / 4)
})

test_that("long format checks its input", {
  d <- moorea_long(ref_species)
  expect_error(do.call(fit_isodars, c(list(d), long_cols, list(hab_cols = "x"))),
               "not both")
  expect_error(fit_isodars(d, species = "species", site = "site",
                           abundance = "abundance"), "'habitat'")
  expect_error(do.call(fit_isodars, c(list(d[, -4]), long_cols)), "not found")
  expect_error(do.call(fit_isodars, c(list(rbind(d, d[1, ])), long_cols)),
               "More than one row")
  neg <- d
  neg$abundance[1] <- -1
  expect_error(do.call(fit_isodars, c(list(neg), long_cols)), "non-negative")
  d$effort <- 1
  expect_error(do.call(fit_isodars, c(list(d), long_cols,
                                      list(effort = "effort", reference_effort = -1))),
               "reference_effort")
  d$effort[1] <- 2
  expect_error(do.call(fit_isodars, c(list(d), long_cols, list(effort = "effort"))),
               "Effort differs")
  expect_error(do.call(fit_isodars, c(list(d[-ncol(d)]), long_cols, list(absent = "reef"))),
               "not among")
})

test_that("species that cannot be fitted give one warning and a note", {
  d <- moorea_long(c(ref_species, "Acanthurus guttatus"))
  expect_warning(fits <- do.call(fit_isodars, c(list(d), long_cols)),
                 "could not be fitted for 1 of 4 species")
  expect_null(fits$fit[[match("Acanthurus guttatus", fits$species)]])
  expect_true(nzchar(fits$note[match("Acanthurus guttatus", fits$species)]))
})
