# Values reported in Granot, Dubiner & Belmaker (Ecology Letters), Table S2
# data. Pairwise isodars use weights = "1/var".
reference <- data.frame(
  species = c("Acanthurus nigrofuscus", "Scarus schlegeli", "Siganus argenteus"),
  pw_di = c(0.468, 0.143, 0.159),
  pw_di_total = c(133, 30, 42),
  pw_dd = c(0.520, 0.414, 0.011),
  ln_di = c(0.387, 0.142, 0.150),
  ln_di_total = c(128, 31, 44),
  ln_dd = c(0.513, 0.412, 0.016),
  var_explained = c(0.666, 0.848, 0.996),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(reference))) {
  ref <- reference[i, ]
  test_that(paste("pairwise reference values:", ref$species), {
    d <- moorea_species(ref$species)
    pw <- fit_isodars(d, hab_cols = moorea_habs)
    di <- iai(pw, weights = "1/var")
    dd <- iai(pw, abundances = 1000, weights = "1/var")
    expect_equal(di$total_abundance, ref$pw_di_total)
    expect_equal(round(di$niche_breadth, 3), ref$pw_di)
    expect_equal(round(dd$niche_breadth, 3), ref$pw_dd)
  })

  test_that(paste("line reference values:", ref$species), {
    d <- moorea_species(ref$species)
    ln <- fit_isodars(d, hab_cols = moorea_habs, dim = "ndim")
    expect_s3_class(ln, "isodar_line")
    expect_equal(attr(ln, "ndim_method"), "pca")
    di <- iai(ln)
    dd <- iai(ln, abundances = 1000)
    expect_equal(di$total_abundance, ref$ln_di_total)
    expect_equal(round(di$niche_breadth, 3), ref$ln_di)
    expect_equal(round(dd$niche_breadth, 3), ref$ln_dd)
    expect_equal(round(attr(ln, "var_explained"), 3), ref$var_explained)
  })
}

test_that("line IAI recovers the true value in a simulation", {
  set.seed(4)
  sim <- simulate_isodars(2, 2, 0, 0, noise = 0.5, n = 60)
  ln <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"), dim = "ndim")
  # with zero intercepts the true allocation is proportional to 1:2:4,
  # whose scaled 1 - Gini is 4/7
  expect_equal(round(iai(ln, abundances = 100)$niche_breadth, 3), 0.571)
  expect_equal(iai(ln, abundances = 100)$niche_breadth, 4 / 7, tolerance = 0.01)
})
