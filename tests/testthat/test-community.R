comm <- matrix(c(5, 0, 2,
                 3, 1, 0,
                 0, 4, 6,
                 1, 2, 2,
                 0, 0, 0), nrow = 5, byrow = TRUE,
               dimnames = list(NULL, c("sp_a", "sp_b", "sp_c")))
samples <- data.frame(plot = c("s1", "s1", "s1", "s2", "s2"),
                      hab = c("reef", "reef", "lagoon", "reef", "lagoon"),
                      area = c(10, 20, 10, 10, 5),
                      lon = c(1, 1, 1, 2, 2))

test_that("community_to_long pools samples within site x habitat", {
  out <- community_to_long(comm, samples, site = "plot", habitat = "hab")
  expect_named(out, c("species", "plot", "hab", "abundance", "effort"))
  # 3 species x 4 sampled site x habitat cells, zeros explicit
  expect_equal(nrow(out), 12)
  a <- out[out$species == "sp_a" & out$plot == "s1" & out$hab == "reef", ]
  expect_equal(a$abundance, 8)
  expect_equal(a$effort, 2)
  b <- out[out$species == "sp_b" & out$plot == "s2" & out$hab == "lagoon", ]
  expect_equal(b$abundance, 0)
  expect_equal(b$effort, 1)
  expect_equal(sum(out$abundance), sum(comm))
})

test_that("community_to_long sums effort and keeps site-level columns", {
  out <- community_to_long(comm, samples, site = "plot", habitat = "hab",
                           effort = "area", keep = "lon")
  expect_named(out, c("species", "plot", "hab", "abundance", "area", "lon"))
  a <- out[out$species == "sp_a" & out$plot == "s1" & out$hab == "reef", ]
  expect_equal(a$area, 30)
  expect_equal(unique(out$lon[out$plot == "s2"]), 2)
})

test_that("community_to_long drops species never recorded", {
  cm <- cbind(comm, sp_d = 0)
  out <- community_to_long(cm, samples, site = "plot", habitat = "hab")
  expect_false("sp_d" %in% out$species)
})

test_that("community_to_long output feeds niche_breadth", {
  d <- moorea_long(ref_species)
  d$sample <- paste(d$site, d$habitat)
  cells <- unique(d[c("sample", "site", "habitat")])
  cm <- sapply(ref_species, function(s) {
    ds <- d[d$species == s, ]
    v <- ds$abundance[match(cells$sample, ds$sample)]
    v[is.na(v)] <- 0
    v
  })
  long <- community_to_long(cm, cells, site = "site", habitat = "habitat")
  nb1 <- niche_breadth(long, species = "species", site = "site",
                       habitat = "habitat", abundance = "abundance",
                       effort = "effort", n_boot = 0,
                       diagnostics = list(n_null = 0))
  nb2 <- niche_breadth(d, species = "species", site = "site",
                       habitat = "habitat", abundance = "abundance",
                       n_boot = 0, diagnostics = list(n_null = 0))
  nb2 <- nb2[match(nb1$species, nb2$species), ]
  expect_equal(nb1$nb_di, nb2$nb_di)
  expect_equal(nb1$nb_dd, nb2$nb_dd)
})

test_that("community_to_long validates its input", {
  expect_error(community_to_long(comm, samples[1:3, ], "plot", "hab"),
               "one row per row")
  expect_error(community_to_long(unname(comm), samples, "plot", "hab"),
               "column names")
  bad <- comm
  bad[1, 1] <- NA
  expect_error(community_to_long(bad, samples, "plot", "hab"), "missing")
  bad[1, 1] <- -1
  expect_error(community_to_long(bad, samples, "plot", "hab"), "non-negative")
  expect_error(community_to_long(comm, samples, "plot", "nope"), "not found")
  expect_error(community_to_long(data.frame(x = "a", sp = 1), samples[1, ],
                                 "plot", "hab"), "numeric")
  s2 <- samples
  s2$lon[2] <- 5
  expect_error(community_to_long(comm, s2, "plot", "hab", keep = "lon"),
               "more than one value")
  expect_error(community_to_long(comm, samples, "plot", "hab", keep = "plot"),
               "must not include")
})
