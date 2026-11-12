test_that("check_coverage() summarises habitats per species", {
  cov <- do.call(check_coverage, c(list(moorea_long(c(ref_species, "Acanthurus guttatus"))),
                                   long_cols))
  expect_equal(nrow(cov), 12)
  expect_equal(unique(cov$n_sites[cov$species == ref_species[1]]), 13)
  expect_equal(cov$status[cov$species %in% ref_species], rep("used", 9))
  d <- moorea_long(ref_species[1])
  expect_equal(cov$total_count[cov$species == ref_species[1]],
               as.vector(tapply(d$abundance, d$habitat, sum)[moorea_habs]))
  all_cov <- do.call(check_coverage, c(list(moorea_fish), long_cols))
  expect_setequal(unique(all_cov$status), c("used", "absent", "undersampled"))
})

# Six sites in two tight clusters 30 km apart, plus a seventh site
geo_data <- function() {
  xy <- data.frame(site = paste0("s", 1:7),
                   lon = c(0, 0.005, 0.01, 0.3, 0.305, 0.31, 0.155),
                   lat = 0)
  d <- expand.grid(site = xy$site, habitat = c("a", "b"), stringsAsFactors = FALSE)
  d$abundance <- seq_len(nrow(d))
  d$effort <- 2
  merge(d, xy)
}

test_that("group_sites() merges sites within the radius and sums counts and effort", {
  d <- geo_data()
  g <- group_sites(d, site = "site", lon = "lon", lat = "lat", radius_km = 2,
                   habitat = "habitat", abundance = "abundance", effort = "effort")
  expect_equal(sort(unique(g$site)), c("s1+s2+s3", "s4+s5+s6", "s7"))
  ga <- g[g$site == "s1+s2+s3" & g$habitat == "a", ]
  expect_equal(ga$abundance, sum(d$abundance[d$site %in% c("s1", "s2", "s3") &
                                               d$habitat == "a"]))
  expect_equal(ga$effort, 6)
  expect_equal(ga$n_sites, 3L)
  expect_equal(sum(g$abundance), sum(d$abundance))
  # an effort column that is not used is pointed out
  expect_warning(group_sites(d, site = "site", lon = "lon", lat = "lat",
                             radius_km = 2, habitat = "habitat",
                             abundance = "abundance"), "effort column")
  # radius 0 keeps every site; no effort column -> number of sites
  d$effort <- NULL
  g0 <- group_sites(d, site = "site", lon = "lon", lat = "lat", radius_km = 0,
                    habitat = "habitat", abundance = "abundance")
  expect_equal(length(unique(g0$site)), 7)
  expect_true(all(g0$effort == 1))
  # complete linkage never joins sites further apart than the radius
  gc <- group_sites(d, site = "site", lon = "lon", lat = "lat", radius_km = 18,
                    habitat = "habitat", abundance = "abundance")
  gs <- group_sites(d, site = "site", lon = "lon", lat = "lat", radius_km = 18,
                    habitat = "habitat", abundance = "abundance", linkage = "single")
  expect_gt(length(unique(gc$site)), length(unique(gs$site)))
  expect_error(group_sites(d, site = "site", lon = "lon", lat = "lat",
                           radius_km = -1, habitat = "habitat",
                           abundance = "abundance"), "radius_km")
  d$lon[1] <- 5
  expect_error(group_sites(d, site = "site", lon = "lon", lat = "lat",
                           radius_km = 1, habitat = "habitat",
                           abundance = "abundance"), "single location")
})

test_that("group_sites() takes a radius per species", {
  d <- geo_data()
  d$effort <- NULL
  d <- rbind(transform(d, sp = "x"), transform(d, sp = "y"))
  g <- group_sites(d, site = "site", lon = "lon", lat = "lat",
                   radius_km = c(y = 0, x = 2), species = "sp",
                   habitat = "habitat", abundance = "abundance")
  expect_equal(length(unique(g$site[g$sp == "x"])), 3)
  expect_equal(length(unique(g$site[g$sp == "y"])), 7)
  expect_error(group_sites(d, site = "site", lon = "lon", lat = "lat",
                           radius_km = c(x = 2), species = "sp",
                           habitat = "habitat", abundance = "abundance"),
               "every species")
})

test_that("radius_sensitivity() returns one row per radius", {
  set.seed(1)
  d <- simulate_isodars(1.5, 1.5, 5, 5, n = 20, p_missing = 0.3, format = "long")
  xy <- data.frame(site = unique(d$site), lon = stats::runif(20, 0, 0.1),
                   lat = stats::runif(20, 0, 0.1))
  d <- merge(d, xy)
  rs <- radius_sensitivity(d, site = "site", lon = "lon", lat = "lat",
                           habitat = "habitat", abundance = "abundance",
                           radii_km = c(0, 3, 6), diagnostics = list(n_null = 0))
  expect_equal(rs$radius_km, c(0, 3, 6))
  expect_equal(rs$n_units[1], 20L)
  expect_true(all(diff(rs$n_units) <= 0))
  expect_true(all(rs$complete_share >= 0 & rs$complete_share <= 1))
})

test_that("simulate_isodars() generalises to any number of habitats", {
  set.seed(4)
  old <- simulate_isodars(2, 2, 0, 0, noise = 0.5, n = 60)
  set.seed(4)
  new <- simulate_isodars(slopes = c(2, 2), intercepts = c(0, 0), noise = 0.5, n = 60)
  expect_equal(old, new)
  set.seed(1)
  s <- simulate_isodars(slopes = c(1, 2, 0.5), intercepts = c(3, 0, 1), n = 50,
                        p_missing = 0.5)
  expect_named(s, c("hab1", "hab2", "hab3", "hab4", "total_abundance"))
  expect_true(all(rowSums(!is.na(s[1:4])) >= 2))
  expect_gt(mean(is.na(as.matrix(s[1:4]))), 0.2)
  long <- simulate_isodars(slopes = c(1, 2, 0.5), intercepts = c(3, 0, 1), n = 5,
                           format = "long")
  expect_named(long, c("site", "habitat", "abundance"))
  expect_equal(nrow(long), 20)
  expect_error(simulate_isodars(slopes = 1:2, intercepts = 1), "same length")
  expect_error(simulate_isodars(p_missing = 1), "p_missing")
})
