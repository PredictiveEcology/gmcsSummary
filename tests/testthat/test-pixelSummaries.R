## Helpers are internal; load them directly so the tests also work outside an installed package
source(file.path("..", "..", "R", "pixelSummaries.R"))

mkCohorts <- function() {
  data.table::data.table(
    pixelGroup = c(1L, 1L, 1L, 2L),
    speciesCode = factor(c("Pice_gla", "Pice_gla", "Popu_tre", "Pice_gla"), levels = c("Pice_gla", "Popu_tre")),
    age = c(30L, 130L, 60L, 30L),
    B = c(100L, 300L, 200L, 50L),
    growthPred = c(80L, 80L, 120L, 50L),
    mortPred = c(10, 30, 20, 5))
}

test_that("age classes include their lower bound and can be customised", {
  expect_equal(.ageClassDefs(c(50, 120))$labels, c("<50", "50-120", "120+"))
  expect_equal(.ageClassDefs(c(20, 60, 100))$names, c("age_lt20", "age_20to60", "age_60to100", "age_100plus"))
  expect_equal(.ageClassDefs(30)$labels, c("<30", "30+"))
  expect_equal(as.character(cut(c(49, 50, 119, 120), .ageClassDefs(c(50, 120))$breaks,
                                labels = .ageClassDefs(c(50, 120))$names, right = FALSE)),
               c("age_lt50", "age_50to120", "age_50to120", "age_120plus"))
})

test_that("pixelGroup statistics: growth is biomass-weighted, mortality is a total", {
  st <- .pixelGroupStats(mkCohorts(), c(50, 120))
  pg1 <- st[pixelGroup == 1]
  expect_equal(pg1$growth, weighted.mean(c(80, 80, 120), c(100, 300, 200)))
  expect_equal(pg1$mortality, 60)
  # per species: the Pice_gla cohorts share one growth value; their mortality is summed
  expect_equal(pg1$growth_Pice_gla, 80)
  expect_equal(pg1$mortality_Pice_gla, 40)
  expect_equal(pg1$mortality_Popu_tre, 20)
  # per age class
  expect_equal(pg1$mortality_age_lt50, 10)
  expect_equal(pg1$mortality_age_50to120, 20)
  expect_equal(pg1$mortality_age_120plus, 30)
  expect_true(is.na(st[pixelGroup == 2, growth_Popu_tre]))
})

test_that("per-pixel means are taken over years with data, and absent species stay NA", {
  pgm <- terra::rast(nrows = 2, ncols = 2, vals = c(1, 1, 2, NA))
  acc <- .initPixelAccum(pgm, c("Pice_gla", "Popu_tre"), c(50, 120))
  st1 <- .pixelGroupStats(mkCohorts(), c(50, 120))
  c2 <- mkCohorts(); c2[, growthPred := growthPred + 10L]
  st2 <- .pixelGroupStats(c2[pixelGroup == 1], c(50, 120)) # pixelGroup 2 has no data in year 2
  acc <- .updatePixelAccum(acc, pgm, st1, 2020)
  acc <- .updatePixelAccum(acc, pgm, st2, 2021)
  maps <- .finalizePixelMaps(acc, pgm)

  expect_equal(names(maps$overall), c("growth", "mortality"))
  expect_equal(names(maps$ageGrowth), c("growth_age_lt50", "growth_age_50to120", "growth_age_120plus"))
  v <- terra::values(maps$overall[["growth"]], mat = FALSE)
  y1 <- weighted.mean(c(80, 80, 120), c(100, 300, 200))
  y2 <- weighted.mean(c(90, 90, 130), c(100, 300, 200))
  expect_equal(v[1:2], rep((y1 + y2) / 2, 2))
  expect_equal(v[3], 50)               # pixelGroup 2: only year 1
  expect_true(is.na(v[4]))             # no pixelGroup
  yn <- terra::values(maps$yearsN[["growth"]], mat = FALSE)
  expect_equal(yn[1:4], c(2, 2, 1, NA))
  # Popu_tre is absent from pixelGroup 2 in both years
  expect_true(is.na(terra::values(maps$speciesGrowth[["growth_Popu_tre"]], mat = FALSE)[3]))
})

test_that("pixels that gain a pixelGroup after the first year are reported", {
  pgm1 <- terra::rast(nrows = 1, ncols = 2, vals = c(1, NA))
  pgm2 <- terra::rast(nrows = 1, ncols = 2, vals = c(1, 2))
  acc <- .initPixelAccum(pgm1, c("Pice_gla", "Popu_tre"), 50)
  expect_warning(.updatePixelAccum(acc, pgm2, .pixelGroupStats(mkCohorts(), 50), 2021), "did not in the first")
})
