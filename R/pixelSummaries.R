## Per-pixel summaries of GMCS predictions.
##
## Each year the cohort-level predictions in `sim$gmcsPredictions` are collapsed to one value per
## pixelGroup, matched to pixels through `pixelGroupMap`, and added to running sums (pixelGroups are
## re-made as the simulation proceeds, so they cannot be averaged across years; pixels can).
## Units: growth is a % of the growth predicted under the climate normal (100 = no change);
## mortality is the change in mortality (g/m2/yr) relative to the climate normal.

#' Biomass-weighted mean
#'
#' Falls back to the unweighted mean when the weights sum to zero.
#'
#' @param x numeric values.
#' @param w non-negative weights.
#' @return a single number.
#' @keywords internal
.wmean <- function(x, w) {
  sw <- sum(w)
  if (is.na(sw) || sw <= 0) mean(x) else sum(x * w) / sw
}

#' Age class definitions
#'
#' @param ageClassBreaks ascending ages (years) that separate age classes. Classes include their lower
#'   bound, so `c(50, 120)` gives `<50`, `50-120` (50 up to, not including, 120) and `120+`.
#' @return a list with `breaks` (for `cut`), `labels` (for plot titles) and `names` (layer name suffixes).
#' @keywords internal
.ageClassDefs <- function(ageClassBreaks) {
  b <- sort(unique(ageClassBreaks))
  stopifnot(length(b) >= 1, all(is.finite(b)), all(b > 0))
  n <- length(b)
  mid <- if (n > 1) list(paste0(b[-n], "-", b[-1]), paste0(b[-n], "to", b[-1])) else list(NULL, NULL)
  list(breaks = c(-Inf, b, Inf),
       labels = c(paste0("<", b[1]), mid[[1]], paste0(b[n], "+")),
       names = paste0("age_", c(paste0("lt", b[1]), mid[[2]], paste0(b[n], "plus"))))
}

#' Names of the layers in the per-pixel rasters
#'
#' @param species character vector of species codes.
#' @inheritParams .ageClassDefs
#' @return a named list of character vectors: `overall`, `speciesGrowth`, `speciesMortality`,
#'   `ageGrowth` and `ageMortality`.
#' @keywords internal
.pixelLayerNames <- function(species, ageClassBreaks) {
  ac <- .ageClassDefs(ageClassBreaks)
  list(overall = c("growth", "mortality"),
       speciesGrowth = paste0("growth_", species),
       speciesMortality = paste0("mortality_", species),
       ageGrowth = paste0("growth_", ac$names),
       ageMortality = paste0("mortality_", ac$names))
}

#' Collapse cohort predictions to one row per pixelGroup
#'
#' Growth is the biomass-weighted mean of `growthPred`. LandR.CS gives every cohort of a species in a
#' pixelGroup the same `growthPred`, so the weighting only matters when combining species or age
#' classes. Mortality is the sum of `mortPred` (g/m2), i.e. the stand-level (or species-level, or
#' age-class-level) total. Only cohorts that were modelled should be supplied.
#'
#' @param dt cohort-level predictions with `pixelGroup`, `speciesCode`, `age`, `B`, `growthPred`, `mortPred`.
#' @inheritParams .ageClassDefs
#' @return a `data.table` with `pixelGroup` and one column per layer (see `.pixelLayerNames()`);
#'   columns are missing for species or age classes that are absent.
#' @keywords internal
.pixelGroupStats <- function(dt, ageClassBreaks = c(50, 120)) {
  ac <- .ageClassDefs(ageClassBreaks)
  dt <- data.table::as.data.table(dt)[, .(pixelGroup, speciesCode, age, B, growthPred, mortPred)]
  if (!nrow(dt)) return(data.table::data.table(pixelGroup = integer()))
  dt[, ageClass := cut(age, breaks = ac$breaks, labels = ac$names, right = FALSE)]
  summ <- function(groupBy) dt[, .(growth = .wmean(growthPred, B), mortality = sum(mortPred)), by = groupBy]

  overall <- summ("pixelGroup")
  sp <- data.table::dcast(summ(c("pixelGroup", "speciesCode")), pixelGroup ~ speciesCode,
                          value.var = c("growth", "mortality"), sep = "_")
  age <- data.table::dcast(summ(c("pixelGroup", "ageClass")), pixelGroup ~ ageClass,
                           value.var = c("growth", "mortality"), sep = "_")
  Reduce(function(a, b) merge(a, b, by = "pixelGroup", all = TRUE), list(overall, sp, age))
}

#' Set up the running per-pixel sums
#'
#' @param pixelGroupMap the `SpatRaster` of pixelGroups. Pixels with a pixelGroup now are tracked.
#' @param species character vector of all species codes (layers are fixed for the whole simulation).
#' @inheritParams .ageClassDefs
#' @return a list holding the tracked `cells`, and per-cell, per-layer running sums (`sum`) and the
#'   number of years with a value (`n`).
#' @keywords internal
.initPixelAccum <- function(pixelGroupMap, species, ageClassBreaks) {
  cells <- which(!is.na(terra::values(pixelGroupMap, mat = FALSE)))
  layers <- unlist(.pixelLayerNames(species, ageClassBreaks), use.names = FALSE)
  dn <- list(NULL, layers)
  list(cells = cells, layers = layers, species = species, ageClassBreaks = ageClassBreaks,
       sum = matrix(0, length(cells), length(layers), dimnames = dn),
       n = matrix(0L, length(cells), length(layers), dimnames = dn),
       years = numeric())
}

#' Add one year of pixelGroup statistics to the running per-pixel sums
#'
#' @param acc the list made by `.initPixelAccum()`.
#' @param pixelGroupMap the current `SpatRaster` of pixelGroups.
#' @param stats output of `.pixelGroupStats()` for the current year.
#' @param year the simulation year.
#' @return `acc`, updated.
#' @keywords internal
.updatePixelAccum <- function(acc, pixelGroupMap, stats, year) {
  pg <- terra::values(pixelGroupMap, mat = FALSE)
  untracked <- sum(!is.na(pg)) - sum(!is.na(pg[acc$cells]))
  if (untracked > 0)
    warning("gmcsSummary: ", untracked, " pixel(s) have a pixelGroup in ", year,
            " but did not in the first summarized year, and are not included in the per-pixel maps",
            call. = FALSE)
  stats <- data.table::copy(stats)
  for (nm in setdiff(acc$layers, names(stats))) data.table::set(stats, j = nm, value = NA_real_)
  vals <- as.matrix(stats[match(pg[acc$cells], stats$pixelGroup), acc$layers, with = FALSE])
  ok <- !is.na(vals)
  vals[!ok] <- 0
  acc$sum <- acc$sum + vals
  acc$n <- acc$n + ok
  acc$years <- c(acc$years, year)
  acc
}

#' Build the per-pixel rasters from the running sums
#'
#' @inheritParams .updatePixelAccum
#' @return a list of `SpatRaster`s: `overall` (growth, mortality), `speciesGrowth`, `speciesMortality`,
#'   `ageGrowth`, `ageMortality` (each the mean over the years with a value; `NA` where there was none)
#'   and `yearsN` (the number of years with a value, for every layer).
#' @keywords internal
.finalizePixelMaps <- function(acc, pixelGroupMap) {
  mkRaster <- function(layers, asMean) {
    out <- lapply(layers, function(l) {
      v <- if (asMean) {
        x <- acc$sum[, l] / acc$n[, l]
        x[acc$n[, l] == 0L] <- NA_real_
        x
      } else as.numeric(acc$n[, l])
      full <- rep(NA_real_, terra::ncell(pixelGroupMap))
      full[acc$cells] <- v
      r <- terra::rast(pixelGroupMap)
      terra::values(r) <- full
      names(r) <- l
      r
    })
    terra::rast(out)
  }
  nm <- .pixelLayerNames(acc$species, acc$ageClassBreaks)
  c(lapply(nm, mkRaster, asMean = TRUE), list(yearsN = mkRaster(acc$layers, asMean = FALSE)))
}

#' Draw maps into the current plot layout
#'
#' The colour scale is centred (growth on 100, mortality change on 0) and shared by all layers, so
#' that layers of one measure are comparable.
#'
#' @param r a `SpatRaster` whose layers share units.
#' @param center value at the middle of the colour scale.
#' @param titles one title per layer, which should include the units (they are not repeated in a legend title,
#'   as a legend title collides with long plot titles).
#' @param rev reverse the palette (so that red means *more* of the measure).
#' @keywords internal
.drawMaps <- function(r, center, titles, rev = FALSE) {
  m <- suppressWarnings(max(abs(unlist(terra::global(r, "range", na.rm = TRUE)) - center), na.rm = TRUE))
  if (!is.finite(m) || m == 0) m <- 1
  for (i in seq_len(terra::nlyr(r)))
    terra::plot(r[[i]], range = center + c(-m, m), col = grDevices::hcl.colors(100, "RdBu", rev = rev),
                main = titles[i], cex.main = 0.9)
}

#' Save the per-pixel maps as figures
#'
#' Growth and mortality have different units, so every panel title states its units.
#' Mortality maps are stand-level totals (g/m2/yr summed over cohorts), not per-cohort means.
#' Layers with no data (e.g., species that never occur in the study area) are left out.
#'
#' @param maps output of `.finalizePixelMaps()`.
#' @inheritParams .updatePixelAccum
#' @param outDir directory for the figures.
#' @return the file paths, invisibly.
#' @keywords internal
.plotPixelMaps <- function(maps, acc, outDir) {
  growthUnits <- "% of climate-normal growth"
  mortUnits <- "g/m2/yr vs climate normal"
  ac <- .ageClassDefs(acc$ageClassBreaks)
  hasData <- function(r) as.numeric(terra::global(r, "notNA")[[1]]) > 0

  panels <- function(file, r, labels, center, rev, ncols = 3) {
    keep <- hasData(r)
    if (!any(keep)) return(NULL)
    r <- r[[which(keep)]]
    nc <- min(ncols, terra::nlyr(r))
    nr <- ceiling(terra::nlyr(r) / nc)
    grDevices::png(file, width = 4.5 * nc, height = 5 * nr, units = "in", res = 120)
    on.exit(grDevices::dev.off())
    graphics::par(mfrow = c(nr, nc))
    .drawMaps(r, center, labels[keep], rev)
    file
  }
  gT <- function(what) paste0(what, "\nGrowth (", growthUnits, ")")
  mT <- function(what) paste0(what, "\nMortality change (", mortUnits, ")")

  files <- character()
  if (any(hasData(maps$overall))) {
    f <- file.path(outDir, "map_overall.png")
    grDevices::png(f, width = 10, height = 5, units = "in", res = 120)
    graphics::par(mfrow = c(1, 2))
    .drawMaps(maps$overall[["growth"]], 100, gT("Mean, all cohorts"))
    .drawMaps(maps$overall[["mortality"]], 0, mT("Mean stand-level total"), rev = TRUE)
    grDevices::dev.off()
    files <- c(files, f)
  }
  files <- c(files,
    panels(file.path(outDir, "map_species_growth.png"), maps$speciesGrowth,
           gT(acc$species), 100, FALSE),
    panels(file.path(outDir, "map_species_mortality.png"), maps$speciesMortality,
           mT(paste(acc$species, "(species total)")), 0, TRUE),
    panels(file.path(outDir, "map_ageClass_growth.png"), maps$ageGrowth,
           gT(paste("Age", ac$labels)), 100, FALSE),
    panels(file.path(outDir, "map_ageClass_mortality.png"), maps$ageMortality,
           mT(paste("Age", ac$labels, "(class total)")), 0, TRUE))
  invisible(files)
}
