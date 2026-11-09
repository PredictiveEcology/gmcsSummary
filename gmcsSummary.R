defineModule(sim, list(
  name = "gmcsSummary",
  description = "Summarize and plot GMCS mortality and growth predictions produced by Biomass_core",
  keywords = c("LandR", "LandR.CS", "GMCS"),
  authors = c(
    person("Jonathan", "Van Elslander", email = "jonathan.vanelslander@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person(c("Ian", "MS"), "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "ctb")
  ),
  childModules = character(0),
  version = list(gmcsSummary = "0.3.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "gmcsSummary.Rmd"),
  reqdPkgs = list(
    "SpaDES.core (>= 3.0.4)",
    "data.table",
    "ggplot2",
    "LandR.CS (>=2.0.0)",
    "terra"
  ),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("plotYears", "numeric", NA, NA, NA,
                    desc = "Years to include in the summary plots. NA uses all available years."),
    defineParameter("ageClassBreaks", "numeric", c(50, 120), 1, NA,
                    desc = paste("Ages (years) that separate the age classes used in the per-pixel age-class maps.",
                                 "Classes include their lower bound, so the default gives <50, 50-120 and 120+."))
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "gmcsPredictions", objectClass = "list",
                 desc = paste("Raw cohort-level GMCS predictions produced by Biomass_core.",
                              "Note that these predictions are capped by default and thus",
                              "likely obfuscate outliers and errors in the data.",
                              "Units: `growthPred` is a % of the growth predicted under the climate normal;",
                              "`mortPred` is the change in mortality (g/m2) per YEAR relative to the climate normal,",
                              "because the models are fitted to annual rates. Run Biomass_core with a 1-year time step.")),
    expectsInput(objectName = "pixelGroupMap", objectClass = "SpatRaster",
                 desc = paste("Raster of pixelGroups from Biomass_core, used to place the cohort-level",
                              "predictions on the landscape."))
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "gmcsSummaryByYear", objectClass = "data.table",
                  desc = paste("Year-level summaries of GMCS predictions.")),
    createsOutput(objectName = "gmcsSummaryByYearSpecies", objectClass = "data.table",
                  desc = paste("Species-by-year summaries of GMCS predictions.")),
    
    createsOutput(objectName = "GMCS summary figures", objectClass = "png",
                  desc = paste("Time-series and species-faceted summary plots of mean and",
                               "median GMCS mortality and growth predictions.")),
    createsOutput(objectName = "gmcsPixelAccum", objectClass = "list",
                  desc = "Running per-pixel sums of the annual predictions (internal)."),
    createsOutput(objectName = "gmcsMapOverall", objectClass = "SpatRaster",
                  desc = paste("Per-pixel mean over years of growth (% of the growth predicted under the climate normal)",
                               "and mortality (change in g/m2/yr relative to the climate normal); layers `growth`, `mortality`.")),
    createsOutput(objectName = "gmcsMapSpeciesGrowth", objectClass = "SpatRaster",
                  desc = "Per-pixel mean annual growth (% of normal), one layer per species."),
    createsOutput(objectName = "gmcsMapSpeciesMortality", objectClass = "SpatRaster",
                  desc = "Per-pixel mean annual mortality change (g/m2/yr vs normal), one layer per species."),
    createsOutput(objectName = "gmcsMapAgeGrowth", objectClass = "SpatRaster",
                  desc = "Per-pixel mean annual growth (% of normal), one layer per age class (see `ageClassBreaks`)."),
    createsOutput(objectName = "gmcsMapAgeMortality", objectClass = "SpatRaster",
                  desc = "Per-pixel mean annual mortality change (g/m2/yr vs normal), one layer per age class."),
    createsOutput(objectName = "gmcsMapYearsN", objectClass = "SpatRaster",
                  desc = "Number of years contributing to each layer of the per-pixel means, same layer names.")
  )
))


doEvent.gmcsSummary <- function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      # mortPred is an annual rate (g/m2/yr); a longer simulation time unit would under-state mortality
      if (!identical(tolower(sub("s$", "", timeunit(sim))), "year"))
        warning("gmcsSummary: the simulation time unit is '", timeunit(sim), "', but GMCS growth and mortality ",
                "predictions are annual rates (mortality in g/m2/yr); the summaries and maps assume 1-year steps.",
                call. = FALSE)
      sim <- scheduleEvent(sim, start(sim), module = "gmcsSummary", eventType = "Summarizing", eventPriority = .last())
      ## later priority than the final Summarizing event, so the last year is included in the plots
      sim <- scheduleEvent(sim, end(sim), module = "gmcsSummary", eventType = "Plotting", eventPriority = .last() + 1)
    },
    
    Summarizing = {
      
      if (is.null(sim$gmcsPredictions) || length(sim$gmcsPredictions) == 0) {
        stop("sim$gmcsPredictions not found or empty.")
      }
      
      gmcs_dt <- data.table::rbindlist(sim$gmcsPredictions,
                                       fill = TRUE, use.names = TRUE, idcol = "predictionGroup")
      
      stopifnot(all(c("pixelGroup", "speciesCode", "age", "B", "growthPred", "mortPred",
                      "pred_currentGrowth") %in% names(gmcs_dt)))
      gmcs_dt[, year := time(sim)]

      # Cohorts younger than `gmcsMinAge` are not modelled (LandR.CS fixes them at growth 100, mortality 0);
      # keeping them would dilute every summary below
      gmcs_dt <- gmcs_dt[!is.na(pred_currentGrowth)]

      # Summary by year (all species together)
      # growth: % of the growth predicted under the climate normal (biomass-weighted mean; median unweighted)
      # mortality: change in mortality relative to the normal, g/m2/yr, per cohort
      year_summary <- gmcs_dt[, .(
        climateYear = climateYear[1],
        n_cohorts = .N,
        mean_mortPred = mean(mortPred, na.rm = TRUE),
        median_mortPred = stats::median(mortPred, na.rm = TRUE),
        mean_growthPred = .wmean(growthPred, B),
        median_growthPred = stats::median(growthPred, na.rm = TRUE)),
        by = year]
      
      if (is.null(sim$gmcsSummaryByYear)) {
        sim$gmcsSummaryByYear <- year_summary
      } else {
        sim$gmcsSummaryByYear <- data.table::rbindlist(list(sim$gmcsSummaryByYear, year_summary), fill = TRUE)
      }
      
      # Summary by species (per year)
      species_summary <- gmcs_dt[, .(
        climateYear = climateYear[1],
        n_cohorts = .N,
        mean_mortPred = mean(mortPred, na.rm = TRUE),
        median_mortPred = stats::median(mortPred, na.rm = TRUE),
        mean_growthPred = .wmean(growthPred, B),
        median_growthPred = stats::median(growthPred, na.rm = TRUE)),
        by = .(speciesCode, year)]
      
      if (is.null(sim$gmcsSummaryByYearSpecies)) {
        sim$gmcsSummaryByYearSpecies <- species_summary
      } else {
        sim$gmcsSummaryByYearSpecies <- data.table::rbindlist(list(sim$gmcsSummaryByYearSpecies, species_summary), fill = TRUE)
      }
      
      # Per-pixel running sums: pixelGroups are re-made as the simulation proceeds, so each year's values
      # are placed on pixels now and averaged over years at the end
      if (is.null(sim$pixelGroupMap)) stop("sim$pixelGroupMap not found.")
      if (is.null(sim$gmcsPixelAccum)) {
        spp <- if (is.factor(gmcs_dt$speciesCode)) levels(gmcs_dt$speciesCode) else sort(unique(gmcs_dt$speciesCode))
        sim$gmcsPixelAccum <- .initPixelAccum(sim$pixelGroupMap, spp, P(sim)$ageClassBreaks)
      }
      sim$gmcsPixelAccum <- .updatePixelAccum(sim$gmcsPixelAccum, sim$pixelGroupMap,
                                              .pixelGroupStats(gmcs_dt, P(sim)$ageClassBreaks), time(sim))

      if (time(sim) < end(sim)) {
        sim <- scheduleEvent(sim, time(sim) + 1, module = "gmcsSummary", eventType = "Summarizing", eventPriority = .last())
      }
    },
    
    Plotting = {
      
      outDir <- figurePath(sim)
      if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
      
      if (is.null(sim$gmcsSummaryByYear) ||
          is.null(sim$gmcsSummaryByYearSpecies)) {
        stop(
          "Required GMCS summary objects not found in simList.\n",
          "Ensure gmcsSummary Summarizing event ran successfully."
        )
      }
      
      yearDT <- data.table::as.data.table(sim$gmcsSummaryByYear)
      spDT   <- data.table::as.data.table(sim$gmcsSummaryByYearSpecies)
      
      if (!is.na(P(sim)$plotYears)[1]) {
        yearDT <- yearDT[year %in% P(sim)$plotYears]
        spDT   <- spDT[year %in% P(sim)$plotYears]
      }
      
      plotDataMean <- data.table::melt(yearDT, id.vars = "year",
                                       measure.vars = c("mean_mortPred", "mean_growthPred"),
                                       variable.name = "Metric", value.name = "Value"
      )
      plotDataMean[, Metric := fifelse(Metric == "mean_mortPred", "Mortality change per cohort (g/m2/yr vs climate normal)",
                                 "Growth (% of climate-normal growth)")]
      
      p_time_mean <- ggplot(plotDataMean, aes(x = year, y = Value, colour = Metric, group = Metric)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        facet_wrap(~Metric, scales = "free_y") +
        labs(title = "Mean Mortality and Growth Predictions",
             x = "Year",
             y = NULL) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position="none")
      
      ggsave(file.path(outDir, "ModelsOverTime_mean.png"),
             p_time_mean, width = 10, height = 6, dpi = 150)
      
      plotDataMedian <- data.table::melt(yearDT, id.vars = "year",
                                         measure.vars = c("median_mortPred", "median_growthPred"),
                                         variable.name = "Metric", value.name = "Value"
      )
      plotDataMedian[, Metric := fifelse(Metric == "median_mortPred", "Mortality change per cohort (g/m2/yr vs climate normal)",
                                   "Growth (% of climate-normal growth)")]
      
      p_time_median <- ggplot(plotDataMedian, aes(x = year, y = Value, colour = Metric, group = Metric)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        facet_wrap(~Metric, scales = "free_y") +
        labs(title = "Median Mortality and Growth Predictions",
             x = "Year",
             y = NULL) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position="none")
      
      ggsave(file.path(outDir, "ModelsOverTime_median.png"),
             p_time_median, width = 10, height = 6, dpi = 150)
      
      ## mean by species
      p_mort_mean <- ggplot(spDT, aes(x = year, y = mean_mortPred, colour = speciesCode, group = speciesCode)) +
        geom_line() +
        geom_point() +
        facet_wrap(~speciesCode, scales = "free_y") +
        labs(title = "Mean GMCS Mortality Predictions by Species",
             x = "Year",
             y = "mean per-cohort mortality change (g/m2/yr vs climate normal)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position="none")
      
      ggsave(file.path(outDir, "species_mean_Mortality.png"),
             p_mort_mean, width = 12, height = 8, dpi = 150)
      
      p_growth_mean <- ggplot(spDT, aes(x = year, y = mean_growthPred, colour = speciesCode, group = speciesCode)) +
        geom_line() +
        geom_point() +
        facet_wrap(~speciesCode, scales = "free_y") +
        labs(title = "Mean GMCS Growth Predictions by Species",
             x = "Year",
             y = "mean growth (% of climate-normal growth)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position="none")
      
      ggsave(file.path(outDir, "species_mean_Growth.png"),
             p_growth_mean, width = 12, height = 8, dpi = 150)
      
      ## median by species
      p_mort_median <- ggplot(spDT, aes(x = year, y = median_mortPred, colour = speciesCode, group = speciesCode)) +
        geom_line() +
        geom_point() +
        facet_wrap(~speciesCode, scales = "free_y") +
        labs(title = "Median GMCS Mortality Predictions by Species",
             x = "Year",
             y = "median per-cohort mortality change (g/m2/yr vs climate normal)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position="none")
      
      ggsave(file.path(outDir, "species_median_Mortality.png"),
             p_mort_median, width = 12, height = 8, dpi = 150)
      
      p_growth_median <- ggplot(spDT, aes(x = year, y = median_growthPred, colour = speciesCode, group = speciesCode)) +
        geom_line() +
        geom_point() +
        facet_wrap(~speciesCode, scales = "free_y") +
        labs(title = "Median GMCS Growth Predictions by Species",
             x = "Year",
             y = "median growth (% of climate-normal growth)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position="none")
      
      ggsave(file.path(outDir, "species_median_Growth.png"),
             p_growth_median, width = 12, height = 8, dpi = 150)

      ## per-pixel maps: mean over years of the annual predictions
      if (is.null(sim$gmcsPixelAccum)) stop("No per-pixel summaries found; the Summarizing event must run first.")
      maps <- .finalizePixelMaps(sim$gmcsPixelAccum, sim$pixelGroupMap)
      sim$gmcsMapOverall <- maps$overall
      sim$gmcsMapSpeciesGrowth <- maps$speciesGrowth
      sim$gmcsMapSpeciesMortality <- maps$speciesMortality
      sim$gmcsMapAgeGrowth <- maps$ageGrowth
      sim$gmcsMapAgeMortality <- maps$ageMortality
      sim$gmcsMapYearsN <- maps$yearsN
      for (nm in c("Overall", "SpeciesGrowth", "SpeciesMortality", "AgeGrowth", "AgeMortality", "YearsN")) {
        terra::writeRaster(sim[[paste0("gmcsMap", nm)]], file.path(outputPath(sim), paste0("gmcsMap", nm, ".tif")),
                           overwrite = TRUE)
      }
      .plotPixelMaps(maps, sim$gmcsPixelAccum, outDir)
    },
    
    warning(paste("Undefined event type:",
                  current(sim)[1, "eventType"],
                  "in module gmcsSummary"))
  )
  invisible(sim)
}