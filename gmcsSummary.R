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
    "LandR.CS (>=2.0.0)"
  ),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("plotYears", "numeric", NA, NA, NA,
                    desc = "Years to include in the summary plots. NA uses all available years.")
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "gmcsPredictions", objectClass = "list",
                 desc = paste("Raw cohort-level GMCS predictions produced by Biomass_core.",
                              "Note that these predictions are capped by default and thus",
                              "likely obfuscate outliers and errors in the data."))
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "gmcsSummaryByYear", objectClass = "data.table",
                  desc = paste("Year-level summaries of GMCS predictions.")),
    createsOutput(objectName = "gmcsSummaryByYearSpecies", objectClass = "data.table",
                  desc = paste("Species-by-year summaries of GMCS predictions.")),
    
    createsOutput(objectName = "GMCS summary figures", objectClass = "png",
                  desc = paste("Time-series and species-faceted summary plots of mean and",
                               "median GMCS mortality and growth predictions."))
  )
))


doEvent.gmcsSummary <- function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
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
      
      gmcs_dt[, year := time(sim)]
      
      #Summary by year (all species together)
      year_summary <- gmcs_dt[, .(
        mean_mortPred = mean(mortPred, na.rm = TRUE),
        median_mortPred = stats::median(mortPred, na.rm = TRUE),
        mean_growthPred = mean(growthPred, na.rm = TRUE),
        median_growthPred = stats::median(growthPred, na.rm = TRUE)),
        by = year]
      
      if (is.null(sim$gmcsSummaryByYear)) {
        sim$gmcsSummaryByYear <- year_summary
      } else {
        sim$gmcsSummaryByYear <- data.table::rbindlist(list(sim$gmcsSummaryByYear, year_summary), fill = TRUE)
      }
      
      #Summary by species (per year)
      species_summary <- gmcs_dt[, .(
        mean_mortPred = mean(mortPred, na.rm = TRUE),
        median_mortPred = stats::median(mortPred, na.rm = TRUE),
        mean_growthPred = mean(growthPred, na.rm = TRUE),
        median_growthPred = stats::median(growthPred, na.rm = TRUE)),
        by = .(speciesCode, year)]
      
      if (is.null(sim$gmcsSummaryByYearSpecies)) {
        sim$gmcsSummaryByYearSpecies <- species_summary
      } else {
        sim$gmcsSummaryByYearSpecies <- data.table::rbindlist(list(sim$gmcsSummaryByYearSpecies, species_summary), fill = TRUE)
      }
      
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
      plotDataMean[, Metric := fifelse(Metric == "mean_mortPred", "Mortality", "Growth")]
      
      p_time_mean <- ggplot(plotDataMean, aes(x = year, y = Value, colour = Metric, group = Metric)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        facet_wrap(~Metric, scales = "free_y") +
        labs(title = "Mean Mortality and Growth Predictions",
             x = "Year",
             y = "mean predicted value") +
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
      plotDataMedian[, Metric := fifelse(Metric == "median_mortPred", "Mortality", "Growth")]
      
      p_time_median <- ggplot(plotDataMedian, aes(x = year, y = Value, colour = Metric, group = Metric)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        facet_wrap(~Metric, scales = "free_y") +
        labs(title = "Median Mortality and Growth Predictions",
             x = "Year",
             y = "median predicted value") +
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
             y = "mean mortality prediction") +
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
             y = "mean growth prediction") +
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
             y = "median mortality prediction") +
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
             y = "median growth prediction") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position="none")
      
      ggsave(file.path(outDir, "species_median_Growth.png"),
             p_growth_median, width = 12, height = 8, dpi = 150)
    },
    
    warning(paste("Undefined event type:",
                  current(sim)[1, "eventType"],
                  "in module gmcsSummary"))
  )
  invisible(sim)
}