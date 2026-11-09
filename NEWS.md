# gmcsSummary 0.0.2

- annual summaries now use only modelled cohorts (those at or above `gmcsMinAge`) and weight growth by biomass; the
  year and species summaries also record `climateYear` and the number of cohorts
- per-pixel maps of the mean (over years) growth and mortality, built each year by matching `gmcsPredictions` to
  `pixelGroupMap`: overall, by species, and by age class (parameter `ageClassBreaks`), plus the number of years behind each
  value; saved as GeoTIFFs and figures
- units are stated in every plot: growth is a % of the growth predicted under the climate normal, mortality is the
  change in g/m2/yr relative to the climate normal
- `Plotting` now runs after the final `Summarizing` event, so the last year is included

# gmcsSummary 0.0.1 (08 April 2026)

- initial module version
