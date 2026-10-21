# Matched Jan1 population, radiance and locked PRC: country contrasts and NL333 locality

**21 October 2026. Analysis prepared; EcoVis results pending.** This extends—not repeats—the 2014–2020 reconstructed-population *country–PRC* sensitivity already reported in `docs/results_population_jan1_fixed_nuts2016_2014_2020.md`. The next question is whether country-specific changes in Jan1 population density coincide with changes in calibrated mean VIIRS radiance and pixel-derived spider/moth polarization receptor contrast, and whether regional and national results differ.

## Locked sources and scope

`src/analyze_jan1_country_radiance_prc_nl333.py` reuses `read_matched_endpoints` and `add_radiance_endpoints` in the existing four-panel plotting script, preserving the **1,273** complete NUTS2016 regions, 2014 and 2020 endpoints, January-1 population divided by frozen region land area and corrected annual VIIRS regional radiance. Both taxa use the **unchanged** spider/moth pixel-derived PRC values stored in the reconstructed endpoint table. The command does not regenerate, recalibrate or overwrite original processed sources.

**Study hypothesis, not assumed outcome:** demographic change is an incomplete predictor of regional/national change in anthropogenic nighttime radiance and inferred biological polarization cue. Some countries may nevertheless exhibit substantial associations, notably Türkiye and Germany. Country associations cannot prove absence of a demographic effect; they can reject oversimplified universal relationships and identify exceptions.

## Comparisons

1. **Within each country**: separate Spearman rho and Pearson r for continuous signed ΔJan1 population density (people/km²) against (a) Δcorrected NUTS3 mean radiance, (b) Δspider PRC, (c) Δmoth PRC, plus country sample sizes, medians and density-growth counts. This includes all eligible countries, not only selected illustrative nations. The earlier Jan1–PRC rho values are benchmarks to reproduce (TR spider −.537, DE −.381, FR +.223, IE +.262, IT −.163, for 2014→2020). The direct population–radiance comparison is the new material.
2. **Between countries**: a single country-level statistic computed as the fixed-land-area-weighted mean of the matched NUTS3 absolute density/radiance/PRC changes; calculate association across countries giving **one point/equal regression weight per country**, rather than pooling all NUTS3 observations (which would overweight Germany). Report inclusion restrictions with at least 1, 5, or 10 matched NUTS3 regions. This is a descriptive ecological country-level relationship, not a country fixed-effects result. Country mean PRC is the regional **pixel-transformed PRC mean**, never PRC evaluated at national mean radiance.
3. **Netherlands**: identify NL333 Delft en Westland in the Dutch matched cohort, report baseline brightness and its ranking in |2014–2020 Δradiance|, |Δspider PRC| and |Δmoth PRC|, recompute Dutch correlations with/without NL333, and compare true NUTS2016 polygons sharing its border or lying within 25 or 50 km measured on EPSG:3035. A four-panel local map highlights NL333 (black outline) and touching Dutch NUTS3 neighbours (orange). This is NOT a municipal or greenhouse-level analysis. Do not import NL333's 2014→2024 radiance anomaly (−179.5 nW) as if it applied in 2014→2020; its 2020 value must be measured independently.

## Reproducible EcoVis execution

From project root, with existing environment activated:

```bash
git switch historical-polarization-visibility
git pull --ff-only
python src/analyze_jan1_country_radiance_prc_nl333.py --self-test
python src/analyze_jan1_country_radiance_prc_nl333.py
```

Expected output directory `data/processed/jan1_population_radiance_prc_country_nl333/`:

- `matched_region_change_2014_2020.csv`
- `within_country_correlations_2014_2020.csv`
- `country_area_weighted_changes_2014_2020.csv`
- `between_country_correlations_2014_2020.csv`
- `nl333_influence_2014_2020.csv`
- `nl333_neighbourhood_nuts3_2014_2020.csv`
- `nl333_neighbourhood_correlations_2014_2020.csv`
- `analysis_manifest.json`

Figures under `figures/jan1_population_radiance_prc_country_nl333/`:

- `selected_country_population_vs_radiance_spider_moth.pdf`: 6 × 3 country-specific population–outcome scatter comparisons; NL333 marked independently
- `between_country_population_vs_radiance_spider_moth.pdf`: country-level scatter for all three outcomes
- `nl333_four_panel_local_80km.pdf` and `.png`: geographical map around Delft en Westland

**Required review:** inspect PDFs and tables before conclusions. Compare existing Jan1 country–PRC rhos to verify the locked pipeline was reproduced. Nominal independent-region significance p-values are deliberately not reported: neighbouring NUTS3 regions are spatially autocorrelated, and NL333 local neighbourhoods can have too few units. A lack of correlation does not prove absence of population effects. In particular the radiance and two PRC responses are not independent because PRC is calculated from the same VIIRS record. Area-weighting by administrative land area is an explicit approximation, not pixel-overlap adjusted inference or human/arthropod exposure weighting.

**Next decision once outputs arrive:** identify empirically defensible national contrasts, determine whether NL333 is influential and/or spatially unusual, and whether a local independent-source/horticultural-lighting investigation is worthwhile. Only then discuss stronger multilevel spatial models or plausible mechanisms.