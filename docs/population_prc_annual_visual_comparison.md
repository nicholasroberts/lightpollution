# Visual investigation: does annual population change coincide spatially with annual spider/moth PRC change?

**Research priority — 21 October 2026.** This is a map- and observation-first descriptive comparison, NOT an explanatory-power model. The previous two-band geography analysis described *where* large PRC change lies; we already know that nonlinear PRC(R) sensitivity increases the importance of modestly lit pixels. Here we ask whether **the NUTS3 regions whose population density changes from year to year are the same regions whose PRC changes**, and whether the sign and magnitude of the changes correlate.

## Existing sources (no user uploads required)

1. `data/processed/nuts3_population_density.parquet`: genuine annual Eurostat density with `nuts_release == 2016`, years **2014–2020**, with `matches_nuts_release` checked. Do NOT merge unharmonised 2021–24 Eurostat NUTS releases into 2016 geography.
2. `data/processed/historical_prc_2014_2024_spider_snow_corrected/prc_nuts3_annual_cumulative.csv` and corresponding `moth_snow_corrected`: **locked** corrected biological records, with `coverage_stable` verified.
3. `data/processed/fixed_nuts2016_viirs_2014_2024.parquet`: annual VIIRS corrected radiance, contextualised with `radiance_2014` and pre-existing 2014–20 population-density categories; **no unfitted LED correction claimed**.
4. NUTS2016 polygon geography via project configuration. Baseline radiance categories remain 0.2–1 and 1–5 nW as before.

Use the **same region cohort with complete seven-year population and PRC records** for all six annual comparisons and for the 2014→2020 endpoint map, ensuring maps are comparable. Regions missing any annual population datum remain grey; never treat missing data as no change.

## Figures in `figures/prc_population_change_2014_2020/`

- **Six three-panel paired maps** `annual_population_spider_moth_2014_2015.pdf` through `annual_population_spider_moth_2019_2020.pdf`. Panels: **population density annual percentage change**, **spider absolute ΔPRC**, **moth absolute ΔPRC**, all on the identical NUTS2016 polygons, colour-centred on zero and **constant per-variable limits across all years**.
- **One full-period three-panel map** `cumulative_population_spider_moth_2014_2020.pdf`, population density percentage change and spider/moth absolute PRC change. All values 2020 minus 2014, not the 2014–2024 PRC endpoints used in earlier figures.
- **Two focused versions**, `cumulative_radiance_02_to_1_population_spider_moth_2014_2020.pdf` and `cumulative_radiance_1_to_5_population_spider_moth_2014_2020.pdf`, masking all NUTS3 regions outside each important **2014 mean VIIRS radiance class**.
- **Six two-species annual scatter pages** `annual_population_vs_spider_moth_prc_scatter_2015_2020.pdf`: Δpopulation (%) x-axis, each taxon's absolute annual ΔPRC y-axis, points coloured by **2014 radiance class**, with Pearson r and Spearman ρ printed.
- **One cumulative two-species scatter sheet** `cumulative_population_vs_spider_moth_prc_scatter_2014_2020.pdf`. This directly answers whether regions with the largest population change also have the greatest cumulative modelled biological contrast changes.

## Data outputs

- `matched_population_and_prc_annual_2014_2020.csv` — seven measured yearly values per common NUTS3 ID including year-on-year and 2014-relative quantities.
- `annual_changes_2015_2020.csv`, `cumulative_change_2014_2020.csv` — direct analysis inputs.
- `correlations_annual_and_cumulative.csv` — raw Pearson and Spearman plus p-values for **each individual year transition** and total 2014–2020; observed population density absolute and percentage changes against spider and moth absolute and percentage PRC change.
- `correlations_cumulative_by_baseline_radiance.csv` — the same endpoint association stratified by 2014 radiance category, so global association is not mistaken for a correlation induced by their shared dependence on pre-existing brightness/density.

The previous nested-model study `src/fit_annual_population_prc_change.py` remains **separate** and optional. Map and visualise the measurements **first**. A correlation is descriptive and may contain confounding from policy, geography, industrial lighting, spectral replacements, baseline brightness and country-level population shifts; it is not yet a causal proportion of biological change attributable to demography.

### EcoVis

```bash
git pull --ff-only
python src/visualize_annual_population_prc_changes.py --self-test
python src/visualize_annual_population_prc_changes.py
```

Upload at first: `cumulative_population_spider_moth_2014_2020.pdf`, `annual_population_vs_spider_moth_prc_scatter_2015_2020.pdf`, `correlations_annual_and_cumulative.csv`. Then we can inspect the geographically aligned annual PDFs, and discuss whether year-on-year population changes correspond to the places with biological PRC changes before doing conditional model attribution.
