# Current research state / restart point

**Updated:** 2026-10-13  
**Branch:** `future-population-models`  
**Current analytical commit before documentation updates:** `2d79e1da`  
**Machine:** EcoVis  
**Repository root:** `/run/media/bznwr/data/lightpollution`

## Immediate task

Rerun the literature-led MOD10CM snow test after the raster-level fix that converts every MOD10CM categorical code >100 to true nodata before NUTS3 extraction.

The previous snow outputs are **invalid** because categorical values leaked into regional snow means, producing physically impossible values >100%.

## Commands

```bash
git pull
python src/run_pipeline.py future-radiance-snow-2021-test
```

Expected startup line:

```text
Snow/VIIRS literature-method integration self-test: PASS
```

The cleaning-version change should force regeneration of the cleaned MOD10CM GeoTIFFs and snow NUTS3 caches. VIIRS literature metrics should remain cached.

## Do not interpret until

- pipeline completes cleanly;
- all finite `snow_mean_pct` values satisfy 0 <= value <= 100;
- February primary, January supportive and March negative-control results are available.

## Return these files

```text
data/processed/fixed_nuts2016_snow_viirs_test_summary.txt
data/processed/fixed_nuts2016_snow_2021_anomaly_tests.csv
data/processed/fixed_nuts2016_snow_2021_anomaly_models.csv
data/processed/fixed_nuts2016_snow_levin_jan_july_summary.csv
data/processed/fixed_nuts2016_snow_levin_percent_lit_thresholds.csv

figures/population_relationship/fixed_nuts2016_snow_viirs_test/february_snow_vs_viirs_anomaly.pdf
figures/population_relationship/fixed_nuts2016_snow_viirs_test/snow_viirs_anomaly_rho_jan_feb_mar.pdf
figures/population_relationship/fixed_nuts2016_snow_viirs_test/levin_style_jan_july_snow_radiance_2021.pdf
figures/population_relationship/fixed_nuts2016_snow_viirs_test/levin_style_percent_lit_thresholds.pdf
```

## Current scientific state

The 2021 Annual VNL pulse is robust to fixed NUTS geography, exists in raw data, survives the annual median product, includes both mask/support and persistent-pixel intensity components, reverses strongly in 2022, and localizes mainly to January-February 2021. The July-specific stray-light-LUT explanation is insufficient.

Snow/albedo is the current leading seasonal mechanism to test, but **there is not yet a valid snow result**. The first apparently strong snow correlation was invalidated by the MOD10CM categorical-code leak.

See:

- `docs/investigations/viirs_2021_anomaly.md`
- `docs/lab_notebook/2026-10-13.md`
