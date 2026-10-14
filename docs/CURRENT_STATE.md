# Current research state / restart point

**Updated:** 2026-10-14  
**Branch:** `future-population-models`  
**Machine:** EcoVis  
**Repository root:** `/run/media/bznwr/data/lightpollution`

## Current conclusion

The corrected MOD10CM snow test has completed and is valid.

The February 2021 anomaly shows:

```text
region Spearman rho = +0.6653
country-median rho = +0.7383
adjusted snow-block partial R2 = 0.2805
adjusted snow-block joint p = 2.48e-14
snow x baseline-brightness p = 1.94e-13
```

January independently supports the relationship; March is a practically null
negative control.

Snow/surface reflectance is therefore strongly supported as a major contributor
to the January-February 2021 VIIRS pulse.

The full 2020-2022 within-region seasonal relationship is not generally
positive, so the conclusion is specific to the winter 2021 event rather than a
universal European snow-brightness rule.

## Next scientific task

Follow the published seasonal-night-light framework one step further by adding:

- MODIS surface albedo, with a visible broadband measure as the primary optical
  variable and white-sky/black-sky variants retained;
- monthly MODIS NDVI as a complementary seasonal surface-state control.

Preserve the existing pre-specification:

```text
February primary
January supportive
March negative control
2021 - mean(2020, 2022)
country + latitude + baseline brightness + cf_cvg controls
surface-state x baseline-brightness interaction
```

Primary mechanistic question:

> Does the albedo anomaly predict the VIIRS anomaly and attenuate the snow term?

This is a mechanistic attenuation test, not formal causal mediation.

## Key validated snow outputs

```text
data/processed/fixed_nuts2016_snow_viirs_test_summary.txt
data/processed/fixed_nuts2016_snow_2021_anomaly_tests.csv
data/processed/fixed_nuts2016_snow_2021_anomaly_models.csv
data/processed/fixed_nuts2016_snow_levin_jan_july_summary.csv
data/processed/fixed_nuts2016_snow_levin_percent_lit_thresholds.csv
data/processed/fixed_nuts2016_snow_levin_seasonal_correlations.csv

figures/population_relationship/fixed_nuts2016_snow_viirs_test/
    february_snow_vs_viirs_anomaly.pdf
    snow_viirs_anomaly_rho_jan_feb_mar.pdf
    levin_style_jan_july_snow_radiance_2021.pdf
    levin_style_percent_lit_thresholds.pdf
```

See:

- `docs/investigations/viirs_2021_anomaly.md`
- `docs/lab_notebook/2026-10-14.md`
