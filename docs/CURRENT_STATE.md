# Current research state / restart point

**Updated:** 2026-10-14  
**Branch:** `future-population-models`  
**Machine:** EcoVis  
**Repository root:** `/run/media/bznwr/data/lightpollution`

## Overarching objective: historical reality and model response must converge

The immediate snow/albedo work is subordinate to the main scientific problem:
the future-radiance model must reproduce the historical density-dependent
radiance behaviour before it is used for projection.

The key diagnostic figure is:

```text
figures/population_relationship/future_radiance_historical_validation/
    historical_nuts3_slope_distributions.pdf
```

It compares the distribution of observed fixed-NUTS2016 2013-2020 NUTS3
radiance trends in each population-density class with the canonical
density-class estimates derived from the longer 2013-2024 modelling framework.

The important feature is not merely whether the canonical diamonds lie within
the broad observed distributions. Their *density-dependent pattern* should agree
with the historical class behaviour. At present the figure shows a systematic
tension: the canonical estimates are more positive in sparse classes and
negative in dense classes, whereas the central tendency of the observed
2013-2020 distributions is approximately flat/negative in sparse classes and
more positive in dense classes.

Therefore the current analysis sequence is:

1. identify and explain non-secular disturbances that can rotate the estimated
   density response, especially the 2021 winter pulse;
2. quantify whether snow/surface-albedo effects account for those disturbances;
3. construct a historically defensible radiance series / model specification;
4. refit the density-dependent temporal response;
5. regenerate `historical_nuts3_slope_distributions.pdf` using the final model;
6. require qualitative and quantitative agreement between observed historical
   class behaviour and model-predicted class trends before using the response
   for future projection.

The albedo attenuation test is therefore not an end in itself. Its role is to
determine whether the 2021 pulse is an environmental observation effect that
should be modelled/controlled rather than interpreted as secular change in
artificial-light emissions.

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
