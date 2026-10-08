# Pipeline design

## Scientific objective

Create a reproducible annual dataset linking VIIRS night-time radiance to
population density for European NUTS3 regions.

The analytical unit is an available `year × NUTS_ID` observation. Within each
country and year, matched NUTS3 radiance/population-density pairs are treated
as the observations used for the log-radiance versus log-population-density
relationship. A NUTS3 region does not need to exist in every year.

## Eurostat-aligned NUTS geography

The workflow uses the NUTS classification represented by the current Eurostat
`DEMO_R_D3DENS` population-density time series. Eurostat retrospectively
recasts historical regional statistics onto newer NUTS classifications, so
the geographical coding in the live population dataset does not always equal
the NUTS release formally applicable in the original reference year.

Direct comparison of the population codes against GISCO NUTS3 releases showed
the following best-matching geography:

- 2013-2020: NUTS 2016
- 2021-2022: NUTS 2021
- 2023-2024: NUTS 2024

VIIRS radiance is therefore aggregated using exactly this mapping. This gives
radiance and population density compatible NUTS identifiers within each year
and avoids large artificial losses of observations caused by joining
retrospectively recoded Eurostat population data to obsolete historical NUTS
codes.

This is a deliberate analytical choice. The objective is the country-level
relationship between population density and radiance, so preserving compatible
regional units is more important here than retaining the exact boundary regime
that was legally in force in each historical year.

The required GISCO Level-3 GeoJSON releases are downloaded and retained
locally. Population observations that still do not match their configured
Eurostat-aligned release are flagged and omitted rather than imputed.

## Pipeline stages

1. **nuts**
   - download NUTS 2016, 2021 and 2024 Level-3 GISCO GeoJSON
   - use 01M resolution and EPSG:4326
   - use `wget` when available

2. **validate**
   - discover local VIIRS annual rasters
   - inspect raster CRS/dimensions
   - confirm all required NUTS releases
   - report the Eurostat-aligned year-to-NUTS-release mapping and missing VIIRS years

3. **weighting-test**
   - compare pixel-centre, exact fractional-overlap and spherical-area-weighted
     means on representative NUTS3 regions

4. **radiance**
   - select the Eurostat-aligned NUTS release configured for each VIIRS year
   - calculate exact-overlap, spherical-area-weighted mean radiance
   - retain valid covered area and provenance
   - write the unmodified raw radiance table

5. **calibration**
   - preserve extracted radiance as `radiance_mean_raw`
   - estimate country-specific additive offsets from the darkest 20% of
     matched 2015-2018 NUTS3 regions
   - apply the offset from 2017 onward for DE, IT, NL, FR and ES
   - retain offset, method, threshold and application flag for every row
   - store the result as `radiance_mean_corrected`
   - retain non-positive corrected values without clipping or pseudocounts

6. **population**
   - retrieve Eurostat `DEMO_R_D3DENS`
   - request `geoLevel=nuts3`
   - retain the raw JSON-stat response
   - use `wget` by default
   - associate each year with the NUTS release used by the harmonized Eurostat population series
   - flag whether each population record matches that release
   - preserve Eurostat status flags

7. **population-fallback-test**
   - diagnostic only; does not alter production data
   - test whether missing `DEMO_R_D3DENS` values can be recovered from
     `demo_r_pjanaggr3` population on 1 January divided by `reg_area3`
     total land area
   - compare reconstructed density with published density wherever both exist
   - quantify the expected difference caused by annual-average versus
     1-January population definitions

8. **population-recover**
   - preserve every published `DEMO_R_D3DENS` value unchanged
   - restrict recovery to NUTS3 codes belonging to the configured geography
   - fill only absent/null density observations using
     `demo_r_pjanaggr3 / reg_area3` total land area
   - restrict production recovery to a validated country allow-list
     (currently Belgium, Estonia, Croatia, Italy and the Netherlands)
   - retain explicit provenance and recovery flags for every derived value
   - Norway is excluded from production recovery because validation overlap is
     sparse and several target NUTS3 regions remain unrecovered

9. **merge**
   - inner join on `year + NUTS_ID + nuts_release`
   - omit unmatched/missing observations rather than harmonising or imputing
   - calculate log10 radiance and log10 population-density fields
   - report the number of paired observations for every country/year

10. **fit**
   - fit a separate OLS line in log10-log10 space for every country/year
   - report n, slope, intercept and R²
   - exclude only non-positive values from the logarithmic fit
   - write a fit-summary CSV and a multi-panel comparison figure

11. **benchmark**
   - compare rebuilt country fits with the archived 2016-2019 analysis
   - report old mean/range and new slope/intercept values
   - flag whether the new value lies inside the legacy four-year range
   - write a comparison CSV and benchmark figure

12. **trends**
   - plot country-specific slope, intercept and R² from 2013-2024
   - mark the 2016/2017 VIIRS calibration transition
   - calculate year-to-year changes in slope, intercept and R²
   - write the change table for later calibration analysis

## Radiance statistic

The production radiance extractor uses exact polygon overlap with
`coverage_weight=area_spherical_m2`. VIIRS is supplied on a geographic
longitude/latitude grid, so physical cell area varies with latitude.

The 2024 validation showed that exact boundary treatment changed mean radiance
by about 2% for some small/coastal regions, while spherical physical-area
weighting changed the result by about 2.5% for a large high-latitude region.
The exact-overlap, area-weighted mean is therefore the production radiance
metric.

The median is not used in the production analysis. Testing showed that the
weighted median/quantile returned by the extraction library did not correspond
to the ordinary pixel-value median expected for zero-dominated dark regions,
whereas the mean passed the direct raster checks and is the relevant metric for
the population-density relationship.

## Local data

Large raw and generated datasets are deliberately excluded from Git. Expected
local locations are defined in `config/pipeline.yaml`.


13. **calibration-offset-test**
   - diagnostic only; does not alter production radiance
   - use matched NUTS3 regions across 2015-2018
   - estimate the excess 2016-to-2017 change in linear radiance
   - compare that estimate across 2016-radiance deciles and countries
   - test whether a simple additive zero-point shift is sufficient


14. **calibration-country-test**
   - diagnostic only; does not alter production radiance
   - estimate a country-specific additive 2017 offset from the darkest 20% of
     matched 2016 NUTS3 regions
   - subtract that constant from 2017+ radiance for the test only
   - re-fit annual country relationships and compare the 2016-to-2017
     slope/intercept discontinuity before and after correction


15. **calibration-sensitivity**
   - diagnostic only; does not alter production radiance
   - repeat the country-specific additive offset estimate using the darkest
     10%, 15%, 20%, 25% and 30% of matched regions
   - quantify residual 2016-to-2017 slope/intercept discontinuities
   - compare corrected 2017-2019 changes relative to 2016 against the archived
     legacy analysis, which used an upstream zero-point correction
   - rank thresholds jointly by discontinuity removal and legacy agreement


## Production radiance calibration

The 2017 VIIRS discontinuity is treated as an additive zero-point shift in
linear radiance. Sensitivity analysis using the darkest 10-30% of matched
NUTS3 regions found the 20% threshold gave the smallest residual 2016-2017
slope/intercept discontinuity and the best agreement with the archived
2016-2019 analysis that had an upstream zero-point correction.

The production estimator uses matched 2015-2018 regions within each configured
country. The darkest 20% are selected by 2016 radiance. The expected natural
2016-2017 change is estimated from the mean of the median 2015-2016 and
2017-2018 changes, and the excess is treated as the additive calibration
offset.

The correction is applied from 2017 onward. Raw and corrected radiance are
both retained. Corrected values are never clipped or replaced with a
pseudocount; non-positive values remain in the linear dataset and are excluded
only from logarithmic fitting.

The production correction is currently validated for DE, IT, NL, FR and ES.
Countries outside that configured set retain their raw radiance unchanged until
their calibration correction is separately validated.


16. **country-facets**
   - generate a multi-page publication PDF with six countries per page
   - order countries alphabetically by country code
   - show all available years in every country panel
   - plot corrected NUTS3 radiance versus population density in log10-log10 space
   - use common axes across all pages for direct comparison

17. **mixed-model**
   - fit a Gaussian linear mixed-effects model to corrected log10 radiance
   - use NUTS3 as a random intercept
   - include fixed effects of log10 population density, centred year and country
   - include country-specific population-density slopes and country-specific
     temporal trends
   - fit nested ML models and likelihood-ratio tests for the added effects
   - report the NUTS3 random-intercept variance and ICC
   - exclude countries represented by fewer than three unique NUTS3 regions
