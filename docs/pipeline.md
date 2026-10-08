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
   - write tidy Parquet and CSV tables

5. **population**
   - retrieve Eurostat `DEMO_R_D3DENS`
   - request `geoLevel=nuts3`
   - retain the raw JSON-stat response
   - use `wget` by default
   - associate each year with the NUTS release used by the harmonized Eurostat population series
   - flag whether each population record matches that release
   - preserve Eurostat status flags

6. **merge**
   - inner join on `year + NUTS_ID + nuts_release`
   - omit unmatched/missing observations rather than harmonising or imputing
   - calculate log10 radiance and log10 population-density fields
   - report the number of paired observations for every country/year

7. **fit**
   - fit a separate OLS line in log10-log10 space for every country/year
   - report n, slope, intercept and R²
   - exclude only non-positive values from the logarithmic fit
   - write a fit-summary CSV and a multi-panel comparison figure

8. **benchmark**
   - compare rebuilt country fits with the archived 2016-2019 analysis
   - report old mean/range and new slope/intercept values
   - flag whether the new value lies inside the legacy four-year range
   - write a comparison CSV and benchmark figure

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
