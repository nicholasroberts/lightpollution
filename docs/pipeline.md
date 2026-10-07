# Pipeline design

## Scientific objective

Create a reproducible annual dataset linking VIIRS night-time radiance to
population density for European NUTS3 regions.

The analytical unit is an available `year × NUTS_ID` observation. Within each
country and year, matched NUTS3 radiance/population-density pairs are treated
as the observations used for the log-radiance versus log-population-density
relationship. A NUTS3 region does not need to exist in every year.

## Year-specific NUTS geography

The workflow deliberately uses the NUTS classification applicable to each
reference year rather than forcing all years onto NUTS 2024:

- 2013-2014: NUTS 2010
- 2015-2017: NUTS 2013
- 2018-2020: NUTS 2016
- 2021-2023: NUTS 2021
- 2024: NUTS 2024

This follows Eurostat's official applicability periods. Boundary changes,
splits, mergers and code changes between NUTS releases therefore do not require
harmonisation for this analysis: each annual country-level relationship uses
the valid NUTS3 observations available for that year.

The required GISCO Level-3 GeoJSON releases are downloaded and retained
locally. Missing Eurostat observations or identifiers that do not match that
year's NUTS release are omitted rather than imputed.

## Pipeline stages

1. **nuts**
   - download NUTS 2010, 2013, 2016, 2021 and 2024 Level-3 GISCO GeoJSON
   - use 01M resolution and EPSG:4326
   - use `wget` when available

2. **validate**
   - discover local VIIRS annual rasters
   - inspect raster CRS/dimensions
   - confirm all required NUTS releases
   - report year-to-NUTS-release mapping and missing VIIRS years

3. **weighting-test**
   - compare pixel-centre, exact fractional-overlap and spherical-area-weighted
     means on representative NUTS3 regions

4. **radiance**
   - select the NUTS release applicable to each VIIRS year
   - calculate exact-overlap, spherical-area-weighted mean radiance
   - retain valid covered area and provenance
   - write tidy Parquet and CSV tables

5. **population**
   - retrieve Eurostat `DEMO_R_D3DENS`
   - request `geoLevel=nuts3`
   - retain the raw JSON-stat response
   - use `wget` by default
   - associate each year with its applicable NUTS release
   - flag whether each population record matches that release
   - preserve Eurostat status flags

6. **merge**
   - inner join on `year + NUTS_ID + nuts_release`
   - omit unmatched/missing observations rather than harmonising or imputing
   - calculate log10 radiance and log10 population-density fields
   - report the number of paired observations for every country/year

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
