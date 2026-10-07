# Pipeline design

## Scientific objective

Create a reproducible annual dataset linking VIIRS night-time radiance to
population density for European NUTS3 regions.

The analytical unit is:

`year × NUTS_ID`

The intended final table contains, at minimum:

- year
- NUTS_ID
- country code
- NUTS3 name
- mean annual VIIRS radiance
- median annual VIIRS radiance
- population density
- provenance fields for the VIIRS product, NUTS release and demographic source

## Fixed geography

The initial workflow uses a single NUTS 2024 level-3 geography for every
analysis year. This prevents administrative boundary changes from being
confounded with temporal changes in radiance or population density.

## Pipeline stages

1. **validate**
   - discover local VIIRS annual rasters
   - inspect CRS and raster dimensions
   - confirm NUTS3 geometry
   - report missing analysis years

2. **radiance**
   - calculate annual NUTS3 radiance summaries
   - use exact polygon overlap
   - use physically appropriate area weighting
   - preserve source nodata/mask semantics
   - write a tidy Parquet table

3. **population**
   - retrieve Eurostat `DEMO_R_D3DENS`
   - retain NUTS3 annual population density
   - preserve Eurostat flags/status where useful
   - write a tidy Parquet table

4. **merge**
   - join on year and NUTS_ID
   - report unmatched regions
   - write Parquet and CSV analysis tables

5. **figures**
   - annual VIIRS small multiples
   - annual radiance change
   - NUTS3 choropleths
   - population-density versus radiance relationships
   - residual maps showing regions brighter/darker than expected for population

## Important methodological decision

The production radiance extractor is intentionally not implemented in the
initial skeleton. VIIRS is supplied on a geographic latitude/longitude grid,
so cell areas vary with latitude. Before production extraction, the workflow
will validate a method that accounts for both fractional polygon overlap and
physical cell area rather than using a simple unweighted mean of raster cells.

## Local data

Large raw and generated datasets are deliberately excluded from Git. Expected
local locations are defined in `config/pipeline.yaml`.
