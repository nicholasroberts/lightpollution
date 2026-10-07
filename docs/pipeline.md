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

The workflow uses a single NUTS 2024 level-3 geography for every analysis
year. This prevents administrative boundary changes from being confounded
with temporal changes in radiance or population density.

## European study area

The GISCO NUTS file includes some geographically distant overseas territories.
The production analysis therefore applies an explicit geographic study-area
filter before extraction.

A NUTS3 region is retained when its representative point lies within the
configured longitude/latitude envelope:

`[-32.0, 27.0, 45.0, 72.5]`

This broad envelope deliberately retains the Azores, Madeira, Canary Islands,
Cyprus and all of Türkiye, while excluding distant overseas regions such as
French Guiana, the Caribbean departments, Réunion and Mayotte. The rule is
stored in `config/pipeline.yaml` and the extractor reports all excluded
regions at runtime.

## Pipeline stages

1. **validate**
   - discover local VIIRS annual rasters
   - inspect CRS and raster dimensions
   - confirm NUTS3 geometry
   - report missing analysis years

2. **weighting-test**
   - select representative NUTS3 regions spanning latitude, area and shape
   - compare pixel-centre, exact fractional-overlap and physical-area-weighted means
   - write a diagnostic CSV and figure
   - use the result to choose the production zonal-statistics method

3. **radiance**
   - calculate annual NUTS3 radiance summaries
   - use exact polygon overlap
   - use physically appropriate area weighting
   - preserve source nodata/mask semantics
   - write a tidy Parquet table

4. **population**
   - retrieve Eurostat `DEMO_R_D3DENS`
   - retain NUTS3 annual population density
   - preserve Eurostat flags/status where useful
   - write a tidy Parquet table

5. **merge**
   - join on year and NUTS_ID
   - report unmatched regions
   - write Parquet and CSV analysis tables

6. **figures**
   - annual VIIRS small multiples
   - annual radiance change
   - NUTS3 choropleths
   - population-density versus radiance relationships
   - residual maps showing regions brighter/darker than expected for population

## Important methodological decision

The production radiance extractor uses exact polygon overlap with
`coverage_weight=area_spherical_m2`. VIIRS is supplied on a geographic
longitude/latitude grid, so physical cell area varies with latitude.

A 2024 validation across representative NUTS3 regions compared:

- a pixel-centre mean;
- an exact fractional-overlap mean; and
- exact fractional overlap weighted by spherical cell area.

Boundary treatment changed mean radiance by about 2% for some small/coastal
regions, while physical-area weighting changed the result by about 2.5% for a
large high-latitude region. The area-weighted exact method was therefore
selected for production extraction.

Both mean and median are calculated using the same spherical-area coverage
weight. The output also records the valid VIIRS-covered area in square
kilometres as a quality-control field.

## Local data

Large raw and generated datasets are deliberately excluded from Git. Expected
local locations are defined in `config/pipeline.yaml`.
