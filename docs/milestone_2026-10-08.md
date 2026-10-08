# Milestone: calibrated European NUTS3 light-pollution pipeline

**Date:** 2026-10-08  
**Branch:** `main`

This document records the first end-to-end production milestone of the European
VIIRS / NUTS3 light-pollution analysis.

## Scope reached

The repository now contains a reproducible Python workflow for annual
2013-2024 VIIRS night-time radiance across European NUTS3 regions, joined to
Eurostat population density and analysed through both annual country fits and
longitudinal mixed-effects models.

The core analytical unit is:

```text
year x NUTS_ID
```

Large raw VIIRS rasters and generated outputs remain local and are excluded
from Git.

## Geography

The analysis follows the NUTS coding represented by the current Eurostat
`DEMO_R_D3DENS` series rather than forcing the legally contemporaneous NUTS
release for each historical year.

Configured mapping:

- 2013-2020 -> NUTS 2016
- 2021-2022 -> NUTS 2021
- 2023-2024 -> NUTS 2024

## VIIRS radiance extraction

Production radiance is derived from EOG Annual VNL v2.2
`average_masked` rasters.

The zonal statistic is the exact polygon-overlap, spherical-cell-area-weighted
mean:

```text
mean(coverage_weight = area_spherical_m2)
```

## Population density

Primary source: Eurostat `DEMO_R_D3DENS`.

Where validated values are missing, production recovery is allowed only for
BE, EE, HR, IT and NL using `demo_r_pjanaggr3 / reg_area3` total land area.
Published Eurostat density is always retained where available and recovered
values carry explicit provenance.

## 2017 VIIRS calibration correction

The 2017 radiance discontinuity is treated as an additive zero-point shift in
linear radiance.

Sensitivity tests used the darkest 10%, 15%, 20%, 25% and 30% of matched
2015-2018 NUTS3 regions. The darkest **20%** gave the best combined performance
for removal of the 2016-2017 discontinuity and agreement with the archived
independently corrected 2016-2019 analysis.

Production offsets are applied from 2017 onward for the currently validated
countries DE, ES, FR, IT and NL.

The raw extraction is preserved as `radiance_mean_raw`; production analyses
use `radiance_mean_corrected`. Corrected zero or negative values are retained
in linear space and excluded only from logarithmic analyses.

## Country-level annual fits

For each country and year:

```text
log10(corrected mean radiance) ~ log10(population density)
```

Outputs include n, slope, intercept and R2, plus annual trend plots and
six-country-per-page all-year facet figures in alphabetical order.

## Longitudinal mixed-effects analysis

The full production model is:

```text
log10(radiance)
~ log10(population density) * year * country
+ (1 | NUTS_ID)
```

The current analysis contains 16,743 observations from 1,588 repeated NUTS3
groups across 31 countries.

Sequential likelihood-ratio tests show highly significant effects of year,
country, population density x country, year x country, population density x
year, and population density x year x country. The dominant incremental effect
is the **year x country** interaction.

Current full-model summaries:

- Nakagawa marginal R2: ~0.912
- Nakagawa conditional R2: ~0.992
- NUTS3 ICC: ~0.914

Incremental term-block effect sizes are reported using likelihood-ratio
pseudo-R2 and the corresponding LR-based f2.

## Publication figures

Current production figures include:

- annual country log-log fits;
- slope/intercept/R2 trends;
- alphabetical six-country-per-page all-year facets;
- annual 600 dpi corrected NUTS3 radiance maps across Europe.

The Europe maps use `radiance_mean_corrected` only, EPSG:3035, one common
logarithmic colour scale across 2013-2024, and a European extent that excludes
overseas territories.

## Main production stages

```text
nuts
validate
radiance
calibration
population
population-recover
merge
fit
trends
country-facets
mixed-model
europe-maps
```

Diagnostic stages for weighting, population recovery and calibration remain in
the repository and are documented in `docs/pipeline.md`.

## Status at this milestone

The project has moved from pipeline construction and validation into a stable
production-analysis phase. This milestone should be treated as the baseline
production implementation against which later methodological changes are
compared.

Likely next steps are detailed interpretation of country-specific temporal
effects, publication-ready mixed-model effect figures, annual spatial-change
analysis, and final methods/reproducibility documentation.
