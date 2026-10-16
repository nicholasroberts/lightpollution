# Historical polarization-cue visibility — proposed analysis pipeline

**Decision:** 16 October 2026. **Status:** research plan only; site-data validation and European mapping have not been run.
**Working geography:** Europe, 2014–2024, with local radiance analysis at approximately 500 m–1 km.
**Branch:** a new research branch is proposed for Monday 19 October; not yet created.

## Pipeline overview

```mermaid
flowchart TD
    A[Original field-site PRC measurements and fitted functions]
    B[Original Light Pollution Map VIIRS 2019 values]
    C[Direct 2019 EOG VIIRS and NASA Black Marble extraction]
    D[Check sites, radiance units, QA and spatial sampling]
    E[Compare annual, monthly and quality-screened daily radiance]
    F[Reproduce and validate the PRC-radiance relationship]
    G{Is biological transfer robust enough?}
    H[Reassess surrounding light, spectra, conditions and uncertainty]
    I[Harmonise 2014–2024 local annual radiance series]
    J[Apply species-specific PRC function for each place and year]
    K[Map annual cue-contrast change and biological exposure]
    L[Compare data products and contrasting habitats including NL333]
    M[Future work: additional taxa, world regions and SSP forecasts]
    A --> F
    B --> D
    C --> D
    D --> E --> F --> G
    G -->|Not yet| H --> F
    G -->|Validated with limits| I --> J --> K --> L --> M
```

## Stage 1: reproduce the original site-level calibration

The original PRC-site predictor was explicitly confirmed to be the **VIIRS 2019 radiance layer** (`viirs_2019`) on lightpollutionmap.info, **not** its modelled sky-brightness layer. The website's present product version may differ from the one used during the original analysis.

Obtain the original site coordinates, field observation dates and times, published 2019 VIIRS values, raw PRC observations and species, measurement uncertainty, and fitted functions. Compare the **same site** against directly extracted 2019 EOG Annual VNL radiance and 2019 NASA Black Marble VNP46A4 snow-free annual radiance, tracking native grid, QA, valid-observation counts, instrument version, radiance units and spatial footprint.

Reproduce the published fit unchanged first. Then compare physically specified pixel and local-neighbourhood predictors, using site-level validation rather than choosing a buffer simply for the best fit. Artificial skyglow can arise from sources beyond the containing VIIRS pixel.

## Stage 2: choose biologically meaningful time support

| Radiance measure | Planned purpose | Limitation |
|---|---|---|
| **2019 annual** | Exact published-layer benchmark | May not represent the experimental date or hour |
| **Monthly / active-season** | Primary alternative matched to field dates | Weather, season and sparse valid observations |
| **Daily, quality screened** | Targeted tests of specific field dates and short-term variability | Overpass timing; gaps and interpolated daily estimates |

Evaluate variability at each original site and whether any alternative improves held-out PRC prediction. Keep processing-product choice separate from averaging-period choice.

For the European retrospective maps, **annual** 2014–2024 composites remain the initial candidate, with a pilot comparison of `f(annual radiance)` against averaged monthly `f(monthly radiance)`. Because PRC is nonlinear, those two are not generally equal.

## Stage 3: construct the historical 2014–2024 radiance fields

Compare EOG Annual VNL V2.2 and NASA Black Marble annual snow-free composites across all 11 years. Their radiance processing is different, but they rely on related VIIRS measurements; this is **not independent sensor replication**. Audit temporal instrument changes, 2017+ EOG pixel calibration, the 2021 snow anomaly, available quality indicators and gaps.

**Critical constraint:** the existing 2021 surface correction was estimated for annual *NUTS3 aggregates*. Do not treat it as a validated, spatially uniform adjustment for 500 m pixels. Develop pixel-scale 2021 handling and quantify the resulting uncertainty.

## Stage 4: calculate local PRC and map historical change

After validating the source radiance–PRC relationship, apply a species-appropriate function at each valid local observation in each year. Propagate radiance and transfer-function uncertainty; map continuous PRC and its change before aggregating to habitat area, NUTS3, countries and Europe. Nonlinearity means `mean(f(R_pixel))` generally differs from `f(mean(R_pixel))`.

Evaluate substantial versus biologically consequential radiance changes, including extremely bright NL333 / Westland and dark/intermediate-radiance sites. Report *continuous receptor contrast* separately from behavioural cue detectability; the latter requires independently supported thresholds.

## Stage 5: limitations and extensions

Explicitly assess blue/LED spectral under-response of VIIRS, non-local atmospheric skyglow, variation in cloud/aerosol and emission geometry, receptor spectral/polarization properties, and behavioural thresholds. Do not assume one species' fitted response generalises globally. Broader taxa/geography and SSP projections remain **subsequent research**, not prerequisites for the historical study.

## Restart instructions

Begin with the original field-site PRC table and fitted function. Do not code a Europe-wide PRC mapper before the site-level predictor is reproduced and evaluated. Do not create the new branch until Monday as agreed.

See [`lab_notebook/2026-10-16.md`](lab_notebook/2026-10-16.md) for the preceding VIIRS and density analyses, scientific reasoning, figure names and data provenance.
