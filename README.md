# European VIIRS Light Pollution Analysis

Reproducible workflow for analysing annual VIIRS night-time light radiance across European NUTS3 administrative regions and relating radiance to population density.

## Objectives

- Process annual VIIRS night-time light radiance data.
- Calculate spatially appropriate radiance statistics for European NUTS3 regions.
- Combine radiance with annual NUTS3 population-density data.
- Analyse changes in the population–radiance relationship through time.
- Produce publication-quality spatial and statistical figures.

## Current design

The pipeline uses a **fixed NUTS 2024 level-3 geography** across all analysis years so that administrative boundary changes are not confused with temporal change.

The target demographic source is Eurostat dataset **`DEMO_R_D3DENS`**.

The production VIIRS zonal-statistics step will use an area-aware method that accounts for both polygon overlap and the changing physical area of geographic-grid cells with latitude. That method is intentionally being validated before it is enabled in the pipeline.

## Pipeline stages

```text
VIIRS annual rasters
        |
        v
validate inputs
        |
        v
NUTS3 radiance extraction
        |
        +----------------------+
        |                      |
        v                      v
radiance table          Eurostat population
        |                      |
        +----------+-----------+
                   |
                   v
              merged table
                   |
                   v
          maps + statistical figures
```

The intended analytical unit is:

```text
year × NUTS_ID
```

## Project structure

- `config/pipeline.yaml` — central workflow configuration
- `src/run_pipeline.py` — pipeline entry point
- `src/validate_inputs.py` — validate local VIIRS and NUTS inputs
- `src/extract_nuts3_radiance.py` — NUTS3 radiance extraction stage
- `src/fetch_eurostat_population.py` — Eurostat population stage
- `src/merge_radiance_population.py` — merge stage
- `src/viirs_annual_visualisation.py` — annual VIIRS visualisation
- `data/raw/viirs/` — local source VIIRS imagery
- `data/raw/nuts/` — local NUTS administrative boundaries
- `data/raw/eurostat/` — local population/demographic source data
- `data/processed/` — derived analytical datasets
- `figures/` — generated figures
- `docs/pipeline.md` — pipeline and methodological notes
- `tests/` — automated tests

Large source datasets and generated outputs are not stored in Git.

## First local check

After installing dependencies and placing source files in the configured data directories:

```bash
python src/run_pipeline.py validate
```

This reports the VIIRS files found, their years and raster metadata, the expected NUTS3 file, and missing analysis years.
