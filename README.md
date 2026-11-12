# European VIIRS Light Pollution Analysis

Reproducible workflow for analysing annual VIIRS night-time light radiance across European NUTS3 administrative regions and relating radiance to population density.

## Objectives

- Process annual VIIRS night-time light radiance data.
- Calculate spatially appropriate radiance statistics for European NUTS3 regions.
- Combine radiance with annual NUTS3 population-density data.
- Analyse changes in the population–radiance relationship through time.
- Produce publication-quality spatial and statistical figures.

## Current production design

The production workflow now covers **2013-2024** and uses the Eurostat-aligned
NUTS3 coding represented by the current `DEMO_R_D3DENS` series:

- 2013-2020: NUTS 2016
- 2021-2022: NUTS 2021
- 2023-2024: NUTS 2024

VIIRS radiance is extracted using exact polygon overlap with spherical physical
cell-area weighting. From 2017 onward, the production calibration applies a
Europe-wide **pixel-level** additive correction of 0.235 nW cm^-2 sr^-1 to
positive Annual VNL V2.2 pixels, clipping corrected positive pixels at zero,
before calculating the NUTS3 mean. Background-zero pixels remain zero and the
raw radiance is retained for provenance. Population-density gaps are recovered
only for validated countries and only where published Eurostat density is
absent.

The pipeline now includes annual country fits, longitudinal mixed-effects
models with NUTS3 random intercepts, all-year country facets, and publication
maps of corrected NUTS3 radiance across Europe.

See [the 2026-10-08 production milestone](docs/milestone_2026-10-08.md) and
[the detailed pipeline notes](docs/pipeline.md).

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
- `docs/lab_notebook/README.md` — chronological research lab notebook, including exploratory and negative results
- `tests/` — automated tests

Large source datasets and generated outputs are not stored in Git.

## First local check

After installing dependencies and placing source files in the configured data directories:

```bash
python src/run_pipeline.py validate
```

This reports the VIIRS files found, their years and raster metadata, the expected NUTS3 file, and missing analysis years.


## Original-sky Drassodes reconstruction (experimental)

The independent Python reconstruction of the original PM-eye contrast has a
read-only EcoVis image-folder discovery command, explicit TIFF/Stokes inputs,
a 60-degree-radius equal-pixel eye integration and analytic maximum/minimum
orthogonal-eye log contrast. It preserves raw images and the locked European
VIIRS/PRC analyses. Synthetic checks pass; original-field reproduction is
pending discovery, image geometry and alignment review on EcoVis.
See [the reconstruction guide](docs/spider_prc_reconstruction.md) for setup,
input conventions, processing alternatives, outputs and run commands, and the
[12 November lab notebook](docs/lab_notebook/2026-11-12.md) for decisions and QA.
