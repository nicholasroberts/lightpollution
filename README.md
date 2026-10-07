# European VIIRS Light Pollution Analysis

Analysis workflow for annual VIIRS night-time light radiance across European NUTS3 administrative regions.

## Objectives

- Process annual VIIRS night-time light radiance data.
- Calculate spatially appropriate mean radiance for European NUTS3 regions.
- Combine radiance with NUTS3 population-density data.
- Analyse changes in the population–radiance relationship through time.
- Produce publication-quality spatial and statistical figures.

## Project structure

- `src/` — analysis and processing code
- `data/raw/viirs/` — source VIIRS imagery
- `data/raw/nuts/` — NUTS administrative boundaries
- `data/raw/eurostat/` — population and demographic data
- `data/processed/` — derived analytical datasets
- `figures/` — generated figures
- `config/` — workflow configuration
- `notebooks/` — exploratory analyses
- `tests/` — automated tests
- `docs/` — documentation

Large source datasets and generated outputs are not stored in Git.
