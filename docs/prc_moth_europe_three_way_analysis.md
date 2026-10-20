# Approved moth PRC: whole-Europe parallel spider, moth and VIIRS comparison

**Date:** 20 October 2026. Branch: \`historical-polarization-visibility\`.

The investigator **approved and locked** the independent moth 16-site observation-year PRC transfer, \`config/prc_transfer_moth.yaml\` with \`status: locked\`, defined as:

\`\`\`text
PRC_moth(R)=0.10128334523974113 + 0.5858115294448945 exp(-0.0782312036944599 R)
PRC_spider(R)=0.1453159020633674 + 0.3108124541438353 exp(-0.1037115607084461 R)
\`\`\`

Both R values are the same appropriately calibrated annual EOG VIIRS VNL DNB radiance, 2017+ offset 0.235 subtracted once per valid positive pixel. PRC is evaluated on **each native VIIRS pixel BEFORE spherical-area averaging** into fixed 2016 NUTS3 regions. Do not equate PRC percentages across taxa with relative behavioural effects: their response ranges are different and both biological transfers are empirical.

## Exact command on EcoVis

\`\`\`bash
git pull --ff-only
python src/run_prc_moth_and_triple_maps.py --self-test
python src/run_prc_moth_and_triple_maps.py
\`\`\`

The first three self-tests are light-weight, including EOG pixel-grouping QA and synthetic original-pixel 2021 patch overlay. The production command:

1. Runs the **entire moth 2014–2024 11-year NUTS3 pixel extraction** using the already established \`src/historical_prc_pixel_nuts3.py\` with \`--transfer-config config/prc_transfer_moth.yaml\`. No changes to the original spider PRC code or source VIIRS.
2. Generates parallel moth annual PRC CSV/Parquet, European time series, ranked NUTS3 table and matching PRC maps/figures to the spider 17-site output.
3. Reuses the previously generated spider+calibrated-radiance 2014/2024 NUTS3 CSV and moth 2014/2024 endpoints to generate **two THREE-PANEL maps** (radiance, spider PRC, moth PRC) for absolute and percentage 2014→2024 change, and a **country-faceted spider/moth paired scatter PDF**, and a merged NUTS3 CSV. It checks the cohort, identical pixel-area coverage and radiance offset. Uses symmetric independent 98th-percentile scales. Dim 2014 radiance <0.2 nW cm-2 sr-1 has undefined radiance percentage displayed grey; spider and moth percentage PRCs are separately defined.
4. Applies existing 2021 native-Europe snow/albedo candidate patch to original global 2021 VIIRS **only while reading** for moth PRC. Writes a *candidate* moth 11-year trajectory and per-NUTS 2021 PRC changes; no raw VNL TIFF or approved production moth atlas overwritten.

### Outputs

**Moth baseline (11 years):**
\`data/processed/historical_prc_2014_2024_locked_moth16/\` — \`prc_nuts3_annual_cumulative.csv\` / \`.parquet\`, \`prc_europe_annual_cumulative.csv\`, \`prc_nuts3_rankings_2014_2024.csv\`; map PDFs in \`figures/historical_prc_2014_2024_locked_moth16/\`.

**Three-way figures:**
\`figures/prc_radiance_spider_moth_2014_2024/\`

- \`compare_absolute_radiance_spider_moth_2014_2024.pdf\` — three absolute change maps for 2014→2024
- \`compare_percent_radiance_spider_moth_2014_2024.pdf\` — three percentage change maps for 2014→2024
- \`scatter_spider_moth_prc_vs_radiance_by_country_2014_2024.pdf\` — by country, overlaid spider and moth ΔPRC vs ΔR
- \`nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv\` — one joined NUTS3 table, 2014/2024 R, two taxon PRCs, absolute and percentage changes, country and QA.

**2021 candidate moth:**
\`data/processed/prc_2021_locked_moth16_combined/\` — 11-year candidate \`prc_nuts3_annual_cumulative_candidate.csv\`, Europe \`prc_europe_annual_cumulative_candidate.csv\`, \`europe_2021_original_vs_candidate.csv\`, \`nuts3_2021_spatial_adjustment.csv\`, \`europe_original_vs_candidate_PRC.pdf\` and \`prc_2021_spatial_overlay_summary.json\`.

The candidate 2021 spatial snow correction remains **NOT LOCKED** (it had removed ~48.5% of the original spider 2021 residual in Europe, without removing all). The source patch is shared by both taxa, but nonlinear responses give potentially different percentage of PRC anomaly explained. Endpoint 2014/2024 maps do **not** depend on the 2021 correction.

**Dependencies already available on EcoVis**: \`figures/prc_radiance_comparison_locked17/nuts3_radiance_prc_cumulative_change_2014_2024.csv\` (calibrated VIIRS radiance and spider PRC endpoints); \`data/processed/prc_2021_spatial_europe_candidate/vnl_2021_snow_albedo_europe_patch.tif\` and \`europe_patch_summary.json\`; source annual VIIRS rasters. If the completed radiance CSV is absent, rerun \`python src/plot_prc_radiance_2014_2024.py\` first; it does not depend on the moth baseline.

**QA criteria**: 1,522 NUTS3 regions/year; 16,742 NUTS-year rows; moth PRC monotonically declining with positive R; identical valid areas between moth and spider at each endpoint, and no alteration to source raster zeros/nodata; year 2021 moth candidate PRC non-decreasing relative to moth original; all years except 2021 unchanged. No unreviewed LED spectral correction enters either response.
