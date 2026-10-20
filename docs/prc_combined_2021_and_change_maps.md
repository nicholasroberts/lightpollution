# Combined PRC correction and radiance-vs-PRC comparisons

## Objective

Use the **locked 17-site PRC receptor-contrast function**, the already generated **2021 native-grid Europe snow/albedo correction patch**, and fixed NUTS2016 NUTS3 geography. Run two independent analyses together:

1. Apply the candidate 2021 snow/albedo patch on top of the original VNL 2021 radiance pixels **only while reading them**; transform each VIIRS pixel using the locked 17-site PRC function and aggregate area-weighted to **all 1,522 NUTS3 regions**. Replace only 2021 rows in a copy of the existing eleven-year PRC table, then report original vs candidate 2014–2024 European trajectory. Keep the original global annual VIIRS and the existing locked-17 series unchanged. The generated 2021 adjusted dataset is **not yet locked**: inspect coverage and overall snow correction signal first.

2. Use the already locked 17-site PRC NUTS3 endpoint values (2014, 2024), and independently extract area-weighted EOG VIIRS radiance at identical NUTS2016 polygon boundaries for both endpoint years. Apply the historical post-2017 **positive pixel −0.235** adjustment per pixel (only once) before averaging radiance, while PRC has already been transformed per original pixel before averaging. Thus the map comparison shows the nonlinear PRC effect honestly: **mean(PRC(R))**, not **PRC(mean(R))**.

## Generated files

In \`data/processed/prc_2021_locked17_combined/\`:

- \`prc_nuts3_annual_cumulative_candidate.csv\` — 11-year NUTS3 panel, only 2021 changed
- \`prc_europe_annual_cumulative_candidate.csv\` — area-weighted European 11-year trajectory
- \`nuts3_2021_spatial_adjustment.csv\` — per-region original and candidate 2021 PRC
- \`europe_2021_original_vs_candidate.csv\` — original and adjusted annual PRC, year by year
- \`europe_original_vs_candidate_PRC.pdf\` — overlay of original vs candidate PRC history
- \`prc_2021_spatial_overlay_summary.json\` — 2021 effect and residual anomaly vs 2020/22; **\`status:candidate_not_locked\`**

In \`figures/prc_radiance_comparison_locked17/\`:

- \`compare_cumulative_absolute_radiance_vs_prc_2014_2024.pdf\` — one pair of maps, ΔR and ΔPRC
- \`compare_cumulative_percent_radiance_vs_prc_2014_2024.pdf\` — second map pair, %ΔR and %ΔPRC
- \`scatter_prc_vs_radiance_change_by_country_2014_2024.pdf\` — country-faceted NUTS3 ΔPRC vs ΔR, 12 countries/page
- \`nuts3_radiance_prc_cumulative_change_2014_2024.csv\` — joined NUTS3 2014/2024 R, PRC, absolute/% changes, country, NUTS ID, baseline and QA

### Visual conventions

The maps use zero-centered diverging colours, the same geographical extent, and separate robust ±98th-percentile limits for each quantity. These limits **are not numerically comparable across PRC and radiance**; the fields have different units. Grey signifies missing/undefined values, notably for radiance percentage changes in sites with very dim 2014 baseline. The default % radiance baseline threshold is **0.2 nW cm⁻² sr⁻¹**. This threshold is explicitly reported in the title and the complete raw values remain in CSV. Both map pairs use changes calculated as 2024 minus 2014 and **do not depend on the 2021 snow correction**.

Facets use absolute change because percentages can explode near the radiance floor; different panels intentionally permit different axes for clarity (numerical slopes should not be compared visually between panels without reading scales). PRC should generally respond negatively to radiance changes, with nonlinear saturation near bright urban sources.

## EcoVis

\`\`\`bash
git pull --ff-only
python src/run_prc_combined_2021_and_maps.py --self-test
python src/run_prc_combined_2021_and_maps.py
\`\`\`

The Europe patch must already exist at:
\`data/processed/prc_2021_spatial_europe_candidate/vnl_2021_snow_albedo_europe_patch.tif\` with \`europe_patch_summary.json\`.

The locked 17-site full baseline table must already exist at:
\`data/processed/historical_prc_2014_2024_locked17/prc_nuts3_annual_cumulative.csv\`.

The original global EOG annual rasters must remain in \`data/raw/viirs\` per config. The patch is a **Europe-only GeoTIFF**; it is *not* copied to the global annual TIFF. A virtual RasterSource overlay in the exactextract extraction preserves original pixels outside the patch and original NoData/zero-background cells.

After running, first share \`prc_2021_spatial_overlay_summary.json\`, \`europe_2021_original_vs_candidate.csv\`, and the pair-map / scatter PDFs for discussion. If geographic results are scientifically plausible and retained-mask/QA tests pass, **then** decide whether to lock the snow anomaly correction. The model fits themselves are not being refitted or selected to smooth 2021.
