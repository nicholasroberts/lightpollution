# 2021 snow/albedo effect on historical PRC: regional-to-pixel bridge test

**Aim:** decide whether the known, localized 2021 snow/surface-albedo disturbance can be corrected adequately at NUTS3 level, or whether it requires a spatially varying pixel correction. This is a **validation stage**, not approval to alter the historical PRC data.

## Two different corrections, kept separate

- **Pre-existing 2017-onward radiance calibration:** project-established subtraction of **0.235 nW cm^-2 sr^-1 from positive EOG VIIRS pixels**, with clipping to zero. All historical PRC data already use this; do not double-apply it. The project's calibration is an empirical solution to the cross-2016/2017 product discontinuity; it should not be represented as a second official NASA correction.
- **2021 surface-state bias:** measured anomaly strongest in January-February, associated with snow cover and albedo, with brightness interaction and geographical variability. It is not a global offset and is treated separately from the instrument calibration. A snow-driven upward satellite-radiance response is not necessarily an increase in luminaires or the animal's actual sky-polarization exposure.

Original regional correction: \`src/build_fixed_nuts2016_2021_surface_correction.py\`, output \`data/processed/fixed_nuts2016_2021_surface_correction/fixed_nuts2016_viirs_2014_2024_surface_corrected.parquet\`, fields \`radiance_mean_corrected\`, \`radiance_mean_corrected_surface_logratio_primary\`, \`radiance_mean_corrected_surface_logratio_snow\`.

**Important:** this original correction is empirical at NUTS2016 regional radiance level. It did not establish the within-region snow footprint or exact 2021 pixel-level counterfactual.

## Test 1 — scientifically controlled bridge (implemented)

For NUTS3 region \`g\`, take its existing annual radiance counterfactual \`R*_g,2021\` and unmodified pixel-calibrated \`R_g,2021\`:

\`\`\`text
q_g = R*_g,2021 / R_g,2021
R*_pixel,2021 = q_g * R_pixel,2021
mean_PRC*_g,2021 = area_weighted_mean( f(R*_pixel,2021) )
\`\`\`

The scaling is performed *after* the existing 0.235 pixel correction, *before* calculating nonlinear PRC and *before* area-weighted NUTS means. A fixed within-region \`q_g\` is a **model assumption**. It is NOT a snow-localized correction and cannot be presented as a corrected 500 m VIIRS layer.

\`src/prc_2021_surface_bridge.py\` uses fixed NUTS2016 polygons and the existing local 2021 annual \`average_masked\` raster. For each region it recomputes unmodified regional mean radiance and PRC first, failing if either differs materially from the previously stored extraction. It then applies the primary existing 2021 NUTS snow-plus-albedo model, and independently evaluates the snow-only alternate. Regions without usable correction or with extreme factors remain explicitly ineligible, not silently counted as corrected.

**Test first (ten regions, smoke test):**

\`\`\`bash
git pull --ff-only
python src/prc_2021_surface_bridge.py --limit 10
\`\`\`

Then, if the baseline validations pass:

\`\`\`bash
python src/prc_2021_surface_bridge.py
\`\`\`

Result folder: \`data/processed/historical_prc_2021_surface_bridge/\`

- \`prc_2021_region_bridge.csv\` — NUTS3 baseline and candidate corrected 2021 mean PRC, radiance factor, eligibility, 2020/22 flanking-year PRC anomalies;
- \`prc_2021_bridge_comparison.csv\` — original/primary/snow-only regional summary for the **same eligible regions**;
- \`prc_2021_bridge_by_country.csv\` — coverage and within-country diagnostics;
- \`prc_2021_bridge_diagnostics.pdf\` — 2021 anomaly before/after regional correction; correction vs radiance;
- \`prc_2021_bridge_summary.txt\` — model and missing-coverage notes.

**Caution:** the 2020/22 flank mean is a diagnostic, *not* ground truth. Product consistency around 2022 is not yet established and real lighting changes can occur. The aim is to assess whether corrections are physically constrained and robust, not to force 2021 onto a straight line. Also inspect country/cold-versus-warm region contrasts and the optical plausibility of fitted factors.

## Test 2 — snow-localized pixel correction (not implemented yet)

To justify native 500 m–1 km correction:

1. Select pilot regions with pronounced 2021 snowfall anomalies (Nordic and Central Europe) and regions without (UK/France/southern Europe) using previously downloaded MOD10CM and VNP43MA3 snow/albedo data.
2. On *identical valid VIIRS pixels* for **the same calendar month** in 2020, 2021 and 2022, calculate excess 2021 observed radiance, controlling for changes in cloud-free counts, lighting trends and 2022 product/platform differences. Winter January and February are prespecified; March and summer are negative controls.
3. Co-register MODIS snow and VIIRS pixels carefully (different resolutions/grid). Use QA-filtered fractional snow/albedo and pre-registered physically interpretable spatial supports and brightness strata.
4. Derive a pixel-level winter surface-state factor (snow/albedo effect conditional on baseline radiance) with independent **spatial hold-out regions**, not by matching 2021 to 2020/22 everywhere. Verify the inferred effect persists on stable, lit pixels and collapses on snow-free controls. Quantify uncertainties.
5. Translate the monthly pixel effect to the annual **actual VIIRS aggregation**, including cloud-free observation weighting. Do **not** assume equal month weighting or that the annual EOG layer is exactly a simple monthly arithmetic mean. Previous NUTS3 reconstruction demonstrated the need to estimate the bridge and validate on control years.
6. Cross-check whether the aggregation of pixel corrections reproduces the NUTS3 annual radiance counterfactual **without forcing it to**. Apply the accepted PRC function pixel-wise; measure discrepancies from regional-only scaling and identify where a NUTS approximation is or is not defendable.

Decision gate: if the region-level approximation agrees with the independently obtained pixel correction within its uncertainty for every tested brightness/snow stratum and NUTS regions, a **NUTS-level reporting correction** may be acceptable. If not, use the spatially resolved estimate or retain raw-2021 PRC with a bias/uncertainty flag. In no case should the old regional correction be advertised as a known pixel-level correction.

## Relationship to 2014–2024 findings

The endpoint PRC difference \`PRC_2024 - PRC_2014\` is **unchanged** by correcting only 2021. This issue directly affects 2020→2021, 2021→2022, cumulative annual sequences that pass through 2021, and any 11-year temporal models, but **not** the 2014→2024 two-endpoint difference. Note: cumulative *level relative to 2014* at 2024 is still unchanged; the path to it is altered in 2021.

Never overwrite \`prc_nuts3_annual_cumulative.csv\` or publish modified historical trends until correction choice and uncertainties have been reviewed.
