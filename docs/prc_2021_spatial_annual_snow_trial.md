# Spatial 2021 annual VIIRS snow/albedo correction: prespecified pilot

**Goal:** after the locked 17-site PRC transfer and the 18-region monthly-to-annual bridge, attempt a scientifically motivated **native VIIRS-pixel spatial correction** to the 2021 annual \`average_masked\` EOG VNL V2 radiance layer. Do not use a single NUTS3/year scalar or adjust snow coefficients to smooth PRC. Previous pilot / PRC outputs stay intact.

## Annual/monthly bridge (completed)

18 independent selected NUTS2016 NUTS3 pilot regions, 36 valid region-years. Official annual lit-pixel average vs twelve-month \`NPP VCMSLCFG\` CF-weighted monthly proxy:

| Statistic | 2020 | 2021 |
|---|---:|---:|
| Annual radiance | 5.715240 | 6.630235 |
| Monthly-weighted proxy | 5.620095 | 6.768495 |
| Monthly proxy − annual | −0.095146 | +0.138259 |
| RMSE proxy − annual | 0.471670 | 0.536027 |
| Pixels >20% proxy error | 0.0013 | 0.0008 |
| Jan–Feb fraction of cloudfree count | 0.1639 | 0.1178 |

Source EOG Annual VNL V2 method (https://eogdata.mines.edu/products/vnl/; Elvidge et al. Remote Sensing 2021 doi:10.3390/rs13050922): 12 monthly cloud-free count weighted average, with median-based outlier rejection and multiyear background masking. January/February observations constitute ~11.8% of the 2021 annual cloud-free count across these lit-pixel samples. The monthly proxy generally resembles annual on lit pixels, but note significant coherent **2020 negative vs 2021 positive monthly–annual source bias**, and higher relative residuals for selected highly lit centres. The official 2021 annual image may have already omitted some anomalous winter months. A reconstruction from monthly radiance is NOT an exact historical source reprocessing. These limitations are material to validation of candidate annual effects.

## Candidate correction (fixed β; spatial VIIRS pixels)

January and February VCMSLCFG NPP source monthly raw radiance \`R_m\` and cloudfree count \`CF_m\` are on a grid mapped to the annual pixels. For each winter month, use MOD10CM snow (0..100%) and VNP43MA3.002 WSA visible QA01 albedo (0..1) in **2020, 2021 and 2022** at that geographic pixel:

\`\`\`text
ΔS_m = max(0, S_2021,m - (S_2020,m+S_2022,m)/2) / 100
ΔA_m = max(0, A_2021,m - (A_2020,m+A_2022,m)/2)
G_m = 1 + 1.626479 ΔS_m + 3.788066 ΔA_m

P = sum(monthly_2021_R_m × CF_m over Jan–Dec)
  / sum(monthly_2021_CF_m over Jan–Dec)
W = sum(R_m × CF_m × (1-1/G_m) over Jan–Feb)
  / sum(monthly_2021_CF_m over Jan–Dec)
f = clamp(W/P,0,1)

R_2021_annual_adj = R_2021_annual_original_raw × (1-f)
\`\`\`

The ratio form rescales the **official annual 2021 raw** radiance using a physically estimated *fraction* from the observed winter monthly data. This avoids directly replacing annual VNL V2 with a different V1 monthly proxy and guarantees \`0<=R_adj<=R_original\`. Preserve annual zero-background mask, original nodata, and unchanged years. Apply 2017+ baseline **−0.235 once only** in the subsequent locked \`PRC(R)\` conversion, as the fitted snow model used **raw VCMSLCFG monthly radiance**. Pixels with missing winter QA for a month receive no correction for that month. Do not impute visible albedo=0 if missing.

**Scientific warning:** \`f\` is the fraction of the raw 12-month monthly proxy radiance attributable to *candidate* snow/albedo excess. It is not necessarily the actual fraction retained in the annual VNL V2 image after median outlier rejection. This is an operational first approximation, requiring physical/spatial and historical validation; not a reimplementation of EOG outlier filtering.

## EcoVis pilot

\`\`\`bash
git pull --ff-only
python src/prc_2021_spatial_annual_snow_trial.py --self-test
python src/prc_2021_spatial_annual_snow_trial.py --limit 6 --export-region-rasters
python src/prc_2021_spatial_annual_snow_trial.py --limit 0 --export-region-rasters
\`\`\`

Outputs under \`data/processed/prc_2021_spatial_annual_snow_trial/\`:

- \`spatial_trial_summary.txt\`: QA, mean proposed annual radiance correction and change in locked17 2021 PRC
- \`spatial_trial_region_diagnostics.csv\`: all pilot NUTS3 region statistics
- \`spatial_trial_source_audit.json\`: annual/moon/monthly/CF/MOD10CM/VNP43MA3 paths and frozen coefficients
- \`pilot_region_pixels/*_2021_snow_albedo_trial.tif\`: five bands (uncorrected annual raw, corrected annual raw, removed winter excess, uncorrected locked17 PRC, candidate corrected locked17 PRC), cropped to each NUTS2016 polygon's bounding box and masked outside polygon

Expected scientific decision: if seasonal signal and QA are spatially coherent, generate the actual **whole-Europe 2021 corrected annual VNL raster** once and run locked PRC 2014–2024. If the predicted annual correction is implausibly high/low or spatially unsupported, first inspect whether official median-based outlier filtering causes mis-specified annual snow contribution. Do not make up an annual correction just to erase the 2021 jump. Fallback (if needed): consistent March–December VIIRS sample for **comparison years too**; excluding only Jan/Feb 2021 biases seasonal comparisons.

No 2014–2024 output is overwritten by this pilot and no annual corrected raster is yet accepted.
