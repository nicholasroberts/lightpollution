# 2021 pixel-snow/albedo pilot: matched VIIRS measurements with spatial holdout

**Status:** proposed observational validation. **Do not alter the calibrated 2014–2024 EOG or PRC history.** This tests whether actual 2021 snow and visible surface-albedo variation spatially explains the 2021 transient VIIRS response, rather than assuming the NUTS3-uniform radiance factor is a true pixel correction.

## Method

1. Select ~18 moderate-size fixed **NUTS2016** polygons objectively using the already-computed **February 2021 MOD10CM snow anomaly**, including snow-increase, snow-decrease and near-zero regional contrasts. Up to three per country. This selection **does not use the VIIRS or PRC response**. The region selection CSV is part of the outputs.
2. For **January, February, March** in **2020, 2021, 2022**, use existing **VCMSLCFG monthly EOG radiance and cloud-free counts**. Require precisely the same VIIRS pixel locations with \`cf_cvg > 0\` and finite, nonnegative radiance in all three years of each month. Reject J01 platform changes for these matched-month tests.
3. Reproject **MOD10CM** 0.05-degree monthly percent snow to the VIIRS grid using **nearest neighbour** without claiming that it becomes an independent 500 m measurement. Preserve original coarse-cell IDs, aggregate VIIRS pixels by that cell and **four fixed radiance classes** to avoid treating all 500 m pixels sharing one coarse snow observation as independent snow observations. Invalid categorical MOD10CM values 211, 250, 253, 254, 255 stay masked.
4. Reproject cached native-sinusoidal **VNP43MA3** visible white-sky albedo \`wsa_qa01_mean\` monthly composites to exactly the same VIIRS grid with nearest neighbour. The existing albedo preprocessing requires original day-level QA 0 or 1, >=25% valid days, and physical 0..1. Albedo remains a *surface measurement*, not nighttime sky polarization. A complete three-year albedo triple is required for the combined model; snow-only remains a separate comparator with broader coverage.
5. Outcome: 2021 minus the mean of 2020 and 2022 radiance **for that same pixel and same month**. Predictors: February/January snow anomaly, albedo anomaly, baseline brightness (including surface-state interactions), and cloud-free count anomaly. Fit weighted models on Jan/Feb; test **whole NUTS3 regions held out**. March is a pre-specified negative control. Report apparent slopes as observational, not proof that snow caused the entire 2021 radiance anomaly.
6. The first-stage output quantifies how spatially informative snow and albedo are. **It does not yet generate a corrected annual EOG 2021 layer or a corrected annual PRC map**: that requires constructing monthly cloud-free-weighted pixel exposure over the whole year, validating the monthly-to-annual mapping against unaffected annual years, and comparing the independently spatial prediction with the previous region-uniform 2021 PRC correction. These are the next design gate, not hidden assumptions.

**Important distinction:** the established annual VNL post-2016 calibration (\`−0.235\` for positive pixels from 2017) remains intact. We do not subtract this offset a second time from monthly 2020–22 pixel contrasts; those are contemporaneous products all within the post-2017 interval.

## EcoVis commands

\`\`\`bash
git pull --ff-only
python src/test_prc_2021_pixel_snow_albedo.py --self-test
python src/test_prc_2021_pixel_snow_albedo.py --limit 6
\`\`\`

If the six-region smoke test has suitable snow/albedo and cloud-free support, use:

\`\`\`bash
python src/test_prc_2021_pixel_snow_albedo.py
\`\`\`

The existing monthly VIIRS inputs, MOD10CM snow archive, and cached VNP43MA3 monthly QA composites must be present; no downloads are performed.

Outputs under \`data/processed/historical_prc_2021_pixel_snow_test/\`:

- \`selected_nuts3_pilot_regions.csv\`
- \`spatial_support_quality.csv\`
- \`matched_snowcell_brightness_groups.csv\`
- \`snow_albedo_model_coefficients.csv\`
- \`nuts_disjoint_holdout_predictions.csv\`
- \`spatial_holdout_diagnostics.csv\`
- \`regional_monthly_spatial_contrast.csv\`
- \`pixel_snow_albedo_spatial_holdout.pdf\`
- \`pixel_snow_albedo_summary.txt\`

**Predeclared decision:** the optical contribution should be directionally plausible in held-out NUTS regions, the February signal should exceed March negative control, observed cloud-free coverage should be sufficient and stable, and the combined snow+albedo model should add useful held-out explanatory value beyond snow alone. Do **not** accept the model based solely on its apparent reduction of the 2021 2020/22 residual: persistent lighting changes and 2022 annual product discontinuities can also affect that comparison. The regional scaling model is a benchmark, not a training target.

No numerical pilot results were observed when this code was written.
