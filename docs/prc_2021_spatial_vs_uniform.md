# PRC sensitivity: spatial vs region-uniform 2021 correction (first test)

**Purpose:** assess whether snow/albedo effects within a region must be spatially resolved to estimate **mean polarization-receptor contrast**. This is a *January–March 2021 matched monthly pixel* experiment. It is not a new annual 2021 correction.

Existing 18-region 2021 spatial pilot output must be available on EcoVis: \`data/processed/historical_prc_2021_pixel_snow_test/nuts_disjoint_holdout_predictions.csv\`. This holds region-disjoint estimates of the surface-correlated radiance component for each native MOD10CM 0.05° snow-cell × VIIRS brightness class × month, from the model trained on January/February data in other NUTS3 regions.

The new stage \`src/prc_2021_spatial_vs_uniform.py\` reconstructs the *same valid monthly VIIRS pixels* in 2020, 2021, 2022 and preserves the snow-valid mask. On each 2021 VIIRS pixel, it subtracts the previously **out-of-region-predicted** snow/albedo-associated monthly radiance component for that coarse cell / brightness group, clipping negative radiance to zero. The method is spatially informed but **does not claim MOD10CM independently observes snow at 500m**.

For comparison, one factor \`q_gm\` is fitted to the **same region, same month, and same eligible pixels** so that
\`\`\`text
mean_area( q_gm * R_2021,pixel ) = mean_area( max(R_2021,pixel - C_cell,brightness,month, 0) ).
\`\`\`
These are *numerically matched corrected mean monthly radiances*. The mean pixelwise PRC then differs solely because the radiance adjustments occur at different locations within that region. This is methodologically distinct from applying the original 2021 annual NUTS3 correction factor (the monthly temporal/processing support is different). Region-uniform comparison uses multiplicative scaling to parallel the earlier regional model's conceptual assumption.

The selected working transfer is \`PRC(R)=0.144845+0.316330 exp(-0.106483R)\`. Both direct monthly EOG radiance \`R\` and an alternative \`max(R - 0.235, 0)\` are evaluated: **the production 0.235 adjustment is established for the annual VNL pipeline, not independently validated for the monthly products**. The latter scenario is a calibration sensitivity, not an assertion that monthly and annual radiance scales are identical.

**Inputs:** monthly EOG VCMSLCFG \`avg_rade9h\` and \`cf_cvg\` rasters from 2020/2021/2022, native MOD10CM snow files/cache for the same months, out-of-region holdout predictions, NUTS2016 polygons. No downloads. No annual PRC outputs overwritten. Albedo used in the existing holdout models; the new test does **not** remap all albedo layers, but excludes model groups with very low albedo QA coverage (\`--min-group-albedo\`, 0.5 default).

Run:
\`\`\`bash
git pull --ff-only
python src/prc_2021_spatial_vs_uniform.py --self-test
python src/prc_2021_spatial_vs_uniform.py --limit 6
\`\`\`

If matched support is sufficient, proceed:
\`\`\`bash
python src/prc_2021_spatial_vs_uniform.py
\`\`\`

Outputs in \`data/processed/historical_prc_2021_prc_sensitivity/\`:
- \`prc_spatial_vs_equal_mean_uniform_regions.csv\` — NUTS3 / month / offset sensitivity: same raw radiance reduction, two spatial allocations, both mean PRCs, differences, matched support
- \`prc_spatial_vs_equal_mean_uniform_brightness.csv\` — independent brightness-bin diagnostics
- \`prc_spatial_vs_uniform_summary.txt\` — numerical comparisons and coverage
- \`prc_spatial_vs_uniform_2021_pilot.pdf\` — mean-PRC effects for spatial vs equal-mean uniform scaling

**Decision gates:** (i) the predicted component must not rely on very low albedo-support cells, (ii) spatial vs uniform PRC difference must be compared to the original 16-site PRC calibration uncertainty (RMSE ~0.027) *and* the scale of actual historical changes, not only a percent deviation from tiny anomalies, (iii) improvement should be evaluated among dark/intermediate radiance and not merely in bright extreme-value regions, (iv) this is *monthly* counterfactual sensitivity; translating it into an **annual** 2021 PRC correction still requires an independently validated cloud-free weighted monthly→annual bridge and a snow-free/control-year assessment.

**Note on brightness-bin table:** each brightness bin is recomputed with its *own* equal-mean uniform factor, so bin-specific comparisons assess heterogeneity *inside that class*; the region-month comparison is the predeclared headline test.

A standalone test does not establish mechanism or validate region-wide generalization. The 18 regions were selected on snow contrast; results are not yet Europe-wide.
