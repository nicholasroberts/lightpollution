# 2021 snow / visible-albedo → VIIRS radiance calibration curve

**Stage:** physically constrained observational calibration, before any snow/albedo correction to annual radiance or historical PRC.

Use already processed matched 2020/2021/2022 same-month VIIRS VCMSLCFG / MOD10CM / VNP43MA3 snow-cell × radiance-stratum aggregates:
\`data/processed/historical_prc_2021_pixel_snow_test/matched_snowcell_brightness_groups.csv\`.

The previously unconstrained regression generated some **negative** predicted snow/albedo effects and therefore *increased* corrected 2021 radiance, and could clip a very large proportion of pixel radiances to zero. Both are inadmissible for the intended operation: remove only additional positive reflectance-related radiance in excess of the local reference surface state, without forcing one calendar year to match a neighbour.

## Prespecified empirical curve

\`\`\`text
G = 1 + beta_s * max(delta_snow_cover_fraction,0)
      + beta_a * max(delta_visible_albedo_fraction,0)
R_2021_snow_albedo_counterfactual = R_2021_observed / G
\`\`\`

Both \`beta_s\` and \`beta_a\` are estimated from data and constrained to \`>=0\`. This transformation guarantees \`0 <= R_counterfactual <= R_observed\` for nonnegative input radiance, never creates negative pixels, and is exactly unchanged when no positive surface-state difference is estimated. A negative snow change relative to the flanking reference is **not corrected upwards** in this deliberately one-sided 2021 brightening adjustment. This does not imply snow can never dim nighttime measurements through other physical mechanisms.

Fit the observed 2021 minus 2020/22 same-month VIIRS radiance on nuisance baseline terms (intercept, February, change in cloudfree coverage, and baseline radiance) **plus** a nonnegative surface component \`baseline_radiance * (G-1)\`. Pixel-albedo quality threshold is 50%, subject to sensitivity. Models: \`snow_only\`, \`albedo_only\`, and \`snow_plus_albedo\`, fitted on *the same QA-complete cell/brightness-class cohort* to avoid misleading differences due to missingness. Fit with robust iterative reweighting and region-month-normalized spatial weights, so that high-coverage very bright regions cannot determine the entire slope by pixel count.

**Holdout and checks:** leave a complete NUTS3 region out of coefficient estimation, then calculate its predicted surface gain without re-fitting. Calculate corrected radiance and deviation relative to 2020/22; **that is validation, not the fitting target used for calibration-parameter selection**. March is withheld from all fitting as a seasonal negative control. Test effects across snow-free, weak snow, strong snow and radiance strata, including the distribution of gains \`G\`, correction magnitudes and whether independent NUTS3 regions improve. A monotone sign constraint alone does not establish transferability. In particular, snow/visible albedo are correlated, the broadband daytime visible albedo product does not directly measure nocturnal VIIRS DNB effective reflectance, and lighting/atmospheric temporal changes are confounders.

The ratio \`R_obs/G\` is *proposed*, not a fully physically identified snow-free light-emission measurement. Recheck gain-vs-observed seasonal/brightness behavior and independent snow-free controls before using this for annual exposure. The monthly EOG products and annual VNL differ in composite construction; cloudfree-weighted annual transfer remains subsequent work. The global post-2016 annual \`-0.235\` pixel correction is retained unchanged and **not re-applied to monthly VIIRS**.

## EcoVis command

\`\`\`bash
git pull --ff-only
python src/calibrate_prc_2021_snow_albedo_gain.py --self-test
python src/calibrate_prc_2021_snow_albedo_gain.py
\`\`\`

New outputs:
\`\`\`text
data/processed/historical_prc_2021_surface_gain_calibration/
  surface_gain_calibration_coefficients.csv
  surface_gain_heldout_groups.csv
  surface_gain_holdout_diagnostics.csv
  surface_gain_region_month_diagnostics.csv
  surface_gain_fitted_curves.pdf
  surface_gain_summary.txt
\`\`\`

Do not overwrite the 2021 EOG raster, old regional corrections or the historical PRC atlas until coefficient uncertainty, controls and transfer to the annual VIIRS composite are established.
