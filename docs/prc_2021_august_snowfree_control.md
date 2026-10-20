# August as the snow-free negative control for the 2021 radiance anomaly

**Decision (19 October 2026):** the investigator rejected March as a snow-free control because substantial European snow may persist in March. Use **August** instead, comparing the same calendar month in 2020, 2021 and (when satellite-compatible) 2022. Classify snow-free areas **from observed MOD10CM in each August**, not a seasonal assumption. Retain January/February as winter calibration months.

**Data availability caveat:** the original 2021 visible-albedo AppEEARS downloads and \`vnp43ma3_002_*_monthly.tif\` cache were specifically assembled for January–March. An August combined snow+albedo placebo cannot be completed unless suitable QA-controlled **August 2020, 2021, 2022 VNP43MA3 visible WSA albedo** exists or is acquired and processed. The new August stage intentionally makes NO assumption that August albedo changes are zero.

**Satellite caveat:** an earlier project inventory flagged the **August 2022 VCMSLCFG product as J01**. The new script inspects the actual \`npp\`/\`j01\` token in each VIIRS file. It uses August 2020 and 2021 only if their satellite platform matches. It allows the additional 2020/2022 flank comparison **only if August 2022 also matches the 2020/2021 platform and grid**. Otherwise it reports August 2021 vs August 2020 and explicitly excludes the 2022 radiance layer from matched temporal contrasts; MOD10CM snow status is still checked across 2020/21/22. Do NOT treat NPP-versus-J01 differences as lighting trends or as a snow-free control result.

## Initial August control (ready on EcoVis)

Stage: \`src/prc_2021_august_snowfree_control.py\`. Uses the same 18 pilot regions previously selected in \`data/processed/historical_prc_2021_pixel_snow_test/selected_nuts3_pilot_regions.csv\`. It reads existing EOG monthly VCMSLCFG August radiance+cf_cvg and original MOD10CM 061 August snow, uses region-windowed reprojection and requires both VIIRS monthly products to have positive cloud-free support at identical pixel coordinates. It labels pixels \`snowfree_all_three\` only if MOD10CM August snow cover is no more than **1% in each year 2020, 2021, 2022**. A sensitivity on 0%/2% can be run; no cells are presumed snow-free.

It retains the native MOD10CM 0.05° snow-cell identifier, and reports 2021–2020 August radiance changes on the matched subset for snow-free and non-snow-free strata, with brightness bins. The resulting August changes are background variability/processing diagnostics, not fitted corrections. **No VNP43MA3 albedo correction can be tested without August albedo QA observations**. Test with:

\`\`\`bash
git pull --ff-only
python src/prc_2021_august_snowfree_control.py --self-test
python src/prc_2021_august_snowfree_control.py --limit 6
\`\`\`

If input source platforms and snow fractions are satisfactory:

\`\`\`bash
python src/prc_2021_august_snowfree_control.py
\`\`\`

Outputs (old winter fit and historical PRC unchanged) under \`data/processed/historical_prc_2021_august_control/\`:

- \`august_source_platforms.json\`: actual file paths, satellite provenance and whether 2022 is usable
- \`august_snowfree_valid_pixel_coverage.csv\`: valid sample area and measured snow-free fraction per NUTS area
- \`august_native_snowcell_brightness_groups.csv\`: native snow-cell × VIIRS brightness group outcomes
- \`august_snowfree_region_comparison.csv\`: region radiance 2021/2020 (and 2022 if valid), stratified by true snow absence
- \`august_snowfree_control_summary.txt\`: QA, methodological limitations, initial findings when run

## Next scientifically meaningful combined-model August test

Acquire/build August **VNP43MA3.002** QA0/1 WSA visible albedo composites for all three years from the same provenance and QA method as the winter series, then evaluate the frozen January/February-fitted snow+albedo gain \`G\` on August matched snow-free cells. Model should infer little/zero *snow* amplification in snow-free locations. However, if positive summertime albedo deviations occur, the **albedo component may still predict a gain**: test explicitly for false-positive VIIRS attenuation, and do not zero out albedo or fit summer coefficients merely to obtain G=1. Contrast albedo-stable and albedo-changing summer groups and include placebo winter snow-free cells. Evaluate August observational prediction with an independent summer/background nuisance baseline; the winter \`T\` parameter is not expected to transfer directly to summer.

August does not automatically control for vegetation, clouds, atmosphere or real lighting changes. It is **an improved test of snow absence**, not an independent physical causal identification of the full calibration. Prior March results remain reported accurately but are **reclassified as an additional late-winter/early-spring sensitivity**, not the primary snow-free negative control.

**No historical annual 2014–2024 VIIRS or PRC values have been modified.**


### Implementation update — 20 October 2026

The August control **source script is now actually committed** to this branch; the earlier version of this document described an intended script before its source existed in the repository. The script uses prior 18 pilot NUTS3 IDs, existing MOD10CM HDF/cached GeoTIFF August products and available EOG August monthly VCMSLCFG rasters. It groups radiance changes at native coarse MOD10CM-cell × VIIRS brightness-category level, requires 2020/21 common unmasked cloud-free VIIRS pixels, checks observed August snow fractions in all three years, and excludes 2022 VIIRS radiance **if sensor platform differs**. This is not a full snow+albedo independent placebo until August visible WSA albedo QA composites are obtained.

The May–October 2026 scientific decision is NOT to force the 2021 anomaly to match 2020/22. Our January–February snow/albedo empirical gain is a candidate; use August to falsify residual year-wide VIIRS shifts under truly snow-free conditions before deciding whether and how to map winter gains to the *annual* 2021 EOG VNL pixel composite. Do not overwrite uncorrected PRC baselines.
