# Historical 2014–2024 pixel-derived PRC mapping — initial production stage

**Purpose:** calculate the provisional biological PRC function on each valid pixel of the local annual EOG VIIRS \`average_masked\` raster before NUTS2016 area-weighted polygon averaging. The resulting yearly table is **not** \`f(mean NUTS3 radiance)\`. Rank NUTS3 areas by mean pixel PRC change and report European-area trends.

Working equation from \`config/prc_transfer.yaml\`:

\`\`\`text
PRC(R) = 0.144845 + 0.316330 * exp(-0.106483 * R)
\`\`\`

where \`R\` is calibrated EOG radiance in nW cm-2 sr-1. As in prior pipeline: 2017+ positive pixels have a 0.235 radiance unit subtraction, clipped to zero. Zero pixels remain valid; negative or nonfinite pixels are masked as physically invalid for this biological function. 2014–2016 pixels are unadjusted.

The geographic anchor is fixed NUTS2016 level 3 for all years, including 2014 and 2024. Therefore rankings are not affected by changing NUTS releases. Geographic *coverage* is what NUTS2016 covers, not every European country. The map display window excludes non-mainland overseas territories, but the tables include the full NUTS2016 inventory.

## Run on EcoVis

From the repo root on branch \`historical-polarization-visibility\`:

\`\`\`bash
git pull --ff-only
python src/historical_prc_pixel_nuts3.py --pilot
\`\`\`

Pilot years 2014, 2018, 2019, 2021, 2024, for testing source-product availability, raster source implementation and trends before committing to full processing. Pilot year-over-year changes only exist where adjacent sampled years are present; 2014-relative changes and endpoint rankings are available. Review the 2021 anomaly against earlier diagnostics.

Full series:

\`\`\`bash
python src/run_pipeline.py historical-prc
\`\`\`

Optional native-resolution Europe GeoTIFF for 2024, containing five bands: PRC 2024, 2023→2024 absolute, 2023→2024 %, 2014→2024 absolute, 2014→2024 %:

\`\`\`bash
python src/historical_prc_pixel_nuts3.py --no-maps --pixel-year 2024
\`\`\`

This optional command re-runs the NUTS3 extraction and can create a large compressed raster; choose appropriate disk space. Additional annual pixel-band exports can be requested with extra \`--pixel-year\` arguments. **Pixel-to-pixel differences** are only written if original datasets share a native grid exactly; otherwise explicit reprojection/registration is needed, never silently assumed.

## Files created (not yet run on the user's local datasets)

\`\`\`text
data/processed/historical_prc_2014_2024/
    prc_nuts3_annual_cumulative.csv
    prc_nuts3_annual_cumulative.parquet
    prc_europe_annual_cumulative.csv
    prc_nuts3_rankings_2014_2024.csv
    prc_pixels_2024_annual_and_2014_cumulative.tif   # optional
figures/historical_prc_2014_2024/
    annual_absolute_prc_nuts3_atlas.pdf
    annual_percent_prc_nuts3_atlas.pdf
    cumulative_absolute_prc_nuts3_atlas.pdf
    cumulative_percent_prc_nuts3_atlas.pdf
    prc_europe_annual_cumulative_changes.pdf
\`\`\`

## Calculation and interpretation

Each pixel is transformed via PRC(R) and the **spherical area-weighted mean of pixel PRC** is computed in a fixed NUTS3 polygon. A region's annual delta is the difference between two yearly mean PRCs, and cumulative delta compares to 2014. Percent changes use the immediately previous year / 2014 reference *PRC*, not a percent change in radiance.

For headline ranks, a region must be present every processed year and have all annual valid VIIRS area totals within approximately **±5% of the 2014 valid area** (threshold set by \`--min-coverage 0.95\`). This is a first data-quality guard, **not identical valid-pixel masking**. NUTS-area masks should be audited in anomalous regions and a strict common-pixel mask considered if annual coverage changes affect results. Area-weighted Europe summary includes all valid NUTS2016 areas, with the aggregate valid km2 provided annually to diagnose differences in coverage.

Absolute positive ΔPRC indicates an **increase** in predicted receptor contrast, not necessarily better navigational performance; no behavioural cue-detection threshold has been established. Percent change is relative to PRC, **not** radiance. PRC is upper-bounded by \`y0=0.461175\` for nonnegative R and approaches \`yf=0.144845\` at extreme brightness. It is inappropriate to substitute PRC of NUTS-mean radiance for the average pixel PRC.

## Outstanding 2021 and provenance checks

The earlier annual 2021 snow/albedo correction is calibrated on **NUTS3 region means**, not an individual-pixel field, so it must not be applied unchanged to this pixel PRC stage. Outputs from 2021, and consequently 2020→2021 and 2021→2022 annual differences, are explicitly **diagnostic/provisional** pending pixel-level snow-state assessment. Likewise do not infer causation for a region's brightening/dimming from the satellite series alone. EOG raster versions may vary (the site 2018/19 rasters were V2.1) and must be reported from actual filenames.

The fitting observations still contain two unidentified PRC datapoints and four provisional site assignments. Main model is the investigator's accepted working function; it is not yet a formally validated extrapolation to all taxa, spectral conditions, geographic regions or behavioural outcomes.
