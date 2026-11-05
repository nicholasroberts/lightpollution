# Annual 2014–2024 spider polarization receptor contrast maps for the talk

**Branch:** \`publication-figures\`. **5 November 2026**. Code and instructions committed; full EcoVis execution and visual QA pending.

## Scope

Plot **11 separate annual maps of spider PRC level**, for each year 2014–2024 inclusive. These are **not** delta-PRC maps and not percentage change maps. Source is the **already LOCKED pixel-derived spider PRC regional annual table**, the identical scientific dataset used by publication Figure 3C. There is no new radiance or PRC calibration, native-pixel extraction, geometry re-aggregation, or interpolation.

The design **mirrors** \`src/plot_publication_annual_nuts3_radiance_maps.py\`: fixed NUTS2016 geometry (1,522 NUTS3 regions), the same continental window (-12°,34°) to (46°,73°) including Türkiye, EPSG:3035 map projection, up to 150 m display-only geometry simplification, one uniform colour bar for **every year**, white 16:9 canvas 12.8 × 7.2 inches with large year and title and **no footer**, country outlines, and Illustrator-compatible SVG/PDF plus talk PNG (240 dpi). There is one combined **11-page 16:9 PDF**.

### Source and 2021 provenance

Required input, read only:

\`\`\`text
data/processed/historical_prc_2014_2024_spider_snow_corrected/prc_nuts3_annual_cumulative.csv
\`\`\`

This is the existing 1,522 × 11 annual NUTS2016 spider PRC product previously created by \`publish_prc_locked_snow.py\`. The \`prc_mean\` in every row is a regional mean of PRC originally **transformed at native VIIRS pixel resolution**, with the approved post-2017 calibration and **locked 2021 native-pixel snow correction already incorporated**. Do **not** redo the 2021 patch or apply the 0.235 pixel offset to the regional PRC.

The new script imports the locked-input validator from Figure 3, requiring all 16,742 rows, every year in 2014–2024, one observation per region/year, finite PRC values in 0–1, positive pixel coverage area and \`snow_status=locked\`. It verifies unchanged region IDs every year and checks each 2014 and 2024 stored cumulative difference against the annual PRC endpoints. When the archived paired radiance/PRC endpoint CSV is present, it **independently compares spider PRC 2014 and 2024 against that source** and stops on material disagreement. Missing country polygons stay grey; neither area weights nor missing-value filling are introduced in the display.

### Scale choice: yearly absolute level, not change

Unlike the radiance atlas, the meaningful spider PRC values occupy a narrow subset of their mathematically possible 0–1 range. Therefore, the shared scale is **linear in original, dimensionless PRC** and is fixed from the **pooled continental minimum and maximum of ALL 11 years** rounded outward to 0.01 PRC; no actual data are clipped. This should reveal variation in the biologically relevant PRC range without exaggerating year-to-year shifts through annual rescaling. Dark purple = lower predicted PRC; brighter yellow = higher PRC. The colour map is \`viridis\`, distinct from the radiance \`magma\` atlas. The exact pooled range and fixed endpoints will be reported **after EcoVis execution** in the provenance JSON; they are not guessed in advance.

If an alternative global scale is desired, use \`--vmin 0.14 --vmax 0.47\` **only if these supplied limits fully cover every continental PRC value**; the script rejects limits that clip the observed data. No per-year colour rescaling is permitted.

The maps are observations of modelled mean receptor contrast, *not* directly measured insect orientation or a signal-to-noise/behavioural threshold; the full range of possible impact of a 0.05 PRC change remains a separate scientific analysis.

## Run on EcoVis

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/plot_publication_annual_nuts3_spider_prc_maps.py --self-test
python src/plot_publication_annual_nuts3_spider_prc_maps.py
\`\`\`

An optional \`--png-dpi 300\` increases PNG resolution. \`--png-dpi 0\` disables PNGs while retaining SVG/PDF and the PDF deck. No global TIFF/raster processing or web downloads are required: this reads only the existing locked annual CSV and existing fixed NUTS2016 geometry.

## File output — independent publication folder

\`\`\`text
figures/publication/annual_nuts3_spider_prc_2014_2024/
    spider_prc_nuts3_2014.svg
    spider_prc_nuts3_2014.pdf
    spider_prc_nuts3_2014.png
    ...
    spider_prc_nuts3_2024.svg
    spider_prc_nuts3_2024.pdf
    spider_prc_nuts3_2024.png
    all_years_nuts3_spider_prc_2014_2024.pdf
    annual_nuts3_spider_prc_2014_2024_map_values.csv
    annual_spider_prc_map_provenance.json
\`\`\`

There are **33 individually exported images + a multipage PDF + two provenance/value files**. All generated figure outputs are ignored by Git, like the radiance atlas; script and instructions are committed.

## QA on the new atlas

1. Confirm \`--self-test\` passes, then full run reports **11 × 1,522** scientific annual values, each \`snow_status=locked\`.
2. Check **the 2021 map** reads the existing snow-corrected biological PRC (does not substitute the uncorrected radiance or redo its native-pixel correction).
3. Verify the colour scale and physical PRC ticks are **identical on all 11 maps** and encompass the pooled source extrema without clipping; do not compare map hue differences as if the colourbar showed ΔPRC.
4. View representative 2014, 2019, 2021, 2024 PNG/PDF/SVG side by side and check continental extent, eastern Türkiye, year labels and absence of footers.
5. Read \`annual_spider_prc_map_provenance.json\` for country geography, endpoint QA, mapped region counts, full-range scale and original source metadata. The CSV remains a map-ready **copy** of existing scientific values rather than a newly estimated response.

The existing annual **radiance** maps and Figure 3A–F plots are not modified. A PRC timelapse can be produced subsequently from these already-rendered consistent-scale PNGs, as in the radiance workflow, if explicitly wanted.
