# Annual NUTS3 radiance atlas, 2014–2024 — publication / talk

**Prepared on `publication-figures`, 3 November 2026. Status: script committed, awaiting EcoVis execution and visual QA.**

## Objective

Create **eleven individual annual radiance maps**, 2014 to 2024 inclusive, with identical map geometry, physical radiance colour scale, legend bounds, title placement and slide layout. The maps display **absolute annual calibrated night-time VIIRS DNB radiance at fixed NUTS2016 NUTS3 level**, not each year's change from 2014 or percentage change.

All maps are rendered as individually editable **vector SVG and PDF** plus a 16:9 high-resolution **PNG** for the talk; one combined **11-page PDF** is also produced. This is a presentation output only. There are **no new light-pollution calibrations or statistical models**.

## Correct LOCKED source values

Uses `src/plot_publication_annual_nuts3_radiance_maps.py`, reading:

- `data/processed/fixed_nuts2016_viirs_2014_2024.parquet` — **calibrated annual NUTS3 mean radiance**, field `radiance_mean_corrected`. This table applies the originally approved 2017-onward positive-native-pixel calibration (−0.235 nW cm⁻² sr⁻¹) and already contains fixed NUTS2016 regional means. Use values **as stored** in 2014–2020 and 2022–2024.
- `data/processed/publication_figures/fig1_2021_locked_snow_radiance.csv` — existing *validated, provenance-checked* 2021 native-pixel snow overlay extraction (and adjoining JSON audit). This cache is created via `plot_publication_figure1_radiance.get_locked_snow_2021` on first use or when source/provenance changes; the source is `config/viirs_2021_snow_correction.yaml` with **status=locked**. The 2021 cache holds both unpatched calibrated NUTS3 means and **locked snow-corrected means**. Original annual parquet 2021 values are checked before replacement. The original raw annual TIFF and locked native patch are not changed. The 2017+ pixel offset is applied **only once, after the raw patch**, in this existing extraction process; the script does not subtract it at NUTS3 level again.
- `figures/prc_radiance_spider_moth_2014_2024/nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv` — frozen 1,522-region **2014 and 2024** endpoints, used for independent cross-checks (the same 99th percentile absolute tolerance used in Figure 1). The earlier annual values are not silently overwritten by endpoint values.
- `data/raw/nuts/NUTS_RG_01M_2016_4326_LEVL_3.geojson` via configured NUTS2016 loader — **1,522** geographic polygons before the display crop. The existing source parquet might have fewer than 1,522 valid observations in a year; those are explicitly **grey and counted** rather than given invented zeros. Actual zero radiance is a valid dark map value.

**Do not substitute** the old 2021 uncorrected VIIRS annual parquet value or a merely candidate 2021 surface-correction model. The Figure 1 locked native-pixel patch is the approved correction.

## Cross-year visual comparability

- **Identical geographic extent for all eleven maps:** longitude −12° to 46°E, latitude 34° to 73°N (continental presentation, deliberately including eastern Türkiye; excludes far-flung overseas NUTS regions from the *display* only).
- Projected **EPSG:3035** and display-only **150 m polygon simplification** to reduce Illustrator path complexity. All actual radiance statistics remain based on the established **unmodified NUTS3 geometry**.
- **Single shared colour scale across all years**: `magma`, dark=low, yellow=bright. Sequential map from **0** to a fixed *rounded-up maximum of all mapped values across 2014–2024*; no annual rescaling and **no clipping** of recorded continental radiance. Physical radiance in **nW cm⁻² sr⁻¹** on the legend.
- **Symmetric-log mathematical scale restricted to non-negative radiance**, with 0–0.2 nW linear and logarithmic above that. This expands detail in dark regions while retaining urban bright areas. It is a colour transformation only and **does not change the data**.
- 16:9 white slide canvas at 12.8 × 7.2 inches; large year labels and the same colourbar position, geographic extent and legend across time. PNG default 240 dpi; slide PNGs and combined PDF preserve exact 16:9 page (vector SVG/PDF individual exports may trim whitespace).
- Each map makes the correction status visible: 2014–2016 original calibrated annual values; 2017–2024 positive-pixel 0.235 correction; 2021 **also** has locked snow correction.

## Run on EcoVis

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/plot_publication_annual_nuts3_radiance_maps.py --self-test
python src/plot_publication_annual_nuts3_radiance_maps.py
```

The first run will inspect and reuse `fig1_2021_locked_snow_radiance.csv` if its JSON audit still matches the original 2021 TIFF, locked patch, model and pixel offset. If it doesn't exist, it must run the **approved native-pixel extraction** using the original 2021 global TIFF; this can take significantly longer than the maps alone. The cache is **not an alternative source** if its provenance fails validation.

Change PNG resolution without modifying science: `--png-dpi 300`. Disable PNG copies (SVG/PDF only): `--png-dpi 0`.

Explicit `--vmax 300` is allowed only if **none** of the continental radiance records in 2014–2024 exceeds 300; otherwise the command fails to avoid silently clipping very bright areas. Omitting `--vmax` automatically determines the single shared upper bound.

## File output

Directory (already inside the main `figures` tree for rsync to the MacBook):

`figures/publication/annual_nuts3_radiance_2014_2024/`

For **each** 2014, 2015, …, 2024:

```
radiance_nuts3_2014.svg
radiance_nuts3_2014.pdf
radiance_nuts3_2014.png
...
radiance_nuts3_2024.svg
radiance_nuts3_2024.pdf
radiance_nuts3_2024.png
```

Plus:

- `all_years_nuts3_radiance_2014_2024.pdf` — eleven-page PDF deck, each page 16:9.
- `annual_nuts3_radiance_2014_2024_map_values.csv` — original corrected value and exact year / NUTS3 identifiers for each included region, including explicit `native_pixel_2021_snow_locked` and per-value source description.
- `annual_radiance_map_provenance.json` — original file paths, locked 2021 cache audit, endpoint mismatch checks, map region coverage for each year and single global colour scale.

Total: **11 maps × 3 formats = 33 individual maps, plus one PDF deck and two QA files**.

## Production QA

Confirm the 2014/2024 annual values pass source endpoint comparisons; 2021 reports the snow-corrected NUTS3 substitution; all eleven years are present; the legend maximum is exactly the same; and the CSV flags 2021. Check PNGs in presentation software and SVG vector text in Illustrator. Source statistics and calibration parameters must remain unchanged. The script prints all output filenames with sizes as the final terminal block. Outputs are ignored by Git and so are generated locally on EcoVis.

To copy to the home Mac (use the same real SSH hostname or configured alias used for login), run from a **local Mac terminal** with university VPN connected:

```bash
rsync -avh --progress \
  bznwr@ECOVIS_HOST:/run/media/bznwr/data/lightpollution/figures/publication/annual_nuts3_radiance_2014_2024/ \
  ~/Desktop/annual_nuts3_radiance_2014_2024/
```


## Talk timelapse and clean frame revision — 3 November 2026

At the investigator's request, **removed both footer lines from all annual maps**:

- `Annual NUTS3 radiance · fixed NUTS2016 geography · common nonlinear colour scale`
- `Annual VIIRS, fixed NUTS2016 regional means` (and the year-specific 2017+/2021 substitute footer captions)

The annual map generator now shows **only the year/title, geographic radiance map and shared radiance colourbar**. Any conditional bottom grey/missing-region legend was also removed for a clean presentation frame; missing NUTS3 polygons remain grey in the visual and are explicitly recorded in `annual_radiance_map_provenance.json`. Footer removal and slightly expanded bottom map margin are **layout-only**. The same source correction details remain in CSV/JSON provenance, not burned into frames. Full 2014–2024 maps are **re-rendered** locally by running the atlas generator again; otherwise existing PNGs retain old text.

Added `src/make_annual_nuts3_radiance_timelapse.py`, which uses FFmpeg to encode the existing **eleven rendered PNGs** (`radiance_nuts3_2014.png` through `radiance_nuts3_2024.png`) at **1 second per year**, total **11 seconds**. Video is **1920×1080**, 30 fps H.264/yuv420p for cross-platform slide playback, 330 frames in total, no interpolation between calendar years. Each input still is held unchanged for one second; all colour/extent/scale comparisons remain as produced by the annual atlas. The original PNGs, SVGs and PDFs are not modified by video creation.

### EcoVis commands (run in this exact order)

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/plot_publication_annual_nuts3_radiance_maps.py --self-test
python src/plot_publication_annual_nuts3_radiance_maps.py

python src/make_annual_nuts3_radiance_timelapse.py --self-test
python src/make_annual_nuts3_radiance_timelapse.py
```

The first full command **replaces the previously rendered annual PNGs/SVGs/PDFs** with cleaned, footer-free versions. The video command then reads those new PNGs. `ffmpeg` and, for duration/frame verification, `ffprobe` must be installed on EcoVis; check `ffmpeg -version` or install with `sudo apt install ffmpeg` as appropriate.

### New outputs (same atlas folder)

```
figures/publication/annual_nuts3_radiance_2014_2024/
  radiance_nuts3_2014.png             # clean footer-free map
  ...
  radiance_nuts3_2024.png
  annual_nuts3_radiance_timelapse_2014_2024.mp4
  annual_nuts3_radiance_timelapse_2014_2024.json
```

The JSON reports all PNG source dimensions, 11 × 1-second hold, frame count, 1080p encode and the actual `ffprobe` duration when available. The MP4 must be visually reviewed before presentation (year labels, colourbar size and absence of footer text); successful encoding cannot be claimed until EcoVis runs it.

### Transfer ONLY the MP4 to the home Mac

From **local Mac Terminal** with the university VPN connected:

```bash
rsync -avh --progress \
  bznwr@137.222.175.129:/run/media/bznwr/data/lightpollution/figures/publication/annual_nuts3_radiance_2014_2024/annual_nuts3_radiance_timelapse_2014_2024.mp4 \
  ~/Desktop/annual_nuts3_radiance_2014_2024/
```

Create the destination folder first if absent, using `mkdir -p ~/Desktop/annual_nuts3_radiance_2014_2024/`.
