# Five local lighting interventions: 2015–2024 EOG V1 monthly Stage B analysis

**Authoring state:** script and this design were committed to `publication-figures` on 9 November 2026. The investigator reports complete extraction of the newly downloaded EOG monthly archives; **the full five-case script has not yet run against EcoVis, so there are no new reported outcomes or inferred effects in this document**.

## Scientific question and scope

Can changes in VIIRS DNB monthly radiance at documented European street-lighting interventions be explained solely by LED spectral response, or are some changes more consistent with actual reductions/increases in **lighting operating hours or output**? We do not claim VIIRS radiance is a direct lumens/flux measurement and cannot uniquely partition LED spectral conversion, operation, atmospheric effects and surface reflectance from VNL V1 alone.

**Source series:** EOG monthly **V1** `VCMSLCFG`, tile **75N060W**, January 2015–December 2024: **120 × 2** monthly `avg_rade9h.tif` and `cf_cvg.tif`. All 84 previously missing month archives were located as of the latest EcoVis inventory; 2020–2022 was already available. Assumption: archive extraction has completed. The script independently requires **all 120 months** and stops on missing input; archive presence alone is not treated as raster validation.

**Source code:** [`src/local_viirs_five_case_analysis.py`](../../../src/local_viirs_five_case_analysis.py)

## Geography and independent treatment hypotheses

| City | Primary spatial unit | Locked NUTS2016 context, not a treatment polygon | Intervention type and date |
|---|---|---|---|
| Échirolles | Official French commune **INSEE 38151**, cached from `geo.api.gouv.fr` | FRK24, Isère | Named street/site switch-offs from July 2022, municipal road-hours switch-off 1 Oct 2022; most LED replacement later |
| Portsmouth | ONS December 2024 local authority **E06000044** | UKJ31 | Staged LED 2018–19; summer 2020 adaptive dimming |
| Coventry | ONS December 2024 local authority **E08000026** | UKG33 | Part-night operation initiated May–July 2024 |
| Bristol | ONS December 2024 local authority **E06000023** | UKK11 | Gradual LED/CMS programme from March 2023, ongoing 2024 |
| Thanet | ONS December 2024 district **E07000114** | UKJ44 **East Kent** (much larger than Thanet) | LED plus more all-night lighting Dec 2016–May 2017; further 2018 stages |

UK boundaries are retrieved **per authority** from ONS 2024 **Local Authority Districts (BGC, 20-m generalised, clipped to coastline)** ArcGIS FeatureServer and cached under `data/raw/local_lighting_boundaries/ons2024_lad_<code>.geojson`. ONS authoritative catalogue: https://www.data.gov.uk/dataset/02fd4563-a97b-47c3-9c06-8e7379d25b06/local-authority-districts-december-2024-boundaries-uk-bfc (full resolution BFC) and https://www.data.gov.uk/dataset/local-authority-districts-december-2024-boundaries-uk-bgc (generalised BGC). The ONS 2024 support is fixed retrospectively across all years: a fixed geographical comparison, not a claim that every street-light circuit was in operation throughout.

**Do not silently use a NUTS3 polygon as a municipality.** The script stops if an official boundary cannot be downloaded/validated. Only the explicit `--allow-nuts3-fallback` diagnostic option permits that proxy, and labels all outputs accordingly. **Never interpret fallback Thanet UKJ44 East Kent as an estimate for Thanet district**.

## Predeclared primary contrasts (avoid choosing years after seeing satellite values)

| Case | Primary outcome | Reference years, same calendar months | Why |
|---|---|---|---|
| **Échirolles** | **Oct–Dec 2022** | **2020 + 2021** Oct–Dec | Late-2022 operation change prior to major 2023 LED rollout; check 2020/21 both rather than 2021 alone |
| **Portsmouth** | **Sep–Dec 2020** | **2018 + 2019** Sep–Dec | Late-2020 adaptive dimming vs prior seasons; 2018–19 are **not pristine sodium-only controls** because LED retrofit was ongoing |
| **Coventry** | **Aug–Dec 2024** | **2022 + 2023** Aug–Dec | Exclude staged May–July 2024 changeover |
| **Bristol** | **Apr–Dec 2023** | **2020–2022** Apr–Dec | Citywide evolution **during**, not after, a rolling LED/CMS programme; LED and dimming remain confounded |
| **Thanet** | **Jun–Nov 2017** | **2015 + 2016** Jun–Nov | After first Dec 2016–May 2017 implementation; 2016 December is excluded as transition; compare sign for possible brightening |

Secondary descriptive follow-ups: Échirolles Oct–Dec 2023 (with LED changes), Portsmouth Sep–Dec 2021 (post/post relative to 2020), Bristol Apr–Dec 2024 (continued rollout), Thanet Jun–Nov 2018 (further stage not individually dated). None is treated as an independent treatment replication.

**Absolute `ΔR` is primary.** Summary percentage is `100 × (mean outcome months − mean reference-month baseline) / mean reference-month baseline`, not the mean of pixel percentages. For each outcome month, reference is the **arithmetic mean of municipality means for the same calendar month in the reference years**. This is a *descriptive seasonal contrast*; no model-inferred causal effects or formal p-values are computed.

## Input validation / cloud-free observations

The script:
1. Indexes EOG `SVDNB_npp|j01` filenames by the **observation start date** and rejects wrong `75N060E` tiles; a later `_cYYYYMMDD` creation date is never mistaken for the observation month.
2. Requires both the primary monthly radiance and `cf_cvg` GeoTIFF for each 2015–2024 month. Missing months cause an explicit error, not silent interpolation.
3. Checks paired radiance/coverage raster dimensions, CRS and affine transform, and fixed-grid alignment across all months in each municipality.
4. Uses pixel centres inside the **fixed municipality polygon** (`all_touched=False`), radiance not nodata, finite values and **`cf_cvg > 0`**. Does not apply the locked **annual V2 −0.235** zero correction to EOG monthly V1.
5. Reports 120 monthly means/medians, interquartile ranges, fraction of static pixel centres observed, median and quartiles of cloud-free observations, and radiance sensitivities requiring **`cf_cvg >=3` and `>=5`**.
6. For each primary spatial contrast, matches **the same native-grid cells** for each treatment month and each corresponding reference year; maps the mean of those month-specific paired radiance differences. A cell enters the aggregate map only when **at least 75%** of the specified months have valid paired observations. Reports count of supported months as a separate map. **Neighbouring VIIRS pixels are spatially correlated, not independent replicates**.
7. Flags **Jan–Mar 2021** as surface reflectance/snow-sensitive months; `cf_cvg` does **not** correct snow, aerosols, cloud optical properties or differences in clear-night timing. Snow in other years remains possible.
8. Keeps 2015–2024 separate from the locked European annual NUTS3/PRC time series. It does not change the historical paper's interpretation or compute a flux partition.

## Outputs

Under `data/processed/local_lighting_interventions/five_case_monthly/`:
- `<city>_2015_2024_monthly.csv` (120 rows for each case)
- `<city>_calendar_matched_contrasts.csv` (monthly outcome-year/reference differences)
- `five_case_monthly_2015_2024.csv` (all selected case rows)
- `five_case_contrast_summary.csv` (one row for each of primary/secondary contrasts; primary includes matched-pixel map median)
- `<city>_<outcome-year>_versus_<reference-years>_delta.tif` (native spatial ΔR primary)
- `<city>_<outcome-year>_versus_<reference-years>_matched_month_count.tif` (number of same-calendar-month paired observations supporting the map)

Under `figures/local_lighting_interventions/five_case_monthly/`, both PNG and PDF:
- `<city>_2015_2024_series`: raw decade trajectory, valid-pixel fraction, cloud-free observations; intervention shaded
- `<city>_seasonal_comparison`: same-month trajectories for selected years, with cf-cvg reference
- `<city>_<year>_primary_contrast` and relevant `secondary_contrast`: year/reference calendar-month series, ΔR, cf-cvg
- `<city>_<year>_versus_<references>_period_maps`: matched baseline, outcome, difference with shared absolute colour scale and zero-centred difference
- `<city>_<year>_versus_<references>_monthly_maps`: per-calendar-month paired differences
- `<city>_<year>_versus_<references>_matched_coverage_map`: supported matched-pixel months
- `<city>_<year>_versus_<references>_cloudfree_observation_maps`: maps of *mean accepted cloud-free observation counts per pixel per month*, the outcome-minus-reference difference and their shared absolute scale; separate from the binary matched-month validity count

## Exact EcoVis command sequence

```bash
cd /run/media/bznwr/data/lightpollution
git status --short
git switch publication-figures
git pull --ff-only origin publication-figures

# 1) Independent code and input checks (no full GIS work)
python src/local_viirs_five_case_analysis.py --self-test
python src/inventory_monthly_viirs_acquisitions.py
python src/local_viirs_five_case_analysis.py --inventory

# 2) Run the previously analysed two cases first
python src/local_viirs_five_case_analysis.py --cases echirolles portsmouth

# 3) Run all five (reruns those two; replaces only this stage's outputs)
python src/local_viirs_five_case_analysis.py
```

If the boundary web service fails, report the error and use/obtain the **official** cached GeoJSON rather than silently switching to NUTS3. Existing EcoVis `.venv` already had the required geospatial dependencies for the pilot. The run may access the ONS and French polygon services but makes no EOG downloads; large raster files stay local.

**Review gate:** First inspect `five_case_contrast_summary.csv`, five `*_seasonal_comparison.png` and the primary `*_period_maps.png`. Evaluate pre-trends and month-specific cloud-free coverage before causal language. Next stage: consider predeclared nearby atmospheric-context polygons (distinct from proven untreated controls), local weather/snow, 2018/19 LED rollout dates and the documented street/column footprints if needed.
