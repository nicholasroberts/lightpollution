# LED-dominated contrasts: same-format annual and municipal monthly VIIRS plots

**Prepared 10 November 2026** for `publication-figures`. This is an **a priori comparison**, not a selection of cases based on their satellite outcomes. A new script has been committed, but **has not yet run on the real EcoVis Roma/Richmond data**, so no VIIRS results are reported here.

## Objective and comparison to Coventry/Échirolles

The switch-off cases show abrupt municipal VNL V1 declines temporally aligned with independently documented changes in operating hours. We now test what the **same annual NUTS2016 VNL V2 + municipal monthly VNL V1** plots look like in two **large LED-dominated replacement programmes**. Crucially the new cases are *not verified undimmed*, and neither is claimed to retain the same total luminous flux for all years.

### Selection (independent historical evidence)

| Case | Fixed NUTS2016 regional comparator | Municipal mask | LED programme window / sources | Focus comparison and interpretation |
|---|---|---|---|---|
| **Rome** | **ITI43, Roma metropolitan region** (much larger than Comune di Roma) | **Union of all Rome municipal subdistrict polygons**, from public ArcGIS river authority layer; validated 900–1700 km² before use | Approx **2015–2019** conversion, ~51k fixtures by Nov 2016, ~170k by Dec 2018, >180k by Dec 2019 | **2018 vs 2017 same calendar months**: independent ACoS data indicate *installed* luminous flux **+1.0%** in this year while energy use **−27.4%**, but night-time dimming unknown. LED programme already under way in 2017. |
| **Richmond upon Thames** | **UKI75, Hounslow and Richmond upon Thames** (contains non-Richmond areas) | ONS December 2024 LAD **E09000027** only | ~11,000 LED replacements financed **Dec 2019**, staged in practice, much progress by 2022; **Dec 2019 is financing, NOT a verified installation start** | **2022 vs 2017/18 same months** as exploratory mixed-rollout comparison. Both earlier years may already contain replacement lamps. Funding notes unchanged operating hours, **but existing dimming settings remain unverified**. |

**These annual and monthly comparisons are distinct products and geographies.** Never subtract or divide annual NUTS3 mean by monthly municipality mean, or read a regional decade decline as the effect of the municipal retrofit.

### Input source: European tile 75N060W — correction, 10 November 2026

**Correction to the original version of this protocol:** EOG's six-tile map labels the Europe-containing tile **Tile 2 (60W75N)**, which is encoded in the files as **`75N060W`**. This tile spans longitudes roughly **60°W to 60°E** and includes **Rome (~12.5°E), Richmond and all five earlier study cities**. **Tile 3 (60E75N), encoded `75N060E`, covers East Asia/Pacific, not Rome.** The initial script/documents wrongly requested Tile 3 for Rome; this error has now been corrected. The user supplied the official EOG tile diagram to confirm.

**No extra EOG data must be downloaded for the current Rome/Richmond comparison.** The pre-existing 2015–2024 extracted `VCMSLCFG` `avg_rade9h` and `cf_cvg` from `75N060W` should cover **Rome 2016–2019 (48 observation months)** and **Richmond 2017–2024 (96 months)**. The original code's search now uses the same tile for both. The `--inventory` mode checks that every actual TIF pair exists before attempting any extraction or plot; if it reports missing months, investigate source files rather than acquiring `75N060E`.

**Remaining separate validation:** the script needs an authoritative **Rome municipality** boundary, which is not the NUTS3 `ITI43` metropolitan boundary. If the official municipal boundary query fails, the script must halt rather than substitute metropolitan Rome. Richmond is clipped to ONS LAD E09000027, not UKI75.

## Execution

<pre>
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures

# Real synthetic rendering test; no actual raw data accessed
python src/plot_led_dominated_contrast_figures.py --self-test

# Check EOG V1 monthly pair availability in the shared European tile
python src/plot_led_dominated_contrast_figures.py --inventory

# Generate the two FULL locked annual NUTS3 histories NOW
python src/plot_led_dominated_contrast_figures.py --annual-only

# Generate the Richmond borough and Rome municipality monthly plots
python src/plot_led_dominated_contrast_figures.py --cases richmond rome
</pre>

To produce both full cases directly from the already extracted European tile:

```bash
python src/plot_led_dominated_contrast_figures.py --cases richmond rome
```

If an error occurs after `--inventory` reports complete raster pairs, the next thing to check is **Rome's municipal polygon lookup**, not whether the imagery was downloaded from an eastern tile. Both cases are coded to use exactly `75N060W`.

## Script and outputs

[`src/plot_led_dominated_contrast_figures.py`](../../../src/plot_led_dominated_contrast_figures.py)

- Read-only locked annual source: `data/processed/fixed_nuts2016_viirs_2014_2024.parquet`, original final PRC endpoint CSV, **audited existing 2021 snow cache** from Figure 1. Validation against the frozen 2014/2024 radiance endpoints is compulsory; no rebuild or recalibration.
- Read-only monthly EOG V1 `VCMSLCFG` `avg_rade9h` and `cf_cvg`, observation-start month, fixed official municipality geometry.
- Annual panel: **2014–2024** NUTS2016 regional. Monthly panel: **2016–2019 Rome** or **2017–2024 Richmond**. Focus anomaly panel: **2018 Rome vs 2017 same-calendar months** or **2022 Richmond vs 2017/18**.
- Output figures (PNG/PDF/editable SVG): `figures/local_lighting_interventions/led_dominated_contrasts/`, named `rome_led_contrast_annual_and_municipal_monthly_2014_2024` and `richmond_led_contrast_annual_and_municipal_monthly_2014_2024`. When monthly inputs absent: `*_led_contrast_nuts3_annual_only_2014_2024`.
- Output CSVs: `data/processed/local_lighting_interventions/led_dominated_contrasts/`: frozen annual regional 22-row extract, each municipality monthly table, and calendar-month anomaly table. Figures are for hypothesis testing, not to infer a causal effect from one year.
- Output explicitly records `cf_cvg` median and `cf_cvg>=3/5` means for planned sensitivity checks.
- A **`--strict`** flag can require all requested monthly inputs instead of skipping missing cases.

## Historical documentary sources

- [Rome official/public-services ACoS LED programme](https://pubblicazioni.agenzia.roma.it/schede-33-piano_led_investimenti_e_obiettivi), and [ACoS 2019 report on installed flux and electricity 2013–2018](https://www.agenzia.roma.it/documenti/schede/relazione_annuale_2019_cap_5_servizi_a_rete.pdf). **2017→2018** total installed streetlighting flux rose from **1,991,000 to 2,010,000 klm (+0.95%)**, with estimated electricity use **115,640 to 84,000 MWh (−27.36%)**. However **2014→2018** flux fell from **3,377,000 to 2,010,000 klm (−40.5%)**, so **Rome's entire conversion was NOT equivalent to unchanged light emission**.
- [2018 ACEA report, >170k fixtures installed](https://sostenibilita2018.gruppo.acea.it/en/relations-stakeholders/customers-and-community/quality-delivered).
- [Richmond 2019 funder description, design brightness and operating hours](https://www.amberinfrastructure.com/sectors/case-studies/london-borough-of-richmond-meef/). The funder says *hours* unchanged, not that late-night dimming is absent.
- [Richmond council's description of LED/CMS operation](https://www.richmond.gov.uk/services/roads_and_transport/street_lights/replacement_street_lights/street_lighting_replacement_programme).
- [ONS Richmond borough E09000027](https://www.ons.gov.uk/explore-local-statistics/areas/E09000027-richmond-upon-thames).
- [Rome municipal subdistrict boundary layer](https://webgis.abdac.it/server/rest/services/T05_P5_Cf/MapServer/1) (query with geographic projection, union all official layer subdistricts, sanity-check geographic area).
- [EOG official monthly VNL V1 product documentation, including tile scheme](https://eogdata.mines.edu/products/vnl/).

## Review gate

The two cities are meaningful **prospective LED-dominated contrasts** regardless of whether VIIRS rises, declines or stays flat. Before describing either as a *true LED-only/no-dimming control*, the municipal lighting department's **actual historical controller set points at ~01:30** need confirmation. For Rome, compare the 2017/18 VIIRS results with the independently documented **installed luminous-flux history** and resist attributing a whole-city nocturnal change to installed flux alone. For Richmond, check conversion commissioning history and whether at-overpass output was controlled. Then retain the same four-point attribution logic as for Coventry/Échirolles.

## Rome municipal-boundary runtime fix — 10 November 2026

After both Rome and Richmond passed the 75N060W raster inventory (**48/48** Rome monthly rad+cf pairs; **96/96** Richmond), the real Rome case aborted at `official_boundary()` with `ValueError: Rome boundary query did not return GeoJSON features`. The old river-authority `T05_P5_Cf/MapServer/1/query` endpoint did not return the assumed 15-feature GeoJSON: it was an **unverified municipality source**, not a raster error. The locked annual VNL V2 history and audited 2021 snow cache had already loaded successfully. Richmond's earlier real plot also completed.

**Correction:** The plotting script now obtains the single **Comune di Roma** polygon from the **ISTAT-derived [geojson-italy Rome-province `P_58` municipal GeoJSON](https://github.com/guglielmo/geojson-italy/blob/main/geojson/limits_P_58_municipalities.geojson)**, accessible at `https://raw.githubusercontent.com/guglielmo/geojson-italy/main/geojson/limits_P_58_municipalities.geojson`. It selects the *exact* feature with both **ISTAT code `058091`** and **cadastral code `H501`**, rejecting a missing or duplicate selection, and validates polygon topology, EPSG:4326 bounds, inclusion of Rome city centre and **municipal area 1,100–1,500 km²** before caching. This is a derived re-publication of ISTAT data under CC BY, not a polygon produced by the VIIRS output. `ITI43` remains the separate **metropolitan NUTS3 annual** geography and is never used as the local municipal mask. The generated vetted municipality is cached in `data/raw/viirs/boundaries/rome_capitale_istat_058091.geojson` or the source boundary folder defined in `local_viirs_five_case_analysis.BOUNDARIES`; the script prints its exact path when validated. The deprecated `rome_capitale_union_15_municipi.geojson` is not read by the new code. No downloads of new satellite archives or adjustments to annual data are required.

**New first-stage check:**

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures

python src/plot_led_dominated_contrast_figures.py --self-test
# Only fetch/validate municipal polygon: no rasters, no annual sources
python src/plot_led_dominated_contrast_figures.py --check-boundaries --cases rome
# Full real 2016–2019 Roma municipal series and plot
python src/plot_led_dominated_contrast_figures.py --cases rome
```

The self-test now verifies **ISTAT feature identity** and rejects non-GeoJSON/error/duplicate responses, alongside the existing synthetic figure rendering; the `--check-boundaries` check isolates upstream GeoJSON/network problems before the monthly processing. The new boundary has **not yet been fetched/verified on EcoVis** in this chat; only the endpoint's documented file/schema and code changes have been inspected. If the primary raw GitHub download fails, do not substitute metropolitan geometry or guess a URL; report the HTTP/error first.
