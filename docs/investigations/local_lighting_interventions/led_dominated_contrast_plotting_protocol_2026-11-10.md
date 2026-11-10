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

### Input source constraints, particularly the EAST tile

The 2015–2024 EOG V1 `VCMSLCFG` monthly TIFs used for Coventry, Échirolles and Richmond are the **western `75N060W`** tile. **Roma at longitude ~12.5° E is in `75N060E`**. The existing western raster coverage does NOT contain Rome.

To avoid an unnecessary 10-year full eastern-tile archive acquisition, the initial Roman monthly plot is **limited to 2016–2019**, which is **48 EOG `75N060E` month pairs** (`avg_rade9h`, `cf_cvg`). That's still a large download and requires authorised EOG archive access. No raw TIFFs are fabricated; the script **reports exactly which months are missing** and can still draw the annual regional panels without Rome monthly data. A past accidental download of individual `75N060E` archives for **July/August 2024** does NOT help the requested 2016–2019 window.

Richmond's 2017–2024 western-tile data should already be available from the prior 120-month 2015–2024 extraction (the script verifies both layers, not archive presence).

The pre-existing [manual archive importer](../../../src/import_eog_manual_monthly_archives.py) **already supports `--tile 75N060E`**, copying the two required layers into the existing per-year raw-data structure without deleting the source archives or mixing the east/west filenames. New code does **not** auto-download from EOG or guess user credentials.

## Execution

<pre>
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures

# Real synthetic rendering test; no actual raw data accessed
python src/plot_led_dominated_contrast_figures.py --self-test

# Check EOG V1 month availability separately for both tile sets
python src/plot_led_dominated_contrast_figures.py --inventory

# Generate the two FULL locked annual NUTS3 histories NOW
python src/plot_led_dominated_contrast_figures.py --annual-only

# Generate the Richmond borough monthly plot (WEST tiles already extracted)
python src/plot_led_dominated_contrast_figures.py --cases richmond
</pre>

Once the **2016–2019 `75N060E` monthly archive files** are downloaded to `~/Downloads`:

<pre>
python src/import_eog_manual_monthly_archives.py \
  --tile 75N060E --product vcmslcfg \
  --years 2016 2017 2018 2019 --dry-run

python src/import_eog_manual_monthly_archives.py \
  --tile 75N060E --product vcmslcfg \
  --years 2016 2017 2018 2019

python src/plot_led_dominated_contrast_figures.py --inventory
python src/plot_led_dominated_contrast_figures.py --cases rome
</pre>

If the boundary source fails, report it and obtain the authoritative local mask; **no automatic substitution of ITI43 or UKI75 for the municipality**. The initial programme window shading represents **multi-year installations**, not a known single switch-on date.

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
