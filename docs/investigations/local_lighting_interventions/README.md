# Local lighting interventions — five-case evidence hub

**Prepared:** 9 November 2026 · **Branch:** `publication-figures` · **Status:** Stage A documentary footprint investigation; NO local VIIRS/ISS effect estimates yet.

This folder consolidates our original five-case register, independently authored research briefs and source-report catalogue. The study tests whether changes in DNB radiance can arise from independently documented operational switching/dimming, spectral LED conversion and/or longer lighting hours. It is supplementary validation of the locked European study; it does **not** refit PRC curves or select treatments by satellite outcome.

## Comparative motivation

| Case | NUTS2016 | Intervention during locked 2014–24 window | Mechanistic comparison | Main strength | Main limitation |
|---|---|---|---|---|---|
| **Échirolles** | FRK24 Isère | Overnight switch-off Oct 2022; significant LED rollout after 2022 | Primarily operating hours **before** LED change | City-named permanently illuminated exceptions; independent municipal feeder electricity graph | Streetlight circuits and exact lamp points not yet geocoded; 2023 changes contaminate later controls |
| **Coventry** | UKG33 | Part-night switch-off implemented May–Jul 2024; main LEDs only Nov 2025 onward | Switch-off without contemporaneous *principal* LED conversion | **51-page switched-street list + XLSX** and council timeline | Lamp-by-lamp exemptions/dates and 2024 LED mix missing; city centre controls not exchangeable |
| **Portsmouth** | UKJ31 | 2018–19 LEDs, 2019 dimming pilot, citywide adaptive dimming summer 2020 | Spectral conversion then subsequent output change | Named trial neighbourhoods with reported unchanged placebo area | Trial placebo identity and dimming logs not yet obtained |
| **Thanet** | UKJ44 East Kent | LEDs from late 2016, move to optimised all-night operation | Potential **brightening** from longer hours despite LED conversion | Early/delayed cohorts; county current lamp inventory + GIS | Exact restoration/commissioning dates absent; delayed columns not random |
| **Bristol** | UKK11 | Rolling LED + CMS started Mar 2023, ongoing after Dec 2024 | Lantern+node vs node-only changes, staggered cohorts | Public geocoded assets, possible half-hour CMS telemetry | Historical 2023–24 schedules and install logs not publicly available |

## Documents here

- [Original five-case research register](./evidence_register_2014_2024.md) — consolidated context, original sources, satellite windows, limitations.
- [Stage A footprint matrix and decision rules](./stage_A_footprints.md) — what treatment and controls can actually be defined, and what is missing.
- [Primary source reports and original-document catalogue](./primary_source_reports.md) — one index of source reports/PDFs/XLSX and council records, with original URLs. **We link original reports, rather than reposting third-party PDF/Excel documents without verified redistribution rights.**
- [Data requests to local authorities](./data_requests.md) — targeted requests for asset/circuit IDs, commissioning and CMS logs.
- Five evidence dossiers: [Échirolles](./echirolles.md) · [Coventry](./coventry.md) · [Portsmouth](./portsmouth.md) · [Thanet](./thanet.md) · [Bristol](./bristol.md).

## Status of source documents

This folder contains **our complete analytical research documents** in GitHub. The external council, supplier and news report *originals* remain at source URLs and are catalogued with their formats, dates and uses. They have **not** been copied as PDFs into this repository: redistribution licences were not established, and the GitHub text-file workflow does not preserve those binaries. Do not mistake links for locally archived reports.

## Further reading

- Previously established ISS datasets: https://zenodo.org/records/7677478
- Existing spectral pipeline: `src/validate_iss_spectral_shift.py`, `src/analyze_iss_spectral_vs_locked_viirs.py`, `docs/results_iss_spectral_vs_locked_viirs_2014_2024.md`


## Analytical safeguards

- **NUTS2016 is a geographical container, not the treatment polygon.** Start from lamp/circuit/municipal documentation, never select treatments based on VIIRS brightening/dimming.
- Distinguish **street, column, lamp head, CMS node, feeder and VIIRS pixel**, including duplicate light heads and partial street coverage; prefer per-column timestamp and geometry.
- Call a source 'verified' only for what it actually establishes; current spatial inventories do not establish **historical** operation.
- Do **not** modify locked 2014–2024 European VIIRS, NUTS2016 cohort, 2017 zero baseline, 2021 snow correction, or fixed spider/moth PRC functions. Maintain standalone local outputs.
- Check actual nocturnal VIIRS overpass vs local switch-off clock time, BST/CET changes, snow/season, moon/cloud and effective observed-pixel coverage. Require pre-trends and compare within treated NUTS3 container, but examine spillover and treatment contamination. Existing 2020–22 monthly VIIRS collection is insufficient for several cases.
- ISS colour mosaics are irregular 2012–13 vs 2014–20 composites. No clean 2023–24 spectral inference from these composites; assess common-valid pixel counts and image dates before B/G and G/R comparisons.
- Energy use is not emitted luminous flux; B/G or G/R is not a unique spectral correction to DNB; no results are claimed before extracting local satellite pixels.

## Verified contacts and exact historical-record requests

See [verified contacts and five ready-to-send, location-specific request drafts](./verified_contacts_and_historical_record_requests_2026-11-09.md). All five city/county contacts were checked against official publications on 9 November 2026; **none of the messages has been sent**.

## Lean monthly event-window alternative (agreed 9 November 2026)

The [12-month VIIRS pilot protocol](./lean_12_month_viirs_event_study_protocol.md) explains why **municipality-scale time series can begin before column-by-column histories arrive**, and specifies cloud-free observation counts (`cf_cvg`), seasonal comparisons, transition-month flags and 2014–24 boundaries. This is the preferred **Stage B first pass**; it does not supersede the full Stage A evidence register, nor claim an observed effect.

## First local monthly extraction (9 November 2026)

**Implementation:** [src/local_viirs_monthly_pilot.py](../../../src/local_viirs_monthly_pilot.py). The script reads only the original EOG VNL V1 monthly `avg_rade9*` and `cf_cvg` rasters and uses: **Portsmouth UKJ31 NUTS2016** (2020–2021), **Échirolles INSEE commune 38151** (2021–January 2023, with French government contour cached separately). It selects pixel centres inside each city polygon and calculates mean/median radiance, valid-pixel proportion, median and total cloud-free observation counts for **both VCMSLCFG and VCMCFG**. It makes annual and calendar-month plots and CSV tables with source raster names. It has been syntax checked and smoke tested with synthetic rasters; **actual EcoVis outputs have not yet been measured**.

On EcoVis, from the existing repository root:
```bash
git fetch origin publication-figures
git show origin/publication-figures:src/local_viirs_monthly_pilot.py > src/local_viirs_monthly_pilot.py
python src/local_viirs_monthly_pilot.py
```
The `git show` path avoids switching branches or pulling unrelated work. Output is saved to `data/processed/local_lighting_interventions/monthly_pilot/` and `figures/local_lighting_interventions/monthly_pilot/`. Upload both city plots and CSVs to assess the first actual signal. This stage is descriptive, neither a controlled event-study nor a replacement for the locked annual NUTS/PRC results. Month-by-month cloud-free coverage and the 2021 snow anomaly require inspection before substantive interpretation.

## First user-run VIIRS outputs and preliminary findings

[Read the 9 Nov 2026 pilot findings](./first_monthly_pilot_findings_2026-11-09.md). Échirolles shows a marked Sep–Dec 2022 change relative to corresponding 2021 months, with full pixel coverage and comparable `cf_cvg`; Portsmouth does not yet show a clean city-scale shift. The pilot script was amended after those first plots to **include Échirolles 2020**, permitting two reference years without downloading anything new. Earlier "2023 Jan" file inventory was probably a **filename creation timestamp** error. These are preliminary descriptive findings, not intervention attribution. No locked annual results were altered.

## Native-grid spatial mapping extension (9 November 2026)

The current `src/local_viirs_monthly_pilot.py` now generates **six additional PNG/PDF maps** (absolute radiance, 4 × matched-month radiance differences, and per-pixel mean cloud-free observation counts for each of Échirolles and Portsmouth) plus two difference GeoTIFFs, without altering any locked outputs. See [spatial map specifications, limits, and the exact EcoVis `git pull` / run commands](./spatial_maps_monthly_pilot_2026-11-09.md). The maps are *defined and committed*; they are not yet a measured result until the script runs against the EcoVis monthly raster files.

## Actual native-grid spatial pilot results — 9 November 2026

[**Spatial results interpretation, paired-cell monthly differences, cloud-free observation counts and limitations**](./spatial_pilot_findings_2026-11-09.md) reviews all **six maps produced by EcoVis**, rather than just the plotting plan. **Échirolles** displays a substantial central/eastern and increasingly widespread late-2022 decline (mean Oct–Dec −38.7% vs 2021; paired-cell map median September −2.99, October −8.76, November −16.80, December −20.18 nW cm⁻² sr⁻¹), without a corresponding large collapse in `cf_cvg`. **Portsmouth** is spatially mixed and remains inconclusive with only 2020–21 monthly data. The uploaded **numeric CSV for Échirolles still covers 2021–22**, so additional 2020 seasonal-reference verification remains outstanding despite the script now supporting 2020. Interpretation remains descriptive pending a comparator/snow/atmosphere check; no locked scientific results are changed.

## EOG V1 monthly acquisition plan for all five local cases (9 Nov 2026)

[**Acquisition priority, exact missing years/months, repository downloader audit, and EcoVis commands**](./monthly_viirs_acquisition_plan_2015_2024.md). Verified from the investigator's EcoVis inventory that both monthly products are present in **2020–2022**, while **2015–2019** and **2023–2024** are missing (a 2023 January filename was a false-positive creation-date match). Download **2019 + 2023 + 2024** first (36 `VCMSLCFG` tiles), followed by **2016–2018** for Thanet (36 more), and optionally **2015** for a pre-conversion year (12 more). All five cities use the same `75N060W` EOG tile. The two existing GitHub auto-download scripts handled historical subsets only; **GitHub does not record how the raw files arrived on EcoVis** because they and their local `manifest.csv` files are gitignored.

New utilities: [inventory/EOG manifest audit](../../../src/inventory_monthly_viirs_acquisitions.py) and [manual archive importer](../../../src/import_eog_manual_monthly_archives.py). Run them only after a normal `git pull --ff-only origin publication-figures`; neither modifies locked annual products.

## Five-case 2015–2024 Stage B monthly analysis (new script, 9 Nov 2026)

New [`src/local_viirs_five_case_analysis.py`](../../../src/local_viirs_five_case_analysis.py) uses all **120 VCMSLCFG 75N060W observation months**, official **municipality/district polygons** (Échirolles INSEE; UK ONS LAD), the independently dated intervention windows, same-calendar-month reference years, source `cf_cvg` sensitivity (`>=1, >=3, >=5` observations), event-aligned monthly/decade plots and common-valid-pixel VIIRS spatial maps/difference GeoTIFFs. It also reports 2023/24 persistence or staggered-rollout descriptive contrasts where interpretable. The **Thanet** mask is now its own LAD E07000114, not the much larger UKJ44 East Kent NUTS3 region.

**[Read full predeclared protocol, data caveats, and exact EcoVis commands](./full_five_case_monthly_analysis_protocol_2026-11-09.md).** All data/figures are generated locally under new `five_case_monthly` output directories and do **not** modify the locked 2014–24 European NUTS3 radiance or spider/moth PRC work. Code committed, but real EcoVis output has **not yet been observed/reviewed**.
