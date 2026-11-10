# Local lighting interventions — five-case evidence hub

**Prepared:** 9 November 2026 · **Branch:** `publication-figures` · **Status:** Stage A documentary evidence plus Stage B **measured 2015–2024 monthly VIIRS pilot results** for all five cases; causal attribution, ISS event-specific validation and snow/atmospheric controls outstanding.

**Primary purpose — not merely spatial mapping:** [Agreed four-point plan for distinguishing real NUTS3 lighting reductions from HPS→LED spectral artefacts](./four_point_attribution_plan_2026-11-10.md). The local analysis is supporting evidence for the **locked 2014–2024 European VIIRS radiance trends and inferred spider/moth PRC changes**. The four steps are **independently documented lighting operations/footprints → locally matched VIIRS trends and NUTS3 bridge → dated ISS B/G and G/R spectral context → triangulation without forcing a universal LED correction**. Do not add further arbitrary subcity analyses unless they discriminate LED-only vs operational dimming. Prioritise **Coventry 2024 and Échirolles 2022** as documented switch-offs, with actual commissioning/LED timing and NUTS3 linkage; Portsmouth 2020 remains inconclusive as an additional VIIRS event.

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

## First full five-case monthly results reviewed — 9 November 2026

[**Read initial observed results, full city table, spatial interpretations and prioritised diagnostic plan**](./five_case_initial_results_2026-11-09.md). These results were provided as actual EcoVis `five_case_contrast_summary.csv` and ten PNG plots. Primary season-matched municipal means: **Échirolles −37.08% (late 2022), Coventry −32.23% (late 2024), Portsmouth −15.11% (late 2020), Bristol +0.02% (Apr–Dec 2023), Thanet +9.40% (Jun–Nov 2017)**. Important independent caveats: Portsmouth's decline predates its 2020 summer dimming; Bristol 2024 November has very low `cf_cvg`; Thanet's brightening map is dominated by a compact unidentified hotspot. These are descriptive satellite contrasts and **not causal estimates or lumen changes**.

[**Standalone sensitivity/earlier-year placebo script**](../../../src/local_viirs_five_case_robustness.py) reads the 120-month city CSVs and generates outcome-year monthly anomalies, explicit pre-versus-post anomaly differentials, prior-year placebo comparisons, cloudfree `>=3/5` thresholds, low-cf month filtering and leave-one-month-out diagnostics. Source code is committed; its output has not yet been observed.

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/local_viirs_five_case_robustness.py --self-test
python src/local_viirs_five_case_robustness.py
```

The five-case monthly Stage B scripts and their standalone outputs do not change the locked European figures/models.

## Observed five-city event anomalies and timing diagnostic (9 November 2026)

[**Reviewed plots and interpretation**](./first_robustness_event_anomaly_plots_2026-11-09.md) record the five actual EcoVis event-anomaly PNGs. Key update: **Portsmouth was already strongly negative versus 2018–19 in Jan–May 2020**, so its Sep–Dec −15.1% season-matched contrast must **not** be attributed to the summer adaptive-dimming programme. **Coventry** shows a cleaner sustained downward departure after May–July 2024 implementation; **Échirolles** develops the largest negative anomalies in Nov–Dec 2022. **Thanet** +9.4% is dominated by October–November and a compact hotspot; **Bristol** 2023 is mixed while rollout continued throughout the year.

Updated [robustness script](../../../src/local_viirs_five_case_robustness.py) shades **March–December 2023** as Bristol's ongoing LED/CMS change (not only March), and adds a labelled **Bristol 2024 secondary sensitivity** including explicit omission of low-observation **November 2024**. The newly requested exact three numerical robustness CSVs remain for review; no causal effect has been inferred.

## Spatial dilution test first — new fixed-baseline analysis (9 November 2026)

**[Full spatial dilution design, actual three robustness CSV findings, output definitions and EcoVis commands](./spatial_dilution_analysis_plan_2026-11-09.md)**. Separate [`src/local_viirs_spatial_dilution.py`](../../../src/local_viirs_spatial_dilution.py) extracts the relevant monthly native VIIRS cells and groups them by **reference-year brightness ranks (bottom 50%, 50–80%, 80–95%, top 5%)** and **fixed NW/NE/SW/SE geographic sectors**. The pixel groups are defined **without inspecting outcome-year changes**; they are descriptive, **not independently verified lamp-treatment polygons**. For all five cities, the script produces strict same-pixel, same-calendar-month pre-event/post-event anomalies where defensible, additive net versus gross subarea contributions, `cf_cvg >=1/3/5` thresholds, CSVs and maps. Secondary follow-ups are explicitly marked as already exposed/staged, not assumed to have untreated pre periods.

**Uploaded robustness CSVs indicate:** Échirolles pre→post anomaly **−14.52** nW cm⁻² sr⁻¹; Coventry **−3.91**; Portsmouth **+1.86** (the pre period was already more negative); Bristol 2023 **+0.37**. **Thanet +9.4%** in 2017 turns **negative (~−1.8%) if October or November is omitted**. Bristol 2024 −18.2% becomes **−11.17%** excluding low-`cf_cvg` November. All are descriptive differences, not fitted causal estimates. Full numerical tables are in the linked design.

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/local_viirs_spatial_dilution.py --self-test
python src/local_viirs_spatial_dilution.py --cases portsmouth bristol
python src/local_viirs_spatial_dilution.py
```

**State:** code and methodology committed, **actual EcoVis spatial-dilution outputs not yet run/reviewed**. Locked annual European VIIRS and PRC results remain unchanged.

## Spatial maps interpreted and Portsmouth actual rollout checked — 10 November 2026

[**Full visual interpretation and date audit**](./spatial_figures_and_portsmouth_date_audit_2026-11-10.md). The investigator identifies Bristol's bright northwest baseline as **Avonmouth industrial docks** and a southern source as **Ashton Gate Stadium surroundings**. Importantly, **the northwestern source shows a positive 2023 post-minus-pre anomaly, while the Ashton Gate-area hotspot shows a NEGATIVE post-minus-pre anomaly**, despite remaining bright relative to past calendar-month references. Thus strong **2023 spatial cancellation** is real descriptively, but it would be erroneous to say both facilities brightened after intervention or that either is definitively the source of the radiance. Bristol's 2024 stage is different: all four pre-event baseline-brightness groups decline, subject to November low-cloudfree coverage caveats.

**Portsmouth dates now audited against original municipal records:** 2019 trial reported underway in February; council approved citywide adaptive dimming on **4 Feb 2020**; crucially, **official June 2022 council guidance states rollout actually occurred in summer 2020**, not merely planned. Individual streetlight activation dates remain unknown, and exceptional routes/heritage lamps were not homogeneous. Extensive 2018/19 LED installation and a 2019 dimming pilot can explain some pre-2020 confounding; 2021 heritage LED work continued. We **retain June–Aug 2020 as an approximate intervention window**, not relocate it to 2021 just to align satellite data. See [Portsmouth dossier date correction](./portsmouth.md) and [full source catalogue](./primary_source_reports.md).

Next useful output: **`portsmouth_primary_2020_spatial_anomalies.png`**, which was not among the most recent set (the user supplied the secondary **2021** map). This would show the actual 2020 before/after seasonal spatial difference, not incorrectly substitute 2021.


## Coventry / Échirolles — regional attribution bridge, 10 November 2026

**[Evidence table: actual council implementation dates, LEDs vs genuine operating-hour switch-offs, and frozen 2014→2024 UKG33/FRK24 radiance+PRC endpoints](./coventry_echirolles_nuts3_radiance_attribution_table_2026-11-10.md).** Coventry's May–July 2024 part-night rollout covered ~60% of eligible lamps, confirmed retrospectively by council; the principal LED retrofit was subsequently approved January 2025. A 2024 council Equality Impact Assessment independently notes some **pre-existing ~55%-output residential dimming and Hillfields LEDs**. Échirolles' municipal roadway switch-offs were **phased October–December 2022**, not all completed exactly 1 October, with earlier August switches and major 2023 LED conversion.

Frozen regional **2014→2024 endpoints**: **UKG33 Coventry 18.676→13.030 nW (−30.23%)**; **FRK24 Isère 2.074→0.880 nW (−57.58%)**. Those are decade-long changes, **not demonstrated local policy effects**. Échirolles is only a portion of Isère. The older source CSV with FRK24 2024 ~0.936 is not the final locked dataset and is deliberately excluded from annual time-series reconstruction.

**[Read-only locked annual bridge exporter](../../../src/export_coventry_echirolles_nuts3_bridge.py)** provides the 11-year 2014–2024 annual radiance histories, validated original 2014/2024 endpoints and existing audited **2021 native-pixel snow** cache, alongside 2023→2024 UKG33 and 2021→2022 / 2022→2023 FRK24 year-on-year tests. It does not alter original annual inputs, snow correction, PRC or figures. **The new script is committed but has not yet been executed on EcoVis**; intermediate event-year annual values are not fabricated in documentation.


## Two documented switch-off figures and two LED contrast audits — 10 November 2026

**Step 1 — plots:** [\`src/plot_coventry_echirolles_annual_monthly_attribution.py\`](../../../src/plot_coventry_echirolles_annual_monthly_attribution.py) makes one three-panel figure each for **Coventry UKG33 and Échirolles/FRK24**: **(A) frozen annual 2014–24 VNL V2 NUTS2016**, **(B) municipal monthly VNL V1 in key surrounding years**, **(C) intervention-year monthly same-calendar-month anomalies**, all with independently documented approval/start/phase-in and later LED dates. Uses original locked Figure 1 2021 snow cache and will abort if source audits are unavailable. **[Exact figure and source protocol](./annual_monthly_attribution_figures_2026-11-10.md).** Code committed; actual figures await EcoVis run.

**Step 2 — historical LED contrast candidates and negative source audit:** [**Rome (ITI43)** and **Richmond upon Thames (UKI75)**](./led_dominant_contrast_candidates_2026-11-10.md) are two prospective LED-*dominated*, **NOT yet confirmed undimmed** cases. Rome is unusually informative because independent official data show **2017→2018 streetlight total lumens approximately constant (+1%) while electricity use fell 27.4%** during major LED conversion; by contrast **total installed lumens actually fell greatly over 2014→2018** so no universal LED-only claim is permissible. Richmond's 2019–2022 scheme pledged the same illumination and unchanged burning hours, **but the council explicitly uses dimming regimes**, so it needs actual 01:30 profiles before identification as an undimmed contrast. Original Plymouth 2014 paperwork confirms dimming in main-road LED installations and shorter switch-on hours; Milan and Turin also integrated dimming. These three are **excluded from strict LED-only controls**.

The relevant **Kyba et al. (2021)** [VIIRS-source partition study](https://doi.org/10.1177/1477153520958463) estimates city-owned streetlights contributed ~13% of Tucson's overpass-time DNB radiance with normal dimming (95% CI 10–16%); it is a measured demonstration why privately operated sources must not be assumed to follow municipal policy.

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/plot_coventry_echirolles_annual_monthly_attribution.py --self-test
python src/plot_coventry_echirolles_annual_monthly_attribution.py
\`\`\`

Both stages are designed to test the primary **genuine operating-hour reduction vs purely LED-spectral explanation** question. No source or frozen PRC arrays are recalculated.

## 10 Nov 2026 — reviewed switch-off figures and prospective LED comparison plots

[**Interpretation notes for the actually generated Coventry and Échirolles 3-panel figures**](./coventry_echirolles_figure_interpretation_2026-11-10.md) now distinguish their **marked post-implementation monthly municipal declines** from the **nonlinear longer-run NUTS3 annual trajectories**. Coventry's new 2024 switch-off follows May–July programming, and its principal LED programme is later than 2024. Échirolles has phased 2022 operations, but is only a **small part of Isère**, and its 2023 LED retrofit complicates attribution of the regional 2022→2023 step. Notes retain the important distinction between real operations and the unknown contribution of non-municipal private sources. **No claim of a quantified causal contribution to the full 2014–24 regional decline.**

[**New Rome and Richmond LED-dominated comparison plotting code**](../../../src/plot_led_dominated_contrast_figures.py) and [**reproducible run/input protocol**](./led_dominated_contrast_plotting_protocol_2026-11-10.md) reproduce the same 3-panel layout: **(A)** frozen 2014–24 NUTS3 annual radiance using the existing 2021 audit-verified snow cache, **(B)** monthly *municipality* VNL V1 radiance and independently dated staged LED periods, and **(C)** a same-calendar-month focus-year comparison. Both are **prospective LED-dominated**, NOT confirmed *without dimming*. Rome has independent 2017→2018 **~+1% installed lumens with −27.4% estimated electricity consumption**, but the *overall 2014–18* installed lumens fell markedly; Richmond's stated unchanged hours do **not** rule out dimmed output at VIIRS overpass.

**Corrected European tile identification (10 Nov 2026):** Both **Rome (~12.5°E)** and **Richmond** are contained in **EOG Tile 2 (60W75N), filename `75N060W`**, exactly the same tile used for Coventry, Échirolles and the other original cities. The earlier instruction to download **`75N060E` for Rome was incorrect**: that is EOG Tile 3 (60E75N), East Asia/Pacific. The already extracted 2015–2024 `75N060W` monthly `VCMSLCFG` rasters should cover Rome 2016–19 and Richmond 2017–24; `--inventory` checks raster-pair completeness. **No additional imagery or eastern-tile acquisition is required** for this comparison. The plotting script is corrected accordingly, but Roman municipality-boundary availability is a separate validation gate. Existing annual/PRC data remain unchanged.

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/plot_led_dominated_contrast_figures.py --self-test
python src/plot_led_dominated_contrast_figures.py --inventory
python src/plot_led_dominated_contrast_figures.py --annual-only
python src/plot_led_dominated_contrast_figures.py --cases richmond
```

**Status:** Scripts/docs committed; new Rome and Richmond real-data plots remain **unrun in EcoVis** and cannot be interpreted until outputs are reviewed.

## Rome municipality boundary repaired — 10 November 2026

On EcoVis, the user ran the corrected **75N060W** LED comparison script. `--self-test` passed and `--inventory` confirmed **48/48 Roma 2016–19** and **96/96 Richmond 2017–24** month pairs, but the Rome runtime failed because the old river-authority ArcGIS endpoint did not supply a valid GeoJSON municipality polygon. **This was not a missing raster or annual VIIRS failure.** The new plotting script uses an **ISTAT-derived GeoJSON subset for Rome province**, matching **municipality ISTAT `058091` / cadastral `H501`**, with strict geometry/location/area validation before caching; the annual NUTS3 ITI43 remains independent of the municipal mask. See [full issue/fix](./led_dominated_contrast_plotting_protocol_2026-11-10.md#rome-municipal-boundary-runtime-fix--10-november-2026).

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/plot_led_dominated_contrast_figures.py --self-test
python src/plot_led_dominated_contrast_figures.py --check-boundaries --cases rome
python src/plot_led_dominated_contrast_figures.py --cases rome
```

The corrected Rome boundary source and tests are committed; **Rome's full monthly result has not yet been rerun successfully on EcoVis**. Richmond's existing completed plot does **not** need regeneration.


## Rome and Richmond: real LED-dominated VIIRS comparison results — 10 November 2026

**[Reviewed actual four figures and five CSVs, with annual/monthly results and full four-case mechanistic synthesis](./rome_richmond_real_viirs_results_2026-11-10.md).** Both municipalities' real EOG VNL V1 plots are successfully generated with the corrected European 75N060W tile and distinct local boundary masks. No new raw imagery, snow cache, PRC or annual inputs modified.

**Rome stands out:** municipal VIIRS monthly **2016→17 −11.6%**, **2017→18 −9.1% (11/12 months lower)**, **2018→19 −0.3%**, corresponding to ITI43 NUTS2016 annual **2016→18 −13.0%** and post-2018 rebound to 2024. Critically, Rome ACoS independently reports *installed streetlighting luminous flux* **2016→17 −27.6%**, but **2017→18 +0.95%** while *estimated energy use* **−27.4%**. Thus **VIIRS dropped 9.1% in 2017→18 despite approximately constant *rated installed streetlighting lumens***. This is an important mechanism contrast **consistent with** LED spectral/directional changes, but *not proof*: actual output at overpass time, non-municipal sources and switching/dimming schedules remain unverified. The 2016→17 interval included a real drop in installed lumens. The whole conversion period is not a no-dimming/constant-lumen experiment.

**Richmond is not a clean LED-only control:** Richmond borough's 2022 municipal VNL V1 is **−19.7%** below the 2017/18 same-month baseline, but radiance was already declining pre-2019 and the funding document's 'unchanged hours' does not rule out dimming. The 2023→24 borough decline is **−13.9%** overall but only **−4.2% in Jan–Oct**: a large portion comes from November–December; UKI75 NUTS3 (includes Hounslow) fell **−11.6% 2023→24**, which cannot automatically be attributed to Richmond's 2019–22 LED rollout.

**Priority next:** a direct Rome 2016–18 **installed-lumens / regional annual / municipal monthly** juxtaposition (different units/scales explicitly marked), independently dated spectral ISS data if available, plus operational control profiles at VIIRS overpass. This preserves the two-sided attribution argument: real operating reductions **can** cause VIIRS declines (Coventry/Échirolles); LED transitions can also change VIIRS in ways that do not straightforwardly track rated luminous flux (Rome).


## Madrid: ground-based sky colour and frozen VIIRS, 10 November 2026

**[Reproducible Madrid dataset/figure protocol, exact paper Table 6 and source links](./madrid_viirs_astmon_comparison_protocol_2026-11-10.md).** A new independent [\`src/plot_madrid_viirs_robles_2021.py\`](../../../src/plot_madrid_viirs_robles_2021.py) assembles: **(A)** locked annual **ES300 NUTS2016 2014–24** with the existing audited 2021 snow patch; **(B)** municipality-only monthly **75N060W VNL V1 VCMSLCFG 2015–20** (or through 2024); **(C)** directly published **Robles et al. (2021) Table 6** after-midnight UCM AstMon Johnson B/V/R medians, normalized **2015=100 within each spectral band** with physically correct \`10^(-0.4 Δm)\`; and **(D)** AstMon B−V / V−R colour. The plotted 2011/2015/2019 AstMon points are the *only* published Table 6 years; no interpolation, invented ground readings or conflation of ground zenith brightness with VIIRS upward emissions.

[The entire original UCM observational archive is also available as AstMon/SQM/TESS CSVs on Zenodo](https://doi.org/10.5281/zenodo.4633001) (about 140 MB total). A future full-year raw reanalysis must reproduce the paper's no-moon/no-cloud, after-midnight filters before comparing annual medians. Broad HPS-to-LED conversion 2015–16 also changed **luminaire wattage, shielding and after-midnight dimming** (HPS 80%, LED 60% according to the paper), making Madrid a **mixed-mechanism**, not LED-spectral-only control.

The city's mask is **official Ayuntamiento Madrid TERM_MPAL FeatureServer layer 4 (~605 km²)**, validated and distinct from **ES300 Comunidad de Madrid**, with NO use of the wrong east Asia tile. This work does not modify locked PRC/annual inputs. **State: script/docs committed, real EcoVis figure NOT yet generated or reviewed.**

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/plot_madrid_viirs_robles_2021.py --self-test
python src/plot_madrid_viirs_robles_2021.py --inventory
python src/plot_madrid_viirs_robles_2021.py --check-boundary
python src/plot_madrid_viirs_robles_2021.py
\`\`\`


## Madrid: SKY BRIGHTNESS, not merely bluer colour — revised second figure

After reviewing the investigator's actual initial four-panel Madrid VIIRS/UCM AstMon figure, the primary mechanistic outcome is **changes in absolute ground-observed zenith sky brightness within spectral bands**, not simply B−V / V−R shifts. Direct 2015→2019 published clear, moonless after-midnight AstMon Table 6 implies **Johnson B +30.6%**, **V +11.7%** (not a significant long-run trend according to authors), **R −1.8%**. [Robles et al. (2021)](https://doi.org/10.3390/rs13081511) separately reports **TESS-W sky brightness +10% over 2016–19**. By contrast **SQM showed −38% 2012–19**, which the authors associate importantly with **aging/filter transmission decline**, *not* a trustworthy physical sky-darkening result. After 2015 little additional Madrid municipal streetlight replacement was reported; neighbouring communities, independent sources and diffuser/operational change may contribute, so no causal LED-only claim.

New [**\`src/plot_madrid_brightness_viirs_direct_contrast.py\`**](../../../src/plot_madrid_brightness_viirs_direct_contrast.py) creates an additional direct **2015→19 within-instrument percentage-change brightness** figure, separates ES300 annual VIIRS and Madrid municipality annualised monthly VIIRS (both indexed to 2015), and shows discrete measured 2011/2015/2019 **actual photometric B/V/R sky brightness**. **Reuses existing real Madrid annual/monthly/published-ground CSVs**; no raster processing, annual source rebuild or PRC change. Complementary earlier figure preserves 2015–20 municipal monthly detail and B−V/V−R colour shifts. [Full source/method notes](./madrid_viirs_astmon_comparison_protocol_2026-11-10.md#ground-sky-brightness-is-a-distinct-central-endpoint--correction-and-second-diagnostic-figure).

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/plot_madrid_brightness_viirs_direct_contrast.py --self-test
python src/plot_madrid_brightness_viirs_direct_contrast.py
\`\`\`

**Interpretation guardrail:** The ground skyglow measurement is zenith and site-specific; satellite measurements are upward/area-based. Relative radiance can diverge for spectrum, direction, weather, timing and other-source reasons. No direct animal polarization contrast can be computed from ground B/V/R brightness alone.


## Madrid UCM fixed-radius satellite/ground study — 10 November 2026

After examining Madrid's real ES300/municipality VIIRS versus UCM AstMon plots, the next study tests the **spatial-scale confound**. The new [**\`src/madrid_ucm_local_ground_bridge.py\`**](../../../src/madrid_ucm_local_ground_bridge.py) calculates 2015–19 original \`75N060W\` EOG VNL V1 monthly VIIRS over **fixed 3 km and 10 km circles** around the precisely published UCM observatory coordinates (IAU I86), compares with the already extracted **municipality-month VIIRS**, and performs **strict same-calendar-month 2015 vs 2019 paired-pixel tests** at \`cf_cvg\` thresholds 1, 3 and 5. **No annual NUTS3 or PRC dataset changes.**

The original Zenodo AstMon CSV was **successfully MD5-verified and inspected on EcoVis**. Its nine fields contain **NO cloud/moon quality flags**; `Filter` identifies Johnson B/V/R and `Magnitude` contains each observation. The authors manually quality-labelled ~28% of the data for their paper, but those flags are absent from this downloaded file. **Do not fabricate clear-sky labels or run `--prepare-ground` on unfiltered records.** The script has been amended to recognise the actual long-form source and require independently reviewed QC annotations. Meanwhile use the paper's published quality-filtered Tables 3–6 and supplementary Tables S1–S3.

[**Complete methodology, data provenance and staged EcoVis commands**](./madrid_ucm_local_astmon_viirs_plan_2026-11-10.md). The script and methods are committed; Python syntax and **synthetic** ground-flag, geographic-area and PNG/PDF/SVG tests were run in the working environment. **Real UCM-local EcoVis satellite results and the raw Zenodo quality flag schema have not yet been inspected or calculated.**

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/madrid_ucm_local_ground_bridge.py --self-test
python src/madrid_ucm_local_ground_bridge.py --viirs
python src/madrid_ucm_local_ground_bridge.py --download-ground --inspect-ground
\`\`\`

**Ground quality blocker confirmed:** The actual source has no clear/moon/cloud flags, so no `astmon_columns_flags.json` can be legitimately derived from it. Author annotations or the independently published filtered statistics are required to extend the ground time series.

## Madrid UCM local 3/10 km results received and validated — 10 November 2026

**[Actual UCM-local VIIRS results, month-by-month pixel support and winter sensitivity](./madrid_ucm_local_results_2026-11-10.md).** The investigator uploaded the successfully generated `madrid_ucm_local_viirs_2015_2019.png` and all three underlying CSVs. In **2015→19**, Madrid municipal VNL V1 **−1.94%**, UCM **3 km −2.60%** and **10 km −3.78%**, whereas published ground UCM after-midnight Johnson B sky brightness **+30.6%**, V **+11.7%**, R **−1.8%**. **The broad blue-skyglow/upward-DNB discrepancy persists at both local satellite scales**; the circles are not full skyglow source footprints and cannot identify the cause.

All 60 monthly VIIRS composites had **172/172 valid cells at 3 km and 1,913/1,913 at 10 km**, unchanged even at `cf_cvg≥5`, so local pixel support is not driving this result. However declines are **concentrated in winter**: March–October 2015→19 changes just **−0.83% at 3 km and −0.88% at 10 km**. The yearly local VIIRS fell most strongly 2015→16 then partly rebounded to 2019, not a uniform sustained decrease. There is still no monthly, independently verified, clear-sky UCM B/V/R ground series to align season-by-season: **the raw AstMon archive has been inspected and has no author quality flags**.

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/madrid_ucm_local_ground_bridge.py --download-ground --inspect-ground
```

There are no quality-code distributions in the released AstMon file. Do not use unfiltered ground data or claim that the measured B brightening proves a uniquely spectral VIIRS bias, rather than a mixture of lamp spectrum, directionality, light output, schedules, scattering and independent sources. Frozen European 2014–24 analyses/PRC are unchanged.


## Madrid AstMon CSV reviewed; published winter 01:30 brightness comparison — 10 November 2026

**Actual EcoVis result:** The Zenodo [AstMon_evol.csv](https://zenodo.org/records/4633001) was MD5-verified (\`634f098973251bf41f7d8a12b0a8ab45\`) and examined. The **nine fields** are \`Julian_Dates\`, \`Datetime_Raw\`, \`Datetime_UTC\`, \`Filter\`, \`Pixel_x\`, \`Pixely\`, \`Pos_0\`, \`Magnitude\`, \`Error\`. This is **long-form B/V/R photometry, with NO cloud/moon quality flags**. The \`Filter\` field is the photometric filter (e.g. Johnson-B), and \`Error\` is not a sky-condition indicator. The authors' manual clear/no-moon selection (119,296 records, 28%) cannot be reconstructed from an unlabelled CSV by relabelling these fields. The existing script's \`--prepare-ground\` still fails closed; parsing is updated for the correct long schema if separately validated quality annotations are obtained.

**Independent source available now:** [Robles et al. (2021), Tables 3–5](https://doi.org/10.3390/rs13081511) give published filtered **winter 01:30** UCM zenith sky brightness, 2015→2019:
- Johnson B **19.30→18.97 mag/arcsec²: +35.5% relative linear brightness**.
- Johnson V **18.07→17.96: +10.7%**.
- Johnson R **17.42→17.45: −2.7%**.

The paper calls **November–January** winter. The directly computed UCM **VIIRS 2015→19** same-calendar-year **January, November, December** means change **−5.9% at 3 km**, **−7.7% at 10 km** (\`cf_cvg>=3\`, all cells valid). **Important:** author winter-year assignment may cross the calendar boundary, and satellite observations are monthly averages, not synchronised with photometry; this is a nearby-season sensitivity contrast, **NOT exactly time-matched or a causal spectral-effect estimate**.

[\`src/plot_madrid_ucm_published_winter_sky_brightness.py\`](../../../src/plot_madrid_ucm_published_winter_sky_brightness.py) produces a PNG, editable SVG and PDF of the independently published filtered 01:30 B/V/R versus actual UCM-local winter satellite percentages. It uses the existing monthly CSV and published table entries, not unfiltered observations or new VNL rasters.

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/madrid_ucm_local_ground_bridge.py --self-test
python src/madrid_ucm_local_ground_bridge.py --inspect-ground
python src/plot_madrid_ucm_published_winter_sky_brightness.py --self-test
python src/plot_madrid_ucm_published_winter_sky_brightness.py
\`\`\`

For complete **published quality-filtered ANNUAL** B/V/R 2014–19, seek the article's **Supplementary Tables S1–S3**. For monthly observational records, seek the original author's per-observation manually reviewed cloud/moon labels; neither exists in the original nine-column Zenodo file. Frozen annual VIIRS/PRC unchanged.
