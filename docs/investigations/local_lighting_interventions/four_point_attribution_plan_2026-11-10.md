# Four-point attribution plan: do European NUTS3 VIIRS declines represent real lighting reduction or LED spectral artefacts?

**Restated 10 November 2026 at investigator's request.** This is the *agreed scientific objective* behind the ancillary local-lighting cases, not a new study goal. **No recalibration, spectral-refit, addition to PRC response curves or modification of the locked European 2014–2024 analysis is authorised.**

## Primary scientific question and threshold of evidence

The locked study shows regional **VIIRS-DNB radiance changes** and calculates the corresponding changes in **spider and moth celestial polarization receptor contrast (PRC)**. Its biological interpretation needs an independent plausibility check: **where VIIRS radiance fell, is there evidence that actual artificial lighting also fell, rather than the satellite seeing a different lamp spectrum after HPS-to-LED conversion?** This test is especially valuable for the subset of NUTS3 regions with declining DNB radiance and consequential PRC differences.

**Modest, defensible evidentiary objective:** demonstrate convincingly that **some** of the historical NUTS3 radiance decreases reflect real operating-hour switch-offs, command-level dimming or reductions in luminous emission, **and therefore cannot universally be explained by spectral substitution alone**. It is **not** necessary to argue that LED conversion always decreases DNB radiance, that LED lighting is uniformly brighter/darker at ground level, or that every declining NUTS3 region was dimmed. If a LED-converted place becomes brighter or unchanged in DNB, that is a useful counterexample, not a failure. Nor should the data be used to estimate a universal LED-vs-operations percentage partition from DNB alone.

**Physical caveat:** The VIIRS DNB response is weak at the ~450 nm blue pump typical of phosphor-converted white LEDs, so HPS→LED technology changes may reduce *satellite-measured DNB radiance* even with unchanged or increased broadband light output; changes in luminaire direction, lamp power, total intensity and atmospheric conditions also matter. The exact effect depends on actual spectra, outputs, emission geometry and instrument response. **ISS RGB ratios constrain relative spectral/colour changes but cannot independently quantify full luminous flux**, and actual ecological polarization cues need independent spectral/biological interpretation if they are to be called physically observed rather than modelled from DNB.

## The agreed four points — preserve this order and scope

### 1. Independently document and spatially define real lighting interventions

Collect council decisions, **actual implementation dates**, switched-street or circuit inventories, historical lamp/CMS schedules (including overnight hours), energy/power or pre/post lamp-flux evidence, and HPS/LED replacement dates. **Separate approval and proposed completion from documented commissioning**. Define intervention polygons/asset points *before examining VIIRS outcome sign*; do not label darkening satellite pixels as «treated» by hindsight. Map each municipal footprint into its fixed **NUTS2016/NUTS3 container** and note whether the municipality occupies most or only a small part of the region.

Examples: Échirolles Oct 2022 street-light shut-offs when 2022 municipal LED adoption was still small; Coventry May–Jul 2024 part-night switch-offs before main 2025 LED upgrade; Portsmouth 2018/19 LED, 2019 pilot, actual summer 2020 citywide adaptive dimming; Thanet LED plus longer all-night operating hours; Bristol staggered 2023/24 LED+CMS.

### 2. Test the documented changes against the actual local VIIRS time series

Use fixed-geometry and same-calendar-month **2015–2024 EOG V1** \`avg_rade9h\` and \`cf_cvg\`, with intervention and pre-trend windows; compare documented treated-area exposure to unaffected/less-exposed geography **only where independent geometry exists**. Check actual overpass vs off/dimming clock hours; cloud-free observation counts, winter snow/ground reflectance, atmospheric effects, isolated hotspots and interannual drift. An observed radiance fall *after* an independently documented switch-off (without major simultaneous replacement) is direct evidence that **some part of the DNB change reflects lamp operation**, rather than a universal spectral explanation.

**Crucial bridge to the main paper:** quantify whether the municipal footprint's contribution to its **parent NUTS3** radiance history is material, and whether the NUTS3 2014–2024 trajectory/timing is consistent. Échirolles commune **is not equivalent** to the much larger **FRK24 Isère**; a strong Échirolles signal alone cannot prove what caused Isère-wide decline. Coventry's local authority largely matches **UKG33** and may be a stronger regional bridge. Preserve the annual locked series; perform this comparison as **separate descriptive validation**, not by altering annual data processing.

**Current results:** Coventry and Échirolles show strong, spatially widespread, event-aligned radiance declines and robustness to \`cf_cvg\` thresholds. Portsmouth's early-2020 pre-window was **already more negative** than Sep–Dec 2020 relative to 2018/19; its 2020 spatial post-minus-pre map is predominantly positive. This **does not** isolate an incremental summer-2020 DNB fall, despite council retrospective confirmation that dimming was rolled out in summer 2020. Bristol shows opposing source trends in 2023 (Avonmouth industrial area vs Ashton Gate surroundings) and 2024 sparse-month sensitivity. Thanet's 2017 positive outcome is mostly highest-brightness-source and Oct/Nov sensitive.

### 3. Independently examine local changes in ISS RGB spectral ratios

Where ISS **date-appropriate, geometrically matched** images are available, examine **B/G and G/R** on fixed common-valid local pixels; independently test whether a move towards white/blue-rich LED spectra accompanies a local VIIRS decline. Use matched source dates rather than assuming an ISS mosaic epoch means "2020". Do not mechanically use RGB ratios as a quantitative correction to DNB without lamp SPD and response modelling.

**Critical temporal limitation:** Existing calibrated ISS European mosaics have broad **2012–2013** and heterogeneous **2014–2020** epochs. They **cannot directly test** Échirolles' **October 2022** switch-offs or Coventry's **2024** switch-offs. For those events, **the historical LED status and verified operation logs** are stronger identification evidence. ISS is more germane to earlier HPS→LED conversions (e.g., Portsmouth 2018/19) and to country/NUTS3-level *spectral context*, not retrospective proof for post-2020 operational events.

**Existing project support, not to be repeated:** \`docs/results_iss_spectral_vs_locked_viirs_2014_2024.md\` and \`src/analyze_iss_spectral_vs_locked_viirs.py\`. Country-scale ISS shifts and DNB declines have already been compared. France has a large DNB decrease but near-zero median ISS Δ(G/R), though the ISS epochs do not align exactly with 2024 and cannot alone certify emitted-flux reduction.

### 4. Triangulate mechanisms, do not force LED-only versus dimming-only

Combine **(a)** outcome timing/location in local VIIRS, **(b)** ISS colour/spectral shifts where dated imagery is informative, and **(c)** independent on/off, CMS command/energy and LED replacement records. Frame each case as operational reduction, spectral change, **both**, increased hours despite LEDs, or unresolved. The desired inference is **a supported counterexample to the claim that all apparent NUTS3 darkening is solely an LED-induced DNB artefact**, not exact wattage or lumens inferred from satellite. LED conversions with stable or **higher** VIIRS radiance and switch-offs with real reductions are complementary evidence.

**Strongest mechanisms to establish first:** Coventry 2024 and Échirolles 2022. **Portsmouth** may separate LED conversion from operational dimming in documentary records, but measured DNB has no demonstrated incremental summer-2020 step; **do not imply otherwise**. **Thanet** is a possible longer-hours/LED brightening contrast, but currently contaminated by two months and an unidentified hotspot. **Bristol** demonstrates why independently operated emitters create different spatial and temporal trends within one municipality, but is not presently evidence for an isolated LED-versus-dimming difference.

## Sentence the investigator wants retained

> Municipal streetlighting policies, industrial lighting, sports facilities and other independently operated sources can produce different spatial and temporal trends within the same city.

Use that sentence as a **specific caution about aggregation**, not a replacement for the four-point causal-mechanism/PRC validation objective.

## Explicit decision gate / efficient remaining work

1. **Do not keep adding arbitrary spatial plots.** The five-city exploratory maps have identified where the next increment of evidence is most valuable.
2. **Coventry:** secure and geocode the council's published 2024 switched-street list, document any 2024 LED subset and actual switch-off hours. Directly test overlap with the negative 2024 native-grid cells without picking target cells from the data. Establish **UKG33 NUTS3** bridge against the locked trend.
3. **Échirolles:** independently verify LED share and September/October 2022 road switch-off schedule/coverage, municipal feeder power reduction and any spring/winter weather contamination. Test what fraction of the **FRK24 Isère** regional radiance trend plausibly reflects the commune and other independently switched-off communes, rather than equating a city with a region.
4. **Portsmouth:** retain historical **2019 pilot vs summer 2020 actual rollout** distinction; get activation schedules, historical LED inventory and identify unchanged placebo pilot neighbourhood before inferring incremental dimming. The 2020 four-panel spatial post-minus-pre map is predominantly positive, not support for a new citywide VIIRS decline.
5. **ISS:** use existing NUTS3 spectral work or source-date-specific local images selectively for **pre-2020 LED** contrast only; do not ask it to verify 2022/24 events.
6. **Write a compact evidence matrix** with columns: site, parent NUTS3, actual implementation and pre/post lamp technology, DNB change and pre-trend, independent actual switching/power reduction, applicable ISS colour change and image dates, strength/limits of inference for physical dimming, spatial share of NUTS3. Only once this is filled should we decide whether a brief supporting validation in the main paper is warranted.

## Key sources and previous work

- [VIIRS DNB limitation under HPS→LED](https://www.mdpi.com/2072-4292/6/12/11915) (2014 radiometry study); LED blue pump contributes light outside DNB sensitivity.
- [Sánchez de Miguel et al. (2022), ISS RGB European lighting colour shifts](https://pmc.ncbi.nlm.nih.gov/articles/PMC9473566/); mosaics' colour ratios and actual capture epochs.
- [Published ISS spectra vs locked VIIRS analysis](../../results_iss_spectral_vs_locked_viirs_2014_2024.md); especially the France example and independent streetlighting operations.
- [Five-case documentary register](./evidence_register_2014_2024.md).
- [Initial measured local VIIRS results](./five_case_initial_results_2026-11-09.md).
- [Portsmouth exact dated council audit and local spatial maps](./spatial_figures_and_portsmouth_date_audit_2026-11-10.md).

The **2017 VIIRS zero correction, 2021 snow handling, population-NUTS2016 cohort and spider/moth PRC curves remain locked**. For country or NUTS3 interpretation, report what could be due to spectrum, what independent evidence supports genuine output/hours changes, and what is not identified by available data.
