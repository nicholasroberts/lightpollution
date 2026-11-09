# Local lighting-intervention evidence register — 2014–2024 VIIRS / ISS

**Consolidated Stage A reports, five detailed case dossiers, source catalogue and data requests:** [Local lighting interventions hub](./local_lighting_interventions/README.md). The material is mirrored there to keep original links intact.\n\n**Status:** Evidence register and *proposed* extraction design, 9 November 2026. No local VIIRS pixels have yet been extracted or analysed in this stage. No local ISS mask or dated photographs have yet been validated. **The five registered cases are Échirolles, Coventry, Portsmouth, Thanet and Bristol.** The original first-four priority order is preserved; Bristol is added as a fifth, complementary mixed LED/dimming case.

## Scientific purpose and constraints

Test the *bounded* claim that observed declines in the VIIRS Day/Night Band (DNB) cannot universally be assigned to HPS→LED spectral sensitivity changes: some local changes should coincide with independently documented switch-offs or dimming. Conversely, seek an increase associated with restored all-night operation despite LED replacement. This is not a country-scale decomposition of optical flux or a correction to the locked biological PRC transfers.

The fixed **NUTS2016** region is a geographical container only; the treated units are municipal or streetlight-programme footprints **determined independently of VIIRS trends**. No fitting, source-product replacement, spectral correction or modification to the locked European population / spider / moth PRC analyses is authorised. New local outputs must be independent and labelled.

## Intervention evidence (pre-analysis)

| Order | Case | Original NUTS2016 container | Independently documented intervention and timing | Inventory / footprint available in public records | Proposed primary test |
|---|---|---|---|---|---|
| 1 | **Échirolles (France)** | **FRK24 — Isère** | Selected parks/buildings/industrial areas from **1 Aug 2022**; **roadway night switch-offs from 1 Oct 2022**, normally midnight–06:00, except major axes. Municipal LED modernisation began **2023**. A 2025 change to roadway 01:00–05:00 is **outside** 2014–24. | 4,900 municipal light points planned; LED share **5% in 2022**, **40% in 2024**. Exceptions: avenue du 8 Mai 1945, cours Jean-Jaurès, Honhoué/Grugliasco/Kimberley axis. Exact switchable circuit/lamp geocodes **not yet obtained**. | **2022 operational switch-off** before most LED conversions; separate 2023–24 mixed replacement period. |
| 2 | **Coventry (UK)** | **UKG33 — Coventry** | Budget decision **20 Feb 2024**; part-night programming **May–July 2024**, affecting approximately **60%** of lights after exclusions (city centre excluded); weekdays **00:00–05:30**, weekends **01:00–05:30**. The *main* LED replacement programme was approved January 2025 and commenced **November 2025**, outside the analysis interval. Streetlighting restored 9 Sep 2025, also outside. | Proposal covers **29,500** PFI streetlights for subsequent LED upgrades; 2024 affected lights ~60% excluding centre. Excluded lights near vulnerable facilities and other safety-sensitive roads, but no authoritative geocoded treated-versus-untreated inventory yet. | **2024 operational reduction** before the main 2025–27 LED retrofit. Do **not** assert all Coventry lamps were non-LED in 2024 without an inventory. |
| 3 | **Portsmouth (UK)** | **UKJ31 — Portsmouth** | Municipal LED replacement programme approved **2017**, conversions underway **2018**; **2019** adaptive lighting pilot; adaptive dimming **citywide summer 2020**. Contemporary report says most LEDs dimmed by around **50% from 22:00–06:00**, with exceptions. | ~**15,000** LED conversions reported by Feb 2020; heritage installations in Old Portsmouth and parts of Southsea not replaced. Precise conversion dates and CMS-dimming profiles by column remain to be obtained. | Compare earlier **2018–19 LED rollout** and subsequent **2020 operational dimming** as two events, rather than assigning both to LED sensitivity. |
| 4 | **Thanet (Kent, UK)** | **UKJ44 — East Kent** | Kent countywide LED/restored all-night lighting policy effective **15 Mar 2016**; Thanet phase 1 **Dec 2016–May 2017**; phase 2 planned from **Apr 2018**, with county programme targeted for **May 2019**. Verify actual completions, rather than treating projected end date as observed. | Thanet inventory reported **9,931** streetlights; **6,922** converted Dec 2016–May 2017; **3,009** yet to convert as of mid-2017, including ~2,500 concrete-column lights; planned second phase. Obtain converted-column list and separate actual date of all-night restoration. | **Candidate brightening** where longer operating hours counteract LED spectral effects. A positive VIIRS change is a hypothesis, **not an established result**. |
| 5 | **Bristol (UK)** | **UKK11 — Bristol, City of** | Citywide rolling LED and central-management-system (CMS) programme **began March 2023** (official council implementation webpage; earlier 2022 planning proposals). Work continued to **2025**, so **2014–24 captures only a partial rollout**. New LEDs are dimmed as road use declines at night; CMS can alter individual or group schedules. | Council's completed-action account: **29,366 lanterns plus nodes** and **4,500 nodes only** (reported after 2025). Council describes approximately **36,000** streetlights in total scope; manufacturer notes substantial LEDs were *already installed* in earlier years, so do **not** equate scope with HPS→LED conversions. Public **geocoded streetlight asset layer** contains unit IDs, positions, road names and `CMS_CONTROLLED`, but does **not** give historical replacement dates, lamp technology or actual operating schedules. Supplier reports **38 schedule variants** and **half-hourly per-lamp energy/activity telemetry**; access to historical 2023–24 logs unverified. | **Mixed spectral+operational intervention**: identify early vs late converted street-level cohorts and separate *lantern replacements* from *node-only CMS activation*; compare local DNB radiance against independently documented installation and midnight-time dimming/energy histories. Do **not** assume current CMS metadata represents 2024 history. |

### Evidence and URLs (primary sources preferred)

**Échirolles**
- Municipal July 2022 first switch-offs: https://www.echirolles.fr/actualites/un-nouveau-plan-de-sobriete-energetique-pour-la-ville
- Municipal notice specifying **1 October 2022** wider streetlighting plan: https://www.echirolles.fr/actualites/le-plan-de-sobriete-partage-avec-les-habitant-es
- Municipal inventory, 2022/2024 LED shares, current schedule and 2025 revision: https://www.echirolles.fr/vie-quotidienne/cadre-de-vie/extinction-de-leclairage-public
- July 2022 street-specific extinguishing near rue Eugène Sue: https://www.echirolles.fr/actualites/extinction-de-leclairage-public-dans-la-rue-eugene-sue-et-dans-le-grand-parking-du

**Coventry**
- Council scrutiny minutes, 2024 decision, May–July programming and 60%: https://edemocracy.coventry.gov.uk/mgAi.aspx?ID=43074
- Council April 2025 minutes, exact part-night hours and exclusions: https://edemocracy.coventry.gov.uk/ieListDocuments.aspx?CId=723&MId=13298
- Council October 2025 LED programme, scheduled Nov 2025–Nov 2027: https://edemocracy.coventry.gov.uk/mgAi.aspx?ID=44974
- Council current streetlighting ward timeline and 9 Sept 2025 reversal: https://www.coventry.gov.uk/street-scene/streetlighting

**Portsmouth**
- Council lighting replacement decision programme: https://portsmouth.moderngov.co.uk/mgIssueHistoryHome.aspx?IId=11501
- Council guidance documenting 2018 retrofit, 2019 pilot, **summer 2020** adaptive rollout: https://www.portsmouth.gov.uk/wp-content/uploads/2022/06/Guidance-for-dealing-with-stray-and-unwanted-light-from-Street-lights-0622.pdf
- Feb 2020 contemporaneous report, 15,000 lamps, reported 50% overnight dimming and pilot neighbourhoods: https://www.portsmouth.co.uk/news/politics/this-is-why-street-lights-across-portsmouth-will-be-dimmed-at-night-1384884
- Sussex/Kent/Hampshire regional dark-skies case study, descriptive map of 2015–19 change: https://storymaps.arcgis.com/stories/1f72ebb75695454691c094c0fd8dae13

**Thanet**
- County streetlighting LED + optimised all-night policy, decision 4 Mar / effective 15 Mar 2016: https://democracy.kent.gov.uk/ieDecisionDetails.aspx?Id=878
- Thanet conversion accounting, Jul 2017: https://theisleofthanetnews.com/2017/07/19/kcc-has-been-unable-to-convert-a-staggering-2500-streetlights-in-thanet-to-led-all-night-lighting/
- March 2018 report of planned April phase 2 and outstanding 3,009: https://theisleofthanetnews.com/2018/03/05/conversion-of-3009-thanet-streetlamps-to-led-all-night-lighting-to-begin-in-april/
- County Sept 2018 committee minutes confirm Thanet concrete-column delays still relevant: https://democracy.kent.gov.uk/ieListDocuments.aspx?CId=831&MID=7976
- Regional dark-skies case-study interpretation: https://storymaps.arcgis.com/stories/1f72ebb75695454691c094c0fd8dae13

**Bristol**
- City council LED upgrade project: https://www.bristol.gov.uk/residents/streets-travel/transport-plans-and-projects/upgrading-to-led-street-lighting (reports actual start **March 2023**, citywide scheme, CMS/dimming and scope; do not substitute an earlier *proposed* November 2022 start).
- City council Climate Emergency Action Plan 2026–2028, p. 11 and completed-action appendix p. 48: https://www.bristol.gov.uk/files/documents/5241-climate-emergency-action-plan/file (documents completed 2022–25 scheme, **29,366 lanterns/nodes and 4,500 nodes-only**, substantial electricity savings; this is **post-2024 contextual evidence**, not a 2024 completion timestamp).
- Schréder supplier account and interviews with Bristol council staff, **17 April 2026**: https://uk.schreder.com/en/insights/blogs/bristol-bridging-smart-city-gap (describes ~35,000 assets, mixed old/new LEDs, half-hourly per-asset energy and activity records, **38 schedule variants**, automated asset commissioning). Supplier-reported reductions in energy/CO2 are not optical-flux or VIIRS results.
- **Verified Bristol** public spatial dataset, publisher: https://opendata.bristol.gov.uk/datasets/bcc::streetlights-3 ; underlying layer https://maps2.bristol.gov.uk/server2/rest/services/ext/ll_transport/MapServer/37 . Fields include `UNITID`, `STREET`, `EASTING`, `NORTHING`, `UNIT_TYPE_DESCRIPTION`, `OWNER_DESCRIPTION`, `CMS_CONTROLLED`. The layer contains other street furniture too: filter to appropriate lighting types and council ownership, inspect duplicates, boundary and completeness; **no historical conversion date or 2024 dimming history is exposed**.
- Council street lighting team contact in its transport design guide: https://services.bristol.gov.uk/files/documents/4350-03-8-2-street-lighting (**lighting@bristol.gov.uk**). Request asset-ID–linked 2023–24 replacement dates, prior/new lamp type and watt/lumen specification, CMS activation dates and time-of-night dimming calendars, plus half-hourly usage/status or aggregated street-level operational histories for the VIIRS overpass period.
- **Inventory caveat:** Council counts and supplier descriptions mix replacement lanterns, existing LEDs fitted with nodes, and whole-asset scope; do not claim that 36,000 formerly HPS lamps were replaced, or that programme-wide emissions reductions demonstrate optical dimming.

**NUTS2016 provenance:** European Commission Regulation 2016/2066, Annex I (FRK24, UKG33, UKJ31, UKJ44, UKK11): https://eur-lex.europa.eu/legal-content/EN/TXT/?uri=CELEX:32016R2066

## VIIRS windows to extract (none have yet been measured in this local stage)

| Case | Annual comparison | Event-focused monthly/daily target | Key timing or control issue |
|---|---|---|---|
| Échirolles | 2019–21 pre; 2022 transition; 2023–24 after switch-off but with increasing LED conversions | Aug–Dec 2022 and into 2023, especially **Oct 2022** | Separate parks from roads and permanently lit major axes. **2021 needs exact approved snow-patch provenance** in native-pixel annual comparisons. |
| Coventry | 2022–23 baseline; **2024 partly treated** | Mar–Dec 2024, with May–July programming transition | Earlier existing monthly collection 2020–22 cannot measure the 2024 event: obtain 2023/24 data. Not all lights switched simultaneously. |
| Portsmouth | 2014–17 before conversion; 2018–19 conversion; 2020 pilot/rollout; 2021–24 after | 2017–21 if feasible, prioritising 2018/19 and **summer 2020** | Separate CMS-dimmed vs exempt and heritage lighting. Post-2020 years may include further adjustments. |
| Thanet | 2014–15 baseline; 2016–17 first retrofit; 2018–19 planned second phase; 2020–24 later period | 2016–19, aligned with confirmed per-column dates | County-wide programme contaminates neighbouring-town controls; select genuinely untreated lighting footprints, not merely adjacent pixels. |
| Bristol | 2014–22 context (some existing LEDs); **2023 start/transition, 2024 partial rollout** | **2022–24** monthlies, especially actual **2023–24 installation and CMS commissioning dates**, ideally daily/overpass-aligned | Within UKK11, compare predetermined **early-converted, late/not-yet-converted and nodes-only** clusters; use public asset IDs but seek *historical* conversion, wattage, lamp technology, CMS activation and dimming schedules. **2025 outcomes are outside the locked European period**, even though later records can validate intervention chronology. |

For all five: match valid fixed-grid DNB pixels to the documented *lighting footprint* and a defensible untreated or differentially treated local comparator within the original NUTS2016 region where possible. Report number of pixels, lighting density/coverage, total/mean/median DNB change and uncertainty; inspect pre-trends, baseline brightness and sensitivity to spatial buffer choice. Avoid pseudoreplication or treating repeated neighbouring DNB pixels as independent municipalities. Keep annual and monthly EOG product definitions, quality/coverage masks and zero-background correction consistent; do not subtract the locked **0.235 nW cm⁻² sr⁻¹** zero-offset twice. The frozen NUTS3 annual parquet's **2021 snow patch provenance is not fully verified**, so do not assume the patched 2021 is already included there. Use the approved native-pixel correction in any full-year comparison including 2021.

**Crucial methodological issue:** a local spatial average of thousands of streetlights can hide or dilute a genuine municipal event, and emission from neighbouring urban areas may contribute light outside an intervention polygon. County-wide policies (Thanet) also make geographically close untreated areas difficult to identify. Predeclare comparisons using streetlight inventory/administrative records, not peaks and troughs in the satellite series.

**Overpass alignment:** VIIRS DNB night-time observation is typically in the early morning, so all four policies have potentially observable nighttime schedules. Check actual acquisition times and local daylight-saving rules for each event, instead of assuming one fixed clock hour. If the switch-off occurs after overpass, a daily/annual composite cannot serve as direct evidence for that night's operation.

## ISS spectral component: what is known / unknown

- **Known product:** Existing calibrated European ISS colour mosaics from *2012–13* and irregular *2014–20* composites: https://zenodo.org/records/7677478. Existing repo pipeline: `src/validate_iss_spectral_shift.py`; broad comparison: `src/analyze_iss_spectral_vs_locked_viirs.py`, with methods and caveats in `docs/results_iss_spectral_vs_locked_viirs_2014_2024.md`.
- **Not yet established:** sufficient **common-valid ISS pixel support within each intervention footprint**; dates of the individual contributing images at each city, and hence whether a given location was photographed before and after the intervention.
- **Plan:** extract local **Δ(B/G)** and **Δ(G/R)** using exactly the established common-valid-pixel masking and calibration rules. Count matching pixels. Reconstruct image dates where possible. Do not assume the heterogeneous post mosaic cleanly observes Portsmouth 2018 vs 2020 separately, or that it captures Échirolles 2022, Coventry 2024 or Bristol 2023–24.
- **Potential later imagery:** NASA dated astronaut photos (https://eol.jsc.nasa.gov/SearchPhotos/) may help, conditional on adequate repeat scenes, registration and photometric calibration. SDGSAT-1 has later RGB colour data but its ~21:30 overpass cannot directly verify midnight switch-offs. **Availability for these five footprints is unverified.**
- No RGB ratio uniquely determines LED fraction, broadband emission, luminous flux, skyglow, spectral response correction or biological PRC change.

## Execution sequence and decision rules

1. **Stage A — independent footprints:** obtain lamp/circuit-level inventories and installation/operating records. Record true affected columns and dates, controls and exceptions; verify coordinates and NUTS2016 intersections. Do not select pixels on the basis of VIIRS outcomes.
2. **Stage B — local VIIRS:** read *original corrected native pixels* 2014–24 and relevant quality-controlled monthlies (acquire missing months only). Extract separate event and comparison series. Inspect and publish raw tables, masks and plots before fitting an event-study.
3. **Stage C — local ISS:** inspect common-valid pixel counts and image dates for all five. Quantify both colour ratios only where observed; explicitly flag temporally unusable mosaics.
4. **Stage D — triangulation:** separately report (a) DNB radiance outcome, (b) ISS spectral outcome, (c) documented engineering/operational intervention. Classify which associations are supported and which remain ambiguous. **Do not estimate unique LED-vs-operation flux percentages or change locked PRC functions.**

**First criterion for a convincing positive result:** a local VIIRS change temporally aligned with independently documented overnight switch-offs **before substantial local LED replacement** (especially Échirolles 2022, Coventry 2024), and absent or weaker in defensible controls. **For Portsmouth:** additional radiance change after 2020 adaptive dimming, distinguishable from 2018–19 LED transition. **For Thanet:** test, do not presume, an increase associated with restoration of all-night lighting. **For Bristol:** test 2023–24 staggered lantern replacement and CMS commissioning as distinct treatment classes, ideally matched to half-hourly operational profiles, but recognise that both spectrum and output may change together and that the rollout was incomplete at the 2024 endpoint.

## Required gaps before any causal interpretation

- Geocoded columns/circuits and actual policy implementation dates for the intervention and control footprints; for Bristol, public geocodes exist but **historical** asset-ID–linked lamp swaps and CMS schedules still need to be obtained.
- Product-consistent monthly/daily DNB around 2016–19, 2022–23 and **2024**, including Bristol's staggered 2023–24 rollout (current repo's 2020–22 monthly atlas is insufficient for most events).
- Native-pixel annual 2014–24 local availability and explicit 2021 correction provenance.
- ISS common-valid pixel counts and **exact contributing acquisition dates** for each locality; the existing post-2013 ISS mosaic cannot measure Bristol's 2023–24 intervention.
- Local lighting point density, lamp wattages/lumen outputs and actual operating schedules if available.
- A comparison design not contaminated by nearby rollouts or pre-existing illumination differences.

**No local effect or trend is asserted in this register; it is documentary evidence plus an analysis preregistration.**
