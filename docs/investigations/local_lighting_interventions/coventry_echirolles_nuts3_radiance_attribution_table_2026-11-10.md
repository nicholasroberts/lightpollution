# Coventry / Échirolles: dated operational changes versus locked NUTS3 VIIRS radiance

**Prepared 10 November 2026, branch \`publication-figures\`.** First table is evidence actually verified from municipal/council primary documents. Regional **2014 and 2024** numbers are from our previous frozen \`nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv\` (the Fig. 1/3 joint dataset), not from the new local-monthly V1 data. The complete **11-year annual NUTS3 trajectory has not been freshly exported** in this document, so no intermediate locked annual numbers are guessed or copied from older unapproved datasets. Read-only annual exporter committed separately. This is historical mechanism validation, not a new correction to annual VIIRS or species PRC.

## Dated operating and LED events, municipality scale

| Site and parent NUTS2016 | Date/period | Independent evidence about operations | LED context and limitation | Local V1 VIIRS observation |
|---|---|---|---|---|
| **Échirolles, FRK24 Isère** | **29 July 2022** | City reported actual extinguishing part of Rue Eugène Sue and its square, and planned August car-park switch-off near the Pathé | Site-specific switching predates widespread LED rollout, but not all lamps were sodium | Part of the **July–Sep 2022 transition**, excluded from clean pre period |
| Échirolles | **1 August 2022** | Municipal lights switched off **00:00–06:00 parks**, **23:00–06:00 municipal buildings and business parks**; park start adjusted to 22:00 in October | Multiple independent mechanisms, including building operation | Excludes Aug/Sep from pure untreated baseline |
| Échirolles | **26 September 2022 decision; October–December 2022 staged road rollout** | Road lights scheduled off **00:00–06:00**, with exemptions for key axes, tram stops and 24-hour healthcare access; official municipal periodical states roadway installation would proceed in phases **through December 2022** | 2025 retrospective municipal figure says ~**5%** LED stock in 2022. This is a **low share, not zero**; a separate 2023 report described a 2022 LED pilot in La Commanderie, so do not claim zero contemporaneous replacement | Oct–Dec 2022 mean **−37.08% vs Oct–Dec 2020/21**; strongest monthly declines **Nov–Dec**. Matched spatial maps show all baseline-brightness groups declining. |
| Échirolles | **June–summer 2023 onward** | Continued overnight switch-offs; later 2025 reduction in roadway hours to **01:00–05:00**, **outside** 2014–24 | City LED programme begins in earnest **June 2023**, with 120 heads replaced on Jean-Jaurès that summer; city says ~**40%** of stock LED by 2024; 2023/24 observations mix **operation and spectrum** | Oct–Dec 2023 secondary **−57.84% vs 2020/21**, **NOT a pure switch-off contrast** |
| **Coventry, UKG33** | **20 February 2024** | Council approved part-night switching at budget meeting | Approval is not actual lamp switch-off | Baseline year, no event assumed yet |
| Coventry | **May–July 2024** | **Retrospective council scrutiny minutes** confirm programming began May and **all eligible lamps were programmed by July**. About **60% of lights** affected after exemptions; city centre excluded. Weeknights **00:00–05:30**, weekends **01:00–05:30** | The Council's contemporaneous Feb 2024 Equality Impact Assessment says many residential lights were **already dimmed to ~55% output** and some Hillfields LED trial lights existed. Thus 2024 is a new **part-night off** intervention, not the first-ever brightness reduction and not a claim of zero pre-existing LEDs | **May–July transition**; post period starts **Aug 2024** |
| Coventry | **August–December 2024** | Part-night off operating across eligible footprint; 2025 council reports electricity-budget savings May–Dec 2024 approximately **£350k** at their earlier tariff, independent of VIIRS | The **major programme to replace up to 29,500 lights with LEDs was approved 14 January 2025 and began November 2025**, after the frozen study interval | Aug–Dec 2024 mean **−32.23% vs Aug–Dec 2022/23**; pre- vs post-seasonal anomaly difference **−3.91 nW cm⁻² sr⁻¹**; broad negative map |
| Coventry | **9 September / November 2025** | Previously off streetlights restored overnight **9 Sep 2025**; citywide LED conversions started **Nov 2025** | Both are **post-2024**, for documentary separation of mechanisms only | Do not use in locked 2014–24 VIIRS/PRC outcomes |

**Dated original evidence:**
- Échirolles [29 July 2022 rue Eugène Sue switch-off](https://www.echirolles.fr/actualites/extinction-de-leclairage-public-dans-la-rue-eugene-sue-et-dans-le-grand-parking-du); [21 July / 1 August 2022 first operations](https://www.echirolles.fr/actualites/un-nouveau-plan-de-sobriete-energetique-pour-la-ville); [23 Sep 2022 municipality/road-hours announcement](https://www.echirolles.fr/actualites/le-plan-de-sobriete-partage-avec-les-habitant-es); [city magazine No.398 Sep/Oct 2022, p.15: phase-by-phase roadway rollout until Dec](https://www.echirolles.fr/sites/default/files/2022-10/CITE%20ECHIROLLES%20SEPT%20OCT%2022.pdf); [city magazine No.404 Sep/Oct 2023, p.27: actual LED conversion and 2022 pilot](https://www.echirolles.fr/sites/default/files/2023-10/Cite%20Echirolles%20N404-web2_compresse.pdf); [municipal retrospective LED shares / 2025 hours](https://www.echirolles.fr/vie-quotidienne/cadre-de-vie/extinction-de-leclairage-public).
- Coventry [30 Jan 2025 scrutiny minutes: actual May–July 2024 completion, ~60%](https://edemocracy.coventry.gov.uk/mgAi.aspx?ID=43074); [Feb 2024 EIA: earlier ~55% dimming and Hillfields LED trial](https://edemocracy.coventry.gov.uk/documents/s59414/Appendix%2019%20-%20EIA%20-%20Street%20Lighting%20Switch%20Off.pdf); [9 Apr 2025 actual operating hours](https://edemocracy.coventry.gov.uk/ieListDocuments.aspx?CId=723&MId=13298); [city council 2025 reinstatement and main LED conversion chronology](https://www.coventry.gov.uk/street-scene/streetlighting).

## Existing frozen NUTS3 endpoint results, not local municipal monthly products

| Fixed NUTS2016 region | 2014 radiance | 2024 radiance | 2014→24 ΔR | ΔR % | Spider PRC 2014→24 | Moth PRC 2014→24 |
|---|---:|---:|---:|---:|---|---|
| **UKG33 Coventry** | **18.675611** | **13.030121** | **−5.645490** | **−30.229%** | 0.222393 → 0.255709 | 0.291554 → 0.358650 |
| **FRK24 Isère** | **2.074264** | **0.879914** | **−1.194349** | **−57.579%** | 0.417721 → 0.437239 | 0.627419 → 0.658224 |

Units radiance nW cm⁻² sr⁻¹; PRC dimensionless. Data: investigator's locked \`nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv\`, which is the source for the frozen Europe Fig. 1/3. These 2014→2024 differences are **NOT** 2023→2024 Coventry or 2021→2022 Isère event effects. The fixed NUTS2016 Isère polygon covers an extensive area including Grenoble and many other municipalities, whereas Échirolles is a small part.

## Proper bridge between the event and the regional time series

| Test from frozen annual parquet | Role | Inference limit |
|---|---|---|
| Coventry **2022→2023** | Pre-2024 regional annual reference change | Context only; other year-specific changes can occur |
| Coventry **2023→2024** | Regional year containing **May–July 2024 programming and Aug–Dec 2024 off period** | If negative, direction/timing compatible with genuine operational change. The 2024 annual V2 composite also contains untreated pre-May months and cannot partition LED/operation/meteorology alone |
| Isère **2021→2022** | Region covering the **July–Dec 2022 Échirolles staged switch-off** | 2021 requires the locked native-pixel snow correction; Échirolles forms only one municipal fraction of FRK24 |
| Isère **2022→2023** | First full post-year in Échirolles, but major **2023 LED conversions** start | Mixed spectra + continued operation; cannot attribute 2023 Isère decline to 2022 municipal switch-offs alone |
| Isère **2023→2024** | More conversion/operation in Échirolles and elsewhere | Changes in wider Isère may have completely different mechanisms |

**Critical source consistency:** The older uploaded \`nuts3_radiance_population_socioeconomic.csv\` shows a **2024 Isère value around 0.93575** whereas the final **locked endpoint is 0.879914**. The older file is not suitable as an unvalidated substitute for the current manuscript's final annual series, particularly in 2021 and 2024. Coventry is only partially represented there because UK population values are missing in later years. **Do not splice these older numbers into this table.**

**Complete read-only annual extraction script:** [\`src/export_coventry_echirolles_nuts3_bridge.py\`](../../../src/export_coventry_echirolles_nuts3_bridge.py). It reads the original frozen \`data/processed/fixed_nuts2016_viirs_2014_2024.parquet\`, validates 2014/2024 against the frozen Fig. 1 endpoints, and splices **only the existing audit-verified 2021 native snow-patch cache** as in publication Figure 1. It will **abort rather than regenerate** if the approved snow cache/audit is unavailable or stale, and it will never modify the source parquet or PRC curves.

Run on EcoVis:

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/export_coventry_echirolles_nuts3_bridge.py --self-test
python src/export_coventry_echirolles_nuts3_bridge.py
\`\`\`

It prints two region-by-year tables and writes:

- \`data/processed/local_lighting_interventions/nuts3_event_bridge/coventry_isere_locked_nuts3_annual_2014_2024.csv\`
- \`data/processed/local_lighting_interventions/nuts3_event_bridge/coventry_echirolles_nuts3_event_match.csv\`

**State: code committed, but real local frozen annual parquet and verified 2021 cache cannot be accessed directly through this chat. Numeric annual 2021/22/23/24 event-year differences are deliberately left pending EcoVis execution rather than made up.**

## Mechanistic interpretation

- **Coventry** is the strongest near-matched administrative bridge (LA approximately UKG33): real off-hours were independently documented before the principal 2025 LED replacement. But **prior 55% lamp dimming and some LEDs** mean this is not a virgin untreated 2023 baseline or proof that no other source changed. The event documents a real operating-hours reduction, and its matched monthly DNB drop is robust, while the *regional annual share* of this event remains to be quantified.
- **Échirolles** is an unusually informative near-low-LED real off-hours example in 2022, but the 2023→24 changes mix LED retrofit and operation. **FRK24** includes a much larger area: municipal evidence alone establishes that genuine reduction occurs **somewhere within Isère**, not what fraction of the **−57.6%** regional decadal fall was due to Échirolles.
- The point of this ancillary work is to reject an *all-spectral-only explanation* of historical changes where supported by independent operations, **not** to estimate a universal LED calibration correction or prove that every falling NUTS3 region reduced luminous output.
- No modifications to the locked European 2014–2024 radiance dataset, the 2017 zero correction, the 2021 snow overlay or the spider/moth PRC transfer functions.
