# Initial event-year anomaly plots — five independent lighting cases

**Observed and recorded:** 9 November 2026. Reviewed the five actual EcoVis PNG figures supplied by the investigator: `bristol_event_anomalies.png`, `coventry_event_anomalies.png`, `echirolles_event_anomalies.png`, `portsmouth_event_anomalies.png`, `thanet_event_anomalies.png`. Their original images are on EcoVis; **these visually inspected data are not the final exact numerical robustness outputs**. The robustness program generates three additional CSVs, but these **were not supplied in this review**, so below we give approximate visual readings and reserve precision for those CSVs.

## Shared interpretation

Each graph contrasts an outcome year's monthly **municipality mean VCMSLCFG EOG V1 radiance** with arithmetic means of the **same calendar month** in specified reference years. The bottom panel shows **median accepted cloud-free observations (`cf_cvg`) per valid pixel**, not cloud fraction. A plotted 12-month anomaly curve is descriptive, not a treatment-control or causal estimate. Outcomes and references span changing lighting spectra, weather, land cover and pre-existing trends; spatially correlated monthly VIIRS pixels are not independent exposure replications.

The earlier first-stage primary municipal mean contrasts were: Échirolles **−37.08%** (Oct–Dec 2022 vs 2020/21), Coventry **−32.23%** (Aug–Dec 2024 vs 2022/23), Portsmouth **−15.11%** (Sep–Dec 2020 vs 2018/19), Bristol **+0.02%** (Apr–Dec 2023 vs 2020–22), Thanet **+9.40%** (Jun–Nov 2017 vs 2015/16). Those are **whole-period seasonal comparisons**, not intervention increments.

## City-by-city diagnostic readings

| City | Before and during intervention | Observed post/late period | Interpretation |
|---|---|---|---|
| **Échirolles** | Apr–Jun 2022 near reference (approximately −1 to −2 nW cm⁻² sr⁻¹), July 2022 strongly positive, September approximately −7. Earlier site-level switch-offs in Jul–Aug precede the formal **1 Oct 2022** road switch-off. | October only modestly below reference; **November and December approximately −24** each. Median cf_cvg in these months is broadly near seasonal reference; August 2022 unusually high cf_cvg is outside the primary late-autumn contrast. | **Substantial and delayed late-autumn departure**, temporally compatible with staged switch-off. Do not call October an instantaneous measured radiance step, or use 2023 persistence as proof independent of 2023 LED conversions. Weather/snow and nearby atmospheric-context municipalities remain essential checks. |
| **Coventry** | Jan–Apr 2024 anomalies mostly near reference (approximately −4 to +1), then May–Jul 2024 rollout window shaded. | **All Aug–Dec months negative**, approximately −4 to −7 nW cm⁻² sr⁻¹, with November especially negative. Outcome cf_cvg is lower in several months but adequate to test with the already coded sensitivities. | **Strongest temporal contrast** of the five, complementary to widespread negative 2024 spatial maps. Confirm early-year vs post-month anomaly differential, pre-existing 2022/23 trends and low-cf sensitivity. The main LED rollout occurred later in 2025. |
| **Portsmouth** | Jan–May 2020 anomaly already **strongly negative** versus 2018/19 (approximately −2 to −7 across months), before summer 2020 adaptive dimming. | Sep–Dec 2020 still negative but **not appreciably more negative** than pre-summer; visual mean of post anomalies is probably *less negative* than pre. | The whole-period **−15.11% is NOT evidence for an incremental summer-2020 dimming effect**. Pre-existing 2018/19 LED conversion and year-to-year conditions complicate it. Strong candidate for a non-identifiable/attribution-limited contrast, **not** proof that dimming had no effect. |
| **Bristol** | Jan–Feb 2023 are the only fully pre-programme months, with modest mixed anomalies. The lighting LED + CMS programme began **Mar 2023 and continued throughout 2023/24**, not merely March. | Apr–Dec 2023 anomalies cross zero repeatedly, with sizable negative midyear and a strong **positive December anomaly (~+5.5)**. Monthly median cf_cvg changes substantially. | The 2023 **+0.02%** primary aggregate conceals strong month-level swings and heterogeneous spatial areas. This is a **partly exposed rolling intervention**, not a clean before/after experiment. The plotting script previously shaded March only, erroneously making the rest of 2023 look fully post-transition; **updated to shade March–December 2023**. |
| **Thanet** | First installation interval is **Dec 2016–May 2017**; no clean within-2017 pre-intervention period exists. Early 2017 January–March differences are mixed/positive, April–May somewhat negative. | Jun–Sep 2017 mostly neutral/negative; **October and November ~+7** each; December negative. The initial 2017 spatial map has a **single compact high-change hotspot** near 1.28 E, 51.35 N. | Aggregate **+9.40%** is influenced by two late-autumn months and a spatially concentrated source. Need hotspot identification / secondary 2018 map and leave-one-month-out tests before attributing brightening to LEDs plus all-night operation. |

**Important caveats:** Values marked approximate were read off plotted lines. Exact primary and robustness results should come from the existing CSV exports, not hand-digitised images. The images alone cannot calculate observation-night support, reliable uncertainty or counterfactual trends.

## Consequences for the study

- Échirolles and Coventry retain the clearest *candidate* associations between independent municipal night-lighting operation changes and reduced VIIRS DNB observations. Their complementary operation-focused policies make them valuable, but **do not yet establish causality**.
- Portsmouth's before-versus-after anomaly pattern is the major counterexample to naive interpretation of a negative season-matched mean contrast; much of the fall had already occurred *before* its putative dimming event. Do not present the primary −15.1% as an estimated dimming effect.
- Bristol shows the limitation of describing an ongoing geographically staggered LED/CMS transformation as an instantaneous event.
- Thanet's pattern illustrates why **typical-cell ΔR** and urban-source/land-use identification matter alongside municipality means.
- Low `cf_cvg` can affect the selection of nights in a monthly composite, even when almost 100% of raster pixel centres have at least one valid observation. It is not itself a snow/atmospheric correction.
- A negative VIIRS radiance difference is not a direct quantitative reduction in lamp luminous flux or sky polarization cues. The frozen European and PRC calculations are not refitted or corrected here.

## Next exact numerical gate

The existing `src/local_viirs_five_case_robustness.py` writes three machine-readable tables. **Review all three before constructing new causal models**:

1. `pre_post_anomaly_differentials.csv` — exact pre-event anomaly, post-event anomaly and post-minus-pre contrast for Échirolles, Coventry, Portsmouth and Bristol; Thanet deliberately reports unavailable pre-event.
2. `earlier_year_placebo_comparisons.csv` — earlier-year same-month shifts, to assess pre-existing city trends and whether similarly sized changes occurred *without* the chosen event.
3. `reference_and_cloudfree_sensitivity.csv` — single-reference-year, cf_cvg>=3/5 spatial-support variations, median cloud-free month filters, and leave-one-month-out results. The robustness script now **explicitly evaluates Bristol April–December 2024** as a separately labelled **secondary_2024** scenario, including the cloud-free thresholds and a complete **leave-one-month-out check**, with **November 2024 omitted** in one row. The primary Bristol scenario remains 2023.

Additional later analyses: fixed external atmospheric-context area comparisons, winter 2022/23 Échirolles surface/snow checks, geographic/lighting footprints for Coventry and the Thanet brightening hotspot, and targeted Portsmouth subcity treatment geography if full municipality remains inconclusive.

## Code and EcoVis commands

The updated [robustness script](../../../src/local_viirs_five_case_robustness.py) has been corrected to **shade the full March–December 2023 ongoing Bristol rollout** rather than just its starting month and to **analyse Bristol 2024 as an explicitly labelled secondary follow-up with leave-November-out sensitivity**. The original primary five-case contrasts are unchanged.

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/local_viirs_five_case_robustness.py --self-test
python src/local_viirs_five_case_robustness.py
```

The outputs remain in `data/processed/local_lighting_interventions/five_case_robustness/` and `figures/local_lighting_interventions/five_case_robustness/`. No re-extraction of monthly EOG archive files, no alteration to European annual processing or PRC model.
