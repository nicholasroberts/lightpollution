# Five-case spatial maps: Bristol competing sources and Portsmouth intervention date audit

**Reviewed 10 November 2026; branch \`publication-figures\`.** Evidence comprises actual user-supplied outputs from \`src/local_viirs_spatial_dilution.py\`, previously reviewed city-level CSV summaries, and independently retrieved council/press documents. **Original PNG figures remain with the investigator, not stored in GitHub.** This is documentation, not a new calculation and not attribution of light flux to individual facilities. No changes to annual European VIIRS/NUTS3 or fixed PRC analyses.

## Bristol: land use and spatial cancellation

The investigator identifies the **northwestern bright zone** (approximately longitude −2.70 to −2.65, latitude 51.50–51.53) as **Avonmouth docks and industry**, and the small southern bright hotspot (near −2.62, 51.44) as **Ashton Gate Stadium surroundings** (football/rugby). Both are geographically plausible identifications, but there has been **no independent source attribution to particular industrial buildings, lamps, floodlights or road lights**. Native VIIRS DNB pixels blur and combine multiple emitters; Suomi NPP passes at approximately 01:30 local solar time, so scheduled match floodlighting is not an established explanation for a stadium-area pixel.

**Crucial distinction between maps:**

- The brightness-rank map reflects **pre-2023 reference-years' absolute brightness**, NOT a change. Both Avonmouth and Ashton Gate surroundings include bright reference pixels.
- The primary 2023 **post-start vs same-calendar-month 2020–22** seasonal anomaly map has positive patches in the northwest and a local southern bright hotspot, even while much of the wider city is grey/blue.
- In the **post-minus-pre seasonal anomaly** map, **northwest/Avonmouth becomes more positive**, whereas **Ashton Gate becomes strongly negative**: that southern location was **even more anomalously bright in Jan–Feb 2023** than across Apr–Dec 2023. It is therefore **incorrect to say both bright centres got brighter after the programme started**, based solely on red 2023 post-year pixels.
- The January–February vs April–December anomaly contrast still compares different calendar-month sets (each benchmarked to corresponding historical calendar months) and may be sensitive to different industrial use, sports schedules and atmospheric conditions; it is not a site-level instrument measuring a lighting retrofit.
- Earlier exact spatial CSV decomposition: **about 74%** of strictly paired Bristol cells have negative primary 2023 ΔR while matched municipal mean was about **+0.08 nW cm⁻² sr⁻¹**. Baseline **top 5%** cells had mean ΔR about **+5.32**, while the bottom half had mean ΔR about **−0.74** (refer to investigator’s \`brightness_strata_summary.csv\` for exact value and signed contributions). This is strong descriptive **positive/negative cancellation**, not proof of which streets received new LED/CMS lighting.
- The investigator's additional **secondary 2024** Bristol map and bar chart are qualitatively *different*: nearly all brightness bands decline, with an approximate mean **−10.4** within the highest 5% of baseline radiance and positive 2024 pockets mainly around the northwest. The 2024 aggregate is sensitive to exceptionally low \`cf_cvg\` in November; do NOT infer the apparent broad decline reflects only the municipal LED/CMS programme.

**Scientific implication:** Municipality-wide changes may be masked by **a few bright independent non-municipal emitters** and/or other localised sources, even if a majority of the mapped city darkens. A future robust test should spatially predeclare masks around Avonmouth industrial-port land and Ashton Gate, with independently georeferenced land-use footprints, separately from any **2023/24 historical streetlight retrofit cohorts**. Compare the municipality net with and without these masks *as observational sensitivity only*, not a data-chosen selection of "treated" pixels. Ask for asset/date/CMS records before interpreting a true programme exposure effect.

## Spatial context from other uploaded figures

- **Coventry primary Aug–Dec 2024**: signed baseline-class contributions are **negative across all four brightness strata**, with geographically widespread post-minus-pre darkening. Brightest pixels have greater absolute per-cell decreases, but the decline is not confined to a handful of highest-output pixels. The May–July citywide part-night implementation period remains a **window**, not an independently mapped lamp activation log.
- **Échirolles primary Oct–Dec 2022 and 2023 follow-up**: all four baseline brightness classes decline, highest-brightness cells lose more in absolute radiance, but **the bottom 50% of baseline-brightness cells also contribute substantially** to the net drop. The post-minus-pre primary map has strong broad negative values, especially northeast/central. The 2023 secondary cannot yield a clean within-year pre period and includes subsequent LED conversion.
- **Thanet primary 2017 and secondary 2018**: both brightness-band charts show the **highest 5%** of baseline cells driving the net positive aggregate, while many intermediate-brightness cells decline. The 2017 spatial map locates a **compact positive hotspot** near 1.28°E, 51.35°N; this has not been identified as a streetlight source. The two bright late-autumn months and spatial source concentration limit interpretation.
- **Portsmouth:** the supplied images here include the **2020 primary brightness-band decomposition** and **2021 secondary spatial map**, not the **2020 primary four-panel spatial anomaly map**. The primary brightness-band changes are negative in all four classes. The secondary 2021 map, which appropriately has no clean within-year pre panel, displays geographically mixed changes; it cannot adjudicate incremental **summer 2020** dimming. Earlier CSV pre/post means indicate Portsmouth's citywide anomaly was **already more negative in Jan–May 2020 than Sep–Dec**, so municipal mean alone does not identify a new 2020 operational effect.

## Portsmouth dates: planned vs independently confirmed historical

| Period/date | Type of evidence | What it actually demonstrates |
|---|---|---|
| **29 June / 11 July 2017** | Council official issue / decision history | LED replacement investment approved; **not commissioning** |
| **2018 through 2019** | Portsmouth council **Authority Monitoring Report 2018 (published Mar 2019)** | LED with CMS rollout underway in tertiary roads; **from Mar 2019** planned on primary/secondary roads; **14,500 by end Sep 2019** was *target*, not completion proof |
| **Late February 2019** | Contemporary 28 February 2019 newspaper interview | 2019 **adaptive-dimming trial had started the prior week** in five selected neighbourhoods, *before* citywide 2020 rollout; one trial area intentionally unchanged |
| **4 February 2020** | Signed Cabinet **record of decisions** | Citywide adaptive streetlighting **approved** following 2019 trial; not itself a switch-on date. Busy routes and some locations exempted |
| **Summer 2020** | Portsmouth City Council's **June 2022 retrospective operational guidance** | Explicitly states adaptive lighting **was rolled out throughout city in summer 2020** after 2019 trial; high confidence broad season, **not exact dates or lamp-by-lamp activation** |
| **March–April 2021 and beyond** | 21 April 2021 reporting of council programme | The final **~1,000 heritage/specialist LED upgrades** were then underway; thus not all LED conversions ended in 2019 |
| **22 March 2022** | Cabinet agenda/minutes listing an **Adaptive Lighting Review** | A review of the adaptive profile existed, but the actual report text/operational changes were **not verified** from the online search (archive report pack inaccessible through screenshot/download); do not assign a 2022 new effective date without further evidence |

**Original dated sources:**

- [Portsmouth 2017 Street Lighting Replacement Programme](https://portsmouth.moderngov.co.uk/mgIssueHistoryHome.aspx?IId=11501)
- [Portsmouth Authority Monitoring Report 2018 (Mar 2019), §3.7.5](https://portsmouth.moderngov.co.uk/documents/s22431/Housing%20Needs%20appendix%20-%20Authority%20Monitoring%20Report%20appendix.pdf)
- [28 Feb 2019 contemporaneous reporting on beginning of dimming trial](https://www.portsmouth.co.uk/news/politics/street-lights-in-portsmouth-set-to-be-dimmed-in-bid-to-save-cash-1310884)
- [Portsmouth Cabinet original decision record, 4 Feb 2020, item 9](https://portsmouth.moderngov.co.uk/documents/s26318/Record%20of%20decisions%20Cabinet%204%20Feb%202020.pdf)
- [5 Feb 2020 reporting on trial neighbourhoods, placebo and citywide plan](https://www.portsmouth.co.uk/news/politics/this-is-why-street-lights-across-portsmouth-will-be-dimmed-at-night-1384884)
- [**Portsmouth City Council, June 2022, published guidance** (p. 1: rollout *was* completed in summer 2020)](https://www.portsmouth.gov.uk/wp-content/uploads/2022/06/Guidance-for-dealing-with-stray-and-unwanted-light-from-Street-lights-0622.pdf)
- [21 Apr 2021 report of remaining heritage LEDs](https://www.hellorayo.co.uk/hits-radio/south-coast/news/concern-over-fate-of-historical-street-lights-in-portsmouth)
- [22 March 2022 Cabinet pack listing adaptive-profile review, pp.403–406](https://portsmouth.moderngov.co.uk/documents/g4931/Public%20reports%20pack%2022nd-Mar-2022%2013.30%20Cabinet.pdf?T=10)

**Critical exclusion:** Hampshire County Council's 1 April 2019 residential part-night switch-off was **not Portsmouth City Council's policy**; Portsmouth manages its own lighting. Do not misclassify a neighbouring authority's change as a Portsmouth intervention.

**Current event coding assessment:** Keep **summer 2020**, with June–August marked **transition/uncertain**, **September–December 2020** an exploratory post-rollout contrast and **January–May 2020** a pre-citywide reference but **not fully untreated**, due to the 2019 pilot and earlier LED conversion. Do **not** recode to 2021 merely to improve the timing match with the VIIRS decline. Without CMS log timestamps, there is still no evidence that all commanded dimming started on 1 September or finished by 31 August.

**Specific additional independent records required:** (i) citywide CMS schedule-history export by luminaire/column ID for **2019–2021**, (ii) activation/commissioning logs and profiles including the exact local 22:00–06:00 settings and exceptions, (iii) original 2019 trial treatment and **identity of placebo area**, (iv) district/street LED replacement completion dates; (v) actual March 2022 review/profile changes. Once obtained, join independent treated fraction to fixed VIIRS grid without selecting on observed changes.

## Next analysis gate

1. Obtain **\`portsmouth_primary_2020_spatial_anomalies.png\`** to inspect the spatial **post-minus-pre** map for the *actual proposed 2020 citywide intervention*. The **secondary_2021** figure cannot replace it.
2. For Bristol, define **geographic source-area sensitivity polygons from independent land-use maps**, with a reproducible radiance contribution decomposition and original citywide contrast side-by-side. Neither red post map pixel nor blue post-minus-pre pixel is direct proof of a particular facility changing lighting.
3. Keep the Coventry and Échirolles strong, broader spatial decreases as **controls for scale of aggregation**, but still add independent neighbouring-area and cloud/snow checks before causal inference.
4. Request Portsmouth CMS history as above before fitting event timing to the 2020 city mean. No automatic redefinition of dates or retrospective treated-pixel labelling.

**No updated numerical scripts are required just to record this evidence audit.** Preserve the original \`local_viirs_five_case_analysis.py\` intervention definitions while the documentary timeline is clarified.
