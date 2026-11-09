# Portsmouth and Échirolles — spatial VIIRS monthly pilot: observed maps

**Review date:** 9 November 2026. **Status: quantitative/visual inspection of six actual EcoVis map PNGs, descriptive only.** User-provided map images (original full-resolution PNGs supplied in conversation; not copied into this GitHub repository): `echirolles_vcmslcfg_matched_month_difference_maps.png`, `echirolles_vcmslcfg_period_cloudfree_maps.png`, `echirolles_vcmslcfg_period_mean_maps.png`, `portsmouth_vcmslcfg_matched_month_difference_maps.png`, `portsmouth_vcmslcfg_period_cloudfree_maps.png`, `portsmouth_vcmslcfg_period_mean_maps.png`. Monthly aggregate figures below are from the user's previously supplied `combined_monthly_radiance_cloudfree.csv`, which has **2021–2022** for Échirolles and **2020–2021** for Portsmouth. The existing script is now configured to include 2020 for Échirolles, but the **2020 results were not included in the CSV reviewed here** and must not be assumed.

## Scientific framing

The bounded aim is to test whether regional VIIRS declines always indicate spectral replacement (HPS→LED) or whether some changes are consistent with independently documented operation changes. Échirolles had an October **2022** street-lighting switch-off, preceded by some site-level earlier interventions, with broad LED rollout increasing later. Portsmouth had an earlier LED retrofit followed by summer **2020** adaptive dimming. The existence of an observed local DNB decline and the documented intervention **remain separate observations** until controls, weather/snow and actual lamp footprints are considered. No locked European, zero-point, snow or biological PRC datasets are changed.

## Échirolles — 54 DNB gridded pixel centres

**Main finding: coherent local decline, strongest in previously brightest central/eastern part of municipality.** The Oct–Dec **2022 minus Oct–Dec 2021** native-grid mean radiance map is dominated by negative values, rather than one disappearing hotspot. The shared-colour-range maps show the pre-existing high-radiance core is markedly dimmer in late 2022. Adjacent grid cells are spatially correlated; do **not** claim 54 independent measurements of a lighting intervention.

| Same-calendar month | Citywide mean ΔR (2022−2021), nW cm⁻² sr⁻¹ | Relative mean change | Spatial median of 54 paired-cell ΔR read from maps |
|---|---:|---:|---:|
| September | −4.18 | −10.1% | −2.99 |
| October | −9.82 | −20.3% | −8.76 |
| November | −18.34 | −41.8% | −16.80 |
| December | −23.98 | −56.4% | −20.18 |

The October–December pooled monthly mean was **38.7% lower** in 2022 than in the corresponding months of 2021; December 2022 municipal spatial median radiance was **19.26** vs **36.47** in December 2021 (−47.2%). The spatial median of paired-pixel **differences**, e.g. December −20.18, is **not** the difference of spatial medians (19.26−36.47=−17.20).

**Temporal nuance:** September 2022 was already darker than September 2021. October 2022 mean radiance (38.55) was *slightly above* September 2022 (37.02); the October decline is relative to the **same-calendar-month 2021 baseline**, not an abrupt within-2022 step. The strongest 2022 decline developed in November and December. This can be consistent with an October implementation/switch-off but also deserves scrutiny for weather, acquisition timing, evolving operation and snow/reflectance.

**Coverage:** all 54 Échirolles cells had `cf_cvg>0` for the relevant months in the VCMSLCFG aggregate. The two period coverage maps are broadly comparable, with local variation but no large patch of missing coverage matching the central/eastern radiance decline. Median per-pixel cloud-free counts for September–December 2021: **11, 13, 7, 5**; for 2022: **10, 11, 7, 6**, respectively. Similar **counts do not** imply the same observed nights, cloud optical properties, visibility, aerosol conditions or snow cover. 2021's early-year anomaly was documented in the broader project; no equivalent 2021/22 late-year local albedo comparison has yet been run.

**Interpretation:** This is stronger evidence for a **real, spatially extensive change in observed VIIRS DNB radiance** that coincides with a independently documented municipal lighting intervention. It is a **promising operation-driven candidate**, but is **not yet proof** that a given share of total emitted optical flux/skyglow was reduced or that spectral change is irrelevant.

## Portsmouth — UKJ31, 297 DNB gridded pixel centres

- Citywide **Jan–May 2020 vs Sep–Dec 2020** native-grid radiance maps retain a broadly similar bright-core structure, with a mixture of local positive and negative differences. This **within-year** comparison is **seasonally confounded**, so no interpretation as a dimming treatment effect should be made.
- Same-calendar September–December **2021 minus 2020** maps are spatially mixed: paired-cell median differences approximately **−0.45 (Sep), +1.28 (Oct), −0.54 (Nov), −4.69 (Dec)** nW cm⁻² sr⁻¹. December 2021 is notably more negative in much of Portsmouth than September–November, but **both years are post-summer-2020**, making this **post/post**, not a valid event discontinuity test.
- The early-vs-late-2020 coverage maps show appreciably **fewer clear-sky observations per pixel in Sep–Dec**. Unequal observation frequency is not itself a radiance correction and may signal seasonal sampling differences.
- The existing CSV records VCMSLCFG valid footprint near 100% most months; good valid-cell coverage is not the same as identical clear-sky conditions or observation dates.

**Interpretation:** No clean large citywide 2020 dimming change is yet detected on these maps; **Portsmouth is inconclusive**, not a null effect. Need pre-2020 monthly reference or better predetermined subcity treatment groups to evaluate staged adaptive dimming.

## Immediate comparison: what the maps add

1. Échirolles' post-October radiance decline is **not confined to a single cell** and becomes more spatially intense from October to December; its apparent concentration in the city's previously bright core is plausible for changes in built-up lighting but does not localise particular lamp circuits.
2. A major missing-data/cloud-free-count collapse does **not** coincide with the Échirolles darkening, which argues against a trivial missing-cell artifact. It does **not** control snow, aerosols, visibility or differing observed nights.
3. Portsmouth's spatial changes are more heterogeneous and its observation frequency more uneven across the seasons; available maps do not isolate the scheduled dimming mechanism.
4. No per-column asset records are required yet to establish these patterns. They become useful for matching affected areas or distinguishing dimming from replacement.

## The next two tests — minimal additional work

**1. Cross-check Échirolles against the already available 2020 same-calendar-month data.** Re-run the updated time series, then compare October–December **2022** to **both 2020 and 2021**. The maps here have a **2021-only** reference and the uploaded historical CSV contains only 2021–2022 Échirolles. Avoid declaring 2022 uniquely anomalous without the extra 2020 comparator.

**2. Compare Échirolles against an independently documented neighbouring municipality or urban area under a different nighttime-lighting schedule.** Similar weather, snow, aerosol and regional scene-level issues could create shared declines. The control must be selected based on lighting policy or geography, not its VIIRS trend; nearby Grenoble is not automatically unaffected. If no clean local control exists, use a regional comparison as **atmospheric context**, not an untreated causal control.

If neither reveals strong counterevidence, targeted scene acquisition time/snow checks and request for the limited municipal switched-road footprint can strengthen the interpretation. No need for another pan-European model, a flux partition, or changes to the PRC functions.

## Execution and storage

Script: [`src/local_viirs_monthly_pilot.py`](../../../src/local_viirs_monthly_pilot.py). Inputs are **local EcoVis EOG monthly GeoTIFFs**. Native-grid PNG/PDF figure locations on EcoVis: `figures/local_lighting_interventions/monthly_pilot/`; numeric local GeoTIFFs and tables: `data/processed/local_lighting_interventions/monthly_pilot/`.

To refresh:
```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/local_viirs_monthly_pilot.py
```
If Git refuses because of local modifications, inspect `git status --short` and preserve any untracked work; never force a checkout. This report is observational documentation, not an independent archive of raw raster files or original PNGs.
