# First monthly local VIIRS pilot — observations from EcoVis exports

**Date:** 9 Nov 2026. **Status: first descriptive inspection**, NOT an event effect estimate and NOT a flux or spectral correction. Source: user's EcoVis outputs from `src/local_viirs_monthly_pilot.py`: `combined_monthly_radiance_cloudfree.csv`, `echirolles_monthly_radiance_and_cloudfree.png`, `echirolles_same_calendar_month_comparison.png`, `portsmouth_monthly_radiance_and_cloudfree.png`. Original CSV/PNG files were supplied in conversation; this note records key computations, does not claim the original binaries are stored in GitHub. Frozen European datasets and PRC functions not changed.

## Observed monthly dataset

**96 records**, two cities × two monthly EOG VIIRS V1 products × 24 actual months each.
- **Échirolles (INSEE municipality 38151):** 2021-01 through 2022-12 (24 observed months/product), **54 polygon grid cells**. Each VCMSLCFG month has valid radiance and `cf_cvg>0` for **54/54 cells**. Note: spatial grid support ≠ 54 independent native 750-m instrument footprints.
- **Portsmouth (NUTS2016 UKJ31):** 2020-01 through 2021-12, **297 polygon grid cells**, with VCMSLCFG monthly valid area generally 297/297 (lowest valid fraction 98.32%).
- The *earlier terminal inventory* showing "2023 January" **appears a false positive**: it searched for `202301` anywhere in a filename, including EOG *creation timestamp* `_c202301...` attached to **December 2022** radiance and cf_cvg files. The actual extraction based on *observation start/end* date ranges has **no 2023 observation** in the exported CSV. This should be verified by listing files rather than assuming a 2023 January product exists.

## Échirolles: 2022 minus same calendar months 2021

**Primary VCMSLCFG product: mean radiance nW cm^-2 sr^-1.**

| Period/month | 2021 | 2022 | 2022 minus 2021 | Relative 2022 change |
|---|---:|---:|---:|---:|
| Jan–Aug monthly mean | 48.99 | 52.83 | +3.84 | +7.8% |
| Sep | 41.20 | 37.02 | -4.18 | -10.1% |
| Oct | 48.37 | 38.55 | -9.82 | -20.3% |
| Nov | 43.86 | 25.52 | -18.34 | -41.8% |
| Dec | 42.55 | 18.57 | -23.98 | -56.4% |
| Oct–Dec monthly mean | 44.92 | 27.55 | -17.38 | -38.7% |

- The 2022 late-year decline coincides approximately with the documented **1 October 2022** wider roadway switch-off and independently described October 2022 municipal power decrease; it is **not an instantaneous October step**. September (before the roadway event, but after some July/August site-level switch-offs) was already below September 2021. October 2022 is **higher** than September 2022 in absolute radiance; the October intervention association appears in the **same-calendar-year contrast**, not an October-only raw discontinuity.
- At December, the **spatial median** of the 54 gridded radiances also declines, **36.47 (2021) to 19.26 (2022), ~-47.2%**, so this is not solely a mean driven by a few bright pixels.
- Both VCMCFG and VCMSLCFG recover essentially identical September–December 2022 magnitudes, including October–December exactly equal where scenes overlap; no evidence that merely choosing stray-light-corrected product creates this step.
- *Coverage check:* `cf_cvg` median count per valid pixel for Sep,Oct,Nov,Dec **2021** is **11,13,7,5**, and **2022** is **10,11,7,6**, with all 54 VCMSLCFG grid cells covered in both years. Hence no evidence for a simple missing-cell or collapsed-count artifact. BUT observed-night selection, weather, snow/reflectance and imperfect instrument/atmospheric correction remain untested. The **2021 winter snow anomaly** in the larger European work was mainly early 2021; do not automatically infer late-2021 reference months were snow free.
- **2022 July caution:** VCMSLCFG median only **2** cloud-free observations/pixel; interpret its apparent high radiance carefully. Stray-light-filtered VCMCFG is absent for several summer months at this latitude.
- This is **a compelling preliminary coincidence, not proof** the whole DNB drop represents operational reduction or unique evidence to assign an LED/flux fraction.

## Portsmouth: no persuasive city-scale breakpoint yet

**Primary VCMSLCFG:** 2020 vs same 2021 months.

| Period | 2020 | 2021 | 2021 minus 2020 |
|---|---:|---:|---:|
| Jan–Sep monthly mean | 22.15 | 21.53 | -2.8% |
| Oct–Dec monthly mean | 20.14 | 19.07 | -5.3% |

- Within **2020**, Jan–Jul average **22.56**, Aug–Dec average **20.38** (about **-9.7%**). Within **2021**, Jan–Jul **22.02** versus Aug–Dec **19.37** (**-12.0%**). This indicates an end-of-year reduction is not unique to the **summer 2020** policy change; however 2021 was post-intervention and is **not a no-treatment control**.
- The lack of a clean sustained citywide break may reflect seasonal variance, smaller effective change, mixed converted/exempt streets, and dilution in a city-wide average. The missing 2019 monthly reference is limiting. Do **not** call Portsmouth's dimming ineffective based on aggregate DNB.
- VCMSLCFG coverage remains near complete; VCMCFG has summer missing months (zero `cf_cvg`), so use VCMSLCFG as initial source and VCMCFG for months where it genuinely has observations.

## Immediate validation with existing inputs, before new documentary searches

1. **Add 2020 for Échirolles.** It is present on EcoVis, yet original script restricted this municipality to 2021 onward. Plot same calendar months in 2020, 2021 and 2022; determine whether late-2022 is exceptional relative to **both** reference years and avoid relying entirely on 2021.
2. **One nearby-comparison trajectory:** select a suitable, independently documented comparable municipality without the same switch-off, or treat Grenoble/other nearby municipality as descriptive geographic context only if its independent treatment status is unknown. Do not assume metropolitan neighbours were untreated.
3. For Échirolles, first verify the fall is spatially widespread by mapping 2022-minus-matched-2021 **per-grid-cell** contrasts in Sep/Oct/Nov/Dec, with cf_cvg sufficiency. This helps identify an actual municipal signal versus one bright hotspot, not 54 independent replications.
4. Snow/surface reflectance and observation timing must be checked before operational attribution. There is no need yet for per-lamp geocoded records as a prerequisite for these **city-level** plots.

## Discussion boundary

The immediate proposition is **limited**: a strong local VIIRS decline is temporally aligned with an independently documented municipal *switch-off* that preceded large-scale LED conversion. Optical/luminous-flux reductions cannot be quantified from DNB and energy evidence alone. This evidence **supports further falsification tests**, but the 12-month data alone do not rule out correlated weather/snow or another contemporaneous change. The existing locked 2014–24 European and spider/moth PRC analyses remain untouched.
