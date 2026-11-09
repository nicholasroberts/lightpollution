# Five-city VIIRS: spatial dilution analysis (independent Stage B extension)

**Prepared 9 November 2026** on GitHub `publication-figures`. The script is committed; **native-grid spatial outputs still have to be calculated on EcoVis**. This document records the study design and the numerical robustness CSV results **actually supplied by the investigator**. No data or figures in this report are presented as results of the new spatial script.

## Question

Does geographic aggregation across an entire municipality conceal a radiance change in smaller subareas, especially in Portsmouth and Bristol? More broadly, do radiance changes follow particular **a priori brightness classes** or **fixed geographical subareas**, and do positive and negative subarea changes cancel in citywide means?

**Important boundary:** With no independently dated lamp/circuit/road polygons, we **cannot** designate pixels with strong negative VIIRS changes as affected streetlight pixels. That would select "treatment" on the outcome and exaggerate any apparent effect. Brightness bands and bounding-box quadrants below are descriptive stratifications, **not estimates of illumination treatment exposure**. Actual spatial dilution by treated lamps can only be established with independently verified treatment geometry.

## Five case-control dates and baseline support

The program reuses the independent lighting dates and boundary catalogue from [the five-case Stage B protocol](./full_five_case_monthly_analysis_protocol_2026-11-09.md):

| Case | Primary outcome months | Same-month reference | Clean pre-event outcome months |
|---|---|---|---|
| Échirolles | October–December **2022** | **2020, 2021** | April–June 2022; July–September already includes some early site-level switching |
| Portsmouth | September–December **2020** | **2018, 2019** | January–May 2020; 2018/19 reference itself overlaps LED transition |
| Coventry | August–December **2024** | **2022, 2023** | January–April 2024; May–July is rollout |
| Bristol | April–December **2023** | **2020–2022** | January–February 2023 only; March–December is ongoing rollout |
| Thanet | June–November **2017** | **2015, 2016** | **None** during 2017; installation ran December 2016–May 2017 |

Also show secondary **Échirolles 2023, Portsmouth 2021, Bristol 2024, Thanet 2018** for persistence/mixed-intervention context, **without claiming an untreated pre period** for those already-treated/ongoing years.

## Before new spatial extraction: exact uploaded robustness findings

These numbers are from user-uploaded `pre_post_anomaly_differentials.csv`, `earlier_year_placebo_comparisons.csv`, and `reference_and_cloudfree_sensitivity.csv`. These are **city-month mean radiance comparisons**, **not spatial decomposition**.

| Place | Pre-event seasonal anomaly | Late/event-season anomaly | Post − pre anomaly | Important robustness |
|---|---:|---:|---:|---|
| Échirolles 2022 | −1.713 | −16.234 | **−14.521** | At >=5 cf observations/pixel −36.20% versus −37.08% all; 2020-only −35.39%, 2021-only −38.68% |
| Coventry 2024 | −1.379 | −5.290 | **−3.911** | At >=5 cf obs/pixel −32.11% vs −32.23%; dropping any single outcome month retains −28.82 to −35.08% |
| Portsmouth 2020 | −5.507 | −3.645 | **+1.863** | At >=5 cf obs/pixel −14.29% vs −15.11%; **2019 vs 2017/18 seasonal placebo −25.87%**, larger than primary 2020 contrast |
| Bristol 2023 | −0.366 | +0.003 | **+0.369** | All months +0.02%; at >=5 cf obs/pixel +2.37%; 2023 mixed changes with rollout ongoing |
| Thanet 2017 | n/a | +1.005 | n/a | Omit October **or** November: −1.84% or −1.81% respectively, versus +9.40% all; brightening dominated by those months |

**Secondary Bristol 2024 (already present in robustness CSV):** Apr–Dec 2024 vs 2020–22 is **−18.20%**, but removing **November** reduces the decline to **−11.17%** (8 months), and retaining only months with outcome/reference median cf_cvg ≥5 yields **−12.25%** (5 months). Thus it is a real descriptive signal but unusually sensitive to sparse-observation November. The cf_cvg>=5 *per-pixel* threshold across the original months produces −27.64% (8 months) because the spatial sample changes; **do not describe this as a correction or evidence the citywide fall was caused by the intervention**.

**Earlier-year placebos:** Échirolles Oct–Dec 2021 vs 2019/20 **+6.80%**; Oct–Dec 2020 vs 2018/19 **−1.33%**. Coventry Aug–Dec 2023 vs 2021/22 **+15.71%**; Aug–Dec 2022 vs 2020/21 **+3.32%**. Portsmouth Sep–Dec 2019 vs 2017/18 **−25.87%**, indicating a large decline predating 2020 dimming. Thanet Jun–Nov 2016 vs 2015 **−18.92%**.

## New predeclared spatial design

Standalone script: [`src/local_viirs_spatial_dilution.py`](../../../src/local_viirs_spatial_dilution.py).

**Input:** Original locally extracted EOG **monthly V1 VCMSLCFG 75N060W** `avg_rade9h` and `cf_cvg` TIFFs, plus existing official/cached **municipal boundaries**. Uses same case references and dates as full five-case analysis; does NOT download or replace EOG data.

- **Reference-year brightness ranking, fixed before outcome:** For each city, compute each native VIIRS cell's mean baseline radiance using only its designated **reference years**, across a prespecified set of event comparison and clean pre-event calendar months. Bin eligible cells into **bottom 50%, 50–80%, 80–95%, highest 5%** by ascending pre-event brightness. The top 5% captures concentrated sources without classifying any city cell by *its later decline*. These groups are not lamp-treated cohorts; noisy baseline ranks also induce regression-to-mean/mathematical coupling.
- **Fixed spatial quadrants:** Assign each pixel centre to **NW, NE, SW, SE**, based only on the municipality's fixed bounding-box midlines. This can identify broad spatial heterogeneity even if area is not bright; it is an arbitrary geographical split, not streetlight treatment geometry.
- **Strictly matched monthly pixels:** For each outcome month, a cell must have valid VIIRS radiance and `cf_cvg >= N` **in the outcome and every referenced earlier year for that same calendar month**. `N=1,3,5` are sensitivity settings. A cell enters a period-average spatial map only if at least **75% of months** have such strict paired support. No missing values are replaced by zero; archive creation timestamps or **75N060E** wrong tiles are not accepted.
- **Spatial anomaly timing:** For the *primary* cases with a genuine pre-event interval, compute pre-event and post-event same-month anomaly **on each common pixel**; map **post minus pre anomaly**. The clean pre-event map is deliberately unavailable for Thanet and for all secondary follow-ups. These are observational changes, not causal differences-in-differences with an external untreated control.
- **Additive citywide decomposition:** Report for each brightness/sector group the **number/fraction of paired pixels**, baseline mean and total radiance share, absolute/percent mean change, signed contribution to the **overall matched-pixel municipal mean** (sum over group ΔR divided by the number of citywide paired cells), proportion of cells declining, and separate shares of gross negative and positive contributions. Group signed contributions **must sum** to matched-pixel municipality mean. Comparing these sums reveals cancellation/dilution **descriptively**, not that an intervention affected particular pixels. Pixels in highly luminous groups can contribute disproportionately even if few.
- **Common pre/post support:** For the timing differential within a group, include **only pixels with valid pre and post anomalies**. The n of common cells is reported separately from the n in the post-season spatial sum. Do not compare estimates if cell support changes without reporting it.
- **Autocorrelation and sensor footprint:** A gridded 15 arc-second cell is not an independent 15 arc-second optical sensor footprint; hundreds of adjacent cells have correlated light response, atmospheric conditions and sensor blur. Do not use naive pixel-level t-tests or infer number of independently switched streets from pixel counts.
- **No LED/flux partition:** VIIRS DNB brightness, city electricity use and luminous flux are not interchangeable, and EOG V1 does not by itself distinguish LED spectral effects from operating-hour changes. These spatial results will identify where **measured radiance** changed, not prove where lamp output changed.

## Outputs

Run separately from the existing `five_case_monthly` and `five_case_robustness` folders. New CSV directory:

`data/processed/local_lighting_interventions/spatial_dilution/`

- `spatial_stratified_contrasts.csv`: primary/secondary cases × baseline brightness / geographic quadrant × `cf_cvg>=1,3,5` group tables, with exact additive net contributions, fractions negative and common pre/post anomalies where meaningful.
- `brightness_strata_summary.csv`: subset for baseline brightness classes, `cf_cvg>=1`, convenient to inspect first.
- `baseline_spatial_strata_metadata.csv`: independent reference-year and baseline-calendar-month definitions; source official boundary.

PNG/PDF figure directory:

`figures/local_lighting_interventions/spatial_dilution/`

- `<city>_<primary|secondary>_<year>_spatial_anomalies`: **2×2** pre-event brightness rank, primary pre-event anomaly, post-event anomaly, post-minus-pre cell-level anomaly differential. For Thanet and secondary cases the pre/post panels are intentionally blank and labelled unavailable.
- `<city>_<primary|secondary>_<year>_geographic_support`: fixed geographical quadrants and count of strictly supported monthly matched pixels.
- `<city>_<primary|secondary>_<year>_brightness_decomposition`: how baseline brightness classes contribute to the municipality-wide net radiance change, with exact-additive columns.

No figures should be interpreted as independently verified treatment footprints. **Nothing will be committed to GitHub as a measured result until the actual EcoVis outputs have been reviewed.**

## EcoVis run commands — established GitHub workflow

```bash
cd /run/media/bznwr/data/lightpollution
git status --short
git switch publication-figures
git pull --ff-only origin publication-figures

# First: synthetic self-test, no VIIRS files read
python src/local_viirs_spatial_dilution.py --self-test

# Then: targeted spatial dilution check on the two inconclusive/heterogeneous cases
python src/local_viirs_spatial_dilution.py --cases portsmouth bristol

# Finally: all five primary cases + four secondary follow-ups
python src/local_viirs_spatial_dilution.py
```

If importing `local_viirs_five_case_analysis.py` fails, verify the script is present in the same `src/` tree. The new spatial analysis reuses its already-tested GeoTIFF indexing and city boundary functions. **Do not overwrite an untracked script during Git pulls**; investigate `git status` before resolving conflicts.

## First review gate

Inspect `brightness_strata_summary.csv` for all five cases, especially Portsmouth and Bristol, together with **primary** spatial-anomaly maps. Questions:

1. Does Portsmouth have a subregion with **an additional negative post-minus-pre anomaly**, despite its citywide post-minus-pre being positive? Are those subregions in independent 2019 trial neighbourhoods or could they reflect non-streetlight sources?
2. Does Bristol's 2023 near-zero aggregate hide **opposing signed contributions** from different baseline brightness bands or quadrants? Is the secondary 2024 decline spatially concentrated in a few baseline bright sources or diffuse?
3. Does Coventry's robust event-aligned decline occur across *many* baseline brightness bands/geographic sectors or mainly high-radiance points?
4. Is Échirolles' substantial 2022 decline broad across its baseline strata and sectors after restricting to strict common-month support?
5. Is Thanet's 2017 +9.4% driven by the **top 5%** bright source pixels and the same two months that dominated the citywide effect?

If the answer shows promising specific subareas, then join **independently documented lighting-treatment streets/circuits** (starting Coventry's official list, or Portsmouth's original 2019 pilot neighbourhoods) to the VIIRS grid; do **not** label low-ΔR pixels as switched off.
