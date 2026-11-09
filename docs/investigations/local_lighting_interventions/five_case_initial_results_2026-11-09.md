# Five city VIIRS VNL V1 monthly case studies — first full results

**Review date:** 9 November 2026. **Status:** grounded in investigator-uploaded `five_case_contrast_summary.csv`, five city seasonal-comparison figures, and five primary period maps generated on EcoVis. These data/figures remain on EcoVis; this GitHub Markdown report records the observed summary, not a new independently executed measurement. No comparison is a causal effect estimate; no formal significance or flux partition is inferred. Locked 2014–2024 European NUTS3, 2017 offset, 2021 snow handling and spider/moth PRC files have **not** been changed.

## The observed municipality-scale contrasts

All values are **arithmetic means of municipality monthly means** for the predeclared same-calendar-month outcomes, minus the corresponding means across reference years. Units for absolute differences: nW cm⁻² sr⁻¹. Percentage is of the reference-year means. Municipal areas use a fixed INSEE commune or ONS 2024 LAD footprint, not fixed NUTS3 for Portsmouth/Coventry/Bristol/Thanet.

| Case and comparison | Reference mean | Outcome mean | ΔR | ΔR (%) | Outcome valid pixels (%) | Outcome / reference median monthly cf_cvg (mean across chosen months) | Valid matched-grid cells | Matched-grid spatial median ΔR |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| Échirolles Oct–Dec 2022 vs 2020–21 | 43.781 | 27.547 | −16.234 | **−37.08%** | 100.0 | 8.00 / 8.67 | 54 | −15.589 |
| Échirolles Oct–Dec 2023 vs 2020–21 (secondary) | 43.781 | 18.458 | −25.322 | **−57.84%** | 100.0 | 7.33 / 8.67 | not mapped in this run | — |
| Portsmouth Sep–Dec 2020 vs 2018–19 | 24.125 | 20.480 | −3.645 | **−15.11%** | 100.0 | 10.50 / 11.13 | 293 | −3.105 |
| Portsmouth Sep–Dec 2021 vs 2018–19 (secondary) | 24.125 | 19.536 | −4.589 | **−19.02%** | 99.8 | 11.50 / 11.13 | not mapped in this run | — |
| Coventry Aug–Dec 2024 vs 2022–23 | 16.412 | 11.122 | −5.290 | **−32.23%** | 100.0 | 10.20 / 12.00 | 751 | −4.619 |
| Bristol Apr–Dec 2023 vs 2020–22 | 19.022 | 19.025 | +0.003 | **+0.02%** | 100.0 | 9.11 / 10.30 | 821 | −0.689 |
| Bristol Apr–Dec 2024 vs 2020–22 (secondary) | 19.022 | 15.560 | −3.461 | **−18.20%** | 99.97 | 6.44 / 10.30 | not mapped in this run | — |
| Thanet Jun–Nov 2017 vs 2015–16 | 10.695 | 11.700 | +1.005 | **+9.40%** | 100.0 | 11.00 / 10.25 | 767 | +0.355 |
| Thanet Jun–Nov 2018 vs 2015–16 (secondary) | 10.695 | 11.832 | +1.137 | **+10.63%** | 100.0 | 10.00 / 10.25 | not mapped in this run | — |

**Beware:** The `matched-grid spatial median ΔR` comes from native-cell, matched-calendar-month differences and is **not** the same estimand as the arithmetic mean of municipality-wide monthly mean radiance. `matched_period_pixels` are spatially autocorrelated, **not independent measurements**. `valid_percent=100%` means some `cf_cvg>0` for nearly every static pixel in each month, not constant night-time weather/acquisition.

## Interpretation of user-provided seasonal and primary spatial plots

### Échirolles: strong, sustained decline; highest priority for independent control

- The **Oct–Dec 2022** decline persists when the reference is the average of **both 2020 and 2021**, rather than the earlier 2021-only pilot. The matched-pixel map shows widespread negative ΔR, largest around the previously brightest central/eastern core, not one isolated hotspot.
- The separately summarised **Oct–Dec 2023** radiance is **57.84% below 2020/21**, stronger than 2022. This is a *follow-up*, but **2023 LED conversion and other changes contaminate attribution**. The gradual 2019–2021 decline visible in the season trajectories also requires explicit pre-trend comparisons.
- The supplied maps do not independently test atmospheric/snow differences or compare an untreated neighbouring area. No stage claims a fraction of reduction attributable to the 2022 operational switch-off.

### Coventry: strong, spatially broad event-compatible decline

- **Aug–Dec 2024** radiance is **32.23% lower** than the pooled 2022/23 seasonal reference, with a predominantly negative native-grid ΔR map and particularly marked reduction around several high-brightness patches of the city. Contrary to the positive centre seen in some scenes, this is not a one-pixel artifact.
- Council reports May–July 2024 part-night switch-off; the main LED programme was subsequently planned to begin in 2025. This is a strong *candidate* for a real change in lamp operation, but a 2024 pre-period and earlier-year placebo are needed to isolate the intervention timing from a pre-existing trend.
- The outcome has **100% mean valid pixel-centre support**, but its median `cf_cvg` averages **10.2** versus **12** in the reference. Run low-observation sensitivity before making a strong quantitative claim.

### Portsmouth: substantial negative 2020 contrast, but pre-existing drift

- **Sep–Dec 2020** is **15.11% lower** than 2018/19, and **Sep–Dec 2021** is **19.02% lower**. The spatial map shows heterogeneous declines with substantial darkening across multiple built-up portions of Portsmouth.
- The seasonal graph also shows **2019 and early-2020 declines relative to 2017/18**, preceding the summer 2020 citywide adaptive dimming window. This challenges attribution of the **entire −15.11%** to 2020 policy; 2018/19 LED conversion, weather and other trends all matter.
- Need a **2020 Jan–May pre-window anomaly vs the 2018/19 same-month reference**, a `2019 vs 2017/18` placebo and column-specific temporal context before concluding any incremental dimming signal.

### Bristol: primary 2023 aggregate near zero, with spatial heterogeneity; 2024 sensitive to coverage

- For **Apr–Dec 2023**, net mean change **+0.02%** hides spatial mixtures. The map shows positive changes in parts of the west/northwest and isolated southern hotspots, alongside negative areas elsewhere. Rolling LED + CMS installation began during 2023 and did not finish at Dec 2024; the municipal mean cannot be interpreted as a completed-treatment effect.
- **Apr–Dec 2024** is **−18.20%** versus the pre-2023 seasonal baseline, but the uploaded seasonal plot shows **November 2024 radiance collapsing to ~5–6 nW cm⁻² sr⁻¹** with the **median cloud-free observation count near 2**. This raises a specific risk of sparse acquisition/selection, not proof of a physically abrupt lighting shutdown.
- The immediate sensitivity is to exclude November 2024, then compare `cf_cvg>=3` and `>=5` mean radiance where enough cells remain, and inspect individual months in the spatial coverage maps.

### Thanet: modest positive mean, but spatially hotspot-dominated

- **Jun–Nov 2017** is **9.40% brighter** than the 2015/16 same-calendar-month mean. The native-grid difference map is predominantly close to zero/heterogeneous, apart from a conspicuous **positive compact hotspot around longitude 1.28, latitude 51.35**. Its source has **not been identified**. The median matched-cell change is only **+0.355**, which demonstrates the difference between a municipality mean and typical spatial change.
- Consequently it is **premature to claim widespread brightening due to longer streetlight operating hours**. That compact source could be unrelated to streetlighting and merits independent land-use/luminaire verification (without selecting a treatment area post hoc).
- **2018 Jun–Nov** is **10.63% above** 2015/16 reference; separately, the seasonal chart shows anomalously high **Nov–Dec 2018** and early **2019** radiances. A secondary spatial map and leave-month-out check are essential before this serves as a supporting example.

## Prioritised diagnostic next step

New standalone script [`src/local_viirs_five_case_robustness.py`](../../../src/local_viirs_five_case_robustness.py) reads **only the five case-level 120-month CSVs** created by the full stage, not the 480-MB archived raster tiles. It computes:

1. **Same-calendar-month anomalies** across all twelve outcome-year months, explicitly shading the transition window and reporting `cf_cvg`.
2. **Pre-versus-post anomaly differential**, with pre-months chosen to avoid the documented switch-over: Échirolles Apr–Jun; Portsmouth Jan–May; Coventry Jan–Apr; Bristol Jan–Feb (only two months; very weak). Thanet has no valid within-2017 pre-period because installation was Dec 2016–May 2017.
3. **Earlier-year placebo contrasts** (previous outcome year against years preceding it) to expose pre-existing downward/upward drift. These are descriptive comparators, NOT certified untreated controls; Portsmouth's placebo includes LED years.
4. **Reference-year sensitivity**, `cf_cvg>=3/5` cell-threshold means, and whole-month exclusion if median `cf_cvg<5` or `<10` in outcomes/reference; thresholding changes the spatial sample and is NOT automatically a correction.
5. **Leave-one-month-out** sensitivity to check November 2024 Bristol and late-2017/2018 Thanet signals.

**Potential next spatial task**: export secondary maps for Échirolles 2023, Bristol 2024 and Thanet 2018 and examine source hotspots, with special attention to months with `cf_cvg` below 5. No claim of causal attribution should precede these tests.

## EcoVis commands (established GitHub workflow)

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures

python src/local_viirs_five_case_robustness.py --self-test
python src/local_viirs_five_case_robustness.py
```

**Outputs:** `data/processed/local_lighting_interventions/five_case_robustness/` for three CSVs:
- `reference_and_cloudfree_sensitivity.csv`
- `pre_post_anomaly_differentials.csv`
- `earlier_year_placebo_comparisons.csv`

And `figures/local_lighting_interventions/five_case_robustness/` for five PNG/PDF monthly anomaly plots.

Review those three CSVs **before** deciding whether to fit a trend model. No changes to annual VIIRS, snow handling, 2017 zero baseline or polarisation receptor contrast.
