# Investigation record — 2021 Annual VNL VIIRS anomaly

**Status:** ACTIVE  
**Started:** 2026-10-12  
**Current restart commit:** `2d79e1da`  
**Primary question:** what caused the large 2021 radiance pulse and 2022 reversal in the fixed-NUTS2016 Annual VNL series?

## Why this investigation exists

A fixed-NUTS2016 2014-2024 re-extraction was designed to remove two sources of ambiguity in the historical future-response analysis:

1. the influential 2013 endpoint;
2. changes in NUTS geography after 2020.

The fixed extraction reproduced production values exactly over 2014-2020, then showed a large 2021 rise followed by a 2022 reversal. Because the same 1,514 NUTS3 regions and the same fixed density classes are used every year, the feature cannot be a boundary/cohort artifact.

The future radiance response must remain unfrozen until this disturbance is understood.

---

## Locked diagnostic substrate

```text
fixed geography: NUTS2016 Level 3
regions: 1,514
countries: 31
years: 2014-2024
density reference: median Eurostat density, 2014-2020
production pixel calibration from 2017:
    positive raw pixel -> max(raw - 0.235, 0)
```

Canonical fixed-panel files:

```text
data/processed/fixed_nuts2016_viirs_2014_2024.parquet
data/processed/fixed_nuts2016_viirs_2014_2024.csv
data/processed/fixed_nuts2016_population_density_reference.csv
```

Internal validation over 2014-2020:

```text
10,192 region-year rows
max absolute corrected difference from production = 0
```

---

## Evidence ledger

### A. First-difference fixed-panel diagnostic

Stage:

```bash
python src/run_pipeline.py future-radiance-2021-diagnostic
```

Key results:

```text
2020 -> 2021 corrected:
    median regional change +10.033%
    78.5% regions positive

2021 -> 2022 corrected:
    median regional change -8.345%
    12.15% regions positive

country-balanced:
    2021 +18.9995%, 83.9% countries positive
    2022 -10.861%, 3.23% countries positive
```

**KEEP:** a common 2021 pulse and 2022 reversal exist on fixed geography. Do not force one 2014-2024 linear slope.

### B. Raw-radiance reversal diagnostic

Stage:

```bash
python src/run_pipeline.py future-radiance-2021-reversal
```

```text
raw 2020 -> 2021  +9.47686%
raw 2021 -> 2022  -7.8411%

corrected 2020 -> 2021 +10.0331%
corrected 2021 -> 2022 -8.34522%

2021 vs 2022 regional-change rho = -0.428565
73.7% opposite sign
70.1% up-then-down
2020 -> 2022 net raw change ~ +0.893%
```

**KEEP:** the event is present before the project's 0.235 correction. It is not created by production calibration.

### C. Annual product-layer audit

Stage:

```bash
python src/run_pipeline.py future-radiance-product-layer-audit
```

```text
2021:
    average +9.4769%
    median  +7.0858%

2022:
    average -7.8411%
    median  -4.1485%

average-vs-median regional change rho in 2021 = 0.890
```

**KEEP:** not a mean-only outlier effect; the robust annual median carries the pulse.

### D. Absolute-offset audit

Stage:

```bash
python src/run_pipeline.py future-radiance-absolute-offset-audit
```

Country-balanced raw median absolute change:

```text
2017 +0.06609 nW
2020 -0.00742
2021 +0.115858  <- largest
2022 -0.097972  <- second largest, opposite sign
2023 -0.015425
2024 +0.005047
```

**KEEP:** 2021 has a common additive-like component but is not well described by one constant annual offset across brightness.

### E. Mask occupancy audit

Stage:

```bash
python src/run_pipeline.py future-radiance-mask-occupancy-audit
```

```text
country-balanced positive-area change:
    2021 +2.2455 percentage points
    2022 -1.328 percentage points

rho(mask-area change, raw radiance change):
    2021 +0.4878
    2022 +0.390

79.9% of regions increased positive occupancy in 2021
```

**KEEP:** changing support is a major contributor, especially in dim regions.

### F. Exact mask/intensity decomposition

Stage:

```bash
python src/run_pipeline.py future-radiance-mask-intensity-decomposition
```

2021 country-balanced exact decomposition:

```text
occupancy +0.065908 nW
intensity +0.048709 nW
total     +0.114617 nW

occupancy share 57.50%
intensity share 42.50%
```

2022:

```text
occupancy -0.035487 nW
intensity -0.056446 nW
```

The 2017 positive control was almost entirely intensity. In density >=1000, the 2021 occupancy term was only about +0.0056 nW versus intensity about +0.2160 nW.

**KEEP:** 2021 is a two-component disturbance. Dense-region behaviour is overwhelmingly intensity, so "mask expansion" is insufficient.

### G. Pixel-support persistence

Stage:

```bash
python src/run_pipeline.py future-radiance-2021-pixel-support
```

2020->2021:

```text
total                    +0.115858 nW
new-positive contribution +0.027913
lost-positive contribution -0.009303
persistent-positive       +0.095567
```

Stable positive 2020-2023 support:

```text
conditional radiance:
    2021 +0.61017 nW
    2022 -0.37607 nW
```

Approximately half of new-positive area was transient.

**KEEP:** reject "mainly transient mask artifact". The same persistently lit pixels rise and reverse.

### H. July VCM/VCMSL pilot

Stage:

```bash
python src/run_pipeline.py future-radiance-july-vcm-vcmsl-pilot
```

The common-support VCM/VCMSL comparison was geographically limited, and July 2021 VCMSLCFG had no usable monthly coverage for several core Nordic countries.

**KEEP:** the documented early-July 2021 stray-light LUT issue cannot by itself explain the Nordic-amplified annual pulse.

### I. Full monthly VCMSLCFG localization

Stage:

```bash
python src/run_pipeline.py future-radiance-monthly-2021-localization
```

Each calendar month uses only pixels with `cf_cvg > 0` in all 2020/2021/2022 versions of that same month.

Country-balanced primary anomalies:

```text
Jan +0.491459 nW
Feb +0.817785 nW
Mar -0.037181 nW
Apr -0.014609 nW
May -0.193368 nW
Jun +0.016569 nW
Jul +0.034040 nW
Sep -0.021495 nW
Oct +0.015005 nW
Nov +0.020344 nW
Dec +0.203571 nW
```

August 2022 is J01/NOAA-20 rather than NPP and is retained only as a mixed-platform sensitivity result.

```text
positive regions:
    Jan 86.8%
    Feb 91.9%

latitude rho:
    Jan +0.210
    Feb +0.616
    Mar -0.319
```

**KEEP:** the annual event is strongly localized to January-February, particularly February. It is not a persistent year-round step.

---

## Paper-facing evidence artifacts and narrative role

The anomaly investigation is likely to generate several supplementary figures.
The table below records the files that currently carry the strongest evidential
weight, the conclusion each supports, and why that conclusion led to the next
analysis. This is intended to prevent later manuscript drafting from selecting
figures retrospectively without preserving the diagnostic logic.

| Diagnostic step | Key figure(s) | Key dataset(s) | Conclusion supported and rationale for the next step |
|---|---|---|---|
| Locked 2014-2024 geometry | `figures/population_relationship/fixed_nuts2016_2014_2024/real_trajectories_2014_2024.pdf`; `figures/population_relationship/fixed_nuts2016_2014_2024/density_trends_2014_2020_vs_2014_2024.pdf` | `data/processed/fixed_nuts2016_viirs_2014_2024.parquet`; `data/processed/fixed_nuts2016_density_trajectories_2014_2024.csv`; `data/processed/fixed_nuts2016_2014_2020_production_validation.csv` | Establishes that the post-2020 feature is present on one fixed NUTS2016 cohort and that the 2014-2020 extraction exactly reproduces production. This justified treating 2021 as a product/time-series problem rather than a geography problem. |
| First-difference 2021 diagnostic | `figures/population_relationship/fixed_nuts2016_2021_diagnostic/annual_within_region_change.pdf`; `figures/population_relationship/fixed_nuts2016_2021_diagnostic/country_change_2020_2021.pdf` | `data/processed/fixed_nuts2016_annual_change_summary_2015_2024.csv`; `data/processed/fixed_nuts2016_country_balanced_annual_change_2015_2024.csv`; `data/processed/fixed_nuts2016_2021_model_tests.csv` | Shows that 2021 is a large Europe-wide pulse rather than an ordinary point on a linear trajectory. The 2022 reversal then had to be tested explicitly. |
| Raw pulse and 2022 reversal | `figures/population_relationship/fixed_nuts2016_2021_reversal/raw_2021_pulse_vs_2022_reversal.pdf`; `figures/population_relationship/fixed_nuts2016_2021_reversal/raw_pulse_reversal_vs_latitude.pdf`; `figures/population_relationship/fixed_nuts2016_2021_reversal/raw_net_2020_2022_by_density_band.pdf` | `data/processed/fixed_nuts2016_2021_2022_region_reversal.csv`; `data/processed/fixed_nuts2016_2021_2022_reversal_diagnostics.csv`; `data/processed/fixed_nuts2016_2021_2022_country_reversal.csv` | Demonstrates that the event exists in raw radiance, is paired with a strong 2022 reversal, and has geographic structure. This ruled out the project's 0.235 correction as the cause and motivated interrogation of the EOG product itself. |
| Annual mean versus median product layers | `figures/population_relationship/viirs_annual_product_layer_audit/average_vs_median_masked_2020_2023.pdf`; `figures/population_relationship/viirs_annual_product_layer_audit/average_vs_median_2021_region_change.pdf`; `figures/population_relationship/viirs_annual_product_layer_audit/radiance_2021_vs_cf_cvg_change.pdf` | `data/processed/viirs_annual_product_layer_annual_summary_2021_2023.csv`; `data/processed/viirs_annual_product_layer_diagnostics_2021_2023.csv`; `data/processed/viirs_annual_product_layer_density_summary_2021_2023.csv` | The pulse survives the annual median product and is not explained by annual cloud-free count alone. This pushed the analysis from summary-statistic robustness toward support/intensity decomposition. |
| Absolute offset test | `figures/population_relationship/fixed_nuts2016_absolute_offset_audit/annual_absolute_change_raw_vs_corrected.pdf`; `figures/population_relationship/fixed_nuts2016_absolute_offset_audit/raw_2021_absolute_change_vs_baseline_brightness.pdf` | `data/processed/fixed_nuts2016_absolute_change_annual_summary_2015_2024.csv`; `data/processed/fixed_nuts2016_absolute_change_by_brightness_2015_2024.csv`; `data/processed/fixed_nuts2016_absolute_change_brightness_slopes_2015_2024.csv` | 2021 contains an additive-like common component but the change is brightness-dependent, unlike a simple universal offset. This made a mask/support versus within-support intensity decomposition necessary. |
| Mask occupancy | `figures/population_relationship/fixed_nuts2016_mask_occupancy_audit/annual_mask_occupancy_change_2015_2024.pdf`; `figures/population_relationship/fixed_nuts2016_mask_occupancy_audit/radiance_vs_mask_change_2021.pdf`; `figures/population_relationship/fixed_nuts2016_mask_occupancy_audit/mask_change_by_baseline_brightness.pdf` | `data/processed/fixed_nuts2016_mask_occupancy_annual_summary_2015_2024.csv`; `data/processed/fixed_nuts2016_mask_occupancy_by_brightness_2015_2024.csv`; `data/processed/fixed_nuts2016_mask_occupancy_by_density_2015_2024.csv` | Shows that more area is retained as positive in 2021, especially in dim regions. Because this did not address radiance changes on already-positive pixels, the next step was an exact decomposition. |
| Exact mask/intensity decomposition | `figures/population_relationship/fixed_nuts2016_mask_intensity_decomposition/annual_occupancy_vs_intensity_components.pdf`; `figures/population_relationship/fixed_nuts2016_mask_intensity_decomposition/decomposition_2021_by_density.pdf`; `figures/population_relationship/fixed_nuts2016_mask_intensity_decomposition/decomposition_2021_by_brightness.pdf` | `data/processed/fixed_nuts2016_mask_intensity_decomposition_annual_2015_2024.csv`; `data/processed/fixed_nuts2016_mask_intensity_decomposition_by_density_2015_2024.csv`; `data/processed/fixed_nuts2016_mask_intensity_decomposition_by_brightness_2015_2024.csv` | Establishes the two-component nature of the event and shows that dense-region 2021 change is overwhelmingly intensity rather than occupancy. This ruled out a pure mask-artifact narrative and motivated direct pixel persistence. |
| Pixel-support persistence | `figures/population_relationship/fixed_nuts2016_2021_pixel_support/new_2021_pixel_persistence_by_density.pdf`; `figures/population_relationship/fixed_nuts2016_2021_pixel_support/stable_pixels_2021_2022_by_density.pdf`; `figures/population_relationship/fixed_nuts2016_2021_pixel_support/pixel_components_2021_by_density.pdf` | `data/processed/fixed_nuts2016_2021_pixel_support_overall.csv`; `data/processed/fixed_nuts2016_2021_pixel_support_by_density.csv`; `data/processed/fixed_nuts2016_2021_pixel_support_by_country.csv` | This is one of the strongest likely supplementary analyses: the same pixels positive throughout 2020-2023 brighten in 2021 and reverse in 2022. It rejects the explanation that the annual pulse is mainly transient support expansion and shifts attention to seasonal/product radiometry on stable support. |
| July VCM versus VCMSL pilot | `figures/population_relationship/fixed_nuts2016_july_vcm_vcmsl_pilot/july_vcmsl_added_coverage.pdf`; `figures/population_relationship/fixed_nuts2016_july_vcm_vcmsl_pilot/july_common_difference_vs_latitude.pdf`; `figures/population_relationship/fixed_nuts2016_july_vcm_vcmsl_pilot/july_2021_country_common_difference.pdf` | `data/processed/fixed_nuts2016_july_vcm_vcmsl_annual_summary.csv`; `data/processed/fixed_nuts2016_july_vcm_vcmsl_by_country.csv`; `data/processed/fixed_nuts2016_july_vcm_vcmsl_regional.csv` | Provides negative evidence against the simple hypothesis that the documented early-July 2021 stray-light LUT problem generated the annual northern pulse. Because July could not account for the affected northern countries, the anomaly had to be localized across the full year. |
| Full monthly localization | `figures/population_relationship/fixed_nuts2016_monthly_2021_localization/monthly_2021_anomaly_vs_flanks.pdf`; `figures/population_relationship/fixed_nuts2016_monthly_2021_localization/monthly_matched_radiance_2020_2022.pdf`; `figures/population_relationship/fixed_nuts2016_monthly_2021_localization/monthly_anomaly_latitude_rho.pdf`; `figures/population_relationship/fixed_nuts2016_monthly_2021_localization/monthly_country_anomaly_heatmap.pdf` | `data/processed/fixed_nuts2016_monthly_2021_localization_summary.csv`; `data/processed/fixed_nuts2016_monthly_2021_localization_by_country.csv`; `data/processed/fixed_nuts2016_monthly_2021_localization_by_density.csv` | Localizes the pulse to January-February, particularly February, with a strong northern structure. This is the key rationale for moving from generic sensor diagnostics to the published snow/albedo mechanism. |
| Snow/albedo test — **clean rerun validated 2026-10-14** | `figures/population_relationship/fixed_nuts2016_snow_viirs_test/february_snow_vs_viirs_anomaly.pdf`; `snow_viirs_anomaly_rho_jan_feb_mar.pdf`; `levin_style_jan_july_snow_radiance_2021.pdf`; `levin_style_percent_lit_thresholds.pdf` | `data/processed/fixed_nuts2016_snow_2021_anomaly_tests.csv`; `fixed_nuts2016_snow_2021_anomaly_models.csv`; `fixed_nuts2016_snow_levin_jan_july_summary.csv`; `fixed_nuts2016_snow_levin_percent_lit_thresholds.csv`; `fixed_nuts2016_mod10cm_snow_monthly_2020_2022.csv` | The corrected MOD10CM extraction obeys the physical 0-100% snow range. February shows a strong snow-radiance anomaly association at region and country scales, survives country/latitude/baseline-brightness/cloud-free-coverage adjustment, and has a strong positive snow x brightness interaction; March is essentially null. This moves snow/surface reflectance from an open hypothesis to a strongly supported explanation for a major component of the February 2021 pulse. The next step is to test the mechanism more directly with MODIS albedo and NDVI rather than treating snow cover itself as the final causal variable. |

### Likely manuscript structure

For a future paper, the diagnostic chain above naturally separates into two
levels. The main narrative can probably remain compact: locked geography shows a
2021/2022 excursion; raw/median/stable-support tests establish that it is real
within the Annual VNL product rather than an artefact of the project's
calibration or changing support alone; monthly localization identifies
January-February; and the snow/albedo analysis tests the resulting mechanism.

The intermediate product audits are particularly well suited to supplementary
material because they document alternative explanations that were actively
tested and rejected. The most likely high-value supplementary figures are:

1. fixed-panel 2014-2024 trajectories;
2. raw 2021 pulse versus 2022 reversal;
3. mean-versus-median product comparison;
4. occupancy/intensity decomposition;
5. persistent-pixel 2021/2022 reversal;
6. monthly January-February localization and latitude structure;
7. the validated snow/anomaly mechanism figure and January/February/March control comparison.

The corresponding CSVs in the table should be treated as the source-of-truth
publication tables. Final figures can be reformatted later without losing the
evidential provenance.

## Current leading hypothesis: winter snow / surface albedo

The timing and latitude pattern make winter surface state a serious candidate.

Published methodological basis:

1. Noam Levin (2017), *The impact of seasonal changes on observed nighttime brightness from 2014 to 2015 monthly VIIRS DNB composites*, Remote Sensing of Environment 193:150-164. DOI 10.1016/j.rse.2017.03.003.
2. Noam Levin & Qingling Zhang (2017), *A global analysis of factors controlling VIIRS nighttime light levels from densely populated areas*, Remote Sensing of Environment 190:366-382. DOI 10.1016/j.rse.2017.01.006.

The project snow test deliberately follows their logic rather than inventing a one-off correlation:

- MOD10CM monthly mean snow cover;
- mean VIIRS radiance;
- cloud-free coverage;
- Spearman rank relationships;
- January-vs-July comparison;
- percent-lit area at 2, 5, 10, 25, 50, 100 and 250 nW cm^-2 sr^-1.

Project-specific anomaly test:

```text
snow anomaly = snow_2021 - mean(snow_2020, snow_2022)
VIIRS anomaly = VIIRS_2021 - mean(VIIRS_2020, VIIRS_2022)
```

Pre-specified roles:

```text
February = primary
January  = supportive
March    = negative control
```

The adjusted sensitivity model controls for country, latitude, baseline brightness and cloud-free coverage anomaly and includes a snow x baseline-brightness interaction.

---

## Snow-run quality-control incident

The first complete MOD10CM run yielded a very strong apparent February relationship, but some calculated regional snow means exceeded 100%.

MOD10CM encodes only 0-100 as snow percentage; >100 values are categorical codes. Therefore those apparent snow results are **invalid**.

The invalid preliminary numbers are preserved only as debugging history and must not be used as evidence.

Commit `2d79e1da` changes the cached raster construction so every source value >100 is physically written as nodata before `exactextract` sees it. A hard 0-100 invariant now stops the pipeline if any categorical code leaks through.

## Clean corrected snow result — 2026-10-14

The corrected rerun completed successfully after two further implementation fixes:
`e5abd14` corrected an accidental literal newline in the source, and
`7674701` distinguished true range violations from roundoff-scale values such
as 100.00000000000001. The stored snow product now obeys the physical ranges:

```text
snow_mean_pct:              0 to 100
snow_valid_area_fraction:   0 to 1
MOD10CM inventory:          12 months 2020, 12 months 2021, 11 months 2022
structural missing month:   October 2022
```

The pre-specified anomaly was:

```text
2021 - mean(2020, 2022)
```

with February primary, January supportive and March the negative control.

### Primary February result

```text
region-level Spearman rho = +0.6653
p = 2.17e-192
n = 1500 NUTS3

country-median Spearman rho = +0.7383
p = 2.13e-06
n = 31 countries

median snow anomaly = +5.865 percentage points
median VIIRS anomaly = +0.615 nW cm^-2 sr^-1
```

In the pre-specified adjusted model controlling for country, latitude, baseline
brightness and cloud-free-coverage anomaly:

```text
base R2 = 0.3021
full R2 = 0.4979
partial R2 of snow block = 0.2805

snow coefficient at mean log-baseline brightness
    +0.08359 nW per snow percentage point
    95% CI +0.05814 to +0.10904
    p = 1.21e-10

snow x baseline-brightness interaction
    coefficient +0.23077
    p = 1.94e-13

joint snow block
    p = 2.48e-14

cloud-free-coverage anomaly
    p = 0.314
```

**KEEP:** the February 2021 radiance anomaly is strongly associated with
anomalous snow cover, the relationship is not removed by country, latitude,
baseline brightness or cloud-free coverage, and the positive interaction
supports the physical prediction that the absolute snow-reflection effect is
larger where more artificial light is available.

### January supportive result

```text
region rho = +0.4875
country-median rho = +0.5407
adjusted snow coefficient = +0.01976 nW per snow percentage point
snow-block partial R2 = 0.1344
joint snow-block p = 9.35e-08
snow x brightness p = 9.00e-05
```

**KEEP:** January independently supports the same winter-surface mechanism,
although less strongly than February.

### March negative control

```text
region rho = -0.0222, p = 0.391
country-median rho = +0.2040, p = 0.271
snow-block partial R2 = 0.00173
snow coefficient p = 0.346
snow x brightness p = 0.558
```

The robust joint Wald p-value for the two snow terms is borderline
(`p = 0.049`), but the individual snow terms are non-significant and the
incremental explanatory power is only ~0.17% of residual variance.

**KEEP:** March is a practically null negative control. The February result is
therefore not simply a generic latitude/snow association that appears whenever
snow is included in a model.

### Levin-style January-minus-July replication

The January-minus-July relationship between snow change and radiance change is:

```text
2020 region rho +0.0795; country rho +0.0470
2021 region rho +0.4644; country rho +0.5000
2022 region rho +0.1679; country rho +0.3983
```

The 2021 relationship is substantially stronger than the flanking years.

The published percent-lit-area concept also shows a strong 2021 response,
especially at low-to-moderate radiance thresholds:

```text
2021 snow change vs percent-lit change:
    >2 nW    rho +0.5519
    >5 nW    rho +0.5708
    >10 nW   rho +0.5372
    >25 nW   rho +0.4508
    >50 nW   rho +0.3147
    >100 nW  rho +0.1791
    >250 nW  rho +0.0552
```

**KEEP:** the event is not confined to the mean-radiance statistic. In 2021,
greater winter snow is also associated with a larger fraction of area crossing
the published low-to-moderate VIIRS brightness thresholds.

### Important limit on generalization

Across all available paired months in 2020-2022, the median within-region
snow-radiance Spearman correlation is `-0.1751`, with only 27.8% of regions
positive.

This does **not** overturn the February anomaly result. It does mean that the
project should not claim a universal rule that more snow always makes every
European NUTS3 brighter in VIIRS. The evidence is strongest for a specific
winter 2021 anomaly and for January-vs-July contrast in that year. This is
consistent with using the Levin framework as a mechanism test while recognizing
that the present NUTS3 sample is much broader than selected urban areas.

### Revised causal interpretation

Snow cover itself is best treated as an indicator of a changed reflecting
surface, not as the final mechanism. The next analysis should therefore test
whether the 2021 snow anomaly produces a corresponding **surface-albedo anomaly**
and whether that albedo anomaly explains/attenuates the snow-radiance
association. NDVI should be added as the complementary vegetation/seasonal
surface control used in the published literature.

The next stage should preserve the same pre-specification:

```text
February = primary
January  = supportive
March    = negative control
2021 - mean(2020, 2022) for every surface variable
country + latitude + baseline brightness + cf_cvg controls
surface-state x baseline-brightness interaction
```

A sequential model comparison should then ask:

1. does albedo anomaly predict VIIRS anomaly?
2. does the snow coefficient attenuate when albedo enters?
3. does albedo x baseline brightness show the predicted stronger effect in brighter regions?
4. does NDVI explain additional seasonal structure?
5. does March remain a negative control?

This should be described as a **mechanistic attenuation test**, not formal causal
mediation, because all variables are observational.

---

## Current KEEP / REJECT / OPEN table

| Status | Conclusion |
|---|---|
| KEEP | 2021 rise and 2022 reversal survive fixed geography. |
| KEEP | The event exists in raw radiance before project calibration. |
| KEEP | Annual median confirms it is not a mean-only outlier phenomenon. |
| KEEP | Mask/support changes contribute substantially in dim regions. |
| KEEP | Persistent identical lit pixels also show the rise and reversal. |
| KEEP | Dense-region 2021 change is predominantly intensity. |
| KEEP | July-specific VCM/VCMSL stray-light issue is insufficient. |
| KEEP | January-February, especially February, dominate the monthly anomaly. |
| REJECT | Changing NUTS geography explains the event. |
| REJECT | The 0.235 calibration creates the event. |
| REJECT | It is mainly transient new positive pixels. |
| REJECT | The first MOD10CM snow-run statistics are valid evidence. |
| KEEP | Correctly masked MOD10CM snow anomaly strongly predicts the February 2021 VIIRS anomaly. |
| KEEP | The February snow association survives country, latitude, baseline brightness and cf_cvg controls. |
| KEEP | The positive snow x brightness interaction supports a reflection mechanism. |
| KEEP | January independently supports the winter-surface mechanism; March is a practically null negative control. |
| KEEP | Levin-style January-minus-July and percent-lit analyses show an unusually strong 2021 snow-linked response. |
| KEEP | The full 2020-2022 within-region seasonal correlation is not generally positive, so the conclusion is event-specific rather than universal. |
| OPEN | How much of the snow association is explained specifically by surface albedo? |
| OPEN | Does NDVI account for additional seasonal surface-state variation? |
| OPEN | How much residual 2021 signal remains after snow/albedo/NDVI, for possible compositing/calibration mechanisms? |

---

## Pre-registered albedo attenuation test — 2026-10-15

The next mechanism stage is now implemented as:

```bash
python src/run_pipeline.py future-radiance-albedo-2021-test
```

It is intentionally narrower than a new exploratory surface-model sweep. Its
single primary question is:

> Does visible surface-albedo anomaly predict the February 2021 VIIRS anomaly,
> and does adding albedo materially attenuate the validated snow block?

Pre-specification:

```text
months:
    February primary
    January supportive
    March negative control

surface variable:
    VNP43MA3.002 Albedo_WSA_vis primary
    Albedo_BSA_vis sensitivity

mandatory-quality handling:
    QA 0/1 primary
    QA 0 only stricter sensitivity

anomaly:
    2021 - mean(2020, 2022)

controls:
    country
    latitude
    baseline VIIRS brightness
    cf_cvg anomaly

interaction:
    snow x baseline brightness
    albedo x baseline brightness
```

The sequential models are frozen before viewing the albedo result:

```text
M0 = controls only
M1 = controls + snow + snow x brightness
M2 = controls + albedo + albedo x brightness
M3 = controls + snow + snow x brightness
              + albedo + albedo x brightness
```

The primary attenuation statistic is the reduction in the snow-block partial
R2 from `M1 vs M0` to `M3 vs M2`. Main-effect coefficient attenuation is
reported as a secondary summary because the snow x brightness interaction means
one snow coefficient alone does not represent the whole snow effect.

Both HC3 and country-clustered uncertainty are retained for the combined model.

### Absolute versus relative radiance is a hard requirement

The primary mechanism outcome is the **absolute** VIIRS anomaly in
nW cm-2 sr-1. The stage also calculates:

```text
relative anomaly (%) =
    100 x [R2021 - mean(R2020,R2022)] / mean(R2020,R2022)
```

and writes a density-class context table and paired figure showing absolute and
relative results.

This distinction is essential for the eventual historical/model reconciliation:
a large percentage change in a very dark NUTS3 can represent much less physical
radiance than a small percentage change in a bright NUTS3. The final
density-dependent historical response must therefore agree with observations in
both proportional and absolute radiance space rather than being tuned to the
percentage-trend figure alone.

### AppEEARS input layout

The stage recursively discovers the three JFM AppEEARS downloads under:

```text
data/raw/modis/VNP43MA3_002/
```

The downloads may be placed in separate year/request subdirectories or left as
ZIP bundles. The required layers are:

```text
Albedo_WSA_vis
Albedo_BSA_vis
BRDF_Albedo_Band_Mandatory_Quality_vis
```

All daily dates from January-March 2020, 2021 and 2022 must be present before
the analysis runs. Missing dates/layers cause a hard stop with an inventory
report rather than a partial analysis.

Daily BRDF/albedo products are aggregated to one monthly surface-state raster
before NUTS3 extraction. The code checks whether AppEEARS has already applied
albedo scale factors and refuses ambiguous scaling instead of silently guessing.

### Outputs reserved for the paper trail

```text
data/processed/fixed_nuts2016_vnp43ma3_manifest.csv
data/processed/fixed_nuts2016_vnp43ma3_monthly_2020_2022.csv
data/processed/fixed_nuts2016_albedo_snow_viirs_anomaly_regional.csv
data/processed/fixed_nuts2016_albedo_2021_anomaly_tests.csv
data/processed/fixed_nuts2016_albedo_snow_attenuation_models.csv
data/processed/fixed_nuts2016_albedo_viirs_density_context.csv
data/processed/fixed_nuts2016_albedo_viirs_test_summary.txt

figures/population_relationship/fixed_nuts2016_albedo_viirs_test/
    february_wsa_albedo_vs_viirs_anomaly_absolute.pdf
    february_wsa_albedo_vs_viirs_anomaly_relative.pdf
    february_snow_albedo_partial_r2.pdf
    albedo_viirs_rho_jan_feb_mar.pdf
    february_viirs_anomaly_absolute_relative_by_density.pdf
```

If the albedo pathway is supported, this stage still does **not** directly
rewrite the annual 2013-2024 series. The following step is to define a
scientifically defensible control/correction for the 2021 surface-reflectance
disturbance and then refit the long-term density response, with the
historical-distribution comparison regenerated in both relative and absolute
radiance space.

## Current next step

The snow hypothesis has now passed the clean pre-specified test. Do **not**
repeat the snow stage unless code/data change.

Next: acquire MODIS surface albedo and NDVI for the same 2020-2022 period and
build a surface-state mechanism stage. Preferred primary albedo variable is a
visible-band broadband MODIS BRDF/albedo measure on a CMG grid; both white-sky
and black-sky albedo should be retained, with one pre-specified as primary and
the other as sensitivity. NDVI should use a monthly CMG product to match the
published Levin framework.

Key validated snow outputs to preserve:

```text
data/processed/fixed_nuts2016_snow_viirs_test_summary.txt
data/processed/fixed_nuts2016_snow_2021_anomaly_tests.csv
data/processed/fixed_nuts2016_snow_2021_anomaly_models.csv
data/processed/fixed_nuts2016_snow_levin_jan_july_summary.csv
data/processed/fixed_nuts2016_snow_levin_percent_lit_thresholds.csv
data/processed/fixed_nuts2016_snow_levin_seasonal_correlations.csv
data/processed/fixed_nuts2016_mod10cm_snow_monthly_2020_2022.csv

figures/population_relationship/fixed_nuts2016_snow_viirs_test/february_snow_vs_viirs_anomaly.pdf
figures/population_relationship/fixed_nuts2016_snow_viirs_test/snow_viirs_anomaly_rho_jan_feb_mar.pdf
figures/population_relationship/fixed_nuts2016_snow_viirs_test/levin_style_jan_july_snow_radiance_2021.pdf
figures/population_relationship/fixed_nuts2016_snow_viirs_test/levin_style_percent_lit_thresholds.pdf
```

## Key commits in this investigation

```text
1c6c63f4  Add fixed NUTS2016 2021 discontinuity diagnostic
6a0348d0  Add 2021 pulse and 2022 reversal diagnostic
6f1cba36  Add VIIRS annual product-layer audit
105221a7  Add fixed-panel absolute radiance offset audit
4ccd5f9d  Add fixed-panel masked-area occupancy audit
932480c2  Add exact mask-intensity decomposition
4ea6ca81  Add 2021 pixel-support persistence audit
676830ae  Add July VCM versus VCMSL pilot diagnostic
f87ec39c  Add monthly localization of 2021 VIIRS anomaly
21d79b1c  Handle NOAA-20 August 2022 as mixed-platform sensitivity
544f6bd2  Add literature-led snow test for 2021 VIIRS anomaly
ee794d6f  Normalize MOD10CM CMG CRS to EPSG:4326 without resampling
f07c53f8  Fix NUTS ID alignment in snow anomaly pivots
2d79e1da  Mask MOD10CM categorical codes before snow extraction
e5abd14   Fix literal newline in snow test constants
7674701   Allow only roundoff-scale MOD10CM boundary overshoot
```

## Validated albedo attenuation result — 2026-10-15

The pre-registered VNP43MA3 visible-albedo attenuation stage completed
successfully after the CRS fix. The primary analysis was February 2021,
white-sky visible albedo (WSA), mandatory QA 0/1, absolute VIIRS anomaly.

### Primary February result

```text
n NUTS3                         1497
region Spearman rho             +0.7198
country-median rho              +0.8440

M0 controls R2                   0.3026
M1 controls + snow R2            0.4981
M2 controls + albedo R2          0.4657
M3 controls + snow + albedo R2   0.5079

snow partial R2 before albedo     0.2803
snow partial R2 after albedo      0.0789
snow-block attenuation           71.8%

albedo partial R2 after snow      0.0196
albedo block p, HC3               0.0073
albedo block p, country cluster   0.0121

snow block p, HC3                 0.0013
snow block p, country cluster     0.00037
```

This meets the pre-specified primary criterion: February albedo anomaly strongly
tracks the absolute VIIRS anomaly, adds explanatory power beyond snow, and
materially attenuates the snow block. Snow remains independently informative,
so the result supports a shared surface-state pathway rather than complete
replacement of snow by albedo.

The February snow and WSA albedo anomalies are themselves strongly correlated
(Spearman rho approximately 0.85), so coefficient-level interpretation is
limited by collinearity. The pre-specified block-level attenuation / partial-R2
comparison is therefore the more appropriate summary.

### Replication and controls

January independently supports the same mechanism. For primary WSA QA0/1,
region rho is +0.5107 and country-median rho +0.7677. The snow block attenuates
by about 90.8% after albedo enters, while albedo retains substantial independent
explanatory power.

March is a useful negative control. The raw WSA QA0/1 rank correlation is small
(+0.1548) and becomes practically null in the adjusted model: albedo partial R2
after snow is about 0.0011, with albedo-block p=0.834 (HC3) and p=0.882
(country-clustered). The very small unadjusted p-value is therefore a sample-size
effect rather than meaningful explanatory power.

Black-sky albedo closely reproduces the primary February result: region rho
+0.7222, snow-block attenuation 71.6%, and albedo partial R2 after snow 0.0210.

### QA sensitivity

The QA0-only sensitivity is much smaller in winter (February n=310, 17
countries, versus 1497 / 31 for QA0/1) and does not independently support the
albedo block after adjustment. This is not treated as equivalent evidence
against the primary result because QA0 requires full BRDF inversion and
therefore conditions on a strongly selective subset of winter observations.
However, it is an important limitation and must be reported: the clean
Europe-wide result depends on including both full and magnitude BRDF inversions
(QA 0/1), as pre-specified.

### Absolute versus relative radiance

The February density-class result demonstrates why percentage trends cannot be
read as physical radiance changes:

```text
density       median absolute anomaly     median relative anomaly
<30                  0.234 nW                    35.4%
30-100               0.314 nW                    37.9%
100-300              0.768 nW                    55.5%
300-1000             1.919 nW                    55.9%
>=1000               5.624 nW                    45.3%
```

Thus the >=1000 class has an absolute February anomaly roughly 24 times the
<30 class despite broadly comparable proportional anomalies. This is now a hard
constraint on the historical-model reconciliation: both absolute and relative
radiance space must be examined.

### Scientific interpretation

The evidence now strongly supports the January-February 2021 disturbance as a
surface-reflectance event. Snow cover and visible albedo are overlapping
indicators of that surface state; albedo substantially explains the component
previously attributed to snow, but snow retains independent information.

This is mechanistic attenuation in observational data, not formal causal
mediation.

The next model-development task should return to the main historical question:
how to control winter surface-state disturbances without tuning specifically to
2021, then refit the density-dependent long-term response and regenerate the
historical validation in both absolute and proportional radiance units.

## Pre-registered historical snow-control test

Following the validated January-February snow and albedo mechanism work, the
next question is whether winter surface state materially biases the inferred
2014-2024 density-dependent annual VIIRS trends.

Stage:

```bash
python src/run_pipeline.py future-radiance-historical-snow-control
```

This stage is diagnostic and uses the locked fixed-NUTS2016 2014-2024 cohort
and density classes. It does not alter annual VIIRS values and does not assign
a special correction to 2021.

### Surface-state covariate

For every NUTS3 region and year, January and February MOD10CM snow are
extracted with the already validated physical-code masking and area-weighted
zonal method. Each monthly value is converted to a region/month-specific
leave-one-out anomaly:

```text
snow anomaly(year, month)
  = snow(year, month) - mean[snow(other years, same region/month)]
```

The winter surface-state index is the mean of the available January/February
month-specific anomalies. This choice is important for the genuine 2016-02
archive gap: January 2016 can contribute as a January *anomaly* without being
treated as a substitute for a raw two-month winter mean.

Primary sample: at least one valid winter month.
Sensitivity sample: both January and February valid.

### Matched model comparison

For each sample and each reporting space, the same observations are used
before and after surface-state control.

Unadjusted:

```text
radiance ~ year x density band + country + NUTS3 random intercept
```

Controlled:

```text
radiance ~ year x density band + country
         + winter snow anomaly x fixed baseline brightness
         + NUTS3 random intercept
```

Reporting spaces:

1. log10 calibrated radiance, converted to % yr-1;
2. calibrated radiance in absolute nW cm-2 sr-1 yr-1.

The snow x baseline-brightness interaction is pre-specified from the validated
2021 mechanism: the same reflecting-surface change can produce a larger
absolute VIIRS response where more artificial light is available to reflect.

### Pre-registered stopping rule

Historical albedo acquisition is considered scientifically warranted only if:

1. the added snow block is significant at alpha=0.05 in at least one primary
   reporting space; **and**
2. at least two density classes move materially after snow control.

A proportional trend shift is material if it is at least 0.25 percentage
points yr-1, at least 25% of a non-trivial unadjusted trend, or causes a
non-trivial sign reversal.

An absolute trend shift is material if it is at least 0.25% of that class's
median baseline radiance per year, with the same relative/sign safeguards.

If this rule is not met, winter surface state will be retained as an
explanation for the 2021 pulse but will not be pursued as the main explanation
for the broader historical/model discrepancy.
