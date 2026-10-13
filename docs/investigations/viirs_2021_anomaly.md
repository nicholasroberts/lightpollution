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
| Snow/albedo test — **pending clean rerun** | After a valid rerun: `figures/population_relationship/fixed_nuts2016_snow_viirs_test/february_snow_vs_viirs_anomaly.pdf`; `snow_viirs_anomaly_rho_jan_feb_mar.pdf`; `levin_style_jan_july_snow_radiance_2021.pdf`; `levin_style_percent_lit_thresholds.pdf` | After a valid rerun: `data/processed/fixed_nuts2016_snow_2021_anomaly_tests.csv`; `fixed_nuts2016_snow_2021_anomaly_models.csv`; `fixed_nuts2016_snow_levin_jan_july_summary.csv`; `fixed_nuts2016_snow_levin_percent_lit_thresholds.csv` | These are **not yet paper-eligible**. The first run was invalid because categorical MOD10CM codes entered the zonal means. If the corrected February association survives the pre-specified controls while March remains null/weak, this will become the direct mechanism test tying the monthly localization to surface snow/albedo. |

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
7. the snow/anomaly mechanism figure, **only after the corrected MOD10CM rerun is validated**.

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
| OPEN | Correctly masked snow explains February 2021. |
| OPEN | Snow effect survives country/latitude/brightness/cf_cvg controls. |
| OPEN | Snow x brightness is positive as reflection physics predicts. |
| OPEN | Residual signal after snow requires additional compositing/calibration explanation. |
| OPEN | NDVI and albedo should be added after the clean snow result. |

---

## Restart instructions

From repository root on EcoVis:

```bash
git pull
python src/run_pipeline.py future-radiance-snow-2021-test
```

Do not interpret snow until:

1. the integration self-test passes;
2. the stage finishes with all snow means in 0-100;
3. the corrected February/January/March results are inspected.

Files to return after the clean run:

```text
data/processed/fixed_nuts2016_snow_viirs_test_summary.txt
data/processed/fixed_nuts2016_snow_2021_anomaly_tests.csv
data/processed/fixed_nuts2016_snow_2021_anomaly_models.csv
data/processed/fixed_nuts2016_snow_levin_jan_july_summary.csv
data/processed/fixed_nuts2016_snow_levin_percent_lit_thresholds.csv

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
```
