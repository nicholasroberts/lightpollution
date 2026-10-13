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
