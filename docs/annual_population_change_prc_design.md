# Annual population change and predicted spider/moth PRC change

Date: 21 October 2026. This supersedes the earlier question about *why* 0.2–1 and 1–5 nW bands contribute most PRC change. That follows substantially from the locked nonlinear PRC response; the question here is whether **measured year-by-year demographic changes help explain within-region year-on-year changes**.

## Existing data (no new input files needed)

- `data/processed/nuts3_population_density.parquet`: original Eurostat annual `DEMO_R_D3DENS` on its matching NUTS-release years. Use **2014–2020 only** with `nuts_release==2016` (and `matches_nuts_release` when present), so consecutive population differences correspond to the same polygons. After 2020 the year-specific Eurostat NUTS release changes, so naïve 2021–24 annual differences would conflate boundary changes and demography.
- `data/processed/fixed_nuts2016_viirs_2014_2024.parquet`: stable NUTS2016 annual calibrated radiance, same geographical support.
- `data/processed/historical_prc_2014_2024_spider_snow_corrected/prc_nuts3_annual_cumulative.csv` and moth counterpart: approved snow-corrected pixel-to-PRC 2014–2024 panels. As first comparison ends in 2020, no 2021 correction is triggered, but **the standard approved datasets remain the sources**.
- `data/processed/fixed_nuts2016_population_density_reference.csv`: fixed category classification based on 2014–2020 median density.

The user-supplied `radiance_density_joint_contributions.csv` is a baseline density×radiance-band group summary and has no annual population changes; it cannot answer this new question by itself.

## Primary first-difference analysis

For each matched fixed NUTS2016 region `i` and annual transition `t-1 → t` from 2014→2015 through 2019→2020, calculate:

- observed population growth `g_it = 100*(D_it/D_i,t-1 - 1)`
- `ΔR_it`, `ΔPRC_spider,it`, `ΔPRC_moth,it` from respective annual values.

For each of three annual change outcomes separately, compare:

- **M0:** `ΔY ~ country×transition-year fixed effects + log1p(R[t-1]) + square(log1p(R[t-1])) + log(D[t-1])`.
- **M1:** `M0 + g_it`, with NUTS3-clustered standard errors for the population coefficient.

First differencing removes time-invariant differences in region means, including long-lived lighting infrastructure and geography. Country×year controls address broad national policy, economics, instrument/calendar and shared seasonal changes. Lagged radiance and density address saturation, regression-to-mean and settlement stage. Do **not** include concurrent Δradiance in the primary PRC regression because the PRC response is generated from the same radiance pixels and it would eliminate the very population→lighting→PRC pathway under examination.

Report for **all three outcomes** the adjusted coefficient and 95% clustered CI, raw total incremental R² and residual variance improvement (the latter is the proportion of M0 residual sum of squares removed by adding observed g). This is **additional association, not a causal attributable fraction of PRC change**. Country-year fixed effects deliberately remove between-country demographic shifts; report this restriction. Biological units differ across taxa.

## Outstanding tests before causal language

- **Confounding/time:** longitudinal NUTS3 effects are differenced, but time-varying LED conversion, street lighting policies, conflict, agricultural/industrial lighting, migration and economic development can still jointly influence population and radiance.
- **Short-lag limitations:** artificial lighting may adjust years after demographic change; repeat with one-/two-year lags of growth when enough observations remain, and a 2014–2020 first-versus-last interval.
- **Model shape:** assess nonlinear population growth × baseline density/radiance interactions, avoiding a broad exploratory fishing expedition. Validate out-of-sample by held-out regions/countries; in-sample R² gain is not forecast skill.
- **Outliers:** NL333 and greenhouse-dominated Dutch regions have exceptional radiance; inspect country and outlier robustness and physical-versus-percentage PRC responses.
- **2021–2024:** avoid claiming year-on-year demographic attribution for these years until independently harmonised population on fixed NUTS2016 polygons is available. Do not silently mix 2021–2024 NUTS2021/2024 rates.

## Code and outputs

Run:

```bash
git pull --ff-only
python src/fit_annual_population_prc_change.py --self-test
python src/fit_annual_population_prc_change.py
```

Output location `data/processed/prc_population_annual_change/`:

- `population_explanatory_power.csv`: M0/M1 fit and cluster-robust annual population-growth association for Δradiance and both ΔPRCs;
- `annual_matched_population_radiance_spider_moth.csv`: matched region/year measurements, annual differences and fitted population-associated statistical terms;
- `model_summary.json`: reproduction and limitations;
- `annual_population_incremental_explanation.pdf`: juxtaposed incremental residual-variance explained.

After inspecting these results, decide if population growth is explanatory enough to merit a main-text panel. If weak, report that directly, which is also biologically interesting given the lighting-policy/spectral drivers.
