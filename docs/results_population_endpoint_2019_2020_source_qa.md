# Eurostat 2020 source audit and fixed-cohort PRC endpoint robustness

**Date:** 21 October 2026. **Status:** EcoVis outputs reviewed; population-source anomaly unresolved. This report is a sensitivity analysis, not a refit or replacement for the locked spider/moth PRC functions, nor a causal assignment of lighting policy. Applies to the *demographic subanalysis*, not a change to the 2014–2024 historical-PRC study period.

## Inputs and reproducibility

Executed `src/audit_population_2020_and_2019_prc_sensitivity.py` against the seven-year complete fixed-NUTS2016 panel (2014–2020). **Same 1,317 NUTS3 regions** for both 2014→2019 and 2014→2020; original spider17/moth16 locked snow-corrected PRC products, original Eurostat `demo_r_d3dens` regional population-density series. Predictor = **absolute change in people/km²**; response = **absolute PRC change**; Spearman within countries; existing descriptive thresholds |Δdensity|>5 people/km², |ΔPRC|>.001.

Verification against existing 2014–2020 outputs: **all 72 old country/taxon correlations reproduced** within tolerance. Thus change in the 2019 endpoint cannot be explained by a silently altered NUTS3 cohort or transfer function.

Uploaded EcoVis results: `country_2014_2019_vs_2020_correlations.csv`, `population_2019_2020_change_qa_by_country.csv`, `eurostat_source_status_by_country_year.csv`, `ireland_region_paired_2014_2019_2020.csv`, `qa_summary.json`. Source outputs remain local to EcoVis at `data/processed/prc_population_endpoint_qa_2014_2019_2020/` and are **not** added to Git (large datasets and generated figures are handled separately).

## 2020 demographic observation problem

- **1,113 / 1,317 (84.51%)** complete-cohort NUTS3 regions report **identical population density in 2019 and 2020**.
- **29 of 35** sampled countries have every matched region's 2020 density exactly equal to 2019 (a *zero* 2019→2020 change); examples **DE 401/401, TR 81/81, IE 8/8**, also BE, CH, PT, PL, RO.
- Other patterns differ: FR 18/101 identical, IT 6/110, NL 38/40, AL 0/12. One Netherlands region has an anomalous **−432.0 people/km²** 2019→2020, and Albania's region-level changes range **−289.6 to +207.3 people/km²**. Audit those source observations rather than immediately imputing/correcting them.
- These are carried in the supplied dataset as original `demo_r_d3dens`, not our population-density fallback: **none** of the matched 2014, 2019, 2020 observations are recorded as `population_density_recovered`. Most 2020 records have no status reported; all 110 Italian 2020 entries are status **b** (Eurostat break in series), Romania's 42 2020 entries status **e** (estimated), Albania's 12 also **e**, and 2 German entries **e**. The absent status flag elsewhere does **not** establish observed independent 2020 measurement.
- Widespread perfectly copied density values cannot reasonably be described as true zero population change. Their provenance—possible retrospective carry-forward, publication practice, NUTS vintage or pipeline mapping—**has not yet been established**. Source records and an independently computed population/land-area check are required.

**Decision for interpretation:** use **2014→2019 as the provisional primary demographic comparison**, retaining 2014→2020 only as a clearly labelled sensitivity. Do not write manuscript text claiming independently observed 2019→2020 demographic stasis. This does NOT invalidate 2020 PRC/radiance values; the satellite series and demographic source are separate.

## Endpoint robustness: within-country Spearman ρ

| Country | Matched n | Spider 2014–19 | Spider 2014–20 | Moth 2014–19 | Moth 2014–20 |
|---|---:|---:|---:|---:|---:|
| All | 1,317 | −0.160 | −0.169 | −0.158 | −0.178 |
| Türkiye TR | 81 | **−0.492** | **−0.566** | **−0.507** | **−0.582** |
| Germany DE | 401 | **−0.394** | **−0.365** | **−0.400** | **−0.375** |
| France FR | 101 | +0.163 | +0.226 | +0.201 | +0.263 |
| Italy IT | 110 | −0.129 | −0.092 | −0.132 | −0.093 |
| Ireland IE | 8 | +0.048 | +0.286 | +0.048 | +0.286 |
| Belgium BE | 44 | −0.055 | +0.064 | −0.032 | +0.068 |
| Portugal PT | 25 | +0.299 | +0.319 | +0.259 | +0.331 |
| Romania RO | 42 | −0.151 | +0.065 | −0.146 | +0.053 |

These show robust negative *descriptive rank* associations in TR and DE, but country-specific heterogeneity elsewhere. France's 2014–19 ρ is positive yet modest; no independent evidence that municipal switch-offs produced the regional association. Italy remains weak. Several countries (e.g. Belgium, Romania) change sign between windows: mechanisms and classification must not be inferred from 2020 alone. Both taxa use the same underlying satellite radiance; agreement across taxa is not independent confirmation.

At the prespecified 2014–2019 main thresholds, spider PRC **loss** occurs in 79/81 TR regions, while PRC **gain** occurs in 94/101 FR regions. DE has 125 PRC-loss and 104 PRC-gain regions at threshold; IT 42 losses and 26 gains. **These counts are not statistical tests or estimates of causality.**

## Ireland: Category B (lighting-transition comparator), with quantified biological contrast

All **8/8** matched Irish NUTS3 regions have a **positive 2014→2019 spider PRC change** and a **positive moth PRC change**, while all eight also show some population-density increase. At |ΔPRC|>.001, 4 spider regions and 6 moth regions have a detectable positive change. **Two** Irish regions have |Δdensity|>5 people/km² (both increase), and just one of these clears the spider PRC +.001 threshold. For 2014→2020, 7/8 have positive spider and moth ΔPRC; the remaining IE053 decrease is tiny (~−0.000043 spider, −0.000054 moth), well below the chosen threshold.

The within-Ireland Spearman ranks for density increase versus absolute ΔPRC are **ρ=+0.0476** (both taxa, 2014–2019, n=8) and **ρ=+0.2857** (2014–2020). **Do not describe these as a demonstrated positive density–PRC association.** 2014–2019 spider Pearson is +0.969 but is heavily dominated by **IE061**, whose density rose 111.8 people/km² and spider PRC +0.005439. Excluding IE061 changes the 2014–19 spider Pearson to **−0.527** and Spearman to **−0.429** (n=7); this illustrates leverage, not a robust inverse relationship. Ireland is useful because modelled PRC improves while density grows, not because regional growth explains improvement.

SEAI documents time-aligned local LED conversion in Monaghan/Cavan (see `docs/country_mechanism_typology_population_prc_2014_2020.md`); **no spatial NUTS3 attribution** or causal effect of these programmes on inferred PRC is established. The VIIRS spectral response caveat applies especially strongly to LED transitions.

## Interpretation and next decision

1. Keep **2014–2024** as the biological PRC/radiance paper horizon. For population-change comparisons provisionally prioritise **2014–2019**, on a fixed complete 2014–2020 cohort; call 2014–2020 a qualified sensitivity until 2020 records are independently corroborated.
2. Audit the Eurostat **2020** source with independent 2019/2020 January-1 population counts (`demo_r_pjanaggr3`) and properly matched land-area denominators; check NUTS 2016 vintage, metadata and publication status before asserting a cause or rewriting the source. Inspect the exceptional **NL** and **AL** deltas separately. Do not silently replace official data.
3. Retain TR/DE negative rank pattern and their **different plausible mechanisms**. France and Ireland show gains without population decline; IT weak rank association; regional NUTS3 patterns and policies must be spatially/time matched before explaining causes.
4. Do not treat unadjusted country correlations as independent significance tests: spatial dependence, different country sizes and baseline differences remain. Do not interpret spider/moth concordance as satellite replication.
