# Quantitative European radiance vs spider and moth PRC (2014–2024)

**Analysis date:** 20 October 2026. **Status:** descriptive results from investigator's seven EcoVis output files; LED spectral correction **not yet applied**. Snow-adjusted 2021 VIIRS pixel patch and both empirical biological transfer functions are **locked**.

## Source files

- `figures/prc_radiance_spider_moth_2014_2024/nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv`
- `data/processed/historical_prc_2014_2024_spider_snow_corrected/prc_europe_annual_cumulative.csv`
- `data/processed/historical_prc_2014_2024_moth_snow_corrected/prc_europe_annual_cumulative.csv`
- both taxa's `prc_nuts3_annual_cumulative.csv` corrected panels
- both taxa's `snow_correction_summary.json`

The submitted files were named with browser suffixes `(2)` or `(3)`. Their identity was **verified from actual response magnitudes and consistent NUTS3 endpoint matching**, rather than guessed by suffix: spider starts near PRC .4381 and moth near .6589.

## Data/QA

Fixed NUTS2016 NUTS3: **1,522 regions**, **11 years**, **16,742 rows per taxon**; both report `snow_status: locked` and all 1,522 regional `coverage_stable=true`. Original and corrected spider/moth 2014 and 2024 values exactly match their values in the combined NUTS3 CSV (max |difference| 0). Identical valid-pixel area per taxon to floating-point rounding. The radiance CSV marks **57** NUTS3 cases with unreported percentage radiance because 2014 regional mean VIIRS radiance below 0.2 nW cm^-2 sr^-1. Region endpoint-area coverages are close but not necessarily exactly the same valid pixels: full common-pixel masking is not guaranteed. Both 2014→2024 comparisons are unchanged by the 2021-only correction.

## NUTS2016 area-weighted Europe

| Metric | 2014 | 2024 | Δ absolute | Δ percent |
|---|---:|---:|---:|---:|
| Corrected-band VIIRS radiance (nW cm^-2 sr^-1) | 1.065270 | 1.079863 | +0.014593 | **+1.370%** |
| Spider mean PRC | 0.438102 | 0.436695 | −0.001407 | **−0.321%** |
| Moth mean PRC | 0.658949 | 0.656985 | −0.001964 | **−0.298%** |

This is the **area-weighted mean radiance**, not the arithmetic mean of region-specific percent changes, nor any measurement of total human lighting or biologically relevant spectrally complete light. Regional mean PRC is first computed from the nonlinear response at **individual native VIIRS pixels**, then weighted by each NUTS3 polygon's valid area. Do not interpret the small net mean PRC decline as proof that all habitats suffered equivalent losses.

Among individual NUTS3 regions, **805 increased radiance, 717 decreased**; **944 spider PRCs decreased, 578 increased**; **922 moth PRCs decreased, 600 increased**. Their regional medians are Δspider PRC −0.000933 and Δmoth PRC −0.001270. This differs from European area-weighted means because highly populated/bright NUTS3 are not geographically equally sized and regional extrema differ.

### 2021 locked snow/albedo adjustment

Spider original European 2021 PRC **0.435918398** → corrected **0.436621483**, Δ **+0.000703085**; flank 2020/22 reference mean **0.437367112**; correction removes **48.53%** of the absolute original deviation. Moth original **0.655603754** → corrected **0.656737846**, Δ **+0.001134091**; 2020/22 reference mean **0.657876097**; correction removes **49.91%** of anomaly. Residual 2021 discontinuity remains. Flank comparison is QA, **not model fit target**. Both taxon's QA summaries identify correction locked and all other years unchanged. Original global annual VIIRS rasters preserved.

## Country contrasts (fixed NUTS2016, area-weighted by covered km²)

| Country | NUTS3 n | Δradiance % | Δspider PRC % | Δmoth PRC % |
|---|---:|---:|---:|---:|
| France | 101 | −47.72 | +1.94 | +2.07 |
| UK | 179 | −23.46 | +0.97 | +1.15 |
| Netherlands | 40 | −39.80 | **+0.034** | **+0.280** |
| Belgium | 44 | −27.09 | +5.36 | +5.80 |
| Germany | 401 | +3.55 | −0.28 | −0.27 |
| Switzerland | 26 | −23.88 | +1.07 | +1.12 |
| Italy | 110 | +1.59 | −0.81 | −0.78 |
| Türkiye | 81 | +79.70 | −2.30 | −2.37 |
| Greece | 52 | +31.73 | −1.25 | −1.28 |

Country percentage PRC is `100*(sum_i area_i*PRC_i,2024/sum_i area_i) / (sum_i area_i*PRC_i,2014/sum_i area_i)-100` using identical 2014 region area weights. Country percentage radiance uses same weighting. The 2016 NUTS geography includes outlying territories; these are *NUTS2016 area-weighted* aggregates, not strictly contemporary EU or mainland-only indicators. The Netherlands' tiny area-averaged net PRC response results from highly heterogeneous spatial and countervailing regional effects and must not be summarized as zero PRC benefits everywhere.

## Regional association is strong in rank, not proportional in magnitude

For all **n=1,522** NUTS3 regions, absolute ΔR vs spider ΔPRC: **Pearson r −0.325**, **Spearman ρ −0.923**; ΔR vs moth ΔPRC: **r −0.375**, **ρ −0.936**. Strong negative rank association is EXPECTED from monotone PRC(R), but magnitude response is nonlinear, and region-mean R changes conceal pixel distributions. Extreme R values attenuate Pearson correlation: excluding top 1% |ΔR| increases absolute Pearson to spider **−0.738**, moth **−0.799**. Spider–moth absolute ΔPRC: Pearson **+0.9903**, Spearman **+0.9985**, so species differ primarily in **magnitude**, much less in regional rank ordering.

Some regions also reverse the simple expected signs because averaging a nonlinear pixel-wise response does not commute with averaging regional radiance: **163/1,522** regions show spider ΔPRC with same sign as regional mean ΔR (151 radiance decreases with spider PRC decreases; 12 radiance increases with spider PRC increases); **139/1,522** analogous moth cases (128,11). Most differences are small: requiring |ΔR|≥0.2 and |ΔPRC|≥0.005 leaves **6 spider** and **9 moth** cases; investigate their source/masking before presenting behavioural inversions.

## Grounded regional exemplars (2014–2024)

| NUTS3 | ΔR (nW) | ΔR % | Δspider PRC | Δspider % | Δmoth PRC | Δmoth % |
|---|---:|---:|---:|---:|---:|---:|
| NL333 Delft en Westland | −179.536 | −65.08 | +0.024667 | +15.72 | +0.060051 | +44.05 |
| UKC23 Sunderland | −15.299 | −53.15 | +0.066936 | +33.49 | +0.139200 | +58.84 |
| TR100 Istanbul | +5.309 | +48.03 | −0.037560 | −10.42 | −0.064370 | −12.28 |

**Do not claim NL333 shows no benefit to moths**: moth model predicts +44% relative PRC, versus +15.7% for spider. It is the greatest numerical ΔR outlier and earlier flagged as anomalous in the radiance project. Its apparent 2014 regional mean raw/corrected R ~275.87 and 2024 ~96.33 deserve historical source/lighting (greenhouse) and spectral QA before being used as a real-world emission claim. NL333 is only **0.414%** of the Dutch covered geographical area, yet contributes **25.7% of the Dutch area-weighted ΔR**. Dutch country aggregate excluding NL333 still has **−35.08% ΔR, +0.008% spider PRC, +0.238% moth PRC**. Netherlands remains a strong example of area-aggregated radiance change not converting linearly into national mean PRC gain, but not a claim about every habitat.

**Jensen's inequality / ecological aggregation:** calculating `PRC(mean(R))` directly at regional mean R is WRONG. At NL333, spider PRC of means would predict ΔPRC about **+0.00001** and moth **+0.00031**, whereas correct pixel-transform→area-aggregate gives **+0.024667** and **+0.060051**. Heterogeneous dimmer pixels continue to affect PRC even when average R is dominated by very bright installations. This example is also susceptible to exceptional 2014/2024 source radiance and must not alone support the biological conclusion.

## Baseline NUTS3 2014 mean-radiance groups (not population density)

| 2014 mean R nW cm^-2 sr^-1 | NUTS3 regions | % land area | share of area-integrated |Δspider PRC| | share of area-integrated |Δmoth PRC| |
|---|---:|---:|---:|---:|
| <0.2 | 57 | 16.09 | 3.31% | 3.22% |
| 0.2–1 | 698 | 57.58 | 50.95% | 49.61% |
| 1–5 | 510 | 23.66 | 37.43% | 37.78% |
| 5–20 | 192 | 2.45 | 7.46% | 8.23% |
| ≥20 | 65 | 0.23 | 0.85% | 1.16% |

**88.4% (spider) and 87.4% (moth)** of `sum_NUTS3(area_km2 × |ΔPRC|)` occurs in regions whose initial mean radiance is **0.2–5 nW**, which cover ~81.2% of the study area. This is a purely **area-integrated modelled contrast-change measure**: not the number of affected animals, an ecological exposure estimate, or a direct link to human population density. Need separate population gridded data before testing proposed density explanation. Largest *per-region* |ΔPRC| occurs mostly in brighter cities where |ΔR| also tends to be large; don't conflate effect size with land-area contribution.

## Interpretation guardrails

1. Positive area-weighted DNB radiance change (+1.37%) plus negative area-weighted PRC contrast change (~−0.3%) is possible because of **nonlinearity, initial radiance and spatial aggregation**. The country pattern includes large west-European DNB declines and southeastern-European increases. No claim of paradox or monotonic transformation failure at the individual pixel: every species PRC function remains decreasing in R.
2. **LED spectral differences have NOT been corrected**: current outputs are post-2017 baseline-calibrated DNB, with locked 2021 snow correction, not total visible/scotopic/skyglow trends. Verify ISS spectra before causal claims.
3. All NUTS3 regions are coverage-stable across annual products but exact common-valid-pixel temporal support is not established. 2014 region annual products may have processing biases. Full regional common-pixel support and NL333 extreme QA must precede submission.
4. European mean of area-weighted region values is an approximation to spatial average of native-pixel contrast in study geography. Region population density does not imply the abundance/distribution of either nocturnal arthropod taxon.
5. Empirical receptor contrast is a predicted cue, **not directly measured behavioural detection or navigation success**. Moth/spider percent changes have different PRC baselines and should not be interpreted as species fitness changes.

**Calculation basis:** uploaded outputs analysed using Pandas/NumPy and Pearson/Spearman correlations, with geographical weights from `valid_area_km2_2014`. Numerical country shares are descriptive and lack sampling/model uncertainty intervals.
