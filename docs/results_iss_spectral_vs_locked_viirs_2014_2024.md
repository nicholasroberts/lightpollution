# ISS colour shifts versus locked corrected VIIRS radiance: first EcoVis results

**Review date: 21 October 2026. Status: quantitative descriptive results, no spectral correction authorised.**  Reuses the locked 2017+ zero-point-calibrated fixed-NUTS2016 VIIRS series and the *previously prepared* ISS paired RGB spectral data. Does **not** reopen the frozen Jan1 population/radiance/spider/moth PRC analysis.

## Investigator-defined scope: bounded evidence against LED-only dimming

**Decision, 21 October 2026:** This ISS analysis is a **supporting validity / interpretation check**, not another mechanistic-model discovery study. Its sole immediate purpose is to establish, with spectral comparison *plus independent operational evidence*, that **European VIIRS radiance declines cannot all be ascribed to changes from sodium-like lamps to white LEDs and VIIRS's sensitivity to source spectra**. Do **not** use this stage to estimate a percent LED-vs-dimming partition, fit a spectral correction, develop new demographic models, or propose a large new set of tests unless required to assess this modest evidentiary claim.

**Central illustrative example — France:**
- Matched regional 2014→2024 median VIIRS log10 endpoint ratio **−0.3275**, equivalent to approximately **−53%** of the median regional ratio. This is **not** the national area-weighted total.
- Paired ISS regional median Δ(G/R) **+0.0080**; within-France association of regional Δ(G/R) with VIIRS 2014→2024 change **Spearman ρ=+0.018**.
- Thus large dimming in the corrected DNB record coexists with very little observed spectral change in the older ISS comparison. Neither Δ(G/R) nor the epochs alone proves an exact amount of true broadband emission reduction, especially because the post ISS epoch is an irregular 2014–2020 mosaic.
- **Independent support:** French municipal overnight streetlight shut-offs and legal/operational controls are documented by Cerema and legislation; Enedis reports a **34.2% decline in mean streetlighting electrical power around 02:00 from May 2022 to May 2023**, as described in `docs/literature_europe_nighttime_dimming_policy_led_2014_2024.md`. This is independent of DNB spectral response, though not a NUTS3-matched optical flux measurement. The 2022–23 evidence helps interpret **2014→2024**, and must not be backdated into claims restricted to 2014→2020.
- Therefore France is a compelling **counterexample to a universal LED-only account**, not proof that the 2014→2024 DNB decline is *entirely* due to switch-offs.

**Secondary supporting contrasts:** UK has larger Δ(G/R) (+.064) alongside VIIRS declines and independent evidence for mixed LED/operational programmes; Belgium has large VIIRS declines but median Δ(G/R) slightly negative (−.016); Italy shows whitening without corresponding median national-region dimming. These reinforce that a single universal spectral explanation is inadequate, without requiring a quantitative apportionment.

**Paper intent:** include a concise paragraph or compact illustrative comparison in the historical-radiance methods/results/discussion. Preserve the locked zero-point and snow corrections and two biological PRC transfers. Return to full figure/table selection only at manuscript-design stage.

## Source output files received

Investigator uploaded `iss_viirs_between_country_associations.csv` (18 rows), `iss_viirs_sign_sensitivity.csv` (36), `iss_viirs_within_country_associations.csv` (312) and `iss_spectral_shift_vs_corrected_viirs_scatter.pdf` (one page; four panels). Calculated with `src/analyze_iss_spectral_vs_locked_viirs.py` and source `data/processed/iss_nuts3_spectral_shift.csv` from `src/validate_iss_spectral_shift.py`. The **1,302 eligible NUTS2016 regional** ISS/VIIRS matches are NOT the same 1,273-region Jan1 population cohort, and their membership depends on ISS common-pixel and VIIRS-positive mask eligibility. Türkiye is absent from the paper-country ISS inclusion list and therefore absent from these spectral associations.

ISS pre colour epoch **2012–2013** versus **heterogeneous composite 2014–2020**; measured G/R and B/G shifts cannot be assigned specific 2014–2020 or 2014–2024 acquisition dates. VIIRS endpoint datasets use calibrated 2014, 2020, 2024; annual full-series 2021 snow patch is not touched by any tested endpoint, so no need to re-extract 2021 for this particular test. No PRC, population, VIIRS source table or model has been modified.

## European NUTS3 pooled geographical associations

Spearman rho between paired ISS spectral change and NUTS3 log10(corrected VIIRS endpoint ratio) (N=1,302):

| Colour ratio | 2014–2020 | 2014–2024 | 2020–2024 |
|---|---:|---:|---:|
| G/R | −0.2156 | **−0.2980** | −0.2116 |
| B/G | −0.1176 | +0.0334 | +0.1581 |

An increasing G/R is consistent with displacement from red-dominated sodium-like sources towards higher green/red light fractions. B/G, the relatively direct blue/green colour indicator, does **not** independently reproduce the decade-long association. Distinguish the colour-ratio estimands; do not call either a measurement of total blue-photon emission, total luminous flux or exact LED conversion fraction. The PDF scatter shows very broad variation around each monotonic trend.

At a stricter minimum of 50 ISS common valid pixels, the 2014–20 G/R rho shifts to −0.245, B/G to −0.113 (n=1,244); the qualitative difference between metrics remains.

## Between countries: stronger association, but not country-level causal identification

Each country receives one equal-weight observation **after taking median regional spectral shift and median regional VIIRS log change**, not a national area-weighted radiance reduction.

| Minimum eligible NUTS3 regions | n countries | rho G/R vs VIIRS 2014–20 | rho G/R vs VIIRS 2014–24 | rho B/G vs VIIRS 2014–24 |
|---|---:|---:|---:|---:|
| 5 | 23 | −0.528 | −0.462 | +0.282 |
| 10 | 20 | −0.493 | −0.516 | +0.386 |
| 20 | 17 | −0.554 | −0.600 | +0.419 |

The cross-country G/R correlation appears resilient to these minimum-size restrictions, and **is stronger than the pooled regional correlation**. But **within-country** 2014–2020 G/R-vs-VIIRS rho is negative for **10** and positive for **10** of the 20 countries with ≥10 valid regions (five ≤−.2 and four ≥+.2). This contrast is consistent with between-country structural differences explaining part of the pooled geographical association. No spatial-dependence-adjusted p-values or confidence intervals are available. Country medians have differing masks and pixel fractions.

## Country contrasts (within-country signed change associations)

| Country | Regions | median ISS ΔG/R | within-country G/R vs VIIRS rho 2014–20 | rho 2014–24 | median NUTS3 VIIRS log change 2014–24 |
|---|---:|---:|---:|---:|---:|
| FR | 96 | +.0080 | +.044 | +.018 | −.3275 |
| UK | 176 | +.0636 | −.202 | −.031 | −.1145 |
| IT | 110 | +.0440 | +.008 | +.072 | +.0046 |
| NL | 40 | −.0133 | −.003 | +.246 | −.0159 |
| DE | 387 | +.0005 | −.081 | −.263 | +.0214 |
| IE | 8 | +.0940 | −.548 | −.762 | +.0237 |
| BE | 44 | −.0161 | +.052 | −.255 | −.1788 |
| RO | 42 | −.0528 | −.088 | −.298 | +.1524 |
| ES | 57 | +.0199 | −.255 | −.227 | −.0137 |

- **France:** enormous regional-median VIIRS decline 2014–24 (log10 ratio −.3275 ≈ **−53%** median region-specific proportional change), but very small median ISS ΔG/R (+.008), almost zero **within-France** 2014–24 rho (+.018). This is compatible with an important non-spectral operational component, independently documented by municipal curfews and electricity records. It does **not** prove optical total emitted output fell by 53% or attribute a regional fraction.
- **UK:** more positive ISS G/R shift (+.0636) and VIIRS decline (2014–24 median-region ≈ −23%); weak within-country relation in 2014–20 (−.202) and near-zero in 2014–24 (−.031). Mixed LED/operation, not a uniform LED-only explanation.
- **Italy:** positive G/R change (+.044) but median regional VIIRS near-flat to slightly higher by 2024; the hypothesised LED-driven apparent dimming cannot simply be imposed on Italian regions as a whole. Milan-specific historical spectral evidence remains plausible and partly predates 2014.
- **Netherlands:** median ISS G/R slightly *decreases* (−.0133), median NUTS3 VIIRS 2014–24 ≈ −3.6%; this **median-region metric is NOT** the previous area-weighted Dutch radiance country change (nor the dramatic NL333 outlier). Dutch greenhouse lighting change has independent sector/operational evidence but not a country-wide observed positive G/R signature in this ISS mask.
- **Germany:** near-zero median G/R change and slight VIIRS increase; heterogeneity within its many regions.
- **Belgium:** negative median G/R but major 2014–24 VIIRS median-region decline (≈−34%): more evidence that the apparent dimming cannot be generally explained by a *positive* G/R colour change; independent Belgian lighting-operation dates/coverage remain necessary.
- **Ireland:** large median positive G/R but n=8; rho −.762 in 2014–24 is not robust national inferential evidence. Median VIIRS *increased* in these eight matched regions by 2024.
- **Romania:** ΔG/R median negative despite earlier literature emphasizing blueward ISS transition in some Romanian locations. This discrepancy requires reviewing regional masking, colour-ratio methods and actual matched 2012–20 geographical support, not overriding either observation with the literature.

The per-country percent conversions above are `100*(10^(median regional Δlog10 VIIRS)−1)` and therefore **median regional endpoint ratios**; they do not equal area-weighted national radiance changes in the earlier 2014–24 report. Do not mix these quantities.

## Signed-region outcome groups: 0 and ±.02 ISS ratio sensitivity

For ΔG/R > +.02, 2014–20 VIIRS dimming occurred in **320/523 = 61.2%** versus **219/497 = 44.1%** for ΔG/R < −.02. For 2014–24 it occurred in **329/523 = 62.9%** versus **181/497 = 36.4%**. Thus signs support an association but not a mechanism assignment. Crucially the *near-zero* G/R group (|ΔG/R| ≤ .02) also exhibits **167/282 = 59.2%** VIIRS dimming over 2014–24; zero G/R change does not exclude true lamp/dimming changes, mixed spectra or ISS measurement uncertainties. B/G sign contrasts are substantially weaker or reversed at 2014–24.

The cutoff .02 is an exploratory ratio tolerance, not an instrument uncertainty-calibrated threshold. Regional counts are spatially dependent and unadjusted for country composition.

## Important timing test

An association between **historic pre/post ISS colour difference** and **2020–2024** VIIRS change remains (pooled G/R rho −.212; between country at n≥10, rho −.466). Because the colour mosaics do **not observe colour changes specifically during 2020–2024**, this warns that persistent geography, common lighting transition trajectories, baseline brightness or policies may confound a naive LED-spectral explanation. It is a **temporal mismatch control**, not a rigorous null control: regions with earlier LED uptake could plausibly have different later lighting practices.

**Do not** convert these associations into a unique DNB→broadband-visible multiplier; two RGB ratios do not identify true source spectra, upward luminous flux, lamp shielding or skyglow. The script's `−log10(1.5)×Δ(G/R)` is illustrative at *fixed green output*, not an observed physical decomposition. The 2012/13 ISS baseline is also not matched to VIIRS's 2014 start date. The ISS raster mask requires raw annual VIIRS in 2013/2020 above 0.5, so selection may vary by region and light-change direction. Independent European NUTS3 spatial dependence and unequal ISS pixel support require further QA.

## Follow-up limited to manuscript evidentiary sufficiency

No additional regression, country-wide causal classification or spectral correction is required for the immediate claim. If a reviewer requests a stronger independent demonstration, use the cited French municipality/Enedis evidence with appropriate date and geography caveats. Previously listed expanded analyses are **deferred**, not active next steps.

**Defensible bounded conclusion:** some historical European DNB radiance decline is consistent with spectral changes, but the French contrast and independently documented reductions in lighting operation show that an LED spectral-response explanation **cannot be presumed to explain 100% of the observed dimming**. Neither the true optical amount nor a region-by-region partition is estimated.
