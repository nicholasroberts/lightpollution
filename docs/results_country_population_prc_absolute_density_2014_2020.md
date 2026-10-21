# Absolute population-density change versus spider/moth PRC: revised country results

**Reviewed 21 October 2026**, directly from EcoVis `absolute_vs_percentage_correlation_sensitivity.csv`, `country_directional_summary_2014_2020.csv`, and `country_annual_correlations_spider.pdf`. This document **SUPERSEDES the percent-population-primary directional interpretation** in `docs/results_country_population_prc_asymmetry_2014_2020.md` (which remains as historical sensitivity). The **primary demographic change variable is Δpopulation density measured in people/km²** over the fixed NUTS2016 NUTS3 area; the response is **absolute ΔPRC**. Percentage population change is retained only as comparative sensitivity.

## 2014→2020 within-country continuous association (Spearman ρ)

| Country | n | Spider absolute-density ρ | Moth absolute-density ρ | Spider percentage-density ρ |
|---|---:|---:|---:|---:|
| ALL | 1,317 | −.169 | −.178 | −.187 |
| Türkiye TR | 81 | **−.566** | **−.582** | −.286 |
| Germany DE | 401 | **−.365** | **−.375** | −.154 |
| France FR | 101 | **+.226** | **+.263** | +.193 |
| Italy IT | 110 | **−.092** | **−.093** | −.369 |
| Netherlands NL | 40 | +.102 | +.129 | +.047 |
| Belgium BE | 44 | +.064 | +.068 | −.079 |
| Switzerland CH | 26 | +.131 | +.141 | +.265 |

**Key finding:** absolute density change substantially strengthens Türkiye and Germany correlations compared with percentage population growth, while Italy's formerly strong −.369 rank correlation nearly disappears (−.092). Percent-based country claims for Italy are not robust. France shows a positive rank association, against a simple growth→PRC loss mechanism. Annual cross-country strength is inconsistent; the same VIIRS radiance is the common input to both taxa, so their concordance is not independent evidence.

## Main signed-change threshold results

Main descriptive thresholds: |Δdensity2014–20| > **5 people/km²**, |absolute ΔPRC| > **0.001**; other cases neutral or mixed.

| Statistic | Spider | Moth |
|---|---:|---:|
| Population growth among 1,317 complete regions | 404 / 1,317 = 30.7% | same |
| Population decline among 1,317 | 103 / 1,317 = 7.8% | same |
| Population growth among PRC-loss regions | 223 / 513 = 43.5% | 244 / 589 = 41.4% |
| Population decline among PRC-gain regions | 21 / 297 = 7.1% | 28 / 384 = 7.3% |
| Both change magnitude thresholds met | 387 | 438 |
| Neutral/mixed (at least one threshold fails) | **930** | **879** |

Raw proportions exhibit enrichment of absolute-density growth among PRC losses but **no corresponding enrichment of absolute-density decline among PRC gains**. HOWEVER: these two descriptive contrasts do not imply a formally strong pooled association once neutral cases are removed: signed 2×2 Fisher odds ratios spider **1.007, p=1.000**, moth **1.119, p=.695**.

A **country-stratified Mantel–Haenszel** test on 8 informative countries (AL, BG, DE, ES, FR, IT, NL, PL), derived from uploaded country 2×2 counts, yields:
- Spider **OR=1.219** (95% CI **0.404–3.685**), p=.751.
- Moth **OR=1.049** (95% CI **0.393–2.799**), p=.929.

This **does NOT support a broadly general within-country 2×2 asymmetry** at the main absolute-density thresholds. It is also underpowered/restricted because 930/1,317 spider regions are neutral/mixed and 8 country strata contribute to the test. Countries such as Türkiye have all PRC observations on the deteriorating side, hence no within-country sign odds ratio; use the continuous correlation there.

Country-specific main-threshold signs:
- **Türkiye TR**: n=81, ALL 81 spider/moth PRC deteriorations >.001; 30 region density increases >5 people/km², none declines below −5; 51 insufficiently large density changes. Thus 30/81 *coincidences* are expected given growth base rate; separate **Spearman −.566/−.582** identifies rank association with magnitude.
- **France FR**: n=101, 86 spider PRC improvements vs 1 loss; 19 significant density increases (>5 people/km²), 4 declines (<−5); **14 density↑ and spider PRC↑**, four density↓ and PRC↑. Strong evidence that French PRC improvement is possible with population growth; not that all changes are LED/policy caused.
- **Italy IT**: n=110, 59 spider PRC losses, 15 gains; 49 density declines, 6 increases; **25 density↓ and PRC↓**, versus five density↑ and PRC↓. The old Italy percent-density correlation of −.369 is not robust to the better density-magnitude metric; avoid emphasizing demographic increases as an Italian explanation.
- **Germany DE**: n=401, 212 spider PRC losses, 49 gains; 224 significant increases in density, five declines; stronger continuous negative rank association but most density declines absent, so 2×2 inference limited.

## Annual correlations and essential 2020 QA

Spider within-country Spearman ρ, annual Δdensity (people/km²) versus absolute ΔPRC:
- TR: −.02, −.30, −.50, −.38, −.17, **undefined** (2014–15 ... 2019–20)
- DE: −.23, −.02, −.29, −.15, −.09, **undefined**
- FR: +.32, +.09, −.09, +.07, −.04, +.19
- IT: −.01, −.11, +.02, −.25, −.08, +.05

Important: **31/35 countries (excluding ALL) have undefined Spearman in 2019→2020**, including DE (n401) and TR (n81); 1,054 matched regions belong to countries with undefined 2019→2020 correlation. This suggests a possible 2020 source/carry-forward issue, or insufficient demographic variation; aggregate files cannot resolve whether individual 2019 and 2020 density values are identical. **Do not claim verified real demographic changes for 2019–2020 without examining the original year-by-year Eurostat density source and fallback status.** Re-run the absolute-density cumulative comparison with **2014→2019** as a prespecified robustness check.

Other caution: NUTS2016 geographic region completeness n1317/1522; pixel-level PRC but NUTS3-level population; ecological receptor prediction not observed animal outcomes; VIIRS LED spectral artefacts uncorrected; spatial autocorrelation makes simple region-independent p-values optimistic; country relationships are observational and cannot establish that demographic change *caused* PRC changes. A +5 people/km² threshold has different relative meaning in cities and countryside—report continuous associations alongside sign thresholds.

## Manuscript status

The robust, defensible conclusion **so far** is **country-specific heterogeneity**, not a pan-European demographic causal law. Particularly TR and DE show stronger negative continuous associations when population changes are measured in people/km²; French contrast recovery clearly does not require demographic decline. The proposed general gain/loss asymmetry may be descriptive but is **not significant** under the updated country-stratified signed analysis. First perform population-source QA and 2014→2019 endpoint sensitivity; only then settle manuscript conclusions.
