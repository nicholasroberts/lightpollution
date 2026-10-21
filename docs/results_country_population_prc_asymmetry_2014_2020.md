# Country-level demographic change vs spider/moth PRC change, 2014–2020

**Analysis reviewed 21 October 2026.** Input files produced on EcoVis and uploaded by investigator:
- \`country_directional_summary_2014_2020.csv\` (72 rows: 36 country labels incl. ALL × 2 taxa)
- \`country_directional_summary_annual_and_cumulative.csv\` (1008 rows, 2014–2020 and six yearly intervals × 2 threshold approaches × 2 taxa × 36 country labels)
- \`france_turkiye_iceland_scatter_2014_2020.pdf\`
- \`map_population_prc_direction_spider_2014_2020.pdf\`
- \`country_annual_correlations_spider.pdf\`.

## Specific hypothesis tested

Whether regions with predicted PRC **deterioration** disproportionately coincide with population growth, while those with PRC **improvement** do not necessarily coincide with population decline. The contrast must account for national base rates: for example most Turkish regions grow regardless of PRC, whereas many Italian regions shrink regardless of lighting conditions. PRC values are **modelled receptor contrasts derived from VIIRS radiance**, not observed animal navigation outcomes. Satellite LED spectral sensitivity is not yet corrected. All demographic and PRC comparisons below are for **2014→2020**, NOT 2014→2024.

**Cohort:** 1,317 complete fixed NUTS2016 NUTS3 regions. Primary quadrant thresholds: **|6-year population change| >1%**, **|absolute PRC change| >0.001**; regions below either threshold are explicitly neutral/unclassified. Magnitude correlations use continuous values without these sign thresholds. Spider and moth responses are two nonlinear transforms of the same radiance record and are not statistically independent biological replications.

## Main pooled finding: descriptive directional asymmetry

| Measure | Spider | Moth |
|---|---:|---:|
| Population growth (>1%) among ALL regions | 672/1317 = 51.0% | same |
| Population decline (<−1%) among ALL regions | 391/1317 = 29.7% | same |
| Regions with PRC deterioration (<−0.001) | 513 | 589 |
| Regions with PRC improvement (>0.001) | 297 | 384 |
| PRC-loss regions also growing (>1%) | **321/513 = 62.6%** | **355/589 = 60.3%** |
| PRC-gain regions also declining (<−1%) | **78/297 = 26.3%** | **107/384 = 27.9%** |
| Enrichment of growth among PRC-loss relative to all-region growth baseline | **+11.5 pp** | **+9.2 pp** |
| Enrichment of decline among PRC-gain relative to all-region decline baseline | **−3.4 pp** | **−1.8 pp** |

The observations are compatible with **stronger direction-specific coincidence of growth with PRC loss than of decline with PRC gain**; not causal attribution.

For regions with both signed changes, pooled (unstratified) **odds ratio** relating population growth (versus decline) to PRC loss (versus gain): spider **1.458** (95% CI 1.031–2.063, Fisher p=.038); moth **1.485** (95% CI 1.091–2.022, Fisher p=.014). All-region Spearman ρ of population growth percentage and continuous ΔPRC remains weak: spider **−.1867**, moth **−.1958**. Those are different estimands: a modest quadrant preference need not generate strong predictive magnitude correlation.

**Country-stratified contingency sensitivity**, computed using the Mantel–Haenszel pooled odds ratio (via statsmodels \`StratifiedTable\`) on country tables with nonzero row/column support:

- At **1%/0.001** thresholds, spider **OR 2.636** (CI 1.285–5.408; p=.014; 13 informative country strata), moth **OR 2.159** (CI 1.172–3.979; p=.021; 15 strata).
- At **zero-only sign** thresholds, spider **OR 1.598** (CI 1.137–2.245; p=.0069; 21 strata), moth **OR 1.410** (CI 1.002–1.984; p=.0488; 21 strata).
- These values are computed from uploaded *aggregate* country 2×2 tables, not underlying NUTS3 data; formal spatial autocorrelation-adjusted inference is **not available**. Underlying countries with nearly all regions having one PRC sign (especially Türkiye, all PRC loss) **do not contribute to a within-country sign odds ratio**. Estimates are sensitive to thresholds and small country cells; do not overstate formal statistical significance.
- Fisher tests *within any individual country* at main thresholds have **no p<.05**; some contingency tables are degenerate because one PRC direction is absent.

**Magnitude-conditioned pooled Spearman** (exploratory, restricted outcomes): for PRC-loss regions, population growth versus loss magnitude spider **+.329**, moth **+.340**; for PRC-gain regions, population **decline magnitude** versus PRC gain spider **−.171**, moth **−.141**. These unadjusted magnitude correlations are not country-specific effects or causal tests and may partly reflect conditioning on outcome sign.

## The illustrative national contrasts

| Country | NUTS3 n | Pop growth >1% | Pop decline <−1% | Spider PRC loss (>0.001 abs) | Spider PRC gain | Spearman ρ (spider, moth) |
|---|---:|---:|---:|---:|---:|---:|
| Türkiye (TR) | 81 | 74 | 3 | **81** | 0 | **−.286, −.295** |
| France (FR) | 101 | 48 | 20 | **1** | **86** | **+.193, +.214** |
| Italy (IT) | 110 | 10 | **87** | **59** | 15 | **−.369, −.348** |
| Germany (DE) | 401 | 304 | 39 | 212 | 49 | −.154, −.163 |
| Netherlands (NL) | 40 | 16 | 9 | 17 | 9 | +.047, +.066 |
| Iceland (IS) | **2** | 2 | 0 | 1 | 0 | **not estimable** |

- **Türkiye:** 74/81 regions both grew and lost PRC (remaining PRC-loss regions: three with population decrease, four with near-stable change); greater growth magnitude correlates moderately with greater losses. However all 81 regions lost PRC, and 74 were growing regardless, so sign coincidence alone is non-informative within Türkiye; use ρ, not the 74/81 count, as evidence of within-country *magnitude* association.
- **France:** 41 regions had **population growth AND spider PRC gain**, 14 had population decline AND PRC gain; moth analogous **44 vs 18**. Robust qualitative example that PRC gain does **not require demographic decline**. Spearman is actually positive; the France scatter displays an isolated extreme +~20% population change / negative PRC case that likely reverses Pearson's sign (spider Pearson −.155 versus rank +.193). Do not interpret FR Pearson without outlier QA.
- **Italy:** 43/110 regions had **population decline AND spider PRC deterioration** (moth 45), whereas just 9 had growth and PRC deterioration. Indicates PRC loss does **not require population growth**, though within-country continuous rank correlation is negative, and why simple national sign stories can mislead.
- **Germany:** widespread growth and PRC loss but more mixed categories, weak rank relation and two-taxon contrast.
- **Iceland:** just two matched NUTS3 regions (the other territory lacking matching population) → cannot estimate meaningful country correlation or use as an inferential example.

## CRITICAL annual Eurostat demographic QA problem

Inspecting \`country_directional_summary_annual_and_cumulative.csv\`: for the **2019→2020** transition, **29 of 35 countries (excluding ALL) have NO NUTS3 population change exceeding ±0.25%** — collectively **1,033 of the 1,317** matched regions. **31/35 countries** have undefined Spearman correlation that transition (includes small country n). This conspicuous departure from 2014→15 through 2018→19 is a **MAJOR potential source discontinuity/carry-forward/fallback/boundary/rounding problem**. The tables do NOT prove the 2020 data were copied from 2019, because values below 0.25% would produce the same categorical result; but \`NaN\` ρ for many countries indicates insufficient numerical variation and warrants inspecting full region-level values, \`population_density_source\` and source release. Example: Germany's 401 regions, Türkiye's 81, Poland's 73 and Romania's 42 have neither positive nor negative >0.25% region-year changes in 2019→20, with undefined correlations.

**Before using any 2014→2020 demographic association as a publishable main result**: examine the 2019 and 2020 NUTS3 yearly population values and fallback/source fields, quantify exact equalities/zero variance, and re-run the exact **2014→2019** cumulative endpoint as a preregistered sensitivity using the same geographic cohort. Do NOT silently infer 2020 population or claim 2014→2020 means true six-year population change if the 2020 data are carried forward. Also inspect scatter outliers +~80% annual (2015–16) and +~300% (2019–20).

Country-annual spider Spearman heatmap is patchy, with signs varying even within the same countries: TR mildly positive 2014–15, negative 2015–16 and 2016–17, near zero subsequently; FR mostly positive or mixed; IT mixed. A single cumulative country result cannot imply a stable year-by-year association.

## Provisional manuscript wording

> Population growth was more frequently observed among regions showing a loss of modelled polarization-receptor contrast than among European NUTS3 regions overall, whereas population decline was not correspondingly overrepresented in regions gaining contrast. However, strong national differences demonstrate that demographic change alone is insufficient to explain contrast trends: widespread deterioration coincided with growth in Türkiye, but PRC recovery often accompanied growth in France, while deterioration also accompanied population decline in Italy. These associations are observational, reflect VIIRS-band light and two calibrated photoreceptor models, and need verification against the underlying 2020 population data and independent lighting changes before claims of causation.

**Next action:** source QA and 2014→2019 sensitivity first; then ascertain whether country-stratified asymmetry and Türkiye/France/Italy contrast survive that test. No new model-fitting needed to answer this.
