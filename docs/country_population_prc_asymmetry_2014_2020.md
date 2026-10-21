# Country-specific asymmetric population growth versus PRC loss/gain

**21 October 2026.** Main question: among the fixed-NUTS2016 regions in which spider/moth PRC **deteriorates**, is measured population **growth** disproportionately frequent and correlated with the loss? Conversely, among regions with PRC **improvement**, is it *not* similarly accompanied by population **decline**? France, Iceland and Türkiye are proposed contrasting examples, but analyse **all countries first**.

## Avoiding a misleading directional interpretation

Two distinct comparisons are needed:

1. Country-specific 2014→2020 **Pearson/Spearman correlations**, comparing **percentage population density change** and **absolute ΔPRC** on exactly the same regions; assess annual transitions separately.
2. Country-specific **four-quadrant sign table**: pop up / PRC down; pop down / PRC up; pop up / PRC up; pop down / PRC down. Report conditional fractions with denominators (e.g. `P(population growth | PRC loss)` and `P(population decline | PRC gain)`) but crucially **also report the baseline proportion growing** and Fisher's exact odds ratio. Otherwise positive associations are inevitable in growing countries like Türkiye even if lighting changed independently of population.

A separate correlation within PRC-loss regions tests **magnitude** of change, using growth (%) against PRC loss magnitude (−ΔPRC); a counterpart in PRC-gain regions uses population decline (−Δpop) against PRC gain (+ΔPRC). These have restricted ranges and must not be mistaken for unbiased causal regressions; small countries often have too few regions for such correlations.

To suppress tiny numerical changes near zero, main quadrant results require **at least 1% cumulative population change and 0.001 absolute PRC**, or **0.25% annual population change and 0.0005 annual absolute PRC**. Raw zero-only sign classifications are also exported to test threshold dependence. Thresholds define descriptive bins, not biological visibility thresholds. Primary separate spider and moth PRC functions remain locked.

## Prior data and provenance

The output is computed on the earlier map analysis's **same complete 1,317 NUTS3-region × 7-year cohort**, using `src/visualize_annual_population_prc_changes.py` as its sole data-loading/matching source. It uses observed NUTS2016 Eurostat population density 2014–2020 and both locked 2021 snow-corrected PRC annual panels (2021 falls outside the observed demographic window). No files need re-uploading. Data/processed artefacts are locally generated on EcoVis and deliberately not committed to Git.

**No country-specific conclusions can be claimed yet:** the national joined NUTS3 observation file is not in GitHub, and the user has uploaded the *Europe-wide summary correlations*, not the underlying region-level 2014→2020 population-change panel. The national asymmetry remains a **hypothesis** until the script runs on EcoVis. Country size/coverage matter; small NUTS3 samples (e.g. Iceland) do not permit confident within-country correlation estimates.

## Run on EcoVis

```bash
git pull --ff-only
python src/analyze_country_population_prc_asymmetry.py --self-test
python src/analyze_country_population_prc_asymmetry.py
```

Produces `figures/prc_population_country_asymmetry_2014_2020/`:

- `france_turkiye_iceland_scatter_2014_2020.pdf`: 3 × 2 country panels, PRC vs population changes, spider and moth; explicit quadrants and actual region counts.
- `map_population_prc_direction_spider_2014_2020.pdf` and moth: Europe coloured by population-up/PRC-down, population-down/PRC-up and discordant sign pairs.
- `country_spearman_by_country_2014_2020.pdf`: side-by-side country-specific spider/moth Spearman correlations.
- `loss_with_growth_by_country_2014_2020.pdf` and `gain_with_decline_by_country_2014_2020.pdf`: conditional proportions with country sample sizes recoverable in tables.
- `country_annual_correlations_spider.pdf` and moth: country × annual-interval signed correlation heatmaps (2014→15 … 2019→20).
- `country_directional_summary_2014_2020.csv`: cumulative country×taxon Pearson/Spearman, 2×2 quadrant counts, rates, Fisher test, threshold sensitivity, conditional magnitude correlations.
- `country_directional_summary_annual_and_cumulative.csv`: same for every annual interval, raw signs and thresholded.
- `extreme_population_change_review.csv`: regions with annual |Δpop%| >20% or six-year |Δpop%|>40%, for source QA rather than invisible trimming.

Inspect scatter and heatmaps BEFORE interpreting national differences or doing regression. A country with high population growth and declining PRC supports *compatibility with demographic effects*, but does not prove causation: economic expansion, new industrial lighting, source spectra, settlement change and demographic registration changes can create confounding. To make a statement that PRC **loss is more growth-associated than PRC gain is decline-associated** requires comparing country-level conditional rates, odds ratios and magnitude trends, not cherry-picking Türkiye vs France.

## File-review limitations

Previous 2014→2020 aggregate: N=1,317, percentage population change vs PRC absolute change, spider Pearson r≈−0.080 and moth r≈−0.081; Spearman ≈−0.187 and −0.196. Per-year European correlations vary in sign. These pooled results do NOT settle the user-specified country-by-direction question.
