# Observed population change vs spider and moth PRC: map/scatter results, 2014–2020

**Review:** 21 October 2026, from three investigator-uploaded outputs of \`src/visualize_annual_population_prc_changes.py\`. Analysis is **observation/map/correlation first**, not a causal or nested-model explanatory power claim.

### Files inspected

- \`figures/prc_population_change_2014_2020/cumulative_population_spider_moth_2014_2020.pdf\` — one three-panel map of percentage population density change and absolute spider/moth PRC changes, 2014→2020, fixed NUTS2016
- \`figures/prc_population_change_2014_2020/annual_population_vs_spider_moth_prc_scatter_2015_2020.pdf\` — six annual comparison pages, spider and moth PRC versus measured population density change; coloured by **2014 mean VIIRS radiance band**
- \`figures/prc_population_change_2014_2020/correlations_annual_and_cumulative.csv\` — 56 rows (7 intervals × 2 population-change metrics × 2 taxa × 2 PRC-change metrics), Pearson/Spearman with p-values

**Matching cohort:** 1,317 complete NUTS2016 regions in each annual/cumulative comparison, 205 fewer than the full 1,522 original PRC geography; grey polygons in map signal missing data/eligibility rather than zero change. Both species use locked snow-corrected PRC source series, though **2021 snow is outside this 2014–2020 interval**.

## Main cumulative result: weak geographical correspondence

Measured **2014→2020 percentage population-density change** vs **absolute 2014→2020 ΔPRC** across n=1,317 regions:

| Statistic | Spider PRC | Moth PRC |
|---|---:|---:|
| Pearson r | −0.079577 | −0.080510 |
| Pearson r² | 0.006332 (0.633%) | 0.006482 (0.648%) |
| Two-sided Pearson p | 0.003856 | 0.003459 |
| Spearman ρ | −0.186695 | −0.195805 |
| Spearman p | 8.55e−12 | 7.58e−13 |

Interpretation: overall direction weakly agrees with increased population → more light → reduced receptor contrast, but **population change alone has very little linear association with regional PRC change**. Small p-values are expected for small correlations with n=1,317 and do not demonstrate explanatory importance or causality.

### Annual pattern: NOT stable

Each row refers to population annual **percentage density change** vs absolute ΔPRC of same annual transition, n=1,317 every year.

| Annual transition | Spider Pearson r | Spider Spearman ρ | Moth Pearson r | Moth Spearman ρ |
|---|---:|---:|---:|---:|
| 2014→2015 | −0.109750 | −0.222064 | −0.107352 | −0.229183 |
| 2015→2016 | **+0.068529** | **+0.069397** | **+0.067428** | **+0.066105** |
| 2016→2017 | −0.086275 | −0.109843 | −0.072959 | −0.090637 |
| 2017→2018 | −0.072222 | −0.158537 | −0.076688 | −0.181178 |
| 2018→2019 | −0.003355 | +0.037531 | +0.008736 | +0.049232 |
| 2019→2020 | +0.009086 | +0.048119 | +0.008855 | +0.035276 |

Annual correlations are **weak and change sign**. Spearman is more negative in 2014→15 and 2017→18, but not consistently so. Annual Pearson r² stays ≤~1.2%. There is no stable large annual population–PRC trend common to all years in these raw comparisons.

## Visual evidence and limitations

The three-map spread shows visibly dissimilar *population* and *PRC* geographical colour fields: strongest PRC improvements in some western European regions, decreases concentrated in parts of southeastern Europe, and substantial mixed regional population trends. Spider and moth PRC maps are much more similar to one another than to population-change geography, as expected from applying two separate nonlinear PRC functions to the same VIIRS history. These geographical impressions should not be overinterpreted from one map alone: maps have separately chosen ±98th-percentile PRC scales and population scales; region/mask weights and missing Britain/NUTS3 coverage matter.

**Important data QA concern:** Annual scatter page 2 (2015→16) shows an apparent **~+80% annual population-density change** outlier, page 4 (2017→18) contains ~−30%, and page 6 (2019→20) contains **>+300% annual population-density change**. These values can arise from annual administrative revisions, release-join mismatch, boundary/geographic error, tiny-denominator change or real exceptional events; the plotting PDFs alone do NOT diagnose which. They stretch x axes and obscure the main 99% of observations. Inspect individual NUTS_ID and Eurostat source/area and verify all year-specific NUTS2016 codes before accepting raw Pearson r or asserting demographic shifts as real. Recompute correlations after *transparent* QA, showing both full sample and defensible outlier-excluded sensitivity. Spearman is less sensitive to extreme magnitude but can still be affected by exceptional ranks.

Potential 2016→17 source/baseline radiance shift is treated by the existing 2017+ calibration; do not assume it can never have residual VIIRS measurement effects. Most of this six-year study ends before the approved snow correction, which is important for future 2014→24 work but irrelevant to this particular endpoint correlation.

## Next *visual* analyses, before more regression

1. **Pair each annual map with a four-colour spatial agreement map:** population growing & PRC decreasing (expected sign); population decreasing & PRC increasing (expected sign); growing & PRC increasing; declining & PRC declining. Use a small symmetric threshold to suppress numerical noise, report the threshold and sensitivity. Both taxa have separate agreement maps; differences are expected modestly in sign.
2. **Target the biologically important baseline radiance bands** 0.2–1 and 1–5 nW within the cumulative endpoint and annual comparisons (already generated cumulative focused maps). Report separately their raw correlations and spatial overlap; do not mistake a pooled correlation driven by regional mix for within-band association.
3. **Read/QA the existing \`annual_changes_2015_2020.csv\` and \`cumulative_change_2014_2020.csv\`** to identify extreme 2015–16 and 2019–20 NUTS3 population values, verify stable boundaries and region area, then plot main-density crop or hexbin with points outside the view transparently counted. The PNG scatter has no readable central cloud when x range stretches to +350%.
4. Only then consider adjusted models or a causally framed discussion. **Do not** report "population change drove X% of predicted biological change" from the current correlation.

**Conclusion to date:** raw geographical and annual correlations offer little evidence that measured regional demographic changes consistently track the changing spatial distribution of predicted polarization-receptor contrast over 2014–2020. The proposed population-density → absolute radiance → PRC nonlinear story describes **baseline sensitivity**, but temporal population-growth association must be tested separately. This null/weak empirical result is scientifically useful and should not be forced into a narrative of strong demographic causation.
