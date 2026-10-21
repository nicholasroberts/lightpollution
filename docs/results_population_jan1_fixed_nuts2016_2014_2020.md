# Fixed-NUTS2016 January population-density reconstruction: results

**Date:** 21 October 2026. **Status:** preliminary scientific interpretation of successful EcoVis reconstruction, source/geometry QA not yet independently finalised. The demographic subanalysis is distinct from the **2014–2024** historical PRC study; the spider and moth response functions remain locked.

## Reconstruction and matched cohort

Executed \`python src/reconstruct_population_density_jan1_nuts2016.py\` using Eurostat \`demo_r_pjanaggr3\` January-1 NUTS3 headcounts and \`reg_area3\` land area fixed by region (prefer 2019, fallback to 2018/2020). Project's original full seven-year complete cohort was 1,317 NUTS2016 regions; **1,273/1,317 (96.66%)** have 2014–2020 Jan1 populations and matched land areas on a consistent region code. The **44 exclusions** are **HR 14, AL 12, BE 8, IT 8, EE 2**. The 2017–2020 Jan1 count coverage was 1,285/1,317; area coverage 1,305/1,317. All areas in the included cohort passed the script's 2018–2020 maximum 1% variation flag (zero flags). Same 1,273-region sample and same PRC values used to contrast old \`demo_r_d3dens\` and new \`Jan1 population / fixed land area\`, separately for 2014→2019 and 2014→2020. January density is **not the same variable as Eurostat's annual-average population density**.

Provenance: \`data/raw/eurostat/demo_r_pjanaggr3_nuts3.json\` (Eurostat last updated 28 Aug 2026), \`data/raw/eurostat/reg_area3_nuts3.json\` (22 Jan 2026). \`data/processed/nuts2016_jan1_population_density_sensitivity/reconstruction_manifest.json\` records hashes, exclusion criteria and locked PRC provenance. Production \`nuts3_population_density.parquet\` remains **unchanged**.

## Comparison: primary population predictor absolute Δdensity in people/km² vs absolute ΔPRC

Below are within-country **Spearman ρ**, not spatial-adjusted inferential tests. All old/new comparisons are for precisely the **same region cohort** in each country.

| Country (matched n) | 2014–19 official spider | 2014–19 Jan1 spider | 2014–20 official spider | 2014–20 Jan1 spider | 2014–19 Jan1 moth | 2014–20 Jan1 moth |
|---|---:|---:|---:|---:|---:|---:|
| ALL (1273) | −.158 | **−.156** | −.165 | **−.130** | −.151 | −.135 |
| Türkiye TR (81) | −.492 | **−.485** | −.566 | **−.537** | −.499 | −.555 |
| Germany DE (401) | −.394 | **−.422** | −.365 | **−.381** | −.428 | −.392 |
| France FR (101) | +.163 | **+.194** | +.226 | **+.223** | +.236 | +.261 |
| Ireland IE (8) | +.048 | **+.214** | +.286 | **+.262** | +.214 | +.262 |
| Italy IT (102) | −.194 | **−.262** | −.173 | **−.163** | −.262 | −.159 |
| Greece EL (52) | −.010 | **+.203** | +.038 | **+.436** | +.261 | +.458 |
| Spain ES (59) | +.115 | **+.308** | −.146 | **−.023** | +.318 | +.019 |
| Portugal PT (25) | +.299 | **+.012** | +.319 | **+.117** | −.054 | +.155 |
| Romania RO (42) | −.151 | **−.225** | +.065 | **+.007** | −.219 | −.005 |

**Stable findings:** Türkiye and Germany retain negative country-level rank relationships using an independently defined population predictor, including for both taxa and both endpoints. Their mechanistic hypotheses remain **distinct**. France's rank correlation stays weakly/moderately positive, and its PRC gains are largely independent of demographic decline. Ireland's signs of PRC gain and growing density remain relevant, but with only 8 regions its continuous rank correlation should not be emphasised; the high Pearson coefficient is heavily influenced by IE061.

**Changed/uncertain findings:** Italy's 2014–19 Jan1 rho (−.262) is more negative than 2014–19 old-data (−.194), although 2014–20 remains weaker (−.163); avoid declaring Italy simply null across endpoints. Greece and Spain become more positive for 2014–19, but 2014–20 Spain returns to almost zero, so a long-term monotonic mechanism is not established. Portugal's 2014–19 rho reduces from +.299 to +.012, and Belgium and Romania also vary; therefore country labels must not be based on the original published density series alone.

**Numerical threshold context:** In 2014–19, spider PRC loss >.001 occurs in 79/81 TR regions (31 jointly exceed density-growth +5 people/km² under the Jan1 measure). Germany's 2014–19 Jan1 >5 density-growth group contains 200 of 401 regions (100 of these lose spider PRC >.001); 125 total have PRC losses and 104 gains. France has 94 spider gains >.001 out of 101 regions, with 15 jointly satisfying growth +5 and PRC gain +.001. Ireland has 4/8 spider gains >.001 and two regions exceeding density growth +5; these signed-threshold data should not replace the underlying rank magnitudes.

**No pooled binary-sign discovery:** On the 1,273-region common cohort, spider 2014–19 thresholded Fisher odds ratio is **1.237, p=.616**; spider 2014–20 is **0.649, p=.215**. Many regions remain neutral or mixed (943/1273 and 904/1273); evidence is about **heterogeneous country-level continuous patterns**, not a universally significant pooled up/down sign law.

## Scientific quality assurance and next step

1. The reconstruction solves the bulk of the *2020 repeated published-density* issue without altering old records. The Jan1 counts are independent published information, but fixed land area and continuity of NUTS2016 codes do not alone certify unaltered actual boundaries in every year. Confirm the NUTS2016 support of each source, Jan1 time reference, land-area definition and any population census/revision breaks before calling this the primary series in the manuscript.
2. Prefer the uniform Jan1 count / fixed land area **for all years**, not an ad hoc replacement of 2020 alone. Retain both 2014→19 and 2014→20 as **prespecified paired endpoint comparisons**, recognising that Jan1 t is not annual-average t and could lag a calendar-year radiance exposure. The 2014→24 PRC study is unaffected.
3. Avoid stronger Italy, Greece, Spain, Portugal or Ireland causal interpretations without checking influence, lighting source transitions and the temporal pattern.
4. Next useful methodological decision: verify independent fixed-NUTS2016 geography and source statuses, then freeze the Jan1 demographic source and prepare publication-quality country contrasts without more repeated tuning of PRC calibration or new regressions merely to chase correlations.

Source files (EcoVis, not added to Git): \`country_prc_correlations_published_vs_jan1.csv\`, \`coverage_and_exclusions_by_country_year.csv\`, \`reconstruction_manifest.json\`.
