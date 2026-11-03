# Figure 1A–C and added NUTS3 histogram: first vector export

## Script

`src/plot_publication_figure1_radiance.py` on branch `publication-figures`.

**All outputs use the pre-existing finished data tables and mapping definitions; this script never reruns radiance extractions, calibration, population models or PRC fits.**

From the repository root on EcoVis, after `git switch publication-figures && git pull --ff-only`:

```bash
source .venv/bin/activate
python src/plot_publication_figure1_radiance.py --self-test
python src/plot_publication_figure1_radiance.py
```

If the environment is already activated, omit the `source` line.

## Exact inputs

1. **Figures 1A, 1B and histogram:**
   `figures/prc_radiance_spider_moth_2014_2024/nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv`.

   This is the archived 1,522-region NUTS2016 source of final corrected-band VIIRS DNB 2014/2024 radiance estimates (2017 onward pixel baseline correction already applied), produced by the existing radiance/spider/moth comparison. The endpoints are unaffected by the separate 2021 snow patch. The script validates 1,522 rows, duplicate-free NUTS IDs, numeric 2014/2024 radiance and equality of the archived absolute and percentage changes to the endpoints. **Do not replace this input with the old population-restricted 1,273-region cohort or the older 1,514-region diagnostic cohort.**

2. **Geographic display:**
   `historical_prc_pixel_nuts3.load_nuts(load_config())` using fixed NUTS2016 level-3 polygons. The continental map window is [-12, 34, 36, 73] in lon/lat (outlying study territories retained in calculation tables but not displayed in mainland map window). The graphics transform to EPSG:3035. Display-only 150 m polygon simplification reduces SVG size; no changes to NUTS boundaries used in any scientific computation.

3. **Figure 1C annual history:**
   `data/processed/fixed_nuts2016_viirs_2014_2024.parquet` using its **already calibrated** `radiance_mean_corrected` column. The script joins the fixed 2014 geographical area weights from the final endpoint table and uses the *balanced intersection of NUTS3 regions available every year*; the resulting cohort n is given in the plot and provenance summary. Before rendering, it checks 2014 and 2024 NUTS values against the final endpoint source and stops for investigation if discrepancies are excessive. The available annual source may have a different complete-coverage cohort from the final 1,522 endpoint regions. This difference is not hidden or reweighted to make it look like the same table.

   **2021 restriction:** the locked 2021 snow/albedo native-pixel patch is applied to the finished spider/moth PRC histories, but the older annual NUTS3 `radiance_mean_corrected` column has not been demonstrated to contain the same corrected 2021 radiance field. Therefore **Figure 1C omits 2021 (a visible gap)** rather than presenting the known uncorrected 2021 anomaly as a final radiance estimate. The QA CSV retains the original 2021 value labelled `figure1c_approved=False`, not for plotting. Producing a genuinely corrected full 2014–2024 radiance curve including 2021 would first require the **already locked 2021 native-pixel patch to be propagated through the radiance zonal statistics**; no alternate model or approximate correction is to be invented for this figure.

## Output

`figures/publication/figure1_radiance/` (ignored by Git):

- `fig1A_absolute_radiance_change_2014_2024.svg`, `.pdf`
- `fig1B_percent_radiance_change_2014_2024.svg`, `.pdf`
- `fig1C_europe_annual_radiance_2014_2024.svg`, `.pdf`
- `fig1D_nuts3_percentage_change_distribution_2014_2024.svg`, `.pdf`
- `fig1C_europe_annual_radiance_support.csv`
- `fig1_provenance_and_summary.json`

Figure 1A and B use signed, zero-centred symmetric-log maps with blue = decreasing and red = increasing radiance, 98th-percentile colour saturation and correctly labelled physical units. The maps are vector (including country borders).

**Histogram:** The existing final 2014 radiance cutoff is **≥ 0.2 nW cm⁻² sr⁻¹**, so 57 of the 1,522 regions have undefined/suppressed endpoint percentages in the archived source and are also excluded from this histogram. **1,465** NUTS3 percentage changes are expected in the histogram. Vertical lines mark the **unweighted arithmetic mean and median across these eligible individual NUTS3 regions**; these statistics must not be confused with the separately reported **European area-weighted mean radiance percentage change** (in the report JSON), or the average over human population. All eligible data points are counted, with a symmetric-log horizontal axis to show the full positive/negative tails.

## Quality and provenance checks

- Confirm the output summary reports n=1,522 endpoints, 1,465 valid percentages and 57 unavailable low-baseline percentages.
- Verify that the endpoint area-weighted European change is consistent with the previously published 2014/2024 result (+1.37% for the locked endpoint table); if it differs substantially, verify area weighting/source definitions before using the output.
- Confirm the annual table's matched cohort and endpoint cross-check, and **do not describe the 2021 gap as a zero radiance value**.
- Open both vector versions in Illustrator: text should remain editable, and no map boundaries or clipped ticks should be missing. Font compatibility is controlled by `src/figure_style.py`.
- Export settings, all figures and descriptive summaries are reproducible; generated images are deliberately not committed.
- Print output list at end of the run; retain source provenance in JSON.

**Current state:** script and instructions committed to GitHub. The large EcoVis processed data are not available in the assistant's container, so the scientific exports must be rendered with the command above on EcoVis before they can be visually inspected.
