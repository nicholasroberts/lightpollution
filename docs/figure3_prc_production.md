# Figure 3A–D: locked spider/moth PRC graphics, 2014–2024

**Branch:** `publication-figures`. **Revision: 5 November 2026 — Figures 3A and 3B linear radiance axes.** **Status:** code committed; new render and Illustrator visual checks on EcoVis pending.

## Inputs: existing final calibrated sources

Plot-only script: `src/plot_publication_figure3_prc.py`.

Two **distinct, locked empirical functions**, both R `SSasymp`-form `f(R)=yf+(y0-yf) exp(-kR)`, with **calibrated** VIIRS DNB radiance `R`:

- `config/prc_transfer.yaml`: spider, 17-site fit, locked.
- `config/prc_transfer_moth.yaml`: moth, 16-site fit, locked.

No refitting, transfer-function adjustments, or application of the 0.235 nW sensor offset to the x-axis data. These curves evaluate the existing fit on an already calibrated radiance axis; the 2017+ zero-point correction is applied to native VIIRS pixels in the historical dataset production, *not a second time when evaluating the illustrative curves*.

**Final native-pixel-derived historical PRC tables:**

```text
data/processed/historical_prc_2014_2024_spider_snow_corrected/prc_nuts3_annual_cumulative.csv
data/processed/historical_prc_2014_2024_moth_snow_corrected/prc_nuts3_annual_cumulative.csv
```

Both were generated from the already locked 2021 native-pixel snow correction by `src/publish_prc_locked_snow.py`. The inputs contain 11 years × **1,522 NUTS2016 NUTS3 regions**, and PRC is transformed **at VIIRS pixel resolution before area-weighted NUTS aggregation**. The figure script requires complete cohorts and `snow_status=locked` and verifies each taxon's reported 2024 cumulative change equals the difference between its stored 2014 and 2024 endpoints. It also verifies both taxa use the same valid-pixel area support (within extraction tolerance).

Optional cross-check:

`figures/prc_radiance_spider_moth_2014_2024/nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv`

When present, the script checks the finished locked panel's 2014 and 2024 PRC values against the earlier published joined endpoint table. The earlier joined source was made before the 2021-only patch was adopted; **it must nonetheless have identical 2014/2024 endpoint PRCs**. Failure stops the figure rather than mixing sources.

## Four independent figures

### 3A — locked empirical PRC response functions

Plot `PRC_spider(R)` and `PRC_moth(R)` in one plot, distinguished by fixed colours from `src/figure_style.py`. Default range: 0–100 nW cm⁻² sr⁻¹ **calibrated radiance**. The horizontal axis is now **fully linear**, with displayed ticks at **0, 20, 40, 60, 80 and 100 nW cm⁻² sr⁻¹** for the default range. This replaces the previous symlog x-axis (whose 0–1 interval was linear); it does **not** change, refit or approximate either locked empirical response function.

### 3B — identical absolute radiance change, different starting levels

For every starting calibrated `R>=0`, compute from the same locked curves:

```text
ΔPRC_taxon = f_taxon(R+ΔR) − f_taxon(R)
```

Default illustrative **ΔR=+1 nW cm⁻² sr⁻¹**; alternative positive step via `--delta-radiance`. Figure 3B now uses **the identical fully linear calibrated-radiance x-axis as Figure 3A** (0–100 and 20-nW tick spacing by default), replacing its former symlog display. The two curves show **signed** ΔPRC (negative for increasing radiance), not percentages or losses inferred from particular NUTS3 averages. No region-specific radiance change has been imposed. The resulting numerical function values are exported to a CSV, and **the change of axis does not modify the calculation**.

### 3C — spider PRC change map, 2014→2024

Fixed 1,522 NUTS2016 cohort, original annual **pixel-derived** spider mean PRC: `PRC2024−PRC2014`. Clip to same continental *display window* as Figure 1, lon/lat `(-12,34,36,73)`, projected to EPSG:3035 for drawing. The full 1,522 source records are retained in the output CSV; the mainland display crop excludes outlying territories and may exclude the eastern part of Türkiye, as in Figure 1.

### 3D — moth PRC change map, 2014→2024

Exact same geography, extent, conventions and period as 3C, but **separate, independent moth PRC**. It has its own zero-centred diverging colour bar and **own ±98th-percentile limits** to preserve the different taxon response scales. Blue means negative PRC change; red positive.

Both maps can use **150 m display-only projected polygon simplification**, matching Figure 1's editing-size control. Neither simplification nor the mapping extent affects original zonal statistics or the exported regional table. No missing regions are invented or zero-filled.

**Why 2021 is not reprocessed here:** These are 2014→2024 endpoint maps. They read the final **already snow-corrected** eleven-year series, but the approved 2021-only patch does not change either endpoint. No reason to run the original global 2021 raster again.

## Running on EcoVis

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/plot_publication_figure3_prc.py --self-test
python src/plot_publication_figure3_prc.py
```

Optional display experiments only, keeping source functions/data unchanged:

```bash
python src/plot_publication_figure3_prc.py --delta-radiance 2 --max-radiance 100
```

**Linear-axis regression test:** `--self-test` now verifies that both publication plots' shared `radiance_axis` helper has `get_xscale()=="linear"`, with 0/20/40/60/80/100 ticks at the default 100 nW extent, and also checks the ticks under an alternate 50 nW maximum. The plot still samples the unchanged 0–`--max-radiance` radiance interval at **1,601 points**. `fig3_provenance.json` records `curve_3a_3b_x_axis="linear (0 to maximum calibrated radiance), plot only"` and the exact tick positions. `fig3AB_locked_function_values.csv` retains **identical values** for the same fitting parameters and arguments; the axis transformation affects presentation only. Run the usual full generator to overwrite 3A/3B SVG and PDFs. It also redraws 3C/3D from the same locked inputs without changing their map design or science.

## File outputs — SVG plus vector PDF

Directory `figures/publication/figure3_prc/` (ignored from Git):

```text
fig3A_locked_spider_moth_prc_functions.svg
fig3A_locked_spider_moth_prc_functions.pdf
fig3B_equal_absolute_radiance_change.svg
fig3B_equal_absolute_radiance_change.pdf
fig3C_spider_absolute_prc_change_2014_2024.svg
fig3C_spider_absolute_prc_change_2014_2024.pdf
fig3D_moth_absolute_prc_change_2014_2024.svg
fig3D_moth_absolute_prc_change_2014_2024.pdf
fig3AB_locked_function_values.csv
fig3CD_locked_prc_endpoint_values.csv
fig3_provenance.json
```

## QA after EcoVis run

Verify that **3A and 3B have uniformly spaced linear x ticks**, that both curves are legible (near-zero radiance features appear more compressed than in the former symlog version), and that the locked curve outputs and units remain unchanged. Verify all eight vector files export, legible colour bars, map file size targets and editable text in Illustrator. Check the JSON for the complete 1,522-region cohort, locked source identifiers, the fixed ΔR parameter, and the source-vs-joined-endpoints comparison. **The maps are regionally averaged predictions from pixel-transformed PRC, not PRC of regional mean radiance.**

No fit parameters, scientific datasets, maps from Figures 1/2, or original analysis scripts are altered. If the final locked annual PRC CSV is missing, inspect the already completed `publish_prc_locked_snow.py` workflow before attempting any rerun; do not replace it with an unapproved earlier PRC table.
