# Figure 3A–F: locked spider/moth PRC graphics, 2014–2024

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

## Six independent figures

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


### 3E — distribution of signed absolute spider PRC change, 2014→2024

Same **1,522 fixed NUTS2016 NUTS3 region** cohort and *exact* locked native-pixel-derived PRC endpoint differences as Figure 3C: \`spider_prc_2024 - spider_prc_2014\`. **Absolute change** means a difference in PRC *units*, not percentage change and **not** \`abs(ΔPRC)\`: positive differences (cue gains) and negative differences (losses) remain on opposite sides of zero.

An equal-region **count** histogram, following Figure **1D's** figure size and style: 55 signed-log-distributed bins, blue bars, vertical dashed orange equal-NUTS3 mean, solid black median and dotted grey zero. The histogram uses a **symmetric-log x-axis** with a default linear interval of ±0.001 PRC (customisable via \`--prc-hist-linthresh\`; display-only). Bins use the same signed \`log1p\` transformation as Figure 1D; labels remain in true dimensionless ΔPRC units. No data clipping, winsorisation, percentage denominator, selection on radiance, area weighting or resampling. A note reports the full region count. The **mean/median are arithmetic statistics of the untransformed 1,522 signed values**, not positions computed in transformed space.

### 3F — distribution of signed absolute moth PRC change, 2014→2024

Identical distribution design and cohort, but the independently locked moth PRC endpoints used for Figure 3D. The same signed-log threshold facilitates comparison, while the bin edges cover each taxon's **entire observed range** independently. Separate taxon summary statistics and 55 histogram-bin counts are recorded; these are not measurements from two independent satellite datasets.

### Figure 1D comparison and scientific units

Figure 1D is a histogram of **percentage radiance change**, which requires strictly positive 2014 radiance and has a very long positive tail. Figures 3E/3F instead show **signed absolute PRC differences** (dimensionless). All 1,522 locked PRC endpoints are finite, so the PRC histograms include **all 1,522** and need **no zero-baseline exclusions**. The symmetric-log scale is solely a presentation device for potentially skewed signed changes; it must not be mistaken for log-transformed statistics or percentage PRC. Separate bars represent counts in equal transformed-width bins. The histograms and maps are NUTS3 regional estimates from native-pixel PRC transformation; a region's ΔPRC is not obtained by evaluating the biological response at the region's average radiance.

## Running on EcoVis

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/plot_publication_figure3_prc.py --self-test
python src/plot_publication_figure3_prc.py
# Optionally change the histogram display-only signed-log threshold:
# python src/plot_publication_figure3_prc.py --prc-hist-linthresh 0.002
```

Optional display experiments only, keeping source functions/data unchanged:

```bash
python src/plot_publication_figure3_prc.py --delta-radiance 2 --max-radiance 100
```

**Linear-axis regression test:** `--self-test` now verifies that both publication plots' shared `radiance_axis` helper has `get_xscale()=="linear"`, with 0/20/40/60/80/100 ticks at the default 100 nW extent, and also checks the ticks under an alternate 50 nW maximum. The plot still samples the unchanged 0–`--max-radiance` radiance interval at **1,601 points**. `fig3_provenance.json` records `curve_3a_3b_x_axis="linear (0 to maximum calibrated radiance), plot only"` and the exact tick positions. `fig3AB_locked_function_values.csv` retains **identical values** for the same fitting parameters and arguments; the axis transformation affects presentation only. Run the usual full generator to overwrite 3A/3B SVG and PDFs. It also redraws 3C/3D from the same locked inputs without changing their map design or science, then adds the new signed-ΔPRC 3E/3F histograms.

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
fig3E_spider_absolute_prc_change_distribution_2014_2024.svg
fig3E_spider_absolute_prc_change_distribution_2014_2024.pdf
fig3F_moth_absolute_prc_change_distribution_2014_2024.svg
fig3F_moth_absolute_prc_change_distribution_2014_2024.pdf
fig3EF_absolute_prc_change_summary.csv
fig3EF_absolute_prc_histogram_bins.csv
fig3AB_locked_function_values.csv
fig3CD_locked_prc_endpoint_values.csv
fig3_provenance.json
```

## QA after EcoVis run

Verify that **3A and 3B have uniformly spaced linear x ticks**, that both curves are legible (near-zero radiance features appear more compressed than in the former symlog version), and that the locked curve outputs and units remain unchanged. Verify all twelve vector files export, legible colour bars, map file size targets and editable text in Illustrator. Check the JSON for the complete 1,522-region cohort, locked source identifiers, the fixed ΔR parameter, and the source-vs-joined-endpoints comparison. For **3E/3F**, verify both histograms have **1,522 counted regions**, that positive/negative PRC changes are represented on their correct sides of zero, that the equal-region mean/median markers agree with `fig3EF_absolute_prc_change_summary.csv`, and that the sum of each taxon's 55 `fig3EF_absolute_prc_histogram_bins.csv` counts is **1,522**. The `--self-test` now checks exact signed differences, preserved tail values, mean/median and bin counts on a toy two-taxon dataset. Figure 3C/3D map colours are unrelated to histogram bar colours; histograms are not area weighted. **The maps are regionally averaged predictions from pixel-transformed PRC, not PRC of regional mean radiance.**

No fit parameters, scientific datasets, maps from Figures 1/2, or original analysis scripts are altered. If the final locked annual PRC CSV is missing, inspect the already completed `publish_prc_locked_snow.py` workflow before attempting any rerun; do not replace it with an unapproved earlier PRC table.


## Companion illustrative Figure 3A: +0.05 PRC from 80 versus 5 nW

**Added 5 November 2026 (presentation example; distinct from the unchanged original Figure 3A/3B).** Script:
`src/plot_publication_fig3a_equal_prc_gain.py`.

This illustration answers the specific inverse question: **starting from a calibrated VIIRS radiance of either 80 or 5 nW cm⁻² sr⁻¹, how far must radiance *decrease* for the empirical spider or moth receptor contrast to *increase by 0.05 absolute units*?** The +0.05 is **not a 5% relative change** or a measured behavioural threshold. The original independently locked spider17 and moth16 response functions are `P(R)=yf+(y0-yf) exp(-kR)` and are loaded/validated from the same original YAMLs used for Figure 3A, without refitting. Solving exactly gives:

```
P_target = P(R_start) + 0.05
R_final = -ln[(P_target - yf)/(y0 - yf)] / k
radiance_reduction = R_start - R_final
```

| Taxon | Starting R | Required final R | Absolute R reduction | Reduction from starting R |
|---|---:|---:|---:|---:|
| Spider | 80 | 17.60284 | 62.39716 | 77.996% |
| Spider | 5 | 2.69388 | 2.30612 | 46.122% |
| Moth | 80 | 31.17422 | 48.82578 | 61.032% |
| Moth | 5 | 3.48070 | 1.51930 | 30.386% |

The absolute radiance decrease needed for spider is approximately **27.06 times** greater at the 80-nW baseline than at the 5-nW baseline; for moth approximately **32.14 times** greater. This reflects the saturation of the empirical PRC response functions at high R. These quantities are an **illustrative inversion of fitted relationships**, not a prediction that removing a given light source will restore the nocturnal sky polarization field by this exact amount. It assumes the original empirical transfer, species and light-spectrum context apply.

**Display:** One 16:9 presentation figure, two species rows: left main **linear 0–100** calibrated VIIRS radiance plot showing 80→final, right **linear 0–7 zoom** highlighting 5→final. Each row uses precisely the original 3A species curve. Orange circles show starting R/PRC and teal diamonds show the lower final R that produces a +0.05 increase; the intervening transfer-curve segment is highlighted with a directional arrow and a vertical +0.05 PRC bracket. *Y axes are specific to each plot (not shared between taxa or zoomed panels)*. The four numerical reductions are printed below their respective panels.

**Files** in `figures/publication/figure3_prc/`, without overwriting any existing Figure 3 files:

```
fig3A_example_equal_0p05_prc_gain_80_vs_5.svg
fig3A_example_equal_0p05_prc_gain_80_vs_5.pdf
fig3A_example_equal_0p05_prc_gain_80_vs_5.png
fig3A_example_equal_0p05_prc_gain_calculations.csv
fig3A_example_equal_0p05_prc_gain_provenance.json
```

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate
python -m py_compile src/plot_publication_fig3a_equal_prc_gain.py
python src/plot_publication_fig3a_equal_prc_gain.py --self-test
python src/plot_publication_fig3a_equal_prc_gain.py
```

**Reproducibility and state:** A preview image, vector SVG and PDF have also been generated and visually inspected in the ChatGPT working container from the same locked numerical coefficients. The **new GitHub script has not yet been executed on EcoVis**. The four inverse values are verified analytically. On EcoVis, the script reads the authoritative config files afresh and records the exact source model IDs and parameters in the new provenance JSON. The original `fig3AB_locked_function_values.csv`, original 3A/3B curves and Figure 3C/3D annual maps remain unchanged.

