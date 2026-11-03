# Figure 1A–C and NUTS3 histogram — illustrator-ready exports

**Current revision: 3 November 2026. Branch: `publication-figures`.**

## Objective and rules

Render **four** standalone vector figures with `src/plot_publication_figure1_radiance.py`, from existing established VIIRS annual data, the frozen NUTS2016 geography, and the already **LOCKED** snow/albedo raster patch for 2021. No re-fitting, spectral adjustment, population analysis, pixel-offset recalibration, or alteration of the original satellite files.

The first version had (1) a gap at 2021 and (2) 57 grey NUTS3 percentage-map polygons caused by a **deliberate `2014 radiance >=0.2 nW cm^-2 sr^-1` threshold in the existing `src/plot_prc_radiance_2014_2024.py` output.** Those were an earlier display/QA rule, not loss of the endpoint measurements. Both are addressed now.

## Sources

1. **Final 2014/2024 endpoints**:
   `figures/prc_radiance_spider_moth_2014_2024/nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv`; full **1,522 fixed NUTS2016 NUTS3**. The calibrated `radiance_2014`, `radiance_2024`, and `radiance_change_abs` columns are reused unchanged.
2. **European annual history 2014–2024**:
   `data/processed/fixed_nuts2016_viirs_2014_2024.parquet`, field `radiance_mean_corrected`. This has the 2017+ 0.235 pixel background correction already applied, but 2021 needs the separate native-pixel snow overlay.
3. **Locked 2021 source**:
   `config/viirs_2021_snow_correction.yaml`; its `patch` and `patch_summary` locate `data/processed/prc_2021_spatial_europe_candidate/vnl_2021_snow_albedo_europe_patch.tif` and the matching manifest. The filename retains `candidate` from its original generation; the **config status is now `locked`**.
4. **Boundaries**: the original fixed NUTS2016 level-3 geometry, with optional 150 m simplification for **vector display only** in EPSG:3035, never for extraction.

## Exactly how 2021 is now included

`src/figure1_2021_locked_snow_radiance.py` uses the same original VIIRS 2021 global annual raster and native-grid patch alignment/validation method as the previously locked spider PRC correction. For each window:

- Read the original 2021 annual **raw** VNL pixel values.
- Substitute band 1 of the already locked Europe-window patch only where valid and positive original pixels overlap. Preserve NoData and zero background, reject any brightening.
- **After** overlay, subtract the existing baseline calibration offset **0.235** from positive pixels **once only**, with lower bound 0. Do not subtract it from the regional mean or from the patch a second time.
- Use `exactextract` with `mean(coverage_weight=area_spherical_m2)` on the original 1,522 NUTS2016 regions. Independently extract unpatched 2021 calibrated radiance in the same pass and verify that the patch never increases any NUTS mean or changes valid coverage.
- Compare original/unpatched 2021 extraction with the corresponding rows in the existing annual history, **before replacing** just that year's radiance. If they disagree outside the specified tolerance, the script stops rather than mixing incompatible sources.
- Apply the corrected 2021 regional means to the existing balanced annual cohort, then calculate the European area-weighted 2014–2024 history using **fixed 2014 area weights**. All other years stay unchanged.
- Cache the 1,522-region overlay extraction under `data/processed/publication_figures/fig1_2021_locked_snow_radiance.csv` with JSON provenance. Cache is reused only if the locked model, original raster and patch paths/size/mtime and baseline offset match; `--rebuild-snow-2021` forces re-extraction.

**The source parquet's complete annual cohort may be smaller than the full 1,522 endpoint cohort** (historically 1,514 in a model-specific matched stage). The Figure 1C series reports the actual cohort; the maps and distribution use all 1,522 endpoints. No additional selection is made for the population dataset.

## All-region percentage-change map and histogram

The older published joined endpoint CSV has `radiance_change_pct` missing for **57 regions** because of the optional 0.2-nW baseline threshold. Figure 1 recalculates from the existing endpoints into the new *display-only* field:

`radiance_change_pct_all = 100 * (radiance_2024 - radiance_2014) / radiance_2014`

- **Every region with genuinely positive 2014 radiance is included**, even when its baseline is below 0.2. Reconstructed percentages must agree with the original archived percentages wherever they existed.
- **If and only if the original 2014 mean is exactly zero**, percentage change cannot be calculated (0 denominator). Such regions remain **visible in purple** in the percentage map, explicitly marked *percentage undefined*. Do **not** give them invented `0%`, `+infinity%` or an arbitrary denominator.
- The histogram includes every mathematically defined regional percentage, even very large increases from extremely dim starting values. It uses a signed-log x axis and shows the **equal-NUTS3 arithmetic mean and median**. The sample count and exactly-zero exceptions are reported dynamically in `fig1_provenance_and_summary.json`.
- **Mean of regional percentage changes is not the same as percentage change of the whole-Europe area-weighted mean radiance**; the latter is included separately in the summary. Expect a highly skewed histogram when very dim baselines are included.

## Production commands on EcoVis

From the root of the `lightpollution` clone:

```bash
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/figure1_2021_locked_snow_radiance.py --self-test
python src/plot_publication_figure1_radiance.py --self-test
python src/plot_publication_figure1_radiance.py
```

The first full run creates the locked 2021 radiance overlay extraction cache; subsequent unchanged runs reuse it. If there is an existing `fig1_2021_locked_snow_radiance.csv` but the original raster or patch changed, it will be re-extracted.

## Outputs

`figures/publication/figure1_radiance/` (ignored by Git):

- `fig1A_absolute_radiance_change_2014_2024.svg` + PDF
- `fig1B_percent_radiance_change_2014_2024.svg` + PDF
- `fig1C_europe_annual_radiance_2014_2024.svg` + PDF, **with corrected 2021 point**
- `fig1D_nuts3_percentage_change_distribution_2014_2024.svg` + PDF
- `fig1C_europe_annual_radiance_support.csv` (includes original 2021 European mean for audit)
- `fig1_provenance_and_summary.json` (includes full 2021 source audit and percentage counts).

The separately ignored `data/processed/publication_figures/` cache contains the 1,522-region **original and snow-corrected 2021** regional radiance table plus its provenance JSON.

## QA notes

Check the new histogram mean and median with the expanded sample. With nearly zero 2014 means, the unweighted mean may be extremely sensitive; don't silently truncate data, winsorise or replace denominators to make the result look reasonable. Confirm every NUTS2016 region is displayed either by its finite percentage or the explicit zero-baseline category; grey should not mean an intentionally dropped nonzero region.

**Status:** source code and documentation committed. These scripts must still be run on EcoVis, where the original VNL TIFFs and locked native 2021 patch are stored; no successful actual render is claimed before that run.
