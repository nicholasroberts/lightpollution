# Spider PRC threshold scenario: DoP = 0.23 and European area series (2014–2024)

**5 November 2026 | publication-figures branch | exploratory talk analysis.**
Script: `src/plot_publication_spider_prc_dop023_threshold.py`.

## Scientific definition

Foster et al. (2019, *Journal of Experimental Biology*, DOI: [10.1242/jeb.188532](https://doi.org/10.1242/jeb.188532)) report **weaker**, but nonzero, dung beetle orientation at skylight **DoP < 0.23** compared with fully polarized cues. Their published behavioural detection interval is broad and explicitly **not** a validated value for our spider `Drassodes`. Do **not** call 0.23 a hard spider signal-detection threshold, individual navigational failure, or a measured DoP at every VIIRS pixel.

For the indicative **opponent receptor** model with orthogonal polarization axes, optimal e-vector alignment and **polarization sensitivity Sp=10**, the assumed absolute dimensionless PRC (not a delta) is:

```text
PRC = ln(R_parallel / R_perp)
PRC_threshold = 2 * atanh(DoP_threshold * (Sp-1)/(Sp+1))
PRC_threshold(0.23, Sp=10) = 0.3809031327477229
```

The formula is contingent on those receptor and alignment assumptions, which are distinct from the empirical locked `VIIRS → spider PRC` transfer. It does not refit or alter `config/prc_transfer.yaml`. The **0.3809 PRC** scenario is an interpretation of the PRC numbers, **not a calibrated spider behavioural threshold**. Lower PRC may weaken navigational signal without making polarization invisible.

## Locked data and maps

Read the **already published, locked** spider `prc_mean` annual NUTS3 table at:

```text
data/processed/historical_prc_2014_2024_spider_snow_corrected/prc_nuts3_annual_cumulative.csv
```

Exactly **1,522 fixed NUTS2016 regions** × all **eleven annual years 2014–2024**. Original PRC was transformed from calibrated VIIRS pixels before NUTS3 aggregation. **2021 uses the locked native-pixel snow-corrected spider PRC** already present; this plotting script performs no correction, sensor offset, radiance recalculation or PRC fit. Calls the same source validator and independent 2014/2024 endpoint QA as the existing annual spider PRC atlas.

Eleven full 16:9 maps follow the already-approved radiance/absolute PRC atlas geography and framing: same NUTS2016 polygons, EPSG:3035 map, 150-m optional geometry simplification **for drawing only**, 34–73°N / -12–46° longitude continental extent (including eastern Türkiye), big year/title, shared colourbar, no footer. Unlike the original annual PRC map palette, this map uses a single **fixed threshold-centred colour scale** computed from pooled 2014–2024 PRC extrema, extended to include the scenario threshold if necessary. Below/at threshold: deep red → orange → yellow (more strongly compromised PRC toward red); strictly above: yellow-green → green → teal. A black colourbar mark identifies **PRC 0.3809 = assumed DoP 0.23**. The colour transform is `TwoSlopeNorm`, piecewise linear around the threshold, **not a change to PRC**. The exact shared display limits and valid regions are written in provenance; no annual rescaling and no clipping.

**Very important:** Every coloured polygon is assigned its already-computed **annual NUTS3 mean** PRC; the map does not claim that *all native pixels* within a below-threshold region fall below the PRC threshold.

## European area time series: primary estimand and limitations

The primary metric is:

```text
100 × Σ_{i in full 1,522 NUTS2016; regional_mean_PRC_i,y <= 0.380903}
        fixed_2014_valid_area_km2_i
      / Σ_{i in full 1,522 NUTS2016} fixed_2014_valid_area_km2_i
```

Weights are each region's **2014 valid PRC pixel-coverage area**, held constant across all 2014–2024 observations. This measures the share of the locked NUTS2016 **study's valid-pixel-supported area** belonging to whole regions whose mean spider PRC is at/below the scenario threshold. It is NOT the true percentage of European **native VIIRS pixels** with PRC below threshold, nor a direct biological performance estimate. `PRC(mean radiance)` is not used; the existing mean of **pixel-derived PRC** is classified. Thresholding NUTS3 mean PRC is not equal to thresholding pixels before averaging, particularly for mixed cities/countryside regions (ecological aggregation).

Additional CSV columns report (a) **the continental display-cohort** sensitivity: all NUTS3 regions intersecting the fixed map bbox, counted with their whole-region 2014 valid area (not partial clipped polygon area); (b) the proportion of region counts below threshold; and (c) an annual-valid-area-weight sensitivity. The primary full-locked and alternate display-cohort area percentages are plotted on one graph; the y-axis is **explicitly zoomed** for readable temporal changes (the original values are available in CSV). Preserve both denominators and do not silently refer to the continental-subset figure as all-Europe area.

Only a genuinely **pixel-level threshold classification on original VIIRS pixel-derived PRC**, with original pixel support and area weights, could estimate the precise physical area of below-threshold PRC. That is a separate heavier analysis, not undertaken in this plot-only initial test.

## Run on EcoVis

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/plot_publication_spider_prc_dop023_threshold.py --self-test
python src/plot_publication_spider_prc_dop023_threshold.py
```

Outputs in the **new sibling folder**:

```text
figures/publication/annual_nuts3_spider_prc_dop023_threshold_2014_2024/
  spider_prc_dop023_threshold_2014.svg/.pdf/.png
  ...
  spider_prc_dop023_threshold_2024.svg/.pdf/.png
  all_years_spider_prc_dop023_threshold_2014_2024.pdf
  spider_prc_dop023_threshold_area_percent_time_series.svg/.pdf/.png
  spider_prc_dop023_threshold_area_percent_2014_2024.csv
  annual_nuts3_spider_prc_threshold_classification.csv
  spider_prc_dop023_threshold_provenance.json
```

The source CSV, original annual PRC map folder and all original Figure 3 outputs are **not modified**.

The script supports a sensitivity check with another DoP level and polarization sensitivity (`--dop-threshold`, `--polarization-sensitivity`), but for comparisons **use a different `--out` folder** because the default output filename stem is intentionally named for DoP=0.23.

## QA before presenting outcomes

- `--self-test` verifies 0.23/Sp=10 → PRC 0.3809031327; inverse formula; below-or-equal classification at the exact threshold; and 2014-reference area-weight arithmetic using synthetic examples.
- Inspect 2014, 2021 and 2024 maps together. Ensure the black colourbar tick and colour split correspond to 0.381 and all years share exactly the same palette limits. Label the maps as an exploratory reduced-performance scenario for a different nocturnal taxon, not established extinction or navigational blindness.
- The area CSV contains 22 observations (eleven years × two defined geographies), exact fixed-area totals, below-threshold km², yearly sample counts, annual area-weight sensitivity and counts of below-threshold administrative regions. Check the full-locked series has exactly 1,522 NUTS3 and unchanged 2014 denominator every year.
- The 2014–2024 percentages and trend **cannot be known until the real locked EcoVis input table is read**; none are guessed or inferred from the locked exponential fit. Examine whether 2021 deviation is present despite the approved snow correction.
- The observed fraction of Europe with low regional mean PRC is **a model-derived exposure indicator**, not observed harm to species, species occupancy, population decline or navigational failure.
