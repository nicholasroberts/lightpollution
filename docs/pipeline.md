# Pipeline design

## Scientific objective

Create a reproducible annual dataset linking VIIRS night-time radiance to
population density for European NUTS3 regions.

The analytical unit is an available `year × NUTS_ID` observation. Within each
country and year, matched NUTS3 radiance/population-density pairs are treated
as the observations used for the log-radiance versus log-population-density
relationship. A NUTS3 region does not need to exist in every year.

## Eurostat-aligned NUTS geography

The workflow uses the NUTS classification represented by the current Eurostat
`DEMO_R_D3DENS` population-density time series. Eurostat retrospectively
recasts historical regional statistics onto newer NUTS classifications, so
the geographical coding in the live population dataset does not always equal
the NUTS release formally applicable in the original reference year.

Direct comparison of the population codes against GISCO NUTS3 releases showed
the following best-matching geography:

- 2013-2020: NUTS 2016
- 2021-2022: NUTS 2021
- 2023-2024: NUTS 2024

VIIRS radiance is therefore aggregated using exactly this mapping. This gives
radiance and population density compatible NUTS identifiers within each year
and avoids large artificial losses of observations caused by joining
retrospectively recoded Eurostat population data to obsolete historical NUTS
codes.

This is a deliberate analytical choice. The objective is the country-level
relationship between population density and radiance, so preserving compatible
regional units is more important here than retaining the exact boundary regime
that was legally in force in each historical year.

The required GISCO Level-3 GeoJSON releases are downloaded and retained
locally. Population observations that still do not match their configured
Eurostat-aligned release are flagged and omitted rather than imputed.

## Pipeline stages

1. **nuts**
   - download NUTS 2016, 2021 and 2024 Level-3 GISCO GeoJSON
   - use 01M resolution and EPSG:4326
   - use `wget` when available

2. **validate**
   - discover local VIIRS annual rasters
   - inspect raster CRS/dimensions
   - confirm all required NUTS releases
   - report the Eurostat-aligned year-to-NUTS-release mapping and missing VIIRS years

3. **weighting-test**
   - compare pixel-centre, exact fractional-overlap and spherical-area-weighted
     means on representative NUTS3 regions

4. **radiance**
   - select the Eurostat-aligned NUTS release configured for each VIIRS year
   - calculate exact-overlap, spherical-area-weighted mean radiance
   - retain valid covered area and provenance
   - write the unmodified raw radiance table

5. **calibration**
   - preserve extracted NUTS3 radiance as `radiance_mean_raw`
   - from 2017 onward apply one 0.235 nW cm^-2 sr^-1 additive correction to
     every positive Annual VNL V2.2 pixel across all analysis countries
   - clip corrected positive pixels at zero; raw zero/background pixels remain zero
   - perform the transformation on-the-fly through an exactextract RasterSource,
     so no duplicate corrected global GeoTIFF is stored
   - calculate the corrected NUTS3 mean with the same exact polygon overlap and
     spherical physical-area weighting used by the raw extraction
   - retain raw-positive area fraction, clipped area fraction, effective
     NUTS3-level offset, method/version and full correction provenance

6. **population**
   - retrieve Eurostat `DEMO_R_D3DENS`
   - request `geoLevel=nuts3`
   - retain the raw JSON-stat response
   - use `wget` by default
   - associate each year with the NUTS release used by the harmonized Eurostat population series
   - flag whether each population record matches that release
   - preserve Eurostat status flags

7. **population-fallback-test**
   - diagnostic only; does not alter production data
   - test whether missing `DEMO_R_D3DENS` values can be recovered from
     `demo_r_pjanaggr3` population on 1 January divided by `reg_area3`
     total land area
   - compare reconstructed density with published density wherever both exist
   - quantify the expected difference caused by annual-average versus
     1-January population definitions

8. **population-recover**
   - preserve every published `DEMO_R_D3DENS` value unchanged
   - restrict recovery to NUTS3 codes belonging to the configured geography
   - fill only absent/null density observations using
     `demo_r_pjanaggr3 / reg_area3` total land area
   - restrict production recovery to a validated country allow-list
     (currently Belgium, Estonia, Croatia, Italy and the Netherlands)
   - retain explicit provenance and recovery flags for every derived value
   - Norway is excluded from production recovery because validation overlap is
     sparse and several target NUTS3 regions remain unrecovered

9. **merge**
   - inner join on `year + NUTS_ID + nuts_release`
   - omit unmatched/missing observations rather than harmonising or imputing
   - calculate log10 radiance and log10 population-density fields
   - report the number of paired observations for every country/year

10. **fit**
   - fit a separate OLS line in log10-log10 space for every country/year
   - report n, slope, intercept and R²
   - exclude only non-positive values from the logarithmic fit
   - write a fit-summary CSV and a multi-panel comparison figure

11. **benchmark**
   - compare rebuilt country fits with the archived 2016-2019 analysis
   - report old mean/range and new slope/intercept values
   - flag whether the new value lies inside the legacy four-year range
   - write a comparison CSV and benchmark figure

12. **trends**
   - plot country-specific slope, intercept and R² from 2013-2024
   - mark the 2016/2017 VIIRS calibration transition
   - calculate year-to-year changes in slope, intercept and R²
   - write the change table for later calibration analysis

## Radiance statistic

The production radiance extractor uses exact polygon overlap with
`coverage_weight=area_spherical_m2`. VIIRS is supplied on a geographic
longitude/latitude grid, so physical cell area varies with latitude.

The 2024 validation showed that exact boundary treatment changed mean radiance
by about 2% for some small/coastal regions, while spherical physical-area
weighting changed the result by about 2.5% for a large high-latitude region.
The exact-overlap, area-weighted mean is therefore the production radiance
metric.

The median is not used in the production analysis. Testing showed that the
weighted median/quantile returned by the extraction library did not correspond
to the ordinary pixel-value median expected for zero-dominated dark regions,
whereas the mean passed the direct raster checks and is the relevant metric for
the population-density relationship.

## Local data

Large raw and generated datasets are deliberately excluded from Git. Expected
local locations are defined in `config/pipeline.yaml`.


13. **calibration-offset-test**
   - diagnostic only; does not alter production radiance
   - use matched NUTS3 regions across 2015-2018
   - estimate the excess 2016-to-2017 change in linear radiance
   - compare that estimate across 2016-radiance deciles and countries
   - test whether a simple additive zero-point shift is sufficient


14. **calibration-country-test**
   - diagnostic only; does not alter production radiance
   - estimate a country-specific additive 2017 offset from the darkest 20% of
     matched 2016 NUTS3 regions
   - subtract that constant from 2017+ radiance for the test only
   - re-fit annual country relationships and compare the 2016-to-2017
     slope/intercept discontinuity before and after correction


15. **calibration-sensitivity**
   - diagnostic only; does not alter production radiance
   - repeat the country-specific additive offset estimate using the darkest
     10%, 15%, 20%, 25% and 30% of matched regions
   - quantify residual 2016-to-2017 slope/intercept discontinuities
   - compare corrected 2017-2019 changes relative to 2016 against the archived
     legacy analysis, which used an upstream zero-point correction
   - rank thresholds jointly by discontinuity removal and legacy agreement


## KNOWN CALIBRATION-SCOPE WARNING

The current branch historically applied the 2017 additive correction only to
DE, IT, NL, FR and ES because those were the original comparison countries used
to validate the diagnostic. This is now recognised as incorrect in scope: the
2017 DNB calibration change is dataset-wide, not country-specific.

Until the global correction diagnostic has been evaluated, the existing
`radiance_mean_corrected` column should be treated as provisional for
cross-2017 temporal analyses. Use:

```bash
python src/run_pipeline.py global-calibration-rotation-test
```

to compare raw radiance, the legacy five-country correction, a published global
0.15 nW cm-2 sr-1 correction, and two Europe-wide empirical global offsets.
The production calibration will be replaced only after this robustness test.

## Production radiance calibration

The production correction was revised on 2026-10-12 after direct tests of the
Annual VNL V2.2 `average_masked` rasters established that the 2016-to-2017
discontinuity remains in pixels that are positive throughout the calibration
interval. The discontinuity therefore cannot be explained by EOG background
zeroing alone.

Independent-brightness and segmented 2013-2019 diagnostics showed that the
effect is predominantly additive rather than multiplicative in the
dim-to-moderate radiance range that is most sensitive to the zero-point change.
A pixel offset of **0.235 nW cm-2 sr-1** was selected: this is the rounded
0.23478 pass-2 estimate and gave the best residual performance through
independent baseline radiance <=10 nW in the pass-3 candidate comparison.

For every 2017+ raster pixel:

```text
raw <= 0 : corrected = raw
raw > 0  : corrected = max(raw - 0.235, 0)
```

This is applied to all NUTS3 countries and territories in the analysis. It is
not country-specific and it is not a subtraction from the NUTS3 mean.

The correction is generated on-the-fly by a custom exactextract RasterSource.
The raw EOG raster is never modified and no duplicate corrected global raster
is stored. The transformed pixels are aggregated using exactly the same
polygon-overlap and spherical physical-area weighting as the raw production
mean.

For each NUTS3/year the calibrated table retains the raw mean, corrected mean,
raw-positive area fraction, area fraction clipped to zero, effective NUTS3
offset, effective percentage correction and calibration provenance. The stage
also re-extracts the unmodified raw mean through the same RasterSource path and
fails if it does not reproduce the existing production raw zonal mean within
numerical tolerance.

The former DE/IT/NL/FR/ES darkest-20% correction remains represented only in
legacy diagnostic configuration so the earlier calibration experiments are
reproducible. It is no longer the production correction.


16. **country-facets**
   - generate a multi-page publication PDF with six countries per page
   - order countries alphabetically by country code
   - show all available years in every country panel
   - plot corrected NUTS3 radiance versus population density in log10-log10 space
   - use common axes across all pages for direct comparison

17. **mixed-model**
   - fit a Gaussian linear mixed-effects model to corrected log10 radiance
   - use NUTS3 as a random intercept
   - include fixed effects of log10 population density, centred year and country
   - include country-specific population-density slopes and country-specific
     temporal trends
   - fit nested ML models and likelihood-ratio tests for the added effects
   - report the NUTS3 random-intercept variance and ICC
   - additionally test population density × year and the
     population density × year × country interaction
   - report Nakagawa-style marginal and conditional R² for every nested model
   - retain Nakagawa marginal and conditional R² as whole-model summaries
   - report likelihood-ratio pseudo-R² for each added effect block using
     R²_LR = 1 - exp(-LR / n)
   - report the corresponding likelihood-ratio f² = R²_LR / (1 - R²_LR)
   - retain Δ marginal R² only as a diagnostic because it is not guaranteed to
     be monotonic across refitted mixed models
   - fit a standardized full model for directly comparable continuous-effect
     coefficients
   - exclude countries represented by fewer than three unique NUTS3 regions


18. **europe-maps**
   - generate one high-resolution PNG for every analysis year
   - map `radiance_mean_corrected` only; raw radiance is never used
   - use the Eurostat-aligned NUTS release configured for each year
   - crop to a European longitude/latitude extent to exclude overseas territories
   - project to ETRS89 / LAEA Europe (EPSG:3035)
   - use one common logarithmic colour scale across 2013-2024
   - draw NUTS3 boundaries and stronger national outlines
   - write the common colour-scale metadata alongside the figures


19. **socioeconomic** *(future-population-models branch)*
   - download Eurostat NUTS3 GDP dataset `nama_10r_3gdp`
   - retain total GDP and GDP per capita in euro and PPS
   - use Eurostat `reg_area3` total land area to derive GDP density
   - create log10 and descriptive year-on-year socioeconomic metrics
   - preserve source/status fields and report annual NUTS3 coverage

20. **socioeconomic-merge** *(future-population-models branch)*
   - merge socioeconomic metrics onto the corrected radiance/population table
   - retain the existing production table unchanged
   - write a separate socioeconomic analysis dataset

21. **socioeconomic-model** *(future-population-models branch)*
   - start with a biologically interpretable maximal fixed-effects model
     containing all main effects, all two-way interactions and all four
     possible three-way interactions among population density, GDP per capita,
     year and country
   - do not include a four-way interaction
   - restrict this model-selection analysis to 2013-2023 because 2024 regional
     GDP coverage is incomplete
   - retain `(1 | NUTS_ID)` throughout
   - simplify top-down from highest-order interactions to the minimum adequate
     model using nested ML likelihood-ratio deletion tests
   - preserve marginality: lower-order components of retained interactions are
     not eligible for deletion
   - write the complete term-by-term simplification history


## Annual VNL V2.2 retained-pixel 2017 discontinuity diagnostic

Before replacing the provisional production calibration, the pipeline now tests
the Annual VNL V2.2 `average_masked` rasters directly.

The key methodological question is whether EOG's zero-valued background mask
already removes the 2017 DNB calibration discontinuity, or whether a step
remains in the radiance of pixels that are persistently retained as lit.

Run:

```bash
python src/run_pipeline.py annual-vnl-2017-pixel-test
```

This diagnostic uses fixed NUTS2016 polygons and the 2015, 2016, 2017 and 2018
raw Annual VNL V2.2 rasters. It performs two separate checks.

**Mask/support check.** Within each NUTS3 it measures the fraction of common
valid pixels that are positive in each year and the fractions switching
zero-to-positive or positive-to-zero across 2015-2016, 2016-2017 and
2017-2018. An unusual 2016-2017 change here would indicate that support/masking
contributes to the discontinuity.

**Persistent-lit radiance check.** The primary analysis then keeps only pixels
that are positive in all four years. For each such pixel:

```text
d15_16 = R2016 - R2015
d16_17 = R2017 - R2016
d17_18 = R2018 - R2017

excess_2017 =
    d16_17 - 0.5 * (d15_16 + d17_18)
```

Because every tested pixel is positive in every year, a positive 2017 excess
cannot be produced merely by pixels moving from the zero/background class into
the lit class.

NUTS3 pixel statistics are summarized by country using medians of NUTS3
statistics. This explicitly tests the previously observed fact that the
2016-to-2017 radiance jump appeared consistently country-to-country rather than
only in the original five comparison countries.

The same retained-pixel diagnostic is stratified by 2016 pixel radiance. If the
2017 excess is approximately constant in absolute linear radiance units across
brightness classes, a simple additive offset remains plausible for retained
pixels. Strong brightness dependence would argue against one universal
additive value.

Finally, the script calculates the same 2015-2018 excess contrast from the
existing raw NUTS3 zonal means and compares it with the direct persistent-lit
pixel contrast. This connects the new product-level diagnostic to the earlier
country/NUTS3 discontinuity result.

Outputs:

```text
data/processed/annual_vnl_2017_pixel_nuts3.csv
data/processed/annual_vnl_2017_pixel_country.csv
data/processed/annual_vnl_2017_pixel_nuts3_brightness.csv
data/processed/annual_vnl_2017_pixel_country_brightness.csv
data/processed/annual_vnl_2017_pixel_brightness_summary.csv
data/processed/annual_vnl_2017_pixel_summary.txt

figures/annual_change/annual_vnl_2017_pixel_test/country_transition_medians.png
figures/annual_change/annual_vnl_2017_pixel_test/country_2017_excess.png
figures/annual_change/annual_vnl_2017_pixel_test/brightness_dependence.png
figures/annual_change/annual_vnl_2017_pixel_test/nuts3_zonal_vs_persistent_pixel_excess.png
```

This stage is diagnostic only. It does not modify
`radiance_mean_corrected` or any production raster.


## Annual VNL V2.2 calibration pass 2: independent brightness and candidate correction

The first retained-pixel test established that the 2016-to-2017 jump remains in
pixels that are positive in all of 2015-2018, and therefore cannot be explained
by EOG background zeroing alone. Its original brightness stratification was,
however, mathematically coupled because 2016 radiance both defined the
brightness bin and appeared in the adjacent transition differences.

The second-pass diagnostic fixes that problem:

```bash
python src/run_pipeline.py annual-vnl-2017-correction-pass2
```

Pixels must be positive in every year from 2013 through 2018. Their independent
brightness reference is:

```text
baseline = mean(R2013, R2014)
```

The 2017 excess remains:

```text
excess_2017 =
    (R2017 - R2016)
    - 0.5 * [(R2016 - R2015) + (R2018 - R2017)]
```

The stage reports country and brightness-stratified offsets without using any of
2015-2018 to define brightness.

It then evaluates candidate **pixel-level** additive corrections. A candidate
offset is subtracted only from positive 2017 and 2018 pixels and values are
clipped at zero for this diagnostic. Candidate values include:

- no correction;
- literature sensitivity values 0.125 and 0.15 nW cm-2 sr-1;
- the median country retained-pixel excess;
- a dim-to-moderate empirical value obtained by first taking each country's
  median independent-bin excess for baseline radiance <=10 nW cm-2 sr-1, then
  taking the median across countries.

For each scenario the 2017 residual excess is recomputed after correction. The
main comparison is the median absolute country residual, supplemented by
country IQR, clipping fraction and residual brightness structure. An empirical
scenario is not accepted merely because its median residual is zero by
construction.

Outputs:

```text
data/processed/annual_vnl_2017_pass2_nuts3.csv
data/processed/annual_vnl_2017_pass2_nuts3_brightness.csv
data/processed/annual_vnl_2017_pass2_country.csv
data/processed/annual_vnl_2017_pass2_country_brightness.csv
data/processed/annual_vnl_2017_pass2_brightness_summary.csv
data/processed/annual_vnl_2017_pass2_candidate_offsets.csv
data/processed/annual_vnl_2017_pass2_scenario_nuts3.csv
data/processed/annual_vnl_2017_pass2_scenario_country.csv
data/processed/annual_vnl_2017_pass2_scenario_summary.csv
data/processed/annual_vnl_2017_pass2_scenario_country_brightness.csv
data/processed/annual_vnl_2017_pass2_scenario_brightness.csv
data/processed/annual_vnl_2017_correction_pass2_summary.txt

figures/annual_change/annual_vnl_2017_correction_pass2/
    independent_brightness_dependence.png
    candidate_country_residuals.png
    candidate_brightness_residuals.png
```

This stage remains diagnostic. If one constant pixel correction performs well,
the next step is to apply it to the full 2017-2024 rasters *diagnostically* and
re-extract NUTS3 means using the production exact-overlap spherical-area
weighting before changing the production calibration.


## Annual VNL V2.2 calibration pass 3: segmented 2013-2019 time series

Pass 2 showed that the independent-brightness 2017 excess is remarkably flat
at about 0.23-0.24 nW cm-2 sr-1 from the faintest retained pixels through a
2013-2014 baseline of roughly 5 nW cm-2 sr-1, but the pattern breaks down at
higher brightness.

The third pass estimates the discontinuity from the full time series rather
than only adjacent year differences:

```bash
python src/run_pipeline.py annual-vnl-2017-correction-pass3
```

The primary interval is 2013-2019 on fixed NUTS2016 geography. 2020 is excluded
from the calibration estimate so the COVID-era year cannot influence the
inferred sensor step. Pixels must be positive in every year of the interval and
are assigned to brightness bins from the independent mean 2013-2014 radiance.

Within every NUTS3 x brightness-bin series the annual median radiance is fitted
as:

```text
R_t = a + b * (year - 2016) + delta * I(year >= 2017)
```

where `delta` is the additive step. A sensitivity model allows a separate
post-2017 slope, and a parallel log10 model estimates the corresponding
multiplicative percentage step. Comparing brightness dependence of the additive
and multiplicative estimates tests which calibration form is more plausible.

Candidate correction shapes are then evaluated against the segmented step:

- no correction;
- universal 0.15 nW;
- universal 0.23478 nW from pass 2;
- universal empirical pass-3 offset estimated from baseline <=5 nW;
- empirical offset only for baseline <=5 nW;
- empirical offset flat through 5 nW and tapering linearly to zero at 10 nW.

The candidate ranking prioritises removal of the residual step through baseline
<=10 nW while also penalising poor residual performance across the full
brightness range. This remains diagnostic only.

Outputs:

```text
data/processed/annual_vnl_2017_pass3_nuts3_brightness.csv
data/processed/annual_vnl_2017_pass3_country_brightness.csv
data/processed/annual_vnl_2017_pass3_brightness_summary.csv
data/processed/annual_vnl_2017_pass3_empirical_country.csv
data/processed/annual_vnl_2017_pass3_scenario_detail.csv
data/processed/annual_vnl_2017_pass3_scenario_summary.csv
data/processed/annual_vnl_2017_correction_pass3_summary.txt

figures/annual_change/annual_vnl_2017_correction_pass3/
    segmented_additive_step_by_brightness.png
    additive_vs_multiplicative_step.png
    candidate_residuals_by_brightness.png
```

If the same simple correction form wins this longer-window test, the next step
is a diagnostic 2017-2024 raster correction followed by full NUTS3 re-extraction
and re-running the population-density rotation analysis before replacing the
production calibration.


## ISS colour ratios and the dense-area LED-switch hypothesis

With the current Europe-wide pixel calibration, the long-term VIIRS trend is
positive below 100 people km-2, near zero at 100-300, and negative above 300.
A plausible mechanism is that denser landscapes have more mature illuminated
extent while conversion from legacy lamps to whiter LEDs reduces the VIIRS DNB
response.

Run:

```bash
python src/run_pipeline.py iss-led-switch-density-test
```

This stage reuses the paired-common-support ISS spectral table from
`iss-spectral-validation` and the current production-calibrated 2013 and 2020
NUTS3 VIIRS values.

The primary mechanistic conversion uses paired Δ(G/R). From Sánchez de Miguel
et al. (2022):

```text
G_ISS / VIIRS = 0.21 * 1.5^(G_ISS / R_ISS)
```

so, at equal green-band output:

```text
predicted spectral Δlog10(VIIRS)
    = -log10(1.5) * Δ(G/R)
```

Positive Δ(G/R), i.e. spectral whitening, therefore predicts a negative VIIRS
component.

The analysis asks whether whitening increases with population density, whether
the predicted spectral component tracks the observed calibrated VIIRS endpoint
change, and whether removing that component reduces the dense-versus-sparse
contrast.

The ISS post mosaic is a 2014-2020 composite, so this remains an endpoint
mechanism test rather than a year-resolved LED adoption model. Equal-green
normalization is not equivalent to constant lumens, electrical power, or
biological exposure.

Outputs:

```text
data/processed/iss_led_switch_density_nuts3.csv
data/processed/iss_led_switch_density_band_summary.csv
data/processed/iss_led_switch_density_models.csv
data/processed/iss_led_switch_density_summary.txt

figures/population_relationship/spectral_validation/led_switch_density/
    density_band_led_switch_decomposition.png
```
