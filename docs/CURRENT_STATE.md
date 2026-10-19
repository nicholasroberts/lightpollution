# Current research state / restart point

**Updated:** 2026-10-19  
**Branch:** `historical-polarization-visibility`  
**Machine:** EcoVis  
**Repository root:** \`/run/media/bznwr/data/lightpollution\`

## Next analysis stage — 2014–2024 VIIRS pixel PRC

`python src/historical_prc_pixel_nuts3.py --pilot` for an initial 2014/2018/2019/2021/2024 check, then `python src/run_pipeline.py historical-prc` for all 11 years. Native EOG pixels are corrected and mapped through the selected PRC curve **before** aggregation to fixed NUTS2016 polygons. Annual and cumulative PRC absolute/% maps, Europe summaries and NUTS3 rankings are generated. Optional native-resolution five-band pixel GeoTIFF: `python src/historical_prc_pixel_nuts3.py --no-maps --pixel-year 2024`. See [analysis instructions](historical_prc_analysis.md). **No real mapping outputs generated yet.** 2021 snow adjustment remains unresolved at pixel level; 2021 annual changes are diagnostic.

---

## LOCKED WORKING PRC FUNCTION — 19 October 2026

**Adopted provisionally for historical 2014–2024 analysis:** `PRC(R) = 0.144845 + 0.316330 * exp(-0.106483 * R)`, where `R` is the *project-calibrated* annual EOG VIIRS pixel radiance (nW cm^-2 sr^-1). Original R `SSasymp` functional form, refitted using **observation-year radiance** (2018 for three sites; 2019 for eleven). `n=16` provisional matches across the 14 Table S1 sites; `R²=0.949279`, `RMSE=0.027146`. Original source values and original fit are secondary comparisons. Parameters: `config/prc_transfer.yaml`. Research context and caveats: [19 October laboratory notebook](lab_notebook/2026-10-19.md).

**Still unresolved:** two original plotted PRC observations have no confirmed site, four assignments remain provisional, and the original R file mapping is pending. Do not interpret PRC as a validated behavioural threshold. European maps have not yet been calculated; the NUTS3-level special 2021 correction is not a per-pixel map correction.

---

## NEW ACTIVE RESTART — 19 October 2026

**Branch now created:** `historical-polarization-visibility`. **Study:** observed 2014–2024 changes in locally modelled celestial polarization-receptor contrast (PRC) across Europe; forecasting remains later work. See [19 October lab notebook](lab_notebook/2026-10-19.md).

**Completed first check:** all 14 experimental sites extracted from local EOG 2018 and 2019 **VNL V2.1** annual `average_masked` VIIRS rasters; the existing post-2017 **0.235 pixel calibration** has been applied. Note the V2.1 filenames despite earlier planning descriptions saying V2.2. Original Light Pollution Map values came from the **2019 VIIRS layer**, `viirs_2019`, which is spatially smoothed relative to pixel samples.

**Historical interim decision (superseded):** 2018 annual EOG at all sites was briefly proposed, then replaced by **observation-year calibrated EOG**: 2018 at the three 2018 field sites and 2019 at the eleven 2019 sites. This updated dataset was used for the provisional 16-observation SSasymp refit adopted above. The 18-entry R `sigg` vector is available but original file-to-site mapping is still pending.

**Next action:** prepare historical per-pixel radiance/PRC reconstruction with the adopted function, while separately recovering the original R filename-to-site mapping and auditing calibration uncertainty. No Black Marble work is needed at this stage. The historical European study remains **2014–2024**.

---

## Updated primary direction — evening 16 October 2026

**Historical polarization-cue visibility across Europe (2014–2024) is now the first-paper objective**, rather than SSP forecasting. Begin by reproducing the original PRC relationship from site-level field data and **Light Pollution Map's `viirs_2019` radiance layer**; do not confuse it with modelled sky brightness. Re-extract 2019 EOG and Black Marble site radiance, then test annual, date-matched monthly and quality-screened daily measures.

**Pipeline and decision diagram:** [`docs/historical_prc_pipeline.md`](historical_prc_pipeline.md). Full reasoning and prior figure/dataset references remain in the [16 October lab book](lab_notebook/2026-10-16.md). Primary historical PRC mapping and a dedicated branch are **not yet implemented**; the latter is proposed for Monday. Preserve the older annual NUTS3 forecast and NL333 diagnostics as previous results.

---

## Previous restart state, retained for provenance

## CURRENT PRIMARY DECISION — 2026-10-16

Read the full dated [2026-10-16 laboratory notebook](lab_notebook/2026-10-16.md) for detailed tables, scientific reasoning, all figure paths, data provenance and commit references.

**Primary modelling input:** 1,514 NUTS2016 regions × every annual observation from 2014 to 2024, with **both** the 2017+ pixel calibration and **2021 log-ratio snow/surface correction**. Specifically:
\`data/processed/fixed_nuts2016_2021_surface_correction/fixed_nuts2016_viirs_2014_2024_surface_corrected.parquet\`, field \`radiance_mean_corrected_surface_logratio_primary\`.

**Working predictor structure:** density-specific longitudinal time slopes, country baseline intercept and NUTS3 random intercept. **H removed from the primary forecasting model** after small out-of-sample gains and failure of national H composition to explain most country trends. Country-specific temporal slopes are diagnostics, not automatically adopted mechanistic predictors.

**NL333:** exceptionally bright Delft en Westland region fell from 275.87 to 96.33 nW cm^-2 sr^-1 (2014→2024). In earlier *pre-snow* holdouts it dominated (~87%) squared absolute errors. Its **year-by-year** cause/timing remains unverified. Do not delete it or claim a greenhouse explanation without a proper audit.

**Next stage committed but not yet run on EcoVis:**

\`\`\`bash
python src/run_pipeline.py future-radiance-corrected-annual-nl333
\`\`\`

**Outputs to inspect first:** \`data/processed/fixed_nuts2016_corrected_annual_nl333/nl333_annual_diagnosis.csv\`, \`corrected_annual_nl333_summary.txt\`, and \`figures/population_relationship/fixed_nuts2016_corrected_annual_nl333/nl333_full_annual_radiance.pdf\`, \`observed_vs_modelled_every_year_by_density.pdf\`.

Endpoint maps/boxplots from \`future-radiance-corrected-endpoint-change\` were completed successfully and remain model-free ground truth, but endpoint comparisons are **not** a substitute for fitting the full annual panel. Pre-2021-correction H/country analyses remain valuable *supplementary negative results*, not the primary fully corrected analysis.

---

## Archived 2026-10-14 restart notes (superseded by the update above)


**Updated:** 2026-10-14  
**Branch:** `future-population-models`  
**Machine:** EcoVis  
**Repository root:** `/run/media/bznwr/data/lightpollution`

## Overarching objective: historical reality and model response must converge

The immediate snow/albedo work is subordinate to the main scientific problem:
the future-radiance model must reproduce the historical density-dependent
radiance behaviour before it is used for projection.

The key diagnostic figure is:

```text
figures/population_relationship/future_radiance_historical_validation/
    historical_nuts3_slope_distributions.pdf
```

It compares the distribution of observed fixed-NUTS2016 2013-2020 NUTS3
radiance trends in each population-density class with the canonical
density-class estimates derived from the longer 2013-2024 modelling framework.

The important feature is not merely whether the canonical diamonds lie within
the broad observed distributions. Their *density-dependent pattern* should agree
with the historical class behaviour. At present the figure shows a systematic
tension: the canonical estimates are more positive in sparse classes and
negative in dense classes, whereas the central tendency of the observed
2013-2020 distributions is approximately flat/negative in sparse classes and
more positive in dense classes.

Therefore the current analysis sequence is:

1. identify and explain non-secular disturbances that can rotate the estimated
   density response, especially the 2021 winter pulse;
2. quantify whether snow/surface-albedo effects account for those disturbances;
3. construct a historically defensible radiance series / model specification;
4. refit the density-dependent temporal response;
5. regenerate `historical_nuts3_slope_distributions.pdf` using the final model;
6. require qualitative and quantitative agreement between observed historical
   class behaviour and model-predicted class trends before using the response
   for future projection.

The albedo attenuation test is therefore not an end in itself. Its role is to
determine whether the 2021 pulse is an environmental observation effect that
should be modelled/controlled rather than interpreted as secular change in
artificial-light emissions.

## Current conclusion

The corrected MOD10CM snow test has completed and is valid.

The February 2021 anomaly shows:

```text
region Spearman rho = +0.6653
country-median rho = +0.7383
adjusted snow-block partial R2 = 0.2805
adjusted snow-block joint p = 2.48e-14
snow x baseline-brightness p = 1.94e-13
```

January independently supports the relationship; March is a practically null
negative control.

Snow/surface reflectance is therefore strongly supported as a major contributor
to the January-February 2021 VIIRS pulse.

The full 2020-2022 within-region seasonal relationship is not generally
positive, so the conclusion is specific to the winter 2021 event rather than a
universal European snow-brightness rule.

## Executable albedo stage

The focused albedo attenuation test is now committed and registered:

```bash
python src/run_pipeline.py future-radiance-albedo-2021-test
```

Place all three AppEEARS January-March bundles (2020, 2021, 2022) somewhere
under:

```text
data/raw/modis/VNP43MA3_002/
```

ZIP bundles can remain zipped. The stage will refuse to run a partial
2020-2022 inventory.

Primary interpretation remains absolute nW radiance; relative percentage
anomaly is reported alongside it specifically to prevent small-denominator
percentage changes in dark regions from being mistaken for large physical
radiance changes.

## Next scientific task

Follow the published seasonal-night-light framework one step further by adding:

- MODIS surface albedo, with a visible broadband measure as the primary optical
  variable and white-sky/black-sky variants retained;
- monthly MODIS NDVI as a complementary seasonal surface-state control.

Preserve the existing pre-specification:

```text
February primary
January supportive
March negative control
2021 - mean(2020, 2022)
country + latitude + baseline brightness + cf_cvg controls
surface-state x baseline-brightness interaction
```

Primary mechanistic question:

> Does the albedo anomaly predict the VIIRS anomaly and attenuate the snow term?

This is a mechanistic attenuation test, not formal causal mediation.

## Key validated snow outputs

```text
data/processed/fixed_nuts2016_snow_viirs_test_summary.txt
data/processed/fixed_nuts2016_snow_2021_anomaly_tests.csv
data/processed/fixed_nuts2016_snow_2021_anomaly_models.csv
data/processed/fixed_nuts2016_snow_levin_jan_july_summary.csv
data/processed/fixed_nuts2016_snow_levin_percent_lit_thresholds.csv
data/processed/fixed_nuts2016_snow_levin_seasonal_correlations.csv

figures/population_relationship/fixed_nuts2016_snow_viirs_test/
    february_snow_vs_viirs_anomaly.pdf
    snow_viirs_anomaly_rho_jan_feb_mar.pdf
    levin_style_jan_july_snow_radiance_2021.pdf
    levin_style_percent_lit_thresholds.pdf
```

See:

- `docs/investigations/viirs_2021_anomaly.md`
- `docs/lab_notebook/2026-10-14.md`

### Albedo first-run implementation note

The first execution of `future-radiance-albedo-2021-test` exposed a CRS bug
before any scientific inference: native sinusoidal VNP43MA3 rasters were zonally
extracted against EPSG:4326 NUTS polygons. This produced no eligible albedo
models and the run stopped at plotting. Commit `64ed94c` fixes the extraction,
invalidates only the bad zonal caches, and preserves the expensive valid monthly
composites. No failed-run numerical result is used.

### Albedo mechanism result

The corrected VNP43MA3 attenuation stage now strongly supports the
January-February 2021 pulse as a surface-reflectance disturbance. Primary
February WSA QA0/1: region rho +0.7198, country-median rho +0.8440; adding
albedo reduces the snow-block partial R2 from 0.2803 to 0.0789 (71.8%
attenuation), while albedo retains independent explanatory power beyond snow
(partial R2 0.0196; p=0.0073 HC3, 0.0121 country-clustered). March is
practically null after adjustment.

This does not yet justify an ad-hoc 2021 correction. The next task is to define
a historically general surface-state control and then revisit the
density-dependent long-term response in both absolute and relative radiance
space.

### Historical winter snow-control stage ready

The complete 2013-2024 MOD10CM archive is now available locally (142/144
monthly granules; documented gaps 2016-02 and 2022-10). A new pre-registered
stage, `future-radiance-historical-snow-control`, is ready to test whether a
general January-February surface-state covariate materially changes the locked
fixed-NUTS2016 2014-2024 density trends.

The stage reports matched before/after trends in both proportional and absolute
radiance space, includes a complete-Jan/Feb sensitivity analysis, and has a
pre-specified stopping rule for whether full historical albedo acquisition is
warranted. No result from this stage has yet been viewed and no annual VIIRS
production value is modified.

### Historical snow-control rerun required after formula audit

The first output from `future-radiance-historical-snow-control` is invalid for
scientific interpretation because the controlled Patsy formula inadvertently
added baseline brightness as part of the tested block. Commit `dcaaca3`
holds baseline brightness constant in both models and tests only snow + snow x
brightness. The pre-registered stopping rule is unchanged. A clean rerun is
required before deciding whether historical albedo acquisition is warranted.

### Corrected historical snow-control primary result

After the formula fix, the intended 2-df snow + snow x baseline-brightness
block is significant in proportional radiance (p=0.000132) but not in absolute
radiance (p=0.151). Despite that within-series explanatory power, the
2014-2024 density-class secular trends move only trivially in proportional
space and the large 2021 annual residual remains mostly intact. Two dark
classes cross the pre-registered normalized absolute materiality threshold only
because very small physical nW shifts are large relative to their low baseline.

The formal pre-registered rule therefore returns
`historical_albedo_warranted=True`, but the corrected complete-Jan/Feb
sensitivity output should be reviewed before initiating a large historical
albedo download.

### Complete-Jan/Feb snow sensitivity

The corrected complete-Jan/Feb sensitivity (14,998 observations; 1,509 NUTS3;
31 countries) strongly confirms the proportional snow/surface-state block
(LR=34.20, df=2, p=3.74e-08; delta AIC=-30.20) but the absolute block is null
(p=0.735; delta AIC=+3.38). No density class shifts materially in either
proportional or absolute trend space under the complete-winter sensitivity.

Thus winter surface state is a real, reproducible source of short-term
proportional VIIRS variation but does not robustly alter the 2014-2024
density-dependent secular trend. Full historical albedo, if pursued, should be
framed as a mechanism/refinement test rather than as the expected explanation
for the historical/model discrepancy.

### Historical surface-state branch closed

The corrected historical snow-control analysis shows that winter surface state
has statistically detectable short-term proportional effects but negligible
practical leverage on the main annual/secular results. For 2021, the
country-balanced proportional residual falls only from 12.7405% to 12.0225%
(~5.6% reduction) after the general snow control, while the long-term
density-class trends are essentially unchanged.

Accordingly, do **not** proceed to a large historical albedo acquisition for the
main trend analysis. Retain the targeted 2021 snow/albedo work as mechanistic
evidence for the winter pulse, but return the main workflow to modelling the
density-dependent anthropogenic radiance trend itself.
