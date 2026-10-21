# Freeze: Jan1 population × VIIRS radiance × locked spider/moth PRC

**Investigator instruction: 21 October 2026.** Freeze the *current* observational multiscale comparison and return later specifically for choices of manuscript panels, layout and table/figure details. Do **not** re-run, refit, revise source products or alter interpretation merely to seek larger correlations.

The frozen analysis means: the fixed NUTS2016 2014–2020 cohort of **1,273** regions using January-1 Eurostat population counts over fixed land area, calibrated mean NUTS3 VIIRS radiance, locked pixel-transformed spider17/moth16 PRCs, country- and between-country correlations, NL333 and near-neighbour tests, and their existing caveats.

Versioned reference documents:
- `docs/results_population_jan1_fixed_nuts2016_2014_2020.md`
- `docs/results_jan1_country_radiance_prc_nl333_2014_2020.md`
- `docs/country_mechanism_typology_population_prc_2014_2020.md`
- `docs/lab_notebook/2026-10-21.md`

Frozen means scientific workflow and datasets; it does **not** mean unresolved NUTS geography provenance, LED spectral sensitivity, spatial autocorrelation or exact source attribution are somehow validated. The physical lighting-operation evidence for greenhouse horticulture in NL333 remains the documented contextual explanation, not a precisely estimated contribution to VIIRS trend. Future spectral tests operate in new files and never overwrite any of these outputs.

# Resumed question: ISS spectral whitening and the apparent 2014–2024 VIIRS decline

**Research question:** Across fixed NUTS2016 regions and countries, do independently measured ISS RGB colour shifts correspond to lower *locked zero-offset-corrected VIIRS DNB radiance*? Which apparent declines plausibly reflect spectral migration beyond VIIRS blue sensitivity, versus independent evidence for true source dimming/switch-off?

## What already exists

- `src/validate_iss_spectral_shift.py` already downloads/validates the original pair of Europe RGB ISS mosaics (Zenodo doi:10.5281/zenodo.7677478), applies Sánchez de Miguel et al. 2022 common-pixel ratio and quality filters, and writes `data/processed/iss_nuts3_spectral_shift.csv`, with `delta_gr_median`, `delta_bg_median` and `n_common_valid_pixels`.
- `src/test_iss_led_switch_density.py` and `src/test_iss_led_switch_annual_slopes.py` analyse the earlier population-dependent slope explanation at 2013→2020; they are not the 2014→2024 geographical dimming test.
- The **locked post-2017 background offset 0.235 nW** is already in `data/processed/fixed_nuts2016_viirs_2014_2024.parquet` as `radiance_mean_corrected`; never subtract it twice.
- **Locked 2021 native-pixel snow correction** `config/viirs_2021_snow_correction.yaml` is applied to the published spider/moth 2014–2024 PRC histories, but the NUTS3 *radiance mean parquet* has **not been verified** to include the same 2021 patch. This matters for full-year or year-by-year radiance analyses; initial spectral comparison deliberately uses ONLY 2014, 2020, and 2024 endpoints, so no 2021 provenance ambiguity affects computed endpoint differences.

## New read-only first analysis

`src/analyze_iss_spectral_vs_locked_viirs.py` joins existing paired, **common-valid-pixel ISS** spectral ratios onto corrected VIIRS NUTS2016 endpoints. It estimates (i) per-country Spearman correlations between ΔISS B/G or G/R and Δlog10 VIIRS for 2014→2020 and 2014→2024, (ii) equal-country-weight comparisons of national **regional-median** spectral and radiance changes, (iii) sign coexistence with transparent near-zero spectral-ratio sensitivity, and (iv) an exploratory 2020→2024 radiance comparison as **negative-timing-control/sensitivity**, since ISS post epoch ends in 2020. Reports per-country eligible-region and minimum-pixel counts; no artificial imputation or remapping.

**Temporal limitation (essential):** The pre-ISS mosaic is 2012–13; post-ISS is a heterogeneous 2014–20 composite. The ΔISS colour shift cannot be dated to 2014→2020, and certainly not 2014→2024. Therefore all spectral-vs-radiance associations are **geographic mechanistic consistency comparisons**, not exact contemporaneous trends. Physical reasons for no strong or weak association include mixed source replacements, operational dimming, sensor masks and different geographical supports.

**Illustrative spectral physics only:** The prior literature-linked empirical response `G_ISS / DNB = 0.21 × 1.5^(G_ISS/R_ISS)` implies Δlog10(DNB) = −log10(1.5) × Δ(G/R) at *equal ISS green-band output*. The script stores this as a labelled illustrative comparison; it is NOT used as a correction to source radiance, emitted broadband light or the locked PRC functions.

**Crucial inference:** strong ISS colour shift with dimming supports a *possible spectral* source of apparent dimming. VIIRS dimming without measured ISS colour shift, especially alongside independent curfew/electricity evidence, supports *possible true operational reduction*. Neither distinguishes lighting power from spectra conclusively. French policy/Enedis, local UK part-night lighting and Dutch greenhouse-energy evidence are potential independent checks. See `docs/literature_europe_nighttime_dimming_policy_led_2014_2024.md`.

## Run on EcoVis

```bash
git switch historical-polarization-visibility
git pull --ff-only
python src/analyze_iss_spectral_vs_locked_viirs.py --self-test
python src/analyze_iss_spectral_vs_locked_viirs.py
```

If `data/processed/iss_nuts3_spectral_shift.csv` does not exist, inspect prior ISS outputs first, then run the existing `python src/run_pipeline.py iss-spectral-validation` pipeline to produce it; this can download large ISS GeoTIFF files. No reason to re-download if files already exist.

New output CSV directory `data/processed/iss_spectral_vs_corrected_viirs_2014_2024/`:
- `iss_viirs_matched_nuts3.csv`
- `iss_viirs_within_country_associations.csv`
- `iss_viirs_country_medians.csv`
- `iss_viirs_between_country_associations.csv`
- `iss_viirs_sign_sensitivity.csv`
- `iss_viirs_analysis_manifest.json`

Figures in `figures/iss_spectral_vs_corrected_viirs_2014_2024/`: `iss_spectral_shift_vs_corrected_viirs_scatter.pdf`, `iss_spectral_shift_vs_viirs_country_medians.pdf`.

**Status: scripted, NOT YET EXECUTED on EcoVis.** Inspect self-test and regional/count coverage before interpreting estimates. Next investigate whether *actual 2021-corrected native VIIRS radiance* is available as a separate annual product and can replace only 2021 in a reference full-series comparison without overwriting historical source parquet. This is not required to answer the endpoint spectral association.

**Manuscript relevance:** separate evidence for a DNB-band change, evidence of source spectral conversion, and independent evidence for actual operational reductions. Never infer a unique spectral-to-total-light conversion from two ISS colour ratios or directly pass 'LED-corrected radiance' through the locked empirical DNB-to-PRC transfer without validation.