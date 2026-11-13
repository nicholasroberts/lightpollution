# Measured lunar UV/visible shape and single green510 pigment

Investigator request: check candidate lunar spectra, extend accepted inputs, remove the illustrative UV peak from the active pigment proxy, rerun levels5/6/7, and label local mean and median DoLP within the60deg full field (zenith30deg radius).

## Sources and acceptance check

- QASUME Davos April2025 workbook: https://zenodo.org/records/21280043 , DOI10.5281/zenodo.21280043. Use QASUME_Irradiance (accepted absolute observations) and five complete scans starting2025-04-10T21:22:33 through2025-04-11T00:52:45 as published. Times are retained as published; no unverified timezone interpretation is needed to use spectral shape. Valid selected irradiance rows352–549nm: the final550nm row has uncertainty fields but missing irradiance.
- air-LUSI2022 release, DOI10.18434/mds2-3397: https://github.com/usnistgov/air-lusi . Select2022-03-16T07:25:15UTC, signed phase−24.7773deg, closest flight to the approximately22deg QASUME selection. File support381.464–1039.833nm. Use standardized-distance irradiance without further distance correction because only relative shape is used. Author issue5 documents an erroneous multiply instruction in Irradiance metadata; no such multiplication is applied. Issue4 concerns zenith metadata with authors reporting no meaningful impact within assigned uncertainty.

Normalize each spectrum by its mean energy over450–480nm. Take the median of five normalized QASUME spectra. Compare integrals in10nm bins from400–540nm against air-LUSI at matched bin support. Acceptance gates: median absolute difference<=5%, maximum<=15%, positivity and coverage checks. Actual median absolute difference0.48%, maximum3.20%. These are engineering shape checks, not statistical equivalence or independent UV calibration. Median between-scan relative standard deviation352–380nm is3.11%; median reported expanded uncertainty for those samples2.47%. Between-scan variation is not a confidence interval.

Composite uses QASUME below400nm, a linear blend400–420nm, air-LUSI above420nm, sampled352–750nm. No output extrapolation below352nm. Both inputs are TOA direct lunar irradiance at different near-full phases, not ground-sky radiance and not exact full-moon observations. Existing natural lambda^-4 and ALAN lambda^-2 transfers remain explicit sensitivity approximations. Wavelength-dependent UV atmospheric transmission has not been independently retrieved for these sites. The model therefore explores consequences of this source-spectrum extension; it does not establish measured ground-level UV or absolute species PRC.

## Pigment and band bookkeeping

Active sensitivity is the Stavenga alpha formula with lambda_max510nm and no added UV Gaussian or separate UV receptor. This is a phenomenological green-sensitive pigment proxy, not evidence of molecular RH1 identity in this species. Original optional UV parameters remain in the helper for backwards compatibility but default amplitude is zero and the production workflow calls the zero-UV default. The ordinary short-wavelength tail of the alpha curve remains.

Preparation now defaults to352–750nm. Legacy filenames HPS_visible_input.csv and LED_mid_visible_input.csv and visible_spectral_nightly_comparison.png are retained for compatibility; plots and provenance report the actual band. DNB-overlap conversion integrates the SAME supplied sky band, rather than retaining the old450–750 numerator. New LED/HPS sky-band energy factor is1.602 at equal DNB-overlap source response. Source emission in weak lamp UV samples may include noise; no new site lamp-mixture claim is made.

Cottage HPS is retargeted to mean local DoLP0.45 at level5 under the revised spectra, preserving the user's scenario request. Same physical site source scaling applies to LED; LED is not separately fitted. Oxford source scale remains1. Atmospheric factor0.8; horizontal posture; heading0–350deg every10deg;20min time sampling and1024field rays remain as before. No raw experimental arrays are read by adjustment or nightly forward models.

## Recalculation

Both nightly models completed, each writing3564 heading-response records for the two sites at levels5/6/7. FoV snapshot statistics below are local DoLP values on an equal-solid-angle display; they are not DoLP of spatially integrated Stokes or the PRC.

| Level | Site | Source | FoV mean DoLP | FoV median DoLP |
|---|---|---|---:|---:|
|5|Cottage|HPS|0.4500|0.4669|
|5|Cottage|LED|0.2870|0.2926|
|5|Oxford business park|HPS|0.2179|0.2317|
|5|Oxford business park|LED|0.1109|0.1171|
|6|Cottage|HPS|0.4212|0.4361|
|6|Cottage|LED|0.2573|0.2611|
|6|Oxford business park|HPS|0.1952|0.2074|
|6|Oxford business park|LED|0.0959|0.1011|
|7|Cottage|HPS|0.3960|0.4089|
|7|Cottage|LED|0.2331|0.2357|
|7|Oxford business park|HPS|0.1769|0.1878|
|7|Oxford business park|LED|0.0845|0.0889|

Plots label mean/median on all modeled AND experimental DoLP panels. Experimental statistics use valid physical camera-channel display pixels in the same FoV; measured camera spectral response remains distinct. JSON records counts and valid fraction. Saved original arrays and registration are unchanged. Actual experimental NPZ files reside on EcoVis, so measured-panel rendering is performed by the local command below; only modeled panels were regenerated in this execution environment. Visual inspection confirmed overlap chart and map label layout.83tests passed.

## EcoVis workflow: one output root

```bash
git pull --ff-only origin publication-figures
python -m pip install -r requirements-spectral.txt
python src/prepare_spider_spectral_pilot.py \
  --output data/processed/spider_uv510_01/inputs
python src/run_independent_site_alan.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --spectra data/processed/spider_uv510_01/inputs \
  --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
  --cottage-target-dop 0.45 --cottage-target-level 5 \
  --levels 5 6 7 \
  --output data/processed/spider_uv510_01/results
xdg-open data/processed/spider_uv510_01/results/comparison_level_5/visible_spectral_nightly_comparison.png
```

Use a fresh root if this one already exists. The same figure is generated under comparison_level_6 and comparison_level_7. Spectral overlap check: inputs/moon_audit/moon_spectrum_overlap.png. Numerical PRCs: results/hps/{heading_responses.csv,prc_ranges.csv} and results/led/{heading_responses.csv,prc_ranges.csv}. Each comparison folder contains sky_snapshot_settings.json with FoV statistics.
