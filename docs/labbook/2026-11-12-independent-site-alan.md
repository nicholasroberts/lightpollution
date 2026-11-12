# Independent Cottage–Oxford ALAN scenarios — 12 November 2026

Investigator correction: do not fit ALAN to measured DoLP or PRC. Establish a site contrast independently, predict HPS/LED responses, then compare experimental maps. No optimizer, DoLP loss or historical PRC target is used in this pipeline.

## Normalization and evidence

The verified reference table has multiyear upward radiance0.36380952381 for Cottage and41.8296153846 for Oxford Business Park. Their ratio114.9766915 is used as a CONDITIONAL skyglow-contrast proxy. This is independent of measured DoLP but is not a calibrated local downward skyglow ratio. Source distribution, spatial convolution, viewing direction, shielding, atmosphere and different observing nights can alter that ratio substantially. A genuinely independent ground-sky ratio can replace it via --skyglow-ratio. The CSV radiance standard deviations are temporal scatter, not error bars for ground-sky conversion; they have not been repurposed as such.

Both sites now share one natural-intensity reference: Cottage's modelled capture-time zenith moonlight. Previously each site's ALAN was separately set equal to its own capture-time natural intensity; that assumption is superseded when --site-alan-config is supplied. The code retains the old mode for reproducibility.

The scenario levels0.3,1,3,10 specify Oxford HPS ALAN450–750nm / Cottage capture-time natural zenith intensity. They are an explicitly uncalibrated absolute-level sensitivity range, not estimates chosen to reproduce Oxford DoLP. Cottage HPS ALAN equals Oxford HPS ALAN divided by114.9766915 in common radiometric units. Changing Moon brightness gives different instantaneous dilution while ALAN remains fixed throughout each night.

## HPS versus LED at matched DNB-overlap source response

NOAA response source: https://ncc.nesdis.noaa.gov/VIIRS/VIIRSSpectralResponseFunctions.php , SNPP Modulated RSR Release1.0, linked file https://ncc.nesdis.noaa.gov/documents/documentation/GT_F1_RSR_DNBLGS_v5_response_orbit6557_ccm.dat . File SHA256, original response samples and negative-response handling are recorded. Negative instrumental response sample(s) are retained in the exported source table and set to zero only for physical integration.

For each author generic lamp source S(lambda), compute K=integral450–750[S(lambda)(550/lambda)^p] / integral_common_support[S(lambda)RSR(lambda)]. Here p=2 is the existing assumed sky-transfer exponent. Source arrays extend300–900nm, so DNB integration uses the overlap with response support,401.99–900nm. Unknown source emission above900nm is not extrapolated. This is a DNB-OVERLAP spectral comparison, not a full operational VIIRS calibration or a verified satellite-to-ground transfer. That limitation is explicit in metadata. Full VIIRS use requires appropriate source spectra beyond900nm and calibration conventions.

LED2958–3969K generic / HPS K ratio is1.462645884 in the retrieved files. Thus at matched DNB-overlap source response the LED scenario has1.463times HPS450–750nm sky energy before applying spider spectral weighting. The resulting receptor contrasts differ; neither spectrum gets a free ALAN amplitude fit. This is not proof of the actual Oxford lamp mix, and does not guarantee the same spectral ordering for other transfer or normalization assumptions.

At scenario level3, local capture-time band ALAN/natural ratios are approximately:

| Spectrum | Cottage | Oxford |
|---|---:|---:|
| HPS |0.02609|1.38285|
| LED group |0.03816|2.02262|

These are scenario settings, not recovered light-pollution measurements. Spectral inputs remain450–750nm; no UV result is claimed.

## Experimental panels

The plotter accepts --batch and adds the saved experimental camera-channel DoLP alongside natural-only, HPS and LED predictions at capture time. It samples the original I,Q,U arrays into a north-up/east-right equisolid display using the reviewed Moon bearing and fitted image centre/radius. It checks the alignment report hash against the saved Moon anchor. No new registration search or rotation fitting is performed. Input arrays are untouched; bilinear resampling is for display only. Source SHA256 hashes are recorded.

Experimental camera-channel DoLP is not spectrally matched to the model's photon-sensitivity-weighted DoLP. Shared0–1colour limits help visual checks but do not turn this into radiometric validation. Terrain, clouds, calibration errors and camera response can explain residual discrepancies. No automatic use of the observed maps to adjust ALAN is implemented.

## EcoVis workflow

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/run_independent_site_alan.py \
 --batch data/processed/spider_prc_all_sites_01 \
 --spectra data/processed/spider_spectral_inputs_01 \
 --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
 --output data/processed/spider_independent_site_alan_01
```

Output folders: alan_configs, hps, led, comparison_level_0.3, comparison_level_1, comparison_level_3, comparison_level_10. Each comparison folder contains visible_spectral_nightly_comparison.png and sky_snapshot_settings.json. The latter records the actual site/spectrum ALAN ratios applied. Raw response tables preserve scenario level separately from applied local ratio. Level1 is not per-site ALAN/natural=1 anymore.

## Verification

38 numerical tests pass across normalization, pigment, lunar brightness, geometry, measured-map reprojection and DNB conversion checks. Both forward sites ran using the new common-reference scenarios. The combined pipeline is end-to-end checked with synthetic measured arrays because the actual EcoVis Stokes arrays are not in this workspace. No claim of agreement with real experimental maps is made before the local run.
