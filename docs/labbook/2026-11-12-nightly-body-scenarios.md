# Nightly spider PRC scenarios — 12 November 2026

Implemented `src/model_spider_nightly_prc.py` for Cottage and Oxford Business Park. It uses anchor EXIF dates and confirmed local time zones plus verified site coordinates to calculate each noon-to-noon lunar track. Capture-date flags, including Cottage August EXIF versus September reference, are carried forward. Moon positions are compass-referenced, so no camera-north inference is needed for the forward model. No photograph fit or historical PRC fit is performed.

The 60-degree full-width field is fixed to body posture. Tilt is from zenith towards body heading; defaults 0,30,60 degrees. Heading is 0–350 degrees from compass north in 10-degree steps. Samples are uniform in solid angle inside the cap. Orthogonal analyzer axes are parallel transported from its centre to every ray, rather than summing image-plane Stokes in incompatible tangent frames. Body longitudinal axis defines the first eye analyzer; anatomical offset is assumed, not established.

Berry eq4.2 unity-peak pattern plus Wang horizon factor is retained as a natural DoLP envelope. Atmospheric depolarization sensitivity factors are 1,0.8,0.6; these multiply Q/U only. They are not measured turbidity, and cannot be independently identified from ALAN with a single monochromatic DoLP snapshot. No Perez spectral turbidity or Wang atmospheric exponential is restored.

ALAN is unpolarized, spatially uniform absolute radiance. Natural radiance is also uniform and constant in this pilot. Ratios 0,0.5,1,3 are illustrative band-integral ALAN/natural ratios, shared between sites for sensitivity, not estimates for those sites. Moon movement changes the polarization geometry only: phase-dependent brightness, angular natural radiance distribution, lunar extinction and below-horizon diffuse lunar light are not modelled. Predictions are restricted to Sun below -6 degrees and Moon above 0.1 and below 89.9 degrees. No twilight polarization or terrain/cloud mask. These restrictions leave gray gaps in heatmaps. Nightly absolute PRC ranges are min/median/max over uniformly sampled headings, not confidence intervals or behavioral probabilities.

## Run on EcoVis

```bash
python src/model_spider_nightly_prc.py \
  --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
  --output data/processed/spider_nightly_scenarios_01
```

Outputs: signed PRC heatmaps per site/atmospheric factor, absolute PRC median/range plots, heading_responses.csv, prc_ranges.csv and provenance.json. This needs Astropy, NumPy, SciPy and matplotlib. No new image registration.

## Spectral input

Add `--spectra path/to/verified_spectra.csv`. Required columns: wavelength_nm,natural,alan,sensitivity. Optional columns: dop_scale (0..1, default1), abs_parallel (default0.4), abs_perpendicular (default0.04). Spectra must share strictly increasing wavelengths. Each radiance spectrum is normalized separately to unit unweighted band integral; ALAN ratios then refer to that wavelength band. Sensitivity can use arbitrary common units because contrast is scale invariant. Positive overlap with the eye response is required. No extrapolation is performed.

The two receptor sums are integrated over wavelength and solid angle before taking the signed natural-log ratio. Spectral weighting uses wavelength-dependent absorptance 1-10^(-absorbance), not an average DoLP. Input spectral ordinates must be energy radiance times a sensitivity with energy-response weighting, or photon radiance times photon-response sensitivity consistently. Convert energy to photons before supplying a photon-based sensitivity; the code does not infer units. Wavelength-dependent natural DoLP is separable as dop_scale and cannot exceed the unity-peak base. Spectra and response functions are constant in time and across sky directions in this pilot.

Without --spectra, outputs explicitly state ACHROMATIC_SCENARIO_NOT_SPECIES_SPECTRAL_MODEL. No invented HPS, LED or spider spectral sensitivity is supplied. Verified natural spectra, HPS/LED skyglow spectra and species sensitivity must be brought together for the spectral comparisons. Raw lamp spectra alone need not equal the skyglow spectrum at the animal.

Validation: 29 numerical tests across Wang/Berry/night model. New checks cover unit tangent axes, fixed cap, rotation covariance, dilution monotonicity and spectral normalization. Both sites completed an end-to-end coarse ephemeris run. No empirical ALAN values inferred.
