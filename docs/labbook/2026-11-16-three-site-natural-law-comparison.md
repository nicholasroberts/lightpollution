# Matched-pixel comparison of two natural DoLP laws at three sites

This is the next controlled check on the angular pattern used for spider PRC.
It changes no production model or empirical radiance–PRC relationship.

## Candidate and controlled quantities

Let W be the selected Berry-normalized complex geometry in the current model,
H the existing Wang horizon factor, and gamma angular distance from the Moon.
The current complex fractional polarization is W H. The candidate is

(W H)/(1+cos²gamma).

The positive real denominator preserves the phase (AoP) and all zeros. With
coincident neutral points and no horizon factor, its magnitude tends to the
ideal Rayleigh DoLP sin²gamma/(1+cos²gamma). For finite neutral splitting it is
an explicit hypothesis, not a derivation from multiple scattering or the full
Wang spectral equation17. It is not an extra division by KS91 brightness.

Both predictions use the identical raw Stokes pixels, ephemeris elevations,
corrected Moon bearings, nominal equisolid geometry and optional terrain mask.
The angular audit domain has radius 60deg; the animal FoV remains 60deg FULL
width, radius30deg. Pixels within15deg of the Moon are excluded. The10deg bins
and their centres are retained from the previous audit. Physical-input exclusions
are shared by all curves. No AoP offset needs estimating to compare amplitudes.

Atmospheric multiplier=1 and added ALAN=0 for both theoretical curves. Those are
controlled natural-envelope predictions, not assertions that the observed skies
have no atmosphere-induced depolarization or ALAN. No site amplitude/width/ALAN
parameter is fitted. Original radiance means and SDs are copied into the results.
Only the separate shape row is peak-normalized; absolute curves and zenith-FoV
means/medians remain unscaled.

## Run on EcoVis

Use the corrected anchors from the latest 180deg Skillinge bearing rerun:

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/compare_spider_natural_dop_laws.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --moon-anchors data/processed/spider_dark_site_validation_02/corrected_moon_anchors.json \
  --output data/processed/spider_natural_law_comparison_01
```

Optional --sky-masks DIR must point to the same reviewed SITE.npy boolean sky
masks used in the preceding analysis. If a mask directory is requested, all
three site masks are required. Without masks, terrain contamination remains
explicitly provisional. The three sites are Cottage, Skillinge, slabhouse.
All source-registration hashes must match the anchors before output is written.

Look in **data/processed/spider_natural_law_comparison_01**:

- three_site_natural_law_comparison.png: absolute DoLP above; peak-normalized
  shape below; measurements/current/candidate on the same axes for each site.
- SITE_matched_dop_maps.png: measured, current and candidate maps with exactly
  the same compared pixels and a white 30deg-radius animal FoV circle.
- matched_profiles.csv: pixel counts and all three bin means.
- comparison.json: FWHM/HWHM and unresolved crossings; peak amplitudes;
  equal-bin normalized-profile RMSE; unscaled mean/median DoLP in the animal
  FoV; measured radiances, source anchors and input hashes.

These are descriptive comparisons. The width interpolation does not resolve
sub-bin confidence intervals. Differences in support are not extrapolated.
Pixel counts are not independent replicates; no pixel-based significance test
is reported. Camera-channel DoLP is not the green510 spider-weighted signal.

## What decides the next step

Check whether the candidate improves width without displacing the band or
requiring site-specific shape adjustment. Amplitude disagreement is kept
separate. The observations are not used to set natural eta or ALAN.
If supported across sites, explicitly define this natural DoLP law and consistent
Stokes amplitudes in a reviewable model variant, then propagate through existing
HPS/LED spectral photon catches, fixed animal field, overnight Moon track and
360deg horizontal headings. Do not silently replace the baseline until that
comparison is reviewed.

## Verification

Three unit tests check the Rayleigh limit, unchanged phase/zeros and shared
pixel selections. A synthetic three-site end-to-end run exercises the maps,
profiles, widths and JSON. These checks validate software, not the real-sky
hypothesis. The EcoVis run is needed for actual measured results.
