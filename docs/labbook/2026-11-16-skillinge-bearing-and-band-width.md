# Skillinge source correction and angular band width

Investigator identified the previously selected lower-left Skillinge bright
source as a streetlight, with the Moon 180deg opposite in image azimuth. The
assistant's prior source identification was wrong; Skillinge results from that
bearing must not be interpreted as Moon-referenced validation.

The override now rotates the vector from the SAVED CAMERA CENTRE to that
rejected streetlight by 180deg. It does not rotate the whole image or reverse
Moon elevation. Source radius still comes from the original 12.59255deg
Moon ephemeris elevation. This is an azimuth reversal at fixed elevation,
NOT a spherical antipode. Slabhouse/Cottage/Oxford source settings are unchanged.
The configuration records the investigator correction and retains the original
streetlight pixel estimate as the input to that rotation. Source review panels
must still be inspected for camera-geometry uncertainty.

## Run

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/run_spider_dark_site_validation.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --moon-anchors \
    data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
    data/processed/moon_skillinge_slabhouse_review_01/moon_ephemeris_anchors.json \
  --fixed-offset-audit data/processed/spider_aop_dark_audit_01/convention_and_atmosphere_audit.json \
  --output data/processed/spider_dark_site_validation_02
```

Outputs in that one directory include the existing review and absolute DoLP
plots plus:

- three_site_normalized_profiles.png: each measured/model profile divided by
  its own peak; horizontal half-maximum line.
- profile_widths.json: peak, half-maximum, interpolated left/right crossings,
  left/right HWHM and FWHM for each measured/model profile. Widths are measured
  only when crossings are resolved in sampled support. Missing crossings are
  explicit nulls, not extrapolated. Linear interpolation of 10deg binned means
  is descriptive; there are no confidence intervals or sub-bin precision claims.
- summary.json includes those width diagnostics and revised frozen-offset AoP
  results. Historical site radiances are retained.

FWHM is measured relative to zero DoLP, without baseline subtraction, smoothing,
peak fitting or width fitting. The model sometimes remains above half maximum
at the end of the sampled angular support; its full width is then not resolved.
Half-widths use the first maximum when a sampled peak is flat.

## Scientific interpretation

Cottage and Slabhouse appear narrower than the prescribed natural envelope.
That deserves quantification after correcting Skillinge geometry. Multiplication
by a constant atmospheric factor eta cannot change normalized width. A constant
ALAN/natural ratio also only rescales the curve. In contrast, a spatially uniform
absolute ALAN intensity can change width if natural sky intensity varies with
viewing direction, because the dilution ratio then varies across the sky.
Angular atmospheric effects, natural-pattern assumptions and camera calibration
may also contribute. These plots do not distinguish those mechanisms.

No angular exponent is fitted and no physical model is changed just to reproduce
width. First inspect corrected source geometry and normalized profiles; then
compare the existing intensity-dependent spectral ALAN forward model and, if
needed, investigate angular atmospheric depolarization. Retain measured camera
versus spider spectral-response differences and sky-mask limitations.

## Software checks

Six runner tests include a 180deg camera-centred bearing reversal, preservation
of elevation/UTC/hash under bearing edits, missing/conflicting anchors, scale
invariance of FWHM, and unresolved crossings. The four-site synthetic workflow
also runs end to end with width outputs. Real results require the EcoVis rerun.
