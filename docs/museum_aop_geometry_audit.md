# Museum AoP geometry diagnostic

This diagnostic uses the checksum-locked linear RAW I/Q/U reconstruction, fixed astronomical Moon coordinates and the previously saved Moon-derived camera-to-north rotation. It never fits a sky-pattern rotation, analyzer zero, source intensity, reflection or ALAN adjustment. The Moon-based north rotation remains provisional.

The unknown constant transmission-axis zero of the rear polarizer prevents an absolute AoP color comparison. The diagnostic therefore measures **AoP changes between pairs of points**, using exp(2i AoP) to respect the 180-degree ambiguity. A constant analyzer offset cancels algebraically. Predicted and observed changes are compared at fixed separations of 16, 32 and 64 original-image pixels by default. I/Q/U are block averaged before calculating angles; AoP itself is not arithmetically averaged. Both horizontal and vertical pairs enter the summary. The picture shows horizontal pairs at the requested primary separation.

Validation metrics use the fixed central 60-degree full-field circle, positive finite intensity, DoLP from 0.05 to 0.95, and exclusions within 10 degrees of each source/antisource singularity. Low polarization makes AoP unstable; these thresholds are diagnostic choices, not biological thresholds. A reversed-U case and a constant-AoP control are reported independently. No best case is automatically selected or applied as a correction.

The primary prediction uses fixed astronomical Moon elevation/azimuth with the saved azimuth rotation and a provisional zenith-axis/equisolid assumption. The geometric separation between its projected Moon and the observed image Moon is recorded, not removed by fitting. A second predictor uses the locked observed image direction as a separately labelled sensitivity check. Ideal rear-analyzer ray transport is used; no empirical fisheye Mueller matrix is claimed.

## Solar control

At the verified museum coordinates (51.758527 N, 1.256267 W), the saved capture UTC is 2019-05-14 21:09:51. An offline Astropy calculation gives geometric Sun elevation **−10.2274 degrees** and azimuth **319.0842 degrees**. The same calculation recovers Moon elevation 41.1983 degrees and azimuth 187.1466 degrees, consistent with the existing anchor. Thus a solar/twilight contribution is a relevant hypothesis to test. Nautical twilight ends at Sun elevation −12 degrees according to the US Naval Observatory: https://aa.usno.navy.mil/faq/RST_defs

The script computes a fixed-Sun scattering-plane control from the verified site coordinates and locked UTC. This is only an ideal angular control: no twilight brightness, atmospheric multiple scattering or Sun/Moon mixture is fitted. Agreement would support further investigation, not prove a unique physical cause. The original lunar model and its nightly curves remain unchanged.

## Run on EcoVis

```bash
git pull origin publication-figures

python src/audit_spider_lunar_aop.py \
  --reference-lock data/processed/spider_registration_museum_linear_01/locked_site_references.json \
  --site museum \
  --output data/processed/museum_aop_geometry_audit_01
```

Use a fresh output directory. Outputs:

- `aop_geometry_audit.png`: measured DoLP/AoP, fixed-Moon prediction, measured/predicted AoP changes and their difference.
- `aop_moon_sun_control.png`: independent Moon/Sun/constant geometry comparisons.
- `aop_geometry_summary.json`: masks, pair counts, median/p95 circular change differences, solar position, input checksums and explicit limitations.

Three tests check exact analyzer-zero cancellation across angle wraps, sensitivity to U reversal, the zenith tangent basis, block averaging and an end-to-end synthetic checksum-locked case with an unknown global analyzer zero and independently computed solar position. Real museum AoP conclusions require running this diagnostic beside the new RAW Stokes file on EcoVis; the supplied DoLP images alone do not contain AoP.
