# Museum registration audit — 16 November 2026

The whole camera was rotated with the polarizer fixed, with three exposures at each nominal angle. The supplied Sigma 8 mm f/3.5 EX DG images have the same dimensions (2000 × 1333). First and last original capture times are 14 May 2019 22:09:51 and 22:14:10 local. Input 1 is the fixed reference; originals are unchanged.

## Implemented refinement

`--registration-model fisheye-shared` adds a bounded shared tangential image-coordinate correction to the existing equisolid SO(3) registration. Both correction coefficients are common to all exposures. Each exposure retains its own attitude and bounded projection-centre offset. The shared correction is estimated from exposures 4 and 7 only; the 135-degree group checks transfer to another angle. No theoretical polarization, Stokes, DoLP, AoP, PRC, ALAN, or astronomical coordinates enter the geometry fit.

Entire 96-pixel reference cells are withheld before both shared and per-exposure fitting. Holdout membership is `(floor(x/96)+2*floor(y/96)) % 4 == 0`. Initializers also exclude held-out points. The old numerical thresholds are retained: median held-out error <2 px, at least 60% within 3 px, minimum match counts, broad spatial coverage, and no optimization boundary hit. The correction is checked for positive Jacobian and inverse roundtrip accuracy throughout the supported field. Final scientific interpolation samples each original exposure once.

This is an empirical registration refinement, not an independent angular or polarimetric lens calibration. It does not demonstrate absence of all residual ghosting. Site status remains `numerically_locked_awaiting_visual_review` after passing checks. Vegetation, stars, saturated edges and nearby structures can still disagree. Nominal analyzer equations are retained, so DoLP calibration is a separate question.

## Full raw-set result

Reference 1 is identity. Other frames all pass the unchanged numerical criteria:

| Exposure | Median held-out error (px) | Held-out landmarks within 3 px |
|---|---:|---:|
| 2 | 0.146 | 99.1% |
| 3 | 0.168 | 99.7% |
| 4 | 0.864 | 96.9% |
| 5 | 0.887 | 96.7% |
| 6 | 0.897 | 98.0% |
| 7 | 1.330 | 92.0% |
| 8 | 1.319 | 91.0% |
| 9 | 1.349 | 91.6% |
| 10 | 1.502 | 90.1% |
| 11 | 1.554 | 90.3% |
| 12 | 1.586 | 90.3% |

Shared coefficients (dimensionless, normalized to the saved rim radius): p1 = 0.0023155306705066876, p2 = 0.0011901354994268156. Minimum correction Jacobian determinant = 0.97896. Reference centre/radius remain the inherited provisional geometry; they were not calibrated against the Moon or polarization pattern.

End-to-end validation rebuilt `linear_stokes.npz`, froze per-exposure mappings and produced a scalar north-up DoLP preview. The image Moon candidate was independently localized at original pixel (654, 610); the fixed astronomical azimuth supplies a provisional display rotation. Moon identity, image handedness and residual geometry still require review; camera-basis Q/U are not relabelled celestial Stokes.

Tests: existing rigid/fisheye tests plus shared-correction inverse, synthetic transfer to an untrained angle, no-reflection and spatial-holdout checks pass (10 tests). Full 12-TIFF integration and CLI/source-lock replay also pass. Museum is the tested real site; no all-site acceptance claim is made.

## Reproduce on EcoVis

```bash
git pull origin publication-figures

python src/lock_spider_site_registration.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --model-run data/processed/spider_heading_all_sites_01/all_sites \
  --sites museum \
  --registration-model fisheye-shared \
  --output data/processed/spider_registration_shared_museum_01

python src/plot_spider_site_details.py \
  --run data/processed/spider_heading_all_sites_01/all_sites \
  --batch data/processed/spider_prc_all_sites_01 \
  --reference-lock data/processed/spider_registration_shared_museum_01/locked_site_references.json \
  --sites museum --alan-level 5 --atmospheric-factor 0.8 \
  --output data/processed/spider_site_details_shared_museum_01
```

Both outputs must be fresh directories. Review `museum/registration_edges.png` and `museum/registered_dop_and_north.png`. The reporting subset flag does not alter the model, 1-degree headings, spectra, or physical scenario.
