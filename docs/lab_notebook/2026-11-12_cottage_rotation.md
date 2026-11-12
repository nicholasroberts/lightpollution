# Cottage rotation and PRC reconstruction — 2026-11-12

Investigator accepted approximate Cottage geometry: movement and individual anatomical variation do not justify pixel-perfect geometry refinement. Use each exposure's middle (5%) fitted centre; reference image 1 radius r60 = rim/sqrt(2). Do not refine geometry to match recorded PRC.

Added src/analyse_spider_sky_rotation.py and tests/test_spider_sky_rotation.py. Original TIFFs on EcoVis have NOT been processed in this workspace. All 28 tests in the current local suite pass. Four new tests cover FFT lag/sign, known coordinate warp, image-based roll recovery and single-quartet end-to-end Stokes/PRC export. Synthetic known Q/I=.2 recovered within .01; both processing previews produced.

## Algorithm

Read numbered TIFFs and verify SHA256 against geometry report; retain blue channel index 2. Use accepted individual centres, without radial scaling. Estimate global scene roll relative to image 1 from polar-sampled local spatial contrast, searching all 360 degrees and refining peak position. Report competing peak, correlation and flags. A scene roll is distinct from eye orientation; warp source coordinates into reference using independently fitted centres and measured roll. No photometric normalization applied to scientific arrays (only diagnostic display scales). Show reference edges red and aligned edges cyan; coincident structure appears white. Residual tilt/parallax/distortion is not corrected.

Nominal four-angle order retained, interleaved groups [1,4,7,10], [2,5,8,11], [3,6,9,12]; four-file sets use [1,2,3,4]. Scene roll sign is opposite physical camera roll. Under nominal 45-degree spacing, reversing the analyzer sign changes U sign and not the maximized global eye contrast; geographic AoP is not established. Measured scene rolls are NOT silently substituted for analyzer calibration. If substantial deviations from nominal spacing appear, review before accepting any PRC.

This implements the historical common image-plane axis approximation, not a full ray-dependent polarimetric lens model. Do not infer exact optical polarimeter calibration from image registration alone. TIFF linearity, exposure matching, export processing and saturation remain possible reproduction limitations.

Two separate outputs: linear_stokes and legacy_dop_aop. Each receives an explicitly UNREVIEWED_PREVIEW PRC result using the existing eye-response code (absorbances .4/.04, natural log ratio after integration). Comparison site ID is s14 (Cottage.csv, rounded PRC .439); reference table notes location alias is not independently confirmed. No model parameter tuned to this number. No historical fit/data changes.

rotation_analysis.json includes analyzer-pair closure, both intensity-retaining and normalized-direction PRC diagnostics, invalid pixel exclusion fractions, and nominal radius factors .95/1/1.05. Their spread is sensitivity, not a confidence interval. Canonical plotted previews retain intensity; alternative processing remains explicit pending identification of the exact historical reduction.

## EcoVis command

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/analyse_spider_sky_rotation.py \
  --geometry data/processed/spider_prc_cottage_geometry_01/geometry_calibration.json \
  --folder data/raw/skyimages/Cottage/OneDrive_10_12-11-2026_extracted \
  --site-id s14 \
  --output data/processed/spider_prc_cottage_rotation_01
```

Output must be new and separate from raw inputs. Raw images remain unchanged. Upload rotation_overlay.png, rotation_analysis.json, and result.json from each of linear_stokes_preview and legacy_dop_aop_preview. These are derived previews and require visual review before scientific acceptance. Do not upload the much larger NPZ arrays at this stage.
