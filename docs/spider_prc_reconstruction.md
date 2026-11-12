# Drassodes sky-image reconstruction

## Current field definition — corrected 12 November 2026

Each PM eye has a **60-degree full field: 30-degree angular radius**. This supersedes earlier 60-degree-radius instructions and manifests. The field is a circular, equally sampled pixel region, with no extra angular weighting. Integrate responses before taking ln(Rmax/Rmin); the axes of the two eyes are globally orthogonal. Use absorbances .4/.04 from supplied code. Retain intensity as the primary explicit reconstruction; normalized-direction results remain a separately reported alternative pending historical provenance.

The core manifest requires `fov.angular_radius_deg: 30` and `fov.radius_px`. Old `radius_60_px` manifests are rejected, not silently reinterpreted. In rotation_analysis.json, `radius_60_px` is retained only as a legacy geometry reference for existing audit scripts; `eye_radius_px`, `eye_angular_radius_deg: 30` and `eye_full_field_deg: 60` identify the actual integration field. Equisolid mapping gives r30/r90=sin(15deg)/sin(45deg).

## All-site run

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/batch_spider_prc_sites.py \
  --root data/raw/skyimages \
  --geometry-reference data/processed/spider_prc_cottage_geometry_01/geometry_calibration.json \
  --output data/processed/spider_prc_all_sites_01
```

The explicit roster contains the 18 historical sites only. Known aliases ChurcHill/ChurchHill, skillinge/Skillinge and WC_pitch/WCpitch are accepted. Bristol, Sandhammeren and Wells are excluded from historical comparison. Waitrose accepts the investigator-confirmed single 1..4 quartet (0/45/90/135). Other complete 1..12 sets retain interleaved triples [1,4,7,10], [2,5,8,11], [3,6,9,12]. Raw filenames/data are not changed.

For Cottage, verify and reuse the accepted per-exposure geometry. Other sites receive explicitly approximate geometry transferred from Cottage image 1, scaling radius with the shorter image dimension and centre by image fractions. This is a practical starting approximation, not independent site calibration. Different aspect ratios/crops receive flags; all transferred fields require visual review. Do not infer reliable geometry from an obstructed mountain skyline. Rotation/translation are then fitted per exposure; weak or non-nominal matches are flagged and retained.

Outputs:

- `all_sites_prc.csv`: site, reference PRC, both methods and response modes, differences, replicate counts, invalid fractions and alignment/geometry flags.
- `batch_status.json`: all site statuses, including failures/missing inputs. Updated after each site.
- `<site>/registration/`: rotations, edge overlay, derived arrays and explicitly unreviewed PRC previews.
- `<site>/polarization/`: matched plots and corrected-field numerical summary.
- `<site>/*.log`: per-stage logs.

Run sequentially to bound memory usage. Failures are reported and the batch continues; final exit status is nonzero if any site failed. `--resume` reuses completed sites only when input, geometry and code hashes match. Partial/stale site outputs are preserved and reported; use a fresh batch output for rerunning those rather than overwriting.

Upload the top-level CSV and status JSON first. We can then select flagged or discrepant sites for image review without uploading the large NPZ arrays. All new site outputs remain provisional until reviewed. Existing fitted PRC/radiance curves and historical reference observations are untouched.

## Existing Cottage arrays

No re-registration required to audit the corrected field:

```bash
python src/audit_spider_polarization.py \
  --run data/processed/spider_prc_cottage_rotation_02 \
  --angular-radius-deg 30 \
  --output data/processed/spider_prc_cottage_fov60_full_01
```

Cottage corrected-field results: legacy intensity-retaining PRC .446818; direct Stokes .454088; reference .439. Near reproduction, not an exact numerical match. Original whole-sky AoP figure provenance remains unresolved. Both DoP/AoP reconstruction methods stay explicit; do not choose parameters to minimize the difference from .439.

## Method limits and verification

Camera and rear polarizer rotate together. Scene image registration is distinct from global eye-axis optimization. Nominal analyzer angles remain 0/45/90/135; measured registration rotations are not silently substituted for a physical analyzer calibration. Global common image-plane axes follow the historical reconstruction model; no full ray-dependent lens Mueller model is implied. Residual tilt, parallax, crop geometry, TIFF radiometric linearity and saturation can affect results. Legacy arithmetic AoP averaging can introduce wrap seams; direct Stokes averaging is retained separately.

Local test suite: `python -m unittest discover -s tests`. As of the all-site addition, 35 tests pass, including coordinate transforms, synthetic polarization recovery, cancellation, wrap handling, unique 18-site roster, Waitrose layout and geometry scaling/flags. Original all-site TIFFs are processed only on EcoVis.
