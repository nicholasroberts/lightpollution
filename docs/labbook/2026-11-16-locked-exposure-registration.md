# Register exposures before locking image coordinates and true-north displays

Ghosted skyline/building structure in museum and Ounsdale means the old
reconstructed experimental intensity and DoLP cannot serve as reviewed alignment
references. Linear/log displays use the same composite. They do not cause it.
The old solver refined mainly inner-sky structure, which could leave outer-field
landmarks misregistered. This is a plausible mechanism, not an established sole
cause at the real sites until their individual exposures are processed.

## Run on the actual raw exposures

```bash
git pull origin publication-figures
source .venv/bin/activate
python -m pip install -r requirements-registration.txt

python src/lock_spider_site_registration.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --model-run data/processed/spider_heading_all_sites_01/all_sites \
  --output data/processed/spider_registration_locked_01
```

The default processes every site with a saved registration report (including
sites missing astronomical timing). To inspect only the clearest failures first,
add `--sites museum Ounsdale` and use a different fresh output directory. Original
TIFFs and prior results are never overwritten. Missing or rejected sites receive
explicit pending entries, never replacement zero DoLP/PRC.

Then generate the same per-site ALAN=5 reports from the frozen products:

```bash
python src/plot_spider_site_details.py \
  --run data/processed/spider_heading_all_sites_01/all_sites \
  --batch data/processed/spider_prc_all_sites_01 \
  --reference-lock data/processed/spider_registration_locked_01/locked_site_references.json \
  --alan-level 5 --atmospheric-factor 0.8 \
  --output data/processed/spider_site_details_locked_01
```

The model's 1 degree PRC curves, actual recorded radiances, spectra, ALAN scales,
atmosphere, Moon azimuth/elevation and UTCs are untouched. Experimental DoLP is
rebuilt; the previously reconstructed experimental PRCs should not be treated as
validated. The report reads new experimental products only for numerically
accepted sites, checks their hashes, and checks that astronomical anchors agree.

## Exposure registration

- Fix raw `1.tif` as the site coordinate reference.
- Verify every TIFF against its original saved SHA256; preserve explicit analyzer
  grouping and channel. No exposure/gamma normalization in scientific output.
- Detect reciprocal SIFT correspondences in the outer field (0.55–1.03 times
  saved rim radius), which includes static trees, buildings and skyline.
- RANSAC fits rotation and translation only: no scale, shear or reflection.
- Withhold every fourth descriptor correspondence before estimating the rigid
  model. At least eight training inliers, four held-out inliers, 60% held-out
  agreement within 3 pixels, median held-out error below 2 pixels, three occupied
  angular sectors and half-rim landmark span are required. These are numerical QC
  criteria, not certification of the biological model or image identity.
- Resample each raw exposure once into image-1 coordinates; unsupported pixels
  remain NaN. Rebuild linear Stokes with the existing nominal analyzer convention.

Inspect **SITE/registration_edges.png**. Reference edges are red, registered
exposure edges cyan, coincident edges white. Inspect buildings and horizon as
well as stars. Parallax, lens distortion, moving vegetation and changing clouds
can prevent a rigid transform from succeeding. Such cases need further review;
the code does not distort a scene to force agreement. The old geometry inherited
from Cottage is frozen and labelled; this change does not independently calibrate
the lens projection or improve the accuracy of a transferred FoV radius.

## Astronomical orientation: rotate experiments, keep theory fixed

The experimentally identified Moon's bearing in raw image 1 sets the image
rotation to the independently calculated compass azimuth. No theoretical curve,
DoLP or AoP is fitted to the images. No global source search or ephemeris-ring
search is used to choose the Moon. Existing audit-bearing pixels seed four sites;
the seven additional clear-source sites use approximate visual seed patches.
Local intensity contrast is searched in a small patch on raw 1.tif. These
identifications remain provisional until the saved raw-reference audit is
inspected. Ambiguous sites stay without a Moon coordinate/north rotation.

Investigators can supply exact identified source pixels with:

```json
{"museum":{"moon_pixel_xy":[630,610],"note":"Example only: replace with identified Moon coordinates in raw 1.tif"}}
```

Pass the actual reviewed file using `--moon-coordinates FILE`. The example numbers
are illustrative, not newly confirmed coordinates. The image-implied elevation
is checked against the fixed ephemeris; discrepancy over 6 degrees leaves north
unresolved. The check does not shift either position or change the camera clock.

North-up export applies a pure, orientation-preserving rotation. The old display
helper's implicit reflection is not inherited. A single Moon fixes bearing only
under the stated level-camera/projection/handedness assumptions; independent
handedness, lens and tilt validation remain separate tasks. North-up products
contain scalar DoLP and intensity only. Camera Q/U remain camera-basis arrays;
they are not relabelled as celestial Stokes without spin-2 calibration.

Inspect **SITE/registered_dop_and_north.png** for the raw reference source, rebuilt
camera DoLP and north-up experimental DoLP where supported.

## Locked reference set

**locked_site_references.json** and **SITE/reference_coords.json** freeze:

- Exact input identities/hashes, reference exposure and pixel convention.
- Reference centre, rim and 60-degree full-FoV radius, with calibration status.
- Per-exposure source-to-reference matrices and source reference-centre pixels.
- All matched landmark coordinates, held-out membership and residuals.
- Moon pixel identity/provisional status, image elevation residual and unchanged
  astronomical anchor.
- Separate camera-to-north rotation matrix, status and display-product paths.

**reference_sha256.json** records output hashes. A change requires a new output
version; downstream plots read these products rather than estimate coordinates
again. Status explicitly remains `numerically_locked_awaiting_visual_review` for
accepted sites: immutable numerical coordinates are not a claim that visual QC
has occurred. Failed/unresolved records remain in the same inventory.

## Validation

Five unit tests check rigid recovery with 30% false correspondences, full SIFT
matching of rotated/translated synthetic texture, no-mirror north rotation, and
Moon-source selection independent of astronomical elevation, and rejection of
changed locked products or astronomical anchors. A four-exposure
synthetic TIFF set exercises the full raw-input, registration, Stokes, north-up
and frozen-manifest path. Real-site alignment still depends on running the above
commands where the original TIFFs reside and inspecting the saved overlays.
