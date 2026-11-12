# Drassodes PM-eye contrast reconstruction

## Scope and current validation

`src/reconstruct_spider_eye_contrast.py` implements the investigator's
12 November specification in Python: two orthogonal PM-eye polarization axes,
a common circular fisheye region extending to angular radius **60 degrees**,
equal contributions from all included pixels, integrated eye responses, and
`PRC = ln(Rmax/Rmin)`. No extra angular, solid-angle or anatomical weights are
introduced. The source image is not manually rotated to maximize eye contrast.
Original raw pixels and files are never rewritten.

This is a new reconstruction from the source equations and investigator's
description, not a recovered verbatim spider program. Synthetic tests have
passed locally. Original images are on **EcoVis**, not available in the coding
workspace: no real image folder has yet been identified or analysed by this
session. Original-site numerical validation is pending. No radiance dataset,
locked spider transfer function, moth result or European map is changed.

## Equations and conventions

The decimal absorbances from the recovered source are 0.4 and 0.04. Define

    r = 10**(-0.4)
    t = 10**(-0.04)
    B = 1 - (r+t)/2
    A = (t-r)/2
    Q_theta = Q*cos(2*theta) + U*sin(2*theta)
    response(theta) = B*I - A*Q_theta

At theta=0 and 90 degrees these match the original R `par`/`perp` expressions.
Here I,Q,U refer to a common analyzer/image-plane reference frame. There is no
spectral integration or implicit conversion to a newly assumed photopigment.

Sum over the included circular FOV BEFORE optimizing orientation:

    M = B * sum(I)
    D = A * hypot(sum(Q), sum(U))
    Rmax = M + D
    Rmin = M - D
    PRC = ln((M+D)/(M-D))

The maximum axis is `0.5*atan2(-sum(U),-sum(Q))` modulo 180 degrees; the other
is 90 degrees away. For zero amplitude the preferred axis is undefined and
reported as null, not a measured heading. The orientation curve is also
exported explicitly. This is a single global eye orientation, never an
independent maximum for every pixel. Oppositely polarized sky patches can
cancel in the integrated response. Ratios of sums and averages agree because
both eyes use identical pixel support. Only the maximum nonnegative log ratio
is reported as PRC; the signed heading curve remains available.

This analytic maximization removes the need to manually align a complete sky
image's polarization pattern to the model eye. It does **not** fix spatial
misregistration of analyzer exposures, arbitrary camera roll conventions,
fisheye geometry, or per-direction polarization reference-frame errors.
The circular FOV is common and centred by the manifest; there is no additional
inferred eye-axis transport over a curved sky.

### Intensity and normalized-direction alternatives

`response_mode: stokes_intensity` sums the actual I/Q/U-derived responses with
equal pixel weights, consistent with the recovered image-response equations.
Spatial intensity variation is part of each response, not an added weighting.
The alternative `normalized_direction` divides Q,U by I and sets I=1 at each
valid pixel before integration. This is a different, explicitly labelled
calculation if historical work normalized sky directions. The investigator's
no-extra-weighting instruction does not require this normalization; default
intensity mode is retained. No choice is made to improve agreement with sigg.

## 1. Locate the EcoVis folders

From the existing repository/venv:

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
source .venv/bin/activate
python -m pip install -r requirements-prc.txt
python src/reconstruct_spider_eye_contrast.py discover \
  --root /run/media/bznwr/data \
  --root "$HOME/Downloads" \
  --output data/processed/spider_prc_discovery_01
```

`discover` uses filenames only; it does not open huge TIFF pixel arrays. It
walks explicitly supplied roots, skips hidden directories and symlinked
subdirectories, deduplicates overlapping roots, reports permission errors and
missing roots, and stops at a directory-count limit instead of silently
truncating. It does not extract archives or guess analyzer angles from arbitrary
filenames. Add another `--root` if the downloaded data are elsewhere on EcoVis.

It writes `inventory.json` with complete/incomplete/duplicate numbered TIFF
sets, other image folders (including unfamiliar TIFF/JPEG/NEF/DNG/NPZ names),
and archive paths not yet extracted. Complete 1–12 folders get candidate JSON
manifests. **Folder names are not assigned to sites automatically.** Inspect
and select one case first. A missing image or duplicate `1.tif`/`01.tiff`
prevents automatic candidate-manifest creation for that folder.

Discovery reports any missing search roots in its stdout and inventory. Zero
candidates is an inventory result, not permission to manufacture a dataset.
Use a fresh output directory name when rerunning: all commands refuse overwrite.

## 2. Configure one selected case

Use the generated candidate manifest, or copy
`config/spider_prc_manifest.example.json` to an ignored location under
`data/processed/`. Keep machine-specific manifests out of Git.

- `input.directory`: actual folder. Relative paths resolve against manifest
  location; input filenames resolve against the input directory.
- `channel`: zero-based colour channel; 2 reproduces R's third array channel.
  Grayscale TIFFs require null. No spectral sensitivity is inferred from this
  selection. TIFFs must already be suitable linear-intensity products; this
  script does not demosaic NEFs, reverse gamma, estimate darks or exposure-match
  images. Establish their processing/exposure identity before scientific use.
- `groups`: one to three explicit replicates at 0/45/90/135 degrees; recovered three-replicate ordering
  is (1,4,7,10), (2,5,8,11), (3,6,9,12). Other filenames can be supplied
  explicitly. No filesystem order is used to establish angle or replicate.
  Waitrose's single replicate is now investigator-confirmed as files
  1,2,3,4 at 0,45,90,135 degrees respectively. Use [["1.tif","2.tif","3.tif","4.tif"]]
  for that site, never repeat those files to manufacture additional replicates.
- `rotation_deg_clockwise`, `shift_xy_px`: rotate about the image centre,
  then translate in displayed x-right/y-down coordinates. Starting -45/-90/-135
  values come from the R source but R/OpenCV/SciPy sign conventions are not
  assumed equivalent. **Inspect their effect on actual images**. If exports
  are already co-registered, use zero transforms instead of rotating twice.
  These fixed transforms are reused across replicates, as in the source.
  Automatic optimization against the target PRC is not performed.
- `fov.center_xy_px`: measured [x,y] centre, zero-based pixel coordinates.
- `fov.radius_60_px`: independently calibrated radius for **60 degrees from
  the centre**, not a 60-degree full diameter. The full circle must fit.
  No image-edge-to-90-degree projection or lens model is guessed. If equidistant
  mapping is independently established, r60=(60/90)*r90; other projections differ.
- `site_id`: explicitly verified s01–s18, or null for an unidentified case.
  Matched values are compared to the existing rounded 3-decimal reference,
  including Slabhouse as a measurement, without refitting or choosing geometry
  to minimize residuals. A rounding match is not validation of image identity.
- `invalid_policy`: `error` by default. If a reviewed mask is necessary,
  `exclude` uses identical support for both eyes and records excluded counts
  and fraction. No zero filling of missing sky or replacement of NaNs with live
  measurements is introduced. Saturation fractions are reported per input;
  saturated values are not automatically corrected or silently excluded.

### Explicit Stokes methods for TIFFs

Select `input.stokes_method`; discovery leaves it unselected:

1. `legacy_dop_aop`: per replicate q=(I0-I90)/(I0+I90),
   u=(I45-I135)/(I45+I135); DoP=hypot(q,u), AoP=atan2(u,q)/2. Preserve the
   supplied R's DoP NaN/>0.95 zeroing and arithmetic AoP averaging across
   replicates. Use average I0+I90 and reconstruct Q=I*DoP*cos(2AoP),
   U=I*DoP*sin(2AoP). This reproduces that source's DoP/AoP intermediates,
   **not proof that these specific intermediates fed the original spider model**.
2. `linear_stokes`: average I=I0+I90, Q=I0-I90, U=I45-I135 across aligned
   replicates. No legacy DoP threshold or arithmetic angle averaging. This is
   an explicitly different method and is never silently substituted for legacy.

Neither method guesses how historical intermediate images were normalized.
Per-replicate analyzer-sum closure `(I0+I90-I45-I135)/(I0+I90+I45+I135)`
provides a diagnostic for intensity mismatch and misalignment. Median/p95
absolute closure within the FOV is recorded, but no invented acceptance
threshold is imposed. Texture, moving clouds, exposure, analyzer throughput,
and registration can all affect it; good closure alone does not prove alignment.

Both methods use derived float64 arrays without changing raw integer pixels
or bit depth. Rotation resampling is bilinear with NaN unsupported borders,
not a new nearest-neighbour or rescaled raw observation. SciPy registration
is not claimed numerically identical to the original RNiftyReg process.

### Saved map inputs

TIFF inputs may use either `numbered_tiffs` (existing manifests) or the alias
`tiff_groups`. Both now accept one to three explicitly ordered four-angle groups.
Discovery also creates unreviewed single-quartet candidates for folders with
exactly files 1–4; this is a layout candidate, not proof that other incomplete
12-image folders are genuine one-replicate observations.

For already processed maps, set `input.kind` to:

- `stokes_npz` with `path` pointing to an NPZ containing 2D arrays I,Q,U; or
- `dop_aop_npz` with arrays DoP,AoP,I and `aop_unit` explicitly degrees/radians.
  I may be omitted only in explicit normalized_direction mode.

All map dimensions must agree. Nonfinite, nonpositive-I and physically invalid
DoLP>1 pixels inside the FOV cause an error unless exclusion is explicitly chosen.

## 3. Inspect a preview, then run

After setting the measured geometry and selecting a processing method:

```bash
python src/reconstruct_spider_eye_contrast.py run \
  --manifest data/processed/spider_prc_discovery_01/candidate_0000.json \
  --preview-only \
  --output data/processed/spider_prc_preview_01
```

The preview is labelled UNREVIEWED_PREVIEW. Check correct sky selection,
geometry, channel/exposure treatment, alignment, analyzer-sum closure, missing
pixels, saturation and processing conventions. Set geometry_reviewed and
registration_reviewed true in the local manifest when those have actually
been checked, then run to a new output directory without --preview-only.
Setting these fields does not make the new algorithm historically validated.

Each output folder includes:

- `result.json`: PRC, integrated responses, maximizing axes, support and
  optional signed difference from the recorded site value.
- `provenance.json`: configuration, input SHA256 hashes/dtypes, code/config
  hashes, Git commit, dependency versions and quality diagnostics.
- `derived_stokes.npz`: separate float64 I,Q,U plus the circular FOV mask.
- `eye_orientation_curve.csv`: signed log contrast vs axis angle, 0–179.5°.
- `qc_preview.png`: intensity, DoLP, AoP and closure images with the FOV outline.

The plots downsample for display only; calculations use all input pixels.
They are initial image QC, not a certificate of successful original registration.
Large input arrays are processed one four-image replicate at a time; the
derived maps and three closure maps are retained. Memory scales with frame size.

## Verification

```bash
python -m unittest discover -s tests -p test_reconstruct_spider_eye_contrast.py -v
```

Tests cover analytic values, original fixed-axis R equations, rotation
invariance, integrate-before-maximize cancellation, unpolarized-light dilution,
sum/mean equivalence, explicitly different normalization methods, bad pixels,
registration sign/borders, legacy processing conventions, folder discovery,
16-bit TIFF loading, unchanged raw file hashes, and overwrite protection.
Use the EcoVis inventory and a selected real case to establish actual input
coverage and compare against the original value before using spectral scenarios.

## TIFF header inspection before choosing the image geometry

The uploaded EcoVis inventory contains 20 numbered 12-image folders and one
four-image Waitrose folder: 244 TIFF filenames, 43 directories scanned,
no duplicate numbered images and no discovery warnings. It is a filename
audit, not a pixel-integrity or scientific-alignment check. The original
`archives_not_extracted` field merely listed ZIP files encountered, even where
their extracted folders already existed. Version 2 calls this field
`archive_files_found`; no new extraction is requested.

The investigator has since renamed Waitrose's DSC_5080–DSC_5083 TIFFs to 1–4
and explicitly confirmed their angle order. The old inventory retains those
historical names, but `inspect` lists current TIFF filenames in that inventory's
other-image folders, so it will inspect the renamed Waitrose files correctly.

To check headers across the inventory and produce an unaligned Gavarnie input
quicklook without guessing a 60-degree pixel radius:

```bash
python src/reconstruct_spider_eye_contrast.py inspect \
  --inventory data/processed/spider_prc_filename_audit_01/inventory.json \
  --preview-folder data/raw/skyimages/gavarnie/OneDrive_12_12-11-2026_extracted \
  --output data/processed/spider_prc_tiff_inspection_01
```

Outputs are `tiff_inspection.json` (page-0 shapes/dtypes, available tags,
compression, page count, file sizes, read errors) and `input_quicklook.png`.
The four raw input previews use the first numbered quartet 1,4,7,10 when
available; filenames remain visible and no exposure registration is performed.
They use individual display scales and the third colour channel (or grayscale
for 2D TIFFs), solely for viewing, not normalization of scientific responses.
All headers are inspected, but not all image pixel data are decoded/tested.
The header report is saved before preview creation so unsupported codecs do
not erase the metadata audit. No exposure tags are assumed to exist.

If site folders are renamed, rerun discovery with a fresh output directory
before inspection because the old inventory contains their former paths.
