# Register the museum linear RAW set

The museum converter produced 2000 × 1332 float32 camera-blue images. Exposure 7 is excluded. The complete reconstruction quartets are `[2,5,8,11]` and `[3,6,9,12]`; reference image 1 and frames 4/10 remain available for geometric registration.

`--linear-raw-manifest` validates conversion completeness, disabled display processing, output checksums and original RAW/TIFF capture timestamps. It reads each prepared TIFF as a single channel (`channel=null`). It requires a single explicit site and prohibits reuse of prior landmarks through `--from-lock`. Old reference centre/radius are transferred only as a provisional canvas-scaled initializer; new landmarks and transformations are fitted. This does not independently calibrate angular projection or establish that the optical axis is zenith.

When frame 7 is absent, shared geometry is trained on frames 4 and 10. Spatial holdout checks remain in force; failure still withholds scientific products. All images are resampled once from the prepared linear inputs. Only complete analyzer quartets enter the Stokes average. Astronomical anchors, theoretical patterns, ALAN and the retained 1-degree model heading grid are unchanged.

```bash
git pull origin publication-figures

python src/lock_spider_site_registration.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --model-run data/processed/spider_heading_all_sites_01/all_sites \
  --sites museum \
  --registration-model fisheye-shared \
  --linear-raw-manifest data/processed/museum_linear_raw_01/linear_raw_manifest.json \
  --output data/processed/spider_registration_museum_linear_01

python src/plot_spider_site_details.py \
  --run data/processed/spider_heading_all_sites_01/all_sites \
  --batch data/processed/spider_prc_all_sites_01 \
  --reference-lock data/processed/spider_registration_museum_linear_01/locked_site_references.json \
  --sites museum --alan-level 5 --atmospheric-factor 0.8 \
  --output data/processed/spider_site_details_museum_linear_01
```

Use fresh output directories. Review the registration overlays and measured Moon localization before interpreting a north-up comparison. Linear RAW extraction is not a complete polarimetric lens calibration. No real-RAW registration result is claimed until this run completes on EcoVis.

Verification: twelve rigid/fisheye/shared-geometry tests pass, including synthetic transfer when frame 7 is absent and an integration fixture with a changed canvas, eleven input files, explicit exclusions, source hashes/timestamps and two complete reconstruction quartets. The fixture checks that reference/registration-only frames do not contribute to the Stokes average.
