# Museum green-channel replication

`prepare_spider_linear_raw.py --channel green` extracts both green sites from each visible Bayer 2x2 cell, subtracts their own reported channel black levels, then averages the two samples. It retains negative values, avoids white balance/rendering/demosaicing, and uses the existing common area downsampling. Each frame records the green positions, black levels, median and p95 absolute fractional green-pair differences, and saturation fractions. The pair comparison includes adjacent-pixel scene structure and is a diagnostic, not a flat-field calibration. Default `--channel blue` and existing blue manifests remain supported.

Green manifests are explicitly identified and require fresh landmark registration. The registration lock records camera_channel and a channel-specific radiometry status. The combined fitter accepts either linear blue or green locks and labels its outputs. Geometry, Moon localization and edge overlays must be reviewed for the new channel. The existing blue result and source NEFs remain unchanged.

On EcoVis:

```bash
git pull origin publication-figures

python src/prepare_spider_linear_raw.py \
  --raw-dir data/raw/skyimages/museum \
  --channel green --exclude-exposures 7 \
  --output data/processed/museum_linear_green_01

python src/lock_spider_site_registration.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --model-run data/processed/spider_heading_all_sites_01/all_sites \
  --sites museum --registration-model fisheye-shared \
  --linear-raw-manifest data/processed/museum_linear_green_01/linear_raw_manifest.json \
  --output data/processed/spider_registration_museum_green_01

python src/fit_spider_twilight_moon.py \
  --reference-lock data/processed/spider_registration_museum_green_01/locked_site_references.json \
  --site museum --handedness reversed \
  --output data/processed/museum_twilight_moon_green_01
```

Use the same explicit handedness hypothesis as the successful blue comparison, without claiming independent verification. If registration fails, do not substitute another channel's arrays into its lock. Return the green fit PNG/JSON, linear_raw_manifest.json, registration_summary.json and the museum registration_edges.png for comparison. The default 96px spatial-block split and central field rule are unchanged; fitted geometry is re-estimated from green landmarks rather than copied silently.
