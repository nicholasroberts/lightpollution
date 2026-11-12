# Moon capture-time anchor — 2026-11-12

Investigator confirmed camera times are local: Oxford BST/Europe-London, Cottage Swedish summer time/Europe-Stockholm. Added anchor_moon_from_local_time.py. Prefers original-exposure EXIF, otherwise explicitly flags TIFF DateTime fallback. Uses zoneinfo for date-aware UTC conversion and existing Astropy built-in lunar ephemeris. BST and Stockholm summer conversions tested successfully (UTC+1/UTC+2). Actual ephemeris/real TIFF execution remains on EcoVis.

Predicted Moon elevation supplies a radius/ring on level equisolid image; compass azimuth does not alone provide pixel angle when camera north unknown. Script finds a provisional brightness candidate within +/-5deg of predicted elevation and shows it, not claiming independent Moon identification. Flags timestamps differing from reference field date, metadata fallback, below-horizon Moon. Original field dates not silently rewritten. Explicit --cottage-local-time/--oxford-local-time allow verified corrections later.

Outputs moon_ephemeris_anchors.png and JSON. No model curve, sky data, historical PRC or camera orientation changed. Do not run old automatic-source comparison as validated; first inspect elevation rings/candidates and timestamp flags.

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python -m pip install -r requirements-wang.txt
python src/anchor_moon_from_local_time.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --output data/processed/moon_cottage_oxford_anchors_01
```

Upload both outputs. If date flags appear, verify actual capture date/time rather than assuming export metadata is exposure time.
