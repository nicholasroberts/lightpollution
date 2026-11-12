# All-site band versus zenith comparison — 2026-11-12

Reviewed Gavarnie view_direction_summary.json: maximum-mean-DoLP candidate at elevation 45deg, image-relative azimuth 120deg. Direct intensity PRC .602928 (band) versus .232896 (zenith), legacy .603623/.234153. Historical .470 lies between them. Viewing direction therefore materially changes the inferred contrast, but band maximum does not exactly reproduce historical value. No inference of exact historical field direction justified.

Investigator requested all sites and site-wise plot of viewing positions versus historical values. Interpreted plot as PRC on y-axis, sites on x-axis, with zenith-centred, band-centred and historical markers. Added batch_spider_view_directions.py. Reuses existing combined arrays; no registration. Explicit 18-site roster checked against reference CSV IDs. For each site run same 30-radius spherical-cap selection and azimuth scan. Selection maximizes mean DoLP, never historical agreement or PRC. Inherited alignment/crop issues and possible terrain contamination require individual overlay review.

Use existing corrected Gavarnie registration override and reuse its completed direction-test result, verifying report and Stokes hashes. Other sites start from all_sites_01 registration directories. Output incremental all_sites_view_directions.csv and view_direction_batch_status.json. Plot sites_zenith_band_historical.png/.svg/.pdf has two panels for direct and legacy processing, each containing zenith/band/historical values; intensity-retaining mode plotted consistently. CSV also preserves normalized-direction alternatives and selected elevation/image azimuth. No method mixing by site.

Two targeted tests pass for known-site overrides and preservation of original/band differences; plot-export smoke check confirms PNG/SVG/PDF produced. Full real-site comparison still requires EcoVis run. --resume checks completed report/array/script identity; failed partial directories are preserved and require fresh output.

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/batch_spider_view_directions.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --site-run Gavarnie=data/processed/spider_prc_gavarnie_alignment_02 \
  --reuse-site Gavarnie=data/processed/spider_prc_gavarnie_band_direction_01 \
  --output data/processed/spider_prc_all_sites_band_01
```

Upload sites_zenith_band_historical.png, all_sites_view_directions.csv and view_direction_batch_status.json. New results do not replace historical reference values, fitted functions or locked products.
