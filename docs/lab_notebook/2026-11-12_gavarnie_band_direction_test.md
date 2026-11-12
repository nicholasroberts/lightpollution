# Gavarnie band-centred field versus azimuth — 2026-11-12

Investigator suggested original field was centred on the polarization band rather than zenith, especially at Gavarnie where high Moon changes band placement, and requested band-centred versus azimuth test. This is distinct from changing global polarization-axis orientation.

Added src/test_spider_view_direction.py. Reuses saved combined Stokes arrays and alignment. Inverse nominal equisolid mapping supplies sky unit vectors. A field contains pixels within a 30-degree spherical cap of its centre, preserving equal pixel contributions; off-axis mask is not a translated fixed pixel circle. Zenith mask recovers original 30-radius circle. Test checks this equality and 90-degree azimuth rotation covariance, both passed.

Select candidate band direction using maximum mean physical local DoLP from direct Stokes on predeclared zenith 0..60deg in 5deg increments and image azimuth 0..350deg in 10deg increments. Position-selection scan downsampled to maximum dimension <=500. Require >=95% physical pixels. This is an operational definition requiring visual review, not evidence of historical field direction. No terrain mask available; terrain could contaminate selected field and must be checked on overlay. No PRC/reference number enters position selection.

Then calculate full-resolution 30-degree-radius eye integration for zenith and a 5-degree azimuth sweep at selected band-centre elevation, for both processing methods and both response modes. Output band_vs_azimuth.png, view_direction_summary.json and azimuth_scan.csv. Image-relative azimuth 0=right,90=down in original TIFF; geographic bearing unknown. Eye-axis angles remain common image-plane model; no ray-dependent anatomical remapping introduced. Geometry and residual alignment reused, no parameter fitting.

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/test_spider_view_direction.py \
  --run data/processed/spider_prc_gavarnie_alignment_02 \
  --output data/processed/spider_prc_gavarnie_band_direction_01
```

Upload view_direction_summary.json and band_vs_azimuth.png; CSV optional. No new registration search. Actual effect on historical .470 unknown pending run.
