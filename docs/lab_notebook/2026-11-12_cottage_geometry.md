# Cottage fisheye geometry audit — 2026-11-12

Investigator confirmed Sigma 8mm f/3.5 EX DG Circular Fisheye on Nikon D850, with camera and rear polarizer rotated together at nominal 0/45/90/135 degrees. Original 60-degree pixel calibration is unavailable. Cottage chosen because its sky boundary is less obstructed than Gavarnie.

Added src/calibrate_spider_fisheye.py. Reads each numbered TIFF independently, uses blue channel (index 2) by default, fits a robust circle to the largest illuminated component's outer boundary at 2%, 5% and 10% of the 5th-to-95th-percentile intensity range. Saves source hashes, pixel geometry, residual diagnostics and overlay. Inputs untouched; outputs require a new directory.

Nominal equisolid assumption: r(theta)=2f sin(theta/2); if fitted rim is 90 degrees, r60/r90=1/sqrt(2). Sigma lens instruction text identifies equisolid projection (manufacturer manual reproduced at https://www.gebruikershandleiding.com/Sigma-8mm-f3.5-EX-DG-Circular-Fisheye/preview-handleiding-734759.html); nominal mapping is not an empirical lens calibration. Optical rim, horizon, vegetation and vignetting must not be conflated. Camera tilt cannot be diagnosed solely from a circular rim. Radius here means 60 degrees from optical axis.

Synthetic displaced uniform disk: centre recovered within 1 px and radius within 2 px at all three thresholds; nominal radius conversion checked. Original Cottage TIFFs exist only on EcoVis and have NOT been fitted here. Threshold spread is a diagnostic, not a confidence interval. All geometry is provisional and not automatically consumed by PRC reconstruction. No registration, Stokes or PRC outputs changed. Camera rotation/analyzer coordinate conventions still require validation during registration.

Run on EcoVis:

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/calibrate_spider_fisheye.py \
  --folder data/raw/skyimages/Cottage/OneDrive_10_12-11-2026_extracted \
  --output data/processed/spider_prc_cottage_geometry_01
```

Review geometry_overlay.png and geometry_calibration.json before accepting any calibration or transferring geometry to another site/crop.
