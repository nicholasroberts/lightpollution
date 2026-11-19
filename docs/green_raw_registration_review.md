# Green RAW registration and visual review

Run from the repository root with the prepared green batch already present:

```bash
python src/register_spider_raw_batch.py \
  --raw-batch data/processed/spider_linear_green_all_sites_01/batch_summary.json \
  --batch data/processed/spider_prc_all_sites_01 \
  --model-run data/processed/spider_heading_all_sites_01/all_sites \
  --output data/processed/spider_registration_green_all_sites_01
xdg-open data/processed/spider_registration_green_all_sites_01/review.html
```

The runner validates prepared-manifest hashes before running fresh fisheye-shared landmark registration per site. Each child process writes its own log and lock under runs/SITE. Missing prior reference geometry is recorded as pending, never borrowed from another location. Failed sites do not stop the remaining attempts; exit status is nonzero if any site is pending or failed. Existing output directories are not overwritten. Optional --sites selects canonical names, and --moon-coordinates forwards investigator Moon pixels to the underlying registration script.

The HTML page links full-resolution red/cyan overlays and reference Moon/DoLP panels, lists numerical checks, and records separate alignment and Moon-identity decisions. Download alignment_visual_review.json before closing the page; restore it using the file selector when returning. There is no automatic persistence. Exported decisions are tied to the summary SHA256 and do not modify numerical locks or authorize model comparisons.

Numerical success reconstructs camera-channel I/Q/U and DoLP using the existing nominal analyzer groups. This batch runner does not compute PRC, fit a sky pattern, validate camera handedness or establish independent angular calibration. North remains provisional until Moon identity and geometry are reviewed. Later AoP/PRC comparisons must retain the calibrated Stokes basis and mask checks. The previously prepared TIFFs are black-subtracted linear RAW samples, not rendered RGB TIFFs.

Review registration_batch_summary.json and export alignment_visual_review.json for the next analysis step. Prior reference report availability, rather than RAW availability alone, controls which sites can be processed by this runner. Wells Cathedral remains its own site and may require a new reference geometry report.
