# Audit Moon identity before further spatial interpretation

The investigator reports incorrect Moon positions at multiple sites. A non-null
saved candidate was insufficient evidence of a verified image bearing. New-site
ephemerides also assume local camera time from country, so timing may require
review independently of bearing. No global 180 degree correction is justified.

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/audit_spider_moon_positions.py \
  --run data/processed/spider_heading_all_sites_01/all_sites \
  --batch data/processed/spider_prc_all_sites_01 \
  --output data/processed/spider_moon_audit_01
```

First inspect `moon_positions_overview.png`. The PDF `moon_position_audit.pdf`
and individual PNGs show native-camera intensity with linear and logarithmic
stretches, experimental DoLP, original pixel axes, saved local/UTC times and
the currently calculated ephemeris. Cyan is the predicted elevation ring,
red is an existing candidate, explicitly requiring visual confirmation.
The audit neither selects a bright source nor fits DoLP/AoP. A ring alone
does not identify the Moon: streetlights and reflections can fall on it.

`moon_audit_details.json` records metadata and report-hash agreement. The
optional `moon_review_template.json` is a blank review record, not an input to
the model. Record actual Moon pixels and independently confirmed timing there
for a subsequent reviewed correction. Missing sites remain labelled.

If the real Moon lies at the expected radius but another bearing, inspect the
camera orientation/source identity. If its radius differs, inspect timestamp,
time zone, geometry and horizon calibration. These are diagnostic possibilities,
not automated attribution of the error. The same equisolid projection as the
current workflow is assumed; this audit does not validate the lens calibration.

The site-detail script now leaves measured maps in camera coordinates by
default even if a candidate exists. `--use-saved-bearings` explicitly opts into
old candidate-based orientation, labelled provisional. Previous reports are not
modified. Retain 1 degree heading sampling; no ALAN, atmosphere or Moon model
parameters are changed here. Incorrect timestamps would affect nightly PRCs;
incorrect camera bearings alone affect image alignment.

Software checks: synthetic native-camera/missing-input audit renders completed;
elevation-radius forward/inverse checks cover 0, 10, 45 and 90 degrees. Actual
Moon identifications require the original site images and investigator review.
