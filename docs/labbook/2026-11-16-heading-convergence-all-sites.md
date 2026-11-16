# Heading-grid sensitivity, followed by the historical-site extension

This executes the investigator's requested sequence: check 10/5/1deg headings
with the current normalized model and frozen spectral/ALAN settings; then add
other historical sites to the forward workflow using 1deg as the numerical
reference. It does not fit a new atmosphere, DoLP amplitude or ALAN level.

## Run on EcoVis

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/run_spider_heading_convergence.py \
  --from-run data/processed/spider_normalized_prc_02 \
  --batch data/processed/spider_prc_all_sites_01 \
  --reviewed-anchors data/processed/spider_dark_site_validation_02/corrected_moon_anchors.json \
  --add-sites \
  --output data/processed/spider_heading_all_sites_01
```

Use the saved normalized run actually containing hps/provenance.json and
led/provenance.json. The _02 folder is the dated-axis run from this conversation.
The corrected four-site anchors preserve the reviewed Cottage/Oxford times and
the current Skillinge/Slabhouse ephemerides. No new source finder is used.

This is a heavier calculation than the old 36-heading runs: the reference grid
has360 headings. Progress prints per spectrum and site. Output must be fresh.
Do not treat an actively computing process as hung merely because one site
takes longer; each field/time sample is still evaluated explicitly.

## 1. Heading check

The original spectrum, natural law, neutral separation, atmospheric factor,
ALAN config, ephemerides, time step, field rays and tilt are read from the saved
provenance. Existing spectrum hashes and parsed ALAN config are checked before
reuse. A single 1deg run per spectrum then provides exact nested 5/10deg samples.
At any selected heading, the field quadrature does not depend on the number of
other headings in the grid, so this is equivalent to separate runs with those
same settings. No interpolation or smoothing produces the statistics.

Look in **data/processed/spider_heading_all_sites_01**:

- hps_heading_convergence.png and led_heading_convergence.png: median |PRC|
  at10/5/1deg headings for the existing sites and level5 (or last supplied level).
- heading_grid_statistics.csv: minimum, median and maximum over headings at
  every valid time/scenario, plus differences from1deg.
- convergence_summary.json: maximum and95th percentile absolute discrepancies
  at10deg and5deg versus1deg, across all retained scenarios. There is no arbitrary
  pass threshold and1deg is not claimed to equal a continuous-heading limit.
- hps/ and led/: actual1deg field/time response CSVs, heatmaps and range plots.

A zero lower envelope can reflect an analyzer/body orientation that cancels the
signal; it is not absence of polarization everywhere. If temporal fluctuations
remain at1deg, heading discretization alone does not explain them. Time/field
quadrature convergence would be a separate check. All such ranges are sampled
scenario ranges, not confidence intervals or behavior probabilities.

## 2. Additional sites

After the grid check, all_site extension proceeds under **all_sites/**.
The nightly model now accepts --sites instead of being hard-coded to two sites.
This is a forward world-coordinate prediction: astronomical Moon azimuth comes
from time/location, so a camera image bearing is unnecessary. Camera alignment
is still needed before overlaying a measured map or interpreting its AoP.

Inputs are prepared for all18 historical reference sites:

- Existing reviewed ephemerides are preserved.
- Other sites use original EXIF capture time from1.tif identified by the saved
  registration report, plus reference latitude/longitude.
- UK sites assume Europe/London, Swedish sites Europe/Stockholm, Gavarnie
  Europe/Paris; date-specific daylight saving is applied. This country-based
  local camera-time convention is recorded as an assumption at new sites.
- EXIF dates must match the reference field date or the following calendar day
  (an observing night may cross midnight). TIFF export/fallback timestamps,
  inconsistent dates, inaccessible raw files and unsupported Moon elevations
  remain pending review. No arbitrary replacement time is assigned.
- Additional anchors contain no selected image source. They are labelled
  ephemeris-only and cannot be used as verified measured-camera overlays.

Look at all_sites/site_input_status.csv FIRST. It lists every reference site,
actual radiance mean/SD, time source/zone, source scale basis and pending reason.
all_sites/site_prediction_status.csv additionally identifies sites with no valid
lunar field/time samples: those are not filled with zero PRC.

ALAN treatment is an explicitly conditional extension:

new_site_scale = saved_Oxford_scale * reference_site_radiance/reference_Oxford_radiance.

Existing scales remain unchanged, including the earlier user-adjusted Cottage
scenario. HPS and LED retain their different spectral/DNB factors. The actual
radiances are historical multiyear upward measurements, NOT capture-time ground
skyglow. New sites are not fitted to measured DoLP/PRC and are not assigned an
observed lamp mixture. Because Cottage retained a prior manual target adjustment,
its scale need not follow the new-site radiance proportionality; the inventory
makes this explicit. This remains a sensitivity scenario needing interpretation.

Outputs:

- all_sites/hps/ and all_sites/led/: each ready site's dated heatmaps, PRC ranges,
  response CSVs and provenance, at1deg headings.
- all_sites/all_sites_prc_overview.png: for each ready site, median over valid
  times of the heading-median |PRC|, at the selected scenario level. Sites have
  different dates/Moon tracks; this is not a simultaneous controlled field trial.
- all_sites/all_site_ephemeris_anchors.json, hps_site_alan.json, led_site_alan.json,
  extension_notes.json: time and source-scaling assumptions, hashes and pending list.

The extension can also be prepared without expensive runs:

```bash
python src/extend_spider_nightly_sites.py \
  --from-run data/processed/spider_normalized_prc_02 \
  --batch data/processed/spider_prc_all_sites_01 \
  --reviewed-anchors data/processed/spider_dark_site_validation_02/corrected_moon_anchors.json \
  --prepare-only \
  --output data/processed/spider_all_site_inputs_01
```

## Checks completed in this environment

Seven unit tests cover nested grids, incomplete-grid rejection, unchanged
physical settings, preservation of existing source scales, conditional radiance
scaling, timezone selection, rejection of export timestamps and EXIF-to-UTC
forward-anchor preparation. Existing nightly tests are rerun.

A two-site coarse-time/small-field fixture executes the entire1deg check,
5/10deg subsampling and subsequent site extension. Only the fixture's two
existing ephemerides are ready; the other16 correctly remain pending because
raw TIFFs are absent here. This is an end-to-end software check, not a result
for the investigator's full spectral inputs. Actual convergence discrepancies
and the number of ready historical sites come from the EcoVis command above.
