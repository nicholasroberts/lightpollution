# Site reports at ALAN level 5 with retained 1 degree heading sampling

Run after the heading convergence and all-site calculation already completed:

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/plot_spider_site_details.py \
  --run data/processed/spider_heading_all_sites_01/all_sites \
  --batch data/processed/spider_prc_all_sites_01 \
  --alan-level 5 \
  --atmospheric-factor 0.8 \
  --output data/processed/spider_site_details_01
```

This reads the saved 1 degree results without recalculating the nightly PRC.
Spectra must match the saved hashes. Each site gets a PNG and a page in
`all_site_details.pdf`, with:

- Natural/HPS/LED spectral shapes normalized by band energy; these are shapes,
  not the absolute source levels used in the model.
- Actual supplied photon sensitivity, without adding a UV peak.
- HPS and LED heading-median absolute PRC through the night at ALAN level 5,
  dated local axes, and a capture-time marker. No smoothing or new fitting.
- Experimental camera DoLP, natural-only, natural + HPS and natural + LED maps
  at the exact capture time, using a common 0–1 color scale.
- A 60 degree full-field circle and mean/median DoLP of valid pixels in that
  zenith cap. No terrain mask is invented; measured contamination remains.
- The recorded historical upward radiance, Moon geometry, atmospheric factor,
  and source-specific ALAN/natural band-energy ratios at capture zenith.

Camera maps with a reviewed image Moon anchor are reprojected north-up using
the existing registration helper and its report-hash check. Other experimental
maps are shown in their native camera frame with north explicitly unknown;
these are not spatially aligned comparisons to the north-up model maps.
Missing or inconsistent measured inputs produce a labelled panel and recorded
reason, not fabricated experimental data. `site_report_details.json` records
map statistics, orientation/missing-data status and actual scenario ratios.

Measured camera-channel DoLP is distinct from the pigment-weighted model DoLP.
The report does not imply a measured spectral response match. New-site ALAN
remains conditional on radiance-ratio scaling to Oxford; existing Cottage and
Oxford settings remain unchanged. The whole-night PRC curve is a forward
prediction, not a measured time series. A site report is made only for sites
present in the supplied all-site run; pending sites remain in its input inventory.

Checks: two tests verify scenario selection/time ordering and known experimental
DoLP reconstruction without inventing north. A two-site synthetic end-to-end
render checks PDF/PNG generation, native-camera and missing-data panels. It is
a software check, not validation against the full experimental data.
