# Trial image-bearing correction after visual Moon audit

Run the existing 1 degree all-site results through an image-alignment-only trial:

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/try_spider_audit_bearings.py \
  --run data/processed/spider_heading_all_sites_01/all_sites \
  --batch data/processed/spider_prc_all_sites_01 \
  --output data/processed/spider_bearing_trial_01
```

Inspect **bearing_audit/moon_positions_overview.png** first, then
**reports/all_site_details.pdf**. Individual figures are in both folders.
The experiment is explicitly provisional; neither brightness nor ring proximity
alone independently establishes Moon identity.

The 16 November 18-page native-camera audit suggested clear sources at
Cottage, minifactory, museum, opensquare, Oxford business park, Skillinge,
slabhouse, southhorrington, stannes, WCpitch and zoology.

The existing four Cottage/Oxford/Skillinge/slabhouse bearings are preserved only
if their report hash and candidate radius remain consistent. Seven additional
sites use approximate image-angle sectors derived visually from that audit:
minifactory -3 degrees, museum -174, opensquare -177, southhorrington 180,
stannes -3, WCpitch 180, zoology 0. Angles are measured from image +x toward +y.
These are search seeds, not measured precision coordinates.

Within a +/-12 degree sector and +/-6 degree elevation band, the code finds
local intensity contrast from saved I only, using Gaussian widths 2 and 15 pixels.
A nonpositive contrast, missing input, hash mismatch or peak close to a search
boundary leaves the site unaligned with an explicit reason. There is no global
brightest-source search, no selection using DoLP/AoP, and no fitting of a Moon
position to a polarization pattern. The candidate's radial residual is recorded
rather than forcing it onto the elevation ring. Tree edges or lamps can still
win locally; inspect the new audit before accepting the trial bearings.

ChurchHill, Ounsdale, WClawn and Waitrose remain in camera coordinates.
EastHorrington, Gavarnie and PennCC are not added to the forward model because
their capture ephemerides remain unresolved.

Only display-bearing metadata in a separate derived run changes. All UTC times,
Moon elevations/azimuths, 1 degree PRCs, radiances, spectra, atmosphere and ALAN
levels remain unchanged. Saved PRC CSVs are copied and their hashes verified.
Original files are not modified. Model north-up maps remain the same;
experimental maps with accepted provisional bearings are reprojected.
The report's orientation labels remain explicitly provisional.

**bearing_review.json** records seeds, selected coordinates, elevation residuals
and rejection reasons. **display_run/** holds derived report inputs and source
provenance hashes. This is not a new physical model run.

Validation limitation: the execution host disconnected during this change, so
this new wrapper has not been executed here. On invocation it first runs a
synthetic check that a brighter opposite source cannot steal the restricted
selection, and then verifies copied PRC hashes. Full validation requires the
local run and inspection of its selected sources.
