# Tilt-aware registration revision using the saved landmark coordinates

The supplied first lock contains hundreds of matches at failed sites, but
two-dimensional rigid transforms leave structured residuals. The missing
experimental panels are intentional rejection of these failed alignments, not
missing original TIFFs. Among the fifteen sites in the spectral report, only
Cottage, opensquare, WCpitch and zoology passed the original registration.

The failure mechanism is consistent with small tilt changes in fisheye camera
orientation, which are not equivalent to planar image rotation/translation.
This is supported by a separate equisolid ray-rotation mapping reducing held-out
landmark residuals. It does not establish a fully calibrated lens model.

## Reuse the saved coordinates and regenerate a fresh reference set

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/lock_spider_site_registration.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --model-run data/processed/spider_heading_all_sites_01/all_sites \
  --from-lock data/processed/spider_registration_locked_01/locked_site_references.json \
  --registration-model fisheye \
  --output data/processed/spider_registration_locked_02
```

No new SIFT detection is needed. Raw TIFF hashes and source-report hashes are
still checked. The registration module retains every previously passing rigid
transform. Only failed exposures get a tilt-aware mapping fitted from the saved
training correspondences. The held-out split remains every fourth match.

The new mapping is an SO(3) ray rotation through the existing equisolid projection,
plus a two-dimensional target projection-centre offset. The input centre and rim
radius remain fixed at their inherited values. Centre offset is bounded to 5%
of rim radius and weakly regularized toward zero with a recorded 0.2 pixel prior
weight. A fit at the bound is rejected. There is no arbitrary elastic warp,
reflection, fitted Moon track, or DoLP/AoP optimization. This is a registration
hypothesis subject to visual and calibration checks, not a physical polarization
or full lens Mueller calibration.

The original median-error, inlier-fraction and landmark-coverage thresholds are
unchanged. This development reused the recorded held-out set when assessing a
new mapping family; it is not a new independent validation dataset. Image-level
review remains essential, particularly outside the outer-field landmark mask.

## Findings from the supplied real-site landmark data

With the final bounded mapping and original passing transforms preserved:

- All exposure checks pass at 13/18 sites: ChurchHill, Cottage, EastHorrington,
  Gavarnie, PennCC, Skillinge, WClawn, WCpitch, opensquare, Oxford business park,
  southhorrington, stannes and zoology.
- Ounsdale still fails exposures 4–6 (held-out medians about 2.1–2.4 pixels).
- museum still fails 9–12 (about 2.0–2.6 pixels), though much improved from the
  original roughly 6–9 pixel residuals.
- minifactory, slabhouse and Waitrose also remain unresolved.

Those counts describe fitting the supplied landmark coordinates here, not yet
the user's re-rendered experimental products. Sites without usable astronomical
ephemerides remain excluded from the spectral model even when exposure
registration succeeds. No failed site is silently included or assigned zero DoLP.

Inspect **SITE/registration_edges.png** and **SITE/registered_dop_and_north.png**.
**registration_summary.json** now provides a compact per-site/per-exposure status
file without the millions of landmark coordinate entries in the full lock.

Then plot:

```bash
python src/plot_spider_site_details.py \
  --run data/processed/spider_heading_all_sites_01/all_sites \
  --batch data/processed/spider_prc_all_sites_01 \
  --reference-lock data/processed/spider_registration_locked_02/locked_site_references.json \
  --alan-level 5 --atmospheric-factor 0.8 \
  --output data/processed/spider_site_details_locked_03
```

All astronomical Moon coordinates, 1-degree PRCs and ALAN/atmosphere/spectral
settings stay fixed. Existing Moon identity limitations remain explicit; a better
exposure registration alone does not confirm the Moon or recover true north.

Three new tests cover projection round trips, tilt-aware recovery with false
matches, and inverse-warp direction/no reflection. A synthetic TIFF rerun exercises
reuse of a saved lock, one failed exposure being refitted, product generation and
the new compact summary. The existing locked-reader/rigid registration tests are
also rerun. Original outputs remain unchanged.
