# Cottage target DoLP0.45 and levels5,6,7 — 13 November2026

Investigator request: raise Cottage ALAN so calculated DoLP is approximately0.45, then examine levels5,6,7. Operational interpretation stated to investigator: target the mean LOCAL photon-sensitivity-weighted450–750nm DoLP inside the zenith-centred60degfull field at Cottage capture time, HPS at level5. This is not the DoLP of spatially integrated Stokes, a peak DoLP, or a direct fit to experimental camera pixels.

The explicit --cottage-target-dop0.45 --cottage-target-level5 options adjust Cottage's common source amplitude once. The same adjusted Cottage source scale passes to the LED scenario with its pre-existing DNB-overlap factor, so LED's DoLP/PRC is not separately forced to0.45. Oxford source scale and existing spectra, geometry, extinction and atmospheric factor remain unchanged. Levels6and7 increase ALAN linearly relative to5, giving lower Cottage DoLP. Atmospheric factor default0.8.

This explicitly supersedes the earlier115-fold site-ratio assumption FOR COTTAGE. The original radiance-based settings remain in provenance but are marked superseded; neither the new effective site ratio nor the resulting ALAN level is an independent real-world measurement. The adjustment is user-authorized scenario design, not an inferred light-pollution estimate. Original experimental maps are never read by the adjustment helper and remain visual checks.

`adjust_cottage_alan.py` uses8192equal-solid-angle cap samples and a bracketed scalar solve. Exact capture-time KS91 relative radiance is used; common phase and distance factors cancel between sky rays and zenith at that same instant. An unreachable target above the no-ALAN field mean is rejected, as is a target requiring a reduction from the existing Cottage ALAN scale. The adjustment metadata is written to alan_configs/cottage_adjustment.json. Oxford remains unchanged. Plot headers and settings record the target rather than claiming wholly independent normalization.

Local numerical run: HPS level5 meanDoLP0.4500; Cottage source scaling increased30.4times from the previous very low Cottage scenario. Three new tests check recovery of a known target, monotonic reduction at levels6/7, and rejection of an unreachable target. Existing numerical suites remain applicable. No actual experimental arrays were refitted.

## EcoVis

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/run_independent_site_alan.py \
 --batch data/processed/spider_prc_all_sites_01 \
 --spectra data/processed/spider_spectral_inputs_01 \
 --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
 --cottage-target-dop 0.45 --cottage-target-level 5 \
 --levels 5 6 7 \
 --output data/processed/spider_cottage045_levels567_01
```

In the command above separate option names from numbers with spaces: --cottage-target-dop 0.45 --cottage-target-level 5 --levels 5 6 7.

Start with comparison_level_5/visible_spectral_nightly_comparison.png beneath that output folder; then compare levels6and7. All three include experimental maps and modelled natural/HPS/LED skies. Previous results are not overwritten.
