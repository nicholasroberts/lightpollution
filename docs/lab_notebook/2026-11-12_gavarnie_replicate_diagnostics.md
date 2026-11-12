# Gavarnie replicate diagnostic — 2026-11-12

Reviewed investigator-uploaded rotation_analysis(3).json and overlay. Multistart recovered image 8 at 90.012deg (correlation .4771), image 12 at 136.596deg (correlation .4052), now matching their groups. Nominal-field direct Stokes intensity PRC .2328957529 and legacy .2341530687, below reference .470 despite correction of gross rotation errors. Residual offsets remain visible. +/-5% radius has very little effect. Do not claim alignment correction reproduced historical PRC.

Added diagnose_spider_replicates.py. Reads saved rotation report, verifies raw image SHA256, reapplies exact saved centres/rotations; NO registration search. For each genuine interleaved quartet and combined result, calculate both processing methods, both response modes, original-reference differences, mean DoLP, integrated DoLP, max eye-axis angle, intensity/Q/U sums, invalid-pixel fractions and analyzer pair closure. Retains corrected 30deg radius/60deg full field. Two PNGs show individual and combined intensity/DoLP/AoP cropped for detail around the field, with fixed polarization scales. Output replicate_summary.json preserves provenance. Raw inputs and old output folders unchanged.

Two targeted tests pass: identical replicates preserve contrast; opposite polarization axes retain individual contrast but cancel in combined linear Stokes. No real replicate diagnostic results claimed before EcoVis execution.

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/diagnose_spider_replicates.py \
  --run data/processed/spider_prc_gavarnie_alignment_02 \
  --output data/processed/spider_prc_gavarnie_replicates_01
```

Upload replicate_summary.json, linear_stokes_replicates.png and legacy_dop_aop_replicates.png. These distinguish low contrast common to individual replicates from cancellation or processing effects during replicate combination. Existing alignment imperfections are reused, not declared resolved.
