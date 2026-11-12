# Cottage matched polarization audit — 2026-11-12

Reviewed uploaded previews (.292714 legacy, .314837 direct Stokes) against original Cottage DoP/AoP PNGs. Broad DoP spatial arrangement appears qualitatively compatible after display-axis transpose; numerical agreement not established. Direct DoP auto colour scale extends above 3 from unphysical samples, obscuring valid sky variation. Legacy AoP shows sharper angular seams, consistent with arithmetic angle-averaging wrap artefacts, not proven solely from PNGs.

Recovered R source explicitly computes AoP=0.5 atan2 and converts arithmetic replicate mean to degrees. Original AoP PNG legend alone does not establish units or confirm the same source version. Do not silently relabel it or infer precise values from its colours.

Added src/audit_spider_polarization.py. Reuses existing linear_stokes.npz and legacy_dop_aop.npz in rotation_02; no raw images, rotation search, or model parameters changed. Outputs matched_polarization.png (DoLP 0..1, AoP -90..90 degrees) with whole-sky and FOV panels; display-only transpose to match R row/column axes, equal pixel aspect ratio. Array values unchanged; out-of-range values merely saturate display and are counted.

polarization_summary.json includes source/script hashes, intensity/DoLP quantiles, invalid/nonfinite counts, physical mean local DoLP, intensity-weighted mean local DoLP, axial orientation coherence, and existing fixed-eye results in both intensity and normalized-direction modes. A separately labelled hypothetical local-alignment benchmark removes spatial angle cancellation. This benchmark is NOT the two-fixed-eye model, must not replace it, and is not a fitted alternative selected for proximity to .439. Its purpose is to distinguish orientation cancellation from intensity differences and identify a possible historical integration distinction.

Full local suite passes 32 tests. New tests verify spatial cancellation for equal opposing Q fields, correct axial handling near +/-90deg, and explicit invalid-pixel counts. Real derived arrays remain on EcoVis; no new real summary results available here.

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/audit_spider_polarization.py \
  --run data/processed/spider_prc_cottage_rotation_02 \
  --output data/processed/spider_prc_cottage_polarization_audit_01
```

Upload matched_polarization.png and polarization_summary.json. No re-registration required. Locked PRC observations and fits unchanged.
