# Cottage and Oxford Business Park Wang comparison — 2026-11-12

Investigator selected Cottage as low-ALAN and Oxford Business Park as high-pollution example. Added compare_wang2016_measured_skies.py to reuse all-sites direct Stokes arrays. Same T3, wavelength450nm and explicit illustrative neutral separation24deg applied to both. Corrected x/y spectral branch and explicit source-meridian spin2 convention used. No model amplitude/ALAN fit; estimates one constant AoP reference offset in the zenith field and reports axial residuals.

Source geometry initially from compact brightness prominence candidate inside nominal sky circle, marked with black cross. It can mistake artificial light for Moon, especially at bright site; source identity MUST be visually confirmed. Optional --cottage-moon-xy x y and --oxford-moon-xy x y allow explicit verified coordinates. No capture UTC/clock accuracy, compass north or ephemeris inferred from candidate. Model is a provisional geometry comparison, not yet a physical ALAN attribution.

PNG shows measured/model DoLP and AoP with fixed scales, both sites; white contour 30deg-radius/60deg-full field. JSON includes source pixel/geometry, DoLP means/error, fitted global AoP offset/residual, nonphysical model counts, inherited alignment flags and input hashes. Broad sky fields include obstructions; numerical comparison restricted to zenith field pending sky masking.450nm is not the camera's broadband channel or spider quantum-catch weighting. Cottage not assumed ALAN-free.

Two targeted tests pass (axial offset without amplitude fitting; compact source candidate). Synthetic two-site workflow completes and writes comparison PNG/JSON. Actual two-site arrays remain on EcoVis; no real fit result claimed here.

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/compare_wang2016_measured_skies.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --neutral-separation-deg 24 \
  --output data/processed/wang_cottage_oxford_comparison_01
```

Upload cottage_oxford_model_comparison.png and cottage_oxford_model_summary.json. Check both source crosses before interpreting residuals. Existing low/high site distinction is investigator-provided, not inferred from current amplitude fit. Locked PRC data and historical fits unchanged.
