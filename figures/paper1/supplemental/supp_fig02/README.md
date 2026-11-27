# Supplementary Figure 2

WClawn and Oxford Business Park: the Figure 1 comparison format.

Generate on EcoVis from the completed fits (no refit):

```bash
python src/plot_paper1_supp_fig02.py
```

Default output: `figures/paper1/supplemental/supp_fig02/v01/`.

Four rows: WClawn DoLP, WClawn AoP, Oxford Business Park DoLP, Oxford
Business Park AoP. Columns: measured green, lunar-only, solar-only, best combined.
Uses the same validated renderer as Figure 1D, linear viridis for DoLP and
cyclic romaO for AoP, north up/east left, the 60° full field circle, and the
same fitted mixture/separation for DoLP and AoP. Single-source panels use
each site's fitted separation. One training-derived global analyzer offset
is used consistently across model panels, as in Figure 1.

PNG/PDF/SVG, figure values and copies of the saved fit reports, summary,
profiles and diagnostic plots are archived together with SHA256 checksums.
Existing version directories are never overwritten. Generated files are
created on EcoVis where the fitted arrays reside.

Reported EcoVis results, 2026-11-27:

| Site | Sun elevation | Solar fraction | Lunar fraction | Separation | Lunar-only held-out RMS | Combined held-out RMS |
|---|---:|---:|---:|---:|---:|---:|
| WClawn | −12.5303° | 0.542604 | 0.457396 | 68.3852° | 13.0659° | 4.8941° |
| oxfordbussinesspark | −18.8369° | 0 | 1 | 66.3410° | 4.2011° | 4.2011° |

RMS values describe spatial AoP-change errors. Oxford Business Park reaches
the lunar-only boundary; no solar component is required by this fit.
Fractions are model mixture weights, not independently measured source radiances.
