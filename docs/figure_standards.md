# Publication figure standards

**Applies to:** `publication-figures` branch of `nicholasroberts/lightpollution`.

**Purpose:** reproducible, compact, **editable Adobe Illustrator** artwork for the scientific presentation and eventual paper. This document specifies export/presentation only. It does **not** authorise changes to frozen data, the VIIRS zero-point and 2021 corrections, fixed NUTS geography, empirical spider/moth PRC functions, or any established scientific analysis.

## Formats and editability

- **Primary master:** SVG. Ordinary labels, ticks, legends and annotations must remain editable text and vector objects.
- **Secondary copy:** vector PDF, with TrueType text embedding (`pdf.fonttype=42`).
- **EPS:** optional on explicit request. EPS transparency is unreliable; inspect alpha/transparency layers after import.
- **No PNG/JPEG as the primary deliverable.** Raster preview copies may be created separately as needed for presentations.
- Keep SVG text editable with `svg.fonttype="none"`. Mathtext or special glyph rendering may still use vector paths: inspect in Illustrator if individual symbols need to be edited.
- Ensure the **same typeface is installed on the plotting workstation and the Illustrator computer**. Default is DejaVu Sans; confirm it in Illustrator, otherwise choose a common installed font before final export. Avoid font substitution.
- All native plot marks, axes, text, legends, scale bars, colour bars and NUTS3 boundaries are **vector by default**.

## Dimensions and typography

- Standard **single-column** width: 90 mm.
- Standard **double-column** width: 180 mm.
- **Talk/slide** working width: 240 mm (may be adjusted for aspect ratio).
- Common baseline font size 9 pt for a double-column figure; titles 10 pt, tick and legend labels 8 pt. Scale text sensibly for slide figures.
- Lines ~1.1 pt; axes/ticks ~0.7 pt. Maintain sufficient contrast and avoid tiny labels.
- Use consistent panel letters, units, caption wording, tick formatting, legend order and colour-bar placement. Keep panel letters editable.
- Design within a clean white background. Avoid unnecessary outlines, shadows, transparency effects and unnecessary high-resolution embedded images.

## Colour and comparable panels

- Use colour-blind-conscious, consistently assigned line/point colours where categorical series recur. Shared defaults:
  - Radiance: `#D55E00`.
  - Population: `#009E73`.
  - Spider PRC: `#0072B2`.
  - Moth PRC: `#CC79A7`.
- **Signed European change maps:** existing project convention `RdBu_r`: blue for decreases and red for increases. Centre at zero.
- Preserve the previously chosen **symmetric-log colour treatment** for radiance changes with wide positive/negative ranges where suitable. Show the **real physical values on colour-bar tick labels**; never present scaled values as raw linear measurements.
- Map panels intended for direct comparison should have common geographic extents and scientifically appropriate, clearly documented scale choices. Do not force identical numerical colour limits across unlike physical quantities.
- For spider and moth PRC, label the independent biological measures and their true units/scales. Do not reuse radiance units.
- State whether a panel displays an absolute difference or percentage change, and specify years in the actual figure.

## Size targets and path management

- Aim for **under 5 MiB** for an ordinary SVG/PDF and **under 10 MiB** for a NUTS3 map; these are **targets, not cut-offs**. Never discard or degrade scientific content merely to hit a threshold.
- Prefer compact vector exports over massive per-point/per-pixel SVG objects. Avoid unnecessary duplicates, repeated heavy transparency layers and excessively dense markers.
- NUTS3 polygons stay vector when reasonably sized. If a map is too large, investigate visual-only polygon-path simplification in a **projected coordinate reference system**, after checking that the displayed coastline/boundaries remain suitable. **Never use simplified geometry in zonal statistics or other calculations.**
- For very dense **native-pixel** backdrops, rasterise **only that artist**, while preserving **all labels, borders, scales and annotations as editable vectors**. This is an explicit exception to the full-vector preference, recorded for each affected figure.
- Embedded raster detail, where unavoidable, should be appropriate to intended physical size (typically 300–600 dpi); no blanket high-DPI setting for pure vectors.
- Both editable SVG and PDF are exported by default. Check output size and inspect the result in Illustrator before declaring production ready.

## Files and provenance

- Common plotting defaults live in `src/figure_style.py`.
- Figure artwork is written under `figures/publication/`, **ignored by Git**. Keep scripts, style settings, documentation and source-data references versioned instead.
- Use stable descriptive filenames: e.g. `figures/publication/radiance_change_2014_2024.svg` and the matching `.pdf`.
- Each producing script must identify its source tables, fixed cohort/geography, time period, output filenames, and any display-only transformations. Prefer using **existing corrected intermediate tables** rather than recalculating analysis values.
- Never silently alter source scientific data or locked calibrations while polishing a figure; those are separate, explicit decisions.

## Minimal usage

When a plotting script is launched as `python src/plot_something.py`:

```python
import matplotlib.pyplot as plt
from figure_style import apply_figure_style, figure_size, save_figure

apply_figure_style()

fig, ax = plt.subplots(figsize=figure_size("double"))
# draw the already-validated results ...
save_figure(fig, "figures/publication/example", map_figure=False)
plt.close(fig)
```

The example creates `example.svg` and `example.pdf` in the publication folder. For an EPS deliverable, pass `formats=("svg", "pdf", "eps")` explicitly.

## Before accepting a figure

1. Verify it uses **the intended frozen input table and cohort** and the correct scientific transformation.
2. Open the SVG in Illustrator; verify labels remain separately editable and fonts are not substituted.
3. Zoom in on maps, lines, colour bars and legends; check no clipping, accidental gaps, missing glyphs or unwanted rasterisation.
4. Check SVG/PDF file sizes and remove unnecessary graphical complexity if needed, without modifying data.
5. Confirm output names, source references and figure-specific settings are reproducible from the producing script.
