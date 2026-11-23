# Figure 1 — Natural solar–lunar polarization

## Agreed scientific progression — 23 November 2026

1A establishes experimental green-channel DoLP and AoP. Preferred contrasting
examples are Cottage (low-ALAN candidate) and museum (stronger ALAN and twilight).
1B shows lunar-only natural predictions; 1C adds the solar contribution to test
whether band position and spatial AoP agreement improve. Add two further
experimental/theoretical comparisons across solar elevations after checks.

The intended endpoint is that spatial structure and AoP are explained, while
natural-model DoLP remains too high. This is a hypothesis to test, not a result
to impose on site selection. Figure 2 introduces the unpolarized ALAN contribution.

## Calculation and validation

Combine solar and lunar I,Q,U in one fixed coordinate frame. Show natural
predictions with unpolarized background B=0; do not display the old diluted
predictions as natural-only. Record how the solar/lunar ratio was estimated.
If reusing a ratio jointly fitted with a nuisance background, disclose that
dependency: setting B=0 for display does not make its estimation independent
of dilution. Assess robustness before treating geometry as explained.

Hold landmark registration, observed Moon anchor, parity and analyzer signs
fixed independently of the sky pattern. No pattern rotation or reflection fit.
Use common clear-sky support for numerical comparisons and spatially withheld
regions. The full visible hemisphere may provide visual context, but quantitative
claims must not extend beyond checked support. Foreground remains identifiable.

Absolute analyzer zero is uncalibrated. Report spatial AoP changes on common
support as the independent orientation check. Displayed absolute AoP maps must
state the unknown constant zero; do not silently fit an offset.

## 1A implementation

`src/plot_paper1_experimental.py` renders paired measured DoLP/AoP maps for any
specified sites, from checksum-verified linear-green Stokes and matching saved
fit reports. Reports supply only the saved display geometry, not fitted
intensities. It transforms the polarization tensor with the image coordinates,
uses a sky-view compass (N up, E left), common scales, observed Moon marker and
central 60 degree full-field circle. No scientific parameter is refitted.
The camera-plane angle convention and uncalibrated off-axis transport remain
explicit limitations.

Cottage is a preferred candidate, not yet confirmed ready in this new pipeline.
Earlier uploaded Cottage comparisons used older reconstruction, provisional
geometry and a fitted constant AoP offset; they must not be relabelled as a
current validated green reconstruction. Museum likewise needs its matching
green lock and report; do not mix its separate lock with the Oxford pilot lock.

Example on EcoVis once both sites are in a common verified lock/report set:

```bash
python src/plot_paper1_experimental.py \
  --reference-lock PATH/locked_site_references.json \
  --fit-root PATH/checked_fits \
  --sites Cottage museum \
  --output data/processed/paper1_fig01a_01
```

The fit root layout is SITE/fit/twilight_moon_fit.json. If locks are separate,
render each site separately with its matching lock and a fresh output directory.
Missing/currently incompatible input fails explicitly; no substitute map is made.

## Deliverables

PNG, editable PDF/SVG and source-hash provenance for 1A. Whole-figure assembly,
natural B=0 predictions and all-site quantitative comparison remain next work.
Solar-only controls stay in companion analyses. No global twilight cutoff is
inferred. Theory: [shared framework](../../theory/theory.md).
