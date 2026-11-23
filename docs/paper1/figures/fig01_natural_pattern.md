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

## Figure 1 — Museum absolute Moon bearing

For the museum reference capture on **14 May 2019 at 22:09:51 BST**
(21:09:51 UTC), the previously calculated fixed astronomical Moon azimuth is
**187.15 degrees**, measured clockwise from **true north**:
north=0, east=90, south=180, west=270 degrees. Thus the Moon is
**7.15 degrees west of due south**. Its corresponding astronomical elevation
is **41.20 degrees**.

In the publication sky view (north up, east left, south down, west right),
the Moon must appear near the bottom, slightly **right of south**. This is
the expected bearing check for Figure 1, independent of the DoLP band or AoP.
The north rotation uses the observed Moon pixel and fixed astronomical
azimuth; it does not move the source radially to force agreement with its
astronomical elevation. The fisheye centre/projection controls that separate
radial/elevation comparison.

These values document the existing capture ephemeris; they were not newly
recomputed for this lab-book entry. The explicit reversed-handedness hypothesis
used for the museum display remains recorded separately from the analyzer
sequence sign and unknown analyzer zero.

## Figure 1 — Museum measured DoLP in the spider field

Investigator-reported EcoVis calculation, 23 November 2026:

| Quantity | Value |
| --- | --- |
| Mean measured camera-green DoLP | 0.139762 |
| Spatial sample standard deviation (ddof=1) | 0.040131 |
| Included physical pixels | 180850 |
| Excluded pixels within the field | 9 |
| Total pixels within the field | 180859 |

Paper display, rounded: **DoLP = 0.140 ± 0.040 (mean ± 1 spatial SD)**,
or **13.98 ± 4.01 percentage points** when expressed as percent polarization.

Definition: equal-pixel mean and spatial SD of per-pixel DoLP within the
zenith-centred **60 degree full field (30 degree angular radius)**, using
the locked reference centre and eye_radius_px. Calculated in the camera frame
before north rotation to avoid additional interpolation. Uses the retained
pair-based intensity I_pair and signed Q,U; DoLP=hypot(Q,U)/I_pair.
Include only pixels inside the field with positive I, finite DoLP and
0 <= DoLP <= 1. No additional clear-sky/foreground mask was applied by this
calculation.

This SD describes spatial variation within this capture, not the standard
error, a confidence interval or independent replicate variation. This is the
mean of pixel DoLP, not DoLP from spatially integrated I,Q,U and not PRC.

Sources on EcoVis:

- data/processed/paper1_museum_stokes_audit_01/museum/measured_stokes_and_polarization.npz
- data/processed/spider_registration_green_all_sites_01/runs/museum/locked_site_references.json

Associated museum reference capture: 14 May 2019, 22:09:51 BST. These numerical
values were supplied from the executed EcoVis calculation; not independently
recomputed in this workspace.
