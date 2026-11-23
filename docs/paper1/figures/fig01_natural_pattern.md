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

## Figure 1B — Explicit lunar-only baseline

Added src/model_paper1_lunar_baseline.py. It consumes the matching locked
reference and north-rotated audited museum arrays. No production data were
run here; EcoVis rendering is required.

| Assumption | 1B implementation |
| --- | --- |
| Natural source | Moon only, fixed capture ephemeris |
| Solar contribution | Zero; introduced in 1C |
| Added unpolarized dilution / ALAN | Zero |
| Intrinsic maximum DoLP | Explicit --pmax; use 0.8 for continuity and 1.0 as ideal-Rayleigh sensitivity |
| Horizon depolarization factor | None |
| Projection | Existing zenith-axis equisolid geometry, rim at zenith angle 90 degrees |
| Atmospheric radiative transfer | Not implemented in this geometric baseline |
| Spectral comparator | Camera-green measurement; achromatic polarization with separable lunar spectrum |
| Coordinate frame | Same north-up/east-left sky frame as measured maps |
| Parameter/rotation/zero-angle fitting | None |

For scattering angle gamma, DoLP=pmax*sin(gamma)^2/(1+cos(gamma)^2).
pmax=0.8 already lowers intrinsic amplitude relative to ideal Rayleigh:
it is not a fitted ALAN background and not an independently calibrated local
atmospheric value. Keep it constant across 1B and 1C. Do not import an older
horizon factor silently. Pure multiplicative attenuation common to a source's
I,Q,U cancels from its DoLP, whereas added light or multiple scattering need not.

Under S(lambda,v)=e_M(lambda) S_geometry(v), spectral integration multiplies
all three components by the same response integral, which cancels from DoLP
and AoP. Thus this achromatic lunar-only panel is comparable to camera-green,
but is not claimed as a calibrated RH1-weighted prediction. RH1/ELF response
weighting remains explicit in later spectral/PRC comparisons; differing source
spectra can change their relative Stokes weights.

Exports measured versus lunar-only DoLP/AoP PNG/PDF/SVG, model TIFFs,
Stokes arrays and source-hash JSON. Report common central-field statistics as
north-resampled diagnostics, not replacements for the published unrotated
museum mean. Full-disc context is not a verified clear-sky mask. Unknown
constant measured analyzer zero prohibits an absolute-color agreement claim.

Validation here: syntax and analytical Rayleigh limiting-angle/amplitude checks.
Earlier phase/projection helpers reused unchanged. End-to-end museum rendering
and independent spatial AoP assessment remain pending. No horizon or pmax
assumption should be described as experimentally established.

## Correction to Figure 1B baseline — no imposed polarization ceiling

Investigator rejected the newly introduced pmax=0.8 Rayleigh-shaped baseline:
Figure 1B must use the earlier basic Berry/Wang calculation. This entry
supersedes the preceding 1B pmax and no-horizon assumptions.

Updated model_paper1_lunar_baseline.py to reuse wang2016_sky.natural_stokes:
Berry-normalized pattern (Berry eq4.2 -4 prefactor in source-meridian convention)
times the existing Wang eq8 horizon factor cos(theta)^(1/10), unity atmospheric
multiplier, uniform natural I=1, no solar contribution and no unpolarized
addition. No separate pmax ceiling or clipping. Nonphysical kernel outputs
raise an error instead of being clipped.

Neutral-point separation must be explicitly supplied from the prior declared
configuration; the script does not invent a default. This is the earlier
reduced angular/horizon implementation, not a claim to full Wang spectral
eq17 or twilight radiative transfer. Its geometric polarization envelope used
as DoLP and its phase-to-camera basis remain declared approximations.
The prior --pmax command is superseded; do not run it.
Validation: syntax and explicit unity-multiplier/no-pmax checks; real EcoVis
rendering and source-basis comparison still required.

## Figure 1B — Two-site theoretical display and statistics

Layout: 2 rows (Cottage, museum) x 2 columns (DoLP, AoP); labels A–D.
Shared DoLP scale 0–1 and cyclic AoP scale -90 to +90 degrees; N up/E left,
ephemeris Moon marker, central zenith-centred 60-degree full-field circle.
Annotate each DoLP map with theoretical mean ± one spatial sample SD.

Added plot_paper1_lunar_two_sites.py to read both saved lunar-only reports and
arrays, require the same kernel hash and 24-degree separation, reconstruct the
complete theoretical central field, and export PNG/PDF/SVG plus
figure_values.csv/json and labbook_values.md. No parameter recalculation/fitting.

Confirmed from investigator-uploaded museum lunar_only.json:
lunar-only theoretical mean **0.5462624388780614**, spatial SD
**0.2258845079975526**, 180854 pixels on the report's *common measured-valid
central support after north resampling*. This is not the same support as the
complete theoretical central circle, so retain it as a diagnostic result.
Measured original-camera value remains **0.139762 ± 0.040131** (180850 valid,
9 excluded), separately recorded above.

Cottage's real-run report has not been supplied in this workspace; do not record
a synthetic smoke-test result as a real output. The new comparison command
calculates both complete-field theoretical statistics, prints them, and writes
the exact values with counts and source hashes for subsequent insertion here.
These are pixel-DoLP spatial statistics, not integrated-Stokes DoLP, PRC,
standard errors or independent replicate estimates.
Validation: syntax plus synthetic two-site rendering, field mean and table export.

## Figure 1B label correction

All four Cottage/museum theoretical DoLP/AoP maps form **one subfigure 1B**.
Removed internal A–D labels from plot_paper1_lunar_two_sites.py and added one
figure-level B. Retain site and DoLP/AoP headings and mean ± spatial SD labels.
This supersedes the earlier internal A–D layout description. Numerical
calculations, field support and source data are unchanged.

## Figure 1C — Museum saved-fit comparison and real fit results

Investigator retains ratio profile and error plots as supplementary information.
Main subfigure 1C is two rows (DoLP, AoP) by four columns: measured north-up
museum, lunar-only, solar-only, best combined. Use batlow/romaO, shared scales,
one C label, observed/ephemeris Moon marks and 60-degree circle. No ALAN
correction has been applied to the measured maps.

Added plot_paper1_solar_lunar_panel.py. It validates the saved fit against source
lock/north hashes and component/kernel versions, recreates components only to
check saved arrays, and renders without refitting. Full-disc theoretical
context outside the central fit mask is not validation against foreground.
Unknown absolute analyzer zero remains explicit.

From investigator-uploaded solar_lunar_fit.json:
solar fraction **0.6793869859731805**, lunar fraction **0.3206130140268195**;
solar:lunar uniform-intensity ratio **2.119024981051921**.
Withheld pooled AoP-change RMS: lunar-only **6.710319578018203 degrees**,
solar-only **3.2023583961994873 degrees**, combined **2.94636942540123 degrees**.
Withheld median absolute errors: **3.4006367127863637**, **1.3964526632290803**,
**1.4001520785539507 degrees**, respectively. Combined improves RMS substantially
over lunar-only and modestly over solar-only; it does not improve the median
relative to solar-only. Avoid a blanket claim that combined improves all metrics.

2784 central mask pixels; 4853 training pairs, 1642 withheld pairs; 21 training
blocks. Leave-one-training-block-out solar fractions span
**0.6554572999541297–0.7779567512322622** (conditional spatial stability, not CI).
The fitted minimum is interior but broad. No independent-pixel significance claim.

Solar geometry: elevation **-10.227399587952283 degrees**, true-north clockwise
azimuth **319.08424097267704 degrees**, at 2019-05-14 21:09:51 UTC.
Statistics are for the stated model and mask, not calibrated atmospheric
solar/lunar radiance fractions. The actual source report is the uploaded
solar_lunar_fit.json, not the synthetic test fixture.

Validation: syntax and synthetic saved-fit component/hash checks plus two-row,
four-column PNG/PDF/SVG export. Actual museum render remains an EcoVis command.
