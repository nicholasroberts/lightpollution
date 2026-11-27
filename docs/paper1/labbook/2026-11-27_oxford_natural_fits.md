# Oxford natural-pattern fits: WClawn and Oxford Business Park

2026-11-27 — extension of the final Figure 1 method to the next two reviewed sites.

## Inputs and method

Reuse the signed, linear camera-green Stokes from
`data/processed/spider_sign_checked_01/registration`.
Both sites have previously reviewed Moon-based orientation. The runner checks
the lock and Stokes against the registration checksum manifest, and checks the
previous twilight fit refers to those same files. It retains that fit's recorded
camera-to-sky handedness; the previous DoLP-based mixture is not reused.

The camera export is an interface adapter, not a new RAW audit: saved I, Q and U
are copied unchanged. In particular, the previously applied analyzer U sign is
not applied a second time. The standard north-rotation script then transforms
Stokes into the fixed Moon-anchored sky frame.

Each site independently uses the current Figure 1 angular-structure fitter.
Start separation at 70 degrees, jointly fit solar fraction and one common
solar/lunar neutral-point separation within that site, and retain fixed source
geometry. pmax=1; existing Wang horizon treatment; no ALAN, fitted dilution or
additional atmospheric depolarization. The same fitted mixture and separation
generate both DoLP and AoP. One global analyzer offset for display is estimated
on training pixels and assessed separately on withheld pixels, as in Figure 1.

## EcoVis

```bash
git pull --ff-only origin paper-1-polarization-model
python src/run_paper1_oxford_natural_fits.py \
  --output data/processed/paper1_natural_fit_oxford_01 \
  --figure-copy figures/paper1/supplementary/solar_lunar_oxford/v01
```

Use fresh output directories. The runner preflights both sites before exporting.
It does not change existing Figure 1 or Figure 2 artifacts.

## Outputs and next decision

Each site has a camera reuse record, north-up Stokes, joint fit report, arrays,
mixture/separation profiles and fit diagnostics. The summary CSV/JSON records UTC,
solar and lunar elevations and azimuths, mixture, separation, boundary flags and
held-out RMS AoP-change errors for lunar-only, solar-only and combined models.

The supplementary PNG/PDF/SVG has two rows per site (DoLP then AoP), with measured,
lunar-only, solar-only and combined columns. DoLP is linear viridis 0–1; AoP is
cyclic romaO. Final fitted separation is used for the single-source panels too.

These are two further empirical points for the proposed solar-elevation
relationship. No elevation law or threshold is fitted from these two points.
Below-horizon solar patterns remain a geometric model approximation. Fit errors
measure spatial agreement within an acquisition, not validation across nights.
Results await execution on EcoVis; no new numerical site results are claimed here.

Validation: three synthetic provenance tests check unchanged signed I/Q/U,
changed-Stokes rejection and mismatched prior-lock rejection.

## Supplementary Figure 2 archive

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
