# Supplementary Figure 4: confirmed Moons and fitting — 2026-11-27

Investigator confirms both saved Moon pixels from the reference review sheet.
Source metadata read from attached registration_readiness(1).json.

| Site | Accepted Moon pixel | Local capture | Moon elevation | Moon azimuth | Image elevation residual | Nominal analyzer sign |
|---|---|---|---:|---:|---:|---:|
| stannes | [1412,631] | 2019-05-14 23:21:05 BST | 37.5681° | 209.1688° | −0.2300° | +1 |
| opensquare | [600,617] | 2019-05-14 23:08:05 BST | 38.5157° | 205.3284° | −0.1766° | −1 |

Exact accepted lock/Stokes hashes and UTC times are recorded in
config/paper1/supp04_reviewed_moons.json. Both registrations pass the saved
numerical checks. This confirmation establishes Moon identity; it does not
invent a new investigator acceptance of every exposure overlay.

```bash
git pull --ff-only origin paper-1-polarization-model
python src/run_paper1_supp04.py
```

Fresh derived registrations preserve original files. Preflight verifies
both sites, checksums, nominal analyzer evidence, Moon coordinates and UTC.
Existing sign helper retains U for stannes and reverses opensquare U once.
Campaign reversed camera-to-sky parity remains a separate fixed convention.
No parity is selected using polarization agreement.

Run current Figure 1 joint solar-fraction/shared-separation fits, starting
at 70 degrees independently per site. pmax=1, existing horizon treatment,
no ALAN/dilution or extra atmospheric depolarization.
Same mixture/separation for DoLP/AoP, and fitted separation in single-source panels.

Outputs:
- data/processed/paper1_natural_fit_supp04_01
- figures/paper1/supplemental/supp_fig04/v01

Four rows, two per site; measured/lunar/solar/combined columns.
Unchanged Figure 1 per-site renderer, linear viridis DoLP, cyclic romaO AoP,
north-up/east-left, 60-degree full field circle. Reports, diagnostics, summary
CSV/JSON and hashes archived beside PNG/PDF/SVG. Existing directories refused.

Validation: script compiles. Actual fits require EcoVis arrays and have not
yet been run here. Assess held-out agreement separately from training optimum.
