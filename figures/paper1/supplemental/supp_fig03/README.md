# Accepted Supplementary Figure 3 workflow — 2026-11-27

Investigator accepts the ChurchHill four-star solution (P1 Arcturus,
P2 Pollux, P3 Regulus, P4 Procyon), RMS 1.004 degrees, reversed parity.
Zoology Moon/registration accepted earlier; pixel [1420,631].

Run:
```bash
git pull --ff-only origin paper-1-polarization-model
python src/run_paper1_supp03.py
```

Requires PyEphem (`python -m pip install ephem` if not already installed).
Outputs:
- data/processed/paper1_natural_fit_supp03_01
- figures/paper1/supplemental/supp_fig03/v01

Exact reviewed source-lock hashes are checked. ChurchHill star solution is
recomputed from saved points and checked against the accepted result.
The existing sign-convention helper checks source checksums and applies the
site-specific signed nominal analyzer convention in fresh derived trees:
ChurchHill U multiplier -1, zoology +1. The source registrations are unchanged.

ChurchHill Moon is predicted from stars and ephemeris, never presented as an
observed detection. Source date/time assumptions are retained. Zoology retains
the campaign reversed parity convention independently of analyzer sequence.
No parity choice is made by polarization fit. The north rotation transforms
I/Q/U together; DoLP and AoP are calculated afterwards.

Per-site fits start at 70 degrees and jointly fit solar fraction and one
solar/lunar separation. pmax=1, existing horizon treatment, no ALAN/dilution
or extra atmospheric depolarization. Uses unchanged Figure 1 numerical fitter
and per-site rendering. Shared layout renderer now accepts optional site names,
title and output stem; old defaults and numerical behavior remain unchanged.

Four rows, two per site; measured/lunar/solar/combined columns. Linear viridis
DoLP, cyclic romaO AoP. Copies of complete fit diagnostics/reports and summary
are included with checksums. Existing output versions are never overwritten.

Validation: Python compilation; accepted four-star solution reproduced with
PyEphem 4.2.1 to numerical precision. Real-image fits run on EcoVis.
