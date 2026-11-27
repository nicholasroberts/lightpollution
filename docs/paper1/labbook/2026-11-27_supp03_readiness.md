# Supplementary Figure 3 — ChurchHill and zoology

Next selected sites, 2026-11-27. Planned layout follows Figure 1D and
Supplementary Figure 2: two rows per site, measured/lunar/solar/combined columns,
linear viridis DoLP and cyclic romaO AoP. Jointly fit solar fraction and a
single solar/lunar neutral separation per site, starting at 70 degrees.
Retain pmax=1 and existing horizon treatment; no ALAN or added atmospheric
depolarization at this stage.

First inspect the current EcoVis registration products:

```bash
python src/check_paper1_supp03_inputs.py
```

The read-only source inspection writes
`data/processed/paper1_supp03_readiness_01/registration_readiness.json`
and, when reference images are available, `reference_moon_check.png`.
It reports all candidate locks rather than silently choosing one; checks
lock/Stokes hashes; reports existing signed analyzer convention or derives the
nominal sequence from saved landmark rotations without changing U.
Numerical registration, visual alignment acceptance, Moon identification,
ephemeris and handedness are separate checks. A saved candidate Moon alone is
not treated as a confirmed celestial bearing. Review records live on EcoVis
and have not yet been inspected here. No new fits or figures are claimed.

Once the accepted source lock and celestial orientation are established, reuse
its Stokes with the current Figure 1 fitting and rendering functions and archive
in a fresh version directory under this folder. Existing Figure 1 and
Supplementary Figure 2 outputs remain unchanged.
