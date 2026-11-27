# Supplementary Figure 4 — stannes and opensquare

Selected by investigator on 2026-11-27 as the next pair. These are further
Oxford sites; minifactory remains another Oxford acquisition outside this pair.

Use the same source-check helper with explicit site names:
```bash
git pull --ff-only origin paper-1-polarization-model
python src/check_paper1_supp03_inputs.py \
  --sites stannes opensquare \
  --output data/processed/paper1_supp04_readiness_01
xdg-open data/processed/paper1_supp04_readiness_01/reference_moon_check.png
```

Inspect/upload registration_readiness.json and reference_moon_check.png from
that directory. The helper checks saved lock/Stokes hashes, reports all
candidate locks, records nominal signed analyzer evidence without applying
it, and produces a reference Moon review sheet. A saved candidate does not
automatically establish Moon identity. If no suitable source exists, retain
the pending status and resolve that input rather than borrow another site.

After source review, follow the established Figure 1 method: fixed celestial
orientation; jointly fit per-site solar fraction and common solar/lunar
neutral-point separation starting at 70 degrees; pmax=1; existing horizon
treatment; no ALAN or additional atmospheric depolarization. Use one mixture
and separation for DoLP and AoP. Compare held-out errors to lunar-only;
a training improvement alone does not establish solar contribution.

Planned outputs: figures/paper1/supplemental/supp_fig04/v01.
Four rows, DoLP/AoP per site; measured/lunar/solar/combined columns, linear
viridis and cyclic romaO, as Figure 1 and previous supplementary figures.
No site fits have yet been run or claimed for this pair.
