# Paper 1 — Polarization modeling

Working branch: `paper-1-polarization-model`.

Start with [the analysis plan](analysis_plan.md), [decisions](decisions.md),
[all-site ledger](site_status.csv), [figure status](figure_status.md) and
[shared theory](../theory/theory.md). This is research infrastructure and
technical documentation, not manuscript text.

## Working layout

- `docs/paper1/labbook/`: dated results, failures and next actions; do not rewrite history.
- `docs/paper1/figures/`: the question, inputs and acceptance criteria for each figure.
- `config/paper1/`: site scope and explicit scientific configuration templates.
- `src/paper1/`: small entry points around existing calculations.
- `src/paper1/figures/`: rendering/audit entry points; do not silently refit models.
- `tests/paper1/`: workflow and inventory checks.

The legacy scripts remain authoritative implementations during this setup.
Configuration templates record decisions; the wrappers forward explicit legacy
CLI arguments and do not yet translate YAML into model arguments automatically.
Do not treat a YAML edit as a changed calculation until that connection exists.

## First action on EcoVis

```bash
python src/paper1/audit_sites.py \
  --raw-root data/raw/skyimages \
  --reference-lock data/processed/spider_sign_checked_01/registration/locked_site_references.json \
  --results-root data/processed/spider_prc_lamp_equivalence_01 \
  --output data/processed/paper1/site_audit_01
```

That inventory does not approve registration or Moon identity. It includes the
configured reference sites, additional named RAW folders and any supplied locks.
Missing evidence remains pending. The tracked ledger is the setup snapshot;
actual run inventories are written under processed outputs.

Use `python src/paper1/run_all.py --help` for explicit stage dispatch. It prints
the command without executing unless `--execute` is supplied. It does not
blindly run five unfinished analyses.

## Reusing the pilot

- `reconstruct_prc.py` forwards to the checked camera-band PRC comparison.
- `run_natural_patterns.py` forwards to the exploratory fixed-geometry fitter.
- `fit_pollution_prc.py` forwards to the PRC-equivalent lamp fitter; it requires
  saved PRC targets and a specific spectral input directory.
- Figure 1 forwards to the saved-fit spatial AoP audit.
- Figures 2 and 3 forward to the existing spatial/energy figure renderer.
  These are working adapters, not yet five independent finished publication figures.
- Equal-luminance and intervention analysis entry points return an explicit
  pending status and nonzero exit code until their inputs/methods are implemented.

Each forwarded run gets an adjacent provenance sidecar with exact command,
configuration hashes, source script hash, code version and outcome. Legacy
source checks remain in force. Large data and figures remain outside Git.

## Output convention

```text
data/processed/paper1/<run_id>/
    provenance.json
    site_status.csv
    fig01/ ... fig05/
figures/paper1/<run_id>/
    fig01/ ... fig05/
```

Use fresh output directories. Export figures as PNG, editable PDF/SVG and the
plotted numerical values. A figure script must read existing fits, not refit
pollution or geometry for visual agreement.

Related strategy: [two-paper plan](../lab_notebook/2026-11-23-two-paper-plan.md).
