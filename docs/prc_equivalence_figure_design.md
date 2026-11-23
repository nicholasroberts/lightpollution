# PRC-equivalent spectral scenarios: publication figure drafts

`plot_spider_prc_equivalence_figure.py` creates two complementary figures.

The summary figure uses only the saved equivalence reports and the published
ELF D850 response. Top panels plot exact PRC-versus-added-energy relationships
for each site, with a dashed measured target and fitted HPS/LED intersections.
Lower panels show required energy with HPS=100 and RH1 versus D850 green
photon sensitivities. These are conditional model curves, not independent
experimental dose-response measurements. The shared LED/HPS percentage is
determined by the assumed spectra and receptor sensitivity.

The optional spatial figure uses checked measured Stokes and saved mixed-model
arrays, removes the old fitted background and adds the new PRC-matched spectral
amount. It shows measured camera-green DoLP (the accepted RH1 approximation),
the RH1-weighted prediction, a residual, and the energy-response curve. HPS and
LED RH1 predictions coincide mathematically under uniform unpolarized addition
at matched PRC; the script checks that equality and shows the model once.
It checks source hashes and does not fit, rotate or calibrate analyzer zero.

Spatial panels use only common central-field support (60-degree full field),
with no extrapolation. A shared north-up/east-left sky-view display is applied
to positions only, using each saved display transform followed by an east-west
display reflection. This display choice does not change model handedness or
Stokes conventions. Rasterized scalar maps keep PDF/SVG files compact; text
and response curves remain editable. Absolute AoP is not plotted or claimed
calibrated. The numerical residual is DoLP, and matching integrated PRC does
not force agreement with the local DoLP pattern.

Outputs: PNG at 300 dpi, PDF and SVG for each figure, numerical CSV and
provenance/diagnostic JSON. Do not interpret schematic/synthetic test images
as experimental data. Synthetic pipeline checks exercise the spatial layout;
actual maps must be generated where the original locked Stokes are present.

```bash
cd /run/media/bznwr/data/lightpollution
git pull origin publication-figures

python src/plot_spider_prc_equivalence_figure.py \
  --results data/processed/spider_prc_lamp_equivalence_01 \
  --reference-lock data/processed/spider_sign_checked_01/registration/locked_site_references.json \
  --fit-root data/processed/spider_sign_checked_01 \
  --prc-root data/processed/spider_prc_checked_01 \
  --output data/processed/spider_prc_publication_01

python -m zipfile -c data/processed/spider_prc_publication_01_results.zip \
  data/processed/spider_prc_publication_01
```

Use a fresh output directory. Without the three optional registration/fit/PRC
arguments, only the summary is generated. No new spectral files are selected
or refitted: reported spectral coefficients define the curves and maps.
