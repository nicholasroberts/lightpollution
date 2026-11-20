# PRC-equivalent HPS and LED spectral scenarios

The investigator accepts reconstructed camera-green PRC as an approximate RH1
target. This is an explicit approximation, not a uniquely recovered animal
spectrum. Use the new checked reconstructions for WClawn and Oxford Business
Park, not historical reference PRCs.

`fit_spider_prc_lamp_equivalence.py` retains the existing mixed fit's Sun/Moon
geometry, polarization and source ratio. It removes the old background and
fits a uniform, unpolarized HPS or LED amount directly to maximum absolute PRC.
No DoLP objective is used for the lamp amount. The source ratio remains the
previous DoLP-derived ratio; this script does not implement a new joint AoP fit.

The assumed lunar and solar spectral template is the same transferred natural
spectrum from the existing spectral pilot. Under this approximation their
camera-band ratio equals their band-energy ratio. Solar twilight colour is
not independently retrieved. HPS and mid-CCT LED spectra are normalized to unit
radiant energy over the same shared band (currently 352–700 nm). The photon
response is wavelength times the receptor's photon sensitivity; hc cancels.
The prior pilot sky-transfer assumptions are retained in its input spectra.

Natural spectral I, Q and U are integrated through the receptor, then spatially
integrated on the exact preceding PRC support. Uniform added light contributes
only to I. With the same absorbances and receptor calculation as the preceding
PRC comparison, solve analytically for the amount at the common target PRC.
An unreachable target is reported, not forced by a negative light amount.
The comparison reports LED/HPS required energy and 100*(1-LED/HPS).
Amounts refer to added sky radiance integrated over the stated wavelength band,
relative to lunar zenith band radiance. They are not electrical power or lux.

The primary run uses the 510 nm pigment proxy. The sensitivity run uses the
published D850 green response, holding the same numerical target fixed. This
assesses the response approximation, not a unique correction of measured data.
For uniform unpolarized addition, the HPS/LED amount ratio is determined by
their receptor spectral overlaps; site geometry and target determine the
absolute required amounts. It should not be interpreted as a site lamp-mix
retrieval or a site-specific measured LED penalty.

## D850 calibration data provenance

`data/reference/d850_elf_response.csv` is a 1 nm sampling of the numerical
RGB curves in ScienceDjinn ELF, `calibration/nikon d850/rgb_calib.mat`.
Original Git blob: `f7227a4e7eff5af7d4d469b76ae15d6291465239`.
The MAT file's r/g/b arrays are tabulated at 0.5 nm from 350 to 700 nm;
every second row is retained without fitting, smoothing or extrapolation.
Values retain their original photon-response scale; spectral calculations
normalize the green curve to its peak because absolute scale cancels.

Sources and attribution:

- Jochen Smolka, ScienceDjinn, Nikon D850 camera calibration report, 24 May 2022:
  https://www.sciencedjinn.com/wp-content/uploads/2022/05/D850-Calibration-Report.pdf
- ELF numerical calibration:
  https://github.com/sciencedjinn/elf/blob/master/calibration/nikon%20d850/rgb_calib.mat
- Nilsson, Smolka & Bok (2022), *The vertical light-gradient and its potential
  impact on animal distribution and behavior*, Front. Ecol. Evol. 10:951328,
  https://doi.org/10.3389/fevo.2022.951328 ; supplementary Figure S1.

The calibration uses the Nikon D850 and Sigma 8 mm f/3.5 fisheye. The added
polarizer used in this project is not part of that calibration. These curves
are used as published system-response approximations; the current RAW
registration, black subtraction and Stokes reconstruction are unchanged.

## Run on EcoVis

```bash
cd /run/media/bznwr/data/lightpollution
git pull origin publication-figures

for site in WClawn oxfordbussinesspark; do
  python src/fit_spider_prc_lamp_equivalence.py \
    --fit "data/processed/spider_sign_checked_01/$site/fit" \
    --prc "data/processed/spider_prc_checked_01/$site" \
    --spectra-root data/processed \
    --output "data/processed/spider_prc_lamp_equivalence_01/$site"
done

python -m zipfile -c data/processed/spider_prc_lamp_equivalence_01_results.zip \
  data/processed/spider_prc_lamp_equivalence_01
```

The script searches beneath `--spectra-root` for `HPS_visible_input.csv` and
`LED_mid_visible_input.csv`. Identical duplicates are allowed; conflicting
versions require a narrower directory containing the intended spectral pilot.
No source spectra are downloaded automatically. Use a fresh output directory.
