# Jones/ESO scattered-moonlight comparison — 27 November 2026

Use the actual published Jones implementation exposed by ESO SkyCalc, not a new
power-law approximation or a handwritten partial implementation.

## Keep the first calculation simple

- One zenith direction per site, using the final Figure 1 capture report.
- Moon altitude from the accepted astronomical anchor; Moon–zenith separation
  is 90 degrees minus that altitude.
- Geocentric Moon–Sun elongation and relative Earth–Moon distance calculated at
  the recorded UTC using offline Astropy ephemerides.
- Fixed ESO Paranal reference atmosphere and PWV=3.5 mm; no fitted atmosphere.
- 400–700 nm, 1 nm sampling, only the scattered-moonlight component enabled.
- No new ALAN fit, no Figure 1 change, no solar/twilight spectrum substitution.

This is a reference-atmosphere spectral comparison at site capture geometry,
not a fully adapted lowland-site prediction. ESO's public service does not expose
arbitrary lowland aerosol/pressure profiles. Do not describe its stock Paranal
atmosphere as the actual Swedish/Oxford atmosphere. Retain that explicit scope
rather than adding unconstrained atmospheric parameters.

## Implementation

src/prepare_paper1_jones_moon_spectra.py calls the official skycalc_cli>=1.6.
It saves all submitted parameters and the returned FITS unchanged. Extract FLUX_SML,
not total FLUX, and verify wavelength and radiance units. The documented output
is photons/s/m2/micron/arcsec2 on a wavelength grid in nm. Convert once using

    E_lambda = photons * h*c/(lambda_nm*1e-9) / 1000 * (180*3600/pi)^2

to W/m2/nm/sr. Then normalize to unit 400–700 nm integrated energy for the
comparison plot. Do not multiply by another lambda^-4 or TRANS: atmospheric
propagation is already included. Preserve physical output alongside normalized
arrays so normalizing the curves does not imply absolute spectral validation.

Jones/ESO uses its own solar-spectrum/lunar-reflectance prescription. The audited
lunar source CSV from our wider input preparation is used only to display
source-context and the prior source-times-lambda^-4 approximation. It is not
injected into SkyCalc or claimed to be its source spectrum.

Outputs: one CSV per site with photon, physical energy and normalized energy
spectra; input JSON; original FITS; combined PNG/PDF/SVG comparison and provenance.
Do not reuse either lunar spectrum for the museum's solar component.

## Run on EcoVis

Use the wider spectral source directory prepared for the 400–700 nm update.

```bash
git pull --ff-only origin paper-1-polarization-model
python -m pip install -r requirements-spectral.txt

python src/prepare_paper1_jones_moon_spectra.py \
  --spectra-root data/processed/paper1_spectral_inputs_02 \
  --output data/processed/paper1_jones_moon_spectra_01 \
  --figure-copy figures/paper1/fig02a/jones_check_v01
```

A fresh output directory is required. Each site makes one online ESO request.
The installed model is not fitted to move the spectral maximum. A blue maximum
is retained if returned by the model.

## Validation

Three synthetic tests pass: energy/photon conversion and energy normalization;
explicit extraction of FLUX_SML with unit checks rather than total FLUX; zenith
geometry and fixed settings. One real official SkyCalc 1.6 request also completed
for Cottage's accepted UTC and Moon altitude. The returned FITS confirms nm
wavelengths and the documented photon units. The sampled normalized-energy peak
for that Paranal-reference/zenith case is 412 nm, so this reference calculation
does not by itself establish agreement with the investigator's longer-wavelength
measurements. No actual measured spectrum was provided for quantitative comparison.
The two-site production comparison and archive run on EcoVis.

## Sources

- Jones et al. 2013, DOI 10.1051/0004-6361/201322433:
  https://arxiv.org/pdf/1310.7030
- Official ESO model, source and site-dependence notice:
  https://www.eso.org/sci/software/pipelines/skytools/skycalc
- Official CLI parameters/output units:
  https://eso.org/observing/etc/doc/skycalc/helpskycalccli.html
