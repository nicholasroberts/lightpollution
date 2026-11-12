# Spectral and nightly brightness pilot — 12 November 2026

## Implemented

* Horizontal body, yaw0–350deg in10deg increments,60deg full field. Uniform solid-angle cap and transported eye axes retained. Cottage date21August2018 confirmed; reference table now agrees with EXIF.
* `spider_spectral_lunar.py`: Stavenga alpha pigment template at510nm plus separately declared Gaussian UV band centred365nm. Curves at365/380nm and amplitudes0,0.29,1 are generated. Width40nm is a scenario, and0.29is an illustrative strength, not fitted to DeVoe. A secondary UV pigment and an intrinsic beta band are not claimed equivalent. Correct template comparison DOI:10.1007/s00359-010-0568-7.
* Photon sensitivity is converted to energy weighting using wavelength before the receptor integrations; no second pigment weight is applied. Polarization absorbances remain0.4/0.04, wavelength independent pending species data. The pigment is a declared effective-sensitivity proxy, not a full self-screening or photochemical model.
* `model_spider_nightly_prc.py` now defaults to `--brightness-model ks91`; `constant` preserves the old intensity assumption. KS1991 uses phase, geometric airmasses, Moon elevation, Moon–view separation, V-band extinction, plus an explicit inverse-square distance correction relative to384400km. Phases and distance are calculated from builtin Astropy ephemerides. `--extinction-v 0.2`is an assumption in mag/airmass; not a fitted atmospheric quantity. Polarization scaling factors remain separate sensitivity parameters.
* Natural brightness is computed for every viewing ray, so integration is intensity weighted. Each site's unit intensity is its predicted capture-time zenith moonlight. Input ALAN ratios multiply that fixed reference; instantaneous ALAN/natural ratios change through the night and are recorded per response. This is NOT equal absolute ALAN across sites. No artificial-light switch-offs, clouds or terrain masks.
* KS91 excludes separations<10deg, and the pilot excludes Sun above-6deg, Moon below0.1deg or above89.9deg. Below-horizon Moon is excluded, not interpreted as zero total natural sky or zero available cue. No separate moonless background. V-band spatial scaling with fixed spectral shape remains a separable approximation rather than full spectral radiative transfer.

## Numerical source retrieval

The author's numerical archive is linked at https://davinci.fmph.uniba.sk/~kundracik1/SpectrumAnalyser/ and associated with Kocifaj,Kundracik&Bily2020 doi10.1093/mnras/stz3260. `prepare_spider_spectral_pilot.py` downloads the linked archive or reads `--archive` and records its SHA256 and member names. It never executes archive binaries. It exports HPS and three generic LED CCT groups, preserving author arrays. These are SOURCE EMISSION spectra, not pure measured skyglow. It also exports measured mixed sky spectra from Bratislava2010 and Vienna2012/2018; negative instrumental/noise samples are retained in these examples, not silently clipped or fed as physical radiance. The archive includes a wider original300–900nm grid; tiny UV lamp signals are not assumed validated emissions.

NIST2013 doi10.6028/jres.118.020 Table1 is transcribed as9tabulated wavelength/irradiance/uncertainty samples from northern-hemisphere near-full Moon measurement2012-11-30. This is TOA source irradiance, not diffuse sky radiance. The tabulated visible samples permit a450–750nm pilot using coarse linear interpolation, with no out-of-range extrapolation. No UV-inclusive natural spectrum is available in this input. Consequently, UV proxy curves are supplied but no UV contribution to PRC is claimed.

For this exploratory visible pilot, natural irradiance is multiplied by(lambda/550nm)^-4 and generic ALAN source spectra by(lambda/550nm)^-2. These are declared scattering-shape scenarios, not a fitted transfer or a substitute for Jones2013. ALAN exponent can be varied0–4. Spectra are normalized to equal450–750nm band-integrated energy before ratios are applied. This is neither equal lux nor equal VIIRS radiance, and no local lamp mix has been inferred.

## Reproduce

```bash
python src/prepare_spider_spectral_pilot.py --output data/processed/spider_spectral_inputs_01
python src/model_spider_nightly_prc.py \
 --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
 --spectra data/processed/spider_spectral_inputs_01/HPS_visible_input.csv \
 --output data/processed/spider_nightly_hps_visible_01
python src/model_spider_nightly_prc.py \
 --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
 --spectra data/processed/spider_spectral_inputs_01/LED_mid_visible_input.csv \
 --output data/processed/spider_nightly_led_visible_01
python src/plot_spider_spectral_pilot.py \
 --spectra data/processed/spider_spectral_inputs_01 \
 --hps data/processed/spider_nightly_hps_visible_01 \
 --led data/processed/spider_nightly_led_visible_01 \
 --output data/processed/spider_spectral_comparison_01
```

## Verification and next boundary

34numerical tests pass, including an independent KS91 implementation's19.855mag/arcsec2example, pigment peak, inverse-square distance, phase dimming, constant-ALAN response scaling and existing geometry checks. Both sites complete at20min/10deg/1024rays:14256responses per spectrum. Input spectra and scenario provenance are bundled with the generated plots, not treated as site observations.

Still required: a defensible UV-inclusive lunar spectral source and local spectral transfer; whole-eye proxy/UV strength sensitivity; independent intensity calibration or explicitly relative ALAN scenarios; realistic moonless background before interpreting the whole night. Jones et al2013 doi10.1051/0004-6361/201322433 is the preferred spectral brightness extension. Atmospheric extinction predicts intensity and wavelength redistribution; phenomenological Q/U depolarization must not duplicate scattering already represented in a coupled model. No biological threshold applied.


## Capture-time sky DoLP sense check

The comparison plotting script now adds six calculated sky panels: natural-only, natural+HPS and natural+LED for each site. These use exact capture-time lunar geometry and KS91 relative sky intensity from the saved run metadata, not the nearest20minute sample. DoLP is normalized after integrating spectral Stokes with the photon-sensitivity weighting in the supplied spectral input. It is labelled effective sky DoLP450–750nm, not raw monochromatic or radiometric-camera DoLP. All panels share0–1limits; white circle30degree radius, crossMoon, northup/eastright, grayKS91circumlunar exclusion. The natural-only column retains the chosen atmospheric factor0.8but has zeroALAN. HPS/LED columns use the same settings as the selected PRC curves.

Defaults are --atmospheric-factor0.8 and --alan-ratio1.0. This ALAN is unpolarized, uniform over sky and constant in time. Its450–750nm energy equals each site's predicted natural zenith intensity at that site's capture time, so this is not equal absolute ALAN between Cottage and Oxford. HPS and LED have equal energy in the specified band within each site; different spectra change eye-weighted dilution. No amount is fitted to observed DoLP orVIIRS. Optional CLIsettings select scenarios already present in saved response tables. No raw-image processing or nightly response recomputation is needed; rerun only the plotter to a new output directory. Snapshot settings are written to sky_snapshot_settings.json.

New numerical check: identical natural/artificial spectra and ALANratio1halve effectiveDoLP in constant-intensity mode; zeroALANpreserves naturalDoLP; atmospheric factor bounds amplitude. Six spectral/lunar tests pass (35combined prior tests including the new check).
