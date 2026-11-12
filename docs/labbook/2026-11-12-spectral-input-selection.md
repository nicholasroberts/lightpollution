# Spectral inputs and next model stage — 12 November 2026

## Confirmed decisions

Cottage capture date is 21 August 2018; September in the reference table is a typo, confirmed by investigator. The reference table is corrected and its date_source records this decision. Existing EXIF-based lunar positions already used August, so this does not change that track. Old anchor files may still carry a superseded date warning; do not edit original saved results silently.

For now model a horizontal body, yaw 0–360 degrees, no tilt sweep. Default --tilts is now 0. Existing assumption: horizontal animal corresponds to zenith-centred field; anatomical offset is not newly measured. All prior fixed shared colour-scale settings retained.

Species spectrum unknown: use a declared wolf-spider-inspired proxy, alpha peak 510nm and secondary UV band centred at365nm (sensitivity scenarios 360–370nm). DeVoe1972 doi10.1085/jgp.59.3.247 / PMID5058960 reports principal-eye peak510nm plus360–370nm or shorter and very variable UV/visible ratios. Secondary anterior-lateral cells peak510nm with shoulder near380nm. Therefore the requested365nm band is not evidence for the target species or posterior median eyes. Do not assume a fixed UV amplitude or interpret it as automatically a beta band of a single pigment. Compare UV-band centres365 and380nm and UV strengths independently.

Stavenga2010 doi10.1007/s00359-010-0560-z (https://pmc.ncbi.nlm.nih.gov/articles/PMC2962788/) compares invertebrate pigment templates. Alpha templates agree well above400nm; below400nm beta-band predictions diverge. A510nm alpha template plus a separately declared UV-band sensitivity scenario is preferable to claiming a validated whole-eye action spectrum. Distinguish absorbance from sensitivity to avoid applying the same pigment weighting twice. Use consistent photon units.

## Primary spectral sources found

* NIST Lykke et al2013, Precise Measurement of Lunar Spectral Irradiance at Visible Wavelengths, doi10.6028/jres.118.020: https://www.nist.gov/publications/precise-measurement-lunar-spectral-irradiance-visible-wavelengths . Calibrated near-full lunar irradiance with visible coverage420–1000nm; source-level shape, not diffuse-sky radiance, and lacks365nm.
* USGS Kieffer & Stone2005, The spectral irradiance of the Moon: https://www.usgs.gov/publications/spectral-irradiance-moon . ROLO observations from northern-hemisphere Flagstaff cover350–2450nm in32bands and provide lunar geometry-dependent spectral source. Interpolation and access to numeric model/data still to establish; not downloaded or ingested yet.
* NIST Woodward et al2022 air-LUSI doi10.1088/1681-7575/ac64dc: https://www.nist.gov/publications/measurements-absolute-si-traceable-lunar-irradiance-airborne-lunar-spectral-irradiance . Northern-hemisphere campaign,350–980nm coverage, phases10–60degrees. Useful UV-inclusive source reference, not an exact-full-moon ground-sky measurement. Do not conflate range of coverage with1percent uncertainty below450nm.
* Bará, Aubé, Barentine & Zamorano2020, Magnitude to luminance conversions and visual brightness of the night sky, doi10.1093/mnras/staa323, https://academic.oup.com/mnras/article/493/2/2429/5722125 .19measured calibrated sky spectra plus separate lamp templates. Select HPS-dominated and LED-dominated cases only after inspecting individual site descriptions. Lamp comparison panel is not pure measured skyglow.
* Emission spectra of light-pollution sources determined from the light-scattering spectrometry of the night sky (2020), doi10.1093/mnras/stz3260, https://academic.oup.com/mnras/article/491/4/5586/5686816 . Night-sky measurements and generic HPS/LED component spectra, including multiple LED CCT groups. Retrieve numerical components or digitize documented figures; distinguish inferred emission from observed skyglow.
* Kolláth et al2020, The Colour of the Night Sky, doi10.3390/jimaging6090090, https://www.mdpi.com/2313-433X/6/9/90 . Natural/artifical spectral decomposition and measured sky examples. Useful cross-check.
* Kocifaj et al2013, Light Pollution in Ultraviolet and Visible Spectrum: Effect on Different Visual Perceptions, https://pmc.ncbi.nlm.nih.gov/articles/PMC3575508/ . UV-inclusive modelled skyglow comparison including HPS and LED; assess as UV scenario support where measured tables truncate at380/400nm. Never fill missing UV with zero without evidence.

These are located literature sources, not yet validated numerical arrays. No spectral rerun has been made.

## Lunar brightness recommendation

First implementation: Krisciunas & Schaefer1991 doi10.1086/132921, https://ntrs.nasa.gov/citations/19910069763 . It calculates scattered moonlight using phase, Moon zenith angle, viewing zenith angle, Moon–view separation and extinction. Suitable transparent relative-intensity pilot. Native empirical V-band model, not a UV/visible spectral model: multiplying a fixed spectrum by its brightness is only an interim separable approximation. Use ephemeris phase and distance as well as position. Add independently declared unpolarized natural background rather than extrapolating a polarized Moon below horizon.

Preferred spectral extension: Jones et al2013, An advanced scattered moonlight model for Cerro Paranal, doi10.1051/0004-6361/201322433, https://www.aanda.org/articles/aa/full_html/2013/12/aa22433-13/aa22433-13.html . Spectral higher-order scattering with geometry and atmospheric extinction; optical evaluation0.36–0.89micron, extended framework0.3–2.5micron. Developed for Paranal, so replace site-specific atmosphere/altitude with relevant northern lowland scenarios. No need to treat the southern-hemisphere source location as a different lunar spectrum; local atmospheric transfer and geometry are what must be changed.

Keep attenuation/scattering used to predict I separate from empirical depolarization of Q/U, but do not double count atmospheric depolarization when combining models. Hold absolute ALAN fixed during each nightly scenario so ALAN/natural varies with lunar brightness. Compare normalization at equal integrated radiance first; equal VIIRS radiance requires detector spectral-response conversion, not equal lux.

## Execution status

GitHub reference date and horizontal default updated. This turn's local compute environment became unavailable; those small edits were not executed locally. Prior generated arrays unchanged. No pigment arrays, measured spectra or lunar-brightness code claimed complete.
