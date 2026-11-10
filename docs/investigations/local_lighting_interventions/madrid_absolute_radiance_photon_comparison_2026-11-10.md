# Madrid: SI radiance and approximate photons, not percentages alone (10 Nov 2026)

An explicit test of the investigator's point: a large percentage change in one spectral band can represent a small absolute photon-radiance difference compared with a small percentage change in a strong broadband detector. **These are still different observational geometries, NOT a conservation-of-photons balance or a directly measurable net change in emissions.**

## Actual input

- **Ground:** Robles et al. (2021), [Remote Sensing 13, 1511](https://doi.org/10.3390/rs13081511), Tables 3–5, quality-filtered *winter* 01:30 AstMon zenith night-sky brightness in Johnson B/V/R. Winter year assignments may span the year boundary.
- **Satellite:** user's real EOG monthly \`VCMSLCFG\` VIIRS VNL V1 file \`data/processed/local_lighting_interventions/madrid_ground_viirs/local_ucm/madrid_ucm_3km_10km_2015_2019_monthly.csv\` from the previously defined **3 km/10 km circles** centred on UCM, \`cf_cvg>=3\`, all pixels present. Calculate equal-weight means of January, November and December separately within 2015 and 2019. These calendar-month periods are a seasonal *proxy* rather than exactly the same winters/nights as the ground paper.
- Read-only workflow: [\`src/plot_madrid_absolute_photon_radiance_2015_2019.py\`](../../../src/plot_madrid_absolute_photon_radiance_2015_2019.py); **does not change** any frozen European NUTS3/VIIRS/PRC or original local raster.

## Calculated results, absolute radiance

All changes within each instrument are translated to the same **SI dimensional unit** (\`µW m⁻² sr⁻¹\`). That does *not* make ground zenith sky radiance equivalent to upward DNB surface radiance.

| Observational series | 2015 µW m⁻² sr⁻¹ | 2019 µW m⁻² sr⁻¹ | Δ µW m⁻² sr⁻¹ | % |
|---|---:|---:|---:|---:|
| Ground UCM Johnson B | 4.553 | 6.170 | **+1.617** | **+35.52** |
| Ground UCM Johnson V | 7.626 | 8.439 | **+0.813** | **+10.66** |
| Ground UCM Johnson R (approximate Johnson-R calibration) | 17.733 | 17.250 | **−0.483** | **−2.73** |
| Orbital VIIRS 3 km (Jan/Nov/Dec calendar mean) | 508.259 | 478.171 | **−30.088** | **−5.92** |
| Orbital VIIRS 10 km | 457.937 | 422.605 | **−35.332** | **−7.72** |

The native VIIRS monthly source is in \`nW cm⁻² sr⁻¹\`; numerically multiply by **10** for \`µW m⁻² sr⁻¹\`. This is a simple conversion of physical **unit dimensions**, not an assertion that DNB acts as a perfect square broadband spectral integrator. NASA/NOAA VIIRS DNB is response-weighted over roughly 500–900 nm and gives NO separate wavelength spectrum ([NOAA FAQ](https://ncc.nesdis.noaa.gov/VIIRS/VIIRSFAQ.php); [NOAA DNB table](https://www.star.nesdis.noaa.gov/jpss/VIIRS.php)).

## Calibrations and photon conversion (approximate)

\`L_lambda = F_lambda,0 × 10^(−0.4m) / [(pi/(180×3600))²]\`. Multiply by nominal photometric effective band width to get an approximate band radiance. Here the adopted zero-point spectral flux density is in **W m⁻² nm⁻¹**, and magnitudes are **mag arcsec⁻²**.

- B: zero \`6.31e-11 W m⁻² nm⁻¹\`, λ \`436nm\`, nominal width \`89nm\` (Vega, [ESO](https://www.hq.eso.org/sci/facilities/eelt/science/drm/tech_data/photdefs/)).
- V: \`3.607e-11\`, λ \`545nm\`, width \`84nm\` ([Oxford calibrations](https://academic.oup.com/view-large/21163628)).
- Johnson **R**: \`1.76e-11\`, λ \`700nm\`, width \`220nm\` (original Johnson R, [Johnson-era calibration](https://dmaitra.webspace.wheatoncollege.edu/conversion.html), [Starlink bandwidth](https://starlink.eao.hawaii.edu/docs/sc6.htx/sc6se7.html)). **Johnson R ≠ Cousins R**: a specific instrument transmission curve and standard-star calibration are needed for accurate absolute R energy.

Estimated photon radiance conversion (per band, using one representative wavelength):

\`N_band = L_band × λ_eff/(h c)\`.

| Series | Change in \`10¹² photons m⁻² s⁻¹ sr⁻¹\` |
|---|---:|
| Ground B, 436 nm | **+3.55** |
| Ground V, 545 nm | **+2.23** |
| Ground R, nominal 700 nm | **−1.70** |
| VIIRS 3 km, **assumed 700nm equivalent** | **−106.0** |
| VIIRS 10 km, **assumed 700nm equivalent** | **−124.5** |

**For VIIRS, these are ILLUSTRATIVE photon equivalents of the instrument's reported energy radiance, NOT measured total photon counts or a spectral inversion.** Assuming a single wavelength anywhere within 500–900 nm spans **−75.7 to −136.3** (3 km) and **−88.9 to −160.1** (10 km) in the same \`10¹²\` unit; whiskers in figure are *not* uncertainty bars. A real VIIRS in-band photon radiance estimate needs the spectral response function and source spectral energy distribution.

**Key caution:** The ground instrument looks at the **zenith sky radiance after atmospheric scattering** from a single site. VIIRS receives **upward radiance** from areas including direct and reflected light, under different atmospheric geometry. Even with dimensional SI conversion they are NOT two parts of one energy/photon budget. No subtracting, no mass balance, no inference that ~30× more photons were extinguished in the environment from the apparent magnitudes. B/V/R passbands overlap; their photon counts **cannot be naively added** to form a total sky flux.

## Interpretation relevant to the original PRC model

The statement **“35.5% ground B increase vs 5.9–7.7% VIIRS fall”** looks dramatically asymmetric, but the measured SI radiance changes have very different baselines. The **numerical** VIIRS energy-equivalent drop, **30–35 µW m⁻² sr⁻¹**, is ~19–22 times the **B-band** sky radiance increase **1.62 µW m⁻² sr⁻¹**. The **illustrative** photon-equivalent magnitudes differ by ~30–35 times (at assumed 700nm for VIIRS). However, comparing these numbers does **not** establish whether the *total* number of atmospheric photons increased or decreased (opposite view geometries, spectra, source footprints). The meaningful finding remains that the ground sky became bluer/brighter in B despite satellite response going the opposite way.

Biologically, receptor-relevant change requires a **sky-radiance spectrum**, angular source/skyglow scattering and, for polarized receptor cues, polarized fraction and angle, with receptor spectral/polarization sensitivity. VIIRS cannot resolve these directly. Unlike some earlier qualitative interpretations, it is **not defensible to say the Madrid % trends represent similarly sized absolute photon changes**, even as a net emission balance.

## Run (EcoVis)

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/plot_madrid_absolute_photon_radiance_2015_2019.py
\`\`\`

Output figure \`figures/local_lighting_interventions/madrid_ground_viirs/local_ucm/madrid_absolute_radiance_photons_2015_2019.{png,pdf,svg}\`; CSV \`data/processed/local_lighting_interventions/madrid_ground_viirs/local_ucm/madrid_2015_2019_absolute_radiance_photon_comparison.csv\`.

This is a **re-plot of actual Madrid observed values using explicit approximate calibrations**, not a recalculation of PRC, VIIRS rasters, or the historical 2014–2024 European dataset.
