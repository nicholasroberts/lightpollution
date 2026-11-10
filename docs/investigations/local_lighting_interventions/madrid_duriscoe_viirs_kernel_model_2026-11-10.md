# Madrid UCM: fixed VIIRS-to-skyglow distance kernel (10 November 2026)

**Scientific status:** **Model implementation is committed to the user's \`publication-figures\` branch; not yet run on the user's actual 300 km monthly VIIRS rasters in EcoVis.** The independent 2015–2019 monthly *3/10 km* VIIRS and published ground brightness results were already observed in earlier work. Do **not** invent the new 25–300 km model outcome.

## Question

Would a spatially explicit, *fixed* propagation relationship applied to original 2015–2019 EOG VNL V1 upward radiance around the UCM observatory reproduce the **direction** of the UCM AstMon Johnson B sky brightening? This isolates the influence of spatial source redistribution under a stationary transfer kernel; it does NOT isolate LED spectral bias or prove an atmospheric cause.

Original directly observed partial comparison:

| Observation | Seasonal 2015→2019 change |
| --- | ---: |
| 3 km UCM VIIRS mean, Jan+Nov+Dec within calendar years | **−5.9%** |
| 10 km UCM VIIRS mean, same months | **−7.7%** |
| UCM Johnson B zenith, published filtered winter 01:30 | **+35.5%** |
| UCM Johnson V, published filtered winter 01:30 | **+10.7%** |
| UCM Johnson R, published filtered winter 01:30 | **−2.7%** |

**These are distinct instruments, wavelengths, viewing directions and seasonal year conventions.**

## Published model closely examined

[Duriscoe, Anderson, Luginbuhl & Baugh (2018), *A simplified model of all-sky artificial sky glow derived from VIIRS Day/Night band data*, JQSRT 214:133–145](https://doi.org/10.1016/j.jqsrt.2018.04.028). [Original NOAA hosted open PDF](https://repository.library.noaa.gov/view/noaa/21007/noaa_21007_DS1.pdf).

The published **Eq. 9** is a sum of VIIRS upward radiance \`r_i\` weighted by a distance-dependent function \`d_i^(-alpha)\` out to **300 km**, followed by a calibration constant \`c\` to predict **ALR = artificial average luminance over the entire sky / natural average sky luminance**. The preferred **Eq. 11** is:

\`\`\`text
alpha(d) = 2.3 * (d/350)**0.28     # d in km
W = sum_i r_i * d_i**(-alpha(d_i))
\`\`\`

The original remote-site calibration **562.72** was based on their historical VIIRS product, specific equal-area raster resolution, **all-sky V-band photometry** and **cloudless, low-aerosol conditions** of predominantly *non-urban* field sites in the southwestern USA. It **must not be copied over** to predict absolute **Johnson B/V/R zenith** brightness *inside Madrid*. The authors explicitly describe the breakdown of their extended source/city-distance relationship **for observer–source distances below 20 km** (paper §2.4). They distinguish their all-sky ALR from zenith ZLR and warn that the method was developed for sites away from major sources; see paper §1, §3.2, §4. We therefore implement **an uncalibrated source-weighted VIIRS index**, not a fully physically validated skyglow prediction.

## Reproducible first-pass transfer experiment

**New code:** [\`src/madrid_ucm_duriscoe_skyglow_kernel.py\`](../../../src/madrid_ucm_duriscoe_skyglow_kernel.py).

- Centre **UCM observatory 40.451111°N, −3.726111°E**, IAU I86, precisely matching previous 3 km and 10 km tests.
- Use the *original already-imported* EOG VNL V1 \`VCMSLCFG\` \`75N060W\` \`avg_rade9h\` and \`cf_cvg\` GeoTIFFs, each year **2015–2019**, default **January, November, December** (15 observations). Full annual calendar or March–October can be requested with \`--months\`.
- Use **one fixed 300 km radius** around UCM; calculate the pixel-centre distances in **EPSG:3035**, consistent with existing UCM circle geometry, and correct original native lon/lat pixel area approximately by \`cos(latitude)/cos(latitude_UCM)\` (the published method resamples to equal-area grid). This last factor cancels in no artificial time-dependence.
- For all pixels with \`cf_cvg>=3\`, clip **original radiance values below 0.5 nW cm⁻² sr⁻¹ TO ZERO** (Duriscoe's original low-background threshold). **Do not subtract 0.5 from every pixel.** Negative original values become zero. All source masks and the 2015→2019 spatial distances are fixed.
- Weight and sum each month over **3, 10, 25, 50, 100, 200, 300 km** circles and **0–3, 3–10, 10–25, 25–50, 50–100, 100–200, 200–300 km** rings.
- Prevent uncontrolled near-field divergences or pretending validation at under 20 km: use a **0.5 km centre-distance floor**, with sensitivity to **0.25, 0.5, 1, 2, 5, 10, 20 km**. Do **not** interpret the floor as an independently measured atmospheric mean free path or physical scattering coefficient.
- Check the **identical spatial support** of the 2015 and 2019 observations for *each calendar month*. Derive sensitivity results using only **cells valid in both years**, record minimum valid fraction, and vary source radius/floor **without selecting best outcome post hoc**. Annual curves use monthly valid radiances separately and may have different supports; compare paired results as the stronger diagnostic. Cloud counts do not guarantee equal aerosols or source spectral mixture.
- Plot (A) each radius's weighted 2015-indexed history, (B) 2015→19 matched-month percent change across distance radius and near-field floor with the observed **Johnson B +35.5%** benchmark, (C) which distance annuli contribute positively/negatively, (D) the 300 km *fixed-kernel model's* 2019/2015 ratio versus *actual* Johnson B/V/R filtered winter measurements.
- The model does **not** estimate flux/sky brightness in absolute SI, scattering fraction, receptor photon catch, a 500–900 nm emission spectrum, or polarization contrast. Predicted ratio remains a **relative VIIRS-based ALR-like source index**, not a measurement of Johnson B.

**Important hypotheses:**

1. **If all radii and near-field assumptions predict a decline while B observed rises**, then spatial redistribution of DNB-measured upward radiance, *under this stationary kernel*, fails to explain the observed B increase. This is compatible with a changed SED/sensor response, optical directionality, aerosol conditions, non-VIIRS blue-band sources, or other phenomena; it does NOT uniquely prove LED spectral bias.
2. **If source radii >10 km predict an increase while inner radii decline**, then surrounding communities and spatial propagation explain at least the sign under this kernel. Investigate which annuli/municipalities contribute, rather than blaming instrument spectral response.
3. **If the sign changes with near-field floor**, the urban-core application is unstable; *no strong skyglow attribution* is warranted. Further work should use a physically based zenith spectral Garstang/Falchi transfer model including emission angular distribution.
4. **If >200 km radius dominates**, audit remote EOG geospatial coverage/cloud conditions and boundary effects before interpreting.

## Run on EcoVis

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures

python src/madrid_ucm_duriscoe_skyglow_kernel.py --self-test
python src/madrid_ucm_duriscoe_skyglow_kernel.py --inventory
python src/madrid_ucm_duriscoe_skyglow_kernel.py

# Optional same instrument, different seasonal sensitivity:
python src/madrid_ucm_duriscoe_skyglow_kernel.py --months 3 4 5 6 7 8 9 10
\`\`\`

No new raw EOG downloads should be required if previous local-month workflow imported the 2015–2019 VCMSLCFG Europe tile and \`--inventory\` shows 15/15 months.

Outputs:
- \`figures/local_lighting_interventions/madrid_ground_viirs/duriscoe_kernel/madrid_ucm_duriscoe_kernel_11112.{png,pdf,svg}\`
- \`data/processed/local_lighting_interventions/madrid_ground_viirs/duriscoe_kernel/madrid_ucm_duriscoe_{monthly,annual,sensitivity,rings}.csv\`
- \`data/processed/local_lighting_interventions/madrid_ground_viirs/duriscoe_kernel/model_audit.json\`

**Validation status:** The full algorithm and rendering route were exercised in a separate current-session synthetic raster prototype, testing additive ring weighting, finite near-field behaviour, exactly matched 2015/19 dates, 300 km masks, 15 calendar-month composites, SVG/PDF/PNG output and the physical requirement that DNB radiance under 0.5 is **set to zero rather than subtracting an offset**. The committed code's compact routine also has a \`--self-test\`, but **the user's real 15-month 300 km model remains unrun** and no model result has been assumed. The existing frozen 2014–24 annual European radiance, snow adjustment, VIIRS zero correction and PRC must not be modified for this local test.
