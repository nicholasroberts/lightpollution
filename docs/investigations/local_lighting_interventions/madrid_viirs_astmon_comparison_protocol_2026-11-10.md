# Madrid comparison: locked NUTS3 annual VIIRS, city monthly VIIRS, and UCM AstMon observations

**Prepared 10 November 2026, branch \`publication-figures\`.** The plotting code and method are committed; **no real Madrid EcoVis annual or monthly VIIRS output has been run/reviewed in this conversation**. Do not report an observed Madrid VIIRS direction or magnitude until that test. This is an independent methodological comparison and **does not alter frozen 2014–2024 European radiance, 2021 snow, 2017 offset, NUTS3 geography or spider/moth PRC**.

## Objective

Evaluate whether the **2015–2016 HPS→LED transformation in Madrid** coincides with a decline, increase or no change in (i) our frozen NUTS3 ES300 annual upward night radiance, (ii) the fixed municipality's monthly EOG VIIRS VNL V1 radiance, and (iii) **independent ground-based zenith sky brightness and colour** reported by [Robles, Zamorano, Pascual, Sánchez de Miguel, Gallego & Gaston (2021), *Evolution of Brightness and Color of the Night Sky in Madrid*, Remote Sensing 13(8):1511](https://doi.org/10.3390/rs13081511).

**Especially relevant original paper results:**
- AstMon UCM full Johnson B, V, R and colour measurements **2010–2020**, with a strong HPS-to-LED expansion **during 2015–16**, continuing older and newer changes. ISS photographs **2011 vs 2017**, SAND spectra for **January 2014 vs January 2019**, show colour changes.
- The authors found after-midnight Johnson B brightness **~32% brighter 2015→2019** after ~14% pre-2015 darkening; Johnson V was less variable, R further darkened. This is real **skyward atmospheric-scattering-sensitive zenith photometry**, not remotely sensed top-of-atmosphere VIIRS DNB.
- **Important combined mechanisms, NOT an undimmed LED retrofit:** lamps changed HPS 250W/150W to LED 150W/70W, diffusers/ULOR changed in some cases, HPS were reportedly set to **80%** after midnight and LEDs **60%**. The reference is consequently **not** a pure test of VIIRS spectral bias independent of dimming/shielding.
- UCM photometer located **40°27′04″N, 3°43′34″W** on UCM campus; its skyglow field of view is not the VIIRS mean over Madrid municipality and certainly not the whole ES300 Madrid region.

## EXACT numerical ground source used in first figure

**Table 6** in the **main 2021 paper**, not a digitised plot or simulated observations, gives **after-midnight yearly median AstMon measurements** for three landmark years, Johnson B/V/R mag arcsec⁻² and colour indices:

| Year | B | V | R | B−V | V−R |
|---|---:|---:|---:|---:|---:|
| **2011** | **19.15 ± 0.17** | **17.82 ± 0.24** | **17.11 ± 0.25** | **1.32 ± 0.20** | **0.78 ± 0.08** |
| **2015** | **19.27 ± 0.18** | **18.07 ± 0.22** | **17.42 ± 0.23** | **1.20 ± 0.11** | **0.66 ± 0.04** |
| **2019** | **18.98 ± 0.14** | **17.95 ± 0.16** | **17.44 ± 0.18** | **1.02 ± 0.08** | **0.51 ± 0.03** |

These are **only the three published anchor years**. The plot never invents or linearly estimates 2012–14 / 2016–18 sky colour. The original ± values are transcribed as published rather than silently interpreted as confidence intervals.

To put bands in intuitive radiance-increase/decrease direction while preserving different arbitrary band zero points, calculate **\`100 × 10^[-0.4 × (m_y − m_2015)]\`**, separately within EACH band. Thus **2015 = 100 per band**; Table 6 yields 2015→2019 **Johnson B approximately +30.6%, V +11.7%, R −1.8%**. This is an *individual-band, within-site relative brightening measure*, not an absolute physical radiance comparable between bands or to VIIRS. The paper's quoted **32% B** uses its own summary and/or rounding rather than the exact rounded three-year table values. **Higher astronomical magnitude means darker sky; lower B−V means bluer sky.**

For the FULL observational time series, the authors separately released [Zenodo 10.5281/zenodo.4633001](https://zenodo.org/records/4633001): \`AstMon_evol.csv\` **45.5 MB**, \`SQM_evol.csv\` **26.8 MB**, \`TESS_evol.csv\` **67.9 MB**. The first comparison uses the already-published audited Table 6 to avoid misparsing the raw observation flag schema and substituting incorrectly filtered records. The future follow-up can use all annual or monthly values from Zenodo **only after inspecting the original raw columns, moon/cloud quality flags and local measurement times**; do not fabricate all-year ground series from Figure 22.

## Original historic VIIRS

The study uses two **separate** VIIRS products/areas:

1. **ES300 Comunidad de Madrid**, **2014–2024 annual VNL V2** fixed NUTS2016 from \`data/processed/fixed_nuts2016_viirs_2014_2024.parquet\`. Annual 2014/2024 values are checked against the locked \`nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv\`, and the **EXISTING audited Figure1 2021 native-pixel snow patch** is merged as in the previous Coventry/Échirolles figures. Source inputs are read-only. No automatic 2021 reconstruction, no double 2017 offset and no changes to PRC.
2. **Madrid MUNICIPALITY**, **2015–2024 monthly VNL V1 VCMSLCFG** \`avg_rade9h\` and \`cf_cvg\` from the existing **75N060W** EOG tile. The original 120-month monthly inventory and imported rasters are reused. Mask is **official City Council municipal boundary**, NOT the broader ES300 NUTS3: [Ayuntamiento de Madrid administrative polygon service, \`TERMINO_MPAL_AYTO_MADRID\` FeatureServer layer 4](https://services-eu1.arcgis.com/UrtTdLP2PYSxwtI6/arcgis/rest/services/L%C3%ADmites_administrativos_oficiales_del_Ayuntamiento_de_Madrid/FeatureServer/4). Requires exactly one valid GeoJSON feature, WGS84 expected location, Puerta del Sol inside and total area **550–665 km²** (Madrid city nominal ~605 km²). The fixed polygon is cached locally after validation. The source municipality boundary is contemporary; this is a constant historical spatial support, not a claim that political boundaries were reissued every year.

The first plot focuses monthly on **2015–2020** while displaying **full 2014–2024 annual NUTS3**. The \`--monthly-end 2024\` option shows 2015–24 monthly. The same European tile includes Madrid; **do not request 75N060E**.

## Figures, outputs and possible interpretation

[\`src/plot_madrid_viirs_robles_2021.py\`](../../../src/plot_madrid_viirs_robles_2021.py) generates **one four-panel figure** in PNG/PDF/Illustrator-editable SVG:

A. **2014–2024 annual locked ES300 VIIRS**. Shaded *broad*, independently documented 2015–2016 major LED/HPS retrofit phase.
B. **Madrid city monthly EOG VNL V1**. Shaded 2015–16 retrofit, no artificial exact “first/last lamp switched” date.
C. **Published 2011, 2015, 2019 after-midnight UCM AstMon B/V/R**, each normalized 2015=100 using correct logarithmic magnitude → linear-radiance conversion. The 2011 AstMon point precedes the VIIRS 2014 study.
D. **Published AstMon B−V and V−R colour index** with the reported spreads. Declining colour index means bluer sky after LEDs.

Output figure:
\`figures/local_lighting_interventions/madrid_ground_viirs/madrid_ES300_annual_municipal_monthly_UCM_AstMon_table6_to_2020.{png,pdf,svg}\`

Output independent CSVs in \`data/processed/local_lighting_interventions/madrid_ground_viirs/\`:
- \`madrid_ES300_locked_annual_vnlv2_2014_2024.csv\`
- \`madrid_municipality_vcmslcfg_2015_2024.csv\` (120 rows, valid pixel/coverage stats)
- \`robles2021_table6_after_midnight_yearly_medians.csv\` (15 published metrics, normalized optical brightness).

**Interpretation:** If VIIRS falls while Johnson B increases, this supports a *disagreement among spectral/spatially different metrics* during a mixed retrofit. It **does not prove** VIIRS spectral bias alone; optical shielding, rated lamp output, night-time dimming and other sources also changed. The ground skyglow is not an observed polarization contrast (PRC), and cannot be used to assert an improved/deteriorated animal cue without the full transfer model.

## EcoVis: GitHub-first workflow

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures

# Synthetic test: photometric conversion and actual PNG/PDF/SVG rendering
python src/plot_madrid_viirs_robles_2021.py --self-test
# Verify the already imported EOG monthly source
python src/plot_madrid_viirs_robles_2021.py --inventory
# Validate municipality before reading 120 large raster tiles
python src/plot_madrid_viirs_robles_2021.py --check-boundary
# 2014–24 NUTS3 annual + city 2015–20 monthly + published ground Table 6
python src/plot_madrid_viirs_robles_2021.py

# Optional, same frozen annual/ground, municipal monthly plot to end-2024
python src/plot_madrid_viirs_robles_2021.py --monthly-end 2024
\`\`\`

If Madrid municipality ArcGIS query fails, use the official council [administrative boundary downloads](https://datos.madrid.es/dataset/900012-0-limites-administrativos-mapas/downloads) to confirm the source, but do not substitute ES300 or make up polygon geometry. The first candidate graph has not yet run on EcoVis: no observed Madrid VIIRS trend is asserted until the user returns its actual output.
