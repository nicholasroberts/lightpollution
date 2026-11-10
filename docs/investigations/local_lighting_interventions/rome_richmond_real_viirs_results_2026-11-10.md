# Rome and Richmond: *actual* LED-dominated VIIRS results and scientific interpretation

**Reviewed 10 November 2026, branch \`publication-figures\`.** Based on investigator-uploaded *actual-run* figures and the following CSVs, not on simulated results:
- \`rome_led_contrast_annual_and_municipal_monthly_2014_2024.png\` and \`rome_led_contrast_nuts3_annual_only_2014_2024.png\`
- \`richmond_led_contrast_annual_and_municipal_monthly_2014_2024.png\` and \`richmond_led_contrast_nuts3_annual_only_2014_2024.png\`
- \`rome_municipality_monthly_vnlv1.csv\`, \`rome_monthly_same_month_anomaly.csv\`
- \`richmond_municipality_monthly_vnlv1.csv\`, \`richmond_monthly_same_month_anomaly.csv\`
- \`rome_richmond_locked_annual_NUTS3_2014_2024.csv\`

These are **interpretation / data audit notes**, not manuscript text. The original uploaded PNG/CSV results are on this conversation, and the EcoVis production outputs are in \`figures/local_lighting_interventions/led_dominated_contrasts/\` and \`data/processed/local_lighting_interventions/led_dominated_contrasts/\`. Values computed by directly grouping user-supplied CSV records by calendar year, with the annual NUTS2016 values taken unmodified from the user's verified 2021-snow-corrected export. Radiance units throughout: nW cm⁻² sr⁻¹.

## Summary table — actual observed outcomes

| Case/geography and independent year comparisons | Baseline → outcome mean | Mean ΔR and % change | Months declining / 12 |
|---|---|---|---|
| **Rome municipality VNL V1 2016→17** | 22.192 → 19.625 | **−2.567 (−11.57%)** | **11** |
| **Rome municipality VNL V1 2017→18** | 19.625 → 17.840 | **−1.785 (−9.10%)** | **11** |
| **Rome municipality VNL V1 2018→19** | 17.840 → 17.779 | **−0.061 (−0.34%)** | **5** |
| **Rome municipality VNL V1 2016→19** | 22.192 → 17.779 | **−4.413 (−19.89%)** | **12** |
| **Richmond borough VNL V1 2018→19** | 17.712 → 16.645 | **−1.067 (−6.02%)** | **9** |
| **Richmond borough VNL V1 2019→20** | 16.645 → 15.446 | **−1.199 (−7.20%)** | **9** |
| **Richmond borough VNL V1 2021→22** | 15.118 → 14.427 | **−0.691 (−4.57%)** | **5** |
| **Richmond borough VNL V1 2022 vs 2017/18 average** | 17.964 → 14.427 | **−3.537 (−19.69%)** | **11** (Jan's anomaly is positive) |
| **Richmond borough VNL V1 2023→24** | 14.685 → 12.637 | **−2.048 (−13.95%)** | **9** |

*Method:* Yearly mean over 12 monthly municipality radiance means, with each month equally weighted. \`2017/18 average\` is the arithmetic average of the matching two prior same-calendar months. Counts for that comparison use \`richmond_monthly_same_month_anomaly.csv\`. All are descriptive differences, **not** estimates of per-lamp or causal treatment effect.

### Frozen, approved annual NUTS2016 radiance

| Region | 2014 | 2016 | 2017 | 2018 | 2019 | 2023 | 2024 | Relevant changes |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| **ITI43 Roma metropolitan** | 9.041 | 9.043 | 8.433 | 7.866 | 7.882 | 8.177 | 8.490 | **2016→18 −13.01%; 2018→24 +7.93%; 2014→24 −6.10%** |
| **UKI75 Hounslow and Richmond** | 23.676 | 22.601 | 21.652 | 21.483 | 20.528 | 19.979 | 17.659 | **2023→24 −11.61%; 2014→24 −25.42%** |

The wide parent NUTS3s do not have the same footprint as municipality monthly masks. Do not compare **absolute** radiance values or imply that an annual NUTS3 delta equals the municipal LED effect.

## Rome — *the strongest counterpoint to a universal 'radiance = emitted lumens' interpretation*

Official audited installed *streetlight* luminous flux from **Rome ACoS Relazione Annuale 2019, table 5.2**:

| Year | Total installed luminous flux, kilolumens | Estimated streetlighting energy, MWh | Municipality VNL V1 annual mean |
|---|---:|---:|---:|
| **2016** | **2,750,046** | **167,852** | **22.192** |
| **2017** | **1,991,000** | **115,640** | **19.625** |
| **2018** | **2,010,000** | **84,000** | **17.840** |
| **2019** | No comparable 2019 flux value in the table 5.2 source | Operator reports further energy savings | **17.779** |

**Mechanism separation at the same city scale:**
- **2016→17:** installed streetlight luminous flux **−27.60%**, while municipal VIIRS **−11.57%**. During the phased LED conversion, an actual fall in rated installed streetlight lumens accompanied the VIIRS decline; spectral replacement and private lighting may also have contributed. No one-to-one radiance/flux response can be inferred.
- **2017→18:** installed streetlight luminous flux **+0.95%** (approximately constant), while municipal VIIRS **−9.10%** in **11/12 matching calendar months**. Council energy use **−27.36%** as LED conversion progressed. This is independently grounded **divergence between installed municipal streetlight lumens and VIIRS radiance**. It is **consistent with** altered LED spectral response, directional shielding, changed night-time operating profiles, and/or independently operated emitting sources. It does **not** *prove* LED spectral bias or unchanged **actual** 01:30 luminosity: official flux is *installed/rated*, NOT measured light emitted at VIIRS acquisition time.
- **2018→19:** municipal monthly VIIRS **−0.34%** (near flat, **5/12** months declining), with annual ITI43 VNL V2 almost flat **+0.20%**; strengthens the description that the major radiance transition was concentrated in **2016–18** while later levels stabilised.
- **Caution on timing:** Jan 2018 was already **−3.87** relative to Jan 2017 and late 2018 differences were generally smaller, not a within-2018 isolated off-switch step. 2017 itself already contained very widespread LED conversion. Avoid labelling **2017–2018** an experiment with identical streetlight inventories or controlled lamp schedules.
- **Quality:** All **12 months in both 2017 and 2018** have **municipality median \`cf_cvg\` ≥8**. Excluding October 2018 by requiring >=10 for both years leaves **11 same-month pairs**, with **−9.35%** instead of **−9.10%**. This rules out a *simple low-count whole-month* explanation for the reported 2018 difference but does **not** guarantee matched-pixel/spatial support, snow/atmospheric invariance, or constant nonmunicipal emissions.
- **Installed flux over entire LED programme was NOT flat:** 2014 **3,377,000 klm** vs 2018 **2,010,000 klm** (~**−40.5%**); the 2016→17 true flux reduction is substantial. Rome is *LED-dominated*, not LED-only/no real dimming for the entire programme.

**Original external evidence:**
- [ACoS Rome, 2019, Section 5.1, Table 5.2](https://www.agenzia.roma.it/documenti/schede/relazione_annuale_2019_cap_5_servizi_a_rete.pdf): 2016–2018 total installed flux and estimated energy, independently monitored programme.
- [ACEA sustainability report, public streetlights 2019](https://reports2019.gruppo.acea.it/sostenibilita/sostenibilita/en/relations-stakeholders/customers-and-community/quality-delivered): by end 2019, 205,670/225,730 lamps LED (91%); >180,000 LEDs installed across programme; not per-night lumen output.
- [Kyba et al. 2021, experimental Tucson source partition](https://doi.org/10.1177/1477153520958463): municipality-operated streetlighting was about **13%** of Tucson VIIRS DNB after midnight, not universally the same share in Rome/Richmond; independently operated light sources can dominate.

## Richmond — a mixed, diffuse signal with a late-2024 anomaly

- **The regional trend begins long before the funded December 2019 LED scheme:** UKI75 falls from **23.676 in 2014** to **21.483 in 2018** (pre-funding, **−9.3%**). Richmond borough municipal monthly also declined **2017→18 (−2.8%)**, **2018→19 (−6.0%)**, **2019→20 (−7.2%)**. A 2019 finance announcement is *not* an exact change point, and the past month-by-month borough LED replacement fraction is still not documented.
- The **2022 against 2017/18 same-calendar-month** average is **−19.69%**, yet this is a multi-year comparison **without a clean untreated 2017/18 baseline** or verified constant 01:30 output. Funding [describes same-design light levels and same operating hours](https://www.amberinfrastructure.com/sectors/case-studies/london-borough-of-richmond-meef/), while potentially programmable dimming existed. It is not a clean 'LED-only undimmed' control.
- The abrupt **2023→24 borough −13.95%** is mostly seasonal: **January–October alone −4.2%**; **November 2024 14.91→7.46** and **December 20.85→9.65** relative to matching 2023 months. The regional UKI75 annual decline is **−11.61%** over 2023→24, but includes Hounslow and unrelated light sources. **Do not attribute this late-2024 step to Richmond's 2019–22 LEDs**.
- Observation support: 2024 Richmond median municipal \`cf_cvg\` = **2** in July and **4** in April; Nov 2024 **7**, Dec 2024 **11**, whereas Nov 2023 **9** and Dec 2023 **6**. Removing only those months with either year's median \`cf_cvg\` <5 yields **−16.35%** across **10 remaining matching months**, still containing November/December; raising median \`cf_cvg\` ≥8 leaves only **5 months, −1.5%** (non-comparable seasonal subset). This supports further month-level and common-pixel checks before interpreting a year-on-year change.
- Pre/post winter seasonal differences can reflect business, sports, traffic, airport lighting, atmospheric/cloud/snow and coverage differences. The 2024 regional step **cannot currently be attributed to LED conversion** on these data.

## Four-case synthesis — operational vs LED-dominated

| Documented intervention | Municipal VNL signal | Independent light-output / operation evidence | What this does and does not demonstrate |
|---|---|---|---|
| **Coventry**, part-night off May–July 2024 | Aug–Dec 2024 **−32.23%** vs 2022/23 same months | Council retrospectively verifies programming completed by July; main LED programme **2025+** | Coherent evidence of actual operating-hours reduction associated with marked VIIRS fall, not full causal proportion of UKG33 decade-long decline |
| **Échirolles**, staged switch-offs July–Dec 2022 | Oct–Dec 2022 **−37.08%** vs 2020/21 | Operating switches verified; low municipal LED stock ~5% in 2022; stronger LED replacement from 2023 | A second real-off-hours mechanism at a small spatial subset of FRK24, not causal attribution of all Isère's 2022→23 annual fall |
| **Rome**, massive 2016–2019 HPS→LED | 2017→18 **−9.10%** | Installed streetlight flux **+0.95%** in 2017→18, energy −27.4%, but no actual per-lamp 01:30 dimming history | **Divergence** between rated installed lumens and orbital VIIRS; consistent with spectral/directional change, NOT proof of pure spectral effect or no operating changes |
| **Richmond**, staged 2019–22 LEDs | 2022 vs 2017/18 **−19.69%** | Financing promised unchanged hours/comparable design brightness, but programmable dimming present and commissioning not verified | Supportive mixed example only; cannot claim LED-only 'no dimming' controls; 2024 decline heavily affected by autumn/winter observations |

**Central conclusion supported, without overclaiming:** The documented Coventry/Échirolles municipal operations and temporally aligned DNB declines contradict a simplistic *all historical European declines must be spectral artefacts* claim. Rome at the same time demonstrates that **a radiance decline can occur despite nearly unchanged installed streetlight luminous flux in a specific comparison year**, so **it is equally wrong to treat every VIIRS decline as a true reduction in lumen output or duration**. The two effects can coexist, and neither permits arbitrary quantitative apportionment across entire NUTS3 regions without independent spatial/time information.

**Priority next validations:**
1. Rome **2017–2018 ISS blue/green spectral ratio** at a tightly localised mosaic/time series, if dated, spectrally comparable ISS imagery exists; distinguish spectral shift from change in amount. No ISS image/year is presumed available.
2. Rome's **operational profile at ~01:30** and matching LED conversion phases / installed luminous flux records, ideally luminaire and road-level.
3. For Rome, **fixed paired-pixel VIIRS change maps** over the actual municipality, and predeclared strata of lamp-bearing vs industrial/commercial areas, rather than a false citywide causal total.
4. Richmond **2023/24 November–December diagnosis** with snow/coverage/spatial support and contemporaneous local changes before using the 2024 annual regional drop.
5. If a second *confirmed* no-dimming LED control is still required, independently source one; Richmond does NOT fulfil that experimental requirement yet.

**Frozen European 2014–2024 NUTS3 VIIRS and polarization receptor contrast unchanged.** This interpretation does not require new model fitting or recalibration.
