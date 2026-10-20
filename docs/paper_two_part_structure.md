# Manuscript structure — a decade of European radiance and biological polarization cues

**Agreed concept (20 October 2026): TWO main empirical results, both biologically motivated.** 2014–2024 represents **eleven annual observations covering a ten-year interval**. The manuscript must include **both spiders and moths** from the original research, with separate species/taxon transfer relationships, and avoid getting lost in minor calibration refinements.

## Working thesis

**Observed night-light radiance trends are not interchangeable with their predicted biological consequences.** The first result creates the most defensible possible radiance history (including known instrument/weather/spectral biases); the second propagates changes through empirically measured **polarization-receptor contrast (PRC)** relationships for nocturnal arthropods, highlighting why equal radiance changes can lead to unequal changes in the celestial polarization cue.

Possible working title: *A decade of changing artificial light across Europe and its unequal consequences for celestial polarization cues*.

## Main result 1 — What happened to Europe's night-light radiance, 2014–2024?

1. Area-weighted **native-pixel VIIRS annual radiance** on a fixed NUTS2016 NUTS3 basis, retaining complete 2014–2024 series, map absolute and percentage radiance changes and verify country trajectories and sharp discontinuities.
2. **Baseline zero/bias shift**: apply existing validated **2017 onward +0.235 nW background subtraction** exactly once to eligible VIIRS pixels; report baseline sensitivity and stability.
3. **2021 snow/surface albedo**: winter 2021 spatial native-pixel gain candidate fitted on 18 heldout NUTS pilot regions, two winter months; do not brighten annual pixels; spatial candidate removes ~48.5% of European 2021 PRC anomaly but **leaves substantial residual**. Status is candidate **not locked**; resolve official annual outlier-treatment/source effects as efficiently as possible and preserve raw comparison.
4. **LED spectral shift**: integrate ISS B/G and G/R changes with VIIRS DNB, spectral-response physics, and independent local infrastructure/lighting policies. Treat correction as **uncertainty-bounded scenario**, not a claimed exact emitted-light reconstruction without actual source spectra/observation-year timeseries. Compare corrected trend and original raw DNB trends for 2014–24.
5. **Country contrasts** as mechanistic cross-checks: France (documented extensive municipal switch-off and Enedis 2022/23 02:00 power reduction), UK (dimming + LED), Germany/Austria (less ISS colour shift in earlier epoch), Switzerland (mixed official VIIRS trend), Italy/Milan (early intense LED conversion). Documented policy changes are *hypotheses to investigate* at regional scale, not per-NUTS scalar corrections.

Proposed Fig. 1: absolute and % DNB radiance maps, Europe annual history, corrected/raw overlays and country contrasts; supplementary 2021 spatial QA/sensor data. Proposed Fig. 2: ISS spectral ratios plus the LED-sensitivity method or scenarios and their impact on interpretation of radiance decrease.

## Main result 2 — Radiance change ≠ predicted biological consequence

**Empirical transfer:** Both **spider PRC(R)** and **moth PRC(R)** from the original R code. Fit/validate each separately using independently verified site/date/VIIRS source radiances. The spider **17-site transfer is already locked** after excluding anomalous-zero annual EOG at Slabhouse; do not reopen it for minor gains. **The moth transfer is NOT yet validated, fitted or locked** and must not reuse spider coefficients. The moth dataset had **17 PRC values** in original code, excludes the first listed CSV and additionally the fifth; its original source order has been captured in \`data/reference/prc_moth_17_observations_from_original_r.csv\`. Confirm exact R file order, label site mapping/years, investigate Slabhouse annual-EOG zero consistently, and retain original moth observations even if one is excluded from a fitted predictor calibration.

Pixel-scale biology: for taxon \`t\`, compute \`PRC_t(R_pixel)\` on each valid satellite pixel (and taxon-relevant light quantity, if spectral model allows) **BEFORE** spherical-area aggregation to fixed NUTS3. Display observed/predicted cue availability **not actual navigation success or a proven behavioural visibility threshold**.

Proposed Fig. 3: **paired maps absolute ΔR vs absolute ΔPRC**, and **paired %ΔR vs %ΔPRC**, for spiders and moths or carefully comparable layouts; NUTS3 by-country scatter with ΔR horizontal and ΔPRC vertical and original regional absolute/percent values. Initial locked-spider 2014→2024 visuals already produced:
- \`figures/prc_radiance_comparison_locked17/compare_cumulative_absolute_radiance_vs_prc_2014_2024.pdf\`
- \`figures/prc_radiance_comparison_locked17/compare_cumulative_percent_radiance_vs_prc_2014_2024.pdf\`
- \`figures/prc_radiance_comparison_locked17/scatter_prc_vs_radiance_change_by_country_2014_2024.pdf\`

**Case studies:** select (a) substantial radiance decrease with large PRC benefit, (b) substantial radiance decrease with **little** PRC benefit (saturation or radiance distribution), (c) radiance increase in previously dark habitats with outsized predicted impact, (d) areas with major colour shift where blue light casts doubt on DNB-only benefit. Specifically investigate **NL333** as user-nominated case of potentially large ΔR with muted ΔPRC, but **NL333 was flagged as an anomalous extreme in earlier 2021/annual radiance work**: check original EOG raw pixels, masking, region validity and source provenance before presenting it as a genuine lighting-change exemplar.

**Population density as an analysis nested within Result 2, not a third major manuscript claim:** Stratify NUTS3 (ideally also pixel/consistent-gridded) baseline population densities, population density changes and initial radiance to discover where **absolute ΔPRC** and/or **relative ΔPRC** are largest for spiders vs moths. Work with exposure-weighted vs geographic area-weighted results separately; population density is explanatory *human activity* proxy, not animal abundance, habitat availability, or biological effect itself. Do not assume maximum biological loss in highest-density cities: saturating PRC can make marginal sensitivity greatest in dim/medium landscapes, while largest changes in R often occur in bright places. Consider an interpretable 2D initial radiance × density binned surface and region-within-country effects. This section has the potential to be the ecological centre of the manuscript.

Proposed Fig. 4: spider/moth PRC versus R curves and derivatives, two-taxon ΔPRC maps or density-stratified biological effect surfaces, density × initial radiance 2D grid, and selected grounded NUTS3 case-study trajectories.

## Methodological and causal guardrails

- Optical estimands are not interchangeable: DNB 500–900 nm radiance, broad optical flux, photopic or scotopic emission, skyglow, and taxon-specific polarization contrast need explicit wavelength sensitivity.
- ISS spectral RGB colour ratios span irregular **2012–13** and **2014–20** observation mosaics, not per-year time-series. Do not imply annual LED conversion years identifiable from the ISS ratios alone.
- Do not automatically feed an estimated broadband LED-corrected R into the empirical DNB-based PRC fit. Test the units and assumptions; prefer DNB-calibrated primary inference with independently bounded spectral scenarios until a receptor-specific spectral mapping is validated.
- 2014 and 2024 endpoint changes **do not depend on 2021-only correction**; the 2021 work impacts within-series consistency, event detection and interpretations involving the anomalous year.
- Don't divide regional mean PRC by population: mean PRC is not a count of individual animals, and region density estimates do not establish where arthropods are located.
- **No causal attribution from national lighting policies purely from co-occurring satellite changes**. Use independent implementation dates/power evidence where possible and identify ecological mechanisms as predictions rather than observed behaviour.
- Document observational uncertainty, 2021 residual and LED correction assumptions; do not suppress original annual observations.
- Do not dilute the paper with another exploratory search for marginally different spider PRC coefficients.

## Concrete next steps

1. Finish one bounded investigation of the leftover 2021 discontinuity and choose a documented annual correction / explicit sensitivity, with QA frozen.
2. Retrieve matching **original moth DRAc vector and filenames** and validate against independently geolocated 2018/2019 sites, then fit a distinct moth model (same biological calibration philosophy, not necessarily same model family if evidence demands otherwise).
3. Overlay published **ISS B/G and G/R ratio products** and quantify country/NUTS3 stratified spectral changes, then compare plausible LED SPD-to-VIIRS spectral-response correction bounds with documented policy changes. Constrain correction, do not assign exact visible light from two colour ratios.
4. Recompute regional population density categories on the SAME fixed NUTS2016 cohort; analyse absolute and percent ΔPRC with two taxa by density and baseline radiance.
5. Choose real NUTS3 biological case studies with explicit pixel-level source QA, first auditing NL333.
6. Write the manuscript's Introduction/Discussion around the dual result and position relative to Li et al. 2026 *Nature*, Sánchez de Miguel et al. 2022 *Science Advances*, Kyba et al. 2017 and Cerema 2025.

**Literature support/provenance:** [European dimming, policy and LED source review](literature_europe_nighttime_dimming_policy_led_2014_2024.md).
