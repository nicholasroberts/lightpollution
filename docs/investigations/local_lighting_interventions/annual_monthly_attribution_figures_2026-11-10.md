# Attributed switch-off cases: paired annual NUTS3 and municipal monthly figures

**Committed 10 November 2026, \`publication-figures\`.** The plotting script was authored and its Matplotlib date-marker transform corrected using an independent synthetic plotting API check. **The full two-city script has not yet run on the actual frozen EcoVis annual parquet / existing audited 2021 cache** in this conversation, so no PNG/PDF/SVG outcome figure is represented as already generated.

## Aim

For the two strongest independently documented **operating-hour switch-off** events, make visually comparable figures with:
- **A. Locked annual 2014–2024 NUTS2016 VNL V2 radiance** of **UKG33 Coventry** or **FRK24 Isère**, including exact previously approved 2021 native-pixel snow patch and 2017 absolute offset *once* through existing frozen inputs.
- **B. Actual EOG VNL V1 2015–2024 MUNICIPAL monthly mean** from existing case-study 120-month CSV, zoomed to **2022–2024 Coventry** or **2020–2024 Échirolles**. No direct subtraction/mixing of VNL V1 municipality monthly and VNL V2 region annual signals.
- **C. Intervention-year monthly municipal anomaly** relative to *same-calendar-month* average of independently selected reference years, Coventry 2022+23 or Échirolles 2020+21, to expose seasonality.
- Exact **documentary labels for council approval vs installation/phase-in vs actual switching** and subsequent LED conversion in a distinct colour. Do **not** choose a change-point by inspecting VIIRS.

## Event-marking specification

| Location | Annual region | Municipal intervention markers | Post-event / subsequent LEDs |
|---|---|---|---|
| **Coventry** | UKG33, full 2014–24 | **20 February 2024** approval; shade **May–July 2024** confirmed part-night programming, mark May as start | First entire post-rollout month **Aug 2024**; shade/plot Aug–Dec anomaly. Large city LED conversion did not start until **Nov 2025** — annotate as outside the 2014–24 record, DO NOT show as if within 2024 |
| **Échirolles** | FRK24 Isère, full 2014–24 | **29 July 2022** first documented specific switch-off; **1 Aug** parks and municipal property; **Oct–Dec 2022** staged roadway switch-offs; shade **Jul–Dec 2022** as broad operating phase-in, mark October start distinctly | Jan 2023 first month after roadway phase ended; shade **from June 2023** subsequent LED rollout as mixed-mechanism period. Do not call all of Oct–Dec 2022 uniformly fully treated |

**Annual points are per-calendar-year observations**; an annual 2022 or 2024 point represents a full year's composite, spanning pre- and post-intervention periods. Annual area in Isère contains other municipalities unaffected by Échirolles switches. Exactly when each lamp switched off is not inferred from satellite imagery.

## Plotting code and sources

[\`src/plot_coventry_echirolles_annual_monthly_attribution.py\`](../../../src/plot_coventry_echirolles_annual_monthly_attribution.py)

- Frozen annual source: \`data/processed/fixed_nuts2016_viirs_2014_2024.parquet\` and \`figures/prc_radiance_spider_moth_2014_2024/nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv\`, validated with existing read-only [\`export_coventry_echirolles_nuts3_bridge.py\`](../../../src/export_coventry_echirolles_nuts3_bridge.py).
- Previously verified 2021 cache and JSON audit must exist at \`data/processed/publication_figures/fig1_2021_locked_snow_radiance.csv\` (and \`.json\`). The script **fails rather than reconstructing or silently substituting** an unapproved value; if absent use the existing publication Figure 1 workflow to build the original audited cache *before* plotting.
- Existing municipal-monthly case tables: \`data/processed/local_lighting_interventions/five_case_monthly/coventry_2015_2024_monthly.csv\` and \`echirolles_2015_2024_monthly.csv\`; these come from \`src/local_viirs_five_case_analysis.py\`, with \`cf_cvg>0\` valid-data masks. They remain separate from annual VNL product.
- Published source documents for both municipality event timelines: [Coventry 30 Jan 2025 post hoc implementation scrutiny](https://edemocracy.coventry.gov.uk/mgAi.aspx?ID=43074), [Coventry principal LEDs from Nov 2025](https://www.coventry.gov.uk/street-scene/coventry-led-streetlighting-upgrade), [Échirolles July/Aug operations](https://www.echirolles.fr/actualites/un-nouveau-plan-de-sobriete-energetique-pour-la-ville), [Échirolles Sep/Oct city periodical](https://www.echirolles.fr/sites/default/files/2022-10/CITE%20ECHIROLLES%20SEPT%20OCT%2022.pdf), [Échirolles summer 2023 LED](https://www.echirolles.fr/sites/default/files/2023-10/Cite%20Echirolles%20N404-web2_compresse.pdf).

## How to run and review

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures

# Test source identities, data structures and actual synthetic figure rendering
python src/plot_coventry_echirolles_annual_monthly_attribution.py --self-test

# Two annual + monthly + monthly-anomaly panels, one figure per place
python src/plot_coventry_echirolles_annual_monthly_attribution.py
\`\`\`

Each figure is saved as **PNG, editable-text SVG and PDF** in \`figures/local_lighting_interventions/nuts3_event_attribution/\`:

- \`coventry_nuts3_annual_and_municipal_monthly_2014_2024.{png,pdf,svg}\`
- \`echirolles_nuts3_annual_and_municipal_monthly_2014_2024.{png,pdf,svg}\`

Tables saved in \`data/processed/local_lighting_interventions/nuts3_event_attribution/\`:

- \`coventry_isere_plot_validated_annual.csv\`: 22 verified annual rows
- \`coventry_event_year_monthly_anomalies.csv\`: 2024 12-month municipal comparisons
- \`echirolles_event_year_monthly_anomalies.csv\`: 2022 12-month municipal comparisons

**Interpretation guardrail:** The source annual radiance may be affected by changes in spectral shape, shielding, luminous flux, lamp hours, other private lighting and atmosphere. We are *not* directly estimating the fraction of regional 2014–24 ΔR attributable to a municipality's policy from these plots.

## The Kyba citation relevant to this design

**Kyba, C. C. M. et al. (2021).** [Direct measurement of the contribution of street lighting to satellite observations of nighttime light emissions from urban areas](https://doi.org/10.1177/1477153520958463). Using controlled Tucson streetlight output experiments with satellite VIIRS, estimated **13% of measured urban DNB radiance** came from city streetlights under normal dimmed post-midnight settings (95% confidence interval 10–16%), increasing to **18% if undimmed** (95% CI 15–23%). Including other operated streetlights raised estimated fractions to 16% dimmed and 21% undimmed. Tucson is not a universal European streetlight fraction; the study **demonstrates that other emissions (building windows, businesses, industry, advertising, vehicles and facilities) can dominate a citywide satellite radiance record**. This supports not treating a municipality-average as the signal of its council-owned street lighting alone.

See [two prospective LED-dominated contrasts and rejected pure LED controls](./led_dominant_contrast_candidates_2026-11-10.md).
