# Figure 2A–D: frozen Jan-1 population and radiance, vector production

**Status: publication plot script committed; actual EcoVis SVG/PDF generation and visual QA pending.**

## Scientific source: no new analysis

Use `src/plot_publication_figure2_population_radiance.py`, based **only** on the previously frozen `1,273` NUTS2016 NUTS3 January-1 population / fixed-land-area cohort for **2014–2020**. This is NOT the older Eurostat `demo_r_d3dens` population-density product and NOT the 1,522-region radiance/PRC cohort used for Figure 1.

Primary frozen input:

```text
data/processed/jan1_population_radiance_prc_country_nl333/
  matched_region_change_2014_2020.csv
```

Created by the existing, previously run `src/analyze_jan1_country_radiance_prc_nl333.py`, which performs the read-only 1,273-region radiance join using `src/plot_jan1_population_prc_change_symlog.py`. It contains `density_jan1_fixed_land_2014`, `density_jan1_fixed_land_2020`, `radiance_2014`, `radiance_2020`, `change_density_jan1_fixed_land_2014_2020`, and `radiance_delta_2014_2020`.

Optional frozen statistics file, if present (not recomputed):

```text
data/processed/jan1_population_radiance_prc_country_nl333/
  within_country_correlations_2014_2020.csv
```

The figure script verifies 1,273 unique regions; that both endpoint differences equal their frozen reported deltas; positive population densities; nonnegative radiance; and original Spearman sample sizes where the frozen table is found. It does **not** regenerate, refit or alter any source datasets. The only plotted regression lines are descriptive OLS fits in log space following the project's original `fit_population_radiance.py` / `plot_country_facets.py` graphics conventions.

## Individual figures (each SVG and vector PDF)

- **2A: Broad spatial log–log relationship (2014).** All matched positive-radiance cases, with each fixed NUTS2016 NUTS3 a point. 2014 descriptive OLS display line and labelled slope.
- **2B: Spatial relationship comparison (2014 versus 2020).** Both annual clouds and their separate descriptive log–log lines, **using exactly the same regions positive in both years**; same axis limits as 2A to avoid artifactual fit rotation due to changing plotting cohorts.
- **2C: Temporal changes across all 1,273 NUTS3.** Absolute change in Jan-1 population density (people/km²) versus change in calibrated regional mean VIIRS DNB radiance (nW cm⁻² sr⁻¹), 2014→2020. Both axes use signed-log display to retain zero, small changes and outliers; no cases removed. If the frozen correlation CSV exists, annotate the pre-existing all-region Spearman statistic; otherwise omit the number rather than invent or recompute it.
- **2D: Country facets, regional temporal changes.** The **same six previously selected countries** `TR, DE, FR, IE, IT, NL` (from the frozen `analyze_jan1_country_radiance_prc_nl333.py` `DISPLAY` constant). Each facet shows regional signed density changes vs VIIRS radiance changes with the **same physical axes and shared limits** as 2C; add original country Spearman values from frozen CSV if available; identify NL333 without excluding it. No new country classifications.

**Important:** Panel 2D presents *within-country NUTS3 regional observations*, not *one point per nation*. A between-country country-aggregate plot would be a different object and is not substituted silently.

### Visual standards

- Editable SVG + vector PDF, shared `src/figure_style.py`, default DejaVu Sans.
- 2A/2B use log10 physical density/radiance as per the existing project graphics; positive radiance at **both** endpoints is needed for paired log comparison. Report any omitted nonpositive-radiance cases in provenance.
- 2C/2D use **linear units on signed-log axes** (`linthresh` defaults: 5 people/km² and 0.1 nW cm⁻² sr⁻¹), not percent population changes or transformations of the raw datasets.
- All points remain vectors; no inference from presentation colours.
- No 2021 snow adjustment is involved because the endpoint is 2020.
- No demographic, radiance, country-correlation, or biological calibration is reopened.

## Run on EcoVis

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/plot_publication_figure2_population_radiance.py --self-test
python src/plot_publication_figure2_population_radiance.py
```

Outputs under `figures/publication/figure2_population_radiance/`:

```text
fig2A_spatial_population_radiance_2014.svg
fig2A_spatial_population_radiance_2014.pdf
fig2B_spatial_relation_2014_2020.svg
fig2B_spatial_relation_2014_2020.pdf
fig2C_population_vs_radiance_change_2014_2020.svg
fig2C_population_vs_radiance_change_2014_2020.pdf
fig2D_country_population_radiance_change_2014_2020.svg
fig2D_country_population_radiance_change_2014_2020.pdf
fig2_descriptive_plot_line_parameters.csv
fig2_provenance.json
```

If working from the project root but the virtual environment is already activated, omit the `source` command.

## Quality checks after rendering

Inspect editable SVG text and file sizes in Illustrator. Compare displayed country ρ to the frozen `within_country_correlations_2014_2020.csv` (expected TR +0.621, DE +0.410, FR −0.410 for Δdensity vs Δradiance). Verify all **1,273** regions are in the 2C point sample and that the 2A/2B paired log-eligible count is explicitly recorded (must not silently log zeros). NL333 should appear as a star in the Netherlands panel. If any source table is absent, locate the existing previously produced EcoVis CSV rather than rerunning the Eurostat reconstruction.

Output files are ignored by Git; figure code and notes are committed. The user will determine the final talk/paper narrative and panel arrangement.
