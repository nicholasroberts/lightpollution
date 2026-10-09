# Socioeconomic extension and Crawley-style model simplification

This branch tests whether annual socioeconomic state can explain the country-
and year-dependent structure in the NUTS3 population-radiance relationship,
with the eventual goal of replacing unexplained calendar-time effects with
variables that can themselves be projected into the future.

## Data

The first socioeconomic source is Eurostat `nama_10r_3gdp`, **Gross domestic
product (GDP) at current market prices by NUTS 3 region**. Eurostat reports
annual data through 2024.

The pipeline downloads and retains:

- total GDP, million euro (`MIO_EUR`);
- total GDP, million PPS (`MIO_PPS_EU27_2020`);
- GDP per inhabitant, euro (`EUR_HAB`);
- GDP per inhabitant, PPS (`PPS_EU27_2020_HAB`);
- GDP per inhabitant relative to the EU27 average
  (`PPS_HAB_EU27_2020`).

Eurostat `reg_area3` total land area is also used to derive:

- GDP density, euro per km2;
- GDP density, PPS per km2.

Log10 forms and descriptive year-on-year changes are written to the processed
table.

The **primary GDP covariate for the first model is GDP per capita in PPS**.
This is preferable to raw euro GDP per capita for cross-country comparison.

GDP density is deliberately not entered in the same model as population
density plus GDP per capita because:

```text
GDP density ~= population density x GDP per capita
```

and including all three would introduce structural collinearity.

## Maximum model

Following Crawley's model-simplification logic, the analysis starts with the
largest fixed-effects structure that remains biologically interpretable rather
than including an automatic four-way interaction.

The maximal model contains all main effects, all six two-way interactions and
all four possible three-way interactions among:

- log10 population density;
- log10 GDP per capita in PPS;
- year;
- country.

The four-way interaction is deliberately omitted.

Conceptually:

```text
log10(corrected radiance)
~ all main effects
+ all two-way interactions
+ all three-way interactions
+ (1 | NUTS_ID)
```

The NUTS3 random intercept is the repeated-measures structure and is held fixed
throughout simplification.

For the model-selection analysis, the working interval is **2013-2023**.
Socioeconomic data for 2024 remain in the downloaded/processed table, but 2024
is excluded from this fit because regional GDP coverage is incomplete.

## Simplification procedure

Fixed-effect comparisons are fitted with maximum likelihood (ML). The
simplification proceeds from highest-order interactions downwards while
preserving marginality:

1. fit the maximal three-way model;
2. test deletion of the three-way currently removable term(s) by
   likelihood-ratio tests against the current model;
3. delete the least significant term when its removal does not cause a
   significant loss of fit (`p >= alpha`);
4. refit and repeat at the same interaction order;
5. once all removable terms at that order are significant, move down one
   interaction order;
6. never delete a lower-order term required by a retained higher-order
   interaction;
7. stop when no further non-significant removable terms remain.

Every deletion test is written to
`data/processed/socioeconomic_model_simplification.csv`, including the model
formula, tested term, LR statistic, degrees of freedom, p value and decision.

This gives an auditable route from the maximal model to the minimum adequate
model rather than selecting among an arbitrary set of candidate formulas.

## Pipeline stages on this branch

```bash
python src/run_pipeline.py socioeconomic
python src/run_pipeline.py socioeconomic-merge
python src/run_pipeline.py socioeconomic-model
```

The GDP stage can be forced to refresh Eurostat source files directly with:

```bash
python src/fetch_eurostat_socioeconomic.py --force-download
```

The existing `main` production outputs are not overwritten by these
experimental stages.


## ISS spectral-shift validation

The persistent country-dependent changes in the radiance-population relationship
may partly reflect lighting technology rather than different human demand for
light. Sánchez de Miguel et al. (2022) provide calibrated ISS RGB mosaics for
Europe in 2012-2013 and 2014-2020.

Run:

```bash
python src/run_pipeline.py iss-spectral-validation
```

The updated validation now follows the published image-processing logic more
closely. It restricts the analysis to the paper's accepted countries, applies
the published colour-ratio quality limits, and masks pixels unless VIIRS-DNB
radiance exceeds 0.5 nW cm-2 sr-1. Because the public ISS mosaics combine many
acquisition dates, the reproducible implementation uses the annual 2013 VIIRS
raster as the pre mask and annual 2020 VIIRS raster as the post mask. This is
an explicit approximation to the paper's per-image temporal matching.

The important extension is that spectral change is now extracted at **NUTS3**,
not only as country medians.

The mechanism analysis now uses **paired common spatial support**. A pixel is
retained only if it is valid in both ISS mosaics, exceeds the VIIRS 0.5
nW cm-2 sr-1 threshold in both 2013 and 2020, and passes the colour-ratio
quality filters in both periods. NUTS3 spectral change is then summarized as
the median of the pixelwise post-minus-pre ratio changes. This prevents a
change in the set of sampled pixels from masquerading as spectral change.

For validation, the pipeline still reports separate-support Europe-wide
pre/post medians because those are the closest reproducible comparison with
the published period summaries. The paired-common-support medians and deltas
are reported alongside them but are the quantities used in the slope/intercept
mechanism tests.

For each NUTS3 region the analysis records pre/post B/G and G/R and their
changes. It then tests:

```text
delta spectral ratio ~ country + within-country log10(population density)
```

and a second model in which the population-density gradient is allowed to vary
by country. The common within-country gradient asks whether LED/spectral
conversion was systematically stronger in dense or sparse regions; the
country interaction asks whether that rollout gradient itself differed among
countries.

This distinction maps directly onto the log-log radiance relationship:

- a spatially uniform spectral change can move the apparent VIIRS level
  (intercept);
- a spectral change that covaries with population density can change the
  apparent VIIRS population slope.

For each country the pipeline therefore also estimates the 2013 and 2020
radiance-population log-log slopes on the same matched NUTS3 set and calculates
their change. The country-specific spectral-change-vs-population gradient is
then compared directly with the observed change in VIIRS slope.

For a level/intercept diagnostic, country spectral change is compared with
change in predicted log10 radiance at a common reference density of
100 people km-2 (log10 density = 2). This is more interpretable than comparing
raw intercepts at a population density of 1 person km-2, although raw intercept
changes are retained in the output.

Finally, a NUTS3 change model tests the mechanism directly:

```text
delta log10(VIIRS radiance)
~ country
+ delta log10(population density)
+ [delta log10(GDP per capita)]
+ within-country baseline log10(population density)
+ within-country spectral change
+ spectral change x population density
```

A nested extension adds country-specific population-density slopes. If the
spectral interaction is supported and that residual country-density block is
no longer required, this is the specific result we are seeking: spectral
conversion has explained the country-dependent movement of the log-log slope.

Key outputs are:

```text
data/processed/iss_country_spectral_shift.csv
data/processed/iss_nuts3_spectral_shift.csv
data/processed/iss_spectral_mask_validation.csv
data/processed/iss_nuts3_spectral_viirs_endpoints.csv
data/processed/iss_spectral_density_models.csv
data/processed/iss_country_spectral_density_gradients.csv
data/processed/iss_country_loglog_change_2013_2020.csv
data/processed/iss_spectral_vs_loglog_change_associations.csv
data/processed/iss_nuts3_viirs_change_models.csv
data/processed/iss_spectral_validation_summary.txt
```

The ISS data are still a historical mechanism test rather than a future
technology trajectory. If the mechanism is supported, the next forecasting
step is to represent spectral/LED conversion with a bounded technology-state
variable rather than linearly extrapolating calendar-year coefficients.


## Bounded LED spectral endpoint for future VIIRS-equivalent radiance

The historical country intercept trajectories are not treated as a common
European trend. The LED endpoint stage isolates only the spectral component of
the VIIRS radiance level and leaves country-specific non-spectral intercept
structure intact.

Run:

```bash
python src/run_pipeline.py led-spectral-plateau
```

The first implementation deliberately uses published empirical relations rather
than assuming that LED colour is a blackbody spectrum.

Sánchez de Miguel et al. (2019; DOI 10.1016/j.rse.2019.01.035) reported the
following DSLR synthetic-photometry relation between CCT and ISS-camera G/R
over the valid range 0.2 <= G/R <= 1:

```text
CCT / 1e4 =
    -3.0 (G/R)^4
    +5.8 (G/R)^3
    -3.2 (G/R)^2
    +1.0 (G/R)
    +0.06
```

The stage numerically inverts this relation for all-LED endpoint scenarios at
2700, 3000, 4000 and 5000 K.

Sánchez de Miguel et al. (2022; DOI 10.1126/sciadv.abl6891, Eq. 2) give the
spectral relationship:

```text
G_ISS / VIIRS = 0.21 * 1.5^(G_ISS / R_ISS)
```

so at equal green-band output the corresponding VIIRS response is proportional
to the reciprocal of this expression.

Two outputs are produced:

```text
data/processed/led_spectral_plateau_scenarios.csv
data/processed/led_country_spectral_plateaus.csv
data/processed/led_spectral_plateau_summary.txt
```

The scenario table reports the all-LED G/R endpoint and its spectral-only VIIRS
response relative to a warm 2200 K legacy reference. The country table instead
uses each country's own post-period ISS G/R as its starting spectral state and
calculates the remaining bounded correction required to reach each LED
endpoint.

This term is intended to enter the future model as a bounded additive
log10-radiance correction. It must **not** replace the country-specific
non-spectral intercept. In a later forecasting stage an adoption variable
`f_LED(t)` constrained to [0,1] can move each country from its anchor toward
the selected all-LED plateau without extrapolating the spectral effect
indefinitely.


## Population-density dependence of temporal radiance change

A dedicated `density-temporal-change` stage tests whether the rotation of the
European population-radiance relationship is driven by faster radiance growth
in low-density NUTS3 regions.

Run:

```bash
python src/run_pipeline.py density-temporal-change
```

NUTS3 regions are assigned to fixed median-population-density bands
(<30, 30-100, 100-300, 300-1000, >=1000 people km-2), then a mixed model tests
`year x density_band` with country fixed intercepts and NUTS3 random
intercepts. A second cubic continuous model estimates the temporal trend as a
function of log10 population density and tests whether country-specific year
trends or country-specific density slopes are still required. A country
density-band composition table quantifies how differently countries sample the
rural-to-urban continuum.


## Decomposing country-specific structure after the rural-urban temporal surface

The next diagnostic quantifies how much of the original country-specific
structure is absorbed by the common nonlinear population-density x time
surface.

Run:

```bash
python src/run_pipeline.py density-country-decomposition
```

Three country blocks are tested with the same number of country-specific
parameters before and after adding the common nonlinear density-time surface:

- country x year;
- country x log10 population density;
- country x log10 population density x year (country-specific rotation).

For each block the stage reports the likelihood-ratio statistic and LR
pseudo-R2 under the simple linear population x year reference and after the
common cubic density-time surface is fitted. The diagnostic effect size is:

```text
fraction removed = 1 - residual block LR / original block LR
```

The same calculation is also reported for LR pseudo-R2. This is a descriptive
decomposition of model structure, not a causal proportion.

Outputs:

```text
data/processed/density_country_block_tests.csv
data/processed/density_country_structure_decomposition.csv
data/processed/density_country_model_comparison.csv
data/processed/density_country_structure_decomposition_summary.txt
figures/population_relationship/density_temporal_change/country_structure_explained_by_density_time_surface.pdf
```
