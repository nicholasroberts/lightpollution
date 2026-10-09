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

A harmonized annual country-by-country dataset of LED street-light penetration
could not be identified at sufficient quality for 2013-2023. Rather than use a
weak sales/trade proxy as if it were installed lighting stock, the branch uses
an independent spectral validation based on Sánchez de Miguel et al. (2022).

The published calibrated ISS mosaics compare Europe in 2012-2013 with
2014-2020 and quantify spectral composition using B/G and G/R ratios. The
processed RGBA mosaics are available openly from Zenodo record 7677478.

Run:

```bash
python src/run_pipeline.py iss-spectral-validation
```

The stage downloads the two ~421 MB calibrated mosaics, extracts country-level
median B/G and G/R before and after the lighting transition, and applies the
same published colour-ratio quality limits (B/G <= 1.2, G/R <= 1.2 and
R/G <= 6).

For 2013-2020 NUTS3 observations it then fits a control mixed model:

```text
log10(corrected VIIRS radiance)
~ log10(population density)
+ log10(GDP per capita PPS)
+ population × GDP
+ country
+ (1 | NUTS_ID)
```

Country-specific temporal slopes are estimated from the residuals. These slopes
are then compared with the independent ISS spectral shifts. The mechanistic
prediction is that stronger shifts toward blue/white lighting should be
associated with more negative VIIRS-band residual trends because VIIRS DNB has
limited sensitivity to the blue LED emission peak.

This is a validation of the lighting-technology explanation, not an annual LED
penetration covariate and not yet part of the future projection model.
