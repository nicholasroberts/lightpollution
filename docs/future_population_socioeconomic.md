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
scientifically plausible maximal fixed-effects structure rather than a sequence
of hand-picked alternative models:

```text
log10(corrected radiance)
~ log10(population density)
* log10(GDP per capita PPS)
* year
* country
+ (1 | NUTS_ID)
```

The NUTS3 random intercept is the repeated-measures structure and is held fixed
throughout simplification.

## Simplification procedure

Fixed-effect comparisons are fitted with maximum likelihood (ML). The
simplification proceeds from highest-order interactions downwards while
preserving marginality:

1. fit the maximal model;
2. test deletion of the highest-order currently removable term(s) by
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
