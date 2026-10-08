# Legacy population-density / VIIRS radiance analysis (2016–2019)

This directory is a preserved snapshot of the original analysis used before the
current `lightpollution` pipeline was rebuilt. It is retained for provenance
and direct comparison with the new year-specific NUTS / VIIRS workflow.

**Archive status:** historical reference only. Do not edit these archived source
files to make them agree with the current pipeline.

## Scientific question

The analysis examined the relationship between mean VIIRS night-time radiance
and population density at NUTS3 level, separately by country and year.

The central country/year model was:

```
log10(mean radiance) ~ log10(population density)
```

Country-specific slopes and intercepts were calculated independently for 2016,
2017, 2018 and 2019. The original script also fitted a mixed-effects model over
all four years with NUTS3 identifier as a random intercept.

## Archived files

### Code

- `code/6. Pop_radiance_fitting.R`
  - original R analysis script;
  - reads the annual radiance table;
  - derives annual population density from population and area tables;
  - merges radiance and demographic observations by NUTS3 identifier/year;
  - fits country/year log-log regressions;
  - writes the country slope/intercept summary;
  - creates the faceted country plots.

### Input data

- `outputs/2.rad_dataframe.csv`
  - NUTS3 mean radiance values for 2016–2019;
  - columns are NUTS3 ID and annual mean radiance;
  - this table is the radiance input used by the archived R script.

- `NUTS3data/NUTS3populations.csv`
  - historical NUTS population table used by the original analysis.

- `NUTS3data/NUTS3area.csv`
  - NUTS area table used to calculate population density.

### Derived outputs

- `outputs/6.rad_density_relationships.csv`
  - country-specific slope/intercept estimates for 2016–2019;
  - also contains the mean slope and mean intercept across those four years.

- `figures/6.fitting_graphs.pdf`
  - original multi-country fitted relationship figure.

## Original processing logic

The script:

1. reads `outputs/2.rad_dataframe.csv`;
2. reshapes annual radiance values for 2016–2019 into long format;
3. reads the archived population and area tables;
4. calculates population density as population divided by area;
5. merges population density and radiance by `NUTS_ID` and year;
6. fits `lm(log10(av_rad) ~ log10(density))` independently for every
   country/year;
7. saves annual slopes/intercepts and their four-year averages.

The script also fits:

```r
lmer(
  log10(av_rad) ~ log10(density) + date + country +
    country:log10(density) + (1 | NUTS_ID),
  REML = FALSE
)
```

to assess the multi-year relationship while treating NUTS3 region as a random
intercept.

## Important provenance note: 2017 calibration shift

The archived radiance values in `2.rad_dataframe.csv` had already had the
known VIIRS 2017 zero-point/calibration shift corrected **before** this R
script was run. The upstream correction code is not included in the files
available for this archive.

Therefore:

- the archived 2016–2019 fitted relationships are useful as a benchmark for the
  rebuilt pipeline;
- this archive should **not** be interpreted as documenting the implementation
  of the zero-point correction itself;
- any future reproduction of that correction should be separately documented
  in the active pipeline.

## Relationship to the current pipeline

The current pipeline differs intentionally. It uses:

- the NUTS release applicable to each analysis year;
- Eurostat population density directly where available;
- EOG Annual VNL v2.2 `average_masked` products;
- exact polygon overlap;
- spherical physical-area weighting;
- explicit provenance/QC fields;
- reproducible country/year fit summaries.

The archived 2016–2019 country slopes and intercepts are therefore a
**historical benchmark**, not values the new pipeline is required to reproduce
exactly.

## Integrity

SHA-256 hashes for the six original uploaded files are recorded in
`SHA256SUMS`.
