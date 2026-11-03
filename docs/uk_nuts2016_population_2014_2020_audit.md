# UK exclusion from Figure 4: source audit before any cohort change

**Prepared 3 November 2026 on `publication-figures`.** This is a diagnostic, not an approved data replacement or a change to the frozen Figure 4 dataset.

## What we can already establish from the repository

The approved population-linked analyses use a balanced **2014–2020 NUTS2016** cohort.

1. Full historical VIIRS/spider/moth NUTS2016 geography: **1,522** regions, including Britain.
2. Original continuous 2014–2020 complete cohort using `demo_r_d3dens`: **1,317** regions. UK already had **0** regional population observations in this *matched* cohort. See `docs/results_population_prc_visual_2014_2020.md` and `docs/country_mechanism_typology_population_prc_2014_2020.md`.
3. Jan-1 population divided by constant land-area reconstruction: **1,273** of the original 1,317. The **44 further excluded** were HR 14, AL 12, BE 8, IT 8, EE 2. UK was **not** among these 44; it had already disappeared at the earlier matching stage (`docs/results_population_jan1_fixed_nuts2016_2014_2020.md`).
4. The `reconstruct_population_density_jan1_nuts2016.py` script calls `calculate(inputs())` *from the earlier 1,317-region original-density complete cohort*, then performs a **left join** of new Jan-1 counts and frozen land area *onto that prior cohort*. Even if the raw Jan-1 dataset contains UK records, this architecture **cannot add any UK region that was missing from the 1,317 source**. This is a **confirmed structural reason** why the Jan-1 reconstruction alone did not restore UK regions.
5. `src/plot_publication_figure4_regional_contrasts.py` then reads the frozen 1,273-region joined output without any other population source. Its absence of UK rows is therefore a downstream consequence, not a map rendering problem.

The **remaining unresolved upstream question** is exactly why UK regions did not enter the original 1,317 NUTS2016 **`DEMO_R_D3DENS` complete-year matched** population/PRC cohort: missing numeric values, NUTS geographic version/codes, processed-data filtering, or a different upstream stage. We should not guess, nor backdate Brexit as a specific technical mechanism without evidence.

## Evidence that population estimates exist from a contemporaneous official source

**ONS direct source matching our period/reference date:**

- ONS, **31 March 2021**: *Population Estimates for 1 January 2020 on NUTS boundaries (Unified Demographic Collection for Eurostat)*, https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/adhocs/13079populationestimatesfor1january2020onnutsboundariesunifieddemographiccollectionforeurostat . This is specifically the 1 January 2020 NUTS geography file supplied to Eurostat, with an original XLSX link on that page.
- ONS, **9 March 2021**: *Annual population figures for England and Wales 1998–2020*, https://www.ons.gov.uk/aboutus/transparencyandgovernance/freedomofinformationfoi/annualpopulationfiguresforenglandandwales1998to2020 , says UK NUTS regional population estimates were available through 1 January 2020 and were supplied to Eurostat.
- Eurostat `demo_r_pjanaggr3` population on 1 January: https://ec.europa.eu/eurostat/databrowser/view/demo_r_pjanaggr3/default/table?lang=en . Our local raw snapshot (August 2026) is `data/raw/eurostat/demo_r_pjanaggr3_nuts3.json`.
- Eurostat `demo_r_d3dens` annual density: https://ec.europa.eu/eurostat/databrowser/view/demo_r_d3dens/default/table?lang=en . Our local raw snapshot is `data/raw/eurostat/demo_r_d3dens_nuts3.json`.
- Eurostat `reg_area3` land area, metadata: https://ec.europa.eu/eurostat/cache/metadata/en/reg_area3_esms.htm (2018–2020 fields correspond to NUTS2016).

**Crucial comparability rule:** ONS `mid-year` (30 June) population estimates are NOT interchangeable with January-1 population from Eurostat's demographic collection. ONS explicitly documents the alternative methodology and reference date used for UK NUTS reporting until 2020. If using ONS source for missing UK regions, start with the dated January-1 NUTS data and verify all seven annual reference dates, geographic release and land-area denominator, not a mid-year fallback.

## EcoVis root-cause diagnostic — NO downloads, no overwrites

Added `src/audit_uk_population_nuts2016_2014_2020.py` to compare the actual archived numerical values, rather than just presence of UK geographic codes:

- UK NUTS2016 level-3 geometries and all 2014–2020 region-years;
- raw `demo_r_d3dens` density vs the already processed NUTS-version-checked density parquet;
- raw `demo_r_pjanaggr3` Jan-1 population total by NUTS code and year;
- raw `reg_area3` fixed 2018–2020 land area and area stability;
- original 1,317 complete-cohort reconstruction audit and final 1,273 output;
- locked spider/moth pixel-derived PRC regions 2014–2020;
- source trace for the example NUTS3 region **UKK11 (Bristol)**.

From EcoVis:

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/audit_uk_population_nuts2016_2014_2020.py --self-test
python src/audit_uk_population_nuts2016_2014_2020.py
```

The read-only report is saved to `data/processed/publication_figures/uk_population_nuts2016_audit/`:

- `uk_per_year_stage_counts.csv` — number of UK NUTS2016 regions with numeric official density, numeric Jan-1 population, valid fixed 2018–2020 area and locked PRC values, **for each year 2014–2020**.
- `uk_region_first_exclusion_reason.csv` — each UK NUTS2016 region and the first stage of exclusion, with separate **complete 2014–2020** coverage flags for every source.
- `uk_nuts3_year_by_year_source_trace.csv` — full UK × year stage trace with original numeric values and flags.
- `uk_population_cohort_audit_summary.json` — UK cohort counts, independent Jan1+fixed-area+both-PRC eligibility count, codes present in Eurostat JSON dimensions but not NUTS2016, and original source paths.

The audit **does not** run Eurostat downloads, alter raw or processed population sources, write a new 2014–2020 cohort, refit anything, or regenerate Figure 4. Diagnostic only. Error if an archived input is absent, instead of silently replacing it with new data.

## Criteria for including UK next

Only after looking at the actual audit results:

1. Establish the numeric **seven-year** UK Jan1 counts, or obtain the official 1 January NUTS geographic totals for any gaps from ONS/Eurostat.
2. Confirm NUTS2016 code + polygon correspondence and a fixed land-area definition consistent with the existing `reg_area3` source. This may need a documented UK area source if the Eurostat snapshot is incomplete.
3. Validate that UK native-pixel PRC values exist for every region-year (independent of demographic data).
4. Prepare a **new explicitly named expanded matched cohort** (e.g. `..._with_UK`), do not silently overwrite the frozen 1,273 source. Maintain a controlled comparison of all old 1,273 rows and recompute distribution denominators transparently.
5. Only then produce a revised Figure 4 set. The existing figures remain valid descriptive products for their declared 1,273-region cohort until then.

**Current state:** structural reason for lack of UK in reconstructed cohort verified; the original raw-Eurostat-vs-processed reason still awaits EcoVis source-trace output.
