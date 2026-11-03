# Figure 4A–D — frozen demographic change / PRC change graphics

**3 November 2026 — branch `publication-figures`.** Read-only figure production from the established, approved datasets, with editable Illustrator SVG and vector PDF. Rendering and actual data QA occur locally on EcoVis.

## Scientific inputs — no reanalysis

This is the **2014→2020** demographic comparison. Do NOT treat it as the 2014→2024 all-Europe radiance/PRC study, where no equivalent complete Jan-1 population-density series is available.

1. **Main 1,273 matched NUTS2016 regions:** `data/processed/jan1_population_radiance_prc_country_nl333/matched_region_change_2014_2020.csv` — previously created by `src/analyze_jan1_country_radiance_prc_nl333.py` and also used by Figure 2. Columns used: signed absolute 2014→2020 `change_density_jan1_fixed_land_2014_2020` (people/km²), `cumulative_delta_spider_prc_2020`, `cumulative_delta_moth_prc_2020`, density endpoints, NUTS ID, country, radiance endpoints, fixed land area. No new fits, country assignments or recalibration.
2. **Figure 4D original Jan-1 annual density** (not legacy Eurostat official-density series): `data/processed/nuts2016_jan1_population_density_sensitivity/regional_annual_density_2014_2020_full_audit.csv`, field `density_jan1_fixed_land`, 2014–2020, including eligibility check.
3. **Figure 4D original final taxon-specific PRC annual series**: `data/processed/historical_prc_2014_2024_spider_snow_corrected/prc_nuts3_annual_cumulative.csv` and moth analogue. Final locked pixel-derived spider/moth PRC annual means, not PRC of NUTS regional mean radiance. 2021-only snow adjustment does not affect 2014–2020.
4. NUTS2016 level 3 geographic boundaries from the already configured `historical_prc_pixel_nuts3.load_nuts`. No geometry affects calculations; optional simplification is vector rendering only.

Program `src/plot_publication_figure4_regional_contrasts.py` refuses a changed 1,273 cohort, duplicate IDs, nonfinite source values, invalid density endpoints or inconsistent previously frozen population endpoint differences. In 4D it additionally verifies that the selected region's **2014 and 2020** annual density and PRC values reproduce the frozen region endpoints. The original CSV inputs are never rewritten.

## Exactly what each plot contains

**4A — two 2×2 classification matrices.** One spider and one moth; population increasing/decreasing on x, PRC increasing/decreasing on y. Each square gives count and percentage of **all 1,273 matched regions**; any exact-zero/neutral regions counted separately rather than silently assigned a sign. Colour assignment used identically in 4B and 4C:

| Category | Meaning | Colour |
|---|---|---|
| `growth_prc_loss` | Population increases, PRC decreases | orange `#D55E00` |
| `growth_prc_gain` | Population increases, PRC increases | green `#009E73` |
| `decline_prc_loss` | Population decreases, PRC decreases | purple `#7B61A8` |
| `decline_prc_gain` | Population decreases, PRC increases | blue `#0072B2` |
| `neutral` | Exactly zero in either variable, or within *explicit chosen* neutral deadband | grey `#989898` |

**Classification default: strictly positive / strictly negative (zero cutoffs).** These are descriptive signs, not claims of significant/meaningful changes. Existing notes had sometimes used a sensitivity convention of **5 people/km² and 0.001 PRC** to suppress small changes. For a separately identified sensitivity run, pass `--pop-deadband 5 --prc-deadband 0.001`. This changes the **display classification**, not original measurements, and will naturally make many regions grey. The exact thresholds, counts and classifications are exported in JSON/CSV. Do not quietly replace default signs with threshold-significance interpretations.

**4B — two matched maps:** spider and moth four-colour sign classes for 2014–2020. Fixed NUTS2016, matched 1,273 source records, continental drawing window `[-12,34,46,73]` (includes most of Türkiye; different from Figure 1's narrower eastern crop), common category colours and shared legend. Polygons outside the matched population cohort are light grey, **distinct from neutral-in-cohort grey**. Display-only 150 m EPSG:3035 simplification; all data remain original.

**4C — selected country comparisons:** for the **same previously used six countries** `TR, DE, FR, IE, IT, NL`, horizontal 100%-stacked shares of *regional NUTS3 classification* for spider and moth, n per country shown. Each country contributes its own NUTS3 observations; the visual is not an area-weighted national change, fitted country effect or causal grouping. The figures are descriptive, and the underlying country counts and percentages are exported.

**4D — three NUTS3 annual examples, 2014–2020:** France, Türkiye and **NL333 Delft en Westland** (already independently studied). Each row contains **two separate graphs**, one for 2014-relative change in Jan-1 population density and one for changes in spider/moth PRC. Same observations and fixed baseline, no dual axes. The default French NUTS3 is chosen among *positive density and positive spider PRC* cases and the Turkish one among *positive density and negative spider PRC* cases, **closest to the conditional bivariate medians by robust IQR-scaled absolute distance** (rather than choosing extreme cases). These are transparent graphical examples only, not claims to be representative of those countries or evidence of mechanism. IDs and full selection criteria saved in `fig4_provenance.json`. Any given NUTS2016 codes may be substituted without modifying source code, e.g. `--fr-id FRxxx --tr-id TRxxx`; real codes must exist in the same frozen 1,273 matched cohort. `--nl-id NL333` is the default, and may also be changed. The optional user-selected IDs take priority.

The 4D trajectories use **all seven measured annual values 2014–2020**, not just two endpoints or inferred/interpolated intermediate years.

## How to run on EcoVis

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python src/plot_publication_figure4_regional_contrasts.py --self-test
python src/plot_publication_figure4_regional_contrasts.py
```

Optional *separately labelled* threshold sensitivity:

```bash
python src/plot_publication_figure4_regional_contrasts.py \
  --pop-deadband 5 --prc-deadband 0.001 \
  --out figures/publication/figure4_regional_contrasts_sensitivity
```

No other analysis should be rerun for Figure 4. If an existing input is absent, inspect where the already-generated EcoVis file is stored before regenerating anything.

## Output list

Directory `figures/publication/figure4_regional_contrasts/` (SVG and PDF for each):

- `fig4A_spider_moth_four_quadrant_counts.svg` / `.pdf`
- `fig4B_europe_nuts3_sign_categories.svg` / `.pdf`
- `fig4C_country_category_proportions.svg` / `.pdf`
- `fig4D_fr_tr_nl333_regional_trajectories.svg` / `.pdf`

Additional audited files:

- `fig4_nuts3_signed_category_audit.csv` — raw population and spider/moth changes plus 1,273 sign classifications.
- `fig4A_category_counts.csv` — every category and taxon, counts and percentage.
- `fig4C_country_category_proportions.csv` — category counts/percentages for each country and taxon.
- `fig4D_selected_regional_annual_values.csv` — real 2014–2020 population density and both PRC series.
- `fig4_provenance.json` — frozen source paths, thresholds, complete matching cohort, country choices, region selection and display properties.

## Quality checks before final production

1. Confirm **1,273** classified NUTS3, five mutually exclusive categories, country stacks adding to 100%, and both taxa retaining matched NUTS IDs. The map itself crops overseas/outlying regions; this is *visual only*.
2. Confirm 4D's region IDs in `fig4_provenance.json` are sensible for the intended presentation; override to investigator-selected examples if desired before the final manuscript version.
3. Confirm the 4D 2014–2020 series matches frozen 2014 and 2020 Jan-1 and PRC endpoints. A mismatch raises an exception and is not reconciled by replacing source data.
4. Inspect the vector map path quality, editability of labels in Illustrator, category contrast and separate absence/neutral legend entries.
5. Confirm 4B includes Türkiye where expected and that its geographic extent deliberately differs from the narrower Figure 1 source crop.
6. All plotting artifacts are ignored in Git; code, notes and sources references are committed.

**Status at commit:** plotting script written and checked against the existing code/column definitions on GitHub, but the user's full EcoVis datasets are not executable here, so rendered figures cannot yet be claimed as inspected.
