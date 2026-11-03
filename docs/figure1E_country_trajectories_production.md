# Figure 1E — seven-country radiance change, 2014–2024

**3 November 2026. Branch `publication-figures`. Prepared script; full scientific run must be done on EcoVis.**

## Task and panel

Add a NEW **Figure 1E** (do not overwrite Figure 1C) showing **seven fixed-country area-weighted NUTS3 radiance change trajectories** on a *single set of axes*, 2014–2024. The exact seven countries are those used in UK-expanded Figure 2D:

**TR Türkiye, DE Germany, FR France, IE Ireland, IT Italy, NL Netherlands, UK United Kingdom.**

Outcome plotted is **signed change in the country area-weighted annual mean calibrated radiance from its own 2014 baseline**:

```
mean_R_country_t = SUM_i(valid_area_km2_2014_i * corrected_R_NUTS3_i_t)
                    / SUM_i(valid_area_km2_2014_i)

Delta_R_country_t = mean_R_country_t - mean_R_country_2014
```

All series start at **zero in 2014** by construction. Vertical axis has **absolute radiance difference** units `nW cm-2 sr-1`, *not percentage*, to avoid instability in nearly dark regions. A companion CSV includes the actual country annual mean R and the derived percentage change for examination. This figure is a radiance map-derived observation; it is not population density, receptor contrast, or a causal test.

**Scientific cohort:** the same available fixed-NUTS2016 annual radiance regions and same fixed 2014 observed area weights as Figure 1C, grouped by `CNTR_CODE`. It is **NOT** restricted to the smaller Jan1 population-matched 1,452-region Figure 2D cohort; the comparison to Figure 2D is by **country names**, not by exactly the same NUTS3 sample. A common region cohort is required across all 11 years, so within-country denominators and area support are time-invariant.

## Why Figure 1C starts near 1.065 rather than zero

The uploaded Figure 1C PDF's first point in 2014 is approximately **1.065 nW cm-2 sr-1**, because Figure 1C plots **absolute area-weighted mean VIIRS radiance**, not change from baseline. It is

`SUM_i(A_i * R_i_2014) / SUM_i(A_i)`,

where `A_i = valid_area_km2_2014` is the 2014 VIIRS-observed area for each of the balanced 1,522 NUTS2016 NUTS3 regions. The exact numeric value will be printed by the new script and should reproduce the original Figure 1C support CSV to within tiny rounding error if the source cohort is unchanged.

Area weighting is preferable for a *geographic surface average* because NUTS3 areas vary widely. Equal-NUTS mean would give one tiny urban area the same vote as a huge rural NUTS3 polygon. Fixed valid-2014-area weights allow a stable average across years, instead of producing spurious movements due to changes in weights. It is **not** population weighting and does **not** represent average exposure for a person. It approximates the VIIRS-land-area mean provided the annual valid-pixel support is sufficiently comparable; the source and coverage QC matter.

### Why 2021 is not dropped or smoothed

Figure 1C has a noticeable 2021 spike despite applying the **approved locked 2021 native-pixel snow correction**. That source correction addresses an observed snow/surface effect, but **is not a licence to smooth, impute or force 2021 between 2020 and 2022**. Figure 1E keeps 2021 and uses the identical locked correction. It will help show which countries contribute to the residual 2021 pattern.

## Correct sources and provenance

Uses `src/plot_publication_fig1e_country_radiance_trajectories.py` (new; only creates Figure 1E outputs). As in existing Figure 1C, it reads

1. `figures/prc_radiance_spider_moth_2014_2024/nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv` — 1,522 original NUTS3 2014/2024 endpoints and `valid_area_km2_2014`; independent endpoint check.
2. `data/processed/fixed_nuts2016_viirs_2014_2024.parquet` — original `radiance_mean_corrected` for all 2014–2024; 2017+ pixel-level baseline offset already included; never subtract a second time.
3. `data/processed/publication_figures/fig1_2021_locked_snow_radiance.csv` and verified source metadata — via the existing `get_locked_snow_2021` routine, under the already **locked native-pixel** 2021 snow correction. The script checks original annual 2021 radiance against the original calibrated value from the same locked source before inserting the patched radiance by NUTS3 ID. If the cache is missing or stale, the **approved** native-pixel extraction will be performed and cached, without changing any original data.
4. Existing NUTS2016 geometry `historical_prc_pixel_nuts3.load_nuts` for region-to-country mapping. It is never reprojected or re-extracted for the statistical aggregate (no new radiance zonal computation).

The script explicitly regenerates the **Europe** weighted 11-year history from the same retained region rows and checks it against the existing `figures/publication/figure1_radiance/fig1C_europe_annual_radiance_support.csv`, if present. It **errors on mismatch** rather than quietly changing source methods. If the original Fig1C CSV is absent, it saves its new independently computed Europe reference and labels that QA as unperformed.

No source files, original figure panels, raw rasters, or model parameters are modified.

## Run on EcoVis

```bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only
source .venv/bin/activate

python -m py_compile src/plot_publication_fig1e_country_radiance_trajectories.py
python src/plot_publication_fig1e_country_radiance_trajectories.py --self-test
python src/plot_publication_fig1e_country_radiance_trajectories.py
```

## Output files

All new outputs are added to the existing Figure 1 publication folder:

`figures/publication/figure1_radiance/`

- `fig1E_seven_country_radiance_change_2014_2024.svg` — Illustrator-editable vector master.
- `fig1E_seven_country_radiance_change_2014_2024.pdf` — vector PDF.
- `fig1E_seven_country_radiance_change_2014_2024.png` — white 16:9 talk image at 240 dpi (can change with `--png-dpi 300`, or disable `0`).
- `fig1E_seven_country_annual_radiance_changes.csv` — all 7 countries × 11 years: fixed area/region count, absolute weighted radiance, signed ΔR and percentage vs 2014.
- `fig1E_europe_reference_area_weighted_qa.csv` — independently recomputed entire-cohort weighted radiance per year, including exact 2014 value.
- `fig1E_seven_country_trajectories_provenance.json` — country cohort, fixed weight, locked correction audit, baseline, endpoint matching and comparison to Figure 1C.

The script prints numeric **2014 Europe mean**, the independent match to the existing Figure 1C, and each country 2014/2020/2024 mean radiance/ΔR at end. Generated plots are ignored by Git; script + production notes are committed on `publication-figures`.

**Render and inspect the plot in EcoVis/Illustrator before using in the talk.** If the independent Fig1C comparator disagrees, do NOT relax QA tolerances or overwrite the original figure; inspect the source cohort and 2021 provenance.


## Percentage-change companion (added 3 November 2026)

The investigator requested Figure 1E **also as a percentage** without replacing the original absolute-change chart. The existing `src/plot_publication_fig1e_country_radiance_trajectories.py` now writes an additional **same-size, same-seven-country, same-colours, same-2021-source** plot:

`fig1E_seven_country_radiance_percent_change_2014_2024.svg/.pdf/.png`

Each point is calculated from already area-weighted country **mean** radiance:

`percent_change_country(year) = 100 × [Rbar_country(year) − Rbar_country(2014)] / Rbar_country(2014)`.

All lines begin at **0% in 2014** and use a **separate country-specific 2014 denominator**, not Europe's 2014 radiance. Crucial estimand distinction: this is the **percentage change of an area-weighted country mean**, NOT the area-weighted mean of individual NUTS3 percentage changes. Since dark NUTS3 regions may have near-zero 2014 baselines, the latter would be unstable and unrepresentative. A country with a zero/nonfinite baseline raises an error instead of defining an arbitrary percentage.

The absolute figure remains `fig1E_seven_country_radiance_change_2014_2024.svg/.pdf/.png`, unchanged in meaning, and the original 7 × 11 values CSV already contains `pct_change_from_2014` and `baseline_2014_mean_nw`, so no new scientific input or sensor/model adjustment is needed.

**Run the same command again on EcoVis after pulling GitHub:**

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only
source .venv/bin/activate
python -m py_compile src/plot_publication_fig1e_country_radiance_trajectories.py
python src/plot_publication_fig1e_country_radiance_trajectories.py --self-test
python src/plot_publication_fig1e_country_radiance_trajectories.py
```

To transfer **only the newly added percentage SVG/PDF/PNG** to the original Mac folder, run from a local Mac terminal (VPN connected, with the same valid SSH host used previously):

```bash
rsync -avh --progress \
  --include='fig1E_seven_country_radiance_percent_change_2014_2024.*' \
  --exclude='*' \
  bznwr@ECOVIS_HOST:/run/media/bznwr/data/lightpollution/figures/publication/figure1_radiance/ \
  ~/Desktop/lightpollution-figures/publication/figure1_radiance/
```

The script still must be run and graphically inspected on EcoVis before those three plot files exist.
