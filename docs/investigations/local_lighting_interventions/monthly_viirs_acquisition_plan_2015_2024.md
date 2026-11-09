# Historical EOG V1 monthly radiance acquisition — five local lighting interventions

**Prepared 9 November 2026.** This is an **acquisition plan and scripts audit** based on GitHub `publication-figures` plus EcoVis file counts supplied by the investigator. No remote EcoVis filesystem inspection has occurred in this turn. The new read-only inventory and manual archive importer have been committed but still need to be **run on EcoVis**. Locked European annual radiance and PRC functions remain unchanged.

## Known existing EcoVis input inventory

User previously counted complete `avg_rade9*` plus `cf_cvg` pairs for **both** `VCMSLCFG` and `VCMCFG` in each month **2020, 2021 and 2022**. A 2023-01 hit was misleading: `202301` was also found in an EOG **creation timestamp** `_c202301...` on a 2022-12 file. The subsequent local pilot extracted **2020–22** Échirolles and **2020–21** Portsmouth from existing files, but no 2023 monthly raster. Annual NUTS2016 parquet and annual `average_masked` are *not* substitutes for the EOG monthly average and cloud-free observation layers.

### Requested windows and why

| Locality | Event / objective | Complete available years | Missing priority months | Optional refinement |
|---|---|---|---|---|
| **Échirolles (FRK24)** | Oct **2022** road switch-off; 12-month centred display Apr 2022–Mar 2023 | 2020–22 | **2023-01, 2023-02, 2023-03** | Rest of 2023 and 2024 to test persistence alongside LED rollouts; important to flag when LED conversion begins |
| **Portsmouth (UKJ31)** | Summer **2020** adaptive dimming, with need for seasonal baseline | 2020–22 | **2019-01 to 2019-12**, same months before intervention | 2017–18 and 2019 for staged LED retrofit |
| **Coventry (UKG33)** | May–Jul **2024** part-night rollout | 2020–22 | **2023-01 to 2024-12** | 2022 already in hand as a second baseline. May–July 2024 is a transition period; compare Jan–Apr / Aug–Dec, seasonally adjusted |
| **Bristol (UKK11)** | LED/CMS programme started Mar **2023**, continued in 2024+ | 2020–22 | **2023-01 to 2024-12** | Already have 2022 baseline. Do not treat 2024 as rollout completion |
| **Thanet (UKJ44)** | First LED/all-night phase Dec **2016**–May **2017**, subsequent 2018+ phase | 2020–22 (too late for first event) | **2016-01 to 2019-12**, with 2015 seasonal baseline desirable | **2015-01 to 2015-12**, especially to anchor 2016 before the conversions; may restrict 2019 if actual phase-2 completion is well dated |

**Efficient shared download coverage:** The EOG **75N060W** northern tile covers all five. Therefore downloading a month **once** serves every city.

- **Priority download batch (36 months):** full **2019, 2023 and 2024** VCMSLCFG. This immediately enables Portsmouth 2019 context, completes Échirolles early-2023 post period, and covers Bristol and Coventry.
- **Second batch (36 months):** full **2016–2018** VCMSLCFG for Thanet and Portsmouth's earlier LED history.
- **Third optional batch (12 months):** **2015** VCMSLCFG as Thanet's pre-conversion seasonal reference.
- Total for all five with 2015 baseline: **84 new VCMSLCFG month archives** (7 years, `2015–2019` plus `2023–2024`), not 84 per city. Approximate historic archive size is several hundred MB per month: plan for **tens of gigabytes** of downloading/extracted data. Check `df -h` before committing to all years.
- **Do not double the download volume yet** by acquiring VCMCFG for these same periods: VCMSLCFG was the more continuously observed series in the 2020–22 pilot, particularly in summer; VCMCFG can be downloaded later for selected seasonal robustness checks. Monthly V1 VCMSLCFG *does* retain stray-light-corrected observations and these may be lower quality; record `cf_cvg` and compare alternate VCMCFG if available.

## Exact EOG files to download

**Source:** EOG monthly VNL **Version 1** (NOT annual VNL V2.2 and NOT nightly data).
- EOG information: https://eogdata.mines.edu/products/vnl/
- EOG archived directory structure established by existing repo scripts:
  `https://eogdata.mines.edu/nighttime_light/monthly/v10/YYYY/YYYYMM/vcmslcfg/`
  (e.g. https://eogdata.mines.edu/nighttime_light/monthly/v10/2023/202301/vcmslcfg/).
- Select the **75N060W** .tgz archive for each observation month, whose actual archive name resembles `SVDNB_npp_20230101-20230131_75N060W_vcmslcfg_v10_c....tgz`. **Use the first date range**, not the `_c...` production timestamp; some months may have unusual end dates and/or satellite designations.
- From each archive import **`avg_rade9h.tif` and `cf_cvg.tif`** (and optionally `cvg.tif` if needed by legacy global diagnostics). `cf_cvg` is the number of accepted cloud-free observations per pixel, *not* cloud fraction. Keep `cf_cvg=0` as missing radiance.
- Source monthly composites must stay consistent across a chosen contrast; **do not use annual `average_masked`** as a missing monthly substitute.

## What the GitHub auto-download history actually proves

The repository has **working code pathways** for EOG downloads, but not a committed audit of whether months were retrieved through direct login versus manual browser downloads.

| GitHub script | Existing scope | Capabilities | Key limitation |
|---|---|---|---|
| `src/fetch_eog_monthly_vcmslcfg_2021.py` | **2021** VCMSLCFG | Direct directory access then OAuth/OpenID fallback; reuse matching manually downloaded `.tgz` in `~/Downloads`; reuse July pilot; extract layers, validate raster CRS/grid, write `manifest.csv` | **Hard-coded to 2021 by config**, not a generic 2015–2024 downloader |
| `src/fetch_eog_monthly_vcmcfg_2020_2022.py` | **2020–2022** VCMCFG | Same EOG auth helper, archive reuse, monthly extraction, validation and per-year `manifest.csv` | Configured only for 2020–22 VCMCFG; not for the newly missing years |
| `src/diagnose_fixed_nuts2016_monthly_product_atlas.py` | **2020–2022** paired monthly product atlas | Local input discovery and validation, not acquisition | Does not download files |
| `src/local_viirs_monthly_pilot.py` | Échirolles and Portsmouth | Reads existing local monthly TIFFs; makes CSV, series and maps | Does not download files; cannot substitute an annual series |

The config sections are `downloads.eog_monthly_vcmslcfg_2021` and `downloads.eog_monthly_vcmcfg_2020_2022` in `config/pipeline.yaml`. Raw data and processed manifests are omitted from Git by `.gitignore`, so **GitHub alone cannot distinguish automated vs manual acquisitions**. Local manifests have `archive_source` values such as `downloaded`, `local_existing`, `existing_extracted`, and are the proper provenance to inspect on EcoVis.

## Scripts added for this acquisition stage

1. [`src/inventory_monthly_viirs_acquisitions.py`](../../../src/inventory_monthly_viirs_acquisitions.py) — **read-only** listing for each month/year/product: both layers extracted? local archive waiting in Downloads? existing automated downloader manifest provenance? Avoids false Jan 2023 matches by parsing the **observation start date**. Saves an inventory CSV in `data/processed/local_lighting_interventions/`.
2. [`src/import_eog_manual_monthly_archives.py`](../../../src/import_eog_manual_monthly_archives.py) — safe local **archive importer**, not a network downloader. Scans EOG `.tgz` already in `~/Downloads`, extracts just `avg_rade9h` and `cf_cvg` to `data/raw/viirs/monthly_vcmslcfg_<year>/tif/`, validates the matching raster grids, and writes `manual_import_manifest.csv`. Existing original `.tgz` archives are **retained**; no source annual data or models are overwritten. `--dry-run` previews what would be imported; `--archive-dir` selects another folder, and `--years` can narrow the selection. **Not yet run on EcoVis.**

### EcoVis workflow (existing GitHub branch and venv)

```bash
cd /run/media/bznwr/data/lightpollution
git status --short
git switch publication-figures
git pull --ff-only origin publication-figures
python src/inventory_monthly_viirs_acquisitions.py
```

If the pull refuses due to a local script or untracked file (as previously happened with `src/local_viirs_monthly_pilot.py`), **do not force it or discard work**; preserve any local file before resolving the conflict.

If monthly archives were downloaded **manually** to `~/Downloads`:

```bash
python src/import_eog_manual_monthly_archives.py --dry-run
python src/import_eog_manual_monthly_archives.py
python src/inventory_monthly_viirs_acquisitions.py
```

The importer **does not fetch files**. It only extracts available matching archives. For another download directory use `--archive-dir /your/path`. For the secondary product, use `--product vcmcfg` when needed.

### Inspect actual provenance

```bash
find data/raw/viirs -name 'manifest.csv' -o -name 'manual_import_manifest.csv'
```

The inventory command summarizes how many manifest entries record `archive_source=downloaded` (the downloader fetched an archive), `local_existing` (previously downloaded archive reused), `existing_extracted`, or `manual_import`. If a month was extracted manually without a manifest it may appear as complete in the inventory but its acquisition method cannot be established reliably.

## Scientific interpretation safeguards

- **Seasonally matching** each intervention month against one or more prior years is more robust than a straight pre/post six-month average; pooled monthly maps still need weather and `cf_cvg` checks.
- Portsmouth's 2020/2021 contrast is **post/post**; 2019 monthly data are an important improvement. Échirolles 2020/2021 are pre-event references, with 2023 months essential to test persistence into winter.
- 2021 snow/albedo and October 2022 missing MODIS monthly snow are known caveats; the 2022 Échirolles spatial decline is promising but must still be checked against nearby controls and atmospheric/surface conditions.
- Once 2023–24 arrive, the case-study script must be extended with **municipality boundaries** for Bristol, Coventry and Thanet. Those follow-up GIS masks are a separate step; downloading data alone does not automatically add their plots.

## Actual EcoVis inventory update — 9 November 2026, 13:26 UTC

**User-supplied inventory output** from `python src/inventory_monthly_viirs_acquisitions.py` (not a remote filesystem scan by the assistant). Previously known `2020–2022` files remain **extracted and ready**, 12/12 monthly pairs in each year for both products.

A further **45 VCMSLCFG V1 `75N060W` monthly archives have been located but are not yet extracted**. They span: 2015 (12), 2016 (Feb/Aug/Sep/Oct/Nov/Dec = 6), 2017 (Apr/May/Jun/Aug/Sep/Oct = 6), 2019 (Jan/Feb/Jun/Aug/Dec = 5), 2023 (Jan–Apr/Sep/Oct = 6), 2024 (Jan–Jun/Sep–Dec = 10). No 2018 archive was found.

**Exactly 39 more VCMSLCFG archives need downloading** to fill the originally targeted 2015–19 + 2023–24 coverage:

| Year | Need to download (archive not located) | Count |
|---|---|---:|
| 2015 | none | 0 |
| 2016 | Jan, Mar, Apr, May, Jun, Jul | 6 |
| 2017 | Jan, Feb, Mar, Jul, Nov, Dec | 6 |
| 2018 | Jan–Dec (all 12) | 12 |
| 2019 | Mar, Apr, May, Jul, Sep, Oct, Nov | 7 |
| 2020–22 | none — extracted/ready | 0 |
| 2023 | May, Jun, Jul, Aug, Nov, Dec | 6 |
| 2024 | Jul, Aug | 2 |
| **Total** | | **39** |

**Immediate acquisition priority: 15 archives** = 2019 missing seven + 2023 missing six + 2024 July/August. These complete the critical Portsmouth/Échirolles 2019/2023 windows and 2023–24 Coventry/Bristol chronology; Thanet and earlier Portsmouth retrofit studies then need the remaining 24 missing 2016–18 archives. All **12 archives for 2015 are already located**, useful as Thanet's pre-conversion seasonal baseline.

**Next workflow: extract before any further downloading.** The inventory counts `archives awaiting extraction`, not ready-to-use `avg_rade9h + cf_cvg` raster pairs. Run:
```bash
cd /run/media/bznwr/data/lightpollution
python src/import_eog_manual_monthly_archives.py --dry-run
python src/import_eog_manual_monthly_archives.py
python src/inventory_monthly_viirs_acquisitions.py
```
The importer defaults to `~/Downloads` only; if some archives found by the inventory actually live under `data/raw/viirs/.../archives/`, use `--archive-dir PATH` separately for each such directory. Do not re-download files already found as archives. No VCMCFG acquisition is needed at this point. Observe `archive_source` in manifests to determine provenance. None of these steps alter the locked annual European data/PRC functions.
