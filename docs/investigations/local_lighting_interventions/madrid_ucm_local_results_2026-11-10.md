# Madrid UCM: observed local-scale VIIRS results — 10 November 2026

**Status:** Investigator executed `src/madrid_ucm_local_ground_bridge.py --viirs` on **actual extracted EOG VIIRS monthly VNL V1 (2015–2019)** and supplied the resulting PNG and three CSVs. This is **observed VIIRS data**, not synthetic plots. The underlying high-temporal-resolution **raw AstMon photometry is not yet processed or quality-validated**; the only ground quantities currently compared are the independently published Robles et al. (2021) Table 6 three-year after-midnight anchors. No new PRC calculation or frozen annual NUTS3 changes.

**Original user-supplied outputs:**
- `madrid_ucm_local_viirs_2015_2019.png`
- `madrid_ucm_3km_10km_2015_2019_monthly.csv`: 360 rows = 2 geographic circles × 60 months × 3 cf thresholds
- `madrid_ucm_fixed_common_pixels_2015_vs_2019.csv`: 72 rows = 2 circles × 12 paired calendar months × 3 cf thresholds
- `madrid_ucm_2015_2019_footprint_ground_comparison.csv`: 3 satellite aggregate % changes (ground photometry not included in this file)

## Principal observational comparison

| Measurement (separate instruments / support) | 2015→19 % change |
|---|---:|
| Frozen ES300 region, annual VNL V2 | **−4.8%** (earlier output; not recalculated in this local CSV) |
| Madrid municipality, average of 12 monthly VNL V1 means | **−1.9431%** |
| UCM circle **10 km**, average of 12 monthly VNL V1 means | **−3.7768%** |
| UCM circle **3 km**, average of 12 monthly VNL V1 means | **−2.6025%** |
| **Ground UCM Johnson B zenith sky brightness** (Robles 2021 Table 6) | **+30.6%** |
| **Ground UCM Johnson V zenith sky brightness** | **+11.7%** |
| Ground UCM Johnson R zenith sky brightness | **−1.8%** |

The satellite decline does not disappear when sampling the immediately surrounding urban area, even though the UCM ground photometer measured markedly greater zenith blue brightness. This is a **strong observational mismatch between measures**, **not proof of a VIIRS spectral error**, nor evidence that *all* Madrid skyglow brightened. VIIRS measures upward DNB radiance from mapped ground areas; AstMon measures wavelength-specific **zenith sky brightness at one ground observatory** with contributions from atmospheric scattering across a much wider, unknown illumination source region.

## Spatial-support integrity

Across **all 60 months 2015–2019**, each **3 km circle had exactly 172 valid/inside VIIRS cells**, and each **10 km circle exactly 1,913 valid/inside cells**. All three tested `cf_cvg` thresholds **1, 3, 5** gave identical selected pixel counts and monthly means; all pixels satisfy `cf_cvg ≥ 5`, including the paired 2015→19 analysis. Median monthly `cf_cvg` over the full sequence did not fall below **8** in either circle. Therefore results cannot be ascribed to changing nominal pixel support or cells only present because a weak cf threshold was used. **Important:** adequate coverage count does not prove atmospheric equivalence; EOG `cf_cvg` is an observation count, not a guarantee of equal overpass timing, aerosol, cloud edge or surface reflectance.

## Longitudinal local monthly means (cf threshold 3)

| Year | 3 km VNL V1 mean | 10 km VNL V1 mean |
|---|---:|---:|
| 2015 | 48.105 | 43.470 |
| 2016 | 45.130 | 40.741 |
| 2017 | 46.843 | 41.712 |
| 2018 | 46.031 | 41.216 |
| 2019 | 46.853 | 41.828 |

Units: nW cm⁻² sr⁻¹. **2015→2016 drop** is visible in both circles (3 km ~−6.2%; 10 km ~−6.3%), followed by partial recovery through 2019. **There is no monotonic 2015–2019 decline.** This matters when comparing to the 2015→19 ground B brightening and independently documented multi-year retrofit; no specific lamp installation date inferred from the satellites.

## Pairwise month-by-month anomalies: season matters

For same-calendar-month, identical-cell comparisons (`cf_cvg≥3`):

- **3 km:** **9/12** months show lower 2019 radiance than 2015; November **−17.32%**, February **−5.39%**, but October **+1.76%**, March **+2.34%** and December **+1.14%**.
- **10 km:** **7/12** months show lower radiance; January **−13.45%**, February **−12.10%** and November **−12.71%** dominate. December **+4.70%**, October **+3.22%**.
- **March–October only** (equal-weight means across same eight months): **3 km −0.83%**, **10 km −0.88%**. **Jan–Oct**: 3 km **−1.51%**, 10 km **−3.72%**. Excluding **November** only: 3 km **−1.27%**, 10 km **−2.99%**.
- The local VIIRS decline is **small and seasonally concentrated**. Do not represent it as a persistent year-round decline or a measured year-by-year change in skyglow. Some same-month differences may reflect changing snow/weather/surface reflectance or other independently operated lighting, even though all pixels have high `cf_cvg`.
- The absolute source means and the published ground B magnitudes have **different units and viewing directions**, so report **relative changes within each band/support** only. **2015 itself overlaps the HPS→LED transition**; it is not an untreated baseline.

## Next validation, not yet completed

**Verified real-data blocker (10 Nov 2026):** EcoVis retrieved the original Zenodo AstMon_evol.csv and verified MD5 \`634f098973251bf41f7d8a12b0a8ab45\`. Actual fields: \`Julian_Dates\`, \`Datetime_Raw\`, \`Datetime_UTC\`, \`Filter\`, \`Pixel_x\`, \`Pixely\`, \`Pos_0\`, \`Magnitude\`, \`Error\`; there is **no published per-observation cloud/moon flag**. \`Filter\` denotes Johnson B/V/R band, NOT sky quality. The original paper manually excluded cloudy and moonlit observations, retaining 119,296 readings (~28%). We cannot truthfully reconstruct their selected monthly photometry from this CSV alone. No quality mapping/config should be fabricated.

**New directly published ground observation comparison, Robles Tables 3–5:** Winter at approximately **01:30**, paper's **Nov–Jan** season, 2015→2019 Johnson B **19.30→18.97 mag arcsec⁻² = +35.52% sky brightness**, V **18.07→17.96 = +10.66%**, R **17.42→17.45 = −2.73%**. Ground observations in these tables were filtered by the authors. In the corresponding local VIIRS CSV, a separately defined **same-calendar-year Jan, Nov, Dec** window yields **UCM 3 km −5.92%** and **UCM 10 km −7.72%**. Both results computed with actual local raster means at \`cf_cvg>=3\`; complete per-month pixel support. **Season year bins are not precisely matched**: the author's winter can span December/January, VIIRS's Jan/Nov/Dec set is within a calendar year. These are independent winter-centred measurement contrasts, **not simultaneous 01:30 satellite observations**, and do not identify a spectral-only causal mechanism.

[New published-winter versus VIIRS plotting script](../../../src/plot_madrid_ucm_published_winter_sky_brightness.py) uses existing UCM-local monthly VIIRS CSV and the authors' Tables 3–5 without ever reading unfiltered AstMon raw data.

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures
python src/madrid_ucm_local_ground_bridge.py --self-test
python src/madrid_ucm_local_ground_bridge.py --inspect-ground
python src/plot_madrid_ucm_published_winter_sky_brightness.py --self-test
python src/plot_madrid_ucm_published_winter_sky_brightness.py
\`\`\`

**Next source priority:** obtain publisher's [Supplementary Tables S1–S3](https://www.mdpi.com/article/10.3390/rs13081511/s1) for the original authors' quality-filtered **annual B/V/R** medians, or the original authors' per-observation manual flag annotations for quality-controlled monthly ground photometry. The latter are NOT present in the Zenodo download, so earlier instructions to create \`astmon_columns_flags.json\` from it were incorrect.


**Scientific scope:** Our key finding so far is that **satellite VIIRS and blue-sensitive ground sky-brightness changes can diverge even at local spatial support**, not a uniquely identified spectral artifact nor a measured polarization receptor contrast. Additional independent spectral observations and emitting-source timelines remain necessary.
