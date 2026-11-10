# Madrid second-stage spatial and temporal comparison — UCM AstMon vs local VIIRS

**Dated 10 November 2026, \`publication-figures\`.** This investigation directly follows the user's actual Madrid four-panel ground/satellite comparison: fixed ES300 annual VNL V2 2015→19 **−4.8%**, municipality monthly VNL V1 mean **−1.9%**, after-midnight UCM Johnson B sky brightness **+30.6%**, V **+11.7%**, R **−1.8%**, with independently documented 2015–16 HPS-to-LED changes. Those aggregate observations do **not** identify an isolated spectral, optical, source or operational cause.

## Research question

Does the discrepancy between Madrid-wide satellite radiance and *single-site* zenith sky brightness persist when VIIRS is measured over **smaller, independently predeclared areas centred on the same observatory**, rather than averaging the whole municipality or ES300?

**Observatory centre from original article:** UCM Observatorio, IAU I86, **40°27′04″N, 03°43′34″W**, i.e. lat **40.451111111**, lon **−3.726111111**. [Robles et al. (2021), §2.2](https://doi.org/10.3390/rs13081511).

**Design fixed before seeing local VIIRS outcomes:**

1. **Three municipality/satellite scales:** municipal mask from the *already verified* original Madrid analysis; **3 km radius** and **10 km radius** full circles centred on UCM (EPSG:3035 metric buffers converted to the EOG geographic grid). Circles are **not clipped to the administrative city boundary** and not physically justified source-contribution kernels; they are explicit sensitivity supports. No substitution of NUTS3 for municipality.
2. **Exactly the same EOG VNL V1 / \`75N060W\` monthly raster and \`cf_cvg\` source**, 2015–2019 only (60 observation months). Source NUTS3 annual 2014–24, 2021 snow cache and PRC remain untouched. Monthly means are computed for each fixed support with *both* \`avg_rade9h\` and \`cf_cvg\` valid and thresholds **≥1, ≥3, ≥5**.
3. **Strict pixel-paired 2015→2019 sensitivity:** for each *same calendar month*, retain only cells valid in both years at chosen threshold. Report mean changes, valid pixel count, baseline radiance and percentages separately for **3 km / 10 km** circles. This is not a causal treatment effect, and the citywide monthly mean remains the original reported value, not a newly refitted layer.
4. **Independent ground data:** [Robles et al. original observation-level dataset, Zenodo v1 DOI 10.5281/zenodo.4633001](https://zenodo.org/records/4633001), file **\`AstMon_evol.csv\` 45.5 MB** with published **MD5 \`634f098973251bf41f7d8a12b0a8ab45\`**. The published article §2.4 says the authors **manually flagged all observations as clear with no moon, cloudy or moon over the horizon**, selecting only 119,296/approximately 425,000 clean data; **it does not specify in the main text the raw CSV field names or integer quality-code semantics**. The file could not be fully read through the present web tools because the web interface enforces a 4MB response limit. Therefore ***no invented column names or flag meanings***.
5. **Fail-closed ground processing:** first download the original *on EcoVis* and print its actual header/quality distributions (\`--download-ground --inspect-ground\`); then create a **reviewed explicit JSON column/quality flag mapping** and call \`--prepare-ground --ground-map <file.json>\`. An empty/unreviewed mapping is refused. Julian date UTC → **Europe/Madrid local time**; consider only **00:00–04:29 local after-midnight** to avoid the known ~04:45–05:30 streetlighting power restoration. For B,V,R separately: take **nightly median magnitude**, followed by **monthly median over nights**, requiring ≥3 nights/month; expose per-month valid observation and night counts. The published paper used a different KDE/nonparametric estimator, so the new monthly series should not be presented as an exact reproduction of Table 6.
6. **Ground brightening transform:** \`100×10^(−0.4×(m_year−m_2015))\`, separately *within each Johnson B, V, R passband*. Compare **2015–19 monthly/annual radiance-index changes** in each satellite spatial support, and independently observed *zenith* sky brightness. Satellite nadir upward radiance ≠ ground zenith sky brightness. No conversion to polarization receptor contrast without separate scattering/polarization modelling.
7. **Period restriction:** Restrict first main test to **2015–2019**; deliberately exclude the anomalous 2021 snow year, so no monthly snow corrections are improvised. A 2014 NUTS3 annual point can be presented as context but not compared directly to the 2015-start V1 monthly product.
8. **Interpretation:** if UCM-local VIIRS decreases or remains flat during ground B brightening, that is a robust observational *spatial-scale* check. It **cannot uniquely separate LED spectral bias, shielding, changes in lamp lumen ratings or time schedules, aerosol/sky conditions, and lights from outside 10 km**. Ground zenith brightness receives scattered light from a much larger area than the chosen circles.

## Published source and acquisition

- [Robles et al. 2021 main article, §2.2 and §2.4 manual quality filtering](https://www.mdpi.com/2072-4292/13/8/1511).
- [Zenodo dataset, exact file record and MD5](https://zenodo.org/records/4633001). Raw download URL: \`https://zenodo.org/records/4633001/files/AstMon_evol.csv?download=1\`. Do not commit the 45.5 MB raw ground CSV to GitHub.
- Article says **~28% of AstMon observation records survive** its visual clear/moonless selection, before the after-midnight restriction. It is inappropriate to parse all ground records as if clear; a 2015/19 discrepancy based on unfiltered raw entries would be unreliable.
- Original paper's timing and ground quality classification: zenith B/V/R measurements at UCM ~5 minute sampling under astronomical conditions; near-01:30 night timing overlaps approximately with VIIRS DNB overpasses, **not matched simultaneous individual nighttime acquisitions**.

## Commands / status

\`\`\`bash
cd /run/media/bznwr/data/lightpollution
git switch publication-figures
git pull --ff-only origin publication-figures

python src/madrid_ucm_local_ground_bridge.py --self-test

# Already extracted EOG 75N060W; read only; compute UCM 3/10 km
python src/madrid_ucm_local_ground_bridge.py --viirs

# Download and verify 45.5 MB original source and inspect actual field schema
python src/madrid_ucm_local_ground_bridge.py --download-ground --inspect-ground

# Once manual, published flag semantics are established from inspection:
python src/madrid_ucm_local_ground_bridge.py --prepare-ground \
  --ground-map data/raw/ground_sky_brightness/robles2021_zenodo4633001/astmon_columns_flags.json

# Generate final local VIIRS + ground brightness figure
python src/madrid_ucm_local_ground_bridge.py --plot
\`\`\`

**The illustrative JSON path is *not* a pre-existing reviewed config.** Do not invent flag values or run the last two commands until the real raw CSV schema and values have been inspected; commit and review the mapping after decoding. The existing published Table 6 figure remains a valid independent comparator meanwhile.

**Projected outputs:** \`data/processed/local_lighting_interventions/madrid_ground_viirs/local_ucm/\` monthly UCM radiance CSV, strict 2015/19 pixel pairs CSV, vetted AstMon monthly/annual brightness CSVs and quality audit JSON. Figures (PNG/PDF/SVG) in \`figures/local_lighting_interventions/madrid_ground_viirs/local_ucm/\`. Local VIIRS-only output is already useful even if the 45MB raw archive is not yet downloaded.

**Status:** analysis design and algorithmic self-tests completed in the current workspace. Real EcoVis VIIRS and raw AstMon file schema still require execution. Avoid stating uncomputed UCM-local percentages or quality-filtered ground correlations as results.
