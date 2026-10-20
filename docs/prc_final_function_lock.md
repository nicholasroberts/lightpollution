# One final check before freezing the empirical PRC transfer

The investigator requested that we stop revisiting tiny differences among biologically equivalent fitted PRC curves. Atmospheric and weather variability, illumination changes and satellite/ground spectral differences dominate such sub-percent parameter perturbations.

**Decision:** adopt the complete, directly geolocated **18-site field-observation-year EOG VIIRS** calibration as the final PRC function, contingent on **one brief sanity check**. Frozen intended coefficients: \`yf=0.1452426021\`, \`y0=0.4558190410\`, \`k=0.1035007080\`, equation \`PRC(R)=yf+(y0-yf)*exp(-kR)\`. Relative to the older 16-site provisional function, the PRC curves are nearly indistinguishable on the biologically relevant scale. This is an empirical radiance-to-receptor contrast relationship, **not** a behavioural detectability threshold.

Do NOT start new iterations on small numerical differences or search model families again. All 18 original spider sample coordinates/dates have now been established and the original spider CSV data have been recovered. Proceed with the long-term **2014–2024 European pixel→PRC→NUTS3 change** study.

**Last sanity check**: \`python src/prc_final_function_sanity_check.py\`

1. Look up Slabhouse original uncalibrated annual EOG 2018/2019 pixel values (corrected EOG may be zero), and \`cf_cvg\` cloud-free observation counts if their companion GeoTIFF is present. Report missing coverage explicitly.
2. Extract PRC at every *original EOG pixel* within all **NUTS2016 NUTS3** polygons for **2014 and 2024 only**, using both the existing provisional function and the confirmed final function. Offset, original raster pixels and geographic area weighting are identical. Summarise continent-wide 2014/2024 PRC, change, changes in NUTS3 deltas and top-gain ranking overlap.

This runs only four NUTS-year extractions, no 11-year rerun, no atlas maps and no behavioural model changes. Outputs are **new files** in \`data/processed/prc_18site_final_sanity_check/\`; the existing \`config/prc_transfer.yaml\` and production history are **not** altered until the results are viewed. Then set \`status: locked\` in \`config/prc_transfer.yaml\`, save the full 18-site fit provenance and coefficients, and rebuild the **2014–2024 PRC atlas once** with the final locked parameters.

Commands:

\`\`\`bash
git pull --ff-only
python src/prc_final_function_sanity_check.py
\`\`\`

Share only \`final_prc_sanity_summary.txt\` and, if anything looks anomalous, \`slabhouse_annual_radiance_and_coverage.csv\`. Don't request further model searches or site-level refits for marginal changes.


## 20 October addendum — exclude Slabhouse and finalise with 17 sites

The investigator has instructed that \`slabhouse.csv\` **not contribute to the EOG observation-year PRC fit** because the EOG original \`average_masked\` annual raster value is exactly 0.0 in both 2018 and 2019 at its recorded coordinates, although the original monthly source has nonzero readings. This is an exclusion from **one calibration regression**, not deletion of the genuine PRC field observation nor retroactive alteration of the published original 18-point fit. No need for further attempts to explain the zero before focusing on the biological study.

Recalculated locally from the investigator's complete uploaded \`spider18_site_eog_corrected_2018_2019.csv\`, using the same multi-start least-squares SSasymp as \`fit_prc_observation_year.py\`: \`n=17\`, \`yf=0.1453159020633674\`, \`y0=0.45612835620720266\`, \`k=0.1037115607084461\`; \`R²=0.938987991\`, residual RMSE \`0.028410631\`, 16 unique measurement-year EOG pixel groups, pixel-blocked CV RMSE \`0.033926526\`. The complete 18 fit was \`yf=0.14524260201531\`, \`y0=0.45581904096102577\`, \`k=0.10350070793972209\`, RMSE \`0.0276111\`. **Maximum absolute difference of the two prediction curves across R=0–70 is 0.00030932 PRC**, negligible relative to other environmental variability.

The **one last native-pixel geographical check** now compares the 18 and 17 response functions directly under the existing EOG 2014 and 2024 rasters using \`src/prc_final_function_sanity_check.py\`. The script self-tests its maximum PRC curve difference, retains the original raw Slabhouse zero audit, writes ONLY new files under \`data/processed/prc_17site_final_sanity_check/\`, and compares area-weighted European PRC endpoints, all NUTS2016 regional absolute change rank agreement and top-20 overlap. Not a new full 11-year reconstruction. Production \`config/prc_transfer.yaml\` and old atlas files stay unchanged until this check is reviewed; then **freeze the 17-site function, rebuild the atlas once, and do not iterate on minute PRC fit differences.**

\`\`\`bash
git pull --ff-only
python src/prc_final_function_sanity_check.py --self-test
python src/prc_final_function_sanity_check.py
\`\`\`

Outputs: \`data/processed/prc_17site_final_sanity_check/final_prc_sanity_summary.txt\` and \`prc_2014_2024_function_sensitivity_nuts3.csv\`.
