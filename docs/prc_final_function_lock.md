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
