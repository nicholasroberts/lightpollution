# Final 2021 snow/albedo correction: monthly-to-annual source bridge

**Status:** PRC receptor-response function **locked** (17 sites). Winter empirical snow/albedo radiance gain fitted on matched January/February 2020–2022 NPP monthly VIIRS observations, 18 pilot NUTS regions, geographically held-out. The independent August snow-free VIIRS test is completed on all 18 pilot regions. The outstanding technical step is whether the **official EOG 2021 annual VNL V2 \`average_masked\` composite** and our V1 NPP monthly VCMSLCFG observations are consistent enough to propagate measured excess snow-associated winter radiance to the annual layer.

**Source methodology:** EOG Annual VNL V2 is derived from monthly cloud-free means weighted by each month's cloud-free observation count. It additionally rejects radiance outliers using a twelve-month median and zero-masks backgrounds according to a multiyear lit threshold. The official annual \`average_masked\` layer **is not necessarily equal** to a simple twelve-month CF-weighted monthly radiance average. Sources: https://eogdata.mines.edu/products/vnl/ and Elvidge et al. (2021), Remote Sensing 13, 922 (https://doi.org/10.3390/rs13050922).

### Completed independent August control

18/18 pilot NUTS3 valid; 139,530/142,549 measured-valid pixels (97.882%) are MOD10CM snow <=1% in August 2020, 2021, 2022. 2,994/3,056 snowcell × brightness groups labelled snow-free. Across these snow-free regions, area-weighted August NPP VCMSLCFG radiance: 2020 3.38054, 2021 3.48968 nW cm^-2 sr^-1, change **+0.10913 (+3.228%)**. Region-level changes: 13/18 brighter, 5/18 dimmer. 2022 August monthly data are **J01** not NPP, and must not enter the 2020–2021 NPP source comparison. Summer albedo QA was not acquired; this is a snow-only control, not a full independent albedo-gain falsification. Do not reinterpret its +0.109 as snow amplification or infer all 2021 seasonal changes were background.

### Stage 1: audit monthly-to-annual consistency, no correction

\`src/prc_2021_annual_monthly_bridge.py\` compares **2020 and 2021** official EOG annual masked radiance on its own native pixel grid against \`sum_month R_m*c_m/sum_month c_m\` from the twelve monthly EOG NPP VCMSLCFG radiance/cloudfree rasters, resampled to the annual grid by nearest neighbour at the exact pixel centre. Restrict raw annual lit >0.2 nW pixels so the official background-zero mask doesn't spuriously dominate goodness-of-fit. Compare 2020 and 2021 separately; record absolute difference, approximate RMSE, frac >20% residual and January/February contribution to overall cloud-free weighting. Works on the **same original 18 selected NUTS3 pilot polygons** and fixed NUTS2016 spatial geography. This is a *source bridge audit* rather than a new model family search.

\`\`\`bash
git pull --ff-only
python src/prc_2021_annual_monthly_bridge.py --self-test
python src/prc_2021_annual_monthly_bridge.py --limit 6
# after the six regions complete:
python src/prc_2021_annual_monthly_bridge.py --limit 0
\`\`\`

Outputs under \`data/processed/prc_2021_annual_monthly_bridge/\`:
- \`annual_monthly_bridge_summary.txt\` — comparison and Jan-Feb cloudfree share
- \`annual_monthly_bridge_by_region.csv\` — detailed 2020 and 2021 NUTS3 sample
- \`source_audit.json\` — all 24 monthly source paths and annual inputs

No snow correction is performed by this stage, and 17-site PRC parameters are unchanged.

### Next stage: conservative January–February snow correction

If official annual and cloudfree-weighted monthly sources agree on lit pixels well enough to interpret a winter-only subtraction, the candidate 2021 annual **per-pixel** excess is

\`DeltaR = [sum_{m=Jan,Feb} cf_2021,m * R_2021,m * (1 - 1/G_m)] / [sum_{m=Jan..Dec} cf_2021,m]\`

with **frozen, heldout-tested snow/albedo** \`G=1+1.626479*max(Δsnow_fraction,0)+3.788066*max(Δvisible_albedo_fraction,0)\`. Require actual January/February MOD10CM and VNP43MA3 QA0/1 visible-albedo data in **all three years** on each cell, matching VIIRS annual native pixels. Do not silently treat unavailable albedo as zero. Modify only annual pixels that are valid lit VIIRS, clip corrected radiance to [0, observed annual], preserve original annual \`average_masked\` zero-background mask and nodata. Avoid double-application of the +2017 0.235 calibration subtraction: the surface correction should act on the same radiance scale used by the fitted G and annual source. Compare "uncorrected" and "Jan-Feb corrected" in separate output directories. Do **not** tweak coefficients to force the 2021 PRC trajectory to meet surrounding years.

**Remaining limitation**: EOG annual V2 rejects outlier months, so simple monthly CF weighting can overstate winter snow's contribution to final annual radiance. If audit differences are substantial, **stop** and establish effective annual month weighting or state that this bridge cannot identify a unique annual correction. Do not lock an unsupported annual-2021 correction by relying only on attractive PRC 2021 continuity.

No existing data or PRC products overwritten.
