# Direct multiplicative snow/albedo VIIRS calibration — 2021 validation

**Purpose:** repeat the constrained snow/albedo fit with a directly multiplicative observation model, resolving the mismatch between the previous additive training equation and multiplicative application. Uses the previously generated matched spatial **18-region** January–March 2020–2022 cell/brightness-group data. Do not modify historical radiance or PRC records.

## Identifiable components and assumptions

The data provide a 2021 calendar-month VIIRS radiance \(R_{21}\), a reference \(B=(R_{20}+R_{22})/2\) for the **same location/month**, and snow and visible-albedo anomalies relative to those years. Fit directly on observed 2021 VIIRS:

\`\`\`text
R_2021 = B × T × G + error
T = exp(γ0 + γFeb*I_Feb + γCF*(delta_cloudfree/10)
        + γbrightness*log1p(B)/log(11))
G = 1 + βsnow*max(deltaSnow_fraction,0)
      + βalbedo*max(deltaAlbedo_fraction,0)
βsnow, βalbedo >= 0
\`\`\`

\`T\` explicitly carries nuisance/background year-on-year differences; it is **not** called snow or albedo, and it is **not removed** from the candidate surface-only correction \`R_2021/G\`. The positive-only rule is deliberate for an optical snow/albedo amplification and forbids increasing corrected radiance. It is a testable constraint, *not evidence all negative snow anomalies have no optical effect*.

Parameters are fitted by directly optimizing observed-radiance prediction (robust \`soft_l1\` loss on unlogged radiance), not by regressing an additive anomaly and treating the fitted terms as a ratio afterwards. The 2020/22 baseline is an observational proxy, not independently validated unchanged lighting. Region-month weights have equal totals, preventing the largest sampled NUTS3 regions from dominating merely by sample size. \`B>=0.2 nW cm^-2 sr^-1\` prevents effectively unidentifiable gains at zero brightness. Snow and albedo measured at different resolutions may remain correlated. Numeric parameter bounds prevent optimizer explosion and **must be reported** if a fit hits them.

### Four comparable families

- \`background_only\`: \`T\` only, no surface correction (\`G=1\`).
- \`snow_only\`: one nonnegative amplification coefficient.
- \`albedo_only\`: one nonnegative amplification coefficient.
- \`snow_plus_albedo\`: two nonnegative amplification coefficients.

The same quality-filtered observations, including VNP43MA3 visible-albedo valid fraction threshold, are used for all four. Fitting is Jan+Feb only. For each region the model is trained on the other regions, then evaluated in held-out January, February and March of that region. March is a seasonally held-out control (the February indicator is zero, so March inherits January background structure; not a fully fitted March background model).

### Diagnostics and decision criteria

Report RMSE of *direct observation prediction* (tests consistency of the multiplicative fit) separately from RMSE of \`R_2021/G\` relative to the 2020/22 flank proxy (tests whether the correction helps historical trajectories). **Do not select/fit coefficients to the latter statistic.** Add per-region results, median/upper-tail gain values and proportion \(G>2\), independent-region coefficient variation, and 4 baseline brightness strata: 0.2–1, 1–5, 5–20 and ≥20 nW cm^-2 sr^-1.

Models/parameter stability should be reviewed, especially the number of regions whose fitted snow/albedo coefficients collapse to zero or hit numerical boundaries. Strong multicollinearity may justify choosing albedo-only rather than a combined law. The genuine question is whether the empirically estimated gain is portable across geographic holdouts and physically reasonable on dark/moderate pixels, and whether March behaves like a sensible seasonal control.

**Critical:** background \`T\` and snow gain \`G\` are not fully causally identifiable from three years of observational satellite composites. Real lighting changes, cloud/weather residuals, sensor version shifts and snow-correlated missingness may be confounded. Changing from an additive to multiplicative fitting equation does not solve those problems automatically. Neither snow albedo nor nighttime sky-polarization changes have been measured directly at those locations.

## Run on EcoVis

\`\`\`bash
git pull --ff-only
python src/fit_prc_2021_multiplicative_gain.py --self-test
python src/fit_prc_2021_multiplicative_gain.py
\`\`\`

Data read (no downloads): \`data/processed/historical_prc_2021_pixel_snow_test/matched_snowcell_brightness_groups.csv\` from completed 18-region pilot.

Results written (old files not overwritten): \`data/processed/historical_prc_2021_multiplicative_calibration/\`
- \`multiplicative_summary.txt\`
- \`multiplicative_coefficients_and_folds.csv\`
- \`multiplicative_coefficient_stability.csv\`
- \`multiplicative_holdout_diagnostics.csv\`
- \`multiplicative_brightness_strata.csv\`
- \`multiplicative_region_month.csv\`
- \`multiplicative_heldout_groups.csv\`
- \`multiplicative_holdout_diagnostics.pdf\`

Compare these held-out results with the *previous* additive method's \`surface_gain_summary.txt\` and \`surface_gain_holdout_diagnostics.csv\` noting that the direct multiplicative test excludes baseline radiance below 0.2 and uses a different robust objective, so raw RMSE is **not an exact same-sample comparison**. For a definitive prediction-performance comparison, recompute both methods on an identical selected cohort; this diagnostic is the next step if differences appear meaningful.

**Publication decision:** no annual 2021 PRC correction until coefficient stability, physical plausibility, effective independent surface support and monthly→annual mapping have been established.
