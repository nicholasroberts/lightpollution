#!/usr/bin/env python3
"""Control-year snow-amplification test using the clean 2020 <-> 2022 contrast.

Why this design
---------------
The 2021 anomaly test uses 2021 - mean(2020, 2022).  Treating 2020 or 2022 as
pseudo-target years against the other two years would put anomalous 2021 into
their reference and mechanically create a mirror image of the 2021 pulse.

The clean available control therefore excludes 2021 entirely:

    radiance contrast = 2020 - 2022
    snow contrast     = 2020 - 2022
    baseline          = mean(2020, 2022)

for Jan, Feb and Mar.  Reversing both differences gives the 2022 - 2020 view
and is algebraically the same independent test, not a second replicate.  Both
orientations are written for audit, but inference is based on one contrast.

Model:
    radiance contrast ~ log1p(baseline radiance) + density band + country
                        + snow contrast
                        + snow contrast x log1p(baseline radiance)

VCMSLCFG is primary; VCMCFG is a product-robustness check.  This stage does not
alter production annual or monthly data.
"""

from __future__ import annotations

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr
import statsmodels.formula.api as smf

from project_config import configured_path, load_config
from test_fixed_nuts2016_monthly_snow_amplification import (
    DENSITY_ORDER,
    MONTH_NAME,
    build_snow_anomaly,
    load_atlas,
)


def spearman_clean(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    z = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(z) < 3 or z["x"].nunique() < 2 or z["y"].nunique() < 2:
        return np.nan, np.nan, len(z)
    r = spearmanr(z["x"], z["y"])
    return float(r.statistic), float(r.pvalue), int(len(z))


def partial_r2(reduced, full) -> float:
    if reduced.ssr <= 0 or not np.isfinite(reduced.ssr):
        return np.nan
    return float(max(0.0, (reduced.ssr - full.ssr) / reduced.ssr))


def build_control_table(
    atlas: pd.DataFrame,
    snow: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    months = tuple(int(x) for x in cfg["months"])
    products = tuple(str(x).lower() for x in cfg["products"])

    required = {"radiance_2020", "radiance_2022"}
    missing = required - set(atlas.columns)
    if missing:
        raise RuntimeError(
            "Monthly atlas is missing raw control-year radiances: "
            + repr(sorted(missing))
        )

    d = atlas[
        atlas["product"].isin(products) & atlas["month"].isin(months)
    ].copy()
    d = d.merge(
        snow[
            [
                "NUTS_ID", "month",
                "snow_2020_pct", "snow_2022_pct",
                "snow_valid_2020", "snow_valid_2022",
            ]
        ],
        on=["NUTS_ID", "month"],
        how="left",
        validate="many_to_one",
    )

    numeric = [
        "radiance_2020", "radiance_2022",
        "snow_2020_pct", "snow_2022_pct",
        "snow_valid_2020", "snow_valid_2022",
    ]
    for col in numeric:
        d[col] = pd.to_numeric(d[col], errors="coerce")

    d["baseline_control_nw"] = (
        d["radiance_2020"] + d["radiance_2022"]
    ) / 2.0
    d["baseline_log1p"] = np.log1p(d["baseline_control_nw"].clip(lower=0))
    d["radiance_diff_2020_minus_2022_nw"] = (
        d["radiance_2020"] - d["radiance_2022"]
    )
    d["snow_diff_2020_minus_2022_pp"] = (
        d["snow_2020_pct"] - d["snow_2022_pct"]
    )
    d["min_control_snow_valid_fraction"] = d[
        ["snow_valid_2020", "snow_valid_2022"]
    ].min(axis=1)
    d["density_band"] = pd.Categorical(
        d["density_band"].astype(str),
        categories=list(DENSITY_ORDER),
        ordered=True,
    )
    d["eligible_control"] = (
        d["radiance_diff_2020_minus_2022_nw"].notna()
        & d["snow_diff_2020_minus_2022_pp"].notna()
        & d["baseline_log1p"].notna()
        & d["density_band"].notna()
        & d["CNTR_CODE"].notna()
        & (
            d["min_control_snow_valid_fraction"]
            >= float(cfg["min_snow_valid_area_fraction"])
        )
    )
    return d


def fit_control_models(
    d: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    models = []
    effects = []
    coefficients = []

    for (product, month), g0 in d.groupby(["product", "month"], observed=True):
        g = g0[g0["eligible_control"]].copy()
        if len(g) < 100:
            raise RuntimeError(
                f"Too few eligible control rows for {product} {month}: {len(g)}"
            )

        mu = float(g["baseline_log1p"].mean())
        sd = float(g["baseline_log1p"].std(ddof=0))
        if sd <= 0:
            raise RuntimeError(f"Degenerate baseline for {product} month {month}.")
        g["baseline_z"] = (g["baseline_log1p"] - mu) / sd

        y = "radiance_diff_2020_minus_2022_nw"
        x = "snow_diff_2020_minus_2022_pp"
        f0 = f"{y} ~ baseline_z + C(density_band) + C(CNTR_CODE)"
        f1 = f0 + f" + {x}"
        f2 = (
            f"{y} ~ {x} * baseline_z "
            "+ C(density_band) + C(CNTR_CODE)"
        )

        m0 = smf.ols(f0, data=g).fit(cov_type="HC3")
        m1 = smf.ols(f1, data=g).fit(cov_type="HC3")
        m2 = smf.ols(f2, data=g).fit(cov_type="HC3")

        int_term = f"{x}:baseline_z"
        joint = m2.wald_test(f"{x} = 0, {int_term} = 0", scalar=True)
        rb, pb, _ = spearman_clean(g["baseline_control_nw"], g[y])
        rs, ps, _ = spearman_clean(g[x], g[y])

        row = {
            "contrast": "2020_minus_2022",
            "product": product,
            "month": int(month),
            "month_name": MONTH_NAME[int(month)],
            "role": "control",
            "n": int(m2.nobs),
            "n_countries": int(g["CNTR_CODE"].nunique()),
            "r2_base": float(m0.rsquared),
            "r2_plus_snow": float(m1.rsquared),
            "r2_full": float(m2.rsquared),
            "partial_r2_snow_main": partial_r2(m0, m1),
            "partial_r2_interaction_after_snow": partial_r2(m1, m2),
            "partial_r2_full_snow_block": partial_r2(m0, m2),
            "snow_coef_at_mean_brightness": float(m2.params[x]),
            "snow_p": float(m2.pvalues[x]),
            "snow_x_logbrightness_coef": float(m2.params[int_term]),
            "snow_x_logbrightness_p": float(m2.pvalues[int_term]),
            "joint_snow_block_wald": float(np.asarray(joint.statistic).squeeze()),
            "joint_snow_block_p": float(np.asarray(joint.pvalue).squeeze()),
            "spearman_baseline_vs_delta": rb,
            "spearman_baseline_vs_delta_p": pb,
            "spearman_snow_vs_delta": rs,
            "spearman_snow_vs_delta_p": ps,
        }
        models.append(row)

        ci = m2.conf_int()
        for term in m2.params.index:
            coefficients.append({
                "contrast": "2020_minus_2022",
                "product": product,
                "month": int(month),
                "month_name": MONTH_NAME[int(month)],
                "role": "control",
                "term": str(term),
                "coef": float(m2.params[term]),
                "std_error_hc3": float(m2.bse[term]),
                "ci_low": float(ci.loc[term, 0]),
                "ci_high": float(ci.loc[term, 1]),
                "p_value": float(m2.pvalues[term]),
            })

        cov = m2.cov_params()
        b0 = float(m2.params[x])
        bi = float(m2.params[int_term])
        for q in (0.10, 0.25, 0.50, 0.75, 0.90):
            base = float(g["baseline_control_nw"].quantile(q))
            z = (math.log1p(max(base, 0.0)) - mu) / sd
            eff = b0 + bi * z
            var = (
                float(cov.loc[x, x])
                + z * z * float(cov.loc[int_term, int_term])
                + 2 * z * float(cov.loc[x, int_term])
            )
            se = math.sqrt(max(var, 0.0))
            zstat = eff / se if se > 0 else np.nan
            p = 2 * norm.sf(abs(zstat)) if np.isfinite(zstat) else np.nan
            effects.append({
                "contrast": "2020_minus_2022",
                "product": product,
                "month": int(month),
                "month_name": MONTH_NAME[int(month)],
                "baseline_quantile": q,
                "baseline_radiance_nw": base,
                "snow_effect_nw_per_pp": eff,
                "ci_low": eff - 1.96 * se,
                "ci_high": eff + 1.96 * se,
                "p_value": p,
            })

    return pd.DataFrame(models), pd.DataFrame(effects), pd.DataFrame(coefficients)


def make_mirror_view(d: pd.DataFrame) -> pd.DataFrame:
    """Write both target-year orientations without pretending they are independent."""
    base_cols = [
        "product", "month", "month_name", "NUTS_ID", "CNTR_CODE",
        "density_band", "baseline_control_nw",
        "radiance_diff_2020_minus_2022_nw",
        "snow_diff_2020_minus_2022_pp",
        "min_control_snow_valid_fraction", "eligible_control",
    ]
    a = d[base_cols].copy()
    a["target_year"] = 2020
    a["reference_year"] = 2022
    a["radiance_target_minus_reference_nw"] = (
        a["radiance_diff_2020_minus_2022_nw"]
    )
    a["snow_target_minus_reference_pp"] = a["snow_diff_2020_minus_2022_pp"]

    b = d[base_cols].copy()
    b["target_year"] = 2022
    b["reference_year"] = 2020
    b["radiance_target_minus_reference_nw"] = (
        -b["radiance_diff_2020_minus_2022_nw"]
    )
    b["snow_target_minus_reference_pp"] = -b["snow_diff_2020_minus_2022_pp"]

    out = pd.concat([a, b], ignore_index=True)
    out["independent_contrast_id"] = "2020_vs_2022"
    out["independent_test"] = out["target_year"].eq(2020)
    return out


def compare_with_2021(
    models: pd.DataFrame,
    processed_root,
) -> pd.DataFrame:
    path = (
        processed_root
        / "fixed_nuts2016_monthly_snow_amplification"
        / "monthly_snow_amplification_models.csv"
    )
    if not path.exists():
        return pd.DataFrame()

    y2021 = pd.read_csv(path)
    y2021 = y2021[y2021["product"].eq("vcmslcfg")][
        [
            "month", "month_name",
            "partial_r2_full_snow_block",
            "snow_x_logbrightness_coef",
            "snow_x_logbrightness_p",
            "joint_snow_block_p",
        ]
    ].rename(columns={
        "partial_r2_full_snow_block": "partial_r2_2021",
        "snow_x_logbrightness_coef": "interaction_2021",
        "snow_x_logbrightness_p": "interaction_p_2021",
        "joint_snow_block_p": "joint_p_2021",
    })

    ctrl = models[models["product"].eq("vcmslcfg")][
        [
            "month", "month_name",
            "partial_r2_full_snow_block",
            "snow_x_logbrightness_coef",
            "snow_x_logbrightness_p",
            "joint_snow_block_p",
        ]
    ].rename(columns={
        "partial_r2_full_snow_block": "partial_r2_control_2020_vs_2022",
        "snow_x_logbrightness_coef": "interaction_control_2020_vs_2022",
        "snow_x_logbrightness_p": "interaction_p_control_2020_vs_2022",
        "joint_snow_block_p": "joint_p_control_2020_vs_2022",
    })

    out = y2021.merge(ctrl, on=["month", "month_name"], validate="one_to_one")
    out["partial_r2_ratio_2021_to_control"] = (
        out["partial_r2_2021"]
        / out["partial_r2_control_2020_vs_2022"].replace(0, np.nan)
    )
    return out


def make_plot(
    models: pd.DataFrame,
    comparison: pd.DataFrame,
    out_dir,
    dpi: int,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    ctrl = models[models["product"].eq("vcmslcfg")].sort_values("month")
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    ax.bar(ctrl["month_name"], 100 * ctrl["partial_r2_full_snow_block"])
    ax.set_xlabel("Control month")
    ax.set_ylabel("Snow block partial R² (%)")
    ax.set_title("Clean 2020 vs 2022 control contrast")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(out_dir / "control_2020_vs_2022_snow_block_partial_r2.pdf", dpi=dpi)
    plt.close(fig)

    if not comparison.empty:
        fig, ax = plt.subplots(figsize=(8.0, 5.5))
        x = np.arange(len(comparison))
        width = 0.36
        ax.bar(
            x - width / 2,
            100 * comparison["partial_r2_2021"],
            width=width,
            label="2021 vs 2020/2022",
        )
        ax.bar(
            x + width / 2,
            100 * comparison["partial_r2_control_2020_vs_2022"],
            width=width,
            label="2020 vs 2022 control",
        )
        ax.set_xticks(x, comparison["month_name"])
        ax.set_ylabel("Snow block partial R² (%)")
        ax.set_title("2021 snow signal versus clean control contrast")
        ax.legend()
        ax.grid(alpha=0.2, axis="y")
        fig.tight_layout()
        fig.savefig(out_dir / "2021_vs_control_partial_r2.pdf", dpi=dpi)
        plt.close(fig)


def write_summary(
    models: pd.DataFrame,
    comparison: pd.DataFrame,
    out,
) -> None:
    with out.open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 2020/2022 MONTHLY SNOW CONTROL TEST\n")
        h.write("=" * 78 + "\n\n")
        h.write(
            "This control deliberately excludes 2021 from both the radiance "
            "and snow contrast. Using 2020 or 2022 against the mean of the "
            "other two years would contaminate the control with the 2021 pulse.\n"
        )
        h.write(
            "Independent contrast: 2020 - 2022. The 2022 - 2020 orientation "
            "is the exact sign-reversed view of the same comparison and is "
            "written for audit only, not counted as a second independent test.\n"
        )
        h.write(
            "Jan, Feb and Mar are all control months; no month is designated "
            "primary or supportive. A statistically non-zero snow response is "
            "not automatically a failure of the control, because snow reflection "
            "is a general physical mechanism. The relevant comparison is effect "
            "magnitude and incremental explanatory power versus the 2021 pulse.\n\n"
        )
        h.write("CONTROL MODELS\n")
        h.write(models.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        if not comparison.empty:
            h.write("\n\n2021 VERSUS CONTROL COMPARISON\n")
            h.write(
                comparison.to_string(
                    index=False, float_format=lambda x: f"{x:.6g}"
                )
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_monthly_snow_control_years"]
    atlas = load_atlas(cfg)

    products = tuple(str(x).lower() for x in cfg["products"])
    months = tuple(int(x) for x in cfg["months"])
    use = atlas[
        atlas["product"].isin(products) & atlas["month"].isin(months)
    ].copy()
    ids = set(use["NUTS_ID"])
    expected_n = int(cfg["expected_monthly_tile_nuts3"])
    if len(ids) != expected_n:
        raise RuntimeError(
            f"Monthly atlas has {len(ids):,} NUTS3; expected {expected_n:,}."
        )

    snow = build_snow_anomaly(config, cfg, ids)
    control = build_control_table(use, snow, cfg)
    models, effects, coefficients = fit_control_models(control)
    mirror = make_mirror_view(control)

    processed_root = configured_path(config, "processed")
    processed = processed_root / "fixed_nuts2016_monthly_snow_control_years"
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_monthly_snow_control_years"
    )
    processed.mkdir(parents=True, exist_ok=True)

    comparison = compare_with_2021(models, processed_root)

    outputs = {
        "regional": processed / "monthly_snow_control_2020_vs_2022_regional.csv",
        "mirror": processed / "monthly_snow_control_target_year_mirror.csv",
        "models": processed / "monthly_snow_control_models.csv",
        "effects": processed / "monthly_snow_control_marginal_effects.csv",
        "coefficients": processed / "monthly_snow_control_coefficients.csv",
        "comparison": processed / "monthly_snow_control_vs_2021.csv",
        "summary": processed / "monthly_snow_control_summary.txt",
    }

    control.to_csv(outputs["regional"], index=False)
    mirror.to_csv(outputs["mirror"], index=False)
    models.to_csv(outputs["models"], index=False)
    effects.to_csv(outputs["effects"], index=False)
    coefficients.to_csv(outputs["coefficients"], index=False)
    comparison.to_csv(outputs["comparison"], index=False)
    write_summary(models, comparison, outputs["summary"])
    make_plot(models, comparison, figures, int(cfg["figure_dpi"]))

    print("\n2020/2022 MONTHLY SNOW CONTROL TEST COMPLETE")
    print(models.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    if not comparison.empty:
        print("\n2021 VERSUS CONTROL")
        print(comparison.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for path in [*outputs.values(), figures]:
        print(path)


if __name__ == "__main__":
    main()
