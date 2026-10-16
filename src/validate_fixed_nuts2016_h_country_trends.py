#!/usr/bin/env python3
"""Test whether the country composition of terrain habitability explains
national differences in VIIRS radiance trends (fixed NUTS2016, 2014-2024).

PRE-SNOW annual VIIRS radiance, with the established 2017+ pixel calibration.
No 2021 surface correction, no density-class redefinition or new raster run.

Biologically motivated question:
Does a country's mix of baseline density and physically habitable terrain
explain its mean radiance trajectory without requiring a separate slope for
each national administrative boundary?

Compare density D, D+H, D+H-by-density, and each with country-specific slopes,
then score on strictly future-year temporal holdouts against persistence.
Interpret the added country terms as residual national circumstances, NOT as
biological mechanisms and NOT automatically justified by improved BIC.

H_core is independent of nominal density; it is not used as density/H.
Longitudinal population change is deferred until harmonised NUTS2016 data
are established. The predictor is a descriptive/short-horizon forecast of
radiance changes from an observed starting radiance, not a causal claim.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_density_temporal_change import lrt
from analyse_fixed_nuts2016_presnow_density_h import (
    load_panel, fit_models, predict_changes, observed_target,
    score_groups,
)
from project_config import configured_path, load_config

STAGE = "fixed_nuts2016_h_country_validation"

D = "log_rad ~ year_c * C(density_band) + C(country)"
H = D + " + h_c * C(density_band) + year_c:h_c"
HX = "log_rad ~ year_c * C(density_band) * h_c + C(country)"
MODELS = {
    "D": D,
    "D_H": H,
    "D_HxD": HX,
    "D_country": D + " + year_c:C(country)",
    "D_H_country": H + " + year_c:C(country)",
    "D_HxD_country": HX + " + year_c:C(country)",
}
# Explicitly non-causal nested comparisons, fit on the *same* 1,514 regions.
COMPARISONS = [
    ("D", "D_H", "Habitability spatial and common temporal block"),
    ("D_H", "D_HxD", "Density-dependent temporal H interaction"),
    ("D", "D_country", "Residual country-specific temporal behaviour"),
    ("D_country", "D_H_country", "Habitability after country-specific trends"),
    ("D_H_country", "D_HxD_country", "H-by-density after country trends"),
    ("D_HxD", "D_HxD_country", "Country trends after full H interaction"),
]

# The two date-anchored validation experiments are identical for every model.
FOLDS = [(2018, [2019, 2020]), (2020, [2022, 2023, 2024])]


def _assert_locked_panel(panel: pd.DataFrame, expected: int) -> None:
    if panel["NUTS_ID"].nunique() != expected:
        raise RuntimeError("Habitability changed the locked 1,514-region cohort.")
    if len(panel) != expected * 11:
        raise RuntimeError("Panel must be a complete 2014-2024 annual series.")
    if panel.duplicated(["NUTS_ID", "year"]).any():
        raise RuntimeError("Duplicate region/year observations.")
    if set(panel["year"].unique()) != set(range(2014, 2025)):
        raise RuntimeError("Unexpected years in the fixed 2014-2024 panel.")
    if not np.isfinite(panel["H_core"].to_numpy(float)).all():
        raise RuntimeError("Missing H_core in locked cohort.")


def nested_evidence(models: dict) -> pd.DataFrame:
    rows = []
    for reduced, full, interpretation in COMPARISONS:
        lr, df, p = lrt(models[reduced], models[full])
        rows.append({
            "reduced": reduced, "full": full, "interpretation": interpretation,
            "lr": lr, "df": df, "p": p,
            "bic_improvement": float(models[reduced].bic - models[full].bic),
        })
    return pd.DataFrame(rows)


def baseline(anchor: pd.DataFrame, years: list[int]) -> pd.DataFrame:
    """Persistence: last observed radiance remains unchanged."""
    base = anchor.copy()
    base["year_c"] = 2019.0 - 2019.0  # overwritten below by caller? no effect
    # Keep model-free baseline; preserve only columns needed by scoring.
    parts = []
    for year in years:
        p = anchor[[
            "NUTS_ID", "CNTR_CODE", "density_band", "H_core",
            "reference_population_density",
        ]].copy()
        p["anchor_year"] = int(anchor["year"].iloc[0])
        p["target_year"] = year
        p["radiance_anchor_nw"] = anchor["radiance_mean_corrected"].to_numpy(float)
        p["predicted_log10_change"] = 0.0
        p["predicted_change_percent"] = 0.0
        p["predicted_change_nw"] = 0.0
        p["predicted_radiance_target_nw"] = p["radiance_anchor_nw"]
        parts.append(p)
    return pd.concat(parts, ignore_index=True)


def predictions_for_fold(panel: pd.DataFrame, models: dict,
                         anchor_year: int, target_years: list[int]) -> pd.DataFrame:
    anchor = panel.loc[panel["year"].eq(anchor_year)].copy()
    parts = []
    for name, model in models.items():
        pred = observed_target(
            panel, predict_changes(model, anchor, target_years)
        )
        pred["model"] = name
        parts.append(pred)
    p0 = observed_target(panel, baseline(anchor, target_years))
    p0["model"] = "persistence"
    parts.append(p0)
    out = pd.concat(parts, ignore_index=True)
    out["train_through"] = anchor_year
    out["horizon_years"] = out["target_year"] - anchor_year
    if out["NUTS_ID"].nunique() != panel["NUTS_ID"].nunique():
        raise RuntimeError("Dropped region(s) in holdout.")
    return out


def test_temporal_predictions(panel: pd.DataFrame):
    """No model is trained using its held-out test years."""
    fit_tables, parts = [], []
    for train_through, years in FOLDS:
        train = panel.loc[panel["year"].le(train_through)].copy()
        models, fit = fit_models(train, MODELS)
        fit["train_through"] = train_through
        fit_tables.append(fit)
        parts.append(predictions_for_fold(
            panel, models, train_through, years
        ))
    preds = pd.concat(parts, ignore_index=True)
    if len(preds) != len(MODELS.keys() | {"persistence"}) * len(panel["NUTS_ID"].unique()) * sum(len(t) for _, t in FOLDS):
        raise RuntimeError("Unexpected holdout prediction count.")
    return (
        pd.concat(fit_tables, ignore_index=True),
        preds,
        score_groups(preds, ["model"]),
        score_groups(preds, ["model", "train_through", "target_year", "horizon_years"]),
        score_groups(preds, ["model", "density_band"]),
        score_groups(preds, ["model", "CNTR_CODE"]),
    )


def endpoint_full_fit(panel: pd.DataFrame, models: dict) -> pd.DataFrame:
    """Full-period fit is an explanatory reconstruction, NOT a holdout."""
    anchor = panel.loc[panel["year"].eq(2014)].copy()
    parts = []
    for name, model in models.items():
        p = observed_target(panel, predict_changes(model, anchor, [2024]))
        p["model"] = name
        parts.append(p)
    return pd.concat(parts, ignore_index=True)


def country_composition(panel: pd.DataFrame, endpoints: pd.DataFrame):
    """Does geographic H composition account for national slope differences?

    Use equal-region country mean log10 changes per year. Compare country
    residuals after D to those after D+H, separately with/without country
    trends. A positive reduction in squared residuals is descriptive
    (in-sample) evidence, not causal mediation or independent validation.
    """
    country_meta = (
        panel.loc[panel["year"].eq(2014)]
        .groupby("CNTR_CODE", observed=True)
        .agg(
            n_nuts3=("NUTS_ID", "nunique"),
            H_mean=("H_core", "mean"),
            H_median=("H_core", "median"),
            H_q25=("H_core", lambda x: float(x.quantile(.25))),
            H_q75=("H_core", lambda x: float(x.quantile(.75))),
            reference_density_median=("reference_population_density", "median"),
            fraction_below_100=("reference_population_density", lambda x: float(np.mean(x < 100))),
        ).reset_index()
    )
    rows = (
        endpoints.groupby(["CNTR_CODE", "model"], observed=True)
        .agg(
            observed_change_log10=("observed_log10_change", "mean"),
            predicted_change_log10=("predicted_log10_change", "mean"),
            median_observed_change_percent=("observed_change_percent", "median"),
            observed_change_nw_mean=("observed_change_nw", "mean"),
            predicted_change_nw_mean=("predicted_change_nw", "mean"),
        ).reset_index()
    )
    rows["observed_slope_log10_per_year"] = rows["observed_change_log10"] / 10.0
    rows["predicted_slope_log10_per_year"] = rows["predicted_change_log10"] / 10.0
    rows["country_residual_slope_log10"] = (
        rows["observed_slope_log10_per_year"] -
        rows["predicted_slope_log10_per_year"]
    )
    rows = rows.merge(country_meta, on="CNTR_CODE", validate="many_to_one")
    country = rows.pivot(
        index="CNTR_CODE", columns="model",
        values="predicted_slope_log10_per_year",
    )
    country.columns = ["pred_slope_" + c for c in country.columns]
    country = country.reset_index()
    obs = (
        rows.loc[rows["model"].eq("D"), [
            "CNTR_CODE", "observed_slope_log10_per_year",
            "median_observed_change_percent",
        ]]
    )
    country = country.merge(
        country_meta, on="CNTR_CODE", validate="one_to_one"
    ).merge(obs, on="CNTR_CODE", validate="one_to_one")
    country["density_unexplained_slope"] = (
        country["observed_slope_log10_per_year"] - country["pred_slope_D"]
    )
    for name in ["D_H", "D_HxD", "D_country",
                 "D_H_country", "D_HxD_country"]:
        country["residual_after_" + name] = (
            country["observed_slope_log10_per_year"] -
            country["pred_slope_" + name]
        )
    country["H_net_added_slope"] = (
        country["pred_slope_D_HxD"] - country["pred_slope_D"]
    )
    country["unexplained_after_H"] = country["residual_after_D_HxD"]

    # Compare national errors under identical region compositions.
    metrics = []
    for name in MODELS:
        g = rows.loc[rows["model"].eq(name)].copy()
        e = g["country_residual_slope_log10"].to_numpy(float)
        weights = g["n_nuts3"].to_numpy(float)
        metrics.append({
            "model": name, "n_countries": len(g),
            "country_rmse_equal_country": float(np.sqrt(np.mean(e ** 2))),
            "country_rmse_weighted_nuts3": float(np.sqrt(np.average(e ** 2, weights=weights))),
            "country_mae_equal_country": float(np.mean(np.abs(e))),
            "country_mae_weighted_nuts3": float(np.average(np.abs(e), weights=weights)),
        })
    m = pd.DataFrame(metrics)
    d = m.set_index("model").loc["D"]
    for col in ("country_rmse_equal_country", "country_rmse_weighted_nuts3"):
        m[col + "_reduction_vs_D_pct"] = (
            100.0 * (1 - m[col] / float(d[col]))
        )
    return rows, country, m


def holdout_error_dominance(pred: pd.DataFrame) -> pd.DataFrame:
    """Expose how a few bright regions can dominate squared nW errors."""
    rows = []
    for model, g in pred.groupby("model", observed=True):
        # A region is represented in each of the 5 test targets: group first
        # to avoid labelling the largest *individual prediction* a region.
        region = (
            g.groupby("NUTS_ID")["residual_nw"]
            .apply(lambda s: float(np.sum(np.square(s))))
            .sort_values(ascending=False)
        )
        total = float(region.sum())
        rows.append({
            "model": model,
            "n_nuts3": len(region),
            "share_nw_squared_error_top_1_region": (
                float(region.head(1).sum() / total) if total > 0 else np.nan
            ),
            "share_nw_squared_error_top_5_regions": (
                float(region.head(5).sum() / total) if total > 0 else np.nan
            ),
            "share_nw_squared_error_top_10_regions": (
                float(region.head(10).sum() / total) if total > 0 else np.nan
            ),
            "worst_nuts3": str(region.index[0]),
        })
    return pd.DataFrame(rows)


def plot_holdout(folds: pd.DataFrame, out: Path, dpi: int) -> None:
    names = ["persistence", *MODELS.keys()]
    targets = [(2018, 2019), (2018, 2020), (2020, 2022), (2020, 2023), (2020, 2024)]
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax, col, ylabel in [
        (axes[0], "rmse_percentage_points", "RMSE in percentage points"),
        (axes[1], "rmse_nw", "RMSE in nW cm$^{-2}$ sr$^{-1}$"),
    ]:
        for model in names:
            t = folds.loc[folds["model"].eq(model)]
            values = []
            for cutoff, year in targets:
                q = t.loc[t["train_through"].eq(cutoff) & t["target_year"].eq(year)]
                values.append(float(q[col].iloc[0]))
            ax.plot(range(len(targets)), values, marker="o", linewidth=1.3, label=model)
        ax.set_xticks(range(len(targets)), [
            f"{s}→{e}" for s, e in targets
        ], rotation=20)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8, ncol=2)
    fig.suptitle("Temporal holdouts: density/H versus country effects and persistence")
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


def plot_composition(country: pd.DataFrame, out: Path, dpi: int) -> None:
    # Zero on the diagonal would mean the H shift perfectly explains
    # the national residual left by the density-only model.
    x = country["density_unexplained_slope"].to_numpy(float)
    y = country["H_net_added_slope"].to_numpy(float)
    weights = country["n_nuts3"].to_numpy(float)
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(x * 1000, y * 1000, s=22 + 2.2 * np.sqrt(weights), alpha=.75)
    for _, row in country.iterrows():
        ax.annotate(
            row["CNTR_CODE"],
            (row["density_unexplained_slope"] * 1000,
             row["H_net_added_slope"] * 1000),
            xytext=(3, 2), textcoords="offset points", fontsize=7.5,
        )
    limit = max(float(np.max(np.abs(np.r_[x, y])) * 1000), .1)
    ax.plot([-limit, limit], [-limit, limit], "--", color=".4", linewidth=.8)
    ax.axhline(0, color=".3", linewidth=.7)
    ax.axvline(0, color=".3", linewidth=.7)
    ax.set_xlim(-limit * 1.1, limit * 1.1)
    ax.set_ylim(-limit * 1.1, limit * 1.1)
    ax.set_xlabel("Observed national slope not explained by density alone\n(log10 radiance/year × 1000)")
    ax.set_ylabel("Net national slope adjustment when H is added\n(log10 radiance/year × 1000)")
    ax.set_title("Can national differences in H composition explain country trends?\nDiagonal: complete explanation; not a causal test")
    ax.grid(alpha=.15)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    cfg = config["analysis"]["presnow_density_h"]
    base, panel = load_panel(config, cfg)
    expected = int(config["analysis"]["fixed_nuts2016_2021_diagnostic"]["expected_n_nuts3"])
    _assert_locked_panel(panel, expected)
    print(f"Locked NUTS2016 cohort: {expected} regions, {panel['CNTR_CODE'].nunique()} countries", flush=True)

    full_models, fit_table = fit_models(panel, MODELS)
    nested = nested_evidence(full_models)
    print("Fitted full-period models; running two temporal holdout experiments.", flush=True)
    training_fit, holdouts, overall, by_fold, by_density, by_country = (
        test_temporal_predictions(panel)
    )
    endpoints = endpoint_full_fit(panel, full_models)
    country_model_rows, country_composition_table, country_errors = (
        country_composition(panel, endpoints)
    )
    error_dominance = holdout_error_dominance(holdouts)

    root = configured_path(config, "processed") / STAGE
    figures = configured_path(config, "figures") / "population_relationship" / STAGE
    root.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    dpi = int(cfg["figure_dpi"])
    outputs = {
        "full_period_model_fit": fit_table,
        "full_period_nested_tests": nested,
        "training_period_model_fit": training_fit,
        "temporal_holdout_predictions": holdouts,
        "temporal_holdout_overall": overall,
        "temporal_holdout_by_fold": by_fold,
        "temporal_holdout_by_density": by_density,
        "temporal_holdout_by_country": by_country,
        "in_sample_endpoint_predictions": endpoints,
        "country_model_residual_slopes": country_model_rows,
        "country_h_composition_explanation": country_composition_table,
        "country_residual_error_comparison": country_errors,
        "absolute_error_concentration": error_dominance,
    }
    for stem, frame in outputs.items():
        frame.to_csv(root / f"{stem}.csv", index=False)
    plot_holdout(
        by_fold, figures / "temporal_holdout_model_comparison.pdf", dpi
    )
    plot_composition(
        country_composition_table,
        figures / "country_h_composition_vs_unexplained_trends.pdf", dpi,
    )
    summary = root / "country_h_temporal_validation_summary.txt"
    with summary.open("w", encoding="utf-8") as f:
        f.write("CAN H COMPOSITION EXPLAIN COUNTRY-SPECIFIC VIIRS TRENDS?\n")
        f.write("=" * 77 + "\n")
        f.write(
            "Original fixed NUTS2016 2014-2024 pixel-calibrated history; "
            "no special 2021 snow correction. 1514 regions, unchanged "
            "2014-2020 reference density classes, independent H_core.\n"
            "Hypothesis: national variation in H and density composition "
            "can explain some national differences in radiance change.\n"
            "Country-specific time slopes represent residual national "
            "circumstances, not biological drivers.\n"
            "No automatic model winner is declared.\n"
            "2021 is part of FULL-PERIOD fits but not any holdout "
            "training period or tested target.\n"
            "Temporal holdouts: train 2014-2018, predict 2019-2020 anchored "
            "on observed 2018; train 2014-2020, predict 2022-2024 anchored "
            "on observed 2020. The anchor itself is observed, not forecast.\n"
            "Country composition analysis uses the IN-SAMPLE 2014-2024 "
            "endpoint fit; it evaluates explanatory plausibility, not "
            "future predictive skill or causal mediation.\n"
            "H uses country/region static values, no density/H substitution.\n"
            "Results include proportional and physical-unit errors and a "
            "no-change (persistence) forecast.\n"
            "We deliberately do NOT introduce a changing population-density "
            "predictor without harmonised annual NUTS2016 data after 2020.\n\n"
        )
        for label, frame in [
            ("FULL MODEL FIT", fit_table),
            ("NESTED TESTS", nested),
            ("HOLDOUT OVERALL", overall),
            ("HOLDOUT BY YEAR", by_fold),
            ("COUNTRY-SLOPE ERROR", country_errors),
            ("ABSOLUTE ERROR CONCENTRATION", error_dominance),
        ]:
            f.write(label + "\n" + frame.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            ) + "\n\n")
    print("\nPREDICTIVE HOLDOUT\n" + overall.to_string(index=False))
    print("\nCOUNTRY EXPLANATION\n" + country_errors.to_string(index=False))
    print("\nWROTE", summary, flush=True)
    for stem in outputs:
        print(root / f"{stem}.csv")
    print(figures / "temporal_holdout_model_comparison.pdf")
    print(figures / "country_h_composition_vs_unexplained_trends.pdf")


if __name__ == "__main__":
    main()
