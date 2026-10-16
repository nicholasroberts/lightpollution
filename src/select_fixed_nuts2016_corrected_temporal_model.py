#!/usr/bin/env python3
"""Select the minimum temporal model for the corrected 2014-2024 VIIRS panel.

This is the model-selection stage that follows the 2021 surface correction.
It uses the multiplicative/log-ratio corrected 2021 value and does not revisit
the snow/albedo mechanism.

Candidate ladder
----------------
M0 common_linear
    one common linear time slope across all fixed density classes.

M1 density_linear
    linear time slope differs by fixed population-density class.

M2 density_linear_common_hinge
    M1 plus one common post-2020 change in slope.

M3 density_specific_hinge
    M2 plus density-specific post-2020 changes in slope.

All models retain country fixed effects and a NUTS3 random intercept.

Selection is strictly sequential. A more complex model is retained only if it:
1. significantly improves the nested likelihood-ratio test;
2. improves BIC by the configured minimum; and
3. improves blocked one-year-ahead temporal holdout RMSE by the configured
   minimum.

The holdout folds train through 2022 and 2023, then predict the within-region
change to 2023 and 2024 respectively. Because the prediction target is change
within the same NUTS3 region, intercept terms cancel and the test directly
assesses the transportability of the temporal terms.

Outputs include observed versus selected-model-predicted 2014->2024 absolute
and percentage changes for each fixed population-density class.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_density_temporal_change import fit_mixedlm, lrt
from project_config import configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import prepare_analysis_panel


MODEL_ORDER = (
    "M0_common_linear",
    "M1_density_linear",
    "M2_density_linear_common_hinge",
    "M3_density_specific_hinge",
)


def formulas() -> dict[str, str]:
    return {
        "M0_common_linear": (
            "log_rad ~ year_c + C(density_band) + C(country)"
        ),
        "M1_density_linear": (
            "log_rad ~ year_c * C(density_band) + C(country)"
        ),
        "M2_density_linear_common_hinge": (
            "log_rad ~ year_c * C(density_band) + hinge + C(country)"
        ),
        "M3_density_specific_hinge": (
            "log_rad ~ year_c * C(density_band) + "
            "hinge * C(density_band) + C(country)"
        ),
    }


def restore_density_categories(
    reference: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    labels = list(
        config["analysis"]["density_temporal_change"][
            "density_band_labels"
        ]
    )
    out = reference.copy()
    out["density_band"] = pd.Categorical(
        out["density_band"].astype(str),
        categories=labels,
        ordered=True,
    )
    if out["density_band"].isna().any():
        raise RuntimeError(
            "Density reference contains labels outside the locked classes."
        )
    return out


def load_corrected_panel(
    config: dict,
    cfg: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    processed = configured_path(config, "processed")
    history_path = (
        processed
        / str(cfg["corrected_history_file"])
    )
    reference_path = (
        processed
        / "fixed_nuts2016_population_density_reference.csv"
    )
    if not history_path.exists():
        raise FileNotFoundError(
            f"{history_path}. Run "
            "'python src/run_pipeline.py future-radiance-2021-surface-correction' first."
        )
    if not reference_path.exists():
        raise FileNotFoundError(reference_path)

    history = pd.read_parquet(history_path)
    reference = restore_density_categories(
        pd.read_csv(reference_path), config
    )
    for d in (history, reference):
        d["NUTS_ID"] = d["NUTS_ID"].astype(str)

    source_col = str(cfg["corrected_radiance_column"])
    if source_col not in history.columns:
        raise KeyError(
            f"Corrected history lacks {source_col!r}. "
            "Rerun the 2021 surface-correction stage after pulling the latest branch."
        )

    history = history.copy()
    history["radiance_mean_corrected"] = pd.to_numeric(
        history[source_col], errors="coerce"
    )
    panel = prepare_analysis_panel(
        history,
        reference,
        config,
    )

    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    panel = panel[
        panel["year"].between(first_year, last_year)
    ].copy()

    expected_years = last_year - first_year + 1
    counts = panel.groupby("NUTS_ID")["year"].nunique()
    keep = counts[counts == expected_years].index
    panel = panel[panel["NUTS_ID"].isin(keep)].copy()

    labels = list(
        config["analysis"]["density_temporal_change"][
            "density_band_labels"
        ]
    )
    panel["density_band"] = pd.Categorical(
        panel["density_band"].astype(str),
        categories=labels,
        ordered=True,
    )
    panel["country"] = panel["CNTR_CODE"].astype("category")
    return panel, reference


def add_time_terms(
    d: pd.DataFrame,
    first_year: int,
    last_year: int,
    hinge_year: int,
) -> pd.DataFrame:
    out = d.copy()
    midpoint = (first_year + last_year) / 2.0
    out["year_c"] = out["year"].astype(float) - midpoint
    out["hinge"] = np.maximum(
        0.0, out["year"].astype(float) - float(hinge_year)
    )
    out.attrs["year_midpoint"] = midpoint
    out.attrs["hinge_year"] = hinge_year
    return out


def fit_candidates(d: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    models = {}
    rows = []
    for name in MODEL_ORDER:
        result, optimizer = fit_mixedlm(formulas()[name], d)
        models[name] = result
        rows.append(
            {
                "model": name,
                "formula": formulas()[name],
                "n_obs": int(result.nobs),
                "n_nuts3": int(d["NUTS_ID"].nunique()),
                "n_countries": int(d["CNTR_CODE"].nunique()),
                "log_likelihood": float(result.llf),
                "aic": float(result.aic),
                "bic": float(result.bic),
                "optimizer": optimizer,
                "random_intercept_variance": float(
                    np.asarray(result.cov_re)[0, 0]
                ),
            }
        )
    return models, pd.DataFrame(rows)


def nested_tests(models: dict) -> pd.DataFrame:
    rows = []
    for reduced_name, full_name in zip(MODEL_ORDER[:-1], MODEL_ORDER[1:]):
        lr, df, p = lrt(
            models[reduced_name],
            models[full_name],
        )
        rows.append(
            {
                "reduced_model": reduced_name,
                "full_model": full_name,
                "lr": lr,
                "df": df,
                "p_value": p,
                "bic_reduced": float(models[reduced_name].bic),
                "bic_full": float(models[full_name].bic),
                "bic_improvement": float(
                    models[reduced_name].bic - models[full_name].bic
                ),
            }
        )
    return pd.DataFrame(rows)


def find_param(result, candidates: list[str]) -> str | None:
    names = set(result.fe_params.index)
    return next((x for x in candidates if x in names), None)


def change_contrast_weights(
    result,
    model_name: str,
    band: str,
    ref_band: str,
    year0: int,
    year1: int,
    hinge_year: int,
) -> dict[str, float]:
    dt = float(year1 - year0)
    dh = (
        max(0.0, float(year1 - hinge_year))
        - max(0.0, float(year0 - hinge_year))
    )
    w: dict[str, float] = {}

    if "year_c" in result.fe_params.index:
        w["year_c"] = dt

    if model_name != "M0_common_linear" and band != ref_band:
        name = find_param(
            result,
            [
                f"year_c:C(density_band)[T.{band}]",
                f"C(density_band)[T.{band}]:year_c",
            ],
        )
        if name is None:
            raise RuntimeError(
                f"Could not find year-by-density term for {band} in {model_name}."
            )
        w[name] = dt

    if model_name in (
        "M2_density_linear_common_hinge",
        "M3_density_specific_hinge",
    ):
        if "hinge" not in result.fe_params.index:
            raise RuntimeError(f"Missing common hinge in {model_name}.")
        w["hinge"] = dh

    if (
        model_name == "M3_density_specific_hinge"
        and band != ref_band
    ):
        name = find_param(
            result,
            [
                f"hinge:C(density_band)[T.{band}]",
                f"C(density_band)[T.{band}]:hinge",
            ],
        )
        if name is None:
            raise RuntimeError(
                f"Could not find hinge-by-density term for {band}."
            )
        w[name] = dh
    return w


def fixed_contrast(
    result,
    weights: dict[str, float],
) -> tuple[float, float]:
    names = list(result.fe_params.index)
    L = np.zeros(len(names), dtype=float)
    for name, value in weights.items():
        L[names.index(name)] = value
    beta = result.fe_params.to_numpy(dtype=float)
    cov = result.cov_params().loc[names, names].to_numpy(dtype=float)
    est = float(L @ beta)
    se = float(np.sqrt(max(0.0, L @ cov @ L)))
    return est, se


def temporal_holdout(
    panel: pd.DataFrame,
    config: dict,
    cfg: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    hinge_year = int(cfg["hinge_year"])
    cutoffs = [int(x) for x in cfg["holdout_cutoffs"]]
    labels = list(
        config["analysis"]["density_temporal_change"][
            "density_band_labels"
        ]
    )
    ref_band = labels[0]

    fold_rows = []
    prediction_rows = []

    for cutoff in cutoffs:
        test_year = cutoff + 1
        if test_year > last_year:
            raise ValueError(
                f"Holdout cutoff {cutoff} has no test year <= {last_year}."
            )

        train = panel[panel["year"].between(first_year, cutoff)].copy()
        train = add_time_terms(
            train, first_year, last_year, hinge_year
        )
        fitted, _ = fit_candidates(train)

        a = panel[panel["year"].eq(cutoff)][
            ["NUTS_ID", "CNTR_CODE", "density_band", "log_rad"]
        ].rename(columns={"log_rad": "log_rad_start"})
        b = panel[panel["year"].eq(test_year)][
            ["NUTS_ID", "CNTR_CODE", "density_band", "log_rad"]
        ].rename(columns={"log_rad": "log_rad_end"})
        test = a.merge(
            b,
            on=["NUTS_ID", "CNTR_CODE", "density_band"],
            how="inner",
            validate="one_to_one",
        )
        test["observed_log10_change"] = (
            test["log_rad_end"] - test["log_rad_start"]
        )
        test["observed_pct_change"] = (
            100.0 * (10.0 ** test["observed_log10_change"] - 1.0)
        )

        for model_name in MODEL_ORDER:
            model = fitted[model_name]
            pred_by_band = {}
            for band in labels:
                w = change_contrast_weights(
                    model,
                    model_name,
                    band,
                    ref_band,
                    cutoff,
                    test_year,
                    hinge_year,
                )
                est, _ = fixed_contrast(model, w)
                pred_by_band[band] = est

            x = test.copy()
            x["model"] = model_name
            x["train_through"] = cutoff
            x["test_year"] = test_year
            x["predicted_log10_change"] = x[
                "density_band"
            ].astype(str).map(pred_by_band)
            x["predicted_pct_change"] = (
                100.0 * (
                    10.0 ** x["predicted_log10_change"] - 1.0
                )
            )
            x["residual_log10"] = (
                x["predicted_log10_change"]
                - x["observed_log10_change"]
            )
            x["residual_pct_points"] = (
                x["predicted_pct_change"]
                - x["observed_pct_change"]
            )
            prediction_rows.append(x)

            fold_rows.append(
                {
                    "model": model_name,
                    "train_through": cutoff,
                    "test_year": test_year,
                    "n_test": len(x),
                    "rmse_log10": float(
                        np.sqrt(
                            np.mean(
                                x["residual_log10"].to_numpy(float) ** 2
                            )
                        )
                    ),
                    "mae_log10": float(
                        np.mean(np.abs(x["residual_log10"]))
                    ),
                    "bias_log10": float(
                        x["residual_log10"].mean()
                    ),
                    "rmse_pct_points": float(
                        np.sqrt(
                            np.mean(
                                x["residual_pct_points"].to_numpy(float)
                                ** 2
                            )
                        )
                    ),
                    "mae_pct_points": float(
                        np.mean(np.abs(x["residual_pct_points"]))
                    ),
                    "bias_pct_points": float(
                        x["residual_pct_points"].mean()
                    ),
                }
            )

    folds = pd.DataFrame(fold_rows)
    preds = pd.concat(prediction_rows, ignore_index=True)

    agg = (
        preds.groupby("model", observed=True)
        .agg(
            n_predictions=("NUTS_ID", "size"),
            rmse_log10=(
                "residual_log10",
                lambda s: float(np.sqrt(np.mean(np.square(s)))),
            ),
            mae_log10=("residual_log10", lambda s: float(np.mean(np.abs(s)))),
            bias_log10=("residual_log10", "mean"),
            rmse_pct_points=(
                "residual_pct_points",
                lambda s: float(np.sqrt(np.mean(np.square(s)))),
            ),
            mae_pct_points=(
                "residual_pct_points",
                lambda s: float(np.mean(np.abs(s)))),
            bias_pct_points=("residual_pct_points", "mean"),
        )
        .reset_index()
    )

    base = agg.set_index("model")["rmse_log10"].to_dict()
    improvements = []
    for i, name in enumerate(MODEL_ORDER):
        prev = MODEL_ORDER[i - 1] if i > 0 else None
        improvements.append(
            np.nan
            if prev is None
            else 100.0 * (
                base[prev] - base[name]
            ) / base[prev]
        )
    agg["rmse_improvement_vs_previous_pct"] = improvements
    return agg, folds


def choose_minimum_model(
    fit_table: pd.DataFrame,
    tests: pd.DataFrame,
    holdout: pd.DataFrame,
    cfg: dict,
) -> tuple[str, pd.DataFrame]:
    alpha = float(cfg["selection_alpha"])
    min_bic = float(cfg["min_bic_improvement"])
    min_rmse = float(cfg["min_holdout_rmse_improvement_pct"])

    fit = fit_table.set_index("model")
    ho = holdout.set_index("model")
    decisions = []
    selected = MODEL_ORDER[0]

    for reduced, full in zip(MODEL_ORDER[:-1], MODEL_ORDER[1:]):
        t = tests[
            (tests["reduced_model"] == reduced)
            & (tests["full_model"] == full)
        ].iloc[0]
        p = float(t["p_value"])
        bic_gain = float(t["bic_improvement"])
        rmse_gain = float(
            ho.loc[full, "rmse_improvement_vs_previous_pct"]
        )
        passes = (
            p < alpha
            and bic_gain >= min_bic
            and rmse_gain >= min_rmse
            and selected == reduced
        )
        decisions.append(
            {
                "reduced_model": reduced,
                "full_model": full,
                "lrt_p": p,
                "bic_improvement": bic_gain,
                "holdout_rmse_improvement_pct": rmse_gain,
                "passes_all_rules": passes,
            }
        )
        if passes:
            selected = full
        else:
            break

    return selected, pd.DataFrame(decisions)


def observed_vs_predicted_change(
    panel: pd.DataFrame,
    model,
    model_name: str,
    config: dict,
    cfg: dict,
) -> pd.DataFrame:
    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    hinge_year = int(cfg["hinge_year"])
    labels = list(
        config["analysis"]["density_temporal_change"][
            "density_band_labels"
        ]
    )
    ref_band = labels[0]

    a = panel[panel["year"].eq(first_year)][
        ["NUTS_ID", "density_band", "radiance_mean_corrected"]
    ].rename(
        columns={"radiance_mean_corrected": "radiance_start"}
    )
    b = panel[panel["year"].eq(last_year)][
        ["NUTS_ID", "density_band", "radiance_mean_corrected"]
    ].rename(
        columns={"radiance_mean_corrected": "radiance_end"}
    )
    paired = a.merge(
        b,
        on=["NUTS_ID", "density_band"],
        how="inner",
        validate="one_to_one",
    )

    rows = []
    years = last_year - first_year
    for band in labels:
        g = paired[paired["density_band"].astype(str).eq(str(band))].copy()
        if g.empty:
            continue
        log_start = np.log10(g["radiance_start"].to_numpy(float))
        log_end = np.log10(g["radiance_end"].to_numpy(float))
        geomean_start = float(10.0 ** log_start.mean())
        geomean_end = float(10.0 ** log_end.mean())
        observed_log_change = float((log_end - log_start).mean())
        observed_pct = 100.0 * (
            10.0 ** observed_log_change - 1.0
        )
        observed_abs = geomean_end - geomean_start

        w = change_contrast_weights(
            model,
            model_name,
            band,
            ref_band,
            first_year,
            last_year,
            hinge_year,
        )
        pred_log, pred_se = fixed_contrast(model, w)
        pred_ratio = 10.0 ** pred_log
        pred_pct = 100.0 * (pred_ratio - 1.0)
        pred_end = geomean_start * pred_ratio
        pred_abs = pred_end - geomean_start

        lo_log = pred_log - 1.96 * pred_se
        hi_log = pred_log + 1.96 * pred_se

        rows.append(
            {
                "density_band": band,
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "geometric_mean_radiance_2014_nw": geomean_start,
                "geometric_mean_radiance_2024_nw": geomean_end,
                "observed_absolute_change_nw": observed_abs,
                "observed_total_pct_change": observed_pct,
                "observed_annualized_pct_change": (
                    100.0 * (
                        10.0 ** (observed_log_change / years) - 1.0
                    )
                ),
                "predicted_radiance_2024_nw": pred_end,
                "predicted_absolute_change_nw": pred_abs,
                "predicted_total_pct_change": pred_pct,
                "predicted_total_pct_ci_low": (
                    100.0 * (10.0 ** lo_log - 1.0)
                ),
                "predicted_total_pct_ci_high": (
                    100.0 * (10.0 ** hi_log - 1.0)
                ),
                "predicted_annualized_pct_change": (
                    100.0 * (
                        10.0 ** (pred_log / years) - 1.0
                    )
                ),
                "absolute_change_error_nw": pred_abs - observed_abs,
                "total_pct_change_error_pp": pred_pct - observed_pct,
                "selected_model": model_name,
            }
        )
    return pd.DataFrame(rows)


def fitted_trajectory_table(
    panel: pd.DataFrame,
    model,
    model_name: str,
    config: dict,
    cfg: dict,
) -> pd.DataFrame:
    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    hinge_year = int(cfg["hinge_year"])
    labels = list(
        config["analysis"]["density_temporal_change"][
            "density_band_labels"
        ]
    )
    ref_band = labels[0]

    baseline = panel[panel["year"].eq(first_year)].copy()
    base_geo = (
        baseline.groupby("density_band", observed=True)[
            "radiance_mean_corrected"
        ]
        .agg(
            lambda s: float(
                10.0 ** np.log10(s.to_numpy(float)).mean()
            )
        )
        .rename("geometric_mean_radiance_2014_nw")
    )

    observed = (
        panel.groupby(["density_band", "year"], observed=True)[
            "radiance_mean_corrected"
        ]
        .agg(
            lambda s: float(
                10.0 ** np.log10(s.to_numpy(float)).mean()
            )
        )
        .rename("observed_geometric_mean_nw")
        .reset_index()
    )

    rows = []
    for band in labels:
        start = float(base_geo.loc[band])
        for year in range(first_year, last_year + 1):
            if year == first_year:
                pred_log = 0.0
            else:
                w = change_contrast_weights(
                    model,
                    model_name,
                    band,
                    ref_band,
                    first_year,
                    year,
                    hinge_year,
                )
                pred_log, _ = fixed_contrast(model, w)
            rows.append(
                {
                    "density_band": band,
                    "year": year,
                    "predicted_geometric_mean_nw": (
                        start * 10.0 ** pred_log
                    ),
                    "predicted_index_2014_100": (
                        100.0 * 10.0 ** pred_log
                    ),
                }
            )
    pred = pd.DataFrame(rows)
    out = observed.merge(
        pred,
        on=["density_band", "year"],
        how="left",
        validate="one_to_one",
    )
    out["observed_index_2014_100"] = (
        100.0
        * out["observed_geometric_mean_nw"]
        / out.groupby("density_band", observed=True)[
            "observed_geometric_mean_nw"
        ].transform("first")
    )
    return out


def plot_selected_trajectories(
    trajectories: pd.DataFrame,
    config: dict,
    output: Path,
) -> None:
    labels = list(
        config["analysis"]["density_temporal_change"][
            "density_band_labels"
        ]
    )
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0), sharex=True)
    axes = axes.ravel()

    for i, band in enumerate(labels):
        ax = axes[i]
        g = trajectories[
            trajectories["density_band"].astype(str).eq(str(band))
        ].sort_values("year")
        ax.plot(
            g["year"],
            g["observed_index_2014_100"],
            marker="o",
            label="Observed corrected",
        )
        ax.plot(
            g["year"],
            g["predicted_index_2014_100"],
            linestyle="--",
            linewidth=1.8,
            label="Selected minimum model",
        )
        ax.axhline(100, linewidth=0.7)
        ax.set_title(str(band))
        ax.grid(alpha=0.2)

    axes[5].axis("off")
    handles, labs = axes[0].get_legend_handles_labels()
    axes[5].legend(handles, labs, loc="center", frameon=False)
    fig.supylabel("Corrected VIIRS radiance index (2014 = 100)")
    fig.supxlabel("Year")
    fig.suptitle(
        "Corrected 2014-2024 history versus selected minimum temporal model"
    )
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.95))
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_corrected_temporal_model"]
    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    hinge_year = int(cfg["hinge_year"])

    panel, _ = load_corrected_panel(config, cfg)
    model_data = add_time_terms(
        panel, first_year, last_year, hinge_year
    )

    models, fit_table = fit_candidates(model_data)
    tests = nested_tests(models)
    holdout, holdout_folds = temporal_holdout(
        panel, config, cfg
    )
    selected, decisions = choose_minimum_model(
        fit_table, tests, holdout, cfg
    )

    change = observed_vs_predicted_change(
        panel,
        models[selected],
        selected,
        config,
        cfg,
    )
    trajectories = fitted_trajectory_table(
        panel,
        models[selected],
        selected,
        config,
        cfg,
    )

    processed = configured_path(config, "processed")
    out_dir = processed / "fixed_nuts2016_corrected_temporal_model"
    fig_dir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_corrected_temporal_model"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "fit": out_dir / "candidate_model_fit.csv",
        "tests": out_dir / "nested_model_tests.csv",
        "holdout": out_dir / "temporal_holdout_summary.csv",
        "holdout_folds": out_dir / "temporal_holdout_folds.csv",
        "decisions": out_dir / "minimum_model_selection.csv",
        "change": out_dir / "observed_vs_predicted_2014_2024_by_density.csv",
        "trajectories": out_dir / "selected_model_trajectories.csv",
        "summary": out_dir / "corrected_minimum_model_summary.txt",
    }

    fit_table.to_csv(outputs["fit"], index=False)
    tests.to_csv(outputs["tests"], index=False)
    holdout.to_csv(outputs["holdout"], index=False)
    holdout_folds.to_csv(outputs["holdout_folds"], index=False)
    decisions.to_csv(outputs["decisions"], index=False)
    change.to_csv(outputs["change"], index=False)
    trajectories.to_csv(outputs["trajectories"], index=False)

    fig = fig_dir / "selected_model_trajectories_2014_2024.pdf"
    plot_selected_trajectories(trajectories, config, fig)

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("MINIMUM TEMPORAL MODEL FOR CORRECTED 2014-2024 VIIRS\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            "Input: multiplicative/log-ratio 2021 surface-corrected fixed-NUTS2016 "
            "annual radiance. The production 2017+ calibration remains embedded "
            "in this series.\n"
        )
        h.write(
            f"Years: {first_year}-{last_year}; fixed hinge for candidate M2/M3: "
            f"{hinge_year}.\n"
        )
        h.write(
            "Model ladder: M0 common linear -> M1 density-specific linear -> "
            "M2 common post-2020 slope change -> M3 density-specific post-2020 "
            "slope changes.\n"
        )
        h.write(
            "Selection rule: advance one step only when LRT p < "
            f"{float(cfg['selection_alpha']):g}, BIC improves by >= "
            f"{float(cfg['min_bic_improvement']):g}, and blocked temporal "
            "holdout RMSE improves by >= "
            f"{float(cfg['min_holdout_rmse_improvement_pct']):g}%.\n\n"
        )

        h.write("CANDIDATE MODEL FIT\n")
        h.write(
            fit_table.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\nNESTED TESTS\n")
        h.write(
            tests.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\nTEMPORAL HOLDOUT\n")
        h.write(
            holdout.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\nSEQUENTIAL DECISIONS\n")
        h.write(
            decisions.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write(f"\n\nSELECTED MODEL: {selected}\n")
        h.write(
            "\nOBSERVED VS SELECTED-MODEL 2014-2024 CHANGE BY DENSITY\n"
        )
        h.write(
            change.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )

    print("\nCORRECTED MINIMUM TEMPORAL MODEL")
    print(f"Selected: {selected}")
    print("\nMODEL SELECTION")
    print(
        decisions.to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print("\nOBSERVED VS PREDICTED 2014-2024")
    print(
        change.to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print("\nWROTE")
    for p in [*outputs.values(), fig]:
        print(p)


if __name__ == "__main__":
    main()
