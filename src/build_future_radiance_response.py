#!/usr/bin/env python3
"""Freeze the historical population -> calibrated-VIIRS response for forecasting.

Purpose
-------
This stage is the end of the historical radiance-model selection step and the
start of the future projection pipeline. It deliberately compares a very small
set of models and chooses the simplest model that is both statistically
supported and predictively better.

The forecast response is the annual change in calibrated VIIRS DNB radiance,
g_R, expressed internally as annual change in log10 radiance and reported as
percent per year.

Candidate models
----------------
All models include country fixed effects and a cubic function of baseline
log10 population density, D.

M0 density_only
    g_R = f(D)

M1 density_plus_growth
    g_R = f(D) + g_D

M2 density_growth_interaction
    g_R = f(D) + g_D + D:g_D

where g_D is annual percentage change in population density.

Model selection is intentionally conservative. A demographic term is retained
only if:
1. its clustered Wald block is significant in the historical 5-year-pair fit;
2. it improves genuinely non-overlapping temporal holdout prediction by at
   least the configured RMSE threshold; and
3. the more complex interaction model must additionally improve on M1.

The 5-year historical fit uses 2013-2018, 2014-2019 and 2015-2020 on fixed
NUTS2016 geography. The independent temporal holdout uses a pre-calibration
training interval (2013-2016) and a post-calibration test interval (2017-2020).
No year belongs to both intervals.

The production calibration is validated before analysis. This stage refuses to
run if the merged table does not carry the configured Europe-wide pixel-level
calibration provenance.

The headline five population-density-class trends remain those from
density-temporal-change. This stage selects only the transfer function used for
future prediction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2
import statsmodels.formula.api as smf

from project_config import configured_path, load_config


MODEL_FORMULAS = {
    "density_only": (
        "annual_log10_radiance_change ~ C(country) + C(window) + "
        "log_density_c + I(log_density_c ** 2) + I(log_density_c ** 3)"
    ),
    "density_plus_growth": (
        "annual_log10_radiance_change ~ C(country) + C(window) + "
        "log_density_c + I(log_density_c ** 2) + I(log_density_c ** 3) + "
        "population_density_change_pct_per_year"
    ),
    "density_growth_interaction": (
        "annual_log10_radiance_change ~ C(country) + C(window) + "
        "log_density_c + I(log_density_c ** 2) + I(log_density_c ** 3) + "
        "population_density_change_pct_per_year + "
        "log_density_c:population_density_change_pct_per_year"
    ),
}


def validate_production_calibration(data: pd.DataFrame, config: dict) -> None:
    cal = config["calibration"]
    required = [
        "radiance_mean_raw",
        "radiance_mean_corrected",
        "calibration_method",
        "calibration_version",
        "calibration_correction_applied",
    ]
    missing = [c for c in required if c not in data.columns]
    if missing:
        raise RuntimeError(
            "Merged production table lacks calibration provenance columns: "
            + ", ".join(missing)
        )

    if cal["method"] != "global_pixel_additive_offset_clipped_zero":
        raise RuntimeError(
            "This stage requires the production Europe-wide pixel correction; "
            f"configured method is {cal['method']!r}."
        )

    start = int(cal["apply_from_year"])
    expected_version = str(cal["version"])

    post = data[data["year"] >= start].copy()
    if post.empty:
        raise RuntimeError("No post-calibration observations found.")

    versions = set(post["calibration_version"].dropna().astype(str).unique())
    if versions != {expected_version}:
        raise RuntimeError(
            "Merged table calibration version does not match config. "
            f"Expected {expected_version!r}; found {sorted(versions)}. "
            "Rerun calibration and merge."
        )

    methods = set(post["calibration_method"].dropna().astype(str).unique())
    if methods != {cal["method"]}:
        raise RuntimeError(
            "Merged table calibration method does not match config. "
            f"Expected {cal['method']!r}; found {sorted(methods)}."
        )

    applied = post["calibration_correction_applied"].fillna(False).astype(bool)
    if not applied.all():
        raise RuntimeError(
            "At least one post-2017 merged row is not flagged as production "
            "calibrated. Rerun calibration and merge."
        )

    if (post["radiance_mean_corrected"] > post["radiance_mean_raw"] + 1e-10).any():
        raise RuntimeError(
            "Production correction increased at least one post-2017 NUTS3 mean."
        )


def load_base_data(config: dict) -> pd.DataFrame:
    processed = configured_path(config, "processed")
    path = processed / config["outputs"]["merged"]
    data = pd.read_parquet(path).copy()
    validate_production_calibration(data, config)

    needed = [
        "year",
        "NUTS_ID",
        "CNTR_CODE",
        "nuts_release",
        "population_density",
        "radiance_mean_corrected",
        "radiance_mean_raw",
        "calibration_method",
        "calibration_version",
        "calibration_correction_applied",
    ]
    missing = [c for c in needed if c not in data.columns]
    if missing:
        raise KeyError(f"Merged table is missing required columns: {missing}")

    for col in ["population_density", "radiance_mean_corrected"]:
        data[col] = pd.to_numeric(data[col], errors="coerce")

    valid = (
        data["NUTS_ID"].notna()
        & data["CNTR_CODE"].notna()
        & np.isfinite(data["population_density"])
        & np.isfinite(data["radiance_mean_corrected"])
        & (data["population_density"] > 0)
        & (data["radiance_mean_corrected"] > 0)
    )
    return data.loc[valid].copy()


def make_interval_pairs(
    data: pd.DataFrame,
    starts: list[int],
    horizon: int,
    nuts_release: int,
    min_country_regions: int,
) -> pd.DataFrame:
    fixed = data[data["nuts_release"] == nuts_release].copy()

    dup = fixed.duplicated(["NUTS_ID", "year"], keep=False)
    if dup.any():
        raise RuntimeError(
            "Duplicate NUTS_ID-year rows found in fixed-geography data."
        )

    rows = []
    for start in starts:
        end = int(start) + int(horizon)
        a = fixed[fixed["year"] == int(start)][
            [
                "NUTS_ID",
                "CNTR_CODE",
                "population_density",
                "radiance_mean_corrected",
            ]
        ].rename(
            columns={
                "population_density": "population_density_start",
                "radiance_mean_corrected": "radiance_start",
            }
        )
        b = fixed[fixed["year"] == end][
            [
                "NUTS_ID",
                "CNTR_CODE",
                "population_density",
                "radiance_mean_corrected",
            ]
        ].rename(
            columns={
                "population_density": "population_density_end",
                "radiance_mean_corrected": "radiance_end",
            }
        )

        p = a.merge(
            b,
            on=["NUTS_ID", "CNTR_CODE"],
            how="inner",
            validate="one_to_one",
        )
        p["start_year"] = int(start)
        p["end_year"] = end
        p["window"] = f"{int(start)}-{end}"
        p["baseline_log10_population_density"] = np.log10(
            p["population_density_start"]
        )
        p["population_density_change_pct_per_year"] = 100.0 * (
            (p["population_density_end"] / p["population_density_start"])
            ** (1.0 / horizon)
            - 1.0
        )
        p["annual_log10_radiance_change"] = (
            np.log10(p["radiance_end"]) - np.log10(p["radiance_start"])
        ) / horizon
        p["radiance_change_pct_per_year"] = 100.0 * (
            10.0 ** p["annual_log10_radiance_change"] - 1.0
        )
        rows.append(p)

    if not rows:
        raise RuntimeError("No interval pairs were constructed.")

    out = pd.concat(rows, ignore_index=True)
    counts = out.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    keep = counts[counts >= min_country_regions].index
    return out[out["CNTR_CODE"].isin(keep)].copy()


def trim_and_center(
    pairs: pd.DataFrame,
    trim_quantiles: tuple[float, float],
) -> pd.DataFrame:
    qlo, qhi = trim_quantiles
    lo = float(
        pairs["population_density_change_pct_per_year"].quantile(qlo)
    )
    hi = float(
        pairs["population_density_change_pct_per_year"].quantile(qhi)
    )
    fit = pairs[
        pairs["population_density_change_pct_per_year"].between(
            lo, hi, inclusive="both"
        )
    ].copy()

    center = float(fit["baseline_log10_population_density"].median())
    fit["log_density_c"] = (
        fit["baseline_log10_population_density"] - center
    )
    fit["country"] = fit["CNTR_CODE"].astype("category")
    fit["window"] = fit["window"].astype("category")
    fit.attrs["log_density_center"] = center
    fit.attrs["trim_low"] = lo
    fit.attrs["trim_high"] = hi
    return fit


def cluster_fit(formula: str, data: pd.DataFrame):
    return smf.ols(formula, data=data).fit(
        cov_type="cluster",
        cov_kwds={"groups": data["NUTS_ID"]},
    )


def wald_block(result, terms: list[str]) -> dict:
    names = list(result.params.index)
    present = [t for t in terms if t in names]
    if not present:
        return {"wald_chi2": np.nan, "df": 0, "p": np.nan}

    idx = [names.index(t) for t in present]
    beta = result.params.iloc[idx].to_numpy(dtype=float)
    cov = result.cov_params().iloc[idx, idx].to_numpy(dtype=float)
    stat = float(beta.T @ np.linalg.pinv(cov) @ beta)
    return {
        "wald_chi2": stat,
        "df": len(idx),
        "p": float(chi2.sf(stat, len(idx))),
    }


def fit_candidate_models(fit: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    models = {}
    rows = []

    for name, formula in MODEL_FORMULAS.items():
        model = cluster_fit(formula, fit)
        models[name] = model

        if name == "density_only":
            block = {"wald_chi2": np.nan, "df": 0, "p": np.nan}
        elif name == "density_plus_growth":
            block = wald_block(
                model,
                ["population_density_change_pct_per_year"],
            )
        else:
            block = wald_block(
                model,
                [
                    "population_density_change_pct_per_year",
                    "log_density_c:population_density_change_pct_per_year",
                ],
            )

        rows.append(
            {
                "model": name,
                "formula": formula,
                "n": int(model.nobs),
                "r2": float(model.rsquared),
                "aic": float(model.aic),
                "bic": float(model.bic),
                "demographic_block_wald_chi2": block["wald_chi2"],
                "demographic_block_df": block["df"],
                "demographic_block_p": block["p"],
            }
        )

    return models, pd.DataFrame(rows)


def build_nonoverlap_holdout(
    data: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = config["analysis"]["future_radiance_response"]
    rel = int(cfg["nuts_release"])
    h = int(cfg["temporal_holdout_horizon_years"])
    train_start = int(cfg["temporal_holdout_train_start"])
    test_start = int(cfg["temporal_holdout_test_start"])
    min_country = int(cfg["min_country_regions"])

    train = make_interval_pairs(
        data, [train_start], h, rel, min_country
    )
    test = make_interval_pairs(
        data, [test_start], h, rel, min_country
    )

    common_countries = sorted(
        set(train["CNTR_CODE"]).intersection(test["CNTR_CODE"])
    )
    train = train[train["CNTR_CODE"].isin(common_countries)].copy()
    test = test[test["CNTR_CODE"].isin(common_countries)].copy()

    qlo, qhi = [
        float(x) for x in cfg["population_change_trim_quantiles"]
    ]
    lo = float(
        train["population_density_change_pct_per_year"].quantile(qlo)
    )
    hi = float(
        train["population_density_change_pct_per_year"].quantile(qhi)
    )

    train = train[
        train["population_density_change_pct_per_year"].between(lo, hi)
    ].copy()
    # Do not trim the holdout by its own distribution. Keep only values inside
    # the historical training support so this is interpolation, not accidental
    # extrapolation.
    test = test[
        test["population_density_change_pct_per_year"].between(lo, hi)
    ].copy()

    center = float(train["baseline_log10_population_density"].median())
    for frame in (train, test):
        frame["log_density_c"] = (
            frame["baseline_log10_population_density"] - center
        )
        frame["country"] = frame["CNTR_CODE"].astype(str)
        frame["window"] = frame["window"].astype(str)

    train.attrs["log_density_center"] = center
    return train, test


def holdout_formula(name: str) -> str:
    # A holdout contains only one interval in train and one in test, so no
    # window fixed effect can be estimated/transported. Country effects remain.
    common = (
        "annual_log10_radiance_change ~ C(country) + "
        "log_density_c + I(log_density_c ** 2) + I(log_density_c ** 3)"
    )
    if name == "density_only":
        return common
    if name == "density_plus_growth":
        return common + " + population_density_change_pct_per_year"
    return (
        common
        + " + population_density_change_pct_per_year"
        + " + log_density_c:population_density_change_pct_per_year"
    )


def temporal_holdout_metrics(
    train: pd.DataFrame,
    test: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    # Restrict to countries actually represented in both frames.
    countries = sorted(
        set(train["country"]).intersection(test["country"])
    )
    train = train[train["country"].isin(countries)].copy()
    test = test[test["country"].isin(countries)].copy()

    for name in MODEL_FORMULAS:
        formula = holdout_formula(name)
        model = smf.ols(formula, data=train).fit()
        pred = np.asarray(model.predict(test), dtype=float)
        obs = test["annual_log10_radiance_change"].to_numpy(dtype=float)
        resid = pred - obs
        rmse = float(np.sqrt(np.mean(resid**2)))
        mae = float(np.mean(np.abs(resid)))
        bias = float(np.mean(resid))
        denom = float(np.sum((obs - np.mean(obs)) ** 2))
        predictive_r2 = (
            1.0 - float(np.sum(resid**2)) / denom
            if denom > 0
            else np.nan
        )

        rows.append(
            {
                "model": name,
                "train_interval": str(train["window"].iloc[0]),
                "test_interval": str(test["window"].iloc[0]),
                "n_train": len(train),
                "n_test": len(test),
                "n_countries": len(countries),
                "rmse_log10_per_year": rmse,
                "mae_log10_per_year": mae,
                "bias_log10_per_year": bias,
                "predictive_r2": predictive_r2,
            }
        )

    out = pd.DataFrame(rows)
    base_rmse = float(
        out.loc[
            out["model"] == "density_only",
            "rmse_log10_per_year",
        ].iloc[0]
    )
    out["rmse_improvement_vs_density_only_pct"] = (
        100.0
        * (base_rmse - out["rmse_log10_per_year"])
        / base_rmse
    )

    m1_rmse = float(
        out.loc[
            out["model"] == "density_plus_growth",
            "rmse_log10_per_year",
        ].iloc[0]
    )
    out["rmse_improvement_vs_growth_model_pct"] = (
        100.0
        * (m1_rmse - out["rmse_log10_per_year"])
        / m1_rmse
    )
    return out


def choose_model(
    model_table: pd.DataFrame,
    holdout: pd.DataFrame,
    config: dict,
) -> tuple[str, list[str]]:
    cfg = config["analysis"]["future_radiance_response"]
    alpha = float(cfg["selection_alpha"])
    min_vs_base = float(cfg["min_holdout_rmse_improvement_pct"])
    min_interaction = float(
        cfg["min_interaction_rmse_improvement_pct"]
    )

    by_model = model_table.set_index("model")
    hv = holdout.set_index("model")
    reasons = []

    m1_ok = (
        float(
            by_model.loc[
                "density_plus_growth", "demographic_block_p"
            ]
        )
        < alpha
        and float(
            hv.loc[
                "density_plus_growth",
                "rmse_improvement_vs_density_only_pct",
            ]
        )
        >= min_vs_base
    )

    m2_ok = (
        float(
            by_model.loc[
                "density_growth_interaction",
                "demographic_block_p",
            ]
        )
        < alpha
        and float(
            hv.loc[
                "density_growth_interaction",
                "rmse_improvement_vs_density_only_pct",
            ]
        )
        >= min_vs_base
        and float(
            hv.loc[
                "density_growth_interaction",
                "rmse_improvement_vs_growth_model_pct",
            ]
        )
        >= min_interaction
    )

    if m2_ok:
        selected = "density_growth_interaction"
        reasons.append(
            "Selected interaction model: demographic block is significant "
            "and temporal holdout RMSE improves over both simpler models."
        )
    elif m1_ok:
        selected = "density_plus_growth"
        reasons.append(
            "Selected additive population-growth model: demographic term is "
            "significant and temporal holdout RMSE improves sufficiently over "
            "density-only."
        )
    else:
        selected = "density_only"
        reasons.append(
            "Selected density-only model: demographic terms did not satisfy "
            "both inferential and temporal-prediction criteria."
        )

    reasons.append(
        f"Selection alpha={alpha:g}; minimum RMSE improvement vs density-only="
        f"{min_vs_base:g}%; interaction must additionally improve over M1 by "
        f"{min_interaction:g}%."
    )
    return selected, reasons


def average_prediction(
    result,
    fit: pd.DataFrame,
    density: float,
    pop_change: float,
) -> tuple[float, float]:
    center = float(fit.attrs["log_density_center"])
    combos = (
        fit.groupby(["country", "window"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    combos["weight"] = combos["n"] / combos["n"].sum()
    pred = combos[["country", "window"]].copy()
    pred["log_density_c"] = np.log10(float(density)) - center
    pred["population_density_change_pct_per_year"] = float(pop_change)

    # Average prediction over the observed country-window mixture. Variance of
    # this population-average linear predictor is obtained from the model
    # design matrix rather than from observation-level prediction intervals.
    from patsy import dmatrix

    rhs = result.model.formula.split("~", 1)[1]
    X = dmatrix(rhs, pred, return_type="dataframe")
    X = X.reindex(columns=result.params.index, fill_value=0.0)
    xbar = np.average(
        X.to_numpy(dtype=float),
        axis=0,
        weights=combos["weight"].to_numpy(dtype=float),
    )
    beta = result.params.to_numpy(dtype=float)
    cov = result.cov_params().to_numpy(dtype=float)
    est = float(xbar @ beta)
    se = float(np.sqrt(max(0.0, xbar @ cov @ xbar)))
    return est, se


def build_transfer_table(
    model,
    fit: pd.DataFrame,
    selected: str,
    config: dict,
) -> pd.DataFrame:
    """Build the canonical transfer function used by future projections.

    If demographic change is not retained, use the full 2013-2024 longitudinal
    mixed-model response from density-temporal-change rather than the shorter
    five-year-pair density-only fit. The five-year pair analysis then serves
    only as the model-selection test for whether g_D is needed.

    If demographic change is retained, use the validated five-year pair model
    because g_D is defined from interval change.
    """
    processed = configured_path(config, "processed")

    if selected == "density_only":
        path = processed / "continuous_density_temporal_trend.csv"
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {path}. Run density-temporal-change first."
            )
        src = pd.read_csv(path)
        required = [
            "population_density",
            "log10_radiance_trend_per_year",
            "trend_se",
            "trend_ci_low",
            "trend_ci_high",
            "pct_radiance_change_per_year",
        ]
        missing = [x for x in required if x not in src.columns]
        if missing:
            raise KeyError(
                "Continuous density temporal trend is missing columns: "
                f"{missing}"
            )

        out = pd.DataFrame(
            {
                "selected_model": "density_only",
                "response_source": (
                    "full_2013_2024_longitudinal_mixed_model"
                ),
                "population_density": src["population_density"].astype(float),
                "population_density_change_pct_per_year": np.nan,
                "annual_log10_calibrated_viirs_change": src[
                    "log10_radiance_trend_per_year"
                ].astype(float),
                "annual_log10_change_se": src["trend_se"].astype(float),
                "annual_log10_change_ci_low": src["trend_ci_low"].astype(float),
                "annual_log10_change_ci_high": src["trend_ci_high"].astype(float),
                "calibrated_viirs_change_pct_per_year": src[
                    "pct_radiance_change_per_year"
                ].astype(float),
            }
        )
        out["calibrated_viirs_change_pct_per_year_ci_low"] = (
            100.0 * (10.0 ** out["annual_log10_change_ci_low"] - 1.0)
        )
        out["calibrated_viirs_change_pct_per_year_ci_high"] = (
            100.0 * (10.0 ** out["annual_log10_change_ci_high"] - 1.0)
        )
        return out.sort_values("population_density").reset_index(drop=True)

    cfg = config["analysis"]["future_radiance_response"]
    lo = float(
        fit["population_density_start"].quantile(
            float(cfg["projection_density_quantiles"][0])
        )
    )
    hi = float(
        fit["population_density_start"].quantile(
            float(cfg["projection_density_quantiles"][1])
        )
    )
    densities = np.geomspace(
        lo, hi, int(cfg["projection_density_grid_points"])
    )
    growth_values = [
        float(x) for x in cfg["projection_population_change_scenarios"]
    ]

    rows = []
    for growth in growth_values:
        for density in densities:
            est, se = average_prediction(
                model, fit, density, growth
            )
            rows.append(
                {
                    "selected_model": selected,
                    "response_source": (
                        "fixed_NUTS2016_five_year_interval_model"
                    ),
                    "population_density": float(density),
                    "population_density_change_pct_per_year": growth,
                    "annual_log10_calibrated_viirs_change": est,
                    "annual_log10_change_se": se,
                    "annual_log10_change_ci_low": est - 1.96 * se,
                    "annual_log10_change_ci_high": est + 1.96 * se,
                    "calibrated_viirs_change_pct_per_year": (
                        100.0 * (10.0**est - 1.0)
                    ),
                    "calibrated_viirs_change_pct_per_year_ci_low": (
                        100.0 * (10.0 ** (est - 1.96 * se) - 1.0)
                    ),
                    "calibrated_viirs_change_pct_per_year_ci_high": (
                        100.0 * (10.0 ** (est + 1.96 * se) - 1.0)
                    ),
                }
            )
    return pd.DataFrame(rows)


def load_headline_density_trends(config: dict) -> pd.DataFrame:
    processed = configured_path(config, "processed")
    path = processed / "density_band_temporal_trends.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run density-temporal-change first."
        )
    out = pd.read_csv(path)
    required = [
        "density_band",
        "log10_radiance_trend_per_year",
        "trend_se",
        "trend_p",
        "pct_radiance_change_per_year",
    ]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise KeyError(
            f"Headline density table missing columns: {missing}"
        )

    years = (
        int(config["analysis"]["density_temporal_change"]["last_year"])
        - int(config["analysis"]["density_temporal_change"]["first_year"])
    )
    out["pct_change_full_interval"] = 100.0 * (
        10.0
        ** (
            out["log10_radiance_trend_per_year"].astype(float)
            * years
        )
        - 1.0
    )
    out["pct_change_per_year_ci_low"] = 100.0 * (
        10.0
        ** (
            out["log10_radiance_trend_per_year"]
            - 1.96 * out["trend_se"]
        )
        - 1.0
    )
    out["pct_change_per_year_ci_high"] = 100.0 * (
        10.0
        ** (
            out["log10_radiance_trend_per_year"]
            + 1.96 * out["trend_se"]
        )
        - 1.0
    )
    return out


def make_figure(
    transfer: pd.DataFrame,
    headline: pd.DataFrame,
    selected: str,
    output: Path,
    dpi: int,
) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 6.0))

    if selected == "density_only":
        g = transfer.sort_values("population_density")
        ax.plot(
            g["population_density"],
            g["calibrated_viirs_change_pct_per_year"],
            linewidth=2.0,
            label="Canonical continuous response",
        )
        ax.fill_between(
            g["population_density"],
            g["calibrated_viirs_change_pct_per_year_ci_low"],
            g["calibrated_viirs_change_pct_per_year_ci_high"],
            alpha=0.2,
        )
    else:
        for growth, g in transfer.groupby(
            "population_density_change_pct_per_year",
            observed=True,
        ):
            g = g.sort_values("population_density")
            ax.plot(
                g["population_density"],
                g["calibrated_viirs_change_pct_per_year"],
                linewidth=1.5,
                label=f"{growth:+g}% yr⁻¹ density change",
            )

    band_x = {
        "<30": 15.0,
        "30-100": np.sqrt(30.0 * 100.0),
        "100-300": np.sqrt(100.0 * 300.0),
        "300-1000": np.sqrt(300.0 * 1000.0),
        ">=1000": 1800.0,
    }
    hx = [
        band_x.get(str(b), np.nan)
        for b in headline["density_band"]
    ]
    ax.scatter(
        hx,
        headline["pct_radiance_change_per_year"],
        s=45,
        marker="o",
        label="Headline density-band estimates",
    )

    ax.axhline(0, linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel("Population density (people km$^{-2}$)")
    ax.set_ylabel(
        "Annual change in calibrated VIIRS DNB radiance (%)"
    )
    ax.set_title(
        "Canonical historical response for future radiance projection"
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["future_radiance_response"]
    data = load_base_data(config)

    five_year_pairs = make_interval_pairs(
        data,
        [int(x) for x in cfg["five_year_start_years"]],
        int(cfg["five_year_horizon_years"]),
        int(cfg["nuts_release"]),
        int(cfg["min_country_regions"]),
    )
    fit = trim_and_center(
        five_year_pairs,
        tuple(
            float(x)
            for x in cfg["population_change_trim_quantiles"]
        ),
    )

    models, model_table = fit_candidate_models(fit)
    train, test = build_nonoverlap_holdout(data, config)
    holdout = temporal_holdout_metrics(train, test)
    selected, reasons = choose_model(
        model_table, holdout, config
    )

    transfer = build_transfer_table(
        models[selected], fit, selected, config
    )
    headline = load_headline_density_trends(config)

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "future_radiance_response"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    paths = {
        "pairs": processed / "future_radiance_historical_5yr_pairs.csv",
        "models": processed / "future_radiance_candidate_models.csv",
        "holdout": processed / "future_radiance_temporal_holdout.csv",
        "decision": processed / "future_radiance_model_decision.json",
        "transfer": processed / "future_radiance_canonical_response.csv",
        "headline": processed / "future_radiance_headline_density_classes.csv",
        "summary": processed / "future_radiance_response_summary.txt",
    }

    five_year_pairs.to_csv(paths["pairs"], index=False)
    model_table.to_csv(paths["models"], index=False)
    holdout.to_csv(paths["holdout"], index=False)
    transfer.to_csv(paths["transfer"], index=False)
    headline.to_csv(paths["headline"], index=False)

    decision = {
        "selected_model": selected,
        "selected_formula": MODEL_FORMULAS[selected],
        "canonical_response_source": (
            "full_2013_2024_longitudinal_mixed_model"
            if selected == "density_only"
            else "fixed_NUTS2016_five_year_interval_model"
        ),
        "calibration_version": config["calibration"]["version"],
        "calibration_method": config["calibration"]["method"],
        "nuts_release_for_model_selection": int(cfg["nuts_release"]),
        "selection_reasons": reasons,
    }
    with paths["decision"].open("w", encoding="utf-8") as h:
        json.dump(decision, h, indent=2)

    fig_path = figures / "canonical_future_radiance_response.pdf"
    make_figure(
        transfer,
        headline,
        selected,
        fig_path,
        int(cfg["figure_dpi"]),
    )

    with paths["summary"].open("w", encoding="utf-8") as h:
        h.write(
            "CANONICAL HISTORICAL RESPONSE FOR FUTURE CALIBRATED VIIRS DNB "
            "RADIANCE\n"
        )
        h.write("=" * 78 + "\n\n")
        h.write(
            f"Production calibration: {config['calibration']['version']}\n"
        )
        h.write(
            f"Calibration method: {config['calibration']['method']}\n"
        )
        h.write(
            f"Model-selection geography: NUTS {cfg['nuts_release']}\n"
        )
        h.write(
            "Five-year fitting windows: "
            + ", ".join(
                f"{y}-{y+int(cfg['five_year_horizon_years'])}"
                for y in cfg["five_year_start_years"]
            )
            + "\n"
        )
        h.write(
            f"Five-year pairs: {len(five_year_pairs)} before trim; "
            f"{len(fit)} after trim; {fit['NUTS_ID'].nunique()} NUTS3; "
            f"{fit['CNTR_CODE'].nunique()} countries.\n"
        )
        h.write(
            f"Population-density change support after trim: "
            f"{fit.attrs['trim_low']:.6g} to "
            f"{fit.attrs['trim_high']:.6g}% yr-1.\n\n"
        )

        h.write("CANDIDATE MODELS\n")
        h.write(
            model_table.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nNON-OVERLAPPING TEMPORAL HOLDOUT\n")
        h.write(
            holdout.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nMODEL DECISION\n")
        h.write(f"Selected model: {selected}\n")
        h.write(
            "Canonical response source: "
            + (
                "full 2013-2024 longitudinal mixed model"
                if selected == "density_only"
                else "fixed-NUTS2016 five-year interval model"
            )
            + "\n"
        )
        for reason in reasons:
            h.write(f"- {reason}\n")

        h.write("\nHEADLINE 2013-2024 DENSITY-CLASS TRENDS\n")
        h.write(
            headline[
                [
                    "density_band",
                    "pct_radiance_change_per_year",
                    "pct_change_per_year_ci_low",
                    "pct_change_per_year_ci_high",
                    "trend_p",
                    "pct_change_full_interval",
                ]
            ].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nINTERPRETATION\n")
        h.write(
            "- These five density-class estimates remain the descriptive "
            "Europe-wide historical result.\n"
        )
        h.write(
            "- The selected canonical model is used only as the continuous "
            "transfer function for future population projections.\n"
        )
        h.write(
            "- Population-density change is retained only if it passes both "
            "clustered inference and non-overlapping temporal prediction "
            "criteria.\n"
        )
        h.write(
            "- Habitability is not fitted here. It enters later as a spatial "
            "constraint on where future population/radiance change can occur.\n"
        )
        h.write(
            "- ISS spectral change is not fitted here. It remains a separate "
            "mechanism/sensitivity analysis and later biological spectral "
            "conversion.\n"
        )

    print("\nMODEL DECISION")
    print(f"Selected: {selected}")
    for reason in reasons:
        print(" ", reason)

    print("\nCANDIDATE MODELS")
    print(
        model_table.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print("\nNON-OVERLAPPING TEMPORAL HOLDOUT")
    print(
        holdout.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print("\nHEADLINE DENSITY-CLASS TRENDS")
    print(
        headline[
            [
                "density_band",
                "pct_radiance_change_per_year",
                "pct_change_per_year_ci_low",
                "pct_change_per_year_ci_high",
                "trend_p",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )

    print("\nWROTE")
    for p in [*paths.values(), fig_path]:
        print(p)


if __name__ == "__main__":
    main()
