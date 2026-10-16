#!/usr/bin/env python3
"""Full annual corrected VIIRS density model and NL333 anomaly audit, 2014-2024.

PRIMARY INPUT: all 11 annual observations for each of the same 1,514 fixed
NUTS2016 regions. The production 2017+ pixel-level 0.235 nW correction and
the 2021 log-ratio snow/surface correction are BOTH included in the model
response. Other 2021 series are only used for an explicitly labelled audit.

Fit annual log10 radiance with a NUTS3 random intercept and country intercept:
  M0: one common European linear time slope
  M1: density-band-specific linear time slopes (the primary scientific model)
Both are fitted to ALL YEARS, not to 2014/2024 endpoint differences.
Observed year-to-year transitions and yearly fitted-vs-observed trajectories
are reported. Chronological holdouts test generalisation to later years,
including corrected 2021 as a test year and then as a training year.

NL333 remains in every fit. A sensitivity removes it ONLY FROM SCORING.
No habitability H, country-specific temporal slopes, or snow-free fitted
data are included. This is not a causal attribution of greenhouse lighting.
"""
from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_density_temporal_change import fit_mixedlm, lrt
from project_config import configured_path, load_config
from select_fixed_nuts2016_corrected_temporal_model import load_corrected_panel

STAGE = "fixed_nuts2016_corrected_annual_nl333"
REGION = "NL333"
FORMULAS = {
    "common_slope": "log_rad ~ year_c + C(density_band) + C(country)",
    "density_slopes": "log_rad ~ year_c * C(density_band) + C(country)",
}
# Cutoffs are strictly earlier than their test years. The later training fold
# includes the CORRECTED 2021 observation.
FOLDS = {
    2018: [2019, 2020],
    2020: [2021, 2022, 2023, 2024],
    2022: [2023, 2024],
}


def load_data(config: dict):
    cfg = config["analysis"]["fixed_nuts2016_corrected_temporal_model"]
    if (int(cfg["first_year"]), int(cfg["last_year"])) != (2014, 2024):
        raise RuntimeError("Expected the locked full 2014-2024 analysis.")
    expected = int(
        config["analysis"]["fixed_nuts2016_2021_diagnostic"]["expected_n_nuts3"]
    )
    panel, _ = load_corrected_panel(config, cfg)
    if panel["NUTS_ID"].nunique() != expected or len(panel) != expected * 11:
        raise RuntimeError("Changed/missing fixed 1,514-region complete panel.")
    if panel.duplicated(["NUTS_ID", "year"]).any():
        raise RuntimeError("Duplicate NUTS3/year in corrected panel.")
    if set(panel["year"].unique()) != set(range(2014, 2025)):
        raise RuntimeError("Incomplete year coverage.")
    panel = panel.sort_values(["NUTS_ID", "year"]).copy()
    panel["year_c"] = panel["year"].astype(float) - 2019.0
    panel["log_rad"] = np.log10(panel["radiance_mean_corrected"].astype(float))
    panel["density_band"] = panel["density_band"].cat.remove_unused_categories()
    if REGION not in set(panel["NUTS_ID"]):
        raise RuntimeError("NL333 not in locked fixed NUTS2016 cohort.")

    history_path = (
        configured_path(config, "processed") / str(cfg["corrected_history_file"])
    )
    history = pd.read_parquet(history_path)
    columns = {"NUTS_ID", "year", "radiance_mean_raw",
               "radiance_mean_corrected", str(cfg["corrected_radiance_column"])}
    if not columns.issubset(history.columns):
        raise KeyError(f"Correction history missing {columns - set(history.columns)}")
    h = history[list(columns)].copy()
    h["NUTS_ID"] = h["NUTS_ID"].astype(str)
    h = h.rename(columns={
        "radiance_mean_corrected": "instrument_calibrated_nw",
        "radiance_mean_raw": "uncalibrated_raw_nw",
        str(cfg["corrected_radiance_column"]): "instrument_and_2021_corrected_nw",
    })
    if h.duplicated(["NUTS_ID", "year"]).any():
        raise RuntimeError("Duplicated original surface-correction history.")
    annual = panel[[
        "NUTS_ID", "CNTR_CODE", "density_band",
        "reference_population_density", "year", "radiance_mean_corrected"
    ]].merge(
        h, on=["NUTS_ID", "year"], how="left", validate="one_to_one"
    )
    fields = [
        "instrument_calibrated_nw",
        "uncalibrated_raw_nw",
        "instrument_and_2021_corrected_nw",
    ]
    if not np.isfinite(annual[fields].to_numpy(float)).all():
        raise RuntimeError("Missing input radiance columns in corrected history.")
    if not np.allclose(
        annual["radiance_mean_corrected"],
        annual["instrument_and_2021_corrected_nw"], atol=1e-9, rtol=1e-9
    ):
        raise RuntimeError("Primary model is NOT using the 2021-corrected input.")
    other_years = annual["year"].ne(2021)
    if not np.allclose(
        annual.loc[other_years, "instrument_calibrated_nw"],
        annual.loc[other_years, "radiance_mean_corrected"], atol=1e-9, rtol=1e-9
    ):
        raise RuntimeError("2021 surface correction has affected non-2021 years.")
    corrected_2021 = annual["year"].eq(2021) & (
        (annual["instrument_calibrated_nw"] -
         annual["radiance_mean_corrected"]).abs() > 1e-8
    )
    if not corrected_2021.any():
        raise RuntimeError("No nonzero 2021 correction found. Wrong dataset?")
    # Ensure the inherited pixel-calibration field is distinguished from raw.
    if not (
        (annual.loc[annual["year"].ge(2017), "uncalibrated_raw_nw"] -
         annual.loc[annual["year"].ge(2017), "instrument_calibrated_nw"]) > 0
    ).any():
        raise RuntimeError("2017+ instrument calibration not detected.")
    counts = {
        "n_nuts3": expected, "n_countries": int(panel["CNTR_CODE"].nunique()),
        "n_rows": len(panel),
        "n_2021_regions_changed": int(corrected_2021.sum()),
    }
    return panel, annual, counts


def fit_annual(panel: pd.DataFrame):
    models, fit_rows = {}, []
    for name, formula in FORMULAS.items():
        model, optimizer = fit_mixedlm(formula, panel)
        models[name] = model
        fit_rows.append({
            "model": name, "formula": formula, "n_obs": int(model.nobs),
            "n_regions": int(panel["NUTS_ID"].nunique()),
            "log_likelihood": float(model.llf),
            "aic": float(model.aic), "bic": float(model.bic),
            "optimizer": optimizer,
        })
    return models, pd.DataFrame(fit_rows)


def endpoint_from_anchor(
    model, anchor: pd.DataFrame, target_year: int, observations: pd.DataFrame
) -> pd.DataFrame:
    first = int(anchor["year"].iloc[0])
    if not anchor["year"].eq(first).all():
        raise RuntimeError("Anchor contains more than one year.")
    pred = anchor.copy()
    pred["year"] = int(target_year)
    pred["year_c"] = float(target_year - 2019)
    ratio_log = (
        np.asarray(model.predict(pred), dtype=float)
        - np.asarray(model.predict(anchor), dtype=float)
    ) if model is not None else np.zeros(len(pred), dtype=float)
    start = anchor["radiance_mean_corrected"].to_numpy(float)
    prediction = start * np.power(10.0, ratio_log)
    out = anchor[[
        "NUTS_ID", "CNTR_CODE", "density_band",
        "reference_population_density",
    ]].copy()
    out["anchor_year"] = first
    out["target_year"] = target_year
    out["start_radiance_nw"] = start
    out["predicted_target_nw"] = prediction
    out["predicted_delta_log10"] = ratio_log
    out["predicted_delta_nw"] = prediction - start
    out["predicted_delta_pct"] = 100.0 * (prediction / start - 1.0)
    actual = observations.loc[
        observations["year"].eq(target_year),
        ["NUTS_ID", "radiance_mean_corrected"],
    ].rename(columns={"radiance_mean_corrected": "observed_target_nw"})
    out = out.merge(actual, on="NUTS_ID", validate="one_to_one")
    obs = out["observed_target_nw"].to_numpy(float)
    out["observed_delta_nw"] = obs - start
    out["observed_delta_pct"] = 100.0 * (obs / start - 1.0)
    out["observed_delta_log10"] = np.log10(obs / start)
    out["error_nw"] = out["predicted_delta_nw"] - out["observed_delta_nw"]
    out["error_pp"] = out["predicted_delta_pct"] - out["observed_delta_pct"]
    out["error_log10"] = (
        out["predicted_delta_log10"] - out["observed_delta_log10"]
    )
    return out


def score(g: pd.DataFrame) -> dict:
    result = {"n_predictions": len(g), "n_nuts3": g["NUTS_ID"].nunique()}
    for key in ("nw", "pp", "log10"):
        values = g[f"error_{key}"].to_numpy(float)
        result[f"rmse_{key}"] = float(np.sqrt(np.mean(values ** 2)))
        result[f"mae_{key}"] = float(np.mean(np.abs(values)))
        result[f"median_ae_{key}"] = float(np.median(np.abs(values)))
        result[f"bias_{key}"] = float(np.mean(values))
    return result


def group_scores(d: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    rows = []
    for idx, g in d.groupby(cols, observed=True):
        key = idx if isinstance(idx, tuple) else (idx,)
        rows.append({**dict(zip(cols, key)), **score(g)})
    return pd.DataFrame(rows)


def validate_forecasts(panel: pd.DataFrame):
    all_results, fits = [], []
    for cutoff, target_years in FOLDS.items():
        train = panel.loc[panel["year"].le(cutoff)].copy()
        fitted, comparison = fit_annual(train)
        comparison.insert(0, "train_through", cutoff)
        fits.append(comparison)
        anchor = panel.loc[panel["year"].eq(cutoff)].copy()
        for target_year in target_years:
            for name in ["persistence", *FORMULAS]:
                p = endpoint_from_anchor(
                    fitted.get(name), anchor, target_year, panel
                )
                p["model"] = name
                p["train_through"] = cutoff
                p["horizon_years"] = target_year - cutoff
                all_results.append(p)
    predictions = pd.concat(all_results, ignore_index=True)
    expected_rows = panel["NUTS_ID"].nunique() * sum(
        len(v) for v in FOLDS.values()
    ) * (len(FORMULAS) + 1)
    if len(predictions) != expected_rows:
        raise RuntimeError("Unexpected number of chronological forecasts.")
    predictions["density_band"] = predictions["density_band"].astype(str)
    by_year = group_scores(
        predictions, ["model", "train_through", "target_year", "horizon_years"]
    )
    by_density = group_scores(predictions, ["model", "density_band"])
    overall = group_scores(predictions, ["model"])
    sensitivity = group_scores(
        predictions.loc[predictions["NUTS_ID"].ne(REGION)], ["model"]
    )
    sensitivity["population"] = "all_except_NL333_scoring_only"
    overall["population"] = "all_1514"
    overall = pd.concat([overall, sensitivity], ignore_index=True)
    return pd.concat(fits, ignore_index=True), predictions, by_year, by_density, overall


def annual_differences(annual: pd.DataFrame):
    d = annual.sort_values(["NUTS_ID", "year"]).copy()
    d["previous_year"] = d.groupby("NUTS_ID")["year"].shift()
    for name, col in [
        ("corrected", "radiance_mean_corrected"),
        ("instrument", "instrument_calibrated_nw"),
        ("raw", "uncalibrated_raw_nw"),
    ]:
        previous = d.groupby("NUTS_ID")[col].shift()
        d[f"annual_delta_{name}_nw"] = d[col] - previous
        d[f"annual_delta_{name}_pct"] = 100.0 * (d[col] / previous - 1.0)
    d["annual_delta_log10"] = (
        np.log10(d["radiance_mean_corrected"])
        - np.log10(d.groupby("NUTS_ID")["radiance_mean_corrected"].shift())
    )
    result = d.loc[d["year"].ge(2015)].copy()
    if not result["previous_year"].add(1).eq(result["year"]).all():
        raise RuntimeError("Gaps in the observed annual transitions.")
    density_year = (
        result.groupby(["density_band", "year"], observed=True)
        .agg(
            n_nuts3=("NUTS_ID", "nunique"),
            mean_log10_annual_change=("annual_delta_log10", "mean"),
            median_annual_change_nw=("annual_delta_corrected_nw", "median"),
            median_annual_change_pct=("annual_delta_corrected_pct", "median"),
        ).reset_index()
    )
    density_year["geometric_mean_annual_change_pct"] = 100.0 * (
        np.power(10.0, density_year["mean_log10_annual_change"]) - 1.0
    )
    return result, density_year


def full_annual_trajectories(panel: pd.DataFrame, models: dict):
    anchor = panel.loc[panel["year"].eq(2014)].copy()
    all_rows = []
    for year in range(2014, 2025):
        observed = panel.loc[panel["year"].eq(year)].copy()
        for name, model in models.items():
            pred = endpoint_from_anchor(model, anchor, year, panel)
            pred["model"] = name
            all_rows.append(pred)
    region = pd.concat(all_rows, ignore_index=True)
    region["density_band"] = region["density_band"].astype(str)
    summary = []
    for (name, band, year), g in region.groupby(
        ["model", "density_band", "target_year"], observed=True
    ):
        summary.append({
            "model": name, "density_band": band, "year": year,
            "n_nuts3": len(g),
            "observed_geometric_mean_nw": float(
                np.power(10.0, np.log10(g["observed_target_nw"]).mean())
            ),
            "predicted_geometric_mean_nw": float(
                np.power(10.0, np.log10(g["predicted_target_nw"]).mean())
            ),
            "observed_mean_log10_delta_from_2014": float(
                g["observed_delta_log10"].mean()
            ),
            "predicted_mean_log10_delta_from_2014": float(
                g["predicted_delta_log10"].mean()
            ),
        })
    tr = pd.DataFrame(summary)
    for metric in ("observed", "predicted"):
        baseline = tr.loc[tr["year"].eq(2014)].set_index(
            ["model", "density_band"]
        )[f"{metric}_geometric_mean_nw"]
        tr[f"{metric}_index_2014_100"] = [
            100.0 * val / baseline.loc[(mod, band)]
            for mod, band, val in zip(
                tr["model"], tr["density_band"],
                tr[f"{metric}_geometric_mean_nw"]
            )
        ]
    return region, tr


def make_nl333_report(annual: pd.DataFrame, changes: pd.DataFrame,
                      full_region: pd.DataFrame, validation: pd.DataFrame):
    nl = annual.loc[annual["NUTS_ID"].eq(REGION)].sort_values("year").copy()
    nl["difference_due_to_2021_surface_nw"] = (
        nl["radiance_mean_corrected"] - nl["instrument_calibrated_nw"]
    )
    nl["year_index_2014_100"] = 100.0 * (
        nl["radiance_mean_corrected"] /
        float(nl.loc[nl["year"].eq(2014), "radiance_mean_corrected"].iloc[0])
    )
    change_cols = [
        "year", "annual_delta_corrected_nw", "annual_delta_corrected_pct",
        "annual_delta_instrument_nw", "annual_delta_raw_nw",
    ]
    nl = nl.merge(
        changes.loc[changes["NUTS_ID"].eq(REGION), change_cols],
        on="year", how="left", validate="one_to_one"
    )
    national = (
        changes.loc[changes["CNTR_CODE"].eq("NL")].groupby("year")
        .agg(
            netherlands_median_annual_change_nw=("annual_delta_corrected_nw", "median"),
            netherlands_median_annual_change_pct=("annual_delta_corrected_pct", "median"),
        ).reset_index()
    )
    peer = (
        changes.loc[
            changes["NUTS_ID"].str.startswith("NL33") &
            changes["NUTS_ID"].ne(REGION)
        ].groupby("year")
        .agg(
            zuid_holland_peer_median_annual_change_nw=("annual_delta_corrected_nw", "median"),
            zuid_holland_peer_median_annual_change_pct=("annual_delta_corrected_pct", "median"),
            n_zuid_holland_peers=("NUTS_ID", "nunique"),
        ).reset_index()
    )
    rank = changes.copy()
    rank["negative_change_rank_europe"] = (
        rank.groupby("year")["annual_delta_corrected_nw"]
        .rank(ascending=True, method="min")
    )
    rank["negative_change_rank_netherlands"] = (
        rank.groupby(["CNTR_CODE", "year"])["annual_delta_corrected_nw"]
        .rank(ascending=True, method="min")
    )
    nl = nl.merge(
        rank.loc[rank["NUTS_ID"].eq(REGION),
                 ["year", "negative_change_rank_europe",
                  "negative_change_rank_netherlands"]],
        on="year", how="left", validate="one_to_one"
    ).merge(national, on="year", how="left").merge(peer, on="year", how="left")

    model_full = full_region.loc[
        full_region["NUTS_ID"].eq(REGION) &
        full_region["model"].eq("density_slopes"),
        ["target_year", "predicted_target_nw", "predicted_delta_nw",
         "error_nw", "predicted_delta_pct"]
    ].rename(columns={"target_year": "year"})
    nl = nl.merge(model_full, on="year", validate="one_to_one")

    residuals = validation.loc[validation["NUTS_ID"].eq(REGION)].copy()
    by_case = []
    for (model, cutoff, year), g in validation.groupby(
        ["model", "train_through", "target_year"], observed=True
    ):
        error_energy = float(np.square(g["error_nw"].to_numpy(float)).sum())
        specific = g.loc[g["NUTS_ID"].eq(REGION)].iloc[0]
        by_case.append({
            "model": model, "train_through": cutoff, "target_year": year,
            "NL333_error_nw": float(specific["error_nw"]),
            "NL333_observed_change_nw": float(specific["observed_delta_nw"]),
            "NL333_predicted_change_nw": float(specific["predicted_delta_nw"]),
            "NL333_share_of_squared_error": (
                float(specific["error_nw"] ** 2 / error_energy)
                if error_energy > 0 else np.nan
            ),
        })
    shares = pd.DataFrame(by_case)
    concentration = []
    for name, g in validation.groupby("model", observed=True):
        per_region = (
            g.groupby("NUTS_ID")["error_nw"]
            .apply(lambda s: float(np.square(s.to_numpy(float)).sum()))
            .sort_values(ascending=False)
        )
        whole = float(per_region.sum())
        concentration.append({
            "model": name,
            "NL333_share_total_squared_error": (
                float(per_region.loc[REGION] / whole) if whole else np.nan
            ),
            "top_region": str(per_region.index[0]),
            "top_5_region_share": (
                float(per_region.head(5).sum() / whole) if whole else np.nan
            ),
            "top_10_region_share": (
                float(per_region.head(10).sum() / whole) if whole else np.nan
            ),
        })
    return nl, residuals, shares, pd.DataFrame(concentration)


def draw_plots(nl: pd.DataFrame, trajectories: pd.DataFrame,
               by_year: pd.DataFrame, labels: list[str],
               fig_dir, dpi: int):
    fig, ax = plt.subplots(figsize=(10.2, 5.8))
    ax.plot(nl["year"], nl["uncalibrated_raw_nw"], ":", marker="o",
            label="Original raw VIIRS mean")
    ax.plot(nl["year"], nl["instrument_calibrated_nw"], "--", marker="s",
            label="2017+ instrument-corrected")
    ax.plot(nl["year"], nl["radiance_mean_corrected"], "-", marker="o",
            label="Both corrections (model input)", linewidth=2.2)
    ax.plot(nl["year"], nl["predicted_target_nw"], "-.", marker="x",
            label="Annual density model (anchored to 2014)")
    ax.set_xticks(range(2014, 2025))
    ax.set_xlabel("Year")
    ax.set_ylabel("NUTS3 mean VIIRS radiance (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("NL333: all 11 annual observations and the fitted trajectory")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(fig_dir / "nl333_full_annual_radiance.pdf", dpi=dpi)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6))
    temp = nl.dropna(subset=["annual_delta_corrected_nw"])
    axes[0].plot(temp["year"], temp["annual_delta_corrected_nw"],
                 "o-", label="NL333")
    axes[0].plot(temp["year"],
                 temp["zuid_holland_peer_median_annual_change_nw"],
                 "--", label="Median other NL33 NUTS3")
    axes[0].axhline(0, color=".3", linewidth=.7)
    axes[0].set_ylabel("Change from previous year (nW cm$^{-2}$ sr$^{-1}$)")
    axes[0].set_xlabel("Ending year")
    axes[0].set_title("When does NL333 change?")
    axes[0].legend()
    axes[0].grid(alpha=0.2)
    report = by_year.loc[
        by_year["model"].isin(["density_slopes", "persistence"])
        & by_year["train_through"].eq(2020)
    ]
    for name, g in report.groupby("model"):
        axes[1].plot(g["target_year"], g["rmse_nw"], marker="o", label=name)
    axes[1].set_xlabel("Held-out year (anchored to 2020)")
    axes[1].set_ylabel("Europe-wide physical-radiance RMSE")
    axes[1].set_title("Forecast error as the horizon grows")
    axes[1].legend()
    axes[1].grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(fig_dir / "nl333_annual_changes_and_temporal_errors.pdf", dpi=dpi)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.5))
    for i, band in enumerate(labels):
        ax = axes.flat[i]
        g = trajectories.loc[trajectories["density_band"].eq(band)]
        for model, ls in [
            ("density_slopes", "--"), ("common_slope", ":")
        ]:
            s = g.loc[g["model"].eq(model)].sort_values("year")
            ax.plot(s["year"], s["predicted_index_2014_100"],
                    ls, label=model)
        obs = g.loc[g["model"].eq("density_slopes")].sort_values("year")
        ax.plot(obs["year"], obs["observed_index_2014_100"],
                "o-", label="Observed; BOTH corrections")
        ax.set_title(f"{band} people km$^{{-2}}$")
        ax.grid(alpha=.2)
    axes.flat[-1].axis("off")
    axes.flat[-1].legend(*axes.flat[0].get_legend_handles_labels(),
                         loc="center", frameon=False)
    fig.supylabel("Geometric-mean radiance index (2014 = 100)")
    fig.supxlabel("Year")
    fig.suptitle("2014–2024: full annual series with 2017+ and 2021 corrections")
    fig.tight_layout(rect=(.02, .02, 1, .96))
    fig.savefig(fig_dir / "observed_vs_modelled_every_year_by_density.pdf", dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()
    config = load_config(args.config)
    panel, annual, evidence = load_data(config)
    model, fit = fit_annual(panel)
    test_lr, test_df, test_p = lrt(
        model["common_slope"], model["density_slopes"]
    )
    training, holdout, by_year, by_density, all_scores = (
        validate_forecasts(panel)
    )
    annual_changes, annual_density = annual_differences(annual)
    full_region, trajectories = full_annual_trajectories(panel, model)
    nl, nl_holdout, nl_case, concentration = make_nl333_report(
        annual, annual_changes, full_region, holdout
    )

    processed = configured_path(config, "processed") / STAGE
    figures = configured_path(config, "figures") / "population_relationship" / STAGE
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    data = {
        "candidate_model_fit": fit,
        "chronological_holdout_training_fit": training,
        "annual_corrected_radiance": annual,
        "all_regions_year_to_year_changes": annual_changes,
        "density_year_to_year_summary": annual_density,
        "observed_vs_fitted_each_year_by_density": trajectories,
        "full_fit_annual_region_predictions": full_region,
        "chronological_holdout_predictions": holdout,
        "chronological_holdout_by_year": by_year,
        "chronological_holdout_by_density": by_density,
        "chronological_holdout_overall_and_NL333_sensitivity": all_scores,
        "nl333_annual_diagnosis": nl,
        "nl333_holdout_predictions": nl_holdout,
        "nl333_holdout_error_share_by_year": nl_case,
        "absolute_squared_error_concentration": concentration,
    }
    for filename, frame in data.items():
        frame.to_csv(processed / f"{filename}.csv", index=False)
    draw_plots(
        nl, trajectories, by_year,
        list(config["analysis"]["density_temporal_change"]["density_band_labels"]),
        figures, 300,
    )
    summary = processed / "corrected_annual_nl333_summary.txt"
    with summary.open("w", encoding="utf-8") as f:
        f.write("FULL 2014-2024 ANNUAL CORRECTED VIIRS MODEL AND NL333 AUDIT\n")
        f.write("=" * 78 + "\n\n")
        f.write(
            "SOURCE: fixed_nuts2016_2021_surface_correction/"
            "fixed_nuts2016_viirs_2014_2024_surface_corrected.parquet\n"
            "RESPONSE: radiance_mean_corrected_surface_logratio_primary\n"
            "INCLUDES: 2017+ pixel correction (0.235 nW cm^-2 sr^-1) AND "
            "2021 logratio snow/surface correction. All 2014-2024 observations "
            "enter the annual longitudinal MixedLM, NOT only endpoints.\n"
            "2014-2020 median density bands are locked; 1,514 fixed NUTS2016 "
            "regions; all observations retained, including NL333.\n"
            "Model M1: annual log10 radiance ~ year_c * fixed density band "
            "+ country intercept + NUTS3 random intercept.\n"
            "M0: common linear temporal slope for comparison.\n"
            "Observed year-on-year data and every annual fitted trajectory "
            "are exported; linear trends are not forced to match endpoints.\n"
            "Holdouts: train to 2018 -> 2019/20, 2020 -> 2021/22/23/24, "
            "2022 -> 2023/24. Corrected 2021 enters the second holdout's "
            "validation and the third holdout's training.\n"
            "The exclusion of NL333 is ONLY a score sensitivity; all fitted "
            "models include that region.\n"
            "No H, no country-specific slopes, no demographic predictor "
            "without harmonised fixed NUTS2016 demographic history.\n"
            "These statistics cannot establish a physical cause for NL333's "
            "change. Follow up with spatial/spectral analysis if indicated.\n\n"
        )
        for k, v in evidence.items():
            f.write(f"{k}={v}\n")
        f.write(
            f"\nM0 versus M1 density slope LRT: LR={test_lr:.8g}, "
            f"df={test_df}, p={test_p:.8g}\n\n"
        )
        for label, df in [
            ("FULL-ANNUAL MODEL FIT", fit),
            ("CHRONOLOGICAL HOLDOUT OVERALL & NL333 SCORING SENSITIVITY", all_scores),
            ("HOLDOUT BY TARGET YEAR", by_year),
            ("NL333 ANNUAL INPUT AND CHANGES", nl),
            ("NL333 ERROR SHARE BY YEAR", nl_case),
            ("ABSOLUTE SQUARED ERROR CONCENTRATION", concentration),
        ]:
            f.write(label + "\n" + df.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            ) + "\n\n")
    print("\n", summary, flush=True)
    print("\nMODEL FIT\n", fit.to_string(index=False), flush=True)
    print("\nHOLDOUT OVERALL\n", all_scores.to_string(index=False), flush=True)
    print("\nNL333 ANNUAL VALUES\n", nl[[
        "year", "radiance_mean_corrected",
        "annual_delta_corrected_nw", "annual_delta_corrected_pct",
        "predicted_target_nw"
    ]].to_string(index=False), flush=True)
    print("\nOUTPUTS\n")
    for name in data:
        print(processed / f"{name}.csv")
    for name in [
        "nl333_full_annual_radiance",
        "nl333_annual_changes_and_temporal_errors",
        "observed_vs_modelled_every_year_by_density",
    ]:
        print(figures / f"{name}.pdf")


if __name__ == "__main__":
    main()
