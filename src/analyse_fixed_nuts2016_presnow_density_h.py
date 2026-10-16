#!/usr/bin/env python3
"""Pre-snow VIIRS density + terrain habitability reconstruction, 2014-2024.

This is a new, non-destructive exploratory stage. It preserves:
 * production 2017+ pixel-level calibration;
 * the *pre-2021-surface-correction* radiance series;
 * the locked fixed-NUTS2016 2014-2024 cohort and 2014-2020 density classes;
 * H_core as an independent geographic attribute, not density/H.

Models are nested in increasing complexity. Selection is NOT automatic:
report likelihood/BIC, out-of-sample prediction in log, percentage points
and absolute nW, and explicitly labelled 2021/country-trend sensitivities.

Population-density change after 2020 cannot be assumed comparable in NUTS2016
from year-specific Eurostat data; it is deliberately left for a later stage.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_density_temporal_change import fit_mixedlm, lrt
from project_config import configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import prepare_analysis_panel


STAGE = "fixed_nuts2016_presnow_density_h"
FORMULAS = {
    "D_density": "log_rad ~ year_c * C(density_band) + C(country)",
    "DH_spatial": (
        "log_rad ~ year_c * C(density_band) + C(country)"
        " + h_c * C(density_band)"
    ),
    "DH_common_trend": (
        "log_rad ~ year_c * C(density_band) + C(country)"
        " + h_c * C(density_band) + year_c:h_c"
    ),
    "DH_density_trend": (
        "log_rad ~ year_c * C(density_band) * h_c + C(country)"
    ),
}
ORDER = list(FORMULAS)


def load_panel(config: dict, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    processed = configured_path(config, "processed")
    fixed_path = processed / "fixed_nuts2016_viirs_2014_2024.parquet"
    reference_path = processed / "fixed_nuts2016_population_density_reference.csv"
    habit_path = processed / config["outputs"]["nuts3_habitability"]
    for p in (fixed_path, reference_path, habit_path):
        if not p.is_file():
            raise FileNotFoundError(f"Missing prerequisite {p}")

    fixed = pd.read_parquet(fixed_path)
    # This specific input must NOT be the later 2021 snow/surface-adjusted file.
    if "radiance_mean_corrected" not in fixed:
        raise KeyError("Pre-snow fixed radiance column is missing.")
    if "nuts_release" in fixed and not fixed["nuts_release"].eq(2016).all():
        raise RuntimeError("Mixed geography in fixed NUTS2016 radiance history.")
    reference = pd.read_csv(reference_path)
    fixed["NUTS_ID"] = fixed["NUTS_ID"].astype(str)
    reference["NUTS_ID"] = reference["NUTS_ID"].astype(str)
    labels = list(config["analysis"]["density_temporal_change"]["density_band_labels"])
    reference["density_band"] = pd.Categorical(
        reference["density_band"].astype(str), categories=labels, ordered=True
    )
    if reference["density_band"].isna().any():
        raise RuntimeError("Unexpected density-band label in locked reference.")
    base = prepare_analysis_panel(fixed, reference, config)
    base = base[base["year"].between(2014, 2024)].copy()
    if base.duplicated(["NUTS_ID", "year"]).any():
        raise RuntimeError("Duplicated NUTS_ID/year in original fixed history.")
    expected = int(
        config["analysis"]["fixed_nuts2016_2021_diagnostic"]["expected_n_nuts3"]
    )
    if base["NUTS_ID"].nunique() != expected or len(base) != 11 * expected:
        raise RuntimeError(
            f"Locked original cohort expected {expected} regions x 11 years; "
            f"got {base['NUTS_ID'].nunique()} x {len(base)} annual rows."
        )
    if (base["radiance_mean_corrected"] <= 0).any():
        raise RuntimeError("Non-positive radiance in original panel.")
    endpoint_path = (
        processed / "fixed_nuts2016_corrected_endpoint_change"
        / "paired_nuts3_radiance_changes_2014_2024.csv"
    )
    if endpoint_path.exists():
        ids = set(pd.read_csv(endpoint_path, usecols=["NUTS_ID"])["NUTS_ID"].astype(str))
        if ids != set(base["NUTS_ID"]):
            raise RuntimeError("Cohort differs from 2014-2024 endpoint ground truth.")

    habit = pd.read_parquet(habit_path)
    needed = {"NUTS_ID", "nuts_release", "habitable_fraction_core"}
    if not needed.issubset(habit.columns):
        raise KeyError(f"H table missing: {sorted(needed - set(habit.columns))}")
    h = habit.loc[
        pd.to_numeric(habit["nuts_release"], errors="coerce").eq(2016),
        ["NUTS_ID", "habitable_fraction_core"],
    ].copy()
    h["NUTS_ID"] = h["NUTS_ID"].astype(str)
    if h["NUTS_ID"].duplicated().any():
        raise RuntimeError("H_core is not unique for fixed NUTS2016 regions.")
    h = h.rename(columns={"habitable_fraction_core": "H_core"})
    h["H_core"] = pd.to_numeric(h["H_core"], errors="coerce")
    d = base.merge(h, on="NUTS_ID", how="left", validate="many_to_one")
    valid_h = np.isfinite(d["H_core"]) & d["H_core"].between(0, 1)
    keep_ids = d.loc[valid_h, "NUTS_ID"].unique()
    n_h = len(keep_ids)
    min_fraction = float(cfg["minimum_h_coverage_fraction"])
    if n_h / expected < min_fraction:
        raise RuntimeError(
            f"H_core coverage {n_h}/{expected} below {min_fraction:.0%} "
            "minimum; investigate missing regions rather than silently "
            "changing the comparison cohort."
        )
    panel = d.loc[d["NUTS_ID"].isin(keep_ids)].copy()
    if len(panel) != 11 * n_h:
        raise RuntimeError("Adding habitability unbalanced the longitudinal cohort.")

    panel["country"] = panel["CNTR_CODE"].astype("category")
    panel["density_band"] = pd.Categorical(
        panel["density_band"], categories=labels, ordered=True
    )
    panel["year_c"] = panel["year"].astype(float) - 2019.0
    h_median = float(panel.drop_duplicates("NUTS_ID")["H_core"].median())
    panel["h_c"] = panel["H_core"] - h_median
    panel.attrs["h_median"] = h_median
    return base, panel


def fit_models(data: pd.DataFrame, formulas: dict[str, str] = FORMULAS):
    models, rows = {}, []
    for name, formula in formulas.items():
        result, optimizer = fit_mixedlm(formula, data)
        models[name] = result
        rows.append({
            "model": name, "formula": formula, "n_obs": int(result.nobs),
            "n_nuts3": int(data["NUTS_ID"].nunique()),
            "log_likelihood": float(result.llf),
            "aic": float(result.aic), "bic": float(result.bic),
            "optimizer": optimizer,
        })
    return models, pd.DataFrame(rows)


def nested_tests(models: dict, order: list[str]) -> pd.DataFrame:
    rows = []
    for reduced, full in zip(order[:-1], order[1:]):
        lr, df, p = lrt(models[reduced], models[full])
        rows.append({
            "reduced": reduced, "full": full,
            "lr": lr, "df": df, "p_value": p,
            "bic_improvement": float(models[reduced].bic - models[full].bic),
        })
    return pd.DataFrame(rows)


def predict_changes(
    result, anchor: pd.DataFrame, target_years: list[int]
) -> pd.DataFrame:
    """Anchor *each* region to its observed start radiance, predict FE change."""
    a = anchor.copy()
    year0 = int(a["year"].iloc[0])
    if not a["year"].eq(year0).all():
        raise RuntimeError("Prediction anchor must refer to one year.")
    r0 = a["radiance_mean_corrected"].to_numpy(float)
    log_fit0 = np.asarray(result.predict(a), dtype=float)
    out = []
    for year in target_years:
        b = a.copy()
        b["year"] = int(year)
        b["year_c"] = float(year - 2019)
        contrast = np.asarray(result.predict(b), dtype=float) - log_fit0
        ratio = np.power(10.0, contrast)
        r1 = r0 * ratio
        predicted = a[[
            "NUTS_ID", "CNTR_CODE", "density_band", "H_core",
            "reference_population_density",
        ]].copy()
        predicted["anchor_year"] = year0
        predicted["target_year"] = int(year)
        predicted["radiance_anchor_nw"] = r0
        predicted["predicted_log10_change"] = contrast
        predicted["predicted_change_percent"] = 100.0 * (ratio - 1.0)
        predicted["predicted_change_nw"] = r1 - r0
        predicted["predicted_radiance_target_nw"] = r1
        out.append(predicted)
    return pd.concat(out, ignore_index=True)


def observed_target(panel: pd.DataFrame, pred: pd.DataFrame) -> pd.DataFrame:
    obs = panel[["NUTS_ID", "year", "radiance_mean_corrected"]].rename(
        columns={
            "year": "target_year",
            "radiance_mean_corrected": "observed_radiance_target_nw",
        }
    )
    out = pred.merge(
        obs, on=["NUTS_ID", "target_year"], how="left", validate="many_to_one"
    )
    if out["observed_radiance_target_nw"].isna().any():
        raise RuntimeError("Predicted target lacks observed same-region target.")
    out["observed_change_nw"] = (
        out["observed_radiance_target_nw"] - out["radiance_anchor_nw"]
    )
    out["observed_change_percent"] = 100.0 * (
        out["observed_radiance_target_nw"] / out["radiance_anchor_nw"] - 1.0
    )
    out["observed_log10_change"] = np.log10(
        out["observed_radiance_target_nw"] / out["radiance_anchor_nw"]
    )
    for metric, suffix in [
        ("log10_change", "log10"),
        ("change_percent", "percentage_points"),
        ("change_nw", "nw"),
    ]:
        out[f"residual_{suffix}"] = out[f"predicted_{metric}"] - out[f"observed_{metric}"]
    return out


def score(group: pd.DataFrame) -> dict:
    metrics = {"n_nuts3": int(group["NUTS_ID"].nunique()), "n": len(group)}
    for suffix in ("log10", "percentage_points", "nw"):
        e = group[f"residual_{suffix}"].to_numpy(float)
        metrics.update({
            f"rmse_{suffix}": float(np.sqrt(np.mean(e ** 2))),
            f"mae_{suffix}": float(np.mean(np.abs(e))),
            f"bias_{suffix}": float(np.mean(e)),
        })
    return metrics


def score_groups(pred: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    rows = []
    for key, group in pred.groupby(keys, observed=True):
        if not isinstance(key, tuple):
            key = (key,)
        rows.append({**dict(zip(keys, key)), **score(group)})
    return pd.DataFrame(rows)


def temporal_holdouts(panel: pd.DataFrame, cfg: dict):
    """Two temporally separated folds; never train on a holdout year.

    Train 2014-2018, anchor 2018, predict 2019-20 (pre-snow test).
    Train 2014-2020, anchor 2020, predict 2022-24 (multi-year test).
    2021 is deliberately not a tested target in the principal comparison.
    """
    fits, predictions = [], []
    for cutoff, target_years in [(2018, [2019, 2020]), (2020, [2022, 2023, 2024])]:
        train = panel.loc[panel["year"].le(cutoff)].copy()
        models, fit_table = fit_models(train)
        fit_table.insert(0, "train_through", cutoff)
        fits.append(fit_table)
        anchor = panel.loc[panel["year"].eq(cutoff)].copy()
        for name, model in models.items():
            p = observed_target(panel, predict_changes(model, anchor, target_years))
            p["model"] = name
            p["train_through"] = cutoff
            p["horizon_years"] = p["target_year"] - cutoff
            predictions.append(p)
    pred = pd.concat(predictions, ignore_index=True)
    return (
        pd.concat(fits, ignore_index=True),
        pred,
        score_groups(pred, ["model", "train_through", "target_year", "horizon_years"]),
        score_groups(pred, ["model"]),
        score_groups(pred, ["model", "density_band"]),
    )


def full_endpoint_predictions(panel: pd.DataFrame, models: dict):
    a = panel.loc[panel["year"].eq(2014)].copy()
    rows = []
    for name, model in models.items():
        p = observed_target(panel, predict_changes(model, a, [2024]))
        p["model"] = name
        rows.append(p)
    all_pred = pd.concat(rows, ignore_index=True)
    scores = score_groups(all_pred, ["model", "density_band"])
    aggregate = []
    for (model, band), g in all_pred.groupby(["model", "density_band"], observed=True):
        aggregate.append({
            "model": model, "density_band": str(band), "n_nuts3": len(g),
            "observed_median_change_percent": float(g["observed_change_percent"].median()),
            "observed_median_change_nw": float(g["observed_change_nw"].median()),
            "observed_geometric_mean_change_percent": float(
                100.0 * (10.0 ** g["observed_log10_change"].mean() - 1.0)
            ),
            "predicted_geometric_mean_change_percent": float(
                100.0 * (10.0 ** g["predicted_log10_change"].mean() - 1.0)
            ),
            "predicted_median_change_nw": float(g["predicted_change_nw"].median()),
        })
    return all_pred, scores, pd.DataFrame(aggregate)


def descriptive_h_groups(panel: pd.DataFrame, full_predictions: pd.DataFrame):
    """Within-density habitability thirds, for description only, not fitting."""
    endpoint = full_predictions.loc[
        full_predictions["model"].eq("D_density"),
        ["NUTS_ID", "density_band", "H_core", "observed_change_percent",
         "observed_change_nw"],
    ].copy()
    summaries = []
    assigned = []
    for band, g in endpoint.groupby("density_band", observed=True):
        q1, q2 = g["H_core"].quantile([1 / 3, 2 / 3]).to_numpy(float)
        if np.isclose(q1, q2):
            # Avoid manufacturing distinctions among tied H values.
            g["H_tertile"] = "unresolved_tied_H"
        else:
            g["H_tertile"] = pd.cut(
                g["H_core"], bins=[-np.inf, q1, q2, np.inf],
                labels=["low_H", "medium_H", "high_H"],
            ).astype(str)
        assigned.append(g)
        for label, x in g.groupby("H_tertile", observed=True):
            summaries.append({
                "density_band": str(band), "H_group": str(label),
                "n_nuts3": len(x),
                "median_H_core": float(x["H_core"].median()),
                "median_observed_percent": float(x["observed_change_percent"].median()),
                "median_observed_nw": float(x["observed_change_nw"].median()),
            })
    return pd.concat(assigned, ignore_index=True), pd.DataFrame(summaries)


def sensitivity_models(panel: pd.DataFrame):
    """2021 omission and country-specific trends are sensitivity analyses."""
    pairs = []
    no21 = panel.loc[~panel["year"].eq(2021)].copy()
    fit_no21, tab_no21 = fit_models(
        no21, {k: FORMULAS[k] for k in ("D_density", "DH_density_trend")}
    )
    tab_no21.insert(0, "sensitivity", "omit_2021_no_refitting_correction")
    pairs.append(tab_no21)
    country_formula = {
        "D_country_trend": FORMULAS["D_density"] + " + year_c:C(country)",
        "DH_country_trend": FORMULAS["DH_density_trend"] + " + year_c:C(country)",
    }
    fit_ct, tab_ct = fit_models(panel, country_formula)
    tab_ct.insert(0, "sensitivity", "country_specific_temporal_trends")
    pairs.append(tab_ct)
    lr, df, p = lrt(fit_ct["D_country_trend"], fit_ct["DH_country_trend"])
    test = pd.DataFrame([{
        "comparison": "D_country_trend to DH_country_trend",
        "lr": lr, "df": df, "p_value": p,
        "bic_improvement": (
            fit_ct["D_country_trend"].bic - fit_ct["DH_country_trend"].bic
        ),
        "note": "Joint spatial and temporal H block after country-specific trends",
    }])
    return pd.concat(pairs, ignore_index=True), test


def plots(aggregate: pd.DataFrame, holdout: pd.DataFrame,
          h_summary: pd.DataFrame, labels: list[str], figures: Path, dpi: int):
    models = list(aggregate["model"].unique())
    fig, ax = plt.subplots(figsize=(10.6, 5.6))
    x = np.arange(len(labels))
    observed = aggregate[aggregate["model"].eq(models[0])].set_index("density_band")
    ax.plot(x, observed.loc[labels, "observed_geometric_mean_change_percent"],
            "ko-", linewidth=2.0, label="Observed 2014–2024 (geometric mean)")
    for name in models:
        g = aggregate[aggregate["model"].eq(name)].set_index("density_band")
        ax.plot(
            x, g.loc[labels, "predicted_geometric_mean_change_percent"],
            "--o", linewidth=1.3, label=name,
        )
    ax.axhline(0, color=".3", linewidth=.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("2014–2024 geometric-mean proportional change (%)")
    ax.set_xlabel("Fixed 2014–2020 population density (people km$^{-2}$)")
    ax.set_title("Pre-snow history: actual versus fitted 2014–2024 changes")
    ax.grid(axis="y", alpha=.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures / "observed_vs_fitted_density_endpoints.pdf", dpi=dpi)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, metric, ylabel in [
        (axes[0], "rmse_percentage_points", "Held-out RMSE (percentage points)"),
        (axes[1], "rmse_nw", "Held-out RMSE (nW cm$^{-2}$ sr$^{-1}$)"),
    ]:
        g = holdout.set_index("model").loc[models]
        ax.bar(np.arange(len(models)), g[metric].to_numpy(float), color=".5")
        ax.set_xticks(np.arange(len(models)), models, rotation=20, ha="right")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=.2)
    fig.suptitle("True temporal holdouts: pooled 2019–20 and 2022–24")
    fig.tight_layout()
    fig.savefig(figures / "temporal_holdout_predictive_error.pdf", dpi=dpi)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4))
    groups = ["low_H", "medium_H", "high_H"]
    for ax, metric, ylabel in [
        (axes[0], "median_observed_percent", "Median 2014–2024 change (%)"),
        (axes[1], "median_observed_nw", "Median 2014–2024 change (nW cm$^{-2}$ sr$^{-1}$)"),
    ]:
        for group in groups:
            g = h_summary.loc[h_summary["H_group"].eq(group)].set_index("density_band")
            vals = [float(g.loc[b, metric]) if b in g.index else np.nan for b in labels]
            ax.plot(np.arange(len(labels)), vals, "o-", label=group)
        ax.axhline(0, color=".3", linewidth=.8)
        ax.set_xticks(np.arange(len(labels)), labels)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=.2)
    axes[0].legend(title="Within-band H tertile", fontsize=9)
    fig.suptitle("Observed endpoints: geography H within fixed density class")
    fig.tight_layout()
    fig.savefig(figures / "observed_density_by_habitability.pdf", dpi=dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    cfg = config["analysis"]["presnow_density_h"]
    labels = list(config["analysis"]["density_temporal_change"]["density_band_labels"])
    processed = configured_path(config, "processed") / STAGE
    figures = configured_path(config, "figures") / "population_relationship" / STAGE
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    dpi = int(cfg["figure_dpi"])

    base, panel = load_panel(config, cfg)
    print(
        f"Pre-snow calibrated 2014-2024: {base['NUTS_ID'].nunique()} "
        f"locked NUTS3, H_core matched {panel['NUTS_ID'].nunique()}",
        flush=True,
    )
    models, model_fit = fit_models(panel)
    nested = nested_tests(models, ORDER)
    folds_fit, predictions, fold_scores, pooled_scores, density_scores = (
        temporal_holdouts(panel, cfg)
    )
    endpoints, endpoint_scores, class_aggregates = (
        full_endpoint_predictions(panel, models)
    )
    h_regions, h_summary = descriptive_h_groups(panel, endpoints)
    sensitivities, country_test = sensitivity_models(panel)

    files = {
        "candidate_fit": model_fit,
        "nested_tests": nested,
        "temporal_holdout_fit": folds_fit,
        "temporal_holdout_predictions": predictions,
        "temporal_holdout_by_fold": fold_scores,
        "temporal_holdout_overall": pooled_scores,
        "temporal_holdout_by_density": density_scores,
        "full_sample_endpoint_predictions": endpoints,
        "full_sample_endpoint_errors_by_density": endpoint_scores,
        "observed_vs_fitted_density_changes": class_aggregates,
        "observed_endpoint_h_region_assignment": h_regions,
        "observed_endpoint_h_summary": h_summary,
        "sensitivity_fit": sensitivities,
        "country_trend_sensitivity_test": country_test,
    }
    for stem, df in files.items():
        df.to_csv(processed / (stem + ".csv"), index=False)

    plots(class_aggregates, pooled_scores, h_summary, labels, figures, dpi)
    summary = processed / "presnow_density_h_summary.txt"
    with summary.open("w", encoding="utf-8") as f:
        f.write("PRESNOW CORRECTED VIIRS: DENSITY + TERRAIN HABITABILITY, 2014-2024\n")
        f.write("=" * 79 + "\n")
        f.write(
            "Source: fixed_nuts2016_viirs_2014_2024.parquet; original "
            "production pixel-calibrated radiance, without the later snow/surface "
            "correction for 2021.\n"
            "Population density: unchanged 2014-2020 NUTS2016 median classes.\n"
            "H predictor: NUTS2016 habitable_fraction_core, separate from density.\n"
            "No automatic winner selected; results are a pre-specified diagnostic.\n"
            "All comparisons involving H use the exact same H-matched regions.\n"
        )
        f.write(
            f"Original locked cohort={base['NUTS_ID'].nunique()}; "
            f"H-matched cohort={panel['NUTS_ID'].nunique()}; "
            f"missing H={base['NUTS_ID'].nunique()-panel['NUTS_ID'].nunique()}; "
            f"H_center={panel.attrs['h_median']:.6g}\n"
        )
        f.write(
            "Predictions are FE slope changes from the observed region-specific "
            "anchor; physical nW forecasts use that observed anchor, rather than "
            "only a percentage-score objective.\n"
            "Holdouts: train 2014-2018 -> 2019/20 (anchored 2018); "
            "train 2014-2020 -> 2022/23/24 (anchored 2020). "
            "2021 is not scored as an out-of-sample test year.\n"
            "Sensitivities: omit 2021 entirely; separately permit country "
            "specific linear time trends. These do not modify the primary data.\n"
            "Caution: full-period fit and endpoints are IN-SAMPLE reconstruction, "
            "not forecast verification. No future population-change predictor "
            "is yet implemented, because comparable post-2020 NUTS2016 "
            "population-density series are not established.\n\n"
        )
        for label, df in [
            ("CANDIDATE FIT", model_fit),
            ("NESTED LIKELIHOOD/BIC", nested),
            ("TEMPORAL HOLDOUT OVERALL", pooled_scores),
            ("TEMPORAL HOLDOUT BY FOLD", fold_scores),
            ("FULL-SAMPLE OBSERVED VS PREDICTED ENDPOINT BY DENSITY", class_aggregates),
            ("SENSITIVITY FIT", sensitivities),
            ("COUNTRY-TREND H TEST", country_test),
        ]:
            f.write("\n" + label + "\n" + df.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            ) + "\n")

    print("\nNESTED MODEL TESTS\n" + nested.to_string(index=False))
    print("\nHELD-OUT PREDICTION\n" + pooled_scores.to_string(index=False))
    print("\nWROTE\n" + str(summary))
    for name in files:
        print(processed / f"{name}.csv")
    for name in [
        "observed_vs_fitted_density_endpoints",
        "temporal_holdout_predictive_error",
        "observed_density_by_habitability",
    ]:
        print(figures / f"{name}.pdf")


if __name__ == "__main__":
    main()
