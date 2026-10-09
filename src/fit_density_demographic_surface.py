#!/usr/bin/env python3
"""Fit a historical population-density x density-change surface for ALAN change.

The response is annualised change in log10 corrected VIIRS radiance over fixed
five-year intervals. Predictors are:

  - baseline log10 population density; and
  - annualised percentage change in population density over the same interval.

To avoid NUTS boundary changes in this first mechanistic test, the default
windows are 2013-2018, 2014-2019 and 2015-2020, all on NUTS 2016 geography.

The central diagnostic is the model-predicted radiance trend when population
density change is exactly zero. This distinguishes a background rural-lighting
trend from change associated with demographic redistribution.

Because the three five-year windows overlap, inference uses NUTS_ID-clustered
robust covariance. Country and start-window fixed effects adjust level
differences; the fitted density x demographic-change surface is common across
Europe in this first pass.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy.stats import chi2
import statsmodels.formula.api as smf

from project_config import configured_path, load_config


def build_pairs(config: dict) -> pd.DataFrame:
    cfg = config["analysis"]["density_demographic_surface"]
    processed = configured_path(config, "processed")
    data = pd.read_parquet(processed / config["outputs"]["merged"]).copy()

    needed = [
        "NUTS_ID",
        "CNTR_CODE",
        "year",
        "population_density",
        "radiance_mean_corrected",
    ]
    missing = [c for c in needed if c not in data.columns]
    if missing:
        raise KeyError(f"Missing required merged columns: {missing}")

    data = data[needed].copy()
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
    data = data.loc[valid].copy()

    dup = data.duplicated(["NUTS_ID", "year"], keep=False)
    if dup.any():
        examples = data.loc[dup, ["NUTS_ID", "year"]].head(10)
        raise ValueError(
            "Merged table has duplicate NUTS_ID-year rows; cannot form unique "
            f"change pairs. Examples:\n{examples.to_string(index=False)}"
        )

    horizon = int(cfg["horizon_years"])
    rows = []
    index = data.set_index(["NUTS_ID", "year"]).sort_index()

    for start_year in cfg["start_years"]:
        start_year = int(start_year)
        end_year = start_year + horizon

        s = data[data["year"] == start_year].copy()
        e = data[data["year"] == end_year].copy()

        pair = s.merge(
            e,
            on="NUTS_ID",
            how="inner",
            suffixes=("_start", "_end"),
            validate="one_to_one",
        )
        pair = pair[
            pair["CNTR_CODE_start"] == pair["CNTR_CODE_end"]
        ].copy()

        pair["start_year"] = start_year
        pair["end_year"] = end_year
        pair["CNTR_CODE"] = pair["CNTR_CODE_start"]

        pair["baseline_log10_population_density"] = np.log10(
            pair["population_density_start"]
        )
        pair["annual_log10_population_density_change"] = (
            np.log10(pair["population_density_end"])
            - np.log10(pair["population_density_start"])
        ) / horizon
        pair["population_density_change_pct_per_year"] = 100.0 * (
            (pair["population_density_end"] / pair["population_density_start"])
            ** (1.0 / horizon)
            - 1.0
        )

        pair["annual_log10_radiance_change"] = (
            np.log10(pair["radiance_mean_corrected_end"])
            - np.log10(pair["radiance_mean_corrected_start"])
        ) / horizon
        pair["radiance_change_pct_per_year"] = 100.0 * (
            (pair["radiance_mean_corrected_end"] / pair["radiance_mean_corrected_start"])
            ** (1.0 / horizon)
            - 1.0
        )

        rows.append(pair)

    if not rows:
        raise RuntimeError("No five-year pairs were constructed.")

    out = pd.concat(rows, ignore_index=True)

    # Restrict country fixed effects to countries represented by enough NUTS3
    # regions in the paired data.
    counts = out.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    keep = counts[counts >= int(cfg["min_country_regions"])].index
    out = out[out["CNTR_CODE"].isin(keep)].copy()

    # Transparent light trimming of extreme demographic changes to prevent a
    # handful of likely revisions/outliers dominating a polynomial surface.
    qlo, qhi = [float(x) for x in cfg["population_change_trim_quantiles"]]
    lo = float(out["population_density_change_pct_per_year"].quantile(qlo))
    hi = float(out["population_density_change_pct_per_year"].quantile(qhi))
    out["included_in_surface_fit"] = out[
        "population_density_change_pct_per_year"
    ].between(lo, hi, inclusive="both")

    fit = out[out["included_in_surface_fit"]].copy()
    x_center = float(fit["baseline_log10_population_density"].median())
    fit["log_density_c"] = (
        fit["baseline_log10_population_density"] - x_center
    )
    fit["country"] = fit["CNTR_CODE"].astype("category")
    fit["window"] = fit["start_year"].astype(str).astype("category")

    # Preserve modelling metadata through attrs only for this in-memory object.
    fit.attrs["log_density_center"] = x_center
    fit.attrs["population_change_trim_low"] = lo
    fit.attrs["population_change_trim_high"] = hi
    return out, fit


def fit_models(fit: pd.DataFrame):
    density_only_formula = (
        "annual_log10_radiance_change ~ "
        "C(country) + C(window) + "
        "log_density_c + I(log_density_c ** 2) + I(log_density_c ** 3)"
    )

    surface_formula = (
        "annual_log10_radiance_change ~ "
        "C(country) + C(window) + "
        "log_density_c + I(log_density_c ** 2) + I(log_density_c ** 3) + "
        "population_density_change_pct_per_year + "
        "I(population_density_change_pct_per_year ** 2) + "
        "log_density_c:population_density_change_pct_per_year + "
        "I(log_density_c ** 2):population_density_change_pct_per_year"
    )

    density_only = smf.ols(density_only_formula, data=fit).fit(
        cov_type="cluster",
        cov_kwds={"groups": fit["NUTS_ID"]},
    )
    surface = smf.ols(surface_formula, data=fit).fit(
        cov_type="cluster",
        cov_kwds={"groups": fit["NUTS_ID"]},
    )

    # Cluster-robust joint Wald test of the demographic-change block.
    block_names = [
        name for name in surface.params.index
        if "population_density_change_pct_per_year" in name
    ]
    indices = [surface.params.index.get_loc(name) for name in block_names]
    beta = surface.params.iloc[indices].to_numpy(float)
    cov = surface.cov_params().iloc[indices, indices].to_numpy(float)
    stat = float(beta.T @ np.linalg.pinv(cov) @ beta)
    df = len(indices)
    p = float(chi2.sf(stat, df))

    test = {
        "test": "population-density-change block",
        "wald_chi2": stat,
        "df": df,
        "p": p,
        "terms": " | ".join(block_names),
        "density_only_r2": float(density_only.rsquared),
        "surface_r2": float(surface.rsquared),
        "delta_r2": float(surface.rsquared - density_only.rsquared),
        "density_only_aic": float(density_only.aic),
        "surface_aic": float(surface.aic),
        "density_only_bic": float(density_only.bic),
        "surface_bic": float(surface.bic),
    }
    return density_only, surface, test


def average_design_vector(result, fit, density, pop_change_pct):
    """Average fixed-effect design over the observed country/window mixture."""
    combos = (
        fit.groupby(["country", "window"], observed=True)
        .size()
        .rename("weight")
        .reset_index()
    )
    combos["weight"] = combos["weight"] / combos["weight"].sum()

    center = float(fit.attrs["log_density_center"])
    pred = combos[["country", "window"]].copy()
    pred["log_density_c"] = np.log10(float(density)) - center
    pred["population_density_change_pct_per_year"] = float(pop_change_pct)

    # Rebuild the RHS design matrix from the fitted formula.  Some recent
    # statsmodels/Python combinations expose model.data as PandasData without
    # a design_info attribute, so relying on result.model.data.design_info is
    # not portable.  Reindexing to the fitted parameter names also guarantees
    # identical coefficient order.
    rhs = result.model.formula.split("~", 1)[1]
    X_df = dmatrix(rhs, pred, return_type="dataframe")
    X_df = X_df.reindex(columns=result.params.index, fill_value=0.0)
    X = X_df.to_numpy(dtype=float)

    w = combos["weight"].to_numpy(float)
    return np.average(X, axis=0, weights=w)


def predict_point(result, fit, density, pop_change_pct):
    xbar = average_design_vector(result, fit, density, pop_change_pct)
    beta = result.params.to_numpy(float)
    cov = result.cov_params().to_numpy(float)
    est = float(xbar @ beta)
    se = float(np.sqrt(max(0.0, xbar @ cov @ xbar)))
    return {
        "population_density": float(density),
        "population_density_change_pct_per_year": float(pop_change_pct),
        "annual_log10_radiance_change": est,
        "annual_log10_radiance_change_se": se,
        "annual_log10_radiance_change_ci_low": est - 1.96 * se,
        "annual_log10_radiance_change_ci_high": est + 1.96 * se,
        "predicted_radiance_change_pct_per_year": 100.0 * (10.0**est - 1.0),
        "predicted_radiance_change_pct_5yr": 100.0 * (10.0 ** (5.0 * est) - 1.0),
        "predicted_radiance_change_pct_per_year_ci_low": (
            100.0 * (10.0 ** (est - 1.96 * se) - 1.0)
        ),
        "predicted_radiance_change_pct_per_year_ci_high": (
            100.0 * (10.0 ** (est + 1.96 * se) - 1.0)
        ),
    }


def make_surface_grid(result, fit, config):
    cfg = config["analysis"]["density_demographic_surface"]
    dlo, dhi = [float(x) for x in cfg["surface_density_quantiles"]]
    plo, phi = [float(x) for x in cfg["surface_population_change_quantiles"]]

    density_lo = float(
        fit["population_density_start"].quantile(dlo)
    )
    density_hi = float(
        fit["population_density_start"].quantile(dhi)
    )
    pop_lo = float(
        fit["population_density_change_pct_per_year"].quantile(plo)
    )
    pop_hi = float(
        fit["population_density_change_pct_per_year"].quantile(phi)
    )

    densities = np.geomspace(
        density_lo,
        density_hi,
        int(cfg["surface_density_points"]),
    )
    pop_changes = np.linspace(
        pop_lo,
        pop_hi,
        int(cfg["surface_population_change_points"]),
    )

    rows = []
    for dp in pop_changes:
        for density in densities:
            rows.append(predict_point(result, fit, density, dp))
    return pd.DataFrame(rows)


def make_reference_scenarios(result, fit, config):
    cfg = config["analysis"]["density_demographic_surface"]
    rows = []
    for dp in cfg["population_change_scenarios_pct_per_year"]:
        for density in cfg["reference_population_densities"]:
            rows.append(predict_point(result, fit, density, dp))
    return pd.DataFrame(rows)


def make_plots(surface_grid, scenarios, config, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)

    # 2-D historical response surface.
    pivot = surface_grid.pivot(
        index="population_density_change_pct_per_year",
        columns="population_density",
        values="predicted_radiance_change_pct_per_year",
    )
    x = pivot.columns.to_numpy(float)
    y = pivot.index.to_numpy(float)
    z = pivot.to_numpy(float)

    fig, ax = plt.subplots(figsize=(9.0, 6.3))
    levels = 18
    cs = ax.contourf(x, y, z, levels=levels)
    ax.contour(x, y, z, levels=10, linewidths=0.6)
    ax.axhline(0.0, linewidth=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("Baseline population density (people km⁻²)")
    ax.set_ylabel("Population-density change (% yr⁻¹)")
    ax.set_title("Historical surface predicting VIIRS radiance change")
    cb = fig.colorbar(cs, ax=ax)
    cb.set_label("Predicted corrected-VIIRS radiance change (% yr⁻¹)")
    fig.tight_layout()
    path = outdir / "radiance_change_density_demographic_surface.pdf"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Slices through the surface at selected demographic trajectories.
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    for dp in config["analysis"]["density_demographic_surface"][
        "population_change_scenarios_pct_per_year"
    ]:
        g = scenarios[
            np.isclose(
                scenarios["population_density_change_pct_per_year"],
                float(dp),
            )
        ].sort_values("population_density")
        ax.plot(
            g["population_density"],
            g["predicted_radiance_change_pct_per_year"],
            marker="o",
            label=f"{float(dp):+g}% yr⁻¹ density change",
        )
    ax.axhline(0.0, linewidth=0.9)
    ax.set_xscale("log")
    ax.set_xlabel("Baseline population density (people km⁻²)")
    ax.set_ylabel("Predicted corrected-VIIRS radiance change (% yr⁻¹)")
    ax.set_title("Demographic trajectories through the historical change surface")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    path = outdir / "radiance_change_surface_slices.pdf"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # The crucial zero-demographic-change slice.
    zero = scenarios[
        np.isclose(
            scenarios["population_density_change_pct_per_year"],
            0.0,
        )
    ].sort_values("population_density")

    fig, ax = plt.subplots(figsize=(8.6, 5.7))
    ax.plot(
        zero["population_density"],
        zero["predicted_radiance_change_pct_per_year"],
        marker="o",
    )
    ax.fill_between(
        zero["population_density"],
        zero["predicted_radiance_change_pct_per_year_ci_low"],
        zero["predicted_radiance_change_pct_per_year_ci_high"],
        alpha=0.2,
    )
    ax.axhline(0.0, linewidth=0.9)
    ax.set_xscale("log")
    ax.set_xlabel("Baseline population density (people km⁻²)")
    ax.set_ylabel("Predicted corrected-VIIRS radiance change (% yr⁻¹)")
    ax.set_title("Predicted radiance trend when population density is stable")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    path = outdir / "zero_population_density_change_slice.pdf"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["density_demographic_surface"]

    all_pairs, fit = build_pairs(config)
    density_only, surface, block_test = fit_models(fit)
    grid = make_surface_grid(surface, fit, config)
    scenarios = make_reference_scenarios(surface, fit, config)

    processed = configured_path(config, "processed")
    outdir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "density_demographic_surface"
    )

    pair_path = processed / "density_demographic_5yr_pairs.csv"
    coef_path = processed / "density_demographic_surface_coefficients.csv"
    test_path = processed / "density_demographic_surface_model_test.csv"
    grid_path = processed / "density_demographic_surface_grid.csv"
    scenario_path = processed / "density_demographic_surface_reference_scenarios.csv"
    zero_path = processed / "density_demographic_zero_change_slice.csv"
    summary_path = processed / "density_demographic_surface_summary.txt"

    all_pairs.to_csv(pair_path, index=False)

    coef = pd.DataFrame({
        "term": surface.params.index,
        "estimate": surface.params.values,
        "cluster_robust_se": surface.bse.values,
        "z": surface.tvalues.values,
        "p": surface.pvalues.values,
    })
    coef.to_csv(coef_path, index=False)
    pd.DataFrame([block_test]).to_csv(test_path, index=False)
    grid.to_csv(grid_path, index=False)
    scenarios.to_csv(scenario_path, index=False)
    zero = scenarios[
        np.isclose(
            scenarios["population_density_change_pct_per_year"],
            0.0,
        )
    ].copy()
    zero.to_csv(zero_path, index=False)

    make_plots(grid, scenarios, config, outdir)

    with summary_path.open("w", encoding="utf-8") as h:
        h.write("POPULATION-DENSITY × DEMOGRAPHIC-CHANGE SURFACE\n")
        h.write("=" * 72 + "\n\n")
        h.write(
            "Population change in this analysis means change in POPULATION "
            "DENSITY. On a fixed NUTS polygon, percentage population-density "
            "change is numerically identical to percentage population change.\n\n"
        )
        h.write(
            f"Five-year windows: "
            + ", ".join(
                f"{int(y)}-{int(y)+int(cfg['horizon_years'])}"
                for y in cfg["start_years"]
            )
            + "\n"
        )
        h.write(
            "These windows remain within NUTS 2016 geography, avoiding the "
            "2021/2024 NUTS boundary switches in this first test.\n"
        )
        h.write(
            f"Pairs available before demographic trimming: {len(all_pairs)}; "
            f"surface-fit rows: {len(fit)}; unique NUTS3: "
            f"{fit['NUTS_ID'].nunique()}; countries: {fit['CNTR_CODE'].nunique()}.\n"
        )
        h.write(
            f"Population-density change fit range after configured trimming: "
            f"{fit.attrs['population_change_trim_low']:.6g} to "
            f"{fit.attrs['population_change_trim_high']:.6g}% yr-1.\n\n"
        )

        h.write("DEMOGRAPHIC-CHANGE BLOCK TEST\n")
        for key, value in block_test.items():
            h.write(f"{key}: {value}\n")

        h.write("\nZERO POPULATION-DENSITY-CHANGE SLICE\n")
        h.write(
            zero[[
                "population_density",
                "predicted_radiance_change_pct_per_year",
                "predicted_radiance_change_pct_per_year_ci_low",
                "predicted_radiance_change_pct_per_year_ci_high",
                "predicted_radiance_change_pct_5yr",
            ]].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nINTERPRETATION\n")
        h.write(
            "- If the zero-density-change slice remains strongly positive at "
            "low baseline density, rural brightening is not simply a by-product "
            "of population growth.\n"
        )
        h.write(
            "- If predicted radiance growth rises strongly with positive "
            "population-density change, demographic redistribution provides a "
            "second forecasting axis beyond baseline density.\n"
        )
        h.write(
            "- This is an empirical historical response surface for predicting "
            "future percentage CHANGE, not absolute future VIIRS radiance.\n"
        )
        h.write(
            "- The first test deliberately uses only 2013-2020 NUTS 2016 "
            "geography. A production forecast should harmonise all years to a "
            "single geography before training on the full 2013-2024 record.\n"
        )

    print("\nDEMOGRAPHIC-CHANGE BLOCK TEST")
    print(pd.DataFrame([block_test]).to_string(
        index=False, float_format=lambda x: f"{x:.6g}"
    ))
    print("\nZERO POPULATION-DENSITY-CHANGE SLICE")
    print(
        zero[[
            "population_density",
            "predicted_radiance_change_pct_per_year",
            "predicted_radiance_change_pct_per_year_ci_low",
            "predicted_radiance_change_pct_per_year_ci_high",
            "predicted_radiance_change_pct_5yr",
        ]].to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print("\nWROTE")
    for p in [
        pair_path,
        coef_path,
        test_path,
        grid_path,
        scenario_path,
        zero_path,
        summary_path,
    ]:
        print(p)


if __name__ == "__main__":
    main()
