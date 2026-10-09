#!/usr/bin/env python3
"""Test terrain habitability H as an independent historical predictor.

H is kept separate from nominal population density.

Three analyses are produced:

1. Longitudinal mixed-model test:
   log10(corrected VIIRS radiance) ~ country x year
   + cubic nominal-density x year + H + H x year
   + (1 | NUTS_ID::NUTS_release)

2. Temporally separated headroom test on fixed NUTS2016 geography:
   2013-2014 define baseline radiance/density; 2015-2020 radiance slope is the
   response. H is tested after baseline density, baseline radiance,
   population-density trend and country.

3. Direct descriptive VIIRS check on matched NUTS2016 regions:
   actual 2013 and 2020 radiance values, absolute change, percent change and
   annualised percent change are compared among low/mid/high H tertiles within
   nominally sparse regions (<100 people km-2 at baseline). Year-by-year
   2013-2020 median radiance trajectories are also written and plotted.

The direct check is deliberately descriptive rather than a replacement for the
models: it lets the expected H pattern be seen explicitly in the VIIRS data.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2
import statsmodels.formula.api as smf

from project_config import configured_path, load_config


def drop_unused_categories(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    for col in out.columns:
        if isinstance(out[col].dtype, pd.CategoricalDtype):
            out[col] = out[col].cat.remove_unused_categories()
    return out


def fit_mixedlm(formula: str, data: pd.DataFrame):
    fit = drop_unused_categories(data)
    model = smf.mixedlm(
        formula,
        data=fit,
        groups=fit["region_key"],
        re_formula="1",
    )
    last = None
    for method in ("lbfgs", "bfgs", "cg"):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = model.fit(
                    reml=False,
                    method=method,
                    maxiter=4000,
                    disp=False,
                )
            if result.converged:
                return result, method
            last = RuntimeError(f"non-converged with {method}")
        except Exception as exc:
            last = exc
    raise RuntimeError(f"MixedLM failed for {formula}: {last}")


def lrt(reduced, full):
    lr = max(0.0, 2.0 * (full.llf - reduced.llf))
    df = int(len(full.fe_params) - len(reduced.fe_params))
    p = float(chi2.sf(lr, df)) if df > 0 else np.nan
    return lr, df, p


def prepare_data(config: dict) -> pd.DataFrame:
    cfg = config["analysis"]["habitability_independent_factor"]
    processed = configured_path(config, "processed")

    merged = pd.read_parquet(processed / config["outputs"]["merged"]).copy()
    habit = pd.read_parquet(
        processed / config["outputs"]["nuts3_habitability"]
    ).copy()

    hcols = [
        "nuts_release",
        "NUTS_ID",
        "habitable_fraction_strict",
        "habitable_fraction_core",
        "habitable_fraction_permissive",
    ]
    habit = habit[hcols].drop_duplicates(["nuts_release", "NUTS_ID"])

    d = merged.merge(
        habit,
        on=["nuts_release", "NUTS_ID"],
        how="inner",
        validate="many_to_one",
    )
    d = d[
        d["year"].between(
            int(cfg["first_year"]),
            int(cfg["last_year"]),
            inclusive="both",
        )
    ].copy()

    good = (
        np.isfinite(pd.to_numeric(d["radiance_mean_corrected"], errors="coerce"))
        & np.isfinite(pd.to_numeric(d["population_density"], errors="coerce"))
        & (d["radiance_mean_corrected"] > 0)
        & (d["population_density"] > 0)
        & d["NUTS_ID"].notna()
        & d["CNTR_CODE"].notna()
    )
    for scenario in ("strict", "core", "permissive"):
        good &= np.isfinite(d[f"habitable_fraction_{scenario}"])

    d = d.loc[good].copy()
    d["region_key"] = (
        d["NUTS_ID"].astype(str) + "::" + d["nuts_release"].astype(str)
    )

    repeated = d.groupby("region_key").size()
    d = d[
        d["region_key"].isin(
            repeated[repeated >= int(cfg["min_repeated_years"])].index
        )
    ].copy()

    counts = d.groupby("CNTR_CODE")["region_key"].nunique()
    keep = counts[counts >= int(cfg["min_country_regions"])].index
    d = d[d["CNTR_CODE"].isin(keep)].copy()

    d["log_rad"] = np.log10(d["radiance_mean_corrected"].astype(float))
    d["log_pop"] = np.log10(d["population_density"].astype(float))
    d["country"] = d["CNTR_CODE"].astype("category")
    return d


def fit_independent_h_temporal(
    data: pd.DataFrame,
    scenario: str,
    first_year: int,
    last_year: int,
    nuts_release: int | None,
):
    fit = data[
        data["year"].between(first_year, last_year, inclusive="both")
    ].copy()
    if nuts_release is not None:
        fit = fit[fit["nuts_release"] == int(nuts_release)].copy()

    repeated = fit.groupby("region_key").size()
    fit = fit[fit["region_key"].isin(repeated[repeated >= 2].index)].copy()

    fit["year_c"] = fit["year"].astype(float) - ((first_year + last_year) / 2.0)
    fit["x"] = fit["log_pop"] - float(fit["log_pop"].median())
    hcol = f"habitable_fraction_{scenario}"
    fit["h_c"] = fit[hcol] - float(fit[hcol].median())
    fit["country"] = fit["CNTR_CODE"].astype("category")

    reduced_formula = (
        "log_rad ~ C(country) * year_c + "
        "year_c * (x + I(x ** 2) + I(x ** 3)) + h_c"
    )
    full_formula = reduced_formula + " + year_c:h_c"

    reduced, opt0 = fit_mixedlm(reduced_formula, fit)
    full, opt1 = fit_mixedlm(full_formula, fit)
    lr, df, p = lrt(reduced, full)

    term = "year_c:h_c"
    beta = float(full.fe_params[term])
    se = float(np.sqrt(full.cov_params().loc[term, term]))

    delta_low_vs_high = beta * (0.1 - 1.0)
    low_vs_high_pct = 100.0 * (10.0 ** delta_low_vs_high - 1.0)

    return {
        "scenario": scenario,
        "interval": f"{first_year}-{last_year}",
        "nuts_release_filter": "all" if nuts_release is None else int(nuts_release),
        "n_observations": len(fit),
        "n_region_keys": fit["region_key"].nunique(),
        "n_countries": fit["CNTR_CODE"].nunique(),
        "H_median": float(fit[hcol].median()),
        "H_x_year_beta_log10_per_year_per_H": beta,
        "H_x_year_se": se,
        "lrt_lr": lr,
        "lrt_df": df,
        "lrt_p": p,
        "relative_annual_growth_lowH_0p1_vs_highH_1_pct": low_vs_high_pct,
        "reduced_aic": float(reduced.aic),
        "full_aic": float(full.aic),
        "reduced_bic": float(reduced.bic),
        "full_bic": float(full.bic),
        "reduced_optimizer": opt0,
        "full_optimizer": opt1,
    }


def linear_slope(group: pd.DataFrame, value_col: str):
    g = group[
        np.isfinite(group[value_col]) & np.isfinite(group["year"])
    ].sort_values("year")
    if len(g) < 2 or g["year"].nunique() < 2:
        return np.nan
    return float(
        np.polyfit(
            g["year"].to_numpy(float),
            g[value_col].to_numpy(float),
            1,
        )[0]
    )


def build_headroom_table(data: pd.DataFrame, config: dict):
    cfg = config["analysis"]["habitability_independent_factor"]
    release = int(cfg["fixed_geography_nuts_release"])
    baseline_years = [int(x) for x in cfg["headroom_baseline_years"]]
    first = int(cfg["headroom_outcome_first_year"])
    last = int(cfg["headroom_outcome_last_year"])

    d = data[data["nuts_release"] == release].copy()
    base = d[d["year"].isin(baseline_years)].copy()
    outcome = d[d["year"].between(first, last, inclusive="both")].copy()

    base_summary = (
        base.groupby(["NUTS_ID", "CNTR_CODE"], as_index=False)
        .agg(
            baseline_log_radiance=("log_rad", "mean"),
            baseline_log_population_density=("log_pop", "mean"),
            H_core=("habitable_fraction_core", "first"),
            n_baseline_years=("year", "nunique"),
        )
    )

    rows = []
    for (nuts_id, country), g in outcome.groupby(
        ["NUTS_ID", "CNTR_CODE"], sort=False
    ):
        rows.append(
            {
                "NUTS_ID": nuts_id,
                "CNTR_CODE": country,
                "radiance_log_slope_outcome": linear_slope(g, "log_rad"),
                "population_log_slope_outcome": linear_slope(g, "log_pop"),
                "n_outcome_years": int(g["year"].nunique()),
            }
        )
    slopes = pd.DataFrame(rows)

    out = base_summary.merge(
        slopes,
        on=["NUTS_ID", "CNTR_CODE"],
        how="inner",
        validate="one_to_one",
    )
    out = out[
        (out["n_baseline_years"] >= len(baseline_years))
        & (out["n_outcome_years"] >= int(cfg["headroom_min_outcome_years"]))
        & np.isfinite(out["radiance_log_slope_outcome"])
        & np.isfinite(out["population_log_slope_outcome"])
    ].copy()

    counts = out.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    out = out[
        out["CNTR_CODE"].isin(
            counts[counts >= int(cfg["min_country_regions"])].index
        )
    ].copy()

    out["country"] = out["CNTR_CODE"].astype("category")
    out["x"] = (
        out["baseline_log_population_density"]
        - float(out["baseline_log_population_density"].median())
    )
    out["baseline_rad_c"] = (
        out["baseline_log_radiance"]
        - float(out["baseline_log_radiance"].median())
    )
    out["pop_slope_c"] = (
        out["population_log_slope_outcome"]
        - float(out["population_log_slope_outcome"].median())
    )
    out["h_c"] = out["H_core"] - float(out["H_core"].median())
    return out


def fit_headroom_models(table: pd.DataFrame):
    reduced_formula = (
        "radiance_log_slope_outcome ~ C(country) + "
        "x + I(x ** 2) + I(x ** 3) + "
        "baseline_rad_c + pop_slope_c"
    )
    h_formula = reduced_formula + " + h_c"
    interaction_formula = h_formula + " + x:h_c"

    d = drop_unused_categories(table)
    reduced = smf.ols(reduced_formula, data=d).fit(cov_type="HC3")
    h_model = smf.ols(h_formula, data=d).fit(cov_type="HC3")
    interaction = smf.ols(interaction_formula, data=d).fit(cov_type="HC3")

    terms = ["h_c", "x:h_c"]
    idx = [interaction.params.index.get_loc(t) for t in terms]
    beta = interaction.params.iloc[idx].to_numpy(float)
    cov = interaction.cov_params().iloc[idx, idx].to_numpy(float)
    stat = float(beta.T @ np.linalg.pinv(cov) @ beta)
    block_p = float(chi2.sf(stat, len(terms)))

    rows = []
    for name, model in [
        ("headroom_reduced", reduced),
        ("headroom_plus_H", h_model),
        ("headroom_plus_H_x_density", interaction),
    ]:
        rows.append(
            {
                "model": name,
                "n": len(table),
                "r2": float(model.rsquared),
                "adj_r2": float(model.rsquared_adj),
                "aic": float(model.aic),
                "bic": float(model.bic),
                "H_beta": float(model.params["h_c"]) if "h_c" in model.params else np.nan,
                "H_se": float(model.bse["h_c"]) if "h_c" in model.params else np.nan,
                "H_p": float(model.pvalues["h_c"]) if "h_c" in model.params else np.nan,
            }
        )

    block = {
        "test": (
            "H_core block after baseline density, baseline radiance, "
            "population-density trend, and country"
        ),
        "wald_chi2": stat,
        "df": len(terms),
        "p": block_p,
        "terms": " + ".join(terms),
        "delta_r2_H_only": float(h_model.rsquared - reduced.rsquared),
        "delta_r2_H_plus_interaction": float(
            interaction.rsquared - reduced.rsquared
        ),
        "H_main_beta": float(h_model.params["h_c"]),
        "H_main_se": float(h_model.bse["h_c"]),
        "H_main_p": float(h_model.pvalues["h_c"]),
    }
    return pd.DataFrame(rows), block


def build_direct_viirs_check(data: pd.DataFrame, config: dict):
    cfg = config["analysis"]["habitability_independent_factor"]
    release = int(cfg["fixed_geography_nuts_release"])
    y0 = int(cfg["direct_start_year"])
    y1 = int(cfg["direct_end_year"])
    cutoff = float(cfg["direct_low_density_cutoff"])

    d = data[
        (data["nuts_release"] == release)
        & data["year"].between(y0, y1, inclusive="both")
    ].copy()

    start_cols = [
        "NUTS_ID",
        "CNTR_CODE",
        "population_density",
        "radiance_mean_corrected",
        "habitable_fraction_core",
    ]
    if "radiance_mean_raw" in d.columns:
        start_cols.append("radiance_mean_raw")

    start = d[d["year"] == y0][start_cols].copy()
    start = start.rename(
        columns={
            "population_density": "population_density_start",
            "radiance_mean_corrected": "radiance_corrected_start",
            "radiance_mean_raw": "radiance_raw_start",
            "habitable_fraction_core": "H_core",
        }
    )

    end_cols = ["NUTS_ID", "radiance_mean_corrected"]
    if "radiance_mean_raw" in d.columns:
        end_cols.append("radiance_mean_raw")
    end = d[d["year"] == y1][end_cols].copy().rename(
        columns={
            "radiance_mean_corrected": "radiance_corrected_end",
            "radiance_mean_raw": "radiance_raw_end",
        }
    )

    pairs = start.merge(
        end,
        on="NUTS_ID",
        how="inner",
        validate="one_to_one",
    )
    pairs = pairs[pairs["population_density_start"] < cutoff].copy()

    q1 = float(pairs["H_core"].quantile(1.0 / 3.0))
    q2 = float(pairs["H_core"].quantile(2.0 / 3.0))
    pairs["H_group"] = pd.cut(
        pairs["H_core"],
        bins=[-np.inf, q1, q2, np.inf],
        labels=["low_H", "mid_H", "high_H"],
        ordered=True,
    )

    years = y1 - y0
    pairs["delta_radiance_corrected"] = (
        pairs["radiance_corrected_end"] - pairs["radiance_corrected_start"]
    )
    pairs["pct_change_corrected"] = 100.0 * (
        pairs["radiance_corrected_end"]
        / pairs["radiance_corrected_start"]
        - 1.0
    )
    pairs["annualised_pct_change_corrected"] = 100.0 * (
        (
            pairs["radiance_corrected_end"]
            / pairs["radiance_corrected_start"]
        ) ** (1.0 / years)
        - 1.0
    )

    if {"radiance_raw_start", "radiance_raw_end"}.issubset(pairs.columns):
        raw_good = (
            np.isfinite(pairs["radiance_raw_start"])
            & np.isfinite(pairs["radiance_raw_end"])
            & (pairs["radiance_raw_start"] > 0)
            & (pairs["radiance_raw_end"] > 0)
        )
        pairs["annualised_pct_change_raw"] = np.nan
        pairs.loc[raw_good, "annualised_pct_change_raw"] = 100.0 * (
            (
                pairs.loc[raw_good, "radiance_raw_end"]
                / pairs.loc[raw_good, "radiance_raw_start"]
            ) ** (1.0 / years)
            - 1.0
        )

    calibration_countries = set(config["calibration"]["countries"])
    pairs["is_calibration_country"] = pairs["CNTR_CODE"].isin(
        calibration_countries
    )

    summaries = []
    for group, g in pairs.groupby("H_group", observed=True):
        row = {
            "H_group": str(group),
            "n_regions": int(len(g)),
            "H_min": float(g["H_core"].min()),
            "H_median": float(g["H_core"].median()),
            "H_max": float(g["H_core"].max()),
            "median_population_density_start": float(
                g["population_density_start"].median()
            ),
            "median_radiance_corrected_start": float(
                g["radiance_corrected_start"].median()
            ),
            "median_radiance_corrected_end": float(
                g["radiance_corrected_end"].median()
            ),
            "median_absolute_change_corrected": float(
                g["delta_radiance_corrected"].median()
            ),
            "median_pct_change_corrected": float(
                g["pct_change_corrected"].median()
            ),
            "median_annualised_pct_change_corrected": float(
                g["annualised_pct_change_corrected"].median()
            ),
            "mean_annualised_pct_change_corrected": float(
                g["annualised_pct_change_corrected"].mean()
            ),
        }
        if "annualised_pct_change_raw" in g.columns:
            row["median_annualised_pct_change_raw"] = float(
                g["annualised_pct_change_raw"].median()
            )
        summaries.append(row)
    summary = pd.DataFrame(summaries)

    noncal_rows = []
    noncal = pairs[~pairs["is_calibration_country"]].copy()
    for group, g in noncal.groupby("H_group", observed=True):
        noncal_rows.append(
            {
                "H_group": str(group),
                "n_regions_noncalibration_countries": int(len(g)),
                "median_annualised_pct_change_corrected_noncalibration": float(
                    g["annualised_pct_change_corrected"].median()
                ),
            }
        )
    noncal_summary = pd.DataFrame(noncal_rows)
    summary = summary.merge(noncal_summary, on="H_group", how="left")

    group_map = pairs[["NUTS_ID", "H_group"]].copy()
    traj = d.merge(group_map, on="NUTS_ID", how="inner", validate="many_to_one")
    trajectory = (
        traj.groupby(["year", "H_group"], observed=True)
        .agg(
            n_regions=("NUTS_ID", "nunique"),
            median_corrected_radiance=("radiance_mean_corrected", "median"),
            mean_corrected_radiance=("radiance_mean_corrected", "mean"),
        )
        .reset_index()
    )

    return pairs, summary, trajectory, q1, q2


def make_direct_plots(pairs, trajectory, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    for group in ["low_H", "mid_H", "high_H"]:
        g = trajectory[
            trajectory["H_group"].astype(str) == group
        ].sort_values("year")
        ax.plot(
            g["year"],
            g["median_corrected_radiance"],
            marker="o",
            label=group.replace("_", " "),
        )
    ax.set_xlabel("Year")
    ax.set_ylabel("Median corrected VIIRS radiance (nW cm⁻² sr⁻¹)")
    ax.set_title("Direct VIIRS trajectories in nominally sparse NUTS3 regions")
    ax.grid(alpha=0.2)
    ax.legend(title="Terrain habitability")
    fig.tight_layout()
    p1 = outdir / "direct_viirs_radiance_by_habitability_2013_2020.pdf"
    fig.savefig(p1, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    groups = ["low_H", "mid_H", "high_H"]
    values = [
        pairs.loc[
            pairs["H_group"].astype(str) == group,
            "annualised_pct_change_corrected",
        ].dropna().to_numpy(float)
        for group in groups
    ]
    ax.boxplot(
        values,
        tick_labels=["Low H", "Mid H", "High H"],
        showfliers=False,
    )
    ax.axhline(0.0, linewidth=0.9)
    ax.set_ylabel("Annualised corrected-VIIRS change (% yr⁻¹)")
    ax.set_xlabel("Terrain habitability tertile")
    ax.set_title("Direct 2013–2020 VIIRS change in sparse NUTS3 regions")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    p2 = outdir / "direct_viirs_change_by_habitability_2013_2020.pdf"
    fig.savefig(p2, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return p1, p2


def partial_effect_plot(results: pd.DataFrame, output: Path):
    fig, ax = plt.subplots(figsize=(8.4, 5.7))
    H = np.linspace(0.1, 1.0, 100)

    for interval in ["2013-2024", "2013-2020"]:
        row = results[
            (results["scenario"] == "core")
            & (results["interval"] == interval)
        ]
        if row.empty:
            continue
        row = row.iloc[0]
        beta = float(row["H_x_year_beta_log10_per_year_per_H"])
        h_med = float(row["H_median"])
        delta = beta * (H - h_med)
        pct = 100.0 * (10.0**delta - 1.0)
        ax.plot(H, pct, label=f"{interval}")

    ax.axhline(0.0, linewidth=0.9)
    ax.set_xlabel("Terrain habitable fraction H")
    ax.set_ylabel(
        "H-dependent annual radiance-growth difference\n"
        "relative to median H (%)"
    )
    ax.set_title("Independent temporal effect of terrain habitability")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["habitability_independent_factor"]
    data = prepare_data(config)

    temporal_rows = []
    for scenario in ("strict", "core", "permissive"):
        temporal_rows.append(
            fit_independent_h_temporal(
                data,
                scenario,
                int(cfg["first_year"]),
                int(cfg["last_year"]),
                nuts_release=None,
            )
        )
        temporal_rows.append(
            fit_independent_h_temporal(
                data,
                scenario,
                int(cfg["fixed_geography_first_year"]),
                int(cfg["fixed_geography_last_year"]),
                nuts_release=int(cfg["fixed_geography_nuts_release"]),
            )
        )
    temporal = pd.DataFrame(temporal_rows)

    headroom_table = build_headroom_table(data, config)
    headroom_models, headroom_test = fit_headroom_models(headroom_table)

    pairs, direct_summary, trajectory, q1, q2 = build_direct_viirs_check(
        data, config
    )

    processed = configured_path(config, "processed")
    outdir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "habitability_independent_factor"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    temporal_path = processed / "habitability_independent_temporal_tests.csv"
    headroom_table_path = processed / "habitability_headroom_region_table.csv"
    headroom_models_path = processed / "habitability_headroom_model_comparison.csv"
    headroom_test_path = processed / "habitability_headroom_H_test.csv"
    direct_pairs_path = processed / "habitability_direct_viirs_endpoint_pairs.csv"
    direct_summary_path = processed / "habitability_direct_viirs_summary.csv"
    trajectory_path = processed / "habitability_direct_viirs_trajectory.csv"
    summary_path = processed / "habitability_independent_factor_summary.txt"

    temporal.to_csv(temporal_path, index=False)
    headroom_table.to_csv(headroom_table_path, index=False)
    headroom_models.to_csv(headroom_models_path, index=False)
    pd.DataFrame([headroom_test]).to_csv(headroom_test_path, index=False)
    pairs.to_csv(direct_pairs_path, index=False)
    direct_summary.to_csv(direct_summary_path, index=False)
    trajectory.to_csv(trajectory_path, index=False)

    partial_path = outdir / "habitability_partial_temporal_effect.pdf"
    partial_effect_plot(temporal, partial_path)
    direct_plot1, direct_plot2 = make_direct_plots(
        pairs, trajectory, outdir
    )

    with summary_path.open("w", encoding="utf-8") as h:
        h.write("TERRAIN HABITABILITY H AS AN INDEPENDENT HISTORICAL FACTOR\n")
        h.write("=" * 76 + "\n\n")
        h.write(
            "H is entered as a separate physical-geography factor. "
            "Population density is NOT divided by H in these models.\n\n"
        )

        h.write("DIRECT MATCHED VIIRS DATA CHECK\n")
        h.write(
            f"Fixed NUTS release: {cfg['fixed_geography_nuts_release']}; "
            f"endpoint years: {cfg['direct_start_year']}-"
            f"{cfg['direct_end_year']}; nominal baseline density "
            f"<{cfg['direct_low_density_cutoff']} people km-2.\n"
        )
        h.write(
            f"H_core tertile cut points among matched sparse regions: "
            f"{q1:.6g}, {q2:.6g}.\n"
        )
        h.write(
            direct_summary.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write(
            "\n\nThese values are direct summaries of the matched VIIRS "
            "observations, not model predictions. Baseline radiance and "
            "absolute change are shown specifically to diagnose headroom.\n"
        )

        h.write("\nLONGITUDINAL H x TIME TESTS\n")
        h.write(
            temporal[
                [
                    "scenario",
                    "interval",
                    "nuts_release_filter",
                    "n_observations",
                    "n_region_keys",
                    "n_countries",
                    "H_x_year_beta_log10_per_year_per_H",
                    "H_x_year_se",
                    "lrt_lr",
                    "lrt_df",
                    "lrt_p",
                    "relative_annual_growth_lowH_0p1_vs_highH_1_pct",
                    "reduced_aic",
                    "full_aic",
                ]
            ].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write(
            "\n\nControls: country-specific temporal trends, cubic nominal "
            "population-density x time, H main effect, and a "
            "NUTS-release-specific random intercept.\n"
        )

        h.write("\nTEMPORALLY SEPARATED HEADROOM TEST\n")
        h.write(
            f"Baseline years: {cfg['headroom_baseline_years']}; "
            f"outcome trend: {cfg['headroom_outcome_first_year']}-"
            f"{cfg['headroom_outcome_last_year']}.\n"
        )
        h.write(
            headroom_models.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nH BLOCK AFTER HEADROOM CONTROLS\n")
        for key, value in headroom_test.items():
            h.write(f"{key}: {value}\n")

        h.write("\nINTERPRETATION GUIDE\n")
        h.write(
            "- First inspect the DIRECT VIIRS table and trajectories: this "
            "shows whether low-H regions visibly brightened more in the raw "
            "matched observations and whether they simply started darker.\n"
        )
        h.write(
            "- A significant longitudinal H x time test means H predicts "
            "change independently of nominal density and country-specific "
            "temporal trends.\n"
        )
        h.write(
            "- Agreement in the fixed NUTS2016 2013-2020 analysis argues "
            "against NUTS boundary-release changes creating the effect.\n"
        )
        h.write(
            "- The temporally separated headroom test uses 2013-2014 "
            "radiance only as a predictor and 2015-2020 change as the outcome. "
            "If H remains supported there, the result is not simply a low "
            "starting-radiance percentage effect.\n"
        )

    print("\nDIRECT MATCHED VIIRS DATA CHECK")
    print(
        direct_summary.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )

    print("\nLONGITUDINAL H x TIME TESTS")
    print(
        temporal[
            [
                "scenario",
                "interval",
                "H_x_year_beta_log10_per_year_per_H",
                "lrt_lr",
                "lrt_df",
                "lrt_p",
                "relative_annual_growth_lowH_0p1_vs_highH_1_pct",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )

    print("\nHEADROOM MODEL COMPARISON")
    print(
        headroom_models.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print("\nH BLOCK AFTER HEADROOM CONTROLS")
    for key, value in headroom_test.items():
        print(f"{key}: {value}")

    print("\nWROTE")
    for p in [
        temporal_path,
        headroom_table_path,
        headroom_models_path,
        headroom_test_path,
        direct_pairs_path,
        direct_summary_path,
        trajectory_path,
        summary_path,
        partial_path,
        direct_plot1,
        direct_plot2,
    ]:
        print(p)


if __name__ == "__main__":
    main()
