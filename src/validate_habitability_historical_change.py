#!/usr/bin/env python3
"""Historically validate the static terrain habitability metric H.

This stage asks whether terrain habitability explains which nominally sparse
NUTS3 regions actually brightened over 2013-2024.

It performs four complementary tests:

1. Refit the headline temporal density-band analysis after replacing nominal
   density with terrain-adjusted density D_H = D / H.
2. Compare common continuous density x time surfaces fitted with nominal versus
   terrain-adjusted density on exactly the same observations.
3. Test whether H modifies the temporal trend after controlling for the full
   cubic nominal-density x time surface.
4. Within nominally sparse regions (<100 people km-2), test whether low-, mid-,
   and high-H regions have different historical radiance trends.

H is joined by NUTS_ID + NUTS release.  A region_key combining both fields is
used as the random-intercept group so that boundary-release changes are not
silently treated as an unchanged repeated spatial unit.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
import statsmodels.formula.api as smf

from project_config import configured_path, load_config


def fit_mixedlm(formula: str, data: pd.DataFrame):
    # Subsetting a pandas Categorical retains unused levels. Patsy then creates
    # all-zero dummy columns for those absent levels, which makes the fixed-
    # effect design matrix singular. This matters especially for the low-density
    # H-tertile subset, where some countries present in the full dataset are
    # absent. Remove unused levels before every model fit.
    fit_data = data.copy()
    for column in fit_data.columns:
        if isinstance(fit_data[column].dtype, pd.CategoricalDtype):
            fit_data[column] = fit_data[column].cat.remove_unused_categories()

    model = smf.mixedlm(
        formula,
        data=fit_data,
        groups=fit_data["region_key"],
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
                    maxiter=3000,
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


def nakagawa_r2(result):
    fixed_pred = np.asarray(
        result.model.exog @ result.fe_params,
        dtype=float,
    )
    var_fixed = float(np.var(fixed_pred, ddof=1))
    var_random = float(np.asarray(result.cov_re)[0, 0])
    var_resid = float(result.scale)
    total = var_fixed + var_random + var_resid
    if total <= 0:
        return np.nan, np.nan
    return (
        var_fixed / total,
        (var_fixed + var_random) / total,
    )


def contrast(result, weights):
    names = list(result.fe_params.index)
    L = np.zeros(len(names), dtype=float)
    for name, value in weights.items():
        if name not in names:
            raise KeyError(
                f"Parameter {name!r} not found. Available: {names}"
            )
        L[names.index(name)] = value
    beta = result.fe_params.to_numpy(dtype=float)
    cov = result.cov_params().loc[names, names].to_numpy(dtype=float)
    est = float(L @ beta)
    var = float(L @ cov @ L)
    se = float(np.sqrt(max(0.0, var)))
    z = est / se if se > 0 else np.nan
    p = float(2.0 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
    return est, se, p


def prepare_data(config: dict) -> pd.DataFrame:
    cfg = config["analysis"]["habitability_historical_validation"]
    processed = configured_path(config, "processed")

    merged = pd.read_parquet(
        processed / config["outputs"]["merged"]
    ).copy()
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
    habit = habit[hcols].drop_duplicates(
        ["nuts_release", "NUTS_ID"]
    )

    data = merged.merge(
        habit,
        on=["nuts_release", "NUTS_ID"],
        how="inner",
        validate="many_to_one",
    )

    data = data[
        data["year"].between(
            int(cfg["first_year"]),
            int(cfg["last_year"]),
            inclusive="both",
        )
    ].copy()

    valid = (
        np.isfinite(
            pd.to_numeric(
                data["radiance_mean_corrected"],
                errors="coerce",
            )
        )
        & np.isfinite(
            pd.to_numeric(
                data["population_density"],
                errors="coerce",
            )
        )
        & (data["radiance_mean_corrected"] > 0)
        & (data["population_density"] > 0)
        & data["CNTR_CODE"].notna()
        & data["NUTS_ID"].notna()
    )
    data = data.loc[valid].copy()

    data["region_key"] = (
        data["NUTS_ID"].astype(str)
        + "::"
        + data["nuts_release"].astype(str)
    )

    repeated = data.groupby("region_key").size()
    data = data[
        data["region_key"].isin(
            repeated[
                repeated >= int(cfg["min_repeated_years"])
            ].index
        )
    ].copy()

    country_counts = data.groupby("CNTR_CODE")["region_key"].nunique()
    keep_countries = country_counts[
        country_counts >= int(cfg["min_country_regions"])
    ].index
    data = data[data["CNTR_CODE"].isin(keep_countries)].copy()

    data["log_rad"] = np.log10(
        data["radiance_mean_corrected"].astype(float)
    )
    data["country"] = data["CNTR_CODE"].astype("category")
    midpoint = (
        int(cfg["first_year"]) + int(cfg["last_year"])
    ) / 2.0
    data["year_c"] = data["year"].astype(float) - midpoint

    min_h = float(
        config["habitability"][
            "min_habitable_fraction_for_adjusted_density"
        ]
    )

    for scenario in ("strict", "core", "permissive"):
        hcol = f"habitable_fraction_{scenario}"
        dcol = f"terrain_adjusted_density_{scenario}"
        good_h = (
            np.isfinite(data[hcol])
            & (data[hcol] >= min_h)
        )
        data[dcol] = np.nan
        data.loc[good_h, dcol] = (
            data.loc[good_h, "population_density"]
            / data.loc[good_h, hcol]
        )

    return data


def add_reference_density_and_band(
    data: pd.DataFrame,
    density_col: str,
    prefix: str,
    edges,
    labels,
):
    ref = (
        data.groupby("region_key", as_index=False)[density_col]
        .median()
        .rename(
            columns={
                density_col: f"reference_density_{prefix}"
            }
        )
    )
    data = data.merge(
        ref,
        on="region_key",
        how="left",
        validate="many_to_one",
    )
    data[f"band_{prefix}"] = pd.cut(
        data[f"reference_density_{prefix}"],
        bins=edges,
        labels=labels,
        right=False,
        ordered=True,
    )
    return data


def fit_band_model(data, band_col, labels, metric_name):
    fit = data[data[band_col].notna()].copy()
    fit["density_band"] = fit[band_col].astype("category")

    reduced_formula = (
        "log_rad ~ year_c + C(density_band) + C(country)"
    )
    full_formula = (
        "log_rad ~ year_c * C(density_band) + C(country)"
    )
    reduced, opt0 = fit_mixedlm(reduced_formula, fit)
    full, opt1 = fit_mixedlm(full_formula, fit)
    lr, df, p = lrt(reduced, full)

    rows = []
    ref = labels[0]
    for band in labels:
        if band not in set(fit["density_band"].astype(str)):
            continue
        weights = {"year_c": 1.0}
        if band != ref:
            candidates = [
                f"year_c:C(density_band)[T.{band}]",
                f"C(density_band)[T.{band}]:year_c",
            ]
            found = next(
                (x for x in candidates if x in full.fe_params.index),
                None,
            )
            if found is None:
                continue
            weights[found] = 1.0

        est, se, pp = contrast(full, weights)
        rows.append(
            {
                "density_metric": metric_name,
                "density_band": band,
                "n_observations": int(
                    (fit["density_band"].astype(str) == band).sum()
                ),
                "n_region_keys": int(
                    fit.loc[
                        fit["density_band"].astype(str) == band,
                        "region_key",
                    ].nunique()
                ),
                "log10_radiance_trend_per_year": est,
                "trend_se": se,
                "trend_p": pp,
                "pct_radiance_change_per_year": (
                    100.0 * (10.0**est - 1.0)
                ),
                "pct_radiance_change_2013_2024_modelled": (
                    100.0 * (10.0 ** (11.0 * est) - 1.0)
                ),
            }
        )

    info = {
        "density_metric": metric_name,
        "interaction_lr": lr,
        "interaction_df": df,
        "interaction_p": p,
        "reduced_optimizer": opt0,
        "full_optimizer": opt1,
        "n_observations": len(fit),
        "n_region_keys": fit["region_key"].nunique(),
        "n_countries": fit["CNTR_CODE"].nunique(),
    }
    return pd.DataFrame(rows), info


def fit_continuous_density_model(
    data,
    density_col,
    model_name,
):
    fit = data[
        np.isfinite(data[density_col])
        & (data[density_col] > 0)
    ].copy()

    fit["log_density"] = np.log10(fit[density_col].astype(float))
    centre = float(fit["log_density"].median())
    fit["x"] = fit["log_density"] - centre

    formula = (
        "log_rad ~ C(country) + "
        "year_c * (x + I(x ** 2) + I(x ** 3))"
    )
    result, optimizer = fit_mixedlm(formula, fit)
    marginal, conditional = nakagawa_r2(result)

    row = {
        "model": model_name,
        "density_column": density_col,
        "formula": formula,
        "n_observations": len(fit),
        "n_region_keys": fit["region_key"].nunique(),
        "n_countries": fit["CNTR_CODE"].nunique(),
        "n_fixed_parameters": len(result.fe_params),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "bic": float(result.bic),
        "marginal_r2": marginal,
        "conditional_r2": conditional,
        "random_intercept_variance": float(
            np.asarray(result.cov_re)[0, 0]
        ),
        "residual_variance": float(result.scale),
        "optimizer": optimizer,
        "log_density_center": centre,
    }
    return result, fit, row


def fit_core_h_temporal_test(data):
    fit = data[
        np.isfinite(data["habitable_fraction_core"])
        & np.isfinite(data["population_density"])
        & (data["population_density"] > 0)
    ].copy()

    fit["log_pop"] = np.log10(
        fit["population_density"].astype(float)
    )
    fit["x"] = fit["log_pop"] - float(fit["log_pop"].median())
    fit["h_c"] = (
        fit["habitable_fraction_core"]
        - float(fit["habitable_fraction_core"].median())
    )

    # Spatial H terms are present in both models.  The full model asks whether
    # H also changes the temporal process after the cubic nominal-density x
    # time surface has already been fitted.
    reduced_formula = (
        "log_rad ~ C(country) + "
        "year_c * (x + I(x ** 2) + I(x ** 3)) + "
        "h_c + x:h_c + I(x ** 2):h_c + I(x ** 3):h_c"
    )
    full_formula = (
        reduced_formula
        + " + year_c:h_c"
        + " + year_c:x:h_c"
        + " + year_c:I(x ** 2):h_c"
        + " + year_c:I(x ** 3):h_c"
    )

    reduced, opt0 = fit_mixedlm(reduced_formula, fit)
    full, opt1 = fit_mixedlm(full_formula, fit)
    lr, df, p = lrt(reduced, full)

    return {
        "test": "Does H_core modify the temporal density-radiance surface?",
        "lr": lr,
        "df": df,
        "p": p,
        "reduced_formula": reduced_formula,
        "full_formula": full_formula,
        "reduced_optimizer": opt0,
        "full_optimizer": opt1,
        "n_observations": len(fit),
        "n_region_keys": fit["region_key"].nunique(),
        "n_countries": fit["CNTR_CODE"].nunique(),
    }


def fit_low_density_h_test(
    data,
    cutoff,
):
    region = (
        data[
            [
                "region_key",
                "NUTS_ID",
                "nuts_release",
                "CNTR_CODE",
                "habitable_fraction_core",
                "reference_density_nominal",
            ]
        ]
        .drop_duplicates("region_key")
        .copy()
    )
    region = region[
        region["reference_density_nominal"] < float(cutoff)
    ].copy()

    q1 = float(region["habitable_fraction_core"].quantile(1.0 / 3.0))
    q2 = float(region["habitable_fraction_core"].quantile(2.0 / 3.0))

    region["H_tertile"] = pd.cut(
        region["habitable_fraction_core"],
        bins=[-np.inf, q1, q2, np.inf],
        labels=["low_H", "mid_H", "high_H"],
        right=True,
        ordered=True,
    )

    fit = data.merge(
        region[["region_key", "H_tertile"]],
        on="region_key",
        how="inner",
        validate="many_to_one",
    )
    fit["H_tertile"] = fit["H_tertile"].astype("category")
    fit["log_pop"] = np.log10(fit["population_density"].astype(float))
    fit["x"] = fit["log_pop"] - float(fit["log_pop"].median())

    reduced_formula = (
        "log_rad ~ year_c + C(H_tertile) + "
        "x + year_c:x + C(country)"
    )
    full_formula = (
        "log_rad ~ year_c * C(H_tertile) + "
        "x + year_c:x + C(country)"
    )

    reduced, opt0 = fit_mixedlm(reduced_formula, fit)
    full, opt1 = fit_mixedlm(full_formula, fit)
    lr, df, p = lrt(reduced, full)

    rows = []
    for group in ["low_H", "mid_H", "high_H"]:
        weights = {"year_c": 1.0}
        if group != "low_H":
            candidates = [
                f"year_c:C(H_tertile)[T.{group}]",
                f"C(H_tertile)[T.{group}]:year_c",
            ]
            found = next(
                (x for x in candidates if x in full.fe_params.index),
                None,
            )
            if found is None:
                raise KeyError(
                    f"Could not find year interaction for {group}"
                )
            weights[found] = 1.0

        est, se, pp = contrast(full, weights)
        rkeys = region.loc[
            region["H_tertile"].astype(str) == group,
            "region_key",
        ]
        rows.append(
            {
                "H_group": group,
                "H_lower_bound": (
                    float(region.loc[
                        region["H_tertile"].astype(str) == group,
                        "habitable_fraction_core",
                    ].min())
                ),
                "H_upper_bound": (
                    float(region.loc[
                        region["H_tertile"].astype(str) == group,
                        "habitable_fraction_core",
                    ].max())
                ),
                "n_region_keys": int(rkeys.nunique()),
                "log10_radiance_trend_per_year": est,
                "trend_se": se,
                "trend_p": pp,
                "pct_radiance_change_per_year": (
                    100.0 * (10.0**est - 1.0)
                ),
                "pct_radiance_change_2013_2024_modelled": (
                    100.0 * (10.0 ** (11.0 * est) - 1.0)
                ),
            }
        )

    info = {
        "test": (
            f"Year x H_core tertile within nominal density <{cutoff:g}"
        ),
        "nominal_density_cutoff": float(cutoff),
        "H_tertile_q1": q1,
        "H_tertile_q2": q2,
        "lr": lr,
        "df": df,
        "p": p,
        "reduced_optimizer": opt0,
        "full_optimizer": opt1,
        "n_observations": len(fit),
        "n_region_keys": fit["region_key"].nunique(),
        "n_countries": fit["CNTR_CODE"].nunique(),
    }
    return pd.DataFrame(rows), info


def plot_band_comparison(bands, outdir: Path):
    fig, ax = plt.subplots(figsize=(9.4, 6.0))
    labels = list(dict.fromkeys(bands["density_band"].tolist()))
    x = np.arange(len(labels), dtype=float)

    metrics = [
        "nominal",
        "terrain_adjusted_core",
        "terrain_adjusted_strict",
        "terrain_adjusted_permissive",
    ]
    for metric in metrics:
        g = bands[bands["density_metric"] == metric].copy()
        if g.empty:
            continue
        g = g.set_index("density_band").reindex(labels)
        ax.plot(
            x,
            g["pct_radiance_change_per_year"],
            marker="o",
            label=metric.replace("_", " "),
        )

    ax.axhline(0.0, linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Fixed population-density band (people km⁻²)")
    ax.set_ylabel("Estimated corrected-VIIRS radiance change (% yr⁻¹)")
    ax.set_title(
        "Historical temporal trends: nominal vs terrain-adjusted density"
    )
    ax.grid(alpha=0.2, axis="y")
    ax.legend()
    fig.tight_layout()
    path = outdir / "nominal_vs_terrain_adjusted_density_band_trends.pdf"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_low_density_h_trends(rows, outdir: Path):
    fig, ax = plt.subplots(figsize=(7.8, 5.5))
    x = np.arange(len(rows), dtype=float)
    y = rows["pct_radiance_change_per_year"].to_numpy(float)
    log_est = rows["log10_radiance_trend_per_year"].to_numpy(float)
    se = rows["trend_se"].to_numpy(float)
    lo = 100.0 * (10.0 ** (log_est - 1.96 * se) - 1.0)
    hi = 100.0 * (10.0 ** (log_est + 1.96 * se) - 1.0)
    ax.errorbar(
        x,
        y,
        yerr=[y - lo, hi - y],
        fmt="o",
        capsize=4,
    )
    ax.axhline(0.0, linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(["Low H", "Mid H", "High H"])
    ax.set_xlabel("Terrain habitability tertile within nominal density <100")
    ax.set_ylabel("Estimated corrected-VIIRS radiance change (% yr⁻¹)")
    ax.set_title(
        "Does terrain habitability explain brightening among sparse regions?"
    )
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    path = outdir / "low_density_radiance_trend_by_habitability.pdf"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["habitability_historical_validation"]

    data = prepare_data(config)
    edges = (
        [-np.inf]
        + [float(x) for x in cfg["density_band_breaks"]]
        + [np.inf]
    )
    labels = list(cfg["density_band_labels"])

    data = add_reference_density_and_band(
        data,
        "population_density",
        "nominal",
        edges,
        labels,
    )

    all_band_rows = []
    band_tests = []

    nominal_rows, nominal_info = fit_band_model(
        data,
        "band_nominal",
        labels,
        "nominal",
    )
    all_band_rows.append(nominal_rows)
    band_tests.append(nominal_info)

    continuous_rows = []

    # Primary nominal model.
    nominal_model, nominal_fit, nominal_row = fit_continuous_density_model(
        data,
        "population_density",
        "nominal_density",
    )
    continuous_rows.append(nominal_row)

    for scenario in ("strict", "core", "permissive"):
        dcol = f"terrain_adjusted_density_{scenario}"
        prefix = f"adjusted_{scenario}"

        scen = data[
            np.isfinite(data[dcol])
            & (data[dcol] > 0)
        ].copy()

        scen = add_reference_density_and_band(
            scen,
            dcol,
            prefix,
            edges,
            labels,
        )

        band_rows, band_info = fit_band_model(
            scen,
            f"band_{prefix}",
            labels,
            f"terrain_adjusted_{scenario}",
        )
        all_band_rows.append(band_rows)
        band_tests.append(band_info)

        # To make nominal-vs-adjusted AIC/BIC comparison fair for each
        # sensitivity scenario, refit the nominal model on exactly the same
        # observations as the adjusted model.
        _, _, paired_nominal_row = fit_continuous_density_model(
            scen,
            "population_density",
            f"nominal_on_{scenario}_valid_sample",
        )
        _, _, adjusted_row = fit_continuous_density_model(
            scen,
            dcol,
            f"terrain_adjusted_{scenario}",
        )
        continuous_rows.extend(
            [paired_nominal_row, adjusted_row]
        )

    bands = pd.concat(all_band_rows, ignore_index=True)
    band_tests_df = pd.DataFrame(band_tests)
    continuous = pd.DataFrame(continuous_rows)

    h_test = fit_core_h_temporal_test(data)
    low_h_rows, low_h_info = fit_low_density_h_test(
        data,
        float(cfg["low_density_cutoff"]),
    )

    processed = configured_path(config, "processed")
    outdir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "habitability_validation"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    bands_path = processed / "habitability_density_band_temporal_trends.csv"
    band_tests_path = processed / "habitability_density_band_model_tests.csv"
    continuous_path = processed / "habitability_density_model_comparison.csv"
    h_test_path = processed / "habitability_temporal_modifier_test.csv"
    low_h_path = processed / "habitability_low_density_h_trends.csv"
    low_h_test_path = processed / "habitability_low_density_h_test.csv"
    summary_path = processed / "habitability_historical_validation_summary.txt"

    bands.to_csv(bands_path, index=False)
    band_tests_df.to_csv(band_tests_path, index=False)
    continuous.to_csv(continuous_path, index=False)
    pd.DataFrame([h_test]).to_csv(h_test_path, index=False)
    low_h_rows.to_csv(low_h_path, index=False)
    pd.DataFrame([low_h_info]).to_csv(low_h_test_path, index=False)

    plot1 = plot_band_comparison(bands, outdir)
    plot2 = plot_low_density_h_trends(low_h_rows, outdir)

    with summary_path.open("w", encoding="utf-8") as h:
        h.write("HISTORICAL VALIDATION OF TERRAIN HABITABILITY H\n")
        h.write("=" * 72 + "\n\n")
        h.write(
            f"Interval: {cfg['first_year']}-{cfg['last_year']}\n"
        )
        h.write(
            f"Observations after joining H: {len(data)}; "
            f"region-release groups: {data['region_key'].nunique()}; "
            f"countries: {data['CNTR_CODE'].nunique()}.\n"
        )
        h.write(
            "Random-intercept group is NUTS_ID::nuts_release so a boundary "
            "release change is not treated as the same repeated geometry.\n\n"
        )

        h.write("DENSITY-BAND TEMPORAL TRENDS\n")
        h.write(
            bands.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nBAND x YEAR TESTS\n")
        h.write(
            band_tests_df.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nCONTINUOUS MODEL COMPARISON\n")
        h.write(
            continuous[
                [
                    "model",
                    "n_observations",
                    "n_region_keys",
                    "log_likelihood",
                    "aic",
                    "bic",
                    "marginal_r2",
                    "conditional_r2",
                    "residual_variance",
                ]
            ].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nDOES H MODIFY THE TEMPORAL SURFACE?\n")
        for key, value in h_test.items():
            if key not in {"reduced_formula", "full_formula"}:
                h.write(f"{key}: {value}\n")

        h.write(
            "\nLOW-DENSITY (<100 people km-2 nominal) H-TERTILE TRENDS\n"
        )
        h.write(
            low_h_rows.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nLOW-DENSITY H-TERTILE INTERACTION TEST\n")
        for key, value in low_h_info.items():
            h.write(f"{key}: {value}\n")

        h.write("\nINTERPRETATION GUIDE\n")
        h.write(
            "- If terrain-adjusted density gives lower AIC/BIC or higher "
            "marginal R2 than nominal density on the same sample, H improves "
            "the historical density-time representation.\n"
        )
        h.write(
            "- If the H temporal-modifier LRT is significant, H explains "
            "historical temporal structure beyond nominal population density.\n"
        )
        h.write(
            "- The low-density H-tertile analysis is the direct Mont-Blanc "
            "versus developable-rural-land test: a stronger trend in high-H "
            "regions would show that nominally sparse but physically habitable "
            "landscapes drive more of the rural brightening.\n"
        )
        h.write(
            "- Strict/core/permissive adjusted-density results are sensitivity "
            "tests of the terrain thresholds; core remains the primary metric.\n"
        )

    print("\nLOW-DENSITY (<100) H-TERTILE TEMPORAL TRENDS")
    print(
        low_h_rows.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print(
        f"\nYear x H-tertile LRT: LR={low_h_info['lr']:.6g}, "
        f"df={low_h_info['df']}, p={low_h_info['p']:.6g}"
    )

    print("\nDOES H MODIFY THE FULL TEMPORAL DENSITY SURFACE?")
    print(
        f"LR={h_test['lr']:.6g}, df={h_test['df']}, "
        f"p={h_test['p']:.6g}"
    )

    print("\nCONTINUOUS MODEL COMPARISON")
    print(
        continuous[
            [
                "model",
                "aic",
                "bic",
                "marginal_r2",
                "conditional_r2",
                "residual_variance",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )

    print("\nWROTE")
    for p in [
        bands_path,
        band_tests_path,
        continuous_path,
        h_test_path,
        low_h_path,
        low_h_test_path,
        summary_path,
        plot1,
        plot2,
    ]:
        print(p)


if __name__ == "__main__":
    main()
