#!/usr/bin/env python3
"""Relate annual calibrated VIIRS slopes to paired ISS colour-ratio change.

Purpose
-------
The endpoint LED-switch analysis showed a very strong increase in ISS G/R
whitening with population density, and a descriptive spectral correction that
made dense-area VIIRS change approximately stable. The endpoint regression,
however, was not itself significant.

This stage aligns the response more directly with the historical result we
want to explain: the annual temporal slope of production-calibrated VIIRS
radiance in each NUTS3 region.

For each configured fixed-NUTS2016 window (2013-2019 and 2013-2020 by default)
we fit, separately for every NUTS3:

    log10(corrected VIIRS) ~ year

and also:

    log10(population density) ~ year

The paired ISS colour shift is then merged onto those regional slopes.

Spectral conversion
-------------------
Using Sanchez de Miguel et al. (2022):

    G_ISS / VIIRS = 0.21 * 1.5 ** (G_ISS / R_ISS)

so at equal green-band output:

    predicted spectral delta log10(VIIRS)
        = -log10(1.5) * delta(G/R)

The ISS pre/post interval is anchored approximately to 2013->2020, so the
spectral term is converted to an endpoint-equivalent annualized log10 slope by
dividing by 7 years. This is a descriptive normalization, NOT a claim that LED
conversion occurred linearly through time.

Formal tests
------------
Within each window we ask:

1. Does ISS whitening increase with baseline population density?
2. Does the spectral term explain additional variation in annual VIIRS slope
   after country fixed effects, baseline density and population-change slope?
3. Is there a density x spectral interaction?
4. Does the population-density coefficient attenuate when the theoretical
   spectral component is removed from each NUTS3 VIIRS slope?
5. Does the dense-vs-sparse contrast attenuate after spectral adjustment?
6. Within dense regions (>= configured cutoff), does the spectral term track
   VIIRS slope?

The ISS post mosaic is a 2014-2020 composite. Accordingly, this remains a
mechanism/sensitivity analysis rather than a year-resolved LED adoption model.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

from project_config import configured_path, load_config


LOG10_1P5 = float(np.log10(1.5))


def fit_simple_slope(years: np.ndarray, values: np.ndarray) -> dict:
    fit = linregress(years.astype(float), values.astype(float))
    return {
        "slope": float(fit.slope),
        "intercept": float(fit.intercept),
        "r2": float(fit.rvalue ** 2),
        "p": float(fit.pvalue),
        "se": float(fit.stderr),
    }


def load_inputs(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    processed = configured_path(config, "processed")

    merged_path = processed / config["outputs"]["merged"]
    merged = pd.read_parquet(merged_path).copy()

    required_merged = {
        "year",
        "NUTS_ID",
        "CNTR_CODE",
        "nuts_release",
        "radiance_mean_corrected",
        "population_density",
    }
    missing = required_merged - set(merged.columns)
    if missing:
        raise KeyError(
            f"Merged production table missing columns: {sorted(missing)}"
        )

    spectral_path = processed / "iss_nuts3_spectral_shift.csv"
    if not spectral_path.exists():
        raise FileNotFoundError(
            f"Missing {spectral_path}. Run: "
            "python src/run_pipeline.py iss-spectral-validation"
        )

    spectral = pd.read_csv(spectral_path)
    required_spectral = {
        "NUTS_ID",
        "CNTR_CODE",
        "n_common_valid_pixels",
        "delta_gr_median",
        "delta_bg_median",
    }
    missing = required_spectral - set(spectral.columns)
    if missing:
        raise KeyError(
            f"ISS spectral table missing columns: {sorted(missing)}"
        )

    return merged, spectral


def build_window_slopes(
    merged: pd.DataFrame,
    spectral: pd.DataFrame,
    config: dict,
    window_cfg: dict,
) -> pd.DataFrame:
    cfg = config["analysis"]["iss_led_switch_annual_slopes"]

    label = str(window_cfg["label"])
    first_year = int(window_cfg["first_year"])
    last_year = int(window_cfg["last_year"])
    nuts_release = int(cfg["nuts_release"])
    required_years = list(range(first_year, last_year + 1))
    require_complete = bool(cfg["require_complete_panel"])

    d = merged[
        merged["year"].between(first_year, last_year, inclusive="both")
        & (merged["nuts_release"] == nuts_release)
    ][
        [
            "year",
            "NUTS_ID",
            "CNTR_CODE",
            "radiance_mean_corrected",
            "population_density",
        ]
    ].copy()

    ok = (
        np.isfinite(pd.to_numeric(d["radiance_mean_corrected"], errors="coerce"))
        & np.isfinite(pd.to_numeric(d["population_density"], errors="coerce"))
        & (d["radiance_mean_corrected"] > 0)
        & (d["population_density"] > 0)
        & d["NUTS_ID"].notna()
        & d["CNTR_CODE"].notna()
    )
    d = d.loc[ok].copy()

    rows = []
    target_set = set(required_years)

    for (nuts_id, country), g in d.groupby(
        ["NUTS_ID", "CNTR_CODE"], sort=False
    ):
        g = g.sort_values("year")
        years_present = set(g["year"].astype(int))

        if require_complete:
            if years_present != target_set:
                continue
        else:
            if len(g) < int(cfg["minimum_panel_years"]):
                continue

        years = g["year"].to_numpy(dtype=float)
        log_rad = np.log10(g["radiance_mean_corrected"].to_numpy(dtype=float))
        log_pop = np.log10(g["population_density"].to_numpy(dtype=float))

        rad_fit = fit_simple_slope(years, log_rad)
        pop_fit = fit_simple_slope(years, log_pop)

        baseline_row = g[g["year"] == first_year]
        if baseline_row.empty:
            continue

        baseline_density = float(
            baseline_row["population_density"].iloc[0]
        )

        rows.append(
            {
                "window": label,
                "first_year": first_year,
                "last_year": last_year,
                "NUTS_ID": nuts_id,
                "CNTR_CODE": country,
                "n_years": int(len(g)),
                "population_density_baseline": baseline_density,
                "log_pop_baseline": float(np.log10(baseline_density)),
                "viirs_log10_slope_per_year": rad_fit["slope"],
                "viirs_slope_se": rad_fit["se"],
                "viirs_slope_p": rad_fit["p"],
                "viirs_slope_r2": rad_fit["r2"],
                "population_log10_slope_per_year": pop_fit["slope"],
                "population_slope_r2": pop_fit["r2"],
            }
        )

    slopes = pd.DataFrame(rows)
    if slopes.empty:
        return slopes

    min_pixels = int(cfg["min_valid_iss_pixels"])
    spec = spectral[
        spectral["n_common_valid_pixels"] >= min_pixels
    ][
        [
            "NUTS_ID",
            "CNTR_CODE",
            "n_common_valid_pixels",
            "delta_gr_median",
            "delta_bg_median",
        ]
    ].copy()

    slopes = slopes.merge(
        spec,
        on=["NUTS_ID", "CNTR_CODE"],
        how="inner",
        validate="one_to_one",
    )

    finite = (
        np.isfinite(pd.to_numeric(slopes["delta_gr_median"], errors="coerce"))
        & np.isfinite(pd.to_numeric(slopes["delta_bg_median"], errors="coerce"))
    )
    slopes = slopes.loc[finite].copy()

    spectral_years = float(cfg["spectral_annualization_years"])
    slopes["predicted_spectral_delta_log10"] = (
        -LOG10_1P5 * slopes["delta_gr_median"]
    )
    slopes["spectral_proxy_log10_per_year"] = (
        slopes["predicted_spectral_delta_log10"] / spectral_years
    )
    slopes["viirs_slope_spectral_adjusted"] = (
        slopes["viirs_log10_slope_per_year"]
        - slopes["spectral_proxy_log10_per_year"]
    )

    slopes["observed_pct_per_year"] = (
        100.0 * (10.0 ** slopes["viirs_log10_slope_per_year"] - 1.0)
    )
    slopes["spectral_proxy_pct_per_year"] = (
        100.0 * (10.0 ** slopes["spectral_proxy_log10_per_year"] - 1.0)
    )
    slopes["spectral_adjusted_pct_per_year"] = (
        100.0 * (10.0 ** slopes["viirs_slope_spectral_adjusted"] - 1.0)
    )

    edges = (
        [-np.inf]
        + [float(x) for x in cfg["density_band_breaks"]]
        + [np.inf]
    )
    labels = list(cfg["density_band_labels"])
    slopes["density_band"] = pd.cut(
        slopes["population_density_baseline"],
        bins=edges,
        labels=labels,
        right=False,
        ordered=True,
    )

    slopes["log_pop_within_country"] = (
        slopes["log_pop_baseline"]
        - slopes.groupby("CNTR_CODE")["log_pop_baseline"].transform("mean")
    )
    slopes["spectral_proxy_within_country"] = (
        slopes["spectral_proxy_log10_per_year"]
        - slopes.groupby("CNTR_CODE")[
            "spectral_proxy_log10_per_year"
        ].transform("mean")
    )

    cutoff = float(cfg["high_density_cutoff"])
    slopes["dense_vs_sparse"] = np.where(
        slopes["population_density_baseline"] >= cutoff,
        1.0,
        np.where(slopes["population_density_baseline"] < 100.0, 0.0, np.nan),
    )

    counts = slopes.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    keep = counts[counts >= int(cfg["min_country_regions"])].index
    slopes = slopes[slopes["CNTR_CODE"].isin(keep)].copy()
    slopes["CNTR_CODE"] = slopes["CNTR_CODE"].astype("category")

    return slopes


def country_balanced_band_summary(
    data: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    cfg = config["analysis"]["iss_led_switch_annual_slopes"]
    labels = list(cfg["density_band_labels"])

    country_band = (
        data.groupby(
            ["window", "CNTR_CODE", "density_band"],
            observed=True,
            as_index=False,
        )
        .agg(
            n_nuts3=("NUTS_ID", "nunique"),
            delta_gr=("delta_gr_median", "mean"),
            observed_log_slope=("viirs_log10_slope_per_year", "mean"),
            spectral_log_slope=("spectral_proxy_log10_per_year", "mean"),
            adjusted_log_slope=("viirs_slope_spectral_adjusted", "mean"),
        )
    )

    rows = []
    for window, w in country_band.groupby("window", observed=True):
        for band in labels:
            g = w[w["density_band"].astype(str) == band]
            if g.empty:
                continue

            row = {
                "window": window,
                "density_band": band,
                "n_countries": int(g["CNTR_CODE"].nunique()),
                "n_nuts3": int(g["n_nuts3"].sum()),
                "mean_country_delta_gr": float(g["delta_gr"].mean()),
                "median_country_delta_gr": float(g["delta_gr"].median()),
            }

            for label, col in [
                ("observed", "observed_log_slope"),
                ("spectral_proxy", "spectral_log_slope"),
                ("spectral_adjusted", "adjusted_log_slope"),
            ]:
                value = float(g[col].mean())
                row[f"{label}_mean_log10_per_year"] = value
                row[f"{label}_pct_per_year"] = 100.0 * (10.0 ** value - 1.0)

            rows.append(row)

    return pd.DataFrame(rows)


def model_row(
    fit,
    *,
    window: str,
    model: str,
    subset: str,
    outcome: str,
    term: str,
) -> dict:
    return {
        "window": window,
        "model": model,
        "subset": subset,
        "outcome": outcome,
        "term": term,
        "estimate": float(fit.params[term]),
        "se_hc3": float(fit.bse[term]),
        "p_hc3": float(fit.pvalues[term]),
        "n": int(fit.nobs),
    }


def nested_test(
    data: pd.DataFrame,
    reduced: str,
    full: str,
    *,
    window: str,
    label: str,
) -> dict:
    m0 = smf.ols(reduced, data=data).fit()
    m1 = smf.ols(full, data=data).fit()
    a = anova_lm(m0, m1).iloc[-1]
    return {
        "window": window,
        "test": label,
        "reduced_formula": reduced,
        "full_formula": full,
        "F": float(a["F"]),
        "df_num": float(a["df_diff"]),
        "p": float(a["Pr(>F)"]),
        "r2_reduced": float(m0.rsquared),
        "r2_full": float(m1.rsquared),
        "delta_r2": float(m1.rsquared - m0.rsquared),
    }


def fit_window_models(
    data: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = config["analysis"]["iss_led_switch_annual_slopes"]

    model_rows = []
    test_rows = []

    for window, d in data.groupby("window", observed=True):
        d = d.copy()

        # Whitening gradient on exactly the NUTS3 support used for this window.
        form = "delta_gr_median ~ C(CNTR_CODE) + log_pop_within_country"
        f = smf.ols(form, data=d).fit(cov_type="HC3")
        model_rows.append(
            model_row(
                f,
                window=window,
                model="whitening_density_gradient",
                subset="all",
                outcome="delta_gr_median",
                term="log_pop_within_country",
            )
        )

        base = (
            "viirs_log10_slope_per_year ~ C(CNTR_CODE) + "
            "population_log10_slope_per_year + log_pop_within_country"
        )
        spectral_main = base + " + spectral_proxy_within_country"
        spectral_interaction = (
            spectral_main
            + " + log_pop_within_country:spectral_proxy_within_country"
        )

        base_hc3 = smf.ols(base, data=d).fit(cov_type="HC3")
        spectral_hc3 = smf.ols(spectral_main, data=d).fit(cov_type="HC3")
        interaction_hc3 = smf.ols(
            spectral_interaction, data=d
        ).fit(cov_type="HC3")

        for model_name, fit, terms in [
            (
                "observed_base",
                base_hc3,
                ["log_pop_within_country"],
            ),
            (
                "observed_plus_spectral",
                spectral_hc3,
                [
                    "log_pop_within_country",
                    "spectral_proxy_within_country",
                ],
            ),
            (
                "observed_plus_spectral_interaction",
                interaction_hc3,
                [
                    "log_pop_within_country",
                    "spectral_proxy_within_country",
                    "log_pop_within_country:spectral_proxy_within_country",
                ],
            ),
        ]:
            for term in terms:
                model_rows.append(
                    model_row(
                        fit,
                        window=window,
                        model=model_name,
                        subset="all",
                        outcome="viirs_log10_slope_per_year",
                        term=term,
                    )
                )

        test_rows.append(
            nested_test(
                d,
                base,
                spectral_main,
                window=window,
                label="add_spectral_main_effect",
            )
        )
        test_rows.append(
            nested_test(
                d,
                spectral_main,
                spectral_interaction,
                window=window,
                label="add_density_x_spectral_interaction",
            )
        )
        test_rows.append(
            nested_test(
                d,
                base,
                spectral_interaction,
                window=window,
                label="add_spectral_main_plus_interaction",
            )
        )

        # Direct attenuation test: same controls, but outcome has the theoretical
        # spectral component removed.
        adjusted_formula = (
            "viirs_slope_spectral_adjusted ~ C(CNTR_CODE) + "
            "population_log10_slope_per_year + log_pop_within_country"
        )
        adjusted_hc3 = smf.ols(adjusted_formula, data=d).fit(cov_type="HC3")
        model_rows.append(
            model_row(
                adjusted_hc3,
                window=window,
                model="spectral_adjusted_density_effect",
                subset="all",
                outcome="viirs_slope_spectral_adjusted",
                term="log_pop_within_country",
            )
        )

        # Dense-vs-sparse contrast, excluding the intermediate 100-300 group.
        contrast = d[np.isfinite(d["dense_vs_sparse"])].copy()
        if len(contrast) >= 20:
            for outcome, model_name in [
                (
                    "viirs_log10_slope_per_year",
                    "dense_vs_sparse_observed",
                ),
                (
                    "viirs_slope_spectral_adjusted",
                    "dense_vs_sparse_spectral_adjusted",
                ),
            ]:
                formula = (
                    f"{outcome} ~ C(CNTR_CODE) + "
                    "population_log10_slope_per_year + dense_vs_sparse"
                )
                fit = smf.ols(formula, data=contrast).fit(cov_type="HC3")
                model_rows.append(
                    model_row(
                        fit,
                        window=window,
                        model=model_name,
                        subset="sparse<100_vs_dense>=cutoff",
                        outcome=outcome,
                        term="dense_vs_sparse",
                    )
                )

        cutoff = float(cfg["high_density_cutoff"])
        dense = d[d["population_density_baseline"] >= cutoff].copy()
        if (
            len(dense) >= int(cfg["minimum_dense_rows"])
            and dense["CNTR_CODE"].nunique() >= 3
        ):
            dense_formula = (
                "viirs_log10_slope_per_year ~ C(CNTR_CODE) + "
                "population_log10_slope_per_year + spectral_proxy_within_country"
            )
            dense_fit = smf.ols(
                dense_formula, data=dense
            ).fit(cov_type="HC3")
            model_rows.append(
                model_row(
                    dense_fit,
                    window=window,
                    model="high_density_spectral_effect",
                    subset=f"density>={cutoff:g}",
                    outcome="viirs_log10_slope_per_year",
                    term="spectral_proxy_within_country",
                )
            )

    return pd.DataFrame(model_rows), pd.DataFrame(test_rows)


def build_attenuation_table(models: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for window, g in models.groupby("window", observed=True):
        def get(model, term):
            z = g[(g["model"] == model) & (g["term"] == term)]
            if z.empty:
                return np.nan, np.nan, np.nan
            r = z.iloc[0]
            return float(r["estimate"]), float(r["se_hc3"]), float(r["p_hc3"])

        obs, obs_se, obs_p = get(
            "observed_base", "log_pop_within_country"
        )
        adj, adj_se, adj_p = get(
            "spectral_adjusted_density_effect",
            "log_pop_within_country",
        )

        ds_obs, ds_obs_se, ds_obs_p = get(
            "dense_vs_sparse_observed", "dense_vs_sparse"
        )
        ds_adj, ds_adj_se, ds_adj_p = get(
            "dense_vs_sparse_spectral_adjusted", "dense_vs_sparse"
        )

        rows.append(
            {
                "window": window,
                "density_effect_observed": obs,
                "density_effect_observed_se": obs_se,
                "density_effect_observed_p": obs_p,
                "density_effect_after_spectral_adjustment": adj,
                "density_effect_after_spectral_adjustment_se": adj_se,
                "density_effect_after_spectral_adjustment_p": adj_p,
                "density_effect_absolute_attenuation": (
                    abs(obs) - abs(adj)
                    if np.isfinite(obs) and np.isfinite(adj)
                    else np.nan
                ),
                "density_effect_fraction_attenuated": (
                    1.0 - abs(adj) / abs(obs)
                    if np.isfinite(obs)
                    and np.isfinite(adj)
                    and obs != 0
                    else np.nan
                ),
                "dense_sparse_observed": ds_obs,
                "dense_sparse_observed_se": ds_obs_se,
                "dense_sparse_observed_p": ds_obs_p,
                "dense_sparse_after_spectral_adjustment": ds_adj,
                "dense_sparse_after_spectral_adjustment_se": ds_adj_se,
                "dense_sparse_after_spectral_adjustment_p": ds_adj_p,
                "dense_sparse_absolute_attenuation": (
                    abs(ds_obs) - abs(ds_adj)
                    if np.isfinite(ds_obs) and np.isfinite(ds_adj)
                    else np.nan
                ),
                "dense_sparse_fraction_attenuated": (
                    1.0 - abs(ds_adj) / abs(ds_obs)
                    if np.isfinite(ds_obs)
                    and np.isfinite(ds_adj)
                    and ds_obs != 0
                    else np.nan
                ),
            }
        )

    return pd.DataFrame(rows)


def plot_band_summary(
    summary: pd.DataFrame,
    window: str,
    output: Path,
    dpi: int,
) -> None:
    s = summary[summary["window"] == window].copy()
    if s.empty:
        return

    x = np.arange(len(s))
    fig, ax = plt.subplots(figsize=(10.2, 6.2))

    ax.plot(
        x,
        s["observed_pct_per_year"],
        marker="o",
        linewidth=2.0,
        label="Observed calibrated VIIRS slope",
    )
    ax.plot(
        x,
        s["spectral_proxy_pct_per_year"],
        marker="s",
        linestyle="--",
        linewidth=1.8,
        label="ISS G/R endpoint-equivalent spectral term",
    )
    ax.plot(
        x,
        s["spectral_adjusted_pct_per_year"],
        marker="D",
        linestyle="-.",
        linewidth=1.8,
        label="After removing spectral term",
    )

    ax.axhline(0, color="0.4", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(
        s["density_band"].str.replace(">=", "≥", regex=False)
    )
    ax.set_xlabel("Baseline population-density band (people km$^{-2}$)")
    ax.set_ylabel("NUTS3 temporal slope (% yr$^{-1}$)")
    ax.set_title(
        f"Annual VIIRS slopes and ISS LED-switch spectral term: {window}"
    )
    ax.grid(axis="y", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_spectral_vs_slope(
    data: pd.DataFrame,
    window: str,
    output: Path,
    dpi: int,
) -> None:
    d = data[data["window"] == window].copy()
    if d.empty:
        return

    fig, ax = plt.subplots(figsize=(8.2, 6.3))
    labels = list(d["density_band"].cat.categories)
    markers = ["o", "s", "D", "^", "X"]

    for i, band in enumerate(labels):
        g = d[d["density_band"] == band]
        if g.empty:
            continue
        ax.scatter(
            g["spectral_proxy_log10_per_year"],
            g["viirs_log10_slope_per_year"],
            s=18,
            alpha=0.45,
            marker=markers[i % len(markers)],
            label=str(band).replace(">=", "≥"),
        )

    ax.axhline(0, color="0.4", linewidth=1)
    ax.axvline(0, color="0.4", linewidth=1)
    ax.set_xlabel("ISS G/R spectral slope proxy (log10 yr$^{-1}$)")
    ax.set_ylabel("Observed calibrated VIIRS slope (log10 yr$^{-1}$)")
    ax.set_title(f"ISS spectral change versus annual VIIRS slope: {window}")
    ax.grid(alpha=0.2)
    ax.legend(title="Density band", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    path: Path,
    slopes: pd.DataFrame,
    bands: pd.DataFrame,
    models: pd.DataFrame,
    tests: pd.DataFrame,
    attenuation: pd.DataFrame,
    config: dict,
) -> None:
    cfg = config["analysis"]["iss_led_switch_annual_slopes"]

    with path.open("w", encoding="utf-8") as h:
        h.write("ISS LED-SWITCH TEST USING ANNUAL CALIBRATED VIIRS SLOPES\n")
        h.write("=" * 78 + "\n\n")
        h.write(
            "Response: NUTS3 OLS slope of log10(production-calibrated VIIRS) "
            "against year.\n"
        )
        h.write(
            "Spectral proxy: -log10(1.5) * paired delta(G/R), divided by "
            f"{cfg['spectral_annualization_years']} years.\n"
        )
        h.write(
            "IMPORTANT: the ISS post mosaic is a 2014-2020 composite; this "
            "annualization is an endpoint-equivalent normalization, not an "
            "assumed linear LED-adoption trajectory.\n\n"
        )

        h.write("WINDOW SUPPORT\n")
        support = (
            slopes.groupby("window", observed=True)
            .agg(
                n_nuts3=("NUTS_ID", "nunique"),
                n_countries=("CNTR_CODE", "nunique"),
                median_viirs_slope_r2=("viirs_slope_r2", "median"),
            )
            .reset_index()
        )
        h.write(
            support.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nCOUNTRY-BALANCED DENSITY-BAND SUMMARY\n")
        h.write(
            bands.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nMODEL COEFFICIENTS (HC3)\n")
        h.write(
            models.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nNESTED SPECTRAL TESTS\n")
        h.write(
            tests.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nDENSITY-EFFECT ATTENUATION\n")
        h.write(
            attenuation.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nINTERPRETATION GUIDE\n")
        h.write(
            "- Positive delta G/R is spectral whitening and gives a negative "
            "VIIRS spectral slope proxy at equal green output.\n"
        )
        h.write(
            "- A significant spectral block indicates that regional ISS "
            "whitening explains additional variation in annual VIIRS slope "
            "after country, baseline density and population change.\n"
        )
        h.write(
            "- If the absolute population-density coefficient or the "
            "dense-vs-sparse coefficient becomes smaller after spectral "
            "adjustment, LED spectral conversion explains part of the observed "
            "density rotation.\n"
        )
        h.write(
            "- The 2013-2019 window is the cleaner pre-COVID sensitivity. "
            "The 2013-2020 window aligns better with the ISS post composite but "
            "includes the unusual 2020 year. Agreement between them is stronger "
            "evidence than either alone.\n"
        )
        h.write(
            "- This test cannot equate the spectral-adjusted residual with true "
            "lighting output: equal-green normalization is not constant lumens, "
            "installed power or biological exposure.\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["iss_led_switch_annual_slopes"]

    merged, spectral = load_inputs(config)

    frames = []
    for window_cfg in cfg["windows"]:
        frame = build_window_slopes(
            merged,
            spectral,
            config,
            window_cfg,
        )
        if not frame.empty:
            frames.append(frame)

    if not frames:
        raise RuntimeError("No annual-slope datasets passed configured filters.")

    slopes = pd.concat(frames, ignore_index=True)
    bands = country_balanced_band_summary(slopes, config)
    models, tests = fit_window_models(slopes, config)
    attenuation = build_attenuation_table(models)

    processed = configured_path(config, "processed")
    figure_dir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "spectral_validation"
        / "led_switch_annual_slopes"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "slopes": processed / "iss_led_switch_annual_nuts3_slopes.csv",
        "bands": processed / "iss_led_switch_annual_density_bands.csv",
        "models": processed / "iss_led_switch_annual_models.csv",
        "tests": processed / "iss_led_switch_annual_tests.csv",
        "attenuation": processed / "iss_led_switch_annual_attenuation.csv",
        "summary": processed / "iss_led_switch_annual_summary.txt",
    }

    slopes.to_csv(paths["slopes"], index=False)
    bands.to_csv(paths["bands"], index=False)
    models.to_csv(paths["models"], index=False)
    tests.to_csv(paths["tests"], index=False)
    attenuation.to_csv(paths["attenuation"], index=False)

    dpi = int(cfg["figure_dpi"])
    figure_paths = []
    for window in slopes["window"].drop_duplicates():
        safe = str(window).replace("/", "_").replace(" ", "_")
        p1 = figure_dir / f"{safe}_density_decomposition.png"
        p2 = figure_dir / f"{safe}_spectral_vs_viirs_slope.png"
        plot_band_summary(bands, window, p1, dpi)
        plot_spectral_vs_slope(slopes, window, p2, dpi)
        figure_paths.extend([p1, p2])

    write_summary(
        paths["summary"],
        slopes,
        bands,
        models,
        tests,
        attenuation,
        config,
    )

    print("\nWINDOW SUPPORT")
    print(
        slopes.groupby("window", observed=True)
        .agg(
            n_nuts3=("NUTS_ID", "nunique"),
            n_countries=("CNTR_CODE", "nunique"),
            median_viirs_slope_r2=("viirs_slope_r2", "median"),
        )
        .reset_index()
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )

    print("\nCOUNTRY-BALANCED DENSITY-BAND SUMMARY")
    show_cols = [
        "window",
        "density_band",
        "n_countries",
        "n_nuts3",
        "mean_country_delta_gr",
        "observed_pct_per_year",
        "spectral_proxy_pct_per_year",
        "spectral_adjusted_pct_per_year",
    ]
    print(
        bands[show_cols].to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )

    print("\nNESTED SPECTRAL TESTS")
    print(
        tests.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )

    print("\nDENSITY-EFFECT ATTENUATION")
    print(
        attenuation.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )

    print("\nWROTE")
    for p in [*paths.values(), *figure_paths]:
        print(p)


if __name__ == "__main__":
    main()
