#!/usr/bin/env python3
"""Sensitivity of population-density VIIRS trends to excluding 2013.

Uses the same fixed NUTS2016 complete-panel regions as the historical
validation stage and changes only the temporal window: 2013-2020 versus
2014-2020. This directly tests whether the 2013 endpoint is driving the
population-density pattern.

Outputs:
- country-adjusted mixed-model density-class trends for both windows;
- real trajectories normalized to 2014 = 100;
- NUTS3 observed trend maps for both windows and their difference.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import linregress

from analyze_density_temporal_change import fit_band_model
from project_config import configured_path, load_config
from validate_future_radiance_history import (
    load_fixed_panel,
    load_geometry,
)


def fit_window(panel, first_year, last_year, config):
    d = panel[panel["year"].between(first_year, last_year)].copy()
    d["log_rad"] = d["log10_radiance"].astype(float)
    d["year_c"] = d["year"] - ((first_year + last_year) / 2.0)
    d["country"] = d["CNTR_CODE"].astype("category")
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    trends, info = fit_band_model(d, labels)
    trends["first_year"] = first_year
    trends["last_year"] = last_year
    trends["interval"] = f"{first_year}-{last_year}"
    trends["pct_change_per_year_ci_low"] = 100.0 * (
        10.0 ** (
            trends["log10_radiance_trend_per_year"] - 1.96 * trends["trend_se"]
        ) - 1.0
    )
    trends["pct_change_per_year_ci_high"] = 100.0 * (
        10.0 ** (
            trends["log10_radiance_trend_per_year"] + 1.96 * trends["trend_se"]
        ) - 1.0
    )
    return trends, info


def region_slopes(panel, first_year, last_year):
    d = panel[panel["year"].between(first_year, last_year)].copy()
    rows = []
    for nuts_id, g in d.groupby("NUTS_ID"):
        g = g.sort_values("year")
        fit = linregress(g["year"], g["log10_radiance"])
        rows.append({
            "NUTS_ID": nuts_id,
            "CNTR_CODE": str(g["CNTR_CODE"].iloc[0]),
            "density_band": str(g["density_band"].iloc[0]),
            "reference_population_density": float(
                g["reference_population_density"].iloc[0]
            ),
            "n_years": int(len(g)),
            "observed_log10_slope_per_year": float(fit.slope),
            "observed_pct_change_per_year": float(
                100.0 * (10.0 ** fit.slope - 1.0)
            ),
        })
    return pd.DataFrame(rows)


def trajectories_2014(panel, baseline_year=2014):
    d = panel[panel["year"] >= baseline_year].copy()
    baseline = d[d["year"] == baseline_year][
        ["NUTS_ID", "log10_radiance"]
    ].rename(columns={"log10_radiance": "baseline_log10_radiance"})
    d = d.merge(baseline, on="NUTS_ID", how="inner", validate="many_to_one")
    d["delta_log10"] = d["log10_radiance"] - d["baseline_log10_radiance"]

    rows = []
    for (band, year), g in d.groupby(["density_band", "year"], observed=True):
        vals = g["delta_log10"].to_numpy(float)
        mean_log = float(vals.mean())
        se = float(vals.std(ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else np.nan
        rows.append({
            "density_band": str(band),
            "year": int(year),
            "n_nuts3": int(g["NUTS_ID"].nunique()),
            "mean_delta_log10": mean_log,
            "se_delta_log10": se,
            "observed_index_2014_100": 100.0 * 10.0 ** mean_log,
            "observed_index_ci_low": (
                100.0 * 10.0 ** (mean_log - 1.96 * se)
                if np.isfinite(se) else np.nan
            ),
            "observed_index_ci_high": (
                100.0 * 10.0 ** (mean_log + 1.96 * se)
                if np.isfinite(se) else np.nan
            ),
        })
    return pd.DataFrame(rows)


def plot_trend_comparison(t13, t14, config, output):
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    x = np.arange(len(labels), dtype=float)
    fig, ax = plt.subplots(figsize=(9.5, 5.8))

    for trends, xpos, label, marker in [
        (t13, x - 0.13, "2013–2020", "o"),
        (t14, x + 0.13, "2014–2020", "s"),
    ]:
        g = trends.set_index("density_band").loc[labels].reset_index()
        est = g["pct_radiance_change_per_year"].to_numpy(float)
        lo = g["pct_change_per_year_ci_low"].to_numpy(float)
        hi = g["pct_change_per_year_ci_high"].to_numpy(float)
        ax.errorbar(
            xpos, est, yerr=[est-lo, hi-est], fmt=marker,
            capsize=4, linewidth=1.2, label=label
        )

    ax.axhline(0, linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Fixed population-density class (people km$^{-2}$)")
    ax.set_ylabel(
        "Country-adjusted change in calibrated VIIRS DNB radiance (% yr$^{-1}$)"
    )
    ax.set_title(
        "Sensitivity of density-dependent VIIRS trends to excluding 2013"
    )
    ax.grid(axis="y", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_trajectories(traj, t14, config, output):
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    trends = t14.set_index("density_band")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0), sharex=True)
    axes = axes.ravel()

    for i, band in enumerate(labels):
        ax = axes[i]
        g = traj[traj["density_band"] == str(band)].sort_values("year")
        years = g["year"].to_numpy(float)
        rate = float(trends.loc[str(band), "pct_radiance_change_per_year"])
        ax.plot(
            years, g["observed_index_2014_100"],
            marker="o", linewidth=1.5, label="Observed data"
        )
        ax.fill_between(
            years, g["observed_index_ci_low"], g["observed_index_ci_high"],
            alpha=0.18
        )
        fitted = 100.0 * (1.0 + rate / 100.0) ** (years - 2014)
        ax.plot(
            years, fitted, linestyle="--", linewidth=1.8,
            label="2014–2020 fitted trend"
        )
        ax.axhline(100, linewidth=0.7)
        ax.set_title(
            f"{band} people km$^{{-2}}$\n{rate:+.2f}% yr$^{{-1}}$"
        )
        ax.grid(alpha=0.2)

    axes[5].axis("off")
    handles, labs = axes[0].get_legend_handles_labels()
    axes[5].legend(handles, labs, loc="center", frameon=False)
    fig.supylabel("Calibrated VIIRS DNB radiance index (2014 = 100)")
    fig.supxlabel("Year")
    fig.suptitle(
        "Historical radiance trajectories after excluding 2013\n"
        "Fixed NUTS2016 complete panel, 2014–2020",
        fontsize=14,
    )
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.94))
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_maps(s13, s14, config, output):
    geometry, europe_box = load_geometry(config)
    id_field = config["nuts"]["id_field"]

    a = s13[["NUTS_ID", "observed_pct_change_per_year"]].rename(
        columns={"observed_pct_change_per_year": "trend_2013_2020"}
    )
    b = s14[["NUTS_ID", "observed_pct_change_per_year"]].rename(
        columns={"observed_pct_change_per_year": "trend_2014_2020"}
    )
    vals = a.merge(b, on="NUTS_ID", validate="one_to_one")
    vals["effect_excluding_2013"] = (
        vals["trend_2014_2020"] - vals["trend_2013_2020"]
    )

    mapped = geometry.merge(
        vals, left_on=id_field, right_on="NUTS_ID",
        how="left", validate="one_to_one"
    ).to_crs("EPSG:3035")

    box3035 = europe_box.to_crs("EPSG:3035")
    xmin, ymin, xmax, ymax = box3035.total_bounds
    q = float(
        config["analysis"]["future_radiance_historical_validation"][
            "map_absolute_quantile"
        ]
    )

    common = np.concatenate([
        vals["trend_2013_2020"].to_numpy(float),
        vals["trend_2014_2020"].to_numpy(float),
    ])
    common = common[np.isfinite(common)]
    lim = max(float(np.quantile(np.abs(common), q)), 0.25)

    diff = vals["effect_excluding_2013"].to_numpy(float)
    diff = diff[np.isfinite(diff)]
    diff_lim = max(float(np.quantile(np.abs(diff), q)), 0.25)

    cmap = "coolwarm"
    trend_norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
    diff_norm = TwoSlopeNorm(vmin=-diff_lim, vcenter=0, vmax=diff_lim)

    panels = [
        ("trend_2013_2020", "Observed trend\n2013–2020", trend_norm),
        ("trend_2014_2020", "Observed trend\n2014–2020", trend_norm),
        (
            "effect_excluding_2013",
            "Effect of excluding 2013\n(2014–2020 minus 2013–2020)",
            diff_norm,
        ),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(17, 7))
    for ax, (col, title, norm) in zip(axes, panels):
        mapped.plot(ax=ax, color="0.92", edgecolor="0.78", linewidth=0.15)
        good = mapped[np.isfinite(mapped[col])].copy()
        good.plot(
            column=col, ax=ax, cmap=cmap, norm=norm,
            edgecolor="white", linewidth=0.08
        )
        mapped.dissolve(by="CNTR_CODE").boundary.plot(
            ax=ax, color="0.18", linewidth=0.35
        )
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_axis_off()
        ax.set_title(title)

    sm = ScalarMappable(norm=trend_norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(
        sm, ax=axes[:2].tolist(), orientation="horizontal",
        fraction=0.045, pad=0.035
    )
    cb.set_label(
        "Observed annual change in calibrated VIIRS DNB radiance (% yr$^{-1}$)"
    )

    sm2 = ScalarMappable(norm=diff_norm, cmap=cmap)
    sm2.set_array([])
    cb2 = fig.colorbar(
        sm2, ax=[axes[2]], orientation="horizontal",
        fraction=0.045, pad=0.035
    )
    cb2.set_label(
        "Change after excluding 2013 (percentage points yr$^{-1}$)"
    )

    fig.suptitle(
        "Spatial sensitivity of historical VIIRS trends to excluding 2013",
        fontsize=15, y=0.97
    )
    fig.savefig(
        output,
        dpi=int(
            config["analysis"]["future_radiance_historical_validation"][
                "figure_dpi"
            ]
        ),
        bbox_inches="tight",
    )
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()
    config = load_config(args.config)

    panel = load_fixed_panel(config)

    t13, info13 = fit_window(panel, 2013, 2020, config)
    t14, info14 = fit_window(panel, 2014, 2020, config)
    s13 = region_slopes(panel, 2013, 2020)
    s14 = region_slopes(panel, 2014, 2020)
    traj14 = trajectories_2014(panel, 2014)

    a = t13[
        ["density_band", "pct_radiance_change_per_year",
         "pct_change_per_year_ci_low", "pct_change_per_year_ci_high", "trend_p"]
    ].rename(columns={
        "pct_radiance_change_per_year": "pct_per_year_2013_2020",
        "pct_change_per_year_ci_low": "ci_low_2013_2020",
        "pct_change_per_year_ci_high": "ci_high_2013_2020",
        "trend_p": "p_2013_2020",
    })
    b = t14[
        ["density_band", "pct_radiance_change_per_year",
         "pct_change_per_year_ci_low", "pct_change_per_year_ci_high", "trend_p"]
    ].rename(columns={
        "pct_radiance_change_per_year": "pct_per_year_2014_2020",
        "pct_change_per_year_ci_low": "ci_low_2014_2020",
        "pct_change_per_year_ci_high": "ci_high_2014_2020",
        "trend_p": "p_2014_2020",
    })
    comparison = a.merge(b, on="density_band", validate="one_to_one")
    comparison["change_after_excluding_2013_pp_per_year"] = (
        comparison["pct_per_year_2014_2020"]
        - comparison["pct_per_year_2013_2020"]
    )

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "future_radiance_historical_validation"
        / "exclude_2013"
    )
    figures.mkdir(parents=True, exist_ok=True)

    comparison_path = processed / "future_radiance_no2013_density_band_comparison.csv"
    slopes_path = processed / "future_radiance_no2013_region_slopes.csv"
    summary_path = processed / "future_radiance_no2013_sensitivity_summary.txt"
    comparison.to_csv(comparison_path, index=False)
    s14.to_csv(slopes_path, index=False)

    fig_compare = figures / "density_trends_2013_2020_vs_2014_2020.pdf"
    fig_traj = figures / "real_trajectories_2014_2020.pdf"
    fig_maps = figures / "observed_maps_2013_2020_vs_2014_2020.pdf"
    plot_trend_comparison(t13, t14, config, fig_compare)
    plot_trajectories(traj14, t14, config, fig_traj)
    plot_maps(s13, s14, config, fig_maps)

    with summary_path.open("w", encoding="utf-8") as h:
        h.write(
            "SENSITIVITY OF POPULATION-DENSITY VIIRS TRENDS TO EXCLUDING 2013\n"
        )
        h.write("=" * 80 + "\n\n")
        h.write(
            "Same fixed NUTS2016 complete-panel regions and same fixed "
            "population-density classes are used in both analyses. The only "
            "change is removal of the 2013 radiance observations.\n\n"
        )
        h.write("COUNTRY-ADJUSTED MIXED-MODEL TRENDS\n")
        h.write(
            comparison.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\nDENSITY-BAND × YEAR TESTS\n")
        h.write(
            f"2013-2020: LR={info13['lr']:.6g}, "
            f"df={info13['df']}, p={info13['p']:.6g}\n"
        )
        h.write(
            f"2014-2020: LR={info14['lr']:.6g}, "
            f"df={info14['df']}, p={info14['p']:.6g}\n"
        )

    print("\n2013 EXCLUSION SENSITIVITY")
    print(
        comparison.to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print("\nWROTE")
    for p in [comparison_path, slopes_path, summary_path,
              fig_compare, fig_traj, fig_maps]:
        print(p)


if __name__ == "__main__":
    main()
