#!/usr/bin/env python3
"""Fit and plot annual country-level log radiance vs log population density.

Each available NUTS3 pair is treated as an independent observation within a
country/year. Fits are ordinary least-squares straight lines in log10-log10
space:

    log10(mean VIIRS radiance) ~ log10(population density)

Rows with non-positive radiance or population density are retained in the
merged data but cannot enter a logarithmic fit, so they are counted and
reported separately.
"""

from __future__ import annotations

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_config import configured_path, load_config


def fit_one(group: pd.DataFrame) -> dict:
    raw_n = len(group)
    fit = group[
        (group["radiance_mean"] > 0)
        & (group["population_density"] > 0)
    ].copy()

    n = len(fit)
    if n < 2:
        return {
            "n_pairs": raw_n,
            "n_fit": n,
            "n_excluded_nonpositive": raw_n - n,
            "slope": np.nan,
            "intercept": np.nan,
            "r2": np.nan,
        }

    x = np.log10(fit["population_density"].to_numpy(dtype=float))
    y = np.log10(fit["radiance_mean"].to_numpy(dtype=float))

    slope, intercept = np.polyfit(x, y, 1)
    predicted = intercept + slope * x

    ss_res = np.sum((y - predicted) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {
        "n_pairs": raw_n,
        "n_fit": n,
        "n_excluded_nonpositive": raw_n - n,
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
    }


def build_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (year, country), group in data.groupby(["year", "CNTR_CODE"], sort=True):
        row = {"year": int(year), "CNTR_CODE": country}
        row.update(fit_one(group))
        rows.append(row)

    return pd.DataFrame(rows).sort_values(["year", "CNTR_CODE"]).reset_index(drop=True)


def plot_countries(
    data: pd.DataFrame,
    summary: pd.DataFrame,
    year: int,
    countries: list[str],
    output,
) -> None:
    subset = data[
        (data["year"] == year)
        & data["CNTR_CODE"].isin(countries)
    ].copy()

    present = [
        country
        for country in countries
        if country in set(subset["CNTR_CODE"])
    ]

    if not present:
        raise ValueError(f"None of the requested countries are present for {year}.")

    ncols = 3
    nrows = math.ceil(len(present) / ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(14, 4.4 * nrows),
        squeeze=False,
        constrained_layout=True,
    )

    for ax, country in zip(axes.flat, present):
        group = subset[subset["CNTR_CODE"] == country].copy()
        group = group[
            (group["radiance_mean"] > 0)
            & (group["population_density"] > 0)
        ]

        x = np.log10(group["population_density"].to_numpy(dtype=float))
        y = np.log10(group["radiance_mean"].to_numpy(dtype=float))

        ax.scatter(x, y, s=22, alpha=0.7)

        stats = summary[
            (summary["year"] == year)
            & (summary["CNTR_CODE"] == country)
        ].iloc[0]

        if len(group) >= 2 and np.isfinite(stats["slope"]):
            xline = np.linspace(x.min(), x.max(), 200)
            yline = stats["intercept"] + stats["slope"] * xline
            ax.plot(xline, yline, linewidth=2)

        ax.set_title(
            f"{country}  n={int(stats['n_fit'])}  "
            f"slope={stats['slope']:.3f}  R²={stats['r2']:.3f}"
        )
        ax.set_xlabel("log10 population density (people km⁻²)")
        ax.set_ylabel("log10 mean VIIRS radiance (nW cm⁻² sr⁻¹)")
        ax.grid(alpha=0.2)

    for ax in axes.flat[len(present):]:
        ax.remove()

    fig.suptitle(
        f"NUTS3 population density vs VIIRS radiance, {year}",
        fontsize=15,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--year", type=int, default=None)
    parser.add_argument(
        "--countries",
        nargs="+",
        default=None,
        help="Country codes for the comparison figure.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")
    figures = configured_path(config, "figures")

    merged_path = processed / config["outputs"]["merged"]
    data = pd.read_parquet(merged_path)

    summary = build_summary(data)

    summary_csv = processed / "country_loglog_fits.csv"
    summary.to_csv(summary_csv, index=False)

    year = args.year or int(data["year"].max())
    countries = (
        args.countries
        if args.countries
        else config["analysis"]["comparison_countries"]
    )

    output = (
        figures
        / "population_relationship"
        / f"population_radiance_loglog_{year}.pdf"
    )

    plot_countries(data, summary, year, countries, output)

    print(f"Fit summary for {year}:")
    display = summary[summary["year"] == year].copy()
    print(display.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    print()
    print(f"Wrote: {summary_csv}")
    print(f"Wrote: {output}")


if __name__ == "__main__":
    main()
