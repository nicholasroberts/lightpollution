#!/usr/bin/env python3
"""Plot country-level fit statistics through time."""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import pandas as pd

from project_config import configured_path, load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument(
        "--countries",
        nargs="+",
        default=None,
        help="Country codes to plot. Defaults to analysis.comparison_countries.",
    )
    parser.add_argument(
        "--all-countries",
        action="store_true",
        help="Plot every country present in country_loglog_fits.csv.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")
    figures = configured_path(config, "figures")

    path = processed / "country_loglog_fits.csv"
    df = pd.read_csv(path)

    if args.all_countries:
        countries = sorted(df["CNTR_CODE"].dropna().unique().tolist())
    else:
        countries = (
            args.countries
            if args.countries
            else config["analysis"]["comparison_countries"]
        )
    d = df[df["CNTR_CODE"].isin(countries)].copy()

    outdir = figures / "population_relationship"
    outdir.mkdir(parents=True, exist_ok=True)

    suffix = "_all_countries" if args.all_countries else ""

    # Slope through time.
    fig_width = 14 if args.all_countries else 10
    fig, ax = plt.subplots(figsize=(fig_width, 7 if args.all_countries else 6))
    cmap = plt.get_cmap("turbo")
    colour = {
        country: cmap(i / max(1, len(countries) - 1))
        for i, country in enumerate(countries)
    }
    for country in countries:
        g = d[d["CNTR_CODE"] == country].sort_values("year")
        if g.empty:
            continue
        ax.plot(
            g["year"],
            g["slope"],
            marker="o",
            markersize=3.5 if args.all_countries else 6,
            linewidth=1.2 if args.all_countries else 1.5,
            color=colour[country],
            label=country,
        )
    ax.axvline(2016.5, linestyle="--", linewidth=1.2)
    ax.text(
        2016.55,
        ax.get_ylim()[1],
        "2017 calibration change",
        va="top",
        ha="left",
        fontsize=9,
    )
    ax.set_xlabel("Year")
    ax.set_ylabel("Slope of log10(radiance) ~ log10(population density)")
    ax.set_title("Country log-log slope through time")
    ax.grid(alpha=0.2)
    if args.all_countries:
        ax.legend(
            ncol=4,
            fontsize=8,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            borderaxespad=0,
        )
        fig.tight_layout(rect=[0, 0, 0.82, 1])
    else:
        ax.legend(ncol=2)
        fig.tight_layout()
    slope_path = outdir / f"country_loglog_slope_by_year{suffix}.pdf"
    fig.savefig(slope_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Intercept through time.
    fig, ax = plt.subplots(figsize=(fig_width, 7 if args.all_countries else 6))
    for country in countries:
        g = d[d["CNTR_CODE"] == country].sort_values("year")
        if g.empty:
            continue
        ax.plot(
            g["year"],
            g["intercept"],
            marker="o",
            markersize=3.5 if args.all_countries else 6,
            linewidth=1.2 if args.all_countries else 1.5,
            color=colour[country],
            label=country,
        )
    ax.axvline(2016.5, linestyle="--", linewidth=1.2)
    ax.text(
        2016.55,
        ax.get_ylim()[1],
        "2017 calibration change",
        va="top",
        ha="left",
        fontsize=9,
    )
    ax.set_xlabel("Year")
    ax.set_ylabel("Intercept")
    ax.set_title("Country log-log intercept through time")
    ax.grid(alpha=0.2)
    if args.all_countries:
        ax.legend(
            ncol=4,
            fontsize=8,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            borderaxespad=0,
        )
        fig.tight_layout(rect=[0, 0, 0.82, 1])
    else:
        ax.legend(ncol=2)
        fig.tight_layout()
    intercept_path = outdir / f"country_loglog_intercept_by_year{suffix}.pdf"
    fig.savefig(intercept_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Fit quality through time.
    fig, ax = plt.subplots(figsize=(fig_width, 7 if args.all_countries else 6))
    for country in countries:
        g = d[d["CNTR_CODE"] == country].sort_values("year")
        if g.empty:
            continue
        ax.plot(
            g["year"],
            g["r2"],
            marker="o",
            markersize=3.5 if args.all_countries else 6,
            linewidth=1.2 if args.all_countries else 1.5,
            color=colour[country],
            label=country,
        )
    ax.axvline(2016.5, linestyle="--", linewidth=1.2)
    ax.set_xlabel("Year")
    ax.set_ylabel("R²")
    ax.set_title("Country log-log fit quality through time")
    ax.grid(alpha=0.2)
    if args.all_countries:
        ax.legend(
            ncol=4,
            fontsize=8,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            borderaxespad=0,
        )
        fig.tight_layout(rect=[0, 0, 0.82, 1])
    else:
        ax.legend(ncol=2)
        fig.tight_layout()
    r2_path = outdir / f"country_loglog_r2_by_year{suffix}.pdf"
    fig.savefig(r2_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Year-to-year changes, useful for detecting step changes.
    changes = d.sort_values(["CNTR_CODE", "year"]).copy()
    changes["delta_slope"] = changes.groupby("CNTR_CODE")["slope"].diff()
    changes["delta_intercept"] = changes.groupby("CNTR_CODE")["intercept"].diff()
    changes["delta_r2"] = changes.groupby("CNTR_CODE")["r2"].diff()

    changes_path = processed / f"country_loglog_fit_year_changes{suffix}.csv"
    changes.to_csv(changes_path, index=False)

    focus = changes[changes["year"].isin([2016, 2017, 2018])]
    print(
        focus[
            [
                "year",
                "CNTR_CODE",
                "n_fit",
                "slope",
                "delta_slope",
                "intercept",
                "delta_intercept",
                "r2",
            ]
        ]
        .sort_values(["CNTR_CODE", "year"])
        .to_string(index=False, float_format=lambda x: f"{x:.4f}")
    )
    print()
    print(f"Wrote: {slope_path}")
    print(f"Wrote: {intercept_path}")
    print(f"Wrote: {r2_path}")
    print(f"Wrote: {changes_path}")


if __name__ == "__main__":
    main()
