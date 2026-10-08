#!/usr/bin/env python3
"""Create production country facets: six alphabetical countries per PDF page.

Each country panel contains all available years of corrected NUTS3 radiance
versus population density in log10-log10 space, with annual fitted lines.
"""

from __future__ import annotations

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from project_config import configured_path, load_config


def fit_line(x: np.ndarray, y: np.ndarray):
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if len(x) < 2:
        return None
    return np.polyfit(x, y, 1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--countries-per-page", type=int, default=6)
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")
    figures = configured_path(config, "figures") / "population_relationship"
    figures.mkdir(parents=True, exist_ok=True)

    data = pd.read_parquet(
        processed / config["outputs"]["merged"]
    ).copy()

    xcol = "log10_population_density"
    ycol = "log10_radiance_mean_corrected"

    data = data[
        np.isfinite(data[xcol]) & np.isfinite(data[ycol])
    ].copy()

    countries = sorted(data["CNTR_CODE"].dropna().unique().tolist())
    years = sorted(data["year"].dropna().astype(int).unique().tolist())

    ncols = 2
    nrows = 3
    if args.countries_per_page != 6:
        nrows = math.ceil(args.countries_per_page / ncols)

    # Use fixed global axes so country panels are directly comparable.
    xmin, xmax = data[xcol].min(), data[xcol].max()
    ymin, ymax = data[ycol].min(), data[ycol].max()
    xpad = 0.03 * (xmax - xmin)
    ypad = 0.03 * (ymax - ymin)

    cmap = plt.get_cmap("viridis")
    colours = {
        year: cmap(i / max(1, len(years) - 1))
        for i, year in enumerate(years)
    }

    output = figures / "country_population_radiance_facets_all_years.pdf"

    with PdfPages(output) as pdf:
        for start in range(0, len(countries), args.countries_per_page):
            page_countries = countries[start:start + args.countries_per_page]

            fig, axes = plt.subplots(
                nrows,
                ncols,
                figsize=(11.69, 16.54),
                squeeze=False,
                sharex=True,
                sharey=True,
            )
            axes = axes.ravel()
            legend_handles = {}

            for ax, country in zip(axes, page_countries):
                c = data[data["CNTR_CODE"] == country].copy()

                for year in years:
                    g = c[c["year"] == year]
                    if g.empty:
                        continue

                    x = g[xcol].to_numpy(dtype=float)
                    y = g[ycol].to_numpy(dtype=float)

                    artist = ax.scatter(
                        x,
                        y,
                        s=11,
                        alpha=0.28,
                        color=colours[year],
                        label=str(year),
                    )
                    legend_handles.setdefault(year, artist)

                    fit = fit_line(x, y)
                    if fit is not None:
                        slope, intercept = fit
                        xs = np.linspace(x.min(), x.max(), 100)
                        ax.plot(
                            xs,
                            intercept + slope * xs,
                            color=colours[year],
                            linewidth=1.45,
                            alpha=0.95,
                        )

                n_regions = c["NUTS_ID"].nunique()
                ax.set_title(
                    f"{country}   NUTS3={n_regions}",
                    fontsize=11,
                    fontweight="bold",
                )
                ax.set_xlim(xmin - xpad, xmax + xpad)
                ax.set_ylim(ymin - ypad, ymax + ypad)
                ax.grid(alpha=0.16)

            for ax in axes[len(page_countries):]:
                ax.axis("off")

            ordered = [year for year in years if year in legend_handles]
            fig.legend(
                [legend_handles[y] for y in ordered],
                [str(y) for y in ordered],
                loc="lower center",
                ncol=6,
                frameon=False,
                bbox_to_anchor=(0.5, 0.018),
                title="Year",
            )

            fig.supxlabel(
                "log10 population density (people km⁻²)",
                fontsize=12,
                y=0.045,
            )
            fig.supylabel(
                "log10 corrected mean VIIRS radiance (nW cm⁻² sr⁻¹)",
                fontsize=12,
                x=0.025,
            )
            fig.suptitle(
                "NUTS3 population density vs corrected VIIRS radiance",
                fontsize=15,
                y=0.985,
            )
            fig.tight_layout(rect=[0.045, 0.065, 0.99, 0.965])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"Countries plotted: {len(countries)}")
    print("Order:", ", ".join(countries))
    print(f"Wrote: {output}")


if __name__ == "__main__":
    main()
