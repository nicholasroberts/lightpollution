#!/usr/bin/env python3
"""Publication-quality annual maps of corrected NUTS3 VIIRS radiance across Europe.

Writes one high-resolution PNG per analysis year using:
- production-corrected mean radiance only;
- the Eurostat-aligned NUTS release configured for each year;
- ETRS89 / LAEA Europe (EPSG:3035);
- a common logarithmic colour scale across all years;
- a mainland-Europe plotting extent that excludes overseas territories.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm
from shapely.geometry import box

from project_config import configured_path, load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")
    nuts_raw = configured_path(config, "nuts_raw")
    figures = (
        configured_path(config, "figures")
        / "annual_radiance"
        / "maps"
    )
    figures.mkdir(parents=True, exist_ok=True)

    input_path = (
        processed / config["outputs"]["radiance_by_nuts3_calibrated"]
    )
    data = pd.read_parquet(input_path).copy()

    value_col = "radiance_mean_corrected"
    if value_col not in data.columns:
        raise ValueError(
            "Production Europe maps require 'radiance_mean_corrected'. "
            "Run the calibration stage first."
        )

    first_year = int(config["analysis"]["first_year"])
    last_year = int(config["analysis"]["last_year"])
    years = list(range(first_year, last_year + 1))

    data = data[data["year"].isin(years)].copy()

    positive = data[
        np.isfinite(data[value_col]) & (data[value_col] > 0)
    ].copy()
    if positive.empty:
        raise ValueError("No positive corrected radiance values available.")

    map_cfg = config["analysis"]["europe_maps"]
    lower_q = float(map_cfg["colour_min_quantile"])
    upper_q = float(map_cfg["colour_max_quantile"])

    vmin = float(positive[value_col].quantile(lower_q))
    vmax = float(positive[value_col].quantile(upper_q))

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin <= 0 or vmax <= vmin:
        raise ValueError(
            f"Invalid common colour limits: vmin={vmin}, vmax={vmax}"
        )

    west, south, east, north = map(
        float, map_cfg["extent_4326"]
    )
    europe_box = gpd.GeoDataFrame(
        geometry=[box(west, south, east, north)],
        crs="EPSG:4326",
    )
    europe_box_3035 = europe_box.to_crs("EPSG:3035")
    xmin, ymin, xmax, ymax = europe_box_3035.total_bounds

    cmap = map_cfg.get("cmap", "inferno")
    norm = LogNorm(vmin=vmin, vmax=vmax)

    scale_table = pd.DataFrame(
        [
            {
                "radiance_field": value_col,
                "vmin": vmin,
                "vmax": vmax,
                "lower_quantile": lower_q,
                "upper_quantile": upper_q,
                "projection": "EPSG:3035",
                "extent_west": west,
                "extent_south": south,
                "extent_east": east,
                "extent_north": north,
                "dpi": args.dpi,
            }
        ]
    )
    scale_path = figures / "europe_nuts3_radiance_map_scale.csv"
    scale_table.to_csv(scale_path, index=False)

    id_field = config["nuts"]["id_field"]
    filename_template = config["nuts"]["filename_template"]
    release_by_year = {
        int(year): int(release)
        for year, release in config["nuts"]["release_by_year"].items()
    }

    for year in years:
        if year not in release_by_year:
            raise KeyError(f"No configured NUTS release for {year}")

        release = release_by_year[year]
        nuts_path = nuts_raw / filename_template.format(release=release)
        if not nuts_path.exists():
            raise FileNotFoundError(nuts_path)

        geometry = gpd.read_file(nuts_path)

        # Drop overseas geometries before projection. Intersections at the
        # European plotting extent are retained.
        geometry = gpd.clip(geometry, europe_box)

        annual = data[data["year"] == year][
            [id_field, value_col]
        ].copy()

        mapped = geometry.merge(
            annual,
            on=id_field,
            how="left",
            validate="one_to_one",
        ).to_crs("EPSG:3035")

        fig, ax = plt.subplots(figsize=(10.6, 9.0))

        # Base layer for NUTS3 regions without a usable positive value.
        mapped.plot(
            ax=ax,
            color="0.93",
            edgecolor="0.78",
            linewidth=0.18,
        )

        positive_map = mapped[
            np.isfinite(mapped[value_col]) & (mapped[value_col] > 0)
        ].copy()

        if not positive_map.empty:
            positive_map.plot(
                column=value_col,
                ax=ax,
                cmap=cmap,
                norm=norm,
                edgecolor="white",
                linewidth=0.12,
            )

        # Stronger national boundaries improve readability without obscuring
        # NUTS3 structure.
        country_outline = mapped.dissolve(by="CNTR_CODE")
        country_outline.boundary.plot(
            ax=ax,
            color="0.15",
            linewidth=0.45,
        )

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_axis_off()

        ax.set_title(
            f"Corrected mean VIIRS night-time radiance by NUTS3, {year}",
            fontsize=16,
            pad=12,
        )

        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(
            sm,
            ax=ax,
            fraction=0.035,
            pad=0.02,
        )
        cbar.set_label(
            "Corrected radiance (nW cm$^{-2}$ sr$^{-1}$)",
            fontsize=11,
        )

        fig.text(
            0.5,
            0.018,
            "Common logarithmic colour scale across 2013–2024; "
            "overseas territories excluded by European map extent.",
            ha="center",
            fontsize=9,
        )

        output = figures / f"europe_nuts3_corrected_radiance_{year}.png"
        fig.savefig(
            output,
            dpi=args.dpi,
            bbox_inches="tight",
            facecolor="white",
        )
        plt.close(fig)

        n_positive = int(len(positive_map))
        print(
            f"{year}: NUTS {release}; mapped positive regions={n_positive:,}; "
            f"wrote {output}"
        )

    print()
    print(
        f"Common corrected-radiance scale: {vmin:.6g} to {vmax:.6g} "
        f"nW cm^-2 sr^-1"
    )
    print(f"Wrote scale metadata: {scale_path}")


if __name__ == "__main__":
    main()
