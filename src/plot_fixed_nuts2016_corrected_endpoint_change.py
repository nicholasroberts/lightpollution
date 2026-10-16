#!/usr/bin/env python3
"""Observed 2014-2024 corrected VIIRS endpoint changes on fixed NUTS2016.

A descriptive, model-free diagnostic. Uses the exact corrected annual input,
balanced NUTS3 cohort and fixed 2014-2020 population-density classes already
used by the corrected minimum temporal-model stage. No refitting, raster
extraction, cohort reselection, or further correction is performed.

One record per NUTS3:
    absolute_change = radiance_2024 - radiance_2014
    percent_change = 100 * (radiance_2024 / radiance_2014 - 1)

All units for absolute radiance changes are nW cm^-2 sr^-1. The map colour
limits are robust display limits ONLY: all values remain in the paired CSV.
The boxplot full-range panels display all points; the companion zoom panels
use the same unmodified data but restrict their visible y axes.
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
from shapely.geometry import box

from project_config import configured_path, load_config
from select_fixed_nuts2016_corrected_temporal_model import load_corrected_panel
from validate_fixed_nuts2016_2014_2024 import load_fixed_nuts


STAGE = "fixed_nuts2016_corrected_endpoint_change"


def paired_endpoint_changes(panel: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Retain the existing locked cohort; compare each region to itself."""
    cfg = config["analysis"]["fixed_nuts2016_corrected_temporal_model"]
    first_year, last_year = int(cfg["first_year"]), int(cfg["last_year"])
    expected_years = last_year - first_year + 1

    if panel.duplicated(["NUTS_ID", "year"]).any():
        raise RuntimeError("Duplicate NUTS3/year rows in the corrected panel.")
    years_per_region = panel.groupby("NUTS_ID")["year"].nunique()
    if years_per_region.empty or not years_per_region.eq(expected_years).all():
        raise RuntimeError("Corrected panel is not balanced over 2014-2024.")

    common = ["NUTS_ID", "CNTR_CODE", "density_band"]
    a = panel.loc[
        panel["year"].eq(first_year),
        common + ["reference_population_density", "radiance_mean_corrected"],
    ].rename(columns={"radiance_mean_corrected": "radiance_2014_nw"})
    b = panel.loc[
        panel["year"].eq(last_year),
        common + ["radiance_mean_corrected"],
    ].rename(columns={"radiance_mean_corrected": "radiance_2024_nw"})
    paired = a.merge(b, on=common, how="inner", validate="one_to_one")

    expected_n = int(
        config["analysis"]["fixed_nuts2016_2021_diagnostic"]["expected_n_nuts3"]
    )
    if len(paired) != expected_n or panel["NUTS_ID"].nunique() != expected_n:
        raise RuntimeError(
            f"Expected the locked {expected_n}-region cohort; found "
            f"{len(paired)} paired endpoints and "
            f"{panel['NUTS_ID'].nunique()} panel regions."
        )

    first = paired["radiance_2014_nw"].to_numpy(dtype=float)
    last = paired["radiance_2024_nw"].to_numpy(dtype=float)
    if (not np.isfinite(first).all() or not np.isfinite(last).all()
            or np.any(first <= 0) or np.any(last <= 0)):
        raise RuntimeError("Non-positive/non-finite endpoint radiance in cohort.")

    paired["absolute_change_nw"] = last - first
    paired["percent_change"] = 100.0 * (last / first - 1.0)
    if not np.isfinite(paired[["absolute_change_nw", "percent_change"]]).all().all():
        raise RuntimeError("Non-finite endpoint change detected.")

    return paired.sort_values(["CNTR_CODE", "NUTS_ID"]).reset_index(drop=True)


def make_summary(paired: pd.DataFrame, labels: list[str]) -> pd.DataFrame:
    """Descriptive regional distribution (one region = one observation)."""
    rows = []
    for band in labels:
        g = paired.loc[paired["density_band"].astype(str).eq(str(band))]
        if g.empty:
            raise RuntimeError(f"No regions in fixed density class {band}.")
        entry = {"density_band": str(band), "n_nuts3": len(g)}
        for metric, prefix in (
            ("absolute_change_nw", "absolute_nw"),
            ("percent_change", "percent"),
        ):
            vals = g[metric].to_numpy(dtype=float)
            entry.update({
                f"{prefix}_mean": float(np.mean(vals)),
                f"{prefix}_median": float(np.median(vals)),
                f"{prefix}_q25": float(np.quantile(vals, 0.25)),
                f"{prefix}_q75": float(np.quantile(vals, 0.75)),
                f"{prefix}_min": float(np.min(vals)),
                f"{prefix}_max": float(np.max(vals)),
            })
        rows.append(entry)
    return pd.DataFrame(rows)


def save_figure(fig, stem: Path, dpi: int) -> None:
    for extension in ("pdf", "png"):
        fig.savefig(
            stem.with_suffix("." + extension),
            dpi=dpi, bbox_inches="tight", facecolor="white",
        )
    plt.close(fig)


def european_map_geometry(config: dict) -> tuple[gpd.GeoDataFrame, tuple]:
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    geometry = load_fixed_nuts(config, release)
    west, south, east, north = map(
        float, config["analysis"]["europe_maps"]["extent_4326"]
    )
    envelope = gpd.GeoDataFrame(
        geometry=[box(west, south, east, north)], crs="EPSG:4326",
    )
    geometry = gpd.clip(geometry, envelope).to_crs("EPSG:3035")
    bounds = tuple(envelope.to_crs("EPSG:3035").total_bounds)
    return geometry[["NUTS_ID", "CNTR_CODE", "geometry"]], bounds


def plot_map(
    geometry: gpd.GeoDataFrame,
    bounds: tuple,
    paired: pd.DataFrame,
    metric: str,
    label: str,
    title: str,
    stem: Path,
    dpi: int,
    quantile: float = 0.98,
) -> None:
    values = paired[metric].to_numpy(dtype=float)
    limit = float(np.quantile(np.abs(values), quantile))
    if not np.isfinite(limit) or limit <= 0:
        limit = 1.0
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    cmap = "RdBu_r"

    mapped = geometry.merge(
        paired[["NUTS_ID", metric]], on="NUTS_ID",
        how="left", validate="one_to_one",
    )
    present = mapped[metric].notna()
    fig, ax = plt.subplots(figsize=(10.6, 9.0))
    # Grey polygons are explicitly outside the locked cohort.
    mapped.plot(
        ax=ax, color="0.91", edgecolor="0.79", linewidth=0.14,
    )
    if present.any():
        mapped.loc[present].plot(
            ax=ax, column=metric, cmap=cmap, norm=norm,
            edgecolor="white", linewidth=0.11,
        )

    mapped.dissolve(by="CNTR_CODE").boundary.plot(
        ax=ax, color="0.28", linewidth=0.42,
    )
    xmin, ymin, xmax, ymax = bounds
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title(title, fontsize=15, pad=14)

    scalar = ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    cbar = fig.colorbar(
        scalar, ax=ax, fraction=0.035, pad=0.018, extend="both",
    )
    cbar.set_label(label, fontsize=11)
    fig.text(
        0.5, 0.025,
        f"Observed paired NUTS2016 endpoints; n={len(paired):,} regions. "
        f"Colour saturation beyond ±{quantile:.0%} quantile of |change|; "
        "grey = region outside cohort.",
        ha="center", va="bottom", fontsize=8.6,
    )
    save_figure(fig, stem, dpi)
    print(
        f"Map {metric}: {int(present.sum())} cohort polygons within extent; "
        f"colour limits ±{limit:.6g}"
    )


def plot_boxplots(
    paired: pd.DataFrame,
    labels: list[str],
    metric: str,
    ylabel: str,
    title: str,
    stem: Path,
    dpi: int,
) -> None:
    """Full and 1st-99th percentile zoom side-by-side; no data trimming."""
    arrays = [
        paired.loc[
            paired["density_band"].astype(str).eq(str(band)), metric
        ].to_numpy(dtype=float)
        for band in labels
    ]
    pooled = np.concatenate(arrays)
    q01, q99 = np.quantile(pooled, [0.01, 0.99])
    span = max(float(q99 - q01), 1e-8)
    zoom_low = min(0.0, float(q01)) - 0.06 * span
    zoom_high = max(0.0, float(q99)) + 0.06 * span
    counts = [len(x) for x in arrays]
    positions = np.arange(1, len(labels) + 1)

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.3))
    for ax, subtitle in zip(axes, ("Full distribution", "Central 1–99% view")):
        bp = ax.boxplot(
            arrays, positions=positions, widths=0.64,
            showfliers=True, patch_artist=True,
            medianprops={"color": "#a33b25", "linewidth": 1.8},
            flierprops={
                "marker": ".", "markersize": 2.4,
                "markerfacecolor": "#303e52", "markeredgecolor": "none",
                "alpha": 0.35,
            },
        )
        for patch in bp["boxes"]:
            patch.set_facecolor("#aac4dc")
            patch.set_alpha(0.60)
        ax.axhline(0, color="0.28", linewidth=0.9, zorder=0)
        ax.set_xticks(positions)
        ax.set_xticklabels(
            [f"{band}\n(n={n})" for band, n in zip(labels, counts)],
            fontsize=9,
        )
        ax.set_xlabel("Fixed population density (people km$^{-2}$)")
        ax.set_ylabel(ylabel)
        ax.set_title(subtitle)
        ax.grid(axis="y", alpha=0.2)
    axes[1].set_ylim(zoom_low, zoom_high)

    fig.suptitle(title, fontsize=14)
    fig.text(
        0.5, 0.01,
        "Each observation is a paired NUTS3 region; boxes show median and IQR, "
        "whiskers extend to 1.5×IQR, dots are outliers. "
        "The right panel changes display limits only.",
        ha="center", fontsize=8.8,
    )
    fig.tight_layout(rect=(0.01, 0.04, 1.0, 0.95))
    save_figure(fig, stem, dpi)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_corrected_temporal_model"]
    if (int(cfg["first_year"]), int(cfg["last_year"])) != (2014, 2024):
        raise RuntimeError("This diagnostic is specifically for 2014-2024.")

    panel, _ = load_corrected_panel(config, cfg)
    paired = paired_endpoint_changes(panel, config)
    labels = list(config["analysis"]["density_temporal_change"]["density_band_labels"])
    summary = make_summary(paired, labels)
    figures = (
        configured_path(config, "figures")
        / "population_relationship" / STAGE
    )
    processed = configured_path(config, "processed") / STAGE
    figures.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)
    dpi = int(config["analysis"]["fixed_nuts2016_2014_2024"]["figure_dpi"])

    paired_file = processed / "paired_nuts3_radiance_changes_2014_2024.csv"
    summary_file = processed / "density_class_endpoint_change_summary.csv"
    paired.to_csv(paired_file, index=False)
    summary.to_csv(summary_file, index=False)

    geometry, bounds = european_map_geometry(config)
    for metric, label, title, filename in (
        (
            "absolute_change_nw",
            "2024 minus 2014 VIIRS radiance (nW cm$^{-2}$ sr$^{-1}$)",
            "Observed absolute change in corrected VIIRS DNB radiance\n2014–2024, fixed NUTS2016",
            "map_absolute_radiance_change_2014_2024",
        ),
        (
            "percent_change",
            "Change in VIIRS radiance, 2014–2024 (%)",
            "Observed percentage change in corrected VIIRS DNB radiance\n2014–2024, fixed NUTS2016",
            "map_percent_radiance_change_2014_2024",
        ),
    ):
        plot_map(
            geometry, bounds, paired, metric, label, title,
            figures / filename, dpi,
        )

    for metric, label, title, filename in (
        (
            "absolute_change_nw",
            "2024 minus 2014 VIIRS radiance (nW cm$^{-2}$ sr$^{-1}$)",
            "Observed 2014–2024 absolute radiance change by population-density class",
            "boxplot_absolute_radiance_change_by_density_2014_2024",
        ),
        (
            "percent_change",
            "Observed 2014–2024 radiance change (%)",
            "Observed 2014–2024 percentage radiance change by population-density class",
            "boxplot_percent_radiance_change_by_density_2014_2024",
        ),
    ):
        plot_boxplots(paired, labels, metric, label, title, figures / filename, dpi)

    print(f"\nCORRECTED 2014–2024 ENDPOINT CHANGES: {len(paired)} NUTS3 regions")
    print(
        summary[
            ["density_band", "n_nuts3", "absolute_nw_median",
             "percent_median", "percent_q25", "percent_q75"]
        ].to_string(index=False, float_format=lambda x: f"{x:.5g}")
    )
    print("\nWROTE")
    print(paired_file)
    print(summary_file)
    for name in (
        "map_absolute_radiance_change_2014_2024",
        "map_percent_radiance_change_2014_2024",
        "boxplot_absolute_radiance_change_by_density_2014_2024",
        "boxplot_percent_radiance_change_by_density_2014_2024",
    ):
        print(figures / f"{name}.pdf")
        print(figures / f"{name}.png")


if __name__ == "__main__":
    main()
