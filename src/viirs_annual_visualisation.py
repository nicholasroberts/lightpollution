#!/usr/bin/env python3
"""
Modern VIIRS annual visualisation test.

Purpose
-------
Create publication-style year-by-year VIIRS night-time-light figures with:
  1. identical spatial extent for every year
  2. one shared radiance scale across all years
  3. a nonlinear display transform that preserves near-zero values
  4. optional NUTS3 boundaries
  5. a second figure showing change from the first year

This script is deliberately plotting-only. It does not yet calculate NUTS3
zonal statistics; that should be a separate pipeline stage.

Example
-------
python viirs_annual_visualisation.py \
    --data-dir VIIRSdata/annual_v22 \
    --nuts geojson/NUTS_RG_01M_2024_4326_LEVL_3.geojson \
    --pattern "*average-masked*.tif" \
    --out-dir figures
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from affine import Affine
from matplotlib.colors import SymLogNorm, TwoSlopeNorm
from rasterio.enums import Resampling
from rasterio.transform import array_bounds
from rasterio.windows import Window, from_bounds
from rasterio.warp import transform_bounds


YEAR_RE = re.compile(r"(20(?:1[2-9]|2[0-9]))")

# Geographic display crop for Europe (lon_min, lat_min, lon_max, lat_max).
# Includes Iceland and Cyprus; adjust later if the analytical study area differs.
EUROPE_BOUNDS_4326 = (-25.0, 34.0, 36.0, 72.5)


def extract_year(path: Path) -> int:
    """Extract the first plausible VIIRS year from a filename."""
    match = YEAR_RE.search(path.name)
    if not match:
        raise ValueError(f"Could not identify year from filename: {path.name}")
    return int(match.group(1))


def discover_rasters(data_dir: Path, pattern: str) -> dict[int, Path]:
    """Find one raster per year and return {year: path}."""
    found: dict[int, Path] = {}

    for path in sorted(data_dir.rglob(pattern)):
        try:
            year = extract_year(path)
        except ValueError:
            continue

        if year in found:
            raise RuntimeError(
                f"More than one file found for {year}:\n"
                f"  {found[year]}\n"
                f"  {path}\n"
                "Use a more specific --pattern."
            )
        found[year] = path

    if not found:
        raise FileNotFoundError(
            f"No rasters matching {pattern!r} found under {data_dir}"
        )

    return dict(sorted(found.items()))


def _clip_window_to_dataset(window: Window, width: int, height: int) -> Window:
    """Clip a raster window to the dataset extent."""
    full = Window(0, 0, width, height)
    clipped = window.intersection(full)
    if clipped.width <= 0 or clipped.height <= 0:
        raise ValueError("Requested map bounds do not overlap raster.")
    return clipped


def read_for_display(
    path: Path,
    bounds_4326: tuple[float, float, float, float],
    max_pixels_across: int = 2200,
) -> tuple[np.ndarray, Affine, rasterio.crs.CRS]:
    """
    Read only the European window and downsample for plotting.

    Downsampling uses average resampling because this is a radiance surface.
    The full-resolution raster remains untouched for later zonal statistics.
    """
    with rasterio.open(path) as src:
        if src.crs is None:
            raise ValueError(f"Raster has no CRS: {path}")

        src_bounds = transform_bounds(
            "EPSG:4326",
            src.crs,
            *bounds_4326,
            densify_pts=21,
        )

        window = from_bounds(*src_bounds, transform=src.transform)
        window = _clip_window_to_dataset(
            window.round_offsets().round_lengths(), src.width, src.height
        )

        scale = max(1.0, window.width / max_pixels_across)
        out_width = max(1, int(round(window.width / scale)))
        out_height = max(1, int(round(window.height / scale)))

        data = src.read(
            1,
            window=window,
            out_shape=(out_height, out_width),
            masked=True,
            resampling=Resampling.average,
        ).astype("float32")

        # Transform corresponding to the downsampled display array.
        native_transform = src.window_transform(window)
        display_transform = native_transform * Affine.scale(
            window.width / out_width,
            window.height / out_height,
        )

        arr = np.asarray(data.filled(np.nan), dtype="float32")

        # For visualisation only: negative annual radiances are not displayed.
        # Do not use this clipping in the later analytical extraction stage.
        arr[arr < 0] = 0.0

        return arr, display_transform, src.crs


def raster_extent(
    array: np.ndarray, transform: Affine
) -> tuple[float, float, float, float]:
    """Return matplotlib extent=(left, right, bottom, top)."""
    bottom, left = None, None
    west, south, east, north = array_bounds(
        array.shape[0], array.shape[1], transform
    )
    return west, east, south, north


def load_nuts(
    nuts_path: Path | None,
    target_crs,
) -> gpd.GeoDataFrame | None:
    """Load NUTS3 boundaries for an unobtrusive line overlay."""
    if nuts_path is None:
        return None

    nuts = gpd.read_file(nuts_path)
    if nuts.crs is None:
        raise ValueError(f"NUTS file has no CRS: {nuts_path}")

    if target_crs is not None and nuts.crs != target_crs:
        nuts = nuts.to_crs(target_crs)

    return nuts


def shared_radiance_norm(
    arrays: dict[int, np.ndarray],
    percentile: float = 99.7,
    linthresh: float = 0.2,
) -> SymLogNorm:
    """
    Build one common radiance norm for all years.

    The linear region close to zero keeps dark/unlit areas interpretable;
    above linthresh the scale behaves logarithmically.
    """
    yearly_highs = []

    for arr in arrays.values():
        vals = arr[np.isfinite(arr) & (arr > 0)]
        if vals.size:
            yearly_highs.append(float(np.percentile(vals, percentile)))

    if not yearly_highs:
        raise ValueError("No positive finite radiance values found.")

    vmax = max(yearly_highs)

    return SymLogNorm(
        linthresh=linthresh,
        linscale=0.8,
        vmin=0.0,
        vmax=vmax,
        base=10,
    )


def panel_layout(n_panels: int, ncols: int = 4) -> tuple[int, int]:
    ncols = min(ncols, n_panels)
    nrows = math.ceil(n_panels / ncols)
    return nrows, ncols


def style_map_axis(ax):
    """Minimal map styling: the raster, boundaries and year carry the figure."""
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_annual_small_multiples(
    arrays: dict[int, np.ndarray],
    transforms: dict[int, Affine],
    nuts: gpd.GeoDataFrame | None,
    output: Path,
    ncols: int = 4,
):
    """Absolute radiance for every year, using one common scale."""
    years = list(arrays)
    norm = shared_radiance_norm(arrays)
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#050505")

    nrows, ncols = panel_layout(len(years), ncols=ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.35 * ncols, 3.4 * nrows),
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes).ravel()

    image = None

    for ax, year in zip(axes, years):
        arr = arrays[year]
        extent = raster_extent(arr, transforms[year])

        image = ax.imshow(
            arr,
            extent=extent,
            origin="upper",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
            rasterized=True,
        )

        if nuts is not None:
            nuts.boundary.plot(
                ax=ax,
                color="white",
                linewidth=0.18,
                alpha=0.42,
                zorder=2,
            )

        ax.set_title(str(year), fontsize=12, weight="semibold", pad=5)
        style_map_axis(ax)

    for ax in axes[len(years):]:
        ax.remove()

    cbar = fig.colorbar(
        image,
        ax=axes[: len(years)].tolist(),
        orientation="horizontal",
        shrink=0.72,
        pad=0.025,
        aspect=45,
    )
    cbar.set_label(
        r"VIIRS annual radiance (nW cm$^{-2}$ sr$^{-1}$)",
        fontsize=10,
    )

    fig.suptitle(
        "Annual VIIRS night-time radiance",
        fontsize=16,
        weight="semibold",
    )

    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_change_from_baseline(
    arrays: dict[int, np.ndarray],
    transforms: dict[int, Affine],
    nuts: gpd.GeoDataFrame | None,
    output: Path,
    ncols: int = 4,
    percentile: float = 99.0,
):
    """
    Plot radiance change relative to the first year.

    This is often much more informative than inspecting annual maps alone.
    """
    years = list(arrays)
    baseline_year = years[0]
    baseline = arrays[baseline_year]

    changes: dict[int, np.ndarray] = {}

    for year in years[1:]:
        if arrays[year].shape != baseline.shape:
            raise ValueError(
                "Display grids differ between years. "
                "For these source files, add a reprojection/alignment stage."
            )
        changes[year] = arrays[year] - baseline

    if not changes:
        return

    # Robust shared symmetric scale.
    highs = []
    for arr in changes.values():
        vals = np.abs(arr[np.isfinite(arr)])
        if vals.size:
            highs.append(float(np.percentile(vals, percentile)))

    limit = max(highs)
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)

    years_to_plot = list(changes)
    nrows, ncols = panel_layout(len(years_to_plot), ncols=ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.35 * ncols, 3.4 * nrows),
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes).ravel()

    image = None

    for ax, year in zip(axes, years_to_plot):
        arr = changes[year]
        extent = raster_extent(arr, transforms[year])

        image = ax.imshow(
            arr,
            extent=extent,
            origin="upper",
            cmap="coolwarm",
            norm=norm,
            interpolation="nearest",
            rasterized=True,
        )

        if nuts is not None:
            nuts.boundary.plot(
                ax=ax,
                color="black",
                linewidth=0.16,
                alpha=0.35,
                zorder=2,
            )

        ax.set_title(
            f"{year} − {baseline_year}",
            fontsize=12,
            weight="semibold",
            pad=5,
        )
        style_map_axis(ax)

    for ax in axes[len(years_to_plot):]:
        ax.remove()

    cbar = fig.colorbar(
        image,
        ax=axes[: len(years_to_plot)].tolist(),
        orientation="horizontal",
        shrink=0.72,
        pad=0.025,
        aspect=45,
    )
    cbar.set_label(
        rf"Radiance change from {baseline_year} "
        r"(nW cm$^{-2}$ sr$^{-1}$)",
        fontsize=10,
    )

    fig.suptitle(
        f"Change in annual VIIRS radiance relative to {baseline_year}",
        fontsize=16,
        weight="semibold",
    )

    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Directory containing annual VIIRS GeoTIFFs.",
    )
    parser.add_argument(
        "--pattern",
        default="*.tif",
        help='Recursive glob, e.g. "*average-masked*.tif".',
    )
    parser.add_argument(
        "--nuts",
        type=Path,
        default=None,
        help="Optional NUTS3 polygon file.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("figures"),
    )
    parser.add_argument(
        "--max-pixels-across",
        type=int,
        default=2200,
        help="Maximum width of each display raster.",
    )
    args = parser.parse_args()

    rasters = discover_rasters(args.data_dir, args.pattern)
    print(f"Found {len(rasters)} annual rasters:")
    for year, path in rasters.items():
        print(f"  {year}: {path}")

    arrays: dict[int, np.ndarray] = {}
    transforms: dict[int, Affine] = {}
    raster_crs = None

    for year, path in rasters.items():
        arr, transform, crs = read_for_display(
            path,
            EUROPE_BOUNDS_4326,
            max_pixels_across=args.max_pixels_across,
        )
        arrays[year] = arr
        transforms[year] = transform

        if raster_crs is None:
            raster_crs = crs
        elif crs != raster_crs:
            raise ValueError("Annual rasters do not share a common CRS.")

    nuts = load_nuts(args.nuts, raster_crs)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    absolute_path = args.out_dir / "viirs_annual_small_multiples.pdf"
    change_path = args.out_dir / "viirs_change_from_baseline.pdf"

    plot_annual_small_multiples(
        arrays,
        transforms,
        nuts,
        absolute_path,
    )
    plot_change_from_baseline(
        arrays,
        transforms,
        nuts,
        change_path,
    )

    print(f"Wrote: {absolute_path}")
    print(f"Wrote: {change_path}")


if __name__ == "__main__":
    main()
