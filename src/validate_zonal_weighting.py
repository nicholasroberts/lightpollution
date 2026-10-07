#!/usr/bin/env python3
"""Compare candidate NUTS3 zonal-mean methods on a representative sample.

This is a validation stage, not the production extractor.

Methods compared
----------------
1. pixel_center_mean
   Mean of valid raster cells whose centres fall inside the polygon
   (rasterio mask with all_touched=False).

2. fractional_overlap_mean
   exactextract mean weighted by the exact fraction of each raster cell
   covered by the NUTS3 polygon.

3. area_weighted_mean
   exactextract mean where coverage is expressed as spherical surface area
   (coverage_weight=area_spherical_m2). This accounts for both exact polygon
   overlap and the declining physical area of longitude/latitude cells toward
   higher latitudes.

The script selects a deterministic sample of European NUTS3 regions spanning
latitude, area, and geometric compactness, then writes a CSV and diagnostic
figure comparing the methods.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from rasterio.mask import mask
from shapely.geometry import mapping

from project_config import configured_path, load_config


YEAR_RE = re.compile(r"(20\d{2})")


def raster_year(path: Path) -> int | None:
    match = YEAR_RE.search(path.name)
    return int(match.group(1)) if match else None


def find_raster(config: dict, year: int) -> Path:
    root = configured_path(config, "viirs_raw")
    pattern = config["viirs"]["raster_pattern"]

    matches = [
        path
        for path in sorted(root.rglob(pattern))
        if raster_year(path) == year
    ]

    if not matches:
        raise FileNotFoundError(
            f"No VIIRS raster for {year} matching {pattern!r} under {root}"
        )
    if len(matches) > 1:
        joined = "\n  ".join(str(path) for path in matches)
        raise RuntimeError(
            f"More than one VIIRS raster found for {year}:\n  {joined}"
        )

    return matches[0]


def load_nuts(config: dict) -> gpd.GeoDataFrame:
    path = (
        configured_path(config, "nuts_raw")
        / config["nuts"]["expected_filename"]
    )
    if not path.exists():
        raise FileNotFoundError(f"NUTS3 file not found: {path}")

    nuts = gpd.read_file(path)
    expected_crs = config["nuts"]["crs"]

    if nuts.crs is None:
        raise ValueError(f"NUTS3 file has no CRS: {path}")
    if nuts.crs.to_string() != expected_crs:
        nuts = nuts.to_crs(expected_crs)

    id_field = config["nuts"]["id_field"]
    if id_field not in nuts.columns:
        raise ValueError(f"Missing NUTS identifier field {id_field!r}")

    return nuts


def select_representative_regions(
    nuts: gpd.GeoDataFrame,
    config: dict,
    sample_size: int,
) -> gpd.GeoDataFrame:
    """Select regions spanning latitude, area, and geometric compactness."""
    west, south, east, north = config["viirs"]["display"]["europe_bounds_4326"]

    work = nuts.copy()
    representative = work.geometry.representative_point()
    work["sample_lon"] = representative.x
    work["sample_lat"] = representative.y

    work = work[
        work["sample_lon"].between(west, east)
        & work["sample_lat"].between(south, north)
        & ~work.geometry.is_empty
        & work.geometry.notna()
    ].copy()

    if work.empty:
        raise ValueError("No NUTS3 regions fall within configured Europe bounds.")

    equal_area = work.to_crs("EPSG:3035")
    work["area_km2"] = equal_area.geometry.area.to_numpy() / 1_000_000.0
    perimeter = equal_area.geometry.length.to_numpy()
    area_m2 = equal_area.geometry.area.to_numpy()
    work["compactness"] = np.where(
        perimeter > 0,
        4.0 * math.pi * area_m2 / np.square(perimeter),
        np.nan,
    )

    selections: list[tuple[int, str]] = []
    seen: set[str] = set()
    id_field = config["nuts"]["id_field"]

    def add(index: int, reason: str) -> None:
        nuts_id = str(work.loc[index, id_field])
        if nuts_id not in seen:
            selections.append((index, reason))
            seen.add(nuts_id)

    # Five latitude positions ensure the test spans southern to northern Europe.
    for quantile in (0.05, 0.25, 0.50, 0.75, 0.95):
        target = float(work["sample_lat"].quantile(quantile))
        index = (work["sample_lat"] - target).abs().idxmin()
        add(index, f"latitude_q{int(quantile * 100):02d}")

    # Add geometric/size extremes to stress boundary treatment.
    add(work["area_km2"].idxmin(), "smallest_area")
    add(work["area_km2"].idxmax(), "largest_area")
    add(work["compactness"].idxmin(), "least_compact")
    add(work["compactness"].idxmax(), "most_compact")

    # If duplicates reduced the sample, fill deterministically by latitude.
    if len(selections) < sample_size:
        ordered = work.sort_values(
            ["sample_lat", id_field],
            kind="stable",
        )
        positions = np.linspace(
            0,
            len(ordered) - 1,
            num=min(sample_size * 2, len(ordered)),
            dtype=int,
        )
        for pos in positions:
            add(ordered.index[pos], "latitude_fill")
            if len(selections) >= sample_size:
                break

    selections = selections[:sample_size]

    selected = work.loc[[idx for idx, _ in selections]].copy()
    selected["selection_reason"] = [reason for _, reason in selections]

    columns = [
        id_field,
        "CNTR_CODE",
        "NUTS_NAME",
        "NAME_LATN",
        "sample_lon",
        "sample_lat",
        "area_km2",
        "compactness",
        "selection_reason",
        "geometry",
    ]
    columns = [column for column in columns if column in selected.columns]

    return selected[columns].reset_index(drop=True)


def pixel_center_mean(
    raster_path: Path,
    geometry,
    geometry_crs,
) -> tuple[float, int]:
    """Mean of valid cells whose centres are inside the polygon."""
    with rasterio.open(raster_path) as src:
        geom = geometry
        if geometry_crs != src.crs:
            geom = gpd.GeoSeries(
                [geometry],
                crs=geometry_crs,
            ).to_crs(src.crs).iloc[0]

        clipped, _ = mask(
            src,
            [mapping(geom)],
            crop=True,
            filled=False,
            all_touched=False,
            indexes=1,
        )

    values = np.ma.asarray(clipped).compressed()
    values = values[np.isfinite(values)]

    if values.size == 0:
        return float("nan"), 0

    return float(values.mean()), int(values.size)


def exact_means(
    raster_path: Path,
    selected: gpd.GeoDataFrame,
    id_field: str,
) -> pd.DataFrame:
    """Calculate exact-overlap and physical-area-weighted means."""
    keep = selected[[id_field, "geometry"]].copy()

    result = exact_extract(
        str(raster_path),
        keep,
        [
            "fractional_overlap_mean=mean",
            "area_weighted_mean=mean(coverage_weight=area_spherical_m2)",
        ],
        include_cols=[id_field],
        output="pandas",
        progress=True,
    )

    return result


def safe_percent_difference(new: pd.Series, reference: pd.Series) -> pd.Series:
    denominator = reference.replace(0, np.nan)
    return 100.0 * (new - reference) / denominator


def make_figure(results: pd.DataFrame, output: Path, id_field: str) -> None:
    labels = (
        results[id_field].astype(str)
        + "  "
        + results["NUTS_NAME"].fillna(results.get("NAME_LATN", "")).astype(str)
    )
    y = np.arange(len(results))

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13, max(5.5, 0.55 * len(results) + 2)),
        constrained_layout=True,
        sharey=True,
    )

    axes[0].barh(y, results["fractional_vs_center_pct"])
    axes[0].axvline(0, color="black", linewidth=0.8)
    axes[0].set_title("Exact overlap vs pixel-centre mean")
    axes[0].set_xlabel("Difference (%)")
    axes[0].set_yticks(y, labels=labels)
    axes[0].invert_yaxis()

    axes[1].barh(y, results["area_vs_fractional_pct"])
    axes[1].axvline(0, color="black", linewidth=0.8)
    axes[1].set_title("Physical-area weighting vs exact overlap")
    axes[1].set_xlabel("Difference (%)")

    fig.suptitle(
        "Sensitivity of NUTS3 mean VIIRS radiance to zonal weighting method",
        fontsize=14,
        weight="semibold",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument(
        "--year",
        type=int,
        default=None,
        help="Override validation year from config.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=None,
        help="Override number of representative NUTS3 regions.",
    )
    parser.add_argument(
        "--nuts-ids",
        nargs="+",
        default=None,
        help="Optional explicit NUTS_ID values instead of automatic selection.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    validation = config["analysis"]["area_weighting"]["validation"]

    year = args.year or int(validation["year"])
    sample_size = args.sample_size or int(validation["sample_size"])
    id_field = config["nuts"]["id_field"]

    raster_path = find_raster(config, year)
    nuts = load_nuts(config)

    if args.nuts_ids:
        selected = nuts[nuts[id_field].isin(args.nuts_ids)].copy()
        missing = sorted(set(args.nuts_ids) - set(selected[id_field]))
        if missing:
            raise ValueError(f"Unknown NUTS_ID values: {missing}")
        selected["selection_reason"] = "explicit"
        representative = selected.geometry.representative_point()
        selected["sample_lon"] = representative.x
        selected["sample_lat"] = representative.y
        equal_area = selected.to_crs("EPSG:3035")
        selected["area_km2"] = equal_area.geometry.area.to_numpy() / 1_000_000.0
        perimeter = equal_area.geometry.length.to_numpy()
        area_m2 = equal_area.geometry.area.to_numpy()
        selected["compactness"] = np.where(
            perimeter > 0,
            4.0 * math.pi * area_m2 / np.square(perimeter),
            np.nan,
        )
    else:
        selected = select_representative_regions(
            nuts,
            config,
            sample_size=sample_size,
        )

    print(f"Raster: {raster_path}")
    print(f"Validation year: {year}")
    print(f"Regions selected: {len(selected)}")
    print()
    print(
        selected[
            [
                column
                for column in [
                    id_field,
                    "CNTR_CODE",
                    "NUTS_NAME",
                    "sample_lat",
                    "area_km2",
                    "selection_reason",
                ]
                if column in selected.columns
            ]
        ].to_string(index=False)
    )
    print()

    exact = exact_means(raster_path, selected, id_field)

    center_rows = []
    for row in selected.itertuples(index=False):
        geometry = row.geometry
        mean_value, pixel_count = pixel_center_mean(
            raster_path,
            geometry,
            selected.crs,
        )
        center_rows.append(
            {
                id_field: getattr(row, id_field),
                "pixel_center_mean": mean_value,
                "center_pixel_count": pixel_count,
            }
        )

    center = pd.DataFrame(center_rows)

    metadata_columns = [
        column
        for column in [
            id_field,
            "CNTR_CODE",
            "NUTS_NAME",
            "NAME_LATN",
            "sample_lon",
            "sample_lat",
            "area_km2",
            "compactness",
            "selection_reason",
        ]
        if column in selected.columns
    ]

    results = (
        selected[metadata_columns]
        .merge(center, on=id_field, how="left", validate="one_to_one")
        .merge(exact, on=id_field, how="left", validate="one_to_one")
    )

    results["fractional_vs_center_abs"] = (
        results["fractional_overlap_mean"] - results["pixel_center_mean"]
    )
    results["fractional_vs_center_pct"] = safe_percent_difference(
        results["fractional_overlap_mean"],
        results["pixel_center_mean"],
    )

    results["area_vs_fractional_abs"] = (
        results["area_weighted_mean"] - results["fractional_overlap_mean"]
    )
    results["area_vs_fractional_pct"] = safe_percent_difference(
        results["area_weighted_mean"],
        results["fractional_overlap_mean"],
    )

    results["year"] = year
    results["viirs_file"] = raster_path.name

    processed = configured_path(config, "processed")
    figures = configured_path(config, "figures")
    processed.mkdir(parents=True, exist_ok=True)

    csv_name = validation["output_csv"].format(year=year)
    figure_name = validation["output_figure"].format(year=year)

    csv_path = processed / csv_name
    figure_path = figures / figure_name

    results.to_csv(csv_path, index=False)
    make_figure(results, figure_path, id_field)

    display_columns = [
        id_field,
        "NUTS_NAME",
        "sample_lat",
        "pixel_center_mean",
        "fractional_overlap_mean",
        "area_weighted_mean",
        "fractional_vs_center_pct",
        "area_vs_fractional_pct",
    ]
    display_columns = [
        column for column in display_columns if column in results.columns
    ]

    print("Method comparison:")
    with pd.option_context(
        "display.max_columns",
        None,
        "display.width",
        180,
        "display.float_format",
        "{:.6g}".format,
    ):
        print(results[display_columns].to_string(index=False))

    print()
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {figure_path}")


if __name__ == "__main__":
    main()
