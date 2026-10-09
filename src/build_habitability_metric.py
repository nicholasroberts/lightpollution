#!/usr/bin/env python3
"""Build a static terrain-constrained NUTS3 habitability fraction H.

H is deliberately independent of current settlement, roads, GDP, VIIRS, and
future projections.  It estimates the proportion of a NUTS3 polygon that is
physically available for human habitation under configurable terrain limits.

For each ~90 m Copernicus DEM cell, three binary terrain-suitability scenarios
are evaluated from elevation, slope, and 3x3 local elevation SD.  Permanent
water and permanent snow/ice are then removed using Copernicus 100 m fractional
cover layers.

Primary metric:
    H_core = effective habitable area / sampled NUTS3 area

Sensitivity variants H_strict and H_permissive are produced so that later
historical validation can test whether conclusions depend on terrain
thresholds.  The metric is aggregated independently onto every configured NUTS
release (2016, 2021, 2024), allowing historical models to use the correct
geography while future work can use NUTS 2024.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.warp import reproject
from scipy.ndimage import uniform_filter
from shapely.geometry import box

from project_config import configured_path, load_config


EARTH_RADIUS_M = 6_371_008.8


def nuts_path(config: dict, release: int) -> Path:
    name = config["nuts"]["filename_template"].format(release=release)
    return configured_path(config, "nuts_raw") / name


def load_nuts_releases(config: dict):
    releases = {}
    for release in config["nuts"]["releases"]:
        release = int(release)
        path = nuts_path(config, release)
        if not path.exists():
            raise FileNotFoundError(
                f"Missing NUTS {release} boundaries: {path}. "
                "Run: python src/run_pipeline.py nuts"
            )
        g = gpd.read_file(path)
        if g.crs is None:
            raise ValueError(f"NUTS file has no CRS: {path}")
        g = g.to_crs("EPSG:4326")
        g = g.reset_index(drop=True)
        g["_row"] = np.arange(len(g), dtype=int)
        releases[release] = g
    return releases


def local_elevation_sd(elevation: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """3x3 local elevation standard deviation with missing-value weighting."""
    x = np.where(valid, elevation, 0.0).astype(np.float64)
    m = valid.astype(np.float64)

    # uniform_filter returns a local mean, so multiply by nine to obtain sums.
    count = uniform_filter(m, size=3, mode="nearest") * 9.0
    sx = uniform_filter(x, size=3, mode="nearest") * 9.0
    sx2 = uniform_filter(x * x, size=3, mode="nearest") * 9.0

    mean = np.divide(sx, count, out=np.zeros_like(sx), where=count > 0)
    mean2 = np.divide(sx2, count, out=np.zeros_like(sx2), where=count > 0)
    var = np.maximum(0.0, mean2 - mean * mean)
    sd = np.sqrt(var)
    sd[count < 5] = np.nan
    return sd.astype(np.float32)


def slope_degrees(elevation: np.ndarray, valid: np.ndarray, transform) -> np.ndarray:
    """Approximate physical slope for a geographic Copernicus DEM tile."""
    z = elevation.astype(np.float64).copy()

    # Copernicus DEM coastal cells are generally populated; this fallback only
    # prevents rare nodata cells from contaminating neighbouring gradients.
    if not np.all(valid):
        fallback = float(np.nanmedian(z[valid])) if np.any(valid) else 0.0
        z[~valid] = fallback

    h, _ = z.shape
    lat_centres = (
        transform.f
        + (np.arange(h, dtype=float) + 0.5) * transform.e
    )
    lat_rad = np.deg2rad(lat_centres)

    # Meridional and zonal metres per degree on WGS84, sufficiently accurate
    # for slope at the scale of one 1-degree DEM tile.
    metres_per_deg_lat = (
        111132.92
        - 559.82 * np.cos(2.0 * lat_rad)
        + 1.175 * np.cos(4.0 * lat_rad)
    )
    metres_per_deg_lon = (
        111412.84 * np.cos(lat_rad)
        - 93.5 * np.cos(3.0 * lat_rad)
    )

    dy = np.abs(transform.e) * metres_per_deg_lat
    dx = np.abs(transform.a) * metres_per_deg_lon

    dz_drow = np.gradient(z, axis=0)
    dz_dcol = np.gradient(z, axis=1)

    dz_dy = dz_drow / dy[:, None]
    dz_dx = dz_dcol / dx[:, None]

    slope = np.degrees(np.arctan(np.sqrt(dz_dx * dz_dx + dz_dy * dz_dy)))
    slope[~valid] = np.nan
    return slope.astype(np.float32)


def pixel_area_km2(transform, height: int, width: int) -> np.ndarray:
    """Spherical physical area for every raster row, broadcast to 2-D."""
    lat_edges = transform.f + np.arange(height + 1, dtype=float) * transform.e
    lat1 = np.deg2rad(lat_edges[:-1])
    lat2 = np.deg2rad(lat_edges[1:])
    dlon = np.deg2rad(abs(transform.a))
    row_area_m2 = (
        EARTH_RADIUS_M**2
        * dlon
        * np.abs(np.sin(lat1) - np.sin(lat2))
    )
    return np.broadcast_to((row_area_m2 / 1e6)[:, None], (height, width))


def resample_fraction(src, dst_shape, dst_transform, dst_crs):
    out = np.full(dst_shape, np.nan, dtype=np.float32)
    reproject(
        source=rasterio.band(src, 1),
        destination=out,
        src_transform=src.transform,
        src_crs=src.crs,
        src_nodata=src.nodata,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        dst_nodata=np.nan,
        resampling=Resampling.average,
    )
    return out


def physical_land_fraction(water: np.ndarray, snow: np.ndarray) -> np.ndarray:
    """Fraction of each DEM cell not permanent water or permanent snow/ice."""
    invalid = (
        ~np.isfinite(water)
        | ~np.isfinite(snow)
        | (water < 0)
        | (snow < 0)
        | (water > 100)
        | (snow > 100)
    )
    land = 1.0 - np.clip(water + snow, 0.0, 100.0) / 100.0
    land[invalid] = 0.0
    return land.astype(np.float32)


def make_accumulator(nuts: gpd.GeoDataFrame):
    n = len(nuts)
    fields = [
        "sampled_area_km2",
        "valid_dem_area_km2",
        "physical_land_area_km2",
        "habitable_area_strict_km2",
        "habitable_area_core_km2",
        "habitable_area_permissive_km2",
        "elevation_area_sum",
        "slope_area_sum",
        "local_sd_area_sum",
    ]
    return {name: np.zeros(n, dtype=np.float64) for name in fields}


def terrain_mask(elev, slope, relief, valid, thresholds):
    return (
        valid
        & np.isfinite(elev)
        & np.isfinite(slope)
        & np.isfinite(relief)
        & (elev <= float(thresholds["max_elevation_m"]))
        & (slope <= float(thresholds["max_slope_deg"]))
        & (relief <= float(thresholds["max_local_elevation_sd_m"]))
    )


def add_bincount(acc, indices, labels, weights_by_metric):
    label_flat = labels.ravel()
    max_label = len(indices)
    for field, weights in weights_by_metric.items():
        sums = np.bincount(
            label_flat,
            weights=weights.ravel(),
            minlength=max_label + 1,
        )
        if len(sums) < max_label + 1:
            sums = np.pad(sums, (0, max_label + 1 - len(sums)))
        acc[field][indices] += sums[1 : max_label + 1]


def process_tile(
    tile_path: Path,
    water_src,
    snow_src,
    releases,
    accumulators,
    cfg,
):
    with rasterio.open(tile_path) as src:
        if src.crs is None:
            raise ValueError(f"DEM tile has no CRS: {tile_path}")
        if src.crs.to_string() != "EPSG:4326":
            raise ValueError(f"Unexpected DEM CRS {src.crs} in {tile_path}")

        raw = src.read(1, masked=True)
        elev = raw.filled(np.nan).astype(np.float32)
        valid = np.isfinite(elev)
        if src.nodata is not None:
            valid &= elev != src.nodata

        slope = slope_degrees(elev, valid, src.transform)
        relief = local_elevation_sd(elev, valid)
        water = resample_fraction(
            water_src, elev.shape, src.transform, src.crs
        )
        snow = resample_fraction(
            snow_src, elev.shape, src.transform, src.crs
        )
        land_fraction = physical_land_fraction(water, snow)

        masks = {}
        for name in ["strict", "core", "permissive"]:
            masks[name] = terrain_mask(
                elev,
                slope,
                relief,
                valid,
                cfg["scenarios"][name],
            ).astype(np.float32) * land_fraction

        area = pixel_area_km2(src.transform, src.height, src.width)
        tile_geom = box(*src.bounds)

        for release, nuts in releases.items():
            idx = nuts.sindex.query(tile_geom, predicate="intersects")
            if len(idx) == 0:
                continue
            idx = np.asarray(idx, dtype=int)
            local = nuts.iloc[idx]

            shapes = [
                (geom, label)
                for label, geom in enumerate(local.geometry, start=1)
                if geom is not None and not geom.is_empty
            ]
            if not shapes:
                continue

            labels = rasterize(
                shapes,
                out_shape=elev.shape,
                transform=src.transform,
                fill=0,
                dtype="int32",
                all_touched=False,
            )

            valid_f = valid.astype(np.float32)
            physical_land = valid_f * land_fraction
            valid_area = area * valid_f

            elev_safe = np.where(valid, elev, 0.0)
            slope_safe = np.where(np.isfinite(slope), slope, 0.0)
            relief_safe = np.where(np.isfinite(relief), relief, 0.0)

            add_bincount(
                accumulators[release],
                idx,
                labels,
                {
                    "sampled_area_km2": area,
                    "valid_dem_area_km2": valid_area,
                    "physical_land_area_km2": area * physical_land,
                    "habitable_area_strict_km2": area * masks["strict"],
                    "habitable_area_core_km2": area * masks["core"],
                    "habitable_area_permissive_km2": area * masks["permissive"],
                    "elevation_area_sum": valid_area * elev_safe,
                    "slope_area_sum": valid_area * slope_safe,
                    "local_sd_area_sum": valid_area * relief_safe,
                },
            )


def finalise_release(nuts, acc, release, config):
    cols = [
        c
        for c in [
            "NUTS_ID",
            "CNTR_CODE",
            "NUTS_NAME",
            "NAME_LATN",
        ]
        if c in nuts.columns
    ]
    out = nuts[cols].copy()
    out.insert(0, "nuts_release", release)

    for key in [
        "sampled_area_km2",
        "valid_dem_area_km2",
        "physical_land_area_km2",
        "habitable_area_strict_km2",
        "habitable_area_core_km2",
        "habitable_area_permissive_km2",
    ]:
        out[key] = acc[key]

    denom = out["sampled_area_km2"].replace(0.0, np.nan)
    for name in ["strict", "core", "permissive"]:
        out[f"habitable_fraction_{name}"] = (
            out[f"habitable_area_{name}_km2"] / denom
        ).clip(0.0, 1.0)

    out["physical_land_fraction"] = (
        out["physical_land_area_km2"] / denom
    ).clip(0.0, 1.0)

    valid_denom = out["valid_dem_area_km2"].replace(0.0, np.nan)
    out["mean_elevation_m"] = acc["elevation_area_sum"] / valid_denom
    out["mean_slope_deg"] = acc["slope_area_sum"] / valid_denom
    out["mean_local_elevation_sd_m"] = acc["local_sd_area_sum"] / valid_denom

    # Independent equal-area polygon area is a useful coverage diagnostic.
    area_gdf = nuts[["geometry"]].to_crs("EPSG:3035")
    out["polygon_area_epsg3035_km2"] = area_gdf.area.to_numpy() / 1e6
    out["rasterised_area_coverage_ratio"] = (
        out["sampled_area_km2"]
        / out["polygon_area_epsg3035_km2"].replace(0.0, np.nan)
    )

    out["terrain_compression_factor_core"] = (
        1.0 / out["habitable_fraction_core"].replace(0.0, np.nan)
    )

    # Join 2024 population density only where geography matches NUTS 2024.
    if release == 2024:
        processed = configured_path(config, "processed")
        merged = pd.read_parquet(
            processed / config["outputs"]["merged"],
            columns=["year", "NUTS_ID", "population_density"],
        )
        pop = (
            merged[merged["year"] == 2024][["NUTS_ID", "population_density"]]
            .drop_duplicates("NUTS_ID")
        )
        out = out.merge(pop, on="NUTS_ID", how="left")
        out = out.rename(
            columns={"population_density": "population_density_2024"}
        )
        min_h = float(
            config["habitability"]["min_habitable_fraction_for_adjusted_density"]
        )
        h = out["habitable_fraction_core"]
        out["terrain_adjusted_population_density_2024"] = np.where(
            h >= min_h,
            out["population_density_2024"] / h,
            np.nan,
        )

    return out


def plot_map(nuts, table, value, title, label, output, cfg):
    merged = nuts.merge(
        table[["NUTS_ID", value]],
        on="NUTS_ID",
        how="left",
    )
    west, south, east, north = cfg["map_extent_4326"]
    extent = box(west, south, east, north)
    merged = merged[merged.geometry.intersects(extent)].copy()
    merged = merged.to_crs("EPSG:3035")

    fig, ax = plt.subplots(figsize=(10.5, 9.0))
    merged.plot(
        column=value,
        ax=ax,
        cmap="viridis",
        linewidth=0.10,
        edgecolor="white",
        missing_kwds={"color": "lightgrey"},
        vmin=0.0 if "fraction" in value else None,
        vmax=1.0 if "fraction" in value else None,
    )
    ax.set_axis_off()
    ax.set_title(title)

    norm = Normalize(
        vmin=0.0 if "fraction" in value else np.nanmin(merged[value]),
        vmax=1.0 if "fraction" in value else np.nanmax(merged[value]),
    )
    sm = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.01)
    cb.set_label(label)
    fig.tight_layout()
    fig.savefig(output, dpi=int(cfg["map_dpi"]), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["habitability"]
    raw_root = configured_path(config, "habitability_raw")
    manifest_path = raw_root / "copernicus_dem_glo90_manifest.csv"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Missing DEM manifest: {manifest_path}. Run: "
            "python src/run_pipeline.py habitability-fetch"
        )

    manifest = pd.read_csv(manifest_path)
    releases = load_nuts_releases(config)
    accumulators = {
        release: make_accumulator(nuts)
        for release, nuts in releases.items()
    }

    lc_dir = raw_root / "copernicus_landcover_2015"
    water_path = lc_dir / cfg["landcover"]["permanent_water"]["filename"]
    snow_path = lc_dir / cfg["landcover"]["snow_ice"]["filename"]
    if not water_path.exists() or not snow_path.exists():
        raise FileNotFoundError(
            "Permanent-water or snow/ice layer missing. Run: "
            "python src/run_pipeline.py habitability-fetch"
        )

    with rasterio.open(water_path) as water_src, rasterio.open(snow_path) as snow_src:
        for i, row in manifest.iterrows():
            tile_path = raw_root / row["relative_path"]
            if not tile_path.exists():
                raise FileNotFoundError(f"Missing DEM tile: {tile_path}")
            print(
                f"[{i+1:,}/{len(manifest):,}] Processing {tile_path.name}",
                flush=True,
            )
            process_tile(
                tile_path,
                water_src,
                snow_src,
                releases,
                accumulators,
                cfg,
            )

    frames = []
    for release, nuts in releases.items():
        frames.append(
            finalise_release(
                nuts,
                accumulators[release],
                release,
                config,
            )
        )
    combined = pd.concat(frames, ignore_index=True)

    processed = configured_path(config, "processed")
    processed.mkdir(parents=True, exist_ok=True)
    parquet_path = processed / config["outputs"]["nuts3_habitability"]
    csv_path = parquet_path.with_suffix(".csv")
    combined.to_parquet(parquet_path, index=False)
    combined.to_csv(csv_path, index=False)

    figure_dir = configured_path(config, "figures") / "habitability"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_2024 = combined[combined["nuts_release"] == 2024].copy()
    nuts_2024 = releases[2024]

    h_map = figure_dir / "nuts3_habitable_fraction_core_2024.pdf"
    plot_map(
        nuts_2024,
        table_2024,
        "habitable_fraction_core",
        "Terrain-constrained habitable fraction (H), NUTS3",
        "Habitable fraction H",
        h_map,
        cfg,
    )

    compression_map = figure_dir / "nuts3_terrain_compression_factor_2024.pdf"
    # Winsorise only for display; the underlying table remains unmodified.
    display = table_2024.copy()
    finite = display["terrain_compression_factor_core"].replace(
        [np.inf, -np.inf], np.nan
    )
    cap = float(finite.quantile(0.98))
    display["terrain_compression_factor_core_display"] = finite.clip(upper=cap)
    plot_map(
        nuts_2024,
        display.rename(
            columns={
                "terrain_compression_factor_core_display":
                "terrain_compression_factor_core_map"
            }
        ),
        "terrain_compression_factor_core_map",
        "Terrain compression of nominal population density, NUTS3",
        f"1 / H (display capped at 98th percentile = {cap:.2f})",
        compression_map,
        cfg,
    )

    summary_path = processed / "nuts3_habitability_summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("STATIC TERRAIN-CONSTRAINED NUTS3 HABITABILITY METRIC H\n")
        handle.write("=" * 72 + "\n\n")
        handle.write(
            "H is independent of population, current settlement, roads, GDP, "
            "VIIRS radiance, and future projections. It uses only DEM-derived "
            "terrain plus permanent water and permanent snow/ice masks.\n\n"
        )
        handle.write("SCENARIOS\n")
        for name in ["strict", "core", "permissive"]:
            s = cfg["scenarios"][name]
            handle.write(
                f"{name}: elevation <= {s['max_elevation_m']} m; "
                f"slope <= {s['max_slope_deg']} deg; "
                f"3x3 elevation SD <= {s['max_local_elevation_sd_m']} m\n"
            )
        handle.write("\n")
        handle.write(
            "Core definition: H = effective physically habitable raster area "
            "/ total rasterised NUTS3 area. Permanent-water and snow/ice cover "
            "fractions reduce the eligible area continuously.\n\n"
        )

        for release in sorted(combined["nuts_release"].unique()):
            g = combined[combined["nuts_release"] == release]
            handle.write(f"NUTS {release}: {len(g)} regions\n")
            handle.write(
                "  H_core median="
                f"{g['habitable_fraction_core'].median():.4f}; "
                "10th="
                f"{g['habitable_fraction_core'].quantile(0.10):.4f}; "
                "90th="
                f"{g['habitable_fraction_core'].quantile(0.90):.4f}\n"
            )

        handle.write("\nLOWEST H_core, NUTS 2024\n")
        cols = [
            c for c in ["NUTS_ID", "NUTS_NAME", "CNTR_CODE",
                        "habitable_fraction_core",
                        "mean_elevation_m", "mean_slope_deg",
                        "mean_local_elevation_sd_m"]
            if c in table_2024.columns
        ]
        handle.write(
            table_2024.nsmallest(20, "habitable_fraction_core")[cols]
            .to_string(index=False, float_format=lambda x: f"{x:.4f}")
        )
        handle.write("\n\nHIGHEST H_core, NUTS 2024\n")
        handle.write(
            table_2024.nlargest(20, "habitable_fraction_core")[cols]
            .to_string(index=False, float_format=lambda x: f"{x:.4f}")
        )
        handle.write("\n\n")
        handle.write(
            "IMPORTANT: the strict/core/permissive thresholds are explicit "
            "model assumptions, not claims of a universal physiological limit. "
            "Their value is testable: historical VIIRS change can determine "
            "whether H improves prediction and whether results are robust to "
            "the threshold scenario.\n"
        )

    print()
    print(f"Wrote: {parquet_path}")
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {h_map}")
    print(f"Wrote: {compression_map}")
    print()
    print("NUTS 2024 H_core summary:")
    print(table_2024["habitable_fraction_core"].describe().to_string())


if __name__ == "__main__":
    main()
