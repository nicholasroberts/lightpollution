#!/usr/bin/env python3
"""Extract annual VIIRS radiance statistics for fixed NUTS3 regions.

Production method
-----------------
For each NUTS3 polygon, radiance summaries are calculated with exactextract.
Coverage is expressed as spherical surface area rather than simple cell
fraction:

    coverage_weight=area_spherical_m2

This means each raster contribution is weighted by:
  1. the exact fraction of the VIIRS cell covered by the NUTS3 polygon, and
  2. the physical area of that longitude/latitude cell.

The method was selected after validation against pixel-centre and uncorrected
fractional-overlap approaches on representative 2024 NUTS3 regions.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import geopandas as gpd
import pandas as pd
import rasterio
from exactextract import exact_extract

from project_config import configured_path, load_config


YEAR_RE = re.compile(r"(20\d{2})")


def raster_year(path: Path) -> int | None:
    match = YEAR_RE.search(path.name)
    return int(match.group(1)) if match else None


def discover_rasters(config: dict, requested_years: list[int] | None) -> dict[int, Path]:
    root = configured_path(config, "viirs_raw")
    pattern = config["viirs"]["raster_pattern"]

    found: dict[int, Path] = {}
    for path in sorted(root.rglob(pattern)):
        year = raster_year(path)
        if year is None:
            continue
        if requested_years is not None and year not in requested_years:
            continue
        if year in found:
            raise RuntimeError(
                f"More than one VIIRS raster found for {year}:\n"
                f"  {found[year]}\n"
                f"  {path}"
            )
        found[year] = path

    if requested_years is not None:
        missing = sorted(set(requested_years) - set(found))
        if missing:
            raise FileNotFoundError(
                "Missing VIIRS rasters for requested years: "
                + ", ".join(map(str, missing))
            )

    if not found:
        raise FileNotFoundError(
            f"No VIIRS rasters matching {pattern!r} found under {root}"
        )

    return dict(sorted(found.items()))


def load_nuts(config: dict) -> gpd.GeoDataFrame:
    path = (
        configured_path(config, "nuts_raw")
        / config["nuts"]["expected_filename"]
    )
    if not path.exists():
        raise FileNotFoundError(f"NUTS3 file not found: {path}")

    nuts = gpd.read_file(path)

    if nuts.crs is None:
        raise ValueError(f"NUTS3 file has no CRS: {path}")

    target_crs = config["nuts"]["crs"]
    if nuts.crs.to_string() != target_crs:
        nuts = nuts.to_crs(target_crs)

    id_field = config["nuts"]["id_field"]
    if id_field not in nuts.columns:
        raise ValueError(f"Missing NUTS identifier field {id_field!r}")

    return nuts



def validate_raster(path: Path, expected_crs: str) -> None:
    with rasterio.open(path) as src:
        if src.crs is None:
            raise ValueError(f"VIIRS raster has no CRS: {path}")
        if src.crs.to_string() != expected_crs:
            raise ValueError(
                f"Unexpected VIIRS CRS for {path.name}: "
                f"{src.crs}; expected {expected_crs}"
            )
        if src.count != 1:
            raise ValueError(
                f"Expected one VIIRS band in {path.name}; found {src.count}"
            )


def extract_year(
    raster_path: Path,
    year: int,
    nuts: gpd.GeoDataFrame,
    config: dict,
) -> pd.DataFrame:
    id_field = config["nuts"]["id_field"]

    include_cols = [
        column
        for column in [
            id_field,
            "CNTR_CODE",
            "NUTS_NAME",
            "NAME_LATN",
            "EU_STAT",
            "EFTA_STAT",
            "CC_STAT",
        ]
        if column in nuts.columns
    ]

    operations = [
        "radiance_mean=mean(coverage_weight=area_spherical_m2)",
        "radiance_median=median(coverage_weight=area_spherical_m2)",
        "valid_area_km2=count(coverage_weight=area_spherical_km2)",
    ]

    print(f"Processing {year}: {raster_path.name}", flush=True)

    result = exact_extract(
        str(raster_path),
        nuts,
        operations,
        include_cols=include_cols,
        output="pandas",
        progress=True,
    )

    result.insert(0, "year", year)
    result["viirs_file"] = raster_path.name
    result["viirs_product"] = config["viirs"]["product"]
    result["radiance_variable"] = config["viirs"]["variable"]
    result["radiance_units"] = config["viirs"]["units"]
    result["nuts_release"] = int(config["nuts"]["release"])
    result["zonal_method"] = config["analysis"]["area_weighting"]["method"]

    return result


def requested_years_from_args(args, config: dict) -> list[int] | None:
    if args.year:
        return sorted(set(args.year))

    if args.all_configured_years:
        first = int(config["analysis"]["first_year"])
        last = int(config["analysis"]["last_year"])
        return list(range(first, last + 1))

    # Default: process every matching raster currently available locally.
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument(
        "--year",
        type=int,
        action="append",
        help="Process a specific year. May be supplied multiple times.",
    )
    parser.add_argument(
        "--all-configured-years",
        action="store_true",
        help="Require and process every year in the configured analysis range.",
    )
    args = parser.parse_args()

    config = load_config(args.config)

    if not config["analysis"]["area_weighting"]["enabled"]:
        raise ValueError("Area weighting is disabled in pipeline configuration.")

    method = config["analysis"]["area_weighting"]["method"]
    if method != "exact_overlap_area_spherical_m2":
        raise ValueError(
            "Production extractor requires "
            "analysis.area_weighting.method=exact_overlap_area_spherical_m2; "
            f"found {method!r}."
        )

    requested = requested_years_from_args(args, config)
    rasters = discover_rasters(config, requested)
    nuts = load_nuts(config)

    expected_crs = config["nuts"]["crs"]
    frames = []

    for year, raster_path in rasters.items():
        validate_raster(raster_path, expected_crs)
        frames.append(extract_year(raster_path, year, nuts, config))

    combined = pd.concat(frames, ignore_index=True)

    id_field = config["nuts"]["id_field"]
    combined = combined.sort_values(["year", id_field]).reset_index(drop=True)

    output_dir = configured_path(config, "processed")
    output_dir.mkdir(parents=True, exist_ok=True)

    parquet_path = output_dir / config["outputs"]["radiance_by_nuts3"]
    csv_name = Path(config["outputs"]["radiance_by_nuts3"]).with_suffix(".csv").name
    csv_path = output_dir / csv_name

    combined.to_parquet(parquet_path, index=False)
    combined.to_csv(csv_path, index=False)

    print()
    print(f"Rows written: {len(combined):,}")
    print(f"Years: {sorted(combined['year'].unique().tolist())}")
    print(f"NUTS3 regions per year: {combined.groupby('year')[id_field].nunique().to_dict()}")
    print(f"Wrote: {parquet_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
