#!/usr/bin/env python3
"""Extract annual VIIRS mean radiance on the configured Eurostat-aligned NUTS geography.

For each NUTS3 polygon, mean radiance is calculated with exact polygon overlap
and spherical physical-area weighting:

    mean(coverage_weight=area_spherical_m2)

The NUTS release is selected from config/pipeline.yaml so that VIIRS radiance
is aggregated on the same retrospectively harmonized NUTS classification used
by the current Eurostat population-density time series.
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


def nuts_release_for_year(config: dict, year: int) -> int:
    mapping = config["nuts"]["release_by_year"]
    release = mapping.get(year, mapping.get(str(year)))
    if release is None:
        raise KeyError(f"No NUTS release configured for analysis year {year}")
    return int(release)


def nuts_path_for_release(config: dict, release: int) -> Path:
    filename = config["nuts"]["filename_template"].format(release=release)
    return configured_path(config, "nuts_raw") / filename


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


def load_nuts(config: dict, year: int) -> tuple[gpd.GeoDataFrame, int]:
    release = nuts_release_for_year(config, year)
    path = nuts_path_for_release(config, release)

    if not path.exists():
        raise FileNotFoundError(
            f"NUTS {release} Level-3 file not found: {path}. "
            "Run: python src/run_pipeline.py nuts"
        )

    nuts = gpd.read_file(path)

    if nuts.crs is None:
        raise ValueError(f"NUTS file has no CRS: {path}")

    target_crs = config["nuts"]["crs"]
    if nuts.crs.to_string() != target_crs:
        nuts = nuts.to_crs(target_crs)

    id_field = config["nuts"]["id_field"]
    if id_field not in nuts.columns:
        raise ValueError(f"Missing NUTS identifier field {id_field!r}")

    return nuts, release


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
    nuts_release: int,
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
        "valid_area_km2=count(coverage_weight=area_spherical_km2)",
    ]

    print(
        f"Processing {year} with NUTS {nuts_release}: {raster_path.name}",
        flush=True,
    )

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
    result["nuts_release"] = nuts_release
    result["zonal_method"] = config["analysis"]["area_weighting"]["method"]

    return result


def requested_years_from_args(args, config: dict) -> list[int] | None:
    if args.year:
        return sorted(set(args.year))

    if args.all_configured_years:
        first = int(config["analysis"]["first_year"])
        last = int(config["analysis"]["last_year"])
        return list(range(first, last + 1))

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
    expected_crs = config["nuts"]["crs"]
    frames = []

    for year, raster_path in rasters.items():
        validate_raster(raster_path, expected_crs)
        nuts, release = load_nuts(config, year)
        frames.append(extract_year(raster_path, year, nuts, release, config))

    combined = pd.concat(frames, ignore_index=True)

    id_field = config["nuts"]["id_field"]
    combined = combined.sort_values(["year", id_field]).reset_index(drop=True)

    output_dir = configured_path(config, "processed")
    output_dir.mkdir(parents=True, exist_ok=True)

    parquet_path = output_dir / config["outputs"]["radiance_by_nuts3"]
    csv_path = parquet_path.with_suffix(".csv")

    combined.to_parquet(parquet_path, index=False)
    combined.to_csv(csv_path, index=False)

    print()
    print(f"Rows written: {len(combined):,}")
    print(f"Years: {sorted(combined['year'].unique().tolist())}")
    print(
        "NUTS3 regions per year: "
        f"{combined.groupby('year')[id_field].nunique().to_dict()}"
    )
    print(
        "NUTS release per year: "
        f"{combined.groupby('year')['nuts_release'].first().to_dict()}"
    )
    print(f"Wrote: {parquet_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
