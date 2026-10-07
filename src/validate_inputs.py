#!/usr/bin/env python3
"""Validate local inputs before running the VIIRS/NUTS3 workflow."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import geopandas as gpd
import rasterio

from project_config import configured_path, load_config


YEAR_RE = re.compile(r"(20\d{2})")


def raster_year(path: Path) -> int | None:
    match = YEAR_RE.search(path.name)
    return int(match.group(1)) if match else None


def nuts_release_for_year(config: dict, year: int) -> int:
    mapping = config["nuts"]["release_by_year"]
    release = mapping.get(year, mapping.get(str(year)))
    if release is None:
        raise KeyError(f"No NUTS release configured for {year}")
    return int(release)


def validate_viirs(config: dict) -> list[Path]:
    root = configured_path(config, "viirs_raw")
    pattern = config["viirs"]["raster_pattern"]
    rasters = sorted(root.rglob(pattern))

    print(f"VIIRS directory: {root}")
    print(f"VIIRS pattern:   {pattern}")
    print(f"VIIRS rasters:   {len(rasters)}")

    for path in rasters:
        year = raster_year(path)
        with rasterio.open(path) as src:
            print(
                f"  {year or '????'}  {path.name}  "
                f"{src.width}x{src.height}  {src.crs}  bands={src.count}"
            )

    return rasters


def validate_nuts(config: dict) -> None:
    root = configured_path(config, "nuts_raw")
    template = config["nuts"]["filename_template"]
    id_field = config["nuts"]["id_field"]

    print("NUTS releases:")
    for release in config["nuts"]["releases"]:
        path = root / template.format(release=release)
        if not path.exists():
            print(f"  {release}: MISSING  {path}")
            continue

        nuts = gpd.read_file(path)
        if id_field not in nuts.columns:
            raise ValueError(f"NUTS ID field {id_field!r} not present in {path}")

        print(
            f"  {release}: {len(nuts):,} features; "
            f"CRS={nuts.crs}; {path.name}"
        )


def validate_year_coverage(config: dict, rasters: list[Path]) -> None:
    first = int(config["analysis"]["first_year"])
    last = int(config["analysis"]["last_year"])
    wanted = set(range(first, last + 1))
    found = {year for path in rasters if (year := raster_year(path)) is not None}

    missing = sorted(wanted - found)
    present = sorted(wanted & found)

    print(f"Analysis years:   {first}-{last}")
    print(f"Years present:    {present or 'none'}")
    print(f"Years missing:    {missing or 'none'}")
    print("Year -> NUTS release:")
    for year in range(first, last + 1):
        print(f"  {year}: NUTS {nuts_release_for_year(config, year)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    rasters = validate_viirs(config)
    validate_nuts(config)
    validate_year_coverage(config, rasters)


if __name__ == "__main__":
    main()
