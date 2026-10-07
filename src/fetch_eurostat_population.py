#!/usr/bin/env python3
"""Download and normalise Eurostat NUTS3 population-density data.

The raw source is the official Eurostat Statistics API:

  DEMO_R_D3DENS
  Population density by NUTS 3 region

By default the script uses wget for a reproducible, inspectable raw download.
The downloaded JSON-stat response is retained under data/raw/eurostat and
normalised to one row per year x NUTS3 identifier.
"""

from __future__ import annotations

import argparse
import itertools
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path

import geopandas as gpd
import pandas as pd

from project_config import configured_path, load_config


EUROSTAT_API = (
    "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/"
    "{dataset}?lang=EN&geoLevel=nuts3"
)


def category_codes(dimension: dict) -> list[str]:
    """Return category codes in JSON-stat positional order."""
    index = dimension["category"]["index"]

    if isinstance(index, dict):
        return [
            code
            for code, _ in sorted(index.items(), key=lambda item: item[1])
        ]

    # JSON-stat also permits an ordered array of category codes.
    return list(index)


def download(url: str, destination: Path, force: bool = False) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists() and not force:
        print(f"Raw Eurostat file already exists: {destination}")
        print("Use --force-download to replace it.")
        return

    if shutil.which("wget"):
        command = [
            "wget",
            "--continue",
            "--output-document",
            str(destination),
            url,
        ]
        print("+", " ".join(command), flush=True)
        subprocess.run(command, check=True)
    else:
        print("wget not found; falling back to Python urllib.", flush=True)
        with urllib.request.urlopen(url) as response, destination.open("wb") as out:
            shutil.copyfileobj(response, out)


def jsonstat_to_frame(payload: dict) -> pd.DataFrame:
    dims = payload["id"]
    sizes = payload["size"]

    if len(dims) != len(sizes):
        raise ValueError("Malformed JSON-stat dataset: id/size lengths differ.")

    categories = {
        dim: category_codes(payload["dimension"][dim])
        for dim in dims
    }

    expected = 1
    for size in sizes:
        expected *= int(size)

    value_obj = payload.get("value", {})
    status_obj = payload.get("status", {})

    def positional(obj, position):
        if isinstance(obj, list):
            return obj[position] if position < len(obj) else None
        if isinstance(obj, dict):
            return obj.get(str(position), obj.get(position))
        return None

    rows = []
    for position, combination in enumerate(
        itertools.product(*(categories[dim] for dim in dims))
    ):
        value = positional(value_obj, position)
        status = positional(status_obj, position)

        # Sparse JSON-stat responses omit absent cells.
        if value is None and status is None:
            continue

        row = dict(zip(dims, combination))
        row["population_density"] = value
        row["status"] = status
        rows.append(row)

    if expected == 0:
        raise ValueError("Eurostat response contains no cells.")

    return pd.DataFrame(rows)


def load_nuts_ids(config: dict) -> set[str]:
    nuts_path = (
        configured_path(config, "nuts_raw")
        / config["nuts"]["expected_filename"]
    )
    nuts = gpd.read_file(nuts_path, ignore_geometry=True)
    id_field = config["nuts"]["id_field"]
    return set(nuts[id_field].astype(str))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Replace the retained raw Eurostat JSON response.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    dataset = config["population"]["dataset"]
    dataset_lower = dataset.lower()

    raw_dir = configured_path(config, "eurostat_raw")
    raw_path = raw_dir / f"{dataset_lower}_nuts3.json"

    url = EUROSTAT_API.format(dataset=dataset)
    download(url, raw_path, force=args.force_download)

    with raw_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    label = payload.get("label", "")
    updated = payload.get("updated", "")
    print(f"Dataset: {dataset}")
    print(f"Label: {label}")
    print(f"Updated: {updated}")

    frame = jsonstat_to_frame(payload)

    # Defensive filtering even though the API query requests NUTS3.
    if "freq" in frame.columns:
        frame = frame[frame["freq"] == config["population"]["frequency"]]

    first_year = int(config["analysis"]["first_year"])
    last_year = int(config["analysis"]["last_year"])

    frame["year"] = pd.to_numeric(frame["time"], errors="coerce")
    frame = frame[
        frame["year"].between(first_year, last_year, inclusive="both")
    ].copy()
    frame["year"] = frame["year"].astype(int)

    frame["population_density"] = pd.to_numeric(
        frame["population_density"],
        errors="coerce",
    )

    frame = frame.rename(columns={"geo": config["nuts"]["id_field"]})
    id_field = config["nuts"]["id_field"]

    nuts2024_ids = load_nuts_ids(config)
    frame["matches_nuts2024"] = frame[id_field].isin(nuts2024_ids)

    keep_columns = [
        "year",
        id_field,
        "population_density",
        "status",
        "unit",
        "freq",
        "matches_nuts2024",
    ]
    keep_columns = [c for c in keep_columns if c in frame.columns]
    frame = frame[keep_columns].sort_values(["year", id_field]).reset_index(drop=True)

    output_dir = configured_path(config, "processed")
    output_dir.mkdir(parents=True, exist_ok=True)

    parquet_path = output_dir / config["outputs"]["population_by_nuts3"]
    csv_path = parquet_path.with_suffix(".csv")

    frame.to_parquet(parquet_path, index=False)
    frame.to_csv(csv_path, index=False)

    print()
    print(f"Rows written: {len(frame):,}")
    print(f"Years: {sorted(frame['year'].unique().tolist())}")

    coverage = (
        frame.groupby("year")
        .agg(
            eurostat_rows=(id_field, "size"),
            non_null_density=("population_density", "count"),
            matching_nuts2024=("matches_nuts2024", "sum"),
        )
    )
    print()
    print("Coverage by year:")
    print(coverage.to_string())

    unmatched = (
        frame.loc[~frame["matches_nuts2024"], id_field]
        .drop_duplicates()
        .sort_values()
    )
    print()
    print(f"Distinct Eurostat NUTS3 codes not in NUTS 2024 geometry: {len(unmatched):,}")
    if len(unmatched):
        print(" ".join(unmatched.tolist()))

    print()
    print(f"Wrote: {raw_path}")
    print(f"Wrote: {parquet_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
