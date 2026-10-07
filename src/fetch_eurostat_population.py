#!/usr/bin/env python3
"""Download and normalise Eurostat NUTS3 population-density data.

The raw source is Eurostat DEMO_R_D3DENS. Each observation is checked against
the NUTS release officially applicable to its reference year.
"""

from __future__ import annotations

import argparse
import itertools
import json
import shutil
import subprocess
import urllib.request

import geopandas as gpd
import pandas as pd

from project_config import configured_path, load_config


EUROSTAT_API = (
    "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/"
    "{dataset}?lang=EN&geoLevel=nuts3"
)


def category_codes(dimension: dict) -> list[str]:
    index = dimension["category"]["index"]
    if isinstance(index, dict):
        return [code for code, _ in sorted(index.items(), key=lambda item: item[1])]
    return list(index)


def download(url: str, destination, force: bool = False) -> None:
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
    categories = {
        dim: category_codes(payload["dimension"][dim])
        for dim in dims
    }

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
        if value is None and status is None:
            continue
        row = dict(zip(dims, combination))
        row["population_density"] = value
        row["status"] = status
        rows.append(row)

    return pd.DataFrame(rows)


def nuts_release_for_year(config: dict, year: int) -> int:
    mapping = config["nuts"]["release_by_year"]
    release = mapping.get(year, mapping.get(str(year)))
    if release is None:
        raise KeyError(f"No NUTS release configured for {year}")
    return int(release)


def load_nuts_ids(config: dict, release: int) -> set[str]:
    filename = config["nuts"]["filename_template"].format(release=release)
    path = configured_path(config, "nuts_raw") / filename

    if not path.exists():
        raise FileNotFoundError(
            f"NUTS {release} file not found: {path}. "
            "Run: python src/run_pipeline.py nuts"
        )

    nuts = gpd.read_file(path, ignore_geometry=True)
    id_field = config["nuts"]["id_field"]
    return set(nuts[id_field].astype(str))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    dataset = config["population"]["dataset"]
    raw_dir = configured_path(config, "eurostat_raw")
    raw_path = raw_dir / f"{dataset.lower()}_nuts3.json"

    download(
        EUROSTAT_API.format(dataset=dataset),
        raw_path,
        force=args.force_download,
    )

    with raw_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    print(f"Dataset: {dataset}")
    print(f"Label: {payload.get('label', '')}")
    print(f"Updated: {payload.get('updated', '')}")

    frame = jsonstat_to_frame(payload)

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
        frame["population_density"], errors="coerce"
    )

    id_field = config["nuts"]["id_field"]
    frame = frame.rename(columns={"geo": id_field})

    id_sets = {
        int(release): load_nuts_ids(config, int(release))
        for release in config["nuts"]["releases"]
    }

    frame["nuts_release"] = frame["year"].map(
        lambda year: nuts_release_for_year(config, int(year))
    )
    frame["matches_nuts_release"] = [
        nuts_id in id_sets[int(release)]
        for nuts_id, release in zip(frame[id_field], frame["nuts_release"])
    ]

    keep_columns = [
        "year",
        id_field,
        "population_density",
        "status",
        "unit",
        "freq",
        "nuts_release",
        "matches_nuts_release",
    ]
    keep_columns = [column for column in keep_columns if column in frame.columns]
    frame = frame[keep_columns].sort_values(["year", id_field]).reset_index(drop=True)

    output_dir = configured_path(config, "processed")
    output_dir.mkdir(parents=True, exist_ok=True)
    parquet_path = output_dir / config["outputs"]["population_by_nuts3"]
    csv_path = parquet_path.with_suffix(".csv")

    frame.to_parquet(parquet_path, index=False)
    frame.to_csv(csv_path, index=False)

    coverage = (
        frame.groupby(["year", "nuts_release"])
        .agg(
            eurostat_rows=(id_field, "size"),
            non_null_density=("population_density", "count"),
            matching_year_nuts=("matches_nuts_release", "sum"),
        )
    )

    print()
    print(f"Rows written: {len(frame):,}")
    print("Coverage by year and applicable NUTS release:")
    print(coverage.to_string())

    print()
    print(f"Wrote: {raw_path}")
    print(f"Wrote: {parquet_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
