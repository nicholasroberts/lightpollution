#!/usr/bin/env python3
"""Test an official-Eurostat fallback for missing NUTS3 population density.

This is a diagnostic stage only. It does NOT alter the production population
table.

Candidate fallback:
    population on 1 January (demo_r_pjanaggr3) / total land area (reg_area3)

The published demo_r_d3dens indicator uses annual-average population divided
by land area, so the reconstructed values are validated against published
density wherever both are available before any production use is considered.
"""

from __future__ import annotations

import argparse
import itertools
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from project_config import configured_path, load_config


EUROSTAT_API = (
    "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/"
    "{dataset}?lang=EN&geoLevel=nuts3"
)

POP_DATASET = "demo_r_pjanaggr3"
AREA_DATASET = "reg_area3"


def category_codes(dimension: dict) -> list[str]:
    index = dimension["category"]["index"]
    if isinstance(index, dict):
        return [code for code, _ in sorted(index.items(), key=lambda item: item[1])]
    return list(index)


def category_labels(dimension: dict) -> dict[str, str]:
    labels = dimension["category"].get("label", {})
    return {str(k): str(v) for k, v in labels.items()} if isinstance(labels, dict) else {}


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
        row["value"] = value
        row["status"] = status
        rows.append(row)
    return pd.DataFrame(rows)


def download(url: str, destination: Path, force: bool = False) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not force:
        print(f"Already present: {destination}")
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


def describe_dimensions(payload: dict, dataset: str) -> None:
    print(f"\n{dataset} dimensions:")
    for dim in payload["id"]:
        d = payload["dimension"][dim]
        labels = category_labels(d)
        codes = category_codes(d)
        preview = ", ".join(
            f"{code}={labels.get(code, '')}" for code in codes[:12]
        )
        if len(codes) > 12:
            preview += ", ..."
        print(f"  {dim}: {len(codes)} codes; {preview}")


def choose_code(
    payload: dict,
    dim: str,
    preferred_codes: list[str],
    label_terms: list[str],
) -> str | None:
    if dim not in payload["dimension"]:
        return None

    d = payload["dimension"][dim]
    codes = category_codes(d)
    labels = category_labels(d)

    for code in preferred_codes:
        if code in codes:
            return code

    for term in label_terms:
        term = term.lower()
        matches = [
            code
            for code in codes
            if term in labels.get(code, "").lower()
        ]
        if len(matches) == 1:
            return matches[0]

    return None


def prepare_population(payload: dict) -> pd.DataFrame:
    frame = jsonstat_to_frame(payload)

    filters = {}
    if "sex" in frame.columns:
        filters["sex"] = choose_code(
            payload, "sex", ["T", "TOTAL"], ["total"]
        )
    if "age" in frame.columns:
        filters["age"] = choose_code(
            payload, "age", ["TOTAL", "Y_TOTAL"], ["total"]
        )
    if "unit" in frame.columns:
        filters["unit"] = choose_code(
            payload, "unit", ["NR"], ["number"]
        )
    if "freq" in frame.columns:
        filters["freq"] = choose_code(
            payload, "freq", ["A"], ["annual"]
        )

    print("\nSelected population dimensions:")
    for dim, code in filters.items():
        print(f"  {dim}={code}")
        if code is None:
            raise RuntimeError(
                f"Could not identify required total population category for {dim}"
            )
        frame = frame[frame[dim] == code]

    frame["year"] = pd.to_numeric(frame["time"], errors="coerce")
    frame["population_jan1"] = pd.to_numeric(frame["value"], errors="coerce")
    frame = frame.rename(columns={"geo": "NUTS_ID"})

    return frame[
        ["year", "NUTS_ID", "population_jan1", "status"]
    ].dropna(subset=["year"]).assign(year=lambda x: x["year"].astype(int))


def prepare_area(payload: dict) -> pd.DataFrame:
    frame = jsonstat_to_frame(payload)

    print("\nSelecting total land area (TLA) from reg_area3.")

    # Find the dimension whose labels/codes contain Total Land Area.
    selector_dim = None
    selector_code = None
    for dim in payload["id"]:
        if dim in {"geo", "time", "freq", "unit"}:
            continue
        code = choose_code(
            payload,
            dim,
            ["TLA", "LAND", "LAND_AREA"],
            ["total land area", "land area"],
        )
        if code is not None:
            selector_dim = dim
            selector_code = code
            break

    if selector_dim is None:
        raise RuntimeError(
            "Could not identify the Total Land Area category in reg_area3. "
            "See the printed dimension listing above."
        )

    print(f"  {selector_dim}={selector_code}")
    frame = frame[frame[selector_dim] == selector_code]

    if "unit" in frame.columns:
        unit = choose_code(
            payload, "unit", ["KM2", "KM2_T"], ["square kilometre", "km"]
        )
        if unit is not None:
            print(f"  unit={unit}")
            frame = frame[frame["unit"] == unit]

    if "freq" in frame.columns:
        freq = choose_code(payload, "freq", ["A"], ["annual"])
        if freq is not None:
            frame = frame[frame["freq"] == freq]

    frame["year"] = pd.to_numeric(frame["time"], errors="coerce")
    frame["land_area_km2"] = pd.to_numeric(frame["value"], errors="coerce")
    frame = frame.rename(columns={"geo": "NUTS_ID"})
    frame = frame.dropna(subset=["year", "land_area_km2"])
    frame["year"] = frame["year"].astype(int)

    return frame[["year", "NUTS_ID", "land_area_km2", "status"]]


def select_area_for_population(
    pop: pd.DataFrame,
    area: pd.DataFrame,
) -> pd.DataFrame:
    """Attach the best available area record for each population NUTS_ID/year.

    Prefer the same reference year. If reg_area3 stores one area per NUTS
    classification rather than every calendar year, use the most recent
    available area record not later than the population year; otherwise the
    nearest available record for that exact NUTS_ID.
    """
    rows = []
    by_code = {
        code: g.sort_values("year")
        for code, g in area.groupby("NUTS_ID")
    }

    for row in pop.itertuples(index=False):
        candidates = by_code.get(row.NUTS_ID)
        if candidates is None or candidates.empty:
            continue

        same = candidates[candidates["year"] == row.year]
        if not same.empty:
            chosen = same.iloc[-1]
        else:
            earlier = candidates[candidates["year"] <= row.year]
            if not earlier.empty:
                chosen = earlier.iloc[-1]
            else:
                chosen = candidates.iloc[0]

        rows.append(
            {
                "year": row.year,
                "NUTS_ID": row.NUTS_ID,
                "population_jan1": row.population_jan1,
                "population_status": row.status,
                "area_reference_year": int(chosen["year"]),
                "land_area_km2": float(chosen["land_area_km2"]),
                "area_status": chosen["status"],
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--country", default="NL")
    parser.add_argument(
        "--years",
        nargs="+",
        type=int,
        default=[2021, 2022, 2023, 2024],
    )
    args = parser.parse_args()

    config = load_config(args.config)
    raw_dir = configured_path(config, "eurostat_raw")

    paths = {}
    for dataset in [POP_DATASET, AREA_DATASET]:
        path = raw_dir / f"{dataset}_nuts3.json"
        download(
            EUROSTAT_API.format(dataset=dataset),
            path,
            force=args.force_download,
        )
        paths[dataset] = path

    payloads = {}
    for dataset, path in paths.items():
        with path.open("r", encoding="utf-8") as handle:
            payloads[dataset] = json.load(handle)
        describe_dimensions(payloads[dataset], dataset)

    pop = prepare_population(payloads[POP_DATASET])
    area = prepare_area(payloads[AREA_DATASET])
    derived = select_area_for_population(pop, area)

    derived["derived_density_jan1"] = (
        derived["population_jan1"] / derived["land_area_km2"]
    )

    official_path = (
        configured_path(config, "processed")
        / config["outputs"]["population_by_nuts3"]
    )
    official = pd.read_parquet(official_path)[
        ["year", "NUTS_ID", "population_density", "status"]
    ].rename(columns={"status": "density_status"})

    comparison = derived.merge(
        official,
        on=["year", "NUTS_ID"],
        how="left",
        validate="one_to_one",
    )

    comparison["density_difference"] = (
        comparison["derived_density_jan1"]
        - comparison["population_density"]
    )
    comparison["density_percent_difference"] = (
        100.0
        * comparison["density_difference"]
        / comparison["population_density"]
    )

    country_prefix = args.country.upper()
    c = comparison[
        comparison["NUTS_ID"].str.startswith(country_prefix)
        & comparison["year"].isin(args.years)
    ].copy()

    output_dir = configured_path(config, "processed")
    output = output_dir / (
        f"population_density_fallback_validation_{country_prefix.lower()}.csv"
    )
    c.to_csv(output, index=False)

    print(
        f"\n{country_prefix}: comparison of Jan-1/land-area density "
        "with official demo_r_d3dens"
    )
    print("=" * 78)

    for year in args.years:
        y = c[c["year"] == year].copy()
        both = y[y["population_density"].notna()].copy()
        missing = y[y["population_density"].isna()].copy()

        print(f"\n{year}")
        print(f"  derived density rows: {len(y)}")
        print(f"  rows also having official density: {len(both)}")
        print(f"  candidate recovered missing density rows: {len(missing)}")

        if not both.empty:
            abs_pct = both["density_percent_difference"].abs()
            print(
                "  validation |% difference|: "
                f"median={abs_pct.median():.3f}%  "
                f"mean={abs_pct.mean():.3f}%  "
                f"max={abs_pct.max():.3f}%"
            )
            corr = both[
                ["derived_density_jan1", "population_density"]
            ].corr().iloc[0, 1]
            print(f"  correlation with official density: {corr:.6f}")

        if not missing.empty:
            print("  recovered codes:")
            for row in missing.sort_values("NUTS_ID").itertuples():
                print(
                    f"    {row.NUTS_ID}: "
                    f"population={row.population_jan1:,.0f}, "
                    f"land_area={row.land_area_km2:,.3f} km2, "
                    f"derived_density={row.derived_density_jan1:,.3f}"
                )

    print(f"\nWrote: {output}")
    print(
        "\nDiagnostic only: production density values have NOT been modified."
    )


if __name__ == "__main__":
    main()
