#!/usr/bin/env python3
"""Recover missing Eurostat NUTS3 population-density observations.

Published DEMO_R_D3DENS values are always retained unchanged. For configured
NUTS3 region/year combinations where published density is absent or null, a
fallback density is derived from:

    demo_r_pjanaggr3 population on 1 January / reg_area3 total land area

Only codes belonging to the configured NUTS release for that year are eligible.
Every recovered value is explicitly marked with provenance fields.
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
    categories = {dim: category_codes(payload["dimension"][dim]) for dim in dims}
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


def choose_code(payload, dim, preferred_codes, label_terms):
    if dim not in payload["dimension"]:
        return None

    d = payload["dimension"][dim]
    codes = category_codes(d)
    labels = category_labels(d)

    for code in preferred_codes:
        if code in codes:
            return code

    for term in label_terms:
        matches = [
            code for code in codes
            if term.lower() in labels.get(code, "").lower()
        ]
        if len(matches) == 1:
            return matches[0]

    return None


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
        with urllib.request.urlopen(url) as response, destination.open("wb") as out:
            shutil.copyfileobj(response, out)


def nuts_release_for_year(config: dict, year: int) -> int:
    mapping = config["nuts"]["release_by_year"]
    release = mapping.get(year, mapping.get(str(year)))
    if release is None:
        raise KeyError(f"No NUTS release configured for {year}")
    return int(release)


def load_target_codes(config: dict) -> pd.DataFrame:
    id_field = config["nuts"]["id_field"]
    root = configured_path(config, "nuts_raw")
    template = config["nuts"]["filename_template"]
    first = int(config["analysis"]["first_year"])
    last = int(config["analysis"]["last_year"])

    frames = []
    cache = {}

    for year in range(first, last + 1):
        release = nuts_release_for_year(config, year)

        if release not in cache:
            path = root / template.format(release=release)
            g = gpd.read_file(path, ignore_geometry=True)
            keep = [id_field]
            if "CNTR_CODE" in g.columns:
                keep.append("CNTR_CODE")
            cache[release] = g[keep].copy()

        y = cache[release].copy()
        y["year"] = year
        y["nuts_release"] = release
        frames.append(y)

    target = pd.concat(frames, ignore_index=True)
    if "CNTR_CODE" not in target.columns:
        target["CNTR_CODE"] = target[id_field].str[:2]

    return target[["year", id_field, "CNTR_CODE", "nuts_release"]]


def prepare_population(payload: dict) -> pd.DataFrame:
    frame = jsonstat_to_frame(payload)

    filters = {
        "sex": choose_code(payload, "sex", ["T", "TOTAL"], ["total"]),
        "age": choose_code(payload, "age", ["TOTAL", "Y_TOTAL"], ["total"]),
        "unit": choose_code(payload, "unit", ["NR"], ["number"]),
        "freq": choose_code(payload, "freq", ["A"], ["annual"]),
    }

    for dim, code in filters.items():
        if dim in frame.columns:
            if code is None:
                raise RuntimeError(
                    f"Could not select required total-population category for {dim}"
                )
            frame = frame[frame[dim] == code]

    frame["year"] = pd.to_numeric(frame["time"], errors="coerce")
    frame["fallback_population_jan1"] = pd.to_numeric(
        frame["value"], errors="coerce"
    )
    frame = frame.rename(columns={"geo": "NUTS_ID"})
    frame = frame.dropna(subset=["year", "fallback_population_jan1"])
    frame["year"] = frame["year"].astype(int)

    return frame[
        ["year", "NUTS_ID", "fallback_population_jan1", "status"]
    ].rename(columns={"status": "fallback_population_status"})


def prepare_area(payload: dict) -> pd.DataFrame:
    frame = jsonstat_to_frame(payload)

    landuse = choose_code(
        payload,
        "landuse",
        ["L0008", "TLA"],
        ["land area - total", "total land area"],
    )
    if landuse is None:
        raise RuntimeError("Could not identify total land area in reg_area3")
    frame = frame[frame["landuse"] == landuse]

    if "unit" in frame.columns:
        unit = choose_code(payload, "unit", ["KM2"], ["square kilometre"])
        if unit is not None:
            frame = frame[frame["unit"] == unit]

    if "freq" in frame.columns:
        freq = choose_code(payload, "freq", ["A"], ["annual"])
        if freq is not None:
            frame = frame[frame["freq"] == freq]

    frame["year"] = pd.to_numeric(frame["time"], errors="coerce")
    frame["fallback_land_area_km2"] = pd.to_numeric(
        frame["value"], errors="coerce"
    )
    frame = frame.rename(columns={"geo": "NUTS_ID"})
    frame = frame.dropna(subset=["year", "fallback_land_area_km2"])
    frame["year"] = frame["year"].astype(int)

    return frame[
        ["year", "NUTS_ID", "fallback_land_area_km2", "status"]
    ].rename(columns={"status": "fallback_area_status"})


def attach_area(pop: pd.DataFrame, area: pd.DataFrame) -> pd.DataFrame:
    """Attach same-year land area or nearest available area for each NUTS_ID."""
    by_code = {
        code: g.sort_values("year")
        for code, g in area.groupby("NUTS_ID")
    }

    rows = []

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
                "fallback_population_jan1": row.fallback_population_jan1,
                "fallback_population_status": row.fallback_population_status,
                "fallback_land_area_km2": float(
                    chosen["fallback_land_area_km2"]
                ),
                "fallback_area_status": chosen["fallback_area_status"],
                "fallback_area_reference_year": int(chosen["year"]),
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    id_field = config["nuts"]["id_field"]
    raw_dir = configured_path(config, "eurostat_raw")
    processed = configured_path(config, "processed")

    raw_paths = {}
    for dataset in [POP_DATASET, AREA_DATASET]:
        path = raw_dir / f"{dataset}_nuts3.json"
        download(
            EUROSTAT_API.format(dataset=dataset),
            path,
            force=args.force_download,
        )
        raw_paths[dataset] = path

    payloads = {}
    for dataset, path in raw_paths.items():
        with path.open("r", encoding="utf-8") as handle:
            payloads[dataset] = json.load(handle)

    fallback_pop = prepare_population(payloads[POP_DATASET])
    fallback_area = prepare_area(payloads[AREA_DATASET])
    fallback = attach_area(fallback_pop, fallback_area)

    fallback["fallback_population_density"] = (
        fallback["fallback_population_jan1"]
        / fallback["fallback_land_area_km2"]
    )

    first = int(config["analysis"]["first_year"])
    last = int(config["analysis"]["last_year"])
    fallback = fallback[
        fallback["year"].between(first, last, inclusive="both")
    ].copy()

    target = load_target_codes(config)

    source_path = processed / config["outputs"]["population_by_nuts3"]
    official = pd.read_parquet(source_path).copy()

    official["population_density_official"] = official["population_density"]
    official["population_density_source"] = "demo_r_d3dens"

    # Restrict the official rows to the configured NUTS geometry. This avoids
    # retaining obsolete codes that coexist in Eurostat historical datasets.
    official_target = target.merge(
        official,
        on=["year", id_field, "nuts_release"],
        how="left",
        validate="one_to_one",
    )

    result = official_target.merge(
        fallback[
            [
                "year",
                id_field,
                "fallback_population_density",
                "fallback_population_jan1",
                "fallback_population_status",
                "fallback_land_area_km2",
                "fallback_area_status",
                "fallback_area_reference_year",
            ]
        ],
        on=["year", id_field],
        how="left",
        validate="one_to_one",
    )

    needs_fill = result["population_density"].isna()

    allowed_countries = set(
        config["population"].get("fallback_recovery_countries", [])
    )
    if allowed_countries:
        eligible_country = result["CNTR_CODE"].isin(allowed_countries)
    else:
        eligible_country = pd.Series(False, index=result.index)

    can_fill = (
        needs_fill
        & result["fallback_population_density"].notna()
        & eligible_country
    )

    result.loc[can_fill, "population_density"] = result.loc[
        can_fill, "fallback_population_density"
    ]
    result.loc[can_fill, "population_density_source"] = (
        "demo_r_pjanaggr3_div_reg_area3_land"
    )

    result["matches_nuts_release"] = True
    result["population_density_recovered"] = can_fill
    result["population_density_missing_after_recovery"] = result[
        "population_density"
    ].isna()

    if "unit" in result.columns:
        result.loc[can_fill & result["unit"].isna(), "unit"] = (
            config["population"]["expected_units"]
        )
    if "freq" in result.columns:
        result.loc[can_fill & result["freq"].isna(), "freq"] = (
            config["population"]["frequency"]
        )

    recovered = result[result["population_density_recovered"]].copy()
    missing = result[result["population_density_missing_after_recovery"]].copy()

    result = result.sort_values(["year", id_field]).reset_index(drop=True)

    result.to_parquet(source_path, index=False)
    result.to_csv(source_path.with_suffix(".csv"), index=False)

    print("\nRecovered population-density observations:")
    if recovered.empty:
        print("  none")
    else:
        print(
            recovered.groupby(["year", "CNTR_CODE"])
            .size()
            .rename("n_recovered")
            .reset_index()
            .to_string(index=False)
        )

    print("\nRemaining missing target NUTS3 densities after recovery:")
    if missing.empty:
        print("  none")
    else:
        print(
            missing.groupby(["year", "CNTR_CODE"])
            .size()
            .rename("n_missing")
            .reset_index()
            .to_string(index=False)
        )

    print(f"\nRows written: {len(result):,}")
    print(f"Wrote: {source_path}")
    print(f"Wrote: {source_path.with_suffix('.csv')}")


if __name__ == "__main__":
    main()
