#!/usr/bin/env python3
"""Download and prepare annual Eurostat NUTS3 socioeconomic metrics.

Primary source:
    nama_10r_3gdp -- GDP at current market prices by NUTS3 region

The script retains several Eurostat GDP representations and derives GDP
density using Eurostat total land area (reg_area3). The primary GDP covariate
for mixed modelling is GDP per capita in PPS, which improves cross-country
comparability.

GDP density is retained for diagnostics and mapping, but is NOT used alongside
population density and GDP per capita in the same regression because, by
construction, it is approximately their product and would be collinear.
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
import numpy as np
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


def category_labels(dimension: dict) -> dict[str, str]:
    labels = dimension["category"].get("label", {})
    if isinstance(labels, dict):
        return {str(k): str(v) for k, v in labels.items()}
    return {}


def jsonstat_to_frame(payload: dict) -> pd.DataFrame:
    dims = payload["id"]
    categories = {
        dim: category_codes(payload["dimension"][dim])
        for dim in dims
    }
    values = payload.get("value", {})
    statuses = payload.get("status", {})

    def positional(obj, position):
        if isinstance(obj, list):
            return obj[position] if position < len(obj) else None
        if isinstance(obj, dict):
            return obj.get(str(position), obj.get(position))
        return None

    rows = []
    for position, combination in enumerate(
        itertools.product(*(categories[d] for d in dims))
    ):
        value = positional(values, position)
        status = positional(statuses, position)
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

    cache = {}
    frames = []
    for year in range(first, last + 1):
        release = nuts_release_for_year(config, year)
        if release not in cache:
            path = root / template.format(release=release)
            if not path.exists():
                raise FileNotFoundError(
                    f"Missing NUTS {release}: {path}. "
                    "Run: python src/run_pipeline.py nuts"
                )
            g = gpd.read_file(path, ignore_geometry=True)
            cols = [id_field]
            if "CNTR_CODE" in g.columns:
                cols.append("CNTR_CODE")
            cache[release] = g[cols].copy()

        y = cache[release].copy()
        y["year"] = year
        y["nuts_release"] = release
        frames.append(y)

    target = pd.concat(frames, ignore_index=True)
    if "CNTR_CODE" not in target.columns:
        target["CNTR_CODE"] = target[id_field].str[:2]

    return target[["year", id_field, "CNTR_CODE", "nuts_release"]]


def prepare_gdp(payload: dict, config: dict) -> pd.DataFrame:
    frame = jsonstat_to_frame(payload)
    socioeconomic = config["socioeconomic"]

    if "freq" in frame.columns:
        frame = frame[frame["freq"] == socioeconomic["frequency"]].copy()

    frame["year"] = pd.to_numeric(frame["time"], errors="coerce")
    frame["value"] = pd.to_numeric(frame["value"], errors="coerce")

    first = int(config["analysis"]["first_year"])
    last = int(config["analysis"]["last_year"])
    frame = frame[
        frame["year"].between(first, last, inclusive="both")
    ].copy()
    frame["year"] = frame["year"].astype(int)
    frame = frame.rename(columns={"geo": "NUTS_ID"})

    wanted = socioeconomic["gdp_units"]
    available = set(frame["unit"].dropna().astype(str))
    missing_units = [
        code for code in wanted.values()
        if code not in available
    ]
    if missing_units:
        labels = category_labels(payload["dimension"]["unit"])
        detail = {code: labels.get(code, "") for code in sorted(available)}
        raise RuntimeError(
            f"Required GDP units not present: {missing_units}. "
            f"Available units: {detail}"
        )

    unit_to_metric = {code: metric for metric, code in wanted.items()}
    subset = frame[frame["unit"].isin(unit_to_metric)].copy()
    subset["metric"] = subset["unit"].map(unit_to_metric)

    values = subset.pivot_table(
        index=["year", "NUTS_ID"],
        columns="metric",
        values="value",
        aggfunc="first",
    ).reset_index()
    values.columns.name = None

    status = subset.pivot_table(
        index=["year", "NUTS_ID"],
        columns="metric",
        values="status",
        aggfunc="first",
    ).reset_index()
    status.columns = [
        "year" if c == "year" else
        "NUTS_ID" if c == "NUTS_ID" else
        f"{c}_status"
        for c in status.columns
    ]

    return values.merge(
        status,
        on=["year", "NUTS_ID"],
        how="left",
        validate="one_to_one",
    )


def prepare_land_area(payload: dict, config: dict) -> pd.DataFrame:
    frame = jsonstat_to_frame(payload)

    labels = category_labels(payload["dimension"]["landuse"])
    landuse_codes = category_codes(payload["dimension"]["landuse"])

    preferred = ["L0008", "TLA"]
    landuse = next((c for c in preferred if c in landuse_codes), None)
    if landuse is None:
        candidates = [
            c for c in landuse_codes
            if "total land area" in labels.get(c, "").lower()
            or "land area - total" in labels.get(c, "").lower()
        ]
        if len(candidates) != 1:
            raise RuntimeError("Could not uniquely identify total land area.")
        landuse = candidates[0]

    frame = frame[frame["landuse"] == landuse].copy()
    if "unit" in frame.columns and "KM2" in set(frame["unit"]):
        frame = frame[frame["unit"] == "KM2"]
    if "freq" in frame.columns and "A" in set(frame["freq"]):
        frame = frame[frame["freq"] == "A"]

    frame["year"] = pd.to_numeric(frame["time"], errors="coerce")
    frame["land_area_km2"] = pd.to_numeric(frame["value"], errors="coerce")
    frame = frame.rename(columns={"geo": "NUTS_ID"})
    frame = frame.dropna(subset=["year", "land_area_km2"]).copy()
    frame["year"] = frame["year"].astype(int)

    return frame[
        ["year", "NUTS_ID", "land_area_km2", "status"]
    ].rename(columns={"status": "land_area_status"})


def attach_nearest_area(
    target: pd.DataFrame,
    area: pd.DataFrame,
) -> pd.DataFrame:
    by_code = {
        code: g.sort_values("year")
        for code, g in area.groupby("NUTS_ID")
    }

    rows = []
    for row in target.itertuples(index=False):
        candidates = by_code.get(row.NUTS_ID)
        chosen = None
        if candidates is not None and not candidates.empty:
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
                "land_area_km2": (
                    float(chosen["land_area_km2"])
                    if chosen is not None else np.nan
                ),
                "land_area_status": (
                    chosen["land_area_status"]
                    if chosen is not None else None
                ),
                "land_area_reference_year": (
                    int(chosen["year"])
                    if chosen is not None else np.nan
                ),
            }
        )

    return pd.DataFrame(rows)


def add_derived_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()

    if {"gdp_total_eur_million", "land_area_km2"} <= set(out.columns):
        out["gdp_density_eur_per_km2"] = (
            out["gdp_total_eur_million"] * 1_000_000
            / out["land_area_km2"]
        )

    if {"gdp_total_pps_million", "land_area_km2"} <= set(out.columns):
        out["gdp_density_pps_per_km2"] = (
            out["gdp_total_pps_million"] * 1_000_000
            / out["land_area_km2"]
        )

    positive_metrics = [
        "gdp_total_eur_million",
        "gdp_total_pps_million",
        "gdp_per_capita_eur",
        "gdp_per_capita_pps",
        "gdp_density_eur_per_km2",
        "gdp_density_pps_per_km2",
    ]
    for col in positive_metrics:
        if col in out.columns:
            log_col = f"log10_{col}"
            out[log_col] = np.nan
            ok = np.isfinite(out[col]) & (out[col] > 0)
            out.loc[ok, log_col] = np.log10(out.loc[ok, col])

    # Descriptive year-on-year changes; not automatically included in the
    # maximal model because they are derived temporal quantities.
    for col in ["gdp_per_capita_pps", "gdp_density_pps_per_km2"]:
        if col in out.columns:
            out = out.sort_values(["NUTS_ID", "year"])
            out[f"{col}_yoy_pct"] = (
                out.groupby("NUTS_ID")[col].pct_change(fill_method=None) * 100
            )

    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    socioeconomic = config["socioeconomic"]
    raw_dir = configured_path(config, "eurostat_raw")
    processed = configured_path(config, "processed")

    datasets = [
        socioeconomic["gdp_dataset"],
        socioeconomic["area_dataset"],
    ]
    raw_paths = {}
    for dataset in datasets:
        path = raw_dir / f"{dataset.lower()}_nuts3.json"
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

    gdp = prepare_gdp(
        payloads[socioeconomic["gdp_dataset"]],
        config,
    )
    area = prepare_land_area(
        payloads[socioeconomic["area_dataset"]],
        config,
    )
    target = load_target_codes(config)

    target_gdp = target.merge(
        gdp,
        on=["year", "NUTS_ID"],
        how="left",
        validate="one_to_one",
    )
    target_area = attach_nearest_area(
        target[["year", "NUTS_ID"]],
        area,
    )
    result = target_gdp.merge(
        target_area,
        on=["year", "NUTS_ID"],
        how="left",
        validate="one_to_one",
    )
    result = add_derived_metrics(result)
    result["gdp_source"] = socioeconomic["gdp_dataset"]
    result["land_area_source"] = socioeconomic["area_dataset"]

    output = processed / config["outputs"]["socioeconomic_by_nuts3"]
    result.to_parquet(output, index=False)
    result.to_csv(output.with_suffix(".csv"), index=False)

    coverage_metrics = [
        "gdp_total_pps_million",
        "gdp_per_capita_pps",
        "gdp_density_pps_per_km2",
    ]
    coverage = result.groupby("year").agg(
        target_nuts3=("NUTS_ID", "size"),
        **{
            f"n_{metric}": (metric, "count")
            for metric in coverage_metrics
            if metric in result.columns
        },
    )

    print(f"GDP dataset: {socioeconomic['gdp_dataset']}")
    print(
        "Primary model GDP covariate: "
        f"{socioeconomic['model_gdp_metric']}"
    )
    print()
    print("Coverage by year:")
    print(coverage.to_string())
    print()
    print(
        "Note: GDP density is retained as a derived metric but should not be "
        "entered in the same model as population density + GDP per capita."
    )
    print(f"Wrote: {output}")
    print(f"Wrote: {output.with_suffix('.csv')}")


if __name__ == "__main__":
    main()
