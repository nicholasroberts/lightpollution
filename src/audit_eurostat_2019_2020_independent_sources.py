#!/usr/bin/env python3
"""Audit Eurostat 2019/2020 NUTS3 population density using independent source series.

READ-ONLY wrt original raw and processed project data. Optional --fetch-current
saves NEW API snapshots in a separate audit directory; never overwrites input
or the locked PRC products. Diagnostic only: no demographic imputation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd

from fetch_eurostat_population import jsonstat_to_frame
from validate_population_fallback import prepare_area, prepare_population
from project_config import PROJECT_ROOT as ROOT

RAW = ROOT / "data/raw/eurostat"
PREVIOUS = ROOT / "data/processed/prc_population_endpoint_qa_2014_2019_2020/region_paired_2014_2019_2020.csv"
OUT = ROOT / "data/processed/eurostat_2019_2020_independent_audit"
DATASETS = ("demo_r_d3dens", "demo_r_pjanaggr3", "reg_area3")
SELECTORS = {
    "demo_r_d3dens": {"freq": "A", "unit": "PER_KM2"},
    "demo_r_pjanaggr3": {"freq": "A", "unit": "NR", "sex": "T", "age": "TOTAL"},
    "reg_area3": {"freq": "A", "unit": "KM2"},
}
YEARS = (2019, 2020)


def input_json(path: Path) -> tuple[dict, dict]:
    data = path.read_bytes()
    parsed = json.loads(data)
    if "dimension" not in parsed or "id" not in parsed:
        raise ValueError(f"Not a Eurostat JSON-stat response: {path}")
    return parsed, {
        "path": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(data).hexdigest(),
        "updated": parsed.get("updated"),
        "label": parsed.get("label"),
    }


def fetch_snapshot(out: Path) -> dict[str, Path]:
    """Fetch current Eurostat JSON-stat selections without touching historic raw files."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    folder = out / "raw_snapshot" / stamp
    folder.mkdir(parents=True, exist_ok=False)
    result = {}
    for name in DATASETS:
        params = dict(lang="EN", geoLevel="nuts3",
                      sinceTimePeriod="2018", untilTimePeriod="2021",
                      **SELECTORS[name])
        url = ("https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/"
               + name + "?" + urlencode(params))
        print(f"Fetching a NEW audit snapshot: {name}", flush=True)
        with urlopen(url, timeout=180) as response:
            raw = response.read()
        payload = json.loads(raw)
        if "dimension" not in payload or "id" not in payload:
            raise ValueError(f"Unexpected Eurostat response for {name}: {str(payload)[:250]}")
        dest = folder / f"{name}_nuts3.json"
        dest.write_bytes(raw)
        result[name] = dest
    return result


def official_density(payload: dict) -> pd.DataFrame:
    p = jsonstat_to_frame(payload)
    for col, value in SELECTORS["demo_r_d3dens"].items():
        if col in p.columns:
            p = p[p[col].eq(value)]
    p["year"] = pd.to_numeric(p["time"], errors="coerce")
    p["official_density"] = pd.to_numeric(p["population_density"], errors="coerce")
    p = p[p.year.isin(YEARS)].rename(columns={"geo": "NUTS_ID"})
    q = p[["NUTS_ID", "year", "official_density", "status"]].copy()
    q["year"] = q.year.astype(int)
    if q.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("More than one official density measurement per NUTS3/year")
    return q


def by_year(df: pd.DataFrame, variable: str, prefix: str) -> pd.DataFrame:
    d = df[df.year.isin(YEARS)].copy()
    if d.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError(f"Duplicate NUTS3/year in {variable}")
    w = d.pivot(index="NUTS_ID", columns="year", values=variable)
    return w.rename(columns={year: f"{prefix}_{year}" for year in YEARS}).reset_index()


def chosen_land_area(areas: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Use most recent land area no later than year; label its reference year."""
    rows = []
    ids = set(cohort.NUTS_ID)
    areas = areas[areas.NUTS_ID.isin(ids) & areas.land_area_km2.gt(0)]
    for year in YEARS:
        d = areas[areas.year.le(year)].sort_values("year").drop_duplicates(
            "NUTS_ID", keep="last").copy()
        d = d[["NUTS_ID", "year", "land_area_km2"]]
        d["target_year"] = year
        rows.append(d.rename(columns={"year": "area_reference_year"}))
    return pd.concat(rows, ignore_index=True)


def build_side_by_side(paths: dict[str, Path], cohort: pd.DataFrame, label: str
                       ) -> tuple[pd.DataFrame, dict]:
    provenance = {}
    data = {}
    for name in DATASETS:
        data[name], provenance[name] = input_json(paths[name])
    official = official_density(data["demo_r_d3dens"])
    jan1 = prepare_population(data["demo_r_pjanaggr3"])
    jan1 = jan1[jan1.year.isin(YEARS)].copy()
    land = prepare_area(data["reg_area3"])
    land = land[land.NUTS_ID.isin(cohort.NUTS_ID)].copy()

    a = cohort.copy()
    for name, d, col in (
        ("published", official, "official_density"),
        ("jan1", jan1, "population_jan1"),
    ):
        a = a.merge(by_year(d, col, name), on="NUTS_ID", how="left",
                    validate="one_to_one")
        for year in YEARS:
            s = d[d.year.eq(year)][["NUTS_ID", "status"]].rename(
                columns={"status": f"{name}_status_{year}"})
            a = a.merge(s, on="NUTS_ID", how="left", validate="one_to_one")

    sel = chosen_land_area(land, cohort)
    for year in YEARS:
        s = sel.loc[sel.target_year.eq(year),
                    ["NUTS_ID", "land_area_km2", "area_reference_year"]].rename(
            columns={"land_area_km2": f"area_km2_{year}",
                     "area_reference_year": f"area_refyear_{year}"})
        a = a.merge(s, on="NUTS_ID", how="left", validate="one_to_one")
        a[f"jan1_per_km2_{year}"] = a[f"jan1_{year}"] / a[f"area_km2_{year}"]

    a["density_2020_eq_2019_processed"] = a.density_2019.eq(a.density_2020)
    a["density_2020_eq_2019_raw"] = (
        a.published_2019.notna() & a.published_2020.notna()
        & a.published_2019.eq(a.published_2020))
    a["jan1_2020_eq_2019"] = (
        a.jan1_2019.notna() & a.jan1_2020.notna()
        & a.jan1_2019.eq(a.jan1_2020))
    a["delta_published"] = a.published_2020 - a.published_2019
    a["delta_processed"] = a.density_2020 - a.density_2019
    a["delta_jan1_people"] = a.jan1_2020 - a.jan1_2019
    a["jan1_change_pct"] = (100 * a.delta_jan1_people
                             / a.jan1_2019.where(a.jan1_2019.gt(0)))
    # Avoid denominator movements masquerading as a population change.
    a["jan1_change_density_on_fixed_2019_land_area"] = (
        a.delta_jan1_people / a.area_km2_2019)
    a["raw_vs_processed_2019"] = a.published_2019 - a.density_2019
    a["raw_vs_processed_2020"] = a.published_2020 - a.density_2020
    a["source_label"] = label
    return a, provenance


def summarise(a: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for code, g in [("ALL", a), *list(a.groupby("country"))]:
        mask = g.published_2019.notna() & g.published_2020.notna()
        jmask = g.jan1_2019.notna() & g.jan1_2020.notna()
        both = mask & jmask
        copied_but_changed = (both & g.density_2020_eq_2019_raw
                              & ~g.jan1_2020_eq_2019)
        rows.append({
            "source_label": a.source_label.iloc[0],
            "country": code,
            "n_regions": len(g),
            "n_official_pairs": int(mask.sum()),
            "n_official_2019_2020_identical": int(g.loc[mask, "density_2020_eq_2019_raw"].sum()),
            "n_processed_2019_2020_identical": int(g.density_2020_eq_2019_processed.sum()),
            "n_jan1_pairs": int(jmask.sum()),
            "n_jan1_values_different": int((jmask & ~g.jan1_2020_eq_2019).sum()),
            "n_official_identical_but_jan1_changed": int(copied_but_changed.sum()),
            "median_jan1_change_pct": float(g.jan1_change_pct.median()) if jmask.any() else None,
            "n_raw_processed_2019_disagree": int(
                (g.raw_vs_processed_2019.abs().gt(1e-8)).sum()),
            "n_raw_processed_2020_disagree": int(
                (g.raw_vs_processed_2020.abs().gt(1e-8)).sum()),
            "n_same_land_area_reference_year": int(
                (g.area_refyear_2019.eq(g.area_refyear_2020)
                 & g.area_refyear_2019.notna()).sum()),
            "n_exact_area_2020_available": int(g.area_refyear_2020.eq(2020).sum()),
        })
    return pd.DataFrame(rows)


def self_test() -> None:
    demo = pd.DataFrame({
        "NUTS_ID": ["DE001", "DE002"], "country": ["DE", "DE"],
        "density_2019": [10., 20.], "density_2020": [10., 20.]})
    a = demo.copy()
    a["published_2019"] = a.density_2019
    a["published_2020"] = a.density_2020
    a["jan1_2019"] = [100., 200.]
    a["jan1_2020"] = [101., 202.]
    a["area_refyear_2019"] = [2018, 2018]
    a["area_refyear_2020"] = [2018, 2018]
    a["density_2020_eq_2019_raw"] = True
    a["density_2020_eq_2019_processed"] = True
    a["jan1_2020_eq_2019"] = False
    a["jan1_change_pct"] = [1., 1.]
    a["raw_vs_processed_2019"] = [0., 0.]
    a["raw_vs_processed_2020"] = [0., 0.]
    a["source_label"] = "synthetic"
    q = summarise(a)
    assert q.loc[q.country.eq("DE"), "n_official_identical_but_jan1_changed"].iloc[0] == 2
    assert q.loc[q.country.eq("ALL"), "n_raw_processed_2020_disagree"].iloc[0] == 0
    print("Eurostat independent-source comparison self-test: PASS")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--fetch-current", action="store_true",
                   help="Download new dated Eurostat JSON snapshots; leave historic raw data untouched")
    p.add_argument("--out", type=Path, default=OUT)
    args = p.parse_args()
    self_test()
    if args.self_test:
        return
    out = args.out if args.out.is_absolute() else ROOT / args.out
    if not PREVIOUS.exists():
        raise FileNotFoundError(f"Run the previous 2019/2020 endpoint QA first: {PREVIOUS}")
    pair = pd.read_csv(PREVIOUS, dtype={"NUTS_ID": str, "country": str})
    need = {"NUTS_ID", "country", "density_2019", "density_2020"}
    if not need.issubset(pair):
        raise ValueError(f"Expected earlier matched panel columns: {need - set(pair.columns)}")
    if len(pair) != 1317 or pair.NUTS_ID.duplicated().any():
        raise ValueError("Reference matched 1,317-region cohort changed")
    cohort = pair[["NUTS_ID", "country", "density_2019", "density_2020"]].copy()
    out.mkdir(parents=True, exist_ok=True)
    snapshots = []
    historic = {name: RAW / f"{name}_nuts3.json" for name in DATASETS}
    if all(path.exists() for path in historic.values()):
        snapshots.append(("archived_local", historic))
    else:
        missing = [name for name, path in historic.items() if not path.exists()]
        print(f"Original raw source incomplete; cannot compare archive: {missing}")
    if args.fetch_current:
        snapshots.append(("fresh_eurostat", fetch_snapshot(out)))
    if not snapshots:
        raise FileNotFoundError(
            "No complete archived raw source. Re-run with --fetch-current to create audit snapshots.")
    provenance = {}
    all_summary = []
    for label, paths in snapshots:
        print(f"\nAnalysing source: {label}", flush=True)
        detail, provenance[label] = build_side_by_side(paths, cohort, label)
        summary = summarise(detail)
        detail.to_csv(out / f"{label}_regional_2019_2020_comparison.csv", index=False)
        summary.to_csv(out / f"{label}_country_2019_2020_comparison.csv", index=False)
        all_summary.append(summary)
        print(summary.loc[summary.country.isin(
            ["ALL", "DE", "TR", "IE", "FR", "IT", "NL", "AL"]),
            ["country", "n_regions", "n_official_pairs",
             "n_official_2019_2020_identical", "n_jan1_pairs",
             "n_official_identical_but_jan1_changed",
             "n_raw_processed_2020_disagree"]].to_string(index=False))
    pd.concat(all_summary, ignore_index=True).to_csv(
        out / "all_sources_country_comparison.csv", index=False)
    manifest = {
        "population_reference": "locked 1317-region complete 2014-2020 NUTS2016 cohort",
        "source_provenance": provenance,
        "statistical_caution": (
            "Eurostat demo_r_d3dens uses annual-average population/land area, "
            "whereas demo_r_pjanaggr3 is January 1 population. "
            "Jan1/area is not an exact reconstruction of the published indicator."),
        "method": "Same geographic NUTS3 codes; no inferred boundary conversion",
        "production_inputs_changed": False,
    }
    (out / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nAll QA outputs in: {out}")


if __name__ == "__main__":
    main()
