#!/usr/bin/env python3
"""Diagnose *exactly* where UK NUTS2016 NUTS3 left the 2014–2020 cohort.

READ ONLY: no downloads, reprocessing, imputation, changes to frozen
1,273-region demographic data, or modifications to Figure 2/4.
Uses original local Eurostat JSON-Stat snapshots, already processed
density/PRC, immutable NUTS2016 geometry and existing cohort audits.
Run on EcoVis after `git pull`, not on absent remote datasets.

Report separately:
  geometry -> raw DEMO_R_D3DENS -> processed density -> old 1317 cohort
  geometry -> raw DEMO_R_PJANAGGR3 -> frozen reg_area3 denominator
  geometry -> locked native-pixel spider/moth PRC -> final 1273 cohort.
Do NOT equate 'code listed in a Eurostat dimension' with a valid
year-specific numeric observation. Never mix Jan1 with ONS mid-year.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from project_config import PROJECT_ROOT as ROOT, load_config
from historical_prc_pixel_nuts3 import load_nuts
from fetch_eurostat_population import jsonstat_to_frame as density_jsonstat
from validate_population_fallback import prepare_population, prepare_area
from reconstruct_population_density_jan1_nuts2016 import fixed_land_areas


YEARS = tuple(range(2014, 2021))
RAW = ROOT / "data/raw/eurostat"
PROC = ROOT / "data/processed"
SOURCES = {
    "density_raw": RAW / "demo_r_d3dens_nuts3.json",
    "population_jan1_raw": RAW / "demo_r_pjanaggr3_nuts3.json",
    "area_raw": RAW / "reg_area3_nuts3.json",
    "density_processed": PROC / "nuts3_population_density.parquet",
    "original_cohort": (
        PROC / "nuts2016_jan1_population_density_sensitivity"
        / "regional_annual_density_2014_2020_full_audit.csv"
    ),
    "final_cohort": (
        PROC / "jan1_population_radiance_prc_country_nl333"
        / "matched_region_change_2014_2020.csv"
    ),
    "spider": (
        PROC / "historical_prc_2014_2024_spider_snow_corrected"
        / "prc_nuts3_annual_cumulative.csv"
    ),
    "moth": (
        PROC / "historical_prc_2014_2024_moth_snow_corrected"
        / "prc_nuts3_annual_cumulative.csv"
    ),
}
DEFAULT_OUT = PROC / "publication_figures/uk_population_nuts2016_audit"


def load_json(path: Path):
    with path.open("r", encoding="utf8") as f:
        return json.load(f)


def geo_codes(payload: dict) -> set[str]:
    index = payload["dimension"]["geo"]["category"]["index"]
    return set(map(str, index if isinstance(index, dict) else index))


def validated_unique(d: pd.DataFrame, key: list[str], name: str) -> pd.DataFrame:
    if d.duplicated(key).any():
        sample = d[d.duplicated(key, keep=False)].head(5)
        raise ValueError(f"{name}: duplicate {key}: {sample.to_dict('records')}")
    return d


def annual_raw_density(payload: dict):
    frame = density_jsonstat(payload)
    needed = {"geo", "time", "population_density"}
    if not needed.issubset(frame):
        raise ValueError(f"Eurostat raw density JSON lacks {sorted(needed - set(frame))}")
    if "freq" in frame:
        frame = frame.loc[frame.freq.eq("A")].copy()
    if "unit" in frame:
        frame = frame.loc[frame.unit.eq("PER_KM2")].copy()
    frame["year"] = pd.to_numeric(frame.time, errors="coerce")
    frame = frame.loc[
        frame.year.isin(YEARS) & frame.geo.astype(str).str.startswith("UK")
    ].copy()
    frame["year"] = frame.year.astype(int)
    frame["population_density"] = pd.to_numeric(
        frame.population_density, errors="coerce"
    )
    frame = frame.rename(columns={
        "geo": "NUTS_ID", "population_density": "raw_density"
    })
    return validated_unique(frame[["NUTS_ID", "year", "raw_density", "status"]],
                            ["NUTS_ID", "year"], "Raw DEMO_R_D3DENS")


def annual_processed_density(path: Path):
    df = pd.read_parquet(path)
    required = {"NUTS_ID", "year", "population_density"}
    if not required.issubset(df):
        raise ValueError(f"Processed density lacks {sorted(required-set(df))}")
    d = df.loc[
        df.year.isin(YEARS) & df.NUTS_ID.astype(str).str.startswith("UK")
    ].copy()
    keep = ["NUTS_ID", "year", "population_density"]
    for name in ("nuts_release", "matches_nuts_release", "status"):
        if name in d:
            keep.append(name)
    d = d[keep].rename(columns={
        "population_density": "processed_density",
        "nuts_release": "processed_nuts_release",
        "matches_nuts_release": "matches_nuts_release",
        "status": "processed_density_status",
    })
    d.NUTS_ID = d.NUTS_ID.astype(str)
    return validated_unique(d, ["NUTS_ID", "year"], "Processed DEMO_R_D3DENS")


def annual_jan1(payload: dict):
    p = prepare_population(payload)
    p["NUTS_ID"] = p.NUTS_ID.astype(str)
    p = p.loc[p.year.isin(YEARS) & p.NUTS_ID.str.startswith("UK")].copy()
    p["population_jan1"] = pd.to_numeric(p.population_jan1, errors="coerce")
    return validated_unique(
        p.rename(columns={"status": "jan1_source_status"}),
        ["NUTS_ID", "year"], "Raw January population"
    )


def annual_taxon(path: Path, name: str):
    a = pd.read_csv(path, dtype={"NUTS_ID": str},
                    usecols=["NUTS_ID", "year", "prc_mean", "snow_status"])
    a = a.loc[a.year.isin(YEARS) & a.NUTS_ID.str.startswith("UK")].copy()
    if not a.snow_status.eq("locked").all():
        raise ValueError(f"{name} input includes a non-locked source")
    return validated_unique(
        a[["NUTS_ID", "year", "prc_mean"]].rename(
            columns={"prc_mean": f"{name}_prc"}
        ), ["NUTS_ID", "year"], f"{name} locked PRC"
    )


def annual_cohort(path: Path):
    a = pd.read_csv(path, dtype={"NUTS_ID": str},
                    usecols=["NUTS_ID", "year"])
    a = a.loc[a.year.isin(YEARS) & a.NUTS_ID.str.startswith("UK")].copy()
    a["in_original_1317"] = True
    return validated_unique(a, ["NUTS_ID", "year"],
                            "Original 1317-region annual reconstruction input")


def annual_skeleton(ids: set[str]):
    return pd.MultiIndex.from_product(
        [sorted(ids), YEARS], names=["NUTS_ID", "year"]
    ).to_frame(index=False)


def classify_reason(g: pd.DataFrame) -> str:
    if not bool(g.raw_density_valid.all()):
        if not bool(g.processed_density_valid.all()):
            return "Original Eurostat population density absent/non-numeric in >=1 year"
        return "Raw snapshot missing but processed population density exists: audit vintage/source"
    if not bool(g.processed_density_valid.all()):
        return "Raw density exists but processed density is missing: pipeline/filter issue"
    if not bool(g.processed_nuts2016_match.all()):
        return "Processed density not matched to locked NUTS2016 code/year"
    if not bool(g.in_original_1317.all()):
        return "Population exists but dropped during original complete-cohort/PRC join"
    if not bool(g.jan1_valid.all()):
        return "Original cohort present but Jan1 population missing/nonpositive"
    if not bool(g.land_valid.all()):
        return "Original cohort present but 2018-2020 fixed land area missing"
    if not bool(g.spider_valid.all()) or not bool(g.moth_valid.all()):
        return "Original cohort present but locked PRC absent/nonfinite"
    if not bool(g.in_final_1273.all()):
        return "All prerequisites present but absent from final 1273 source: investigate join"
    return "Included in frozen 1273 cohort"


def audit(config: dict, paths: dict[str, Path]):
    nuts = load_nuts(config)
    if len(nuts) != 1522:
        raise ValueError("Fixed NUTS2016 geographic cohort is not 1522")
    uk = nuts.loc[
        nuts.CNTR_CODE.eq("UK"),
        ["NUTS_ID", "NUTS_NAME", "CNTR_CODE"]
    ].copy()
    if len(uk) == 0 or not uk.NUTS_ID.is_unique:
        raise ValueError("No distinct UK regions in locked NUTS2016 geometry")
    ids = set(uk.NUTS_ID.astype(str))
    raw_density = load_json(paths["density_raw"])
    raw_jan1 = load_json(paths["population_jan1_raw"])
    raw_area = load_json(paths["area_raw"])
    jsonstat_membership = {
        "density_geo_category_has_UK": len({
            c for c in geo_codes(raw_density) if c.startswith("UK")
        }),
        "jan1_geo_category_has_UK": len({
            c for c in geo_codes(raw_jan1) if c.startswith("UK")
        }),
        "area_geo_category_has_UK": len({
            c for c in geo_codes(raw_area) if c.startswith("UK")
        }),
        "density_uk_codes_outside_frozen_nuts2016": sorted(
            c for c in geo_codes(raw_density) if c.startswith("UK")
            and len(c)==5 and c not in ids
        ),
        "jan1_uk_codes_outside_frozen_nuts2016": sorted(
            c for c in geo_codes(raw_jan1) if c.startswith("UK")
            and len(c)==5 and c not in ids
        ),
    }
    raw_d = annual_raw_density(raw_density)
    processed = annual_processed_density(paths["density_processed"])
    pjan = annual_jan1(raw_jan1)
    area = prepare_area(raw_area)
    land = fixed_land_areas(area)
    land = land.loc[land.NUTS_ID.isin(ids),
                    ["NUTS_ID", "fixed_land_area_km2",
                     "land_area_reference_year", "area_geography_warning"]
                   if "area_geography_warning" in land else
                   ["NUTS_ID", "fixed_land_area_km2",
                    "land_area_reference_year", "area_change_over_1pct"]
                   ].copy()
    if "area_change_over_1pct" in land:
        land = land.rename(columns={
            "area_change_over_1pct": "area_geography_warning"
        })
    original = annual_cohort(paths["original_cohort"])
    final = pd.read_csv(paths["final_cohort"], dtype={"NUTS_ID": str},
                        usecols=["NUTS_ID"])
    if final.NUTS_ID.duplicated().any():
        raise ValueError("Frozen 1273 source has duplicate IDs")
    frozen_uk = set(final.loc[final.NUTS_ID.str.startswith("UK"), "NUTS_ID"])
    d = annual_skeleton(ids)
    for extra in (
        raw_d, processed, pjan, annual_taxon(paths["spider"], "spider"),
        annual_taxon(paths["moth"], "moth"), original,
    ):
        d = d.merge(extra, on=["NUTS_ID", "year"],
                    how="left", validate="one_to_one")
    d = d.merge(land, on="NUTS_ID", how="left", validate="many_to_one")
    d = d.merge(uk[["NUTS_ID", "NUTS_NAME"]], on="NUTS_ID",
                how="left", validate="many_to_one")
    d["in_final_1273"] = d.NUTS_ID.isin(frozen_uk)
    d["raw_density_valid"] = pd.to_numeric(
        d.raw_density, errors="coerce"
    ).gt(0) & np.isfinite(pd.to_numeric(d.raw_density, errors="coerce"))
    d["processed_density_valid"] = pd.to_numeric(
        d.processed_density, errors="coerce"
    ).gt(0) & np.isfinite(pd.to_numeric(d.processed_density, errors="coerce"))
    d["processed_nuts2016_match"] = (
        d.get("processed_nuts_release", pd.Series(np.nan,index=d.index)).eq(2016)
        & d.get("matches_nuts_release", pd.Series(False,index=d.index))
         .fillna(False).astype(bool)
    )
    d["jan1_valid"] = pd.to_numeric(
        d.population_jan1, errors="coerce"
    ).gt(0) & np.isfinite(pd.to_numeric(d.population_jan1, errors="coerce"))
    d["land_valid"] = (
        pd.to_numeric(d.fixed_land_area_km2, errors="coerce").gt(0)
        & ~d.area_geography_warning.fillna(False).astype(bool)
    )
    for taxon in ("spider", "moth"):
        d[f"{taxon}_valid"] = pd.to_numeric(
            d[f"{taxon}_prc"], errors="coerce"
        ).between(0., 1., inclusive="both")
    d["in_original_1317"] = d.in_original_1317.fillna(False).astype(bool)
    if len(d) != len(ids)*len(YEARS):
        raise ValueError("Audit join changed the frozen UK spatial/year skeleton")

    indicators = [
        "raw_density_valid", "processed_density_valid",
        "processed_nuts2016_match", "jan1_valid", "land_valid",
        "spider_valid", "moth_valid", "in_original_1317", "in_final_1273"
    ]
    per_year = d.groupby("year")[indicators].sum().astype(int).reset_index()
    per_year.insert(1, "nuts2016_uk_regions", len(ids))
    reasons = pd.DataFrame([
        {"NUTS_ID": name, "first_exclusion_reason": classify_reason(g),
         **{f"complete_{field}": bool(g[field].all()) for field in indicators}}
        for name, g in d.groupby("NUTS_ID")
    ]).sort_values(["first_exclusion_reason","NUTS_ID"]).reset_index(drop=True)
    summary = {
        "years": list(YEARS),
        "locked_nuts2016_uk_regions": len(ids),
        "uk_2020_region_known_example": "UKK11",
        "UKK11_in_locked_geometry": "UKK11" in ids,
        "frozen_1273_uk_regions": len(frozen_uk),
        "uk_regions_complete_raw_density": int(reasons.complete_raw_density_valid.sum()),
        "uk_regions_complete_processed_density": int(reasons.complete_processed_density_valid.sum()),
        "uk_regions_complete_Jan1_population": int(reasons.complete_jan1_valid.sum()),
        "uk_regions_with_fixed_2018_20_land_area": int(reasons.complete_land_valid.sum()),
        "uk_regions_complete_spider": int(reasons.complete_spider_valid.sum()),
        "uk_regions_complete_moth": int(reasons.complete_moth_valid.sum()),
        "uk_regions_potentially_eligible_Jan1_land_PRC_independent_of_original": int(
            (reasons.complete_jan1_valid &
             reasons.complete_land_valid &
             reasons.complete_spider_valid &
             reasons.complete_moth_valid).sum()
        ),
        "first_exclusion_reason_counts": reasons.first_exclusion_reason.value_counts().to_dict(),
        "eurostat_jsonstat_geo_dimension": jsonstat_membership,
        "input_files": {name:str(p) for name,p in paths.items()},
        "IMPORTANT": (
            "A present geo code is not evidence of available annual counts; "
            "the 2014-2020 numerical value must be checked. "
            "Jan1 and ONS midyear population are NOT interchangeable."
        ),
    }
    return d, per_year, reasons, summary


def self_test():
    s = pd.DataFrame({
        "raw_density_valid": [False]*7,
        "processed_density_valid": [False]*7,
        "processed_nuts2016_match": [False]*7,
        "jan1_valid": [True]*7,
        "land_valid": [True]*7,
        "spider_valid": [True]*7,
        "moth_valid": [True]*7,
        "in_original_1317": [False]*7,
        "in_final_1273": [False]*7,
    })
    assert "absent/non-numeric" in classify_reason(s)
    s["raw_density_valid"] = True
    assert "pipeline/filter" in classify_reason(s)
    s["processed_density_valid"] = True
    s["processed_nuts2016_match"] = True
    assert "original complete-cohort" in classify_reason(s)
    s["in_original_1317"] = True
    s["in_final_1273"] = True
    assert classify_reason(s) == "Included in frozen 1273 cohort"
    print("UK cohort audit logical stage tests: PASS")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    self_test()
    if args.self_test:
        return
    missing = [f"{label}: {p}" for label,p in SOURCES.items() if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "Input sources not found on EcoVis; no replacement downloaded:\n"+
            "\n".join(missing)
        )
    d, years, reasons, summary = audit(load_config(), SOURCES)
    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out/"uk_nuts3_year_by_year_source_trace.csv", index=False)
    years.to_csv(out/"uk_per_year_stage_counts.csv", index=False)
    reasons.to_csv(out/"uk_region_first_exclusion_reason.csv", index=False)
    (out/"uk_population_cohort_audit_summary.json").write_text(
        json.dumps(summary, indent=2)+"\n", encoding="utf-8"
    )
    print("\nUK NUTS2016 stage-by-stage, each year:", flush=True)
    print(years.to_string(index=False), flush=True)
    print("\nUK region-level first exclusion reasons:", flush=True)
    print(reasons.first_exclusion_reason.value_counts().to_string(), flush=True)
    print("\nUKK11 (Bristol) case:", flush=True)
    print(d.loc[d.NUTS_ID.eq("UKK11"),[
        "NUTS_ID","year","raw_density","processed_density",
        "population_jan1","fixed_land_area_km2",
        "in_original_1317","in_final_1273"
    ]].to_string(index=False),flush=True)
    print("\nPotential UK regions eligible for Jan1+land+both locked PRCs "
          "WITHOUT the old density-stage dependency:",
          summary["uk_regions_potentially_eligible_Jan1_land_PRC_independent_of_original"],
          flush=True)
    print("\nAUDIT OUTPUTS",flush=True)
    for file in sorted(out.glob("uk_*")):
        print(f"  {file} ({file.stat().st_size/1024:.1f} KiB)",flush=True)


if __name__ == "__main__":
    main()
