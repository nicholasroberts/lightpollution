#!/usr/bin/env python3
"""Reconstruct Jan-1 population density on FIXED NUTS2016 NUTS3 land area.

Scientific sensitivity analysis; does NOT change production Eurostat density,
VIIRS radiance, spider/moth PRC or their calibrations. It compares published
density and population-on-1-January / constant land area on exactly the SAME
complete cohort. All unavailable or suspect geographical observations are
flagged rather than filled from the problematic density input.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from project_config import PROJECT_ROOT as ROOT
from validate_population_fallback import prepare_area, prepare_population
from analyze_country_population_prc_asymmetry import (
    countries_by_year, directional_summary
)
from visualize_annual_population_prc_changes import inputs, calculate

YEARS = tuple(range(2014, 2021))
TAXA = ("spider", "moth")
RAW = ROOT / "data/raw/eurostat"
DEFAULT_OUT = ROOT / "data/processed/nuts2016_jan1_population_density_sensitivity"


def load_json(path: Path) -> tuple[dict, dict]:
    data = path.read_bytes()
    payload = json.loads(data)
    if "dimension" not in payload or "id" not in payload:
        raise ValueError(f"Unexpected Eurostat JSON-stat structure: {path}")
    return payload, {
        "source": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(data).hexdigest(),
        "eurostat_updated": payload.get("updated"),
        "dataset_label": payload.get("label"),
    }


def fixed_land_areas(area: pd.DataFrame) -> pd.DataFrame:
    """Prefer 2019 reg_area3 land area (NUTS2016 era), hold it fixed for all years."""
    a = area.copy()
    a = a[np.isfinite(a.land_area_km2) & a.land_area_km2.gt(0)]
    if a.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Area dataset has ambiguous NUTS3/year values")
    # 2019 is an appropriate reference year for the frozen NUTS2016 cohort;
    # 2018/2020 are backups, not an unrecorded time-varying denominator.
    eligible = a[a.year.isin((2018, 2019, 2020))].copy()
    if eligible.empty:
        raise ValueError("No 2018-2020 land area observations in reg_area3")
    eligible["year_priority"] = eligible.year.map({2019: 0, 2018: 1, 2020: 2})
    selected = (eligible.sort_values(["NUTS_ID", "year_priority"])
                .drop_duplicates("NUTS_ID", keep="first"))
    selected = selected.rename(columns={
        "land_area_km2": "fixed_land_area_km2",
        "year": "land_area_reference_year",
        "status": "land_area_status",
    })[["NUTS_ID", "fixed_land_area_km2", "land_area_reference_year",
        "land_area_status"]]
    spread = eligible.groupby("NUTS_ID").land_area_km2.agg(["min", "max", "count"])
    spread["area_spread_fraction"] = (spread["max"] - spread["min"]) / spread["min"]
    spread = spread[["count", "area_spread_fraction"]].rename(
        columns={"count": "n_2018_2020_area_years"}).reset_index()
    selected = selected.merge(spread, on="NUTS_ID", validate="one_to_one")
    selected["area_change_over_1pct"] = selected.area_spread_fraction.gt(.01)
    return selected


def reconstruct(full: pd.DataFrame, population: pd.DataFrame,
                land: pd.DataFrame) -> pd.DataFrame:
    """Append derived density and transparent coverage flags to the old panel."""
    required = {"NUTS_ID", "year", "population_jan1", "status"}
    if not required.issubset(population):
        raise ValueError(f"Missing Jan-1 population columns: {required-set(population.columns)}")
    pop = population[population.year.isin(YEARS)].copy()
    if pop.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Jan-1 population ambiguous by NUTS3 and year")
    pop = pop.rename(columns={"status": "jan1_population_status"})
    panel = full.merge(
        pop[["NUTS_ID", "year", "population_jan1", "jan1_population_status"]],
        on=["NUTS_ID", "year"], how="left", validate="one_to_one")
    panel = panel.merge(land, on="NUTS_ID", how="left", validate="many_to_one")
    panel["original_eurostat_density"] = panel["population_density"]
    panel["density_jan1_fixed_land"] = (
        panel.population_jan1 / panel.fixed_land_area_km2)
    panel["jan1_available"] = (panel.population_jan1.gt(0)
                               & np.isfinite(panel.population_jan1))
    panel["land_area_available"] = (
        panel.fixed_land_area_km2.gt(0)
        & np.isfinite(panel.fixed_land_area_km2))
    panel["area_geography_warning"] = panel.area_change_over_1pct.fillna(False)
    panel["reconstruction_eligible_year"] = (
        panel.jan1_available & panel.land_area_available
        & ~panel.area_geography_warning
        & np.isfinite(panel.density_jan1_fixed_land))
    return panel


def same_complete_cohort(panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """All analyses use identical NUTS3 regions for ALL 2014..2020 years."""
    if panel.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate region/year in original source")
    if not panel.year.isin(YEARS).all():
        raise ValueError("Unexpected year in original complete-cohort input")
    n = panel.groupby("NUTS_ID").year.nunique()
    if not n.eq(len(YEARS)).all():
        raise ValueError("Baseline population/PRC cohort is not complete")
    eligible = panel.groupby("NUTS_ID").reconstruction_eligible_year.all()
    ids = eligible.index[eligible]
    matched = panel[panel.NUTS_ID.isin(ids)].copy()
    if matched.empty:
        raise ValueError("No complete reconstructed Jan-1 density cohort")
    if len(matched) != len(ids) * len(YEARS):
        raise ValueError("Matched cohort differs among endpoint years")
    return matched, pd.DataFrame({
        "NUTS_ID": eligible.index,
        "complete_reconstructed_2014_2020": eligible.values,
    })


def compare_countries(matched: pd.DataFrame) -> pd.DataFrame:
    """Existing country correlation/threshold functions, unchanged."""
    rows = []
    base = matched[matched.year.eq(2014)].set_index("NUTS_ID")
    for source, col in [
        ("published_demo_r_d3dens", "original_eurostat_density"),
        ("jan1_population_div_fixed_land", "density_jan1_fixed_land"),
    ]:
        for end_year in (2019, 2020):
            d = matched[matched.year.eq(end_year)].copy()
            d["cumulative_delta_population_density"] = (
                d[col] - d.NUTS_ID.map(base[col]))
            d["cumulative_pct_population_density"] = (
                100 * d["cumulative_delta_population_density"]
                / d.NUTS_ID.map(base[col]))
            for country, g in [("ALL", d), *list(d.groupby("country"))]:
                for taxon in TAXA:
                    s = directional_summary(
                        g, country, taxon, "cumulative", 5.0, .001)
                    s["density_source"] = source
                    s["start_year"] = 2014
                    s["end_year"] = end_year
                    rows.append(s)
    return pd.DataFrame(rows).sort_values(
        ["country", "taxon", "end_year", "density_source"]).reset_index(drop=True)


def coverage_table(panel: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    """Expose all excluded country/year records explicitly."""
    complete_ids = set(matched.NUTS_ID)
    rows = []
    for (country, year), g in panel.groupby(["country", "year"]):
        rows.append({
            "country": country, "year": int(year),
            "original_cohort_regions": g.NUTS_ID.nunique(),
            "n_jan1_available": int(g.jan1_available.sum()),
            "n_area_available": int(g.land_area_available.sum()),
            "n_area_geometry_warning": int(g.area_geography_warning.sum()),
            "n_complete_reconstructed_2014_2020": int(
                g.NUTS_ID.isin(complete_ids).sum()),
            "n_excluded_missing_jan1": int((~g.jan1_available).sum()),
            "n_excluded_missing_area": int((~g.land_area_available).sum()),
        })
    return pd.DataFrame(rows)


def paired_region_endpoints(matched: pd.DataFrame) -> pd.DataFrame:
    """Region-level audit table, each region once, with two density measures."""
    columns = ["NUTS_ID", "country", "year", "population_jan1",
               "fixed_land_area_km2", "land_area_reference_year",
               "original_eurostat_density", "density_jan1_fixed_land",
               "cumulative_delta_spider_prc", "cumulative_delta_moth_prc"]
    d = matched[columns].copy()
    index = ["NUTS_ID", "country", "fixed_land_area_km2",
             "land_area_reference_year"]
    names = ["population_jan1", "original_eurostat_density",
             "density_jan1_fixed_land", "cumulative_delta_spider_prc",
             "cumulative_delta_moth_prc"]
    p = d.pivot(index=index, columns="year", values=names)
    p.columns = [f"{field}_{year}" for field, year in p.columns]
    p = p.reset_index()
    for source in ("original_eurostat_density", "density_jan1_fixed_land"):
        for year in (2019, 2020):
            p[f"change_{source}_2014_{year}"] = (
                p[f"{source}_{year}"] - p[f"{source}_2014"])
        p[f"change_{source}_2019_2020"] = (
            p[f"{source}_2020"] - p[f"{source}_2019"])
    return p.sort_values(["country", "NUTS_ID"]).reset_index(drop=True)


def synthetic_test() -> None:
    base = pd.DataFrame([
        {"NUTS_ID": nid, "year": year, "country": "DE",
         "population_density": 10., "spider_prc": .2, "moth_prc": .4}
        for nid in ("DE001", "DE002") for year in YEARS])
    p = pd.DataFrame([
        {"NUTS_ID": nid, "year": year,
         "population_jan1": (100 if nid == "DE001" else 200) + year - 2014,
         "status": None}
        for nid in ("DE001", "DE002") for year in YEARS])
    a = pd.DataFrame([
        {"NUTS_ID": nid, "year": year, "land_area_km2": 10., "status": None}
        for nid in ("DE001", "DE002") for year in (2018, 2019, 2020)])
    area = fixed_land_areas(a)
    assert area.land_area_reference_year.eq(2019).all()
    joined = reconstruct(base, p, area)
    good, _ = same_complete_cohort(joined)
    assert good.NUTS_ID.nunique() == 2
    assert np.isclose(good.loc[
        (good.NUTS_ID.eq("DE001")) & good.year.eq(2020),
        "density_jan1_fixed_land"].iloc[0], 10.6)
    assert np.isclose(good.loc[
        (good.NUTS_ID.eq("DE001")) & good.year.eq(2014),
        "density_jan1_fixed_land"].iloc[0], 10.)
    print("Jan-1 fixed-land reconstruction/paired-cohort self-test: PASS")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    synthetic_test()
    if args.self_test:
        return
    output = args.out if args.out.is_absolute() else ROOT / args.out

    pjan, pm = load_json(RAW / "demo_r_pjanaggr3_nuts3.json")
    area_json, am = load_json(RAW / "reg_area3_nuts3.json")
    pjan_frame = prepare_population(pjan)
    area_frame = prepare_area(area_json)
    land = fixed_land_areas(area_frame)
    original, annual, cumulative = calculate(inputs())
    _, cm, _ = countries_by_year(annual, cumulative)
    orig_regions = set(original.NUTS_ID)
    if len(orig_regions) != 1317:
        raise ValueError("Expected 1,317-region original country-analysis cohort")
    meta = cm[["NUTS_ID", "country"]].drop_duplicates()
    original = original.merge(meta, on="NUTS_ID", validate="many_to_one")
    panel = reconstruct(original, pjan_frame, land)
    matched, flags = same_complete_cohort(panel)
    summary = compare_countries(matched)
    region_endpoints = paired_region_endpoints(matched)
    coverage = coverage_table(panel, matched)
    original_count = len(orig_regions)
    matched_count = matched.NUTS_ID.nunique()

    output.mkdir(parents=True, exist_ok=True)
    panel[["NUTS_ID", "country", "year", "original_eurostat_density",
           "population_jan1", "jan1_population_status",
           "fixed_land_area_km2", "land_area_reference_year",
           "land_area_status", "area_spread_fraction",
           "area_geography_warning", "density_jan1_fixed_land",
           "reconstruction_eligible_year"]].to_csv(
        output / "regional_annual_density_2014_2020_full_audit.csv", index=False)
    region_endpoints.to_csv(
        output / "matched_region_density_2014_2019_2020.csv", index=False)
    summary.to_csv(
        output / "country_prc_correlations_published_vs_jan1.csv", index=False)
    coverage.to_csv(
        output / "coverage_and_exclusions_by_country_year.csv", index=False)
    flags.to_csv(
        output / "reconstructed_complete_cohort_flags.csv", index=False)
    land.to_csv(
        output / "fixed_land_area_source_2019_nuts2016.csv", index=False)
    manifest = {
        "population_source": pm, "land_area_source": am,
        "population_definition": "January-1 population / FIXED NUTS3 reg_area3 land km2",
        "land_area_reference": "2019 preferred (NUTS2016 era); 2018/2020 fallback labelled",
        "area_quality": "Exclude regions with >1 percent range in 2018-2020 land-area series; report exclusions",
        "original_complete_cohort_regions": original_count,
        "reconstructed_common_complete_cohort_regions": matched_count,
        "comparison_years": [2014, 2019, 2020],
        "PRC_source": "Original locked spider/moth snow-corrected annual panel; unchanged",
        "warnings": [
            "January-1 density is not Eurostat annual-average density",
            "NUTS2016 code continuity does not independently certify unchanged polygon borders",
            "Missing January counts or matched land areas are not filled",
            "Same reconstructed NUTS3 cohort is used for both old/new correlations and both endpoints",
            "Reg_area3 land area, not EPSG4326 polygon area or total area including water",
            "No causal interpretation of population-PRC association",
        ],
        "production_inputs_modified": False,
    }
    (output / "reconstruction_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"\nCommon complete cohort: {matched_count} / {original_count} NUTS2016 NUTS3 regions")
    print("\nCountries of interest; comparable Spearman rank correlations:")
    focus = summary[(summary.country.isin(
        ["ALL", "TR", "DE", "FR", "IE", "IT", "BE", "NL", "AL"]))
                    & summary.taxon.eq("spider")]
    print(focus[["country", "density_source", "end_year", "n_regions",
                 "spearman_rho_pop_vs_delta_prc"]].to_string(index=False))
    print(f"\nOutputs: {output}")


if __name__ == "__main__":
    main()
