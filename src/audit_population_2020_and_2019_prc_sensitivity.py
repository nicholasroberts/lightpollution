#!/usr/bin/env python3
"""Fixed-cohort 2014-2019 vs 2014-2020 country PRC/population sensitivity.

Runs only on the existing, approved NUTS2016 2014-2020 seven-year panel
and independently locked spider/moth PRC fields. Does not refit PRC,
modify input data, or claim population causes lighting change.

Outputs under data/processed/prc_population_endpoint_qa_2014_2019_2020.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_country_population_prc_asymmetry import (
    countries_by_year, directional_summary, exact_corr,
)
from project_config import PROJECT_ROOT as ROOT
from visualize_annual_population_prc_changes import inputs, calculate

DEFAULT_OUT = ROOT / "data/processed/prc_population_endpoint_qa_2014_2019_2020"
TAXA = ("spider", "moth")
YEARS = (2014, 2019, 2020)


def region_endpoints(full: pd.DataFrame) -> pd.DataFrame:
    """Construct paired regions without admitting a different 2019 cohort."""
    if full.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate NUTS3/year observations")
    counts = full.groupby("NUTS_ID").year.nunique()
    if not counts.eq(7).all() or not set(full.year.unique()) == set(range(2014, 2021)):
        raise ValueError("Expected complete fixed seven-year panel (2014-2020)")

    out = full[["NUTS_ID", "country"]].drop_duplicates()
    for year in YEARS:
        d = full.loc[full.year.eq(year),
            ["NUTS_ID", "population_density", "spider_prc", "moth_prc"]]
        d = d.rename(columns={
            "population_density": f"density_{year}",
            "spider_prc": f"spider_prc_{year}",
            "moth_prc": f"moth_prc_{year}",
        })
        out = out.merge(d, on="NUTS_ID", validate="one_to_one")
    if len(out) != len(counts):
        raise ValueError("Endpoint intersection altered the reference cohort")

    for end in (2019, 2020):
        out[f"delta_density_2014_{end}"] = out[f"density_{end}"] - out["density_2014"]
        for taxon in TAXA:
            out[f"delta_{taxon}_prc_2014_{end}"] = (
                out[f"{taxon}_prc_{end}"] - out[f"{taxon}_prc_2014"])
    out["delta_density_2019_2020"] = out["density_2020"] - out["density_2019"]
    out["identical_density_2019_2020"] = out["density_2019"].eq(out["density_2020"])
    return out.sort_values(["country", "NUTS_ID"]).reset_index(drop=True)


def country_correlations(full: pd.DataFrame) -> pd.DataFrame:
    """Reuse the original absolute-density rank analysis and sign thresholds."""
    rows = []
    for end in (2019, 2020):
        d = full[full.year.eq(end)].copy()
        for country, g in [("ALL", d), *list(d.groupby("country"))]:
            for taxon in TAXA:
                s = directional_summary(
                    g, country, taxon, "cumulative", 5.0, 0.001)
                _, pct_rho = exact_corr(
                    g["cumulative_pct_population_density"],
                    g[f"cumulative_delta_{taxon}_prc"])
                s.update(start_year=2014, end_year=end,
                         percentage_density_spearman_rho=pct_rho)
                rows.append(s)
    return pd.DataFrame(rows).sort_values(
        ["country", "taxon", "end_year"]).reset_index(drop=True)


def change_qa(wide: pd.DataFrame) -> pd.DataFrame:
    """Identify zero/constant 2019-20 changes, not merely undefined rho."""
    rows = []
    for country, g in [("ALL", wide), *list(wide.groupby("country"))]:
        v = g["delta_density_2019_2020"].astype(float)
        rows.append({
            "country": country,
            "n_matched_regions": len(g),
            "n_exactly_identical_density": int(g["identical_density_2019_2020"].sum()),
            "fraction_exactly_identical_density": float(g["identical_density_2019_2020"].mean()),
            "n_distinct_delta_density_rounded_6dp": int(v.round(6).nunique()),
            "min_delta_density_people_per_km2": float(v.min()),
            "median_delta_density_people_per_km2": float(v.median()),
            "max_delta_density_people_per_km2": float(v.max()),
            "sd_delta_density_people_per_km2": float(v.std(ddof=0)),
            "constant_change_2019_2020": bool(v.round(6).nunique() <= 1),
        })
    return pd.DataFrame(rows).sort_values("country").reset_index(drop=True)


def population_source_qa(full: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Audit original Eurostat/recovery provenance without replacing values."""
    raw = pd.read_parquet(path)
    need = {"NUTS_ID", "year", "population_density", "nuts_release"}
    if not need.issubset(raw.columns):
        raise ValueError(f"Missing source columns: {need - set(raw.columns)}")
    raw = raw[raw.year.isin(YEARS) & raw.nuts_release.eq(2016)].copy()
    if "matches_nuts_release" in raw.columns:
        raw = raw[raw.matches_nuts_release.fillna(False)]
    if raw.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate original Eurostat NUTS3/year observations")
    for col in ("status", "population_density_source",
                "population_density_recovered"):
        if col not in raw:
            raw[col] = pd.NA
    raw = raw[["NUTS_ID", "year", "population_density", "status",
               "population_density_source", "population_density_recovered"]]
    panel = full.loc[full.year.isin(YEARS),
        ["NUTS_ID", "country", "year", "population_density"]].merge(
            raw, on=["NUTS_ID", "year"], how="left",
            validate="one_to_one", indicator=True,
            suffixes=("_panel", "_source"))
    if not panel["_merge"].eq("both").all():
        raise ValueError("Missing original-source row in matched population panel")
    if not np.allclose(panel.population_density_panel,
                       panel.population_density_source,
                       rtol=0, atol=1e-9, equal_nan=False):
        raise ValueError("Current panel density differs from original source")
    panel["status"] = panel["status"].fillna("(not reported)").astype(str)
    panel["population_density_source_name"] = (
        panel["population_density_source"].fillna("(not reported)").astype(str))
    panel["population_density_recovered"] = (
        panel["population_density_recovered"].fillna(False).astype(bool))
    return panel


def verify_previous_2020(now: pd.DataFrame, old_path: Path) -> str:
    """Where available, require exact reproduction of the existing 2020 table."""
    if not old_path.is_file():
        return f"Previous output not present (comparison skipped): {old_path}"
    old = pd.read_csv(old_path)
    old = old[old.taxon.isin(TAXA)]
    new = now[now.end_year.eq(2020)]
    check = old[["country", "taxon", "n_regions",
                 "spearman_rho_pop_vs_delta_prc"]].merge(
        new[["country", "taxon", "n_regions",
             "spearman_rho_pop_vs_delta_prc"]],
        on=["country", "taxon"], how="outer", validate="one_to_one",
        suffixes=("_previous", "_new"), indicator=True)
    if not check["_merge"].eq("both").all():
        raise ValueError("2014-2020 country set differs from earlier analysis")
    if not check.n_regions_previous.eq(check.n_regions_new).all():
        raise ValueError("2014-2020 reference cohort size changed")
    a = check.spearman_rho_pop_vs_delta_prc_previous.to_numpy(float)
    b = check.spearman_rho_pop_vs_delta_prc_new.to_numpy(float)
    if not np.allclose(a, b, atol=1e-9, rtol=0, equal_nan=True):
        raise ValueError("2014-2020 Spearman results differ from original outputs")
    return f"2014-2020 reference correlations replicated: {len(check)} rows"


def synthetic_test() -> None:
    records = []
    for idx in range(10):
        for year in range(2014, 2021):
            d = year - 2014
            records.append({
                "NUTS_ID": f"IE{idx:03d}", "country": "IE", "year": year,
                "population_density": 20 + idx + d * (idx - 4) / 10,
                "spider_prc": .30 - .0001 * d * idx,
                "moth_prc": .40 - .0002 * d * idx,
                "cumulative_delta_population_density": d * (idx - 4) / 10,
                "cumulative_pct_population_density":
                    100 * d * (idx - 4) / 10 / (20 + idx),
                "cumulative_delta_spider_prc": -.0001 * d * idx,
                "cumulative_delta_moth_prc": -.0002 * d * idx,
            })
    panel = pd.DataFrame(records)
    w = region_endpoints(panel)
    assert len(w) == 10 and len(w.columns) > 15
    assert np.allclose(w.delta_density_2014_2019,
                       5 * (np.arange(10) - 4) / 10)
    assert not w.identical_density_2019_2020.all()
    q = change_qa(w)
    assert int(q.loc[q.country.eq("IE"), "n_matched_regions"].iloc[0]) == 10
    c = country_correlations(panel)
    assert set(c.end_year) == {2019, 2020}
    assert set(c.taxon) == set(TAXA)
    assert c.n_regions.eq(10).all()
    print("Fixed-cohort endpoint QA self-test: PASS")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    synthetic_test()
    if args.self_test:
        return

    full, annual, cumulative = calculate(inputs())
    _, previous, _ = countries_by_year(annual, cumulative)
    country_codes = previous[["NUTS_ID", "country"]]
    full = full.merge(country_codes, on="NUTS_ID", validate="many_to_one")
    if full.NUTS_ID.nunique() < 1000:
        raise ValueError("Unexpectedly small seven-year matched NUTS3 cohort")
    wide = region_endpoints(full)
    country = country_correlations(full)
    changes = change_qa(wide)
    prior = (ROOT / "figures/prc_population_country_absolute_density_2014_2020"
             / "country_directional_summary_2014_2020.csv")
    baseline_check = verify_previous_2020(country, prior)

    pop_source = population_source_qa(
        full, ROOT / "data/processed/nuts3_population_density.parquet")
    source_summary = pop_source.groupby(
        ["country", "year", "population_density_source_name",
         "status", "population_density_recovered"], dropna=False
    ).size().reset_index(name="n_region_years").sort_values(
        ["country", "year", "population_density_source_name", "status"])

    out = args.out if args.out.is_absolute() else ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    wide.to_csv(out / "region_paired_2014_2019_2020.csv", index=False)
    wide[wide.country.eq("IE")].to_csv(
        out / "ireland_region_paired_2014_2019_2020.csv", index=False)
    country.to_csv(out / "country_2014_2019_vs_2020_correlations.csv",
                   index=False)
    changes.to_csv(out / "population_2019_2020_change_qa_by_country.csv",
                   index=False)
    source_summary.to_csv(out / "eurostat_source_status_by_country_year.csv",
                          index=False)
    pop_source.to_csv(out / "eurostat_matched_population_provenance.csv",
                      index=False)
    summary = {
        "reference_cohort": "same complete 2014-2020 seven-year NUTS2016 panel",
        "n_regions": len(wide),
        "years_compared": ["2014-2019", "2014-2020"],
        "primary_predictor": "absolute change in people per km2",
        "primary_response": "absolute spider/moth PRC change",
        "fixed_sign_thresholds": {
            "density_people_per_km2": 5.0, "absolute_prc": 0.001},
        "2014_2020_reproduction": baseline_check,
        "population_2019_2020_constant_change_countries": changes.loc[
            changes.country.ne("ALL") & changes.constant_change_2019_2020,
            "country"].tolist(),
        "ireland_n_matched_regions": int(wide.country.eq("IE").sum()),
        "caveats": [
            "Population values are not independently validated against external yearly records",
            "A constant 2019-2020 difference may reflect source construction or actual low variation",
            "Country Spearman values are descriptive, not spatially adjusted p-values",
            "Same NUTS2016 regions used for both endpoints; thresholds unchanged",
            "No causal or LED-spectral attribution is established",
        ],
    }
    (out / "qa_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("\n" + baseline_check)
    print(f"Common seven-year cohort: {len(wide)} regions")
    print("\nSelected within-country absolute-density vs absolute PRC Spearman:")
    show = country[country.country.isin(["ALL", "TR", "DE", "FR", "IT", "IE"])]
    print(show[["country", "taxon", "end_year", "n_regions",
                "spearman_rho_pop_vs_delta_prc",
                "percentage_density_spearman_rho"]].to_string(index=False))
    print("\n2019-2020 density diagnostics for selected countries:")
    print(changes[changes.country.isin(["ALL", "TR", "DE", "FR", "IT", "IE"])].to_string(index=False))
    print(f"\nIreland matched regions: {summary['ireland_n_matched_regions']}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
