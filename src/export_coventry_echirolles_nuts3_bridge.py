#!/usr/bin/env python3
"""Read-only NUTS2016 annual radiance export for Coventry and Isère, 2014–24.

Compares independently dated *municipal* lighting changes with the FROZEN
European annual NUTS3 VNL V2 radiance history (not monthly VNL V1 averages).

Scientific safeguards:
- use fixed NUTS2016 annual parquet, same as publication Figure 1
- reuse Figure 1's verified locked 2021 native-pixel snow overlay cache,
  abort if its identity/audit is absent or stale, never rebuild it here
- validate both endpoints against locked 2014/2024 radiance+PRC CSV
- never apply the 2017 zero offset twice
- annual 2022/2024 scenes contain pre- and post-intervention months
- compare municipality-only EOG V1 contrasts SEPARATELY (not fused means)
- do not change models, PRC transfer, source data, locked figures or cohorts

Run:
  python src/export_coventry_echirolles_nuts3_bridge.py --self-test
  python src/export_coventry_echirolles_nuts3_bridge.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from project_config import PROJECT_ROOT, load_config
from extract_nuts3_radiance import discover_rasters
from figure1_2021_locked_snow_radiance import _lock_info
from plot_publication_figure1_radiance import get_locked_snow_2021

ROOT = PROJECT_ROOT
ANNUAL = ROOT / "data/processed/fixed_nuts2016_viirs_2014_2024.parquet"
ENDPOINT = (ROOT / "figures/prc_radiance_spider_moth_2014_2024/"
            "nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv")
SNOW_CACHE = (ROOT / "data/processed/publication_figures/"
              "fig1_2021_locked_snow_radiance.csv")
MUNICIPAL = (ROOT / "data/processed/local_lighting_interventions/"
             "five_case_monthly/five_case_contrast_summary.csv")
OUT = ROOT / "data/processed/local_lighting_interventions/nuts3_event_bridge"
YEARS = tuple(range(2014, 2025))
REGIONS = {
    "UKG33": {
        "label": "Coventry", "municipality_case": "coventry",
        "event": "part-night switch-offs", "decision": "2024-02-20",
        "rollout_start": "2024-05", "rollout_end": "2024-07",
        "first_full_post_month": "2024-08",
        "affected_share": "about 60% eligible city streetlights; not all sources",
    },
    "FRK24": {
        "label": "Isère (includes Échirolles)", "municipality_case": "echirolles",
        "event": "Échirolles partial night switch-offs",
        "decision": "2022-09-26",
        "rollout_start": "2022-07", "rollout_end": "2022-12",
        "first_full_post_month": "2023-01",
        "affected_share": "only the Échirolles municipal subset of Isère",
    },
}


def strict_snow_cache():
    """Validate cache *before* calling shared reader; never rebuild sources."""
    if not SNOW_CACHE.is_file() or not SNOW_CACHE.with_suffix(".json").is_file():
        raise FileNotFoundError(
            f"Missing verified Figure 1 2021 snow cache/audit: {SNOW_CACHE}. "
            "Regenerate only via the separate established Figure 1 workflow; "
            "this read-only exporter will not create or recalibrate it."
        )
    cfg = load_config()
    lock, offset = _lock_info(cfg)
    patch = ROOT / lock["patch"]
    original = discover_rasters(cfg, [2021])[2021]
    expected = {
        "status": "locked",
        "correction_model": lock["model_id"],
        "original_annual_vnl": str(original.resolve()),
        "patch": str(patch.resolve()),
        "patch_bytes": patch.stat().st_size,
        "patch_mtime_ns": patch.stat().st_mtime_ns,
        "original_annual_bytes": original.stat().st_size,
        "original_annual_mtime_ns": original.stat().st_mtime_ns,
        "pixel_offset_nw_cm2_sr": offset,
    }
    audit = json.loads(
        SNOW_CACHE.with_suffix(".json").read_text(encoding="utf-8"))
    bad = [k for k, v in expected.items() if audit.get(k) != v]
    if bad:
        raise RuntimeError(
            "Figure 1 2021 snow cache changed or is unverified: "
            + ",".join(bad) +
            ". Do not use an unvalidated 2021 value or regenerate silently."
        )
    snow, _ = get_locked_snow_2021(SNOW_CACHE)
    return snow


def endpoint_series(frame):
    required = ("NUTS_ID", "radiance_2014", "radiance_2024",
                "radiance_change_abs", "radiance_change_pct",
                "spider_prc_2014", "spider_prc_2024",
                "moth_prc_2014", "moth_prc_2024")
    missing = set(required) - set(frame)
    if missing:
        raise ValueError(f"Missing frozen endpoint columns: {sorted(missing)}")
    e = frame.loc[frame.NUTS_ID.isin(REGIONS), required].copy()
    if len(e) != len(REGIONS) or e.NUTS_ID.duplicated().any():
        raise ValueError("Missing/duplicate NUTS3 endpoint region")
    for row in e.itertuples(index=False):
        if not np.isclose(row.radiance_2024 - row.radiance_2014,
                          row.radiance_change_abs, atol=1e-6):
            raise ValueError(f"Endpoint identity inconsistent: {row.NUTS_ID}")
    return e


def annual_series(source, endpoints, snow):
    needed = ("NUTS_ID", "year", "radiance_mean_corrected")
    if set(needed) - set(source):
        raise ValueError("Annual parquet does not match locked Figure 1 schema")
    a = source.loc[source.NUTS_ID.isin(REGIONS), needed].copy()
    a["year"] = a.year.astype(int)
    if (len(a) != len(REGIONS) * len(YEARS) or
            a.duplicated(["NUTS_ID", "year"]).any() or
            set(a.year) != set(YEARS) or
            a.groupby("NUTS_ID").year.nunique().ne(len(YEARS)).any()):
        raise ValueError("Both NUTS2016 regions must have all 11 annual rows")
    a["radiance_mean_corrected"] = pd.to_numeric(
        a.radiance_mean_corrected, errors="raise")
    if (not np.isfinite(a.radiance_mean_corrected).all() or
            (a.radiance_mean_corrected < 0).any()):
        raise ValueError("Non-finite/negative annual calibrated radiance")
    z = snow.loc[snow.NUTS_ID.isin(REGIONS), [
        "NUTS_ID", "radiance_2021_original_calibrated",
        "radiance_2021_snow_corrected"
    ]]
    a = a.merge(z, on="NUTS_ID", validate="many_to_one")
    idx = a.year.eq(2021)
    delta = (a.loc[idx, "radiance_mean_corrected"].to_numpy(float) -
             a.loc[idx, "radiance_2021_original_calibrated"].to_numpy(float))
    allowable = (2e-5 * np.maximum(
        1., np.abs(a.loc[idx, "radiance_mean_corrected"].to_numpy(float))))
    if not np.all(np.abs(delta) <= allowable):
        raise ValueError("2021 original annual radiance disagrees with verified snow cache")
    a["radiance_locked_nw_cm2_sr"] = a.radiance_mean_corrected.copy()
    a.loc[idx, "radiance_locked_nw_cm2_sr"] = (
        a.loc[idx, "radiance_2021_snow_corrected"].to_numpy(float))
    a["snow_corrected_2021"] = idx
    a = a.merge(endpoints, on="NUTS_ID", validate="many_to_one")
    for year, col in ((2014, "radiance_2014"), (2024, "radiance_2024")):
        rows = a[a.year.eq(year)]
        bad = np.abs(rows.radiance_locked_nw_cm2_sr - rows[col])
        if (bad > 0.025).any():
            raise ValueError(f"{year} annual parquet vs published endpoint "
                             f"mismatch: {rows.loc[bad>0.025, 'NUTS_ID'].tolist()}")
    a["region_name"] = a.NUTS_ID.map(
        {key: value["label"] for key, value in REGIONS.items()})
    a.sort_values(["NUTS_ID", "year"], inplace=True)
    a["delta_previous_year_nw"] = (
        a.groupby("NUTS_ID").radiance_locked_nw_cm2_sr.diff())
    a["delta_previous_year_pct"] = (
        a.groupby("NUTS_ID").radiance_locked_nw_cm2_sr.pct_change(
            fill_method=None) * 100)
    a["change_since_2014_nw"] = (
        a.radiance_locked_nw_cm2_sr - a.radiance_2014)
    a["change_since_2014_pct"] = (
        a.change_since_2014_nw / a.radiance_2014 * 100)
    return a


def event_summary(a, municipal):
    # The 2022 and 2024 intervention years are mixes of pre/post months.
    events = [
        ("UKG33", 2022, 2023, "pre-2024 regional comparison"),
        ("UKG33", 2023, 2024, "includes May–July rollout in outcome year"),
        ("FRK24", 2021, 2022, "2021 snow patched; 2022 partly treated in Échirolles"),
        ("FRK24", 2022, 2023, "Échirolles operational hours + 2023 LEDs, mixed"),
        ("FRK24", 2023, 2024, "further municipality LED changes, mixed"),
    ]
    x = a.set_index(["NUTS_ID", "year"])
    rows = []
    for code, before, after, status in events:
        source = float(x.loc[(code, before), "radiance_locked_nw_cm2_sr"])
        target = float(x.loc[(code, after), "radiance_locked_nw_cm2_sr"])
        meta = REGIONS[code]
        rows.append({
            "NUTS_ID": code,
            "region_name": meta["label"],
            "comparison": f"{before}->{after}",
            "annual_before_nw": source,
            "annual_after_nw": target,
            "delta_nw": target-source,
            "delta_pct": 100*(target/source-1),
            "local_intervention": meta["event"],
            "municipal_rollout_start": meta["rollout_start"],
            "municipal_rollout_end": meta["rollout_end"],
            "first_full_month_after_rollout": meta["first_full_post_month"],
            "status": status,
            "municipality_share_warning": meta["affected_share"],
            "annual_data_kind": "LOCKED EOG annual VNL V2; 2021 snow corrected",
            "not_causal": True,
        })
    out = pd.DataFrame(rows)
    if municipal is not None:
        # Keep monthly VNL V1 in separate fields with distinct units,
        # geography and baseline; DO NOT subtract it from annual VNL V2.
        m = municipal.copy()
        if {"case","kind","period_delta_percent","period_delta_absolute"}.issubset(m):
            m = m[m.kind.eq("primary")].set_index("case")
            out["municipal_v1_primary_delta_pct"] = [
                m.loc[REGIONS[code]["municipality_case"],
                      "period_delta_percent"] for code in out.NUTS_ID]
            out["municipal_v1_primary_delta_nw"] = [
                m.loc[REGIONS[code]["municipality_case"],
                      "period_delta_absolute"] for code in out.NUTS_ID]
    return out


def self_test():
    endpoints = pd.DataFrame({
        "NUTS_ID": list(REGIONS), "radiance_2014": [20., 4.],
        "radiance_2024": [15., 3.], "radiance_change_abs": [-5., -1.],
        "radiance_change_pct": [-25., -25.],
        "spider_prc_2014": [.2, .4], "spider_prc_2024": [.3, .5],
        "moth_prc_2014": [.3, .6], "moth_prc_2024": [.4, .7]})
    raw = pd.DataFrame([
        {"NUTS_ID": code, "year": year,
         "radiance_mean_corrected": (20 if code == "UKG33" else 4)
         + (year-2014)*(-.5 if code == "UKG33" else -.1)}
        for code in REGIONS for year in YEARS])
    snow = pd.DataFrame({
        "NUTS_ID": list(REGIONS),
        "radiance_2021_original_calibrated": [16.5, 3.3],
        "radiance_2021_snow_corrected": [16.4, 3.2]})
    result = annual_series(raw, endpoints, snow)
    assert len(result) == 22
    assert np.isclose(result.loc[
        result.NUTS_ID.eq("UKG33") & result.year.eq(2021),
        "radiance_locked_nw_cm2_sr"].iloc[0], 16.4)
    assert np.isclose(result.loc[
        result.NUTS_ID.eq("FRK24") & result.year.eq(2024),
        "change_since_2014_pct"].iloc[0], -25)
    events = event_summary(result, None)
    assert len(events) == 5
    print("Self-test passed: 11 years, endpoint identity, locked 2021 splice, "
          "year-on-year and pre/post output")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return
    if not ANNUAL.is_file() or not ENDPOINT.is_file():
        raise FileNotFoundError(
            f"Need locked Figure 1 source files: {ANNUAL}, {ENDPOINT}")
    snow = strict_snow_cache()
    source = pd.read_parquet(
        ANNUAL, columns=["NUTS_ID","year","radiance_mean_corrected"])
    endpoints = endpoint_series(pd.read_csv(
        ENDPOINT, dtype={"NUTS_ID": str}))
    annual = annual_series(source, endpoints, snow)
    municipal = (pd.read_csv(MUNICIPAL) if MUNICIPAL.is_file() else None)
    bridge = event_summary(annual, municipal)
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["NUTS_ID","region_name","year",
            "radiance_locked_nw_cm2_sr","delta_previous_year_nw",
            "delta_previous_year_pct","change_since_2014_nw",
            "change_since_2014_pct","snow_corrected_2021"]
    annual[cols].to_csv(
        OUT / "coventry_isere_locked_nuts3_annual_2014_2024.csv",
        index=False)
    bridge.to_csv(OUT / "coventry_echirolles_nuts3_event_match.csv",
                  index=False)
    print("\nLOCKED ANNUAL NUTS3 HISTORY (nW cm^-2 sr^-1):")
    print(annual.pivot(index="year", columns="NUTS_ID",
                       values="radiance_locked_nw_cm2_sr").round(3).to_string())
    print("\nEVENT-YEAR NUTS3 COMPARISONS:")
    print(bridge[["region_name","comparison","annual_before_nw",
                  "annual_after_nw","delta_pct","status"]].
          round(3).to_string(index=False))
    print("\nOutputs:", OUT)
    print("Original annual, 2021 patch, endpoint and municipal VIIRS tables "
          "were READ ONLY and unchanged.")


if __name__ == "__main__":
    main()
