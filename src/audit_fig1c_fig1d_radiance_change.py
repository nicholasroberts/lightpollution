#!/usr/bin/env python3
"""Read-only audit of Figure 1C Europe trajectory vs Figure 1D histogram.

These figures answer different questions:
  Fig1C: percentage change of fixed-2014-area-weighted European mean radiance.
  Fig1D: equal-region distribution, mean and median of NUTS3 percentages.

Reads the exact two source exports of src/plot_publication_figure1_radiance.py,
plus the endpoint summary JSON, and reports the mathematical reconciliation.
Does not regenerate, correct, impute or change any scientific table/figure.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from project_config import PROJECT_ROOT as ROOT

ENDPOINT = ROOT / (
    "figures/prc_radiance_spider_moth_2014_2024/"
    "nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv"
)
FIGDIR = ROOT / "figures/publication/figure1_radiance"


def result_from_files(endpoint_path: Path, figdir: Path) -> str:
    series_path = figdir / "fig1C_europe_annual_radiance_support.csv"
    summary_path = figdir / "fig1_provenance_and_summary.json"
    endpoints = pd.read_csv(endpoint_path, dtype={"NUTS_ID": str})
    annual = pd.read_csv(series_path)
    summary = json.loads(summary_path.read_text(encoding="utf8"))
    cols = {"NUTS_ID", "radiance_2014", "radiance_2024",
            "valid_area_km2_2014"}
    if not cols.issubset(endpoints):
        raise ValueError(f"Endpoint fields missing: {cols-set(endpoints)}")
    if not {"year", "area_weighted_radiance_nw"}.issubset(annual):
        raise ValueError("Figure 1C annual support CSV lacks expected columns")
    if endpoints.NUTS_ID.isna().any() or endpoints.NUTS_ID.duplicated().any():
        raise ValueError("Endpoint NUTS3 keys missing/duplicated")
    if len(endpoints)!=1522:
        raise ValueError(f"Expected 1,522 NUTS3 endpoints, got {len(endpoints)}")
    if annual.year.duplicated().any() or set(annual.year)!={*range(2014,2025)}:
        raise ValueError("Figure 1C must have 2014–2024 exactly once")
    e=endpoints[["NUTS_ID","radiance_2014","radiance_2024",
                 "valid_area_km2_2014"]].copy()
    for k in ["radiance_2014","radiance_2024","valid_area_km2_2014"]:
        e[k]=pd.to_numeric(e[k],errors="raise")
        if not np.isfinite(e[k]).all():
            raise ValueError(f"Nonfinite value {k}")
    if (e[["radiance_2014","radiance_2024"]]<0).any().any():
        raise ValueError("Negative radiance in endpoint data")
    if (e.valid_area_km2_2014<=0).any():
        raise ValueError("Zero/negative area weights")
    r0=e.radiance_2014.to_numpy(float)
    r1=e.radiance_2024.to_numpy(float)
    area=e.valid_area_km2_2014.to_numpy(float)
    positive=r0>0.
    pct=np.full(len(e),np.nan)
    pct[positive]=100*(r1[positive]/r0[positive]-1)
    finite=pct[np.isfinite(pct)]
    mean=float(finite.mean())
    median=float(np.median(finite))
    er0=float(np.average(r0,weights=area))
    er1=float(np.average(r1,weights=area))
    endpoint_pct=100*(er1/er0-1)
    # Algebraic identity for positive 2014 baseline regions; exact-zero
    # 2014 baselines still add new 2024 light to the European aggregate.
    denom=float(np.sum(area*r0))
    nonzero_contribution=float(np.sum(
        area[positive]*r0[positive]*pct[positive]/100.
    ))
    zero_contribution=float(np.sum(area[~positive]*r1[~positive]))
    identity=100*(nonzero_contribution+zero_contribution)/denom
    if not np.isclose(identity,endpoint_pct,atol=1e-9,rtol=1e-9):
        raise ValueError("Europe / weighted NUTS3 percentage identity failed")
    a=annual.set_index("year")
    c14=float(a.loc[2014,"area_weighted_radiance_nw"])
    c24=float(a.loc[2024,"area_weighted_radiance_nw"])
    chart_pct=100*(c24/c14-1)
    diff14=c14-er0
    diff24=c24-er1
    if "mean_nuts3_percent_change_equal_region_weight" in summary:
        if not np.isclose(summary["mean_nuts3_percent_change_equal_region_weight"],
                          mean,rtol=1e-8,atol=1e-8):
            raise ValueError("Output JSON NUTS3 mean contradicts original endpoints")
    if "median_nuts3_percent_change_equal_region_weight" in summary:
        if not np.isclose(summary["median_nuts3_percent_change_equal_region_weight"],
                          median,rtol=1e-8,atol=1e-8):
            raise ValueError("Output JSON NUTS3 median contradicts original endpoints")
    if "europe_change_of_area_weighted_radiance_percent" in summary:
        if not np.isclose(summary["europe_change_of_area_weighted_radiance_percent"],
                          endpoint_pct,rtol=1e-8,atol=1e-8):
            raise ValueError("Output JSON Europe percentage contradicts endpoint data")
    lines=[
        "FIGURE 1C VERSUS FIGURE 1D — READ-ONLY RECONCILIATION",
        "="*66,
        f"Figure 1C 2014 Europe weighted radiance: {c14:.10f}",
        f"Figure 1C 2024 Europe weighted radiance: {c24:.10f}",
        f"Figure 1C Europe 2014-to-2024 change: {chart_pct:+.6f}%",
        "",
        f"Figure 1D endpoints: {len(e)} NUTS3",
        f"Figure 1D percentage-defined (R2014 > 0): {len(finite)}",
        f"Figure 1D exact-zero baselines (undefined regional %): {int((~positive).sum())}",
        f"Figure 1D equal-NUTS3 arithmetic mean percentage: {mean:+.6f}%",
        f"Figure 1D median regional percentage: {median:+.6f}%",
        "",
        f"Endpoints 2014 weighted Europe mean: {er0:.10f}",
        f"Endpoints 2024 weighted Europe mean: {er1:.10f}",
        f"Endpoints percentage change: {endpoint_pct:+.6f}%",
        f"Area- and starting-radiance-weighted percentage contribution: "
        f"{100*nonzero_contribution/denom:+.6f} percentage points",
        f"Extra contribution from zero-2014 baseline regions: "
        f"{100*zero_contribution/denom:+.6f} percentage points",
        f"Exact weighted identity check: {identity:+.6f}%",
        "",
        "Annual Figure1C table vs independent endpoint CSV:",
        f"2014 absolute difference: {diff14:+.10f} nW cm^-2 sr^-1",
        f"2024 absolute difference: {diff24:+.10f} nW cm^-2 sr^-1",
        f"Change difference (Fig1C - endpoints): {chart_pct-endpoint_pct:+.6f} pp",
        "",
        "The histogram is an equal-region distribution, not an estimate",
        "of percentage change in Europe's aggregate mean. A histogram",
        "mode is not necessarily its arithmetic mean or median.",
        "The 2021 snow adjustment cannot change either 2014/2024 endpoint.",
        "All sources remain unchanged."
    ]
    return "\n".join(lines)+"\n"


def self_test() -> None:
    # Fixed area example: one dim region +100%, one bright region -10%
    # => equal-region mean +45% but Europe (-8/11)=-72.73%.
    a=np.array([1.,10.])
    r0=np.array([1.,10.])
    r1=np.array([2.,9.])
    percent=100*(r1/r0-1)
    aggregate=100*((a*r1).sum()/(a*r0).sum()-1)
    weight=(a*r0)/(a*r0).sum()
    assert np.isclose(aggregate,(weight*percent).sum())
    assert not np.isclose(float(percent.mean()),aggregate)
    print("Figure 1C vs Figure 1D estimand self-test: PASS")


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--endpoints",type=Path,default=ENDPOINT)
    p.add_argument("--figdir",type=Path,default=FIGDIR)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    self_test()
    if args.self_test:
        return
    e=args.endpoints if args.endpoints.is_absolute() else ROOT/args.endpoints
    d=args.figdir if args.figdir.is_absolute() else ROOT/args.figdir
    print(result_from_files(e,d))


if __name__=="__main__":
    main()
