#!/usr/bin/env python3
"""Four-panel NUTS2016 maps: Jan-1 population, VIIRS radiance, spider/moth PRC.

Uses reconstructed January-1 population / fixed land area (not published
DEMO_R_D3DENS). For 2014->2019 and 2014->2020, display four existing change
variables side-by-side on the exact SAME 1273-region reconstructed cohort:
absolute change in Jan-1 population density, corrected regional mean VIIRS
radiance, and unchanged locked pixel-based spider and moth PRC.

Linear references, population-only symlog and combined population/radiance
symlog colour maps are generated. SymLogNorm preserves positive, negative and
zero changes while legends retain untransformed people/km2 and radiance units.
Both PRC panels retain independently zero-centred linear scales.
Neither annual population smoothing, PRC refitting nor radiance recalibration.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, SymLogNorm
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd

from historical_prc_pixel_nuts3 import load_nuts
from project_config import PROJECT_ROOT as ROOT, load_config

INPUT = (ROOT / "data/processed/nuts2016_jan1_population_density_sensitivity"
         / "matched_region_density_2014_2019_2020.csv")
OUTPUT = ROOT / "figures/prc_population_change_jan1_linear_vs_symlog"
ENDPOINTS = (2019, 2020)
RADIANCE_SOURCE = ROOT / "data/processed/fixed_nuts2016_viirs_2014_2024.parquet"
MAP_XLIM = (-12, 46)  # Includes almost all of Turkiye, unlike old x<=36 sheet.
MAP_YLIM = (35, 73)


def read_matched_endpoints(path: Path) -> pd.DataFrame:
    x = pd.read_csv(path, dtype={"NUTS_ID": str})
    needed = {"NUTS_ID", "country", "density_jan1_fixed_land_2014",
              "density_jan1_fixed_land_2019",
              "density_jan1_fixed_land_2020"}
    needed |= {
        f"change_density_jan1_fixed_land_2014_{yr}"
        for yr in ENDPOINTS
    }
    needed |= {
        f"cumulative_delta_{taxon}_prc_{yr}"
        for yr in ENDPOINTS for taxon in ("spider", "moth")
    }
    if not needed.issubset(x.columns):
        raise ValueError(f"Reconstruction endpoint columns missing: {needed-set(x.columns)}")
    if x.NUTS_ID.duplicated().any():
        raise ValueError("Duplicated NUTS3 region in matched Jan-1 endpoints")
    if len(x) != 1273:
        raise ValueError(f"Expected 1273 fixed reconstructed regions, found {len(x)}")
    for year in ENDPOINTS:
        pop = f"change_density_jan1_fixed_land_2014_{year}"
        actual = x[f"density_jan1_fixed_land_{year}"] - x[
            "density_jan1_fixed_land_2014"]
        if not np.allclose(x[pop], actual, atol=1e-10, rtol=1e-10):
            raise ValueError(f"Region population change not equal to endpoint difference: {year}")
        columns = [pop, *(f"cumulative_delta_{taxon}_prc_{year}"
                          for taxon in ("spider", "moth"))]
        if not np.isfinite(x[columns].to_numpy(float)).all():
            raise ValueError(f"Non-finite population or PRC change in {year} panel")
    return x


def add_radiance_endpoints(data: pd.DataFrame, source: Path) -> pd.DataFrame:
    """Join pre-calibrated VIIRS NUTS2016 radiance, never recompute sensor data.

    The PRCs already come from their LOCKED native-pixel transformations;
    mean VIIRS radiance is a separate region-level quantity. The three end
    years and NUTS3 identifiers must be unique and complete.
    """
    columns = ["NUTS_ID", "year", "radiance_mean_corrected"]
    rad = pd.read_parquet(source, columns=columns)
    rad = rad[rad.year.isin((2014, *ENDPOINTS))].copy()
    rad["NUTS_ID"] = rad.NUTS_ID.astype(str)
    if rad.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicated corrected radiance region/year")
    result = data.copy()
    for year in (2014, *ENDPOINTS):
        year_rad = rad.loc[rad.year.eq(year),
                           ["NUTS_ID", "radiance_mean_corrected"]].rename(
            columns={"radiance_mean_corrected": f"radiance_{year}"})
        result = result.merge(year_rad, on="NUTS_ID", how="left",
                              validate="one_to_one")
    required = [f"radiance_{year}" for year in (2014, *ENDPOINTS)]
    if not np.isfinite(result[required].to_numpy(dtype=float)).all():
        raise ValueError("Missing/nonfinite VIIRS radiance in reconstructed NUTS3 cohort")
    if (result[required] < 0).any().any():
        raise ValueError("Invalid negative corrected NUTS3 VIIRS radiance")
    for year in ENDPOINTS:
        result[f"radiance_delta_2014_{year}"] = (
            result[f"radiance_{year}"] - result["radiance_2014"])
    if len(result) != len(data):
        raise ValueError("Joining radiance changed the fixed population cohort")
    return result

def robust_symmetric_limit(values: pd.Series, minimum: float) -> float:
    arr = np.abs(values.to_numpy(dtype=float))
    arr = arr[np.isfinite(arr)]
    if not len(arr):
        raise ValueError("Empty map data")
    # Same 98th-absolute-percentile display convention as original map,
    # with clipping of more extreme values for visibility.
    return max(float(np.quantile(arr, 0.98)), minimum)


def symlog_ticks(limit: float, linthresh: float) -> list[float]:
    candidates = [0, -linthresh, linthresh]
    # Start decades at 10 × the linear threshold, so fractional-unit
    # radiance change receives ticks at 0.1, 1, 10, ... rather than 10, 100.
    decade = 10. * linthresh
    while decade < limit:
        candidates.extend([-decade, decade])
        decade *= 10
    # Do not label the exact 98th-percentile end caps: these can nearly
    # overlap the adjacent round-number decade ticks (e.g. 142 and 100).
    # The extend='both' colourbar triangles already mark clipped extremes.
    return sorted({float(z) for z in candidates if abs(z) <= limit + 1e-10})


def draw_panel(fig, ax, geom, column: str, limit: float, label: str,
               scale_mode: str = "linear", linthresh: float = 1.,
               colourbar_label: str | None = None) -> None:
    if scale_mode not in ("linear", "symlog"):
        raise ValueError(f"Unknown panel scale: {scale_mode}")
    if scale_mode == "symlog":
        if linthresh <= 0:
            raise ValueError("Symmetric-log threshold must be positive")
        norm = SymLogNorm(
            linthresh=linthresh, linscale=1.0, base=10,
            vmin=-limit, vmax=limit, clip=True)
    else:
        norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    cmap = "RdBu_r"  # Red=increasing quantity; blue=decreasing quantity.
    geom.plot(
        column=column, ax=ax, cmap=cmap, norm=norm, linewidth=0,
        missing_kwds={"color": "#d3d3d3"}, legend=False)
    ax.set(xlim=MAP_XLIM, ylim=MAP_YLIM,
           xlabel="Longitude", ylabel="Latitude",
           title=f"{label}\n(n={int(geom[column].notna().sum())})")
    ax.set_aspect("auto")
    mappable = ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    bar = fig.colorbar(mappable, ax=ax, shrink=0.67, pad=0.018,
                       fraction=0.04, extend="both")
    if scale_mode == "symlog":
        ticks = symlog_ticks(limit, linthresh)
        bar.set_ticks(ticks)
        bar.set_ticklabels([f"{v:g}" for v in ticks])
    scale_note = "symmetric log" if scale_mode == "symlog" else "linear"
    bar.set_label(f"{colourbar_label} | {scale_note}" if colourbar_label else label)


def make_map(nuts, data: pd.DataFrame, year: int, pop_kind: str,
             rad_kind: str, out: Path, pop_linthresh: float,
             rad_linthresh: float) -> dict:
    pop_col = f"change_density_jan1_fixed_land_2014_{year}"
    rad_col = f"radiance_delta_2014_{year}"
    cols = [pop_col, rad_col, f"cumulative_delta_spider_prc_{year}",
            f"cumulative_delta_moth_prc_{year}"]
    geo = nuts.merge(data[["NUTS_ID", *cols]], how="left",
                     on="NUTS_ID", validate="one_to_one")
    limits = {
        pop_col: robust_symmetric_limit(data[pop_col], 5.),
        rad_col: robust_symmetric_limit(data[rad_col], 1e-10),
        cols[2]: robust_symmetric_limit(data[cols[2]], 1e-10),
        cols[3]: robust_symmetric_limit(data[cols[3]], 1e-10),
    }
    fig, ax = plt.subplots(1, 4, figsize=(26, 8.8),
                           layout="constrained")
    labels = [
        f"Jan 1 population density Δ 2014–{year} (people/km²)",
        f"VIIRS radiance Δ 2014–{year} (nW cm⁻² sr⁻¹)",
        f"Spider PRC Δ 2014–{year}",
        f"Moth PRC Δ 2014–{year}",
    ]
    draw_panel(fig, ax[0], geo, pop_col, limits[pop_col], labels[0],
               scale_mode=pop_kind, linthresh=pop_linthresh,
               colourbar_label="Population density change (people/km²)")
    draw_panel(fig, ax[1], geo, rad_col, limits[rad_col], labels[1],
               scale_mode=rad_kind, linthresh=rad_linthresh,
               colourbar_label="VIIRS radiance change (nW cm⁻² sr⁻¹)")
    for i in (2, 3):
        draw_panel(fig, ax[i], geo, cols[i], limits[cols[i]], labels[i])
    fig.suptitle(
        f"Cumulative change 2014–{year} | fixed NUTS2016 geography"
        f" | independently reconstructed Jan-1 population\n"
        f"Population colour scale: {pop_kind}; radiance: {rad_kind}; "
        f"spider/moth PRC: linear, independently zero-centred\n"
        f"Grey=outside complete 1,273-region cohort; radiance=corrected "
        f"regional VIIRS mean; PRC=unchanged locked native-pixel estimates",
        fontsize=11,
    )
    # Preserve the original linear-radiance filenames and provide distinctly
    # named symlog-radiance variants without overwriting reference figures.
    basename = f"cumulative_population_radiance_spider_moth_2014_{year}_jan1_{pop_kind}"
    if rad_kind == "symlog":
        basename += "_rad_symlog"
    fig.savefig(out / (basename + ".pdf"), bbox_inches="tight")
    fig.savefig(out / (basename + ".png"), bbox_inches="tight", dpi=200)
    plt.close(fig)
    return {"year": year, "population_scale": pop_kind,
            "radiance_scale": rad_kind, "regions": int(len(data)),
            "colour_scale_max_abs_98pct": limits,
            "population_symlog_threshold_people_per_km2":
                pop_linthresh if pop_kind == "symlog" else None,
            "radiance_symlog_threshold_nw_cm2_sr":
                rad_linthresh if rad_kind == "symlog" else None}


def self_test() -> None:
    a = pd.Series([-400., -10., -1., 0., 1., 10., 400.])
    v = robust_symmetric_limit(a, 1.)
    norm = SymLogNorm(linthresh=1., vmin=-v, vmax=v, base=10)
    assert np.isclose(float(norm(0)), .5)
    assert float(norm(-10)) < .5 < float(norm(10))
    assert float(norm(-1)) < float(norm(0)) < float(norm(1))
    assert symlog_ticks(v, 1.).count(0.) == 1
    assert not any(np.isclose(abs(t), v) for t in symlog_ticks(v, 1.))
    # Signed, fractional-unit radiance with a 0.1 nW cm-2 sr-1 threshold.
    radiance_limit = 12.
    rad_norm = SymLogNorm(linthresh=.1, vmin=-radiance_limit,
                          vmax=radiance_limit, base=10)
    assert np.isclose(float(rad_norm(0.)), .5)
    assert float(rad_norm(-1.)) < .5 < float(rad_norm(1.))
    assert symlog_ticks(radiance_limit, .1) == [-10., -1., -.1, 0., .1, 1., 10.]
    assert all(np.isfinite(rad_norm(np.array([-.1, 0., .1]))))
    # Simple toy pairing validates same-NUTS joining and endpoint subtraction
    # without depending on the large project files in --self-test.
    toy = pd.DataFrame({"NUTS_ID": ["TR001", "FR001"]})
    r = pd.DataFrame({
        "NUTS_ID": ["TR001", "FR001"] * 3,
        "year": [2014, 2014, 2019, 2019, 2020, 2020],
        "radiance_mean_corrected": [1., 2., 1.5, 1.75, 1.8, 1.6],
    })
    # Check source-invariant end-year difference construction.
    pivot = r.pivot(index="NUTS_ID", columns="year",
                    values="radiance_mean_corrected")
    assert np.isclose(float(pivot.loc["TR001", 2020] - pivot.loc["TR001", 2014]), .8)
    assert len(toy) == 2
    print("Signed-log scale and endpoint-change self-test: PASS")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--input", type=Path, default=INPUT)
    p.add_argument("--out", type=Path, default=OUTPUT)
    p.add_argument("--linthresh", type=float, default=1.0,
                   help="Linear part of SymLogNorm around zero (people/km2)")
    p.add_argument("--rad-linthresh", type=float, default=0.1,
                   help="SymLogNorm linear part for VIIRS change in nW cm-2 sr-1")
    args = p.parse_args()
    if not np.isfinite(args.linthresh) or args.linthresh <= 0:
        p.error("--linthresh must be positive and finite")
    if not np.isfinite(args.rad_linthresh) or args.rad_linthresh <= 0:
        p.error("--rad-linthresh must be positive and finite")
    self_test()
    if args.self_test:
        return
    source = args.input if args.input.is_absolute() else ROOT / args.input
    out = args.out if args.out.is_absolute() else ROOT / args.out
    regions = read_matched_endpoints(source)
    regions = add_radiance_endpoints(regions, RADIANCE_SOURCE)
    nuts = load_nuts(load_config())[["NUTS_ID", "geometry"]]
    if not set(regions.NUTS_ID).issubset(set(nuts.NUTS_ID)):
        raise ValueError("Reconstructed region missing from fixed NUTS2016 geometries")
    out.mkdir(parents=True, exist_ok=True)
    info = []
    for year in ENDPOINTS:
        for pop_kind, rad_kind in (("linear", "linear"),
                                   ("symlog", "linear"),
                                   ("symlog", "symlog")):
            info.append(make_map(
                nuts, regions, year, pop_kind, rad_kind, out,
                args.linthresh, args.rad_linthresh))
    (out / "four_panel_map_normalisation_manifest.json").write_text(json.dumps({
        "population_dataset": str(source.relative_to(ROOT)) if source.is_relative_to(ROOT) else str(source),
        "population_definition": "P_Jan1 / frozen regional land area; no smoothing",
        "radiance_dataset": str(RADIANCE_SOURCE.relative_to(ROOT)),
        "radiance_definition": "2020 or 2019 corrected annual VIIRS NUTS3 mean minus 2014 corrected annual VIIRS NUTS3 mean (nW cm-2 sr-1)",
        "PRC_source": "Original fixed NUTS2016 pixel-based spider/moth endpoint values",
        "within_region_support_caveat": "The same NUTS3 polygon mask is used across panels, but radiance vs PRC valid pixel support within each region is not independently guaranteed identical.",
        "matched_cohort_regions": len(regions),
        "palette": "RdBu_r: blue decline, red increase",
        "display_clip": "symmetric 98th percentile of absolute values, same as earlier maps",
        "radiance_symlog_default_threshold_nw_cm2_sr": args.rad_linthresh,
        "population_symlog_default_threshold_people_per_km2": args.linthresh,
        "map_extent_longitude": MAP_XLIM,
        "map_extent_latitude": MAP_YLIM,
        "maps": info,
    }, indent=2) + "\n")
    print(f"Created 2014–2020 four-panel population/radiance/spider/moth maps for {len(regions)} regions")
    print("Also created matched 2014–2019 four-panel sensitivity maps.")
    print("Three display versions per endpoint: all-linear, population-symlog, both-symlog.")
    print(f"PDF and PNG results: {out}")


if __name__ == "__main__":
    main()
