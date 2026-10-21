#!/usr/bin/env python3
"""Replot original cumulative NUTS2016 population / spider / moth PRC maps.

Uses the separately reconstructed January-1 population / FIXED land-area
dataset (not published DEMO_R_D3DENS). For each 2014->2019 and 2014->2020,
produce original-style three-panel maps with BOTH linear and symmetric-log
population colour normalisations. Spider/moth PRC data and calculations are
unaltered; all maps share the EXACT same 1273-region reconstructed cohort.

SymLogNorm maps signed change in people/km2, allowing declines, growth and
true zeros, unlike log10(delta). The colour bar remains in original units.
Neither temporal smoothing nor averaging of Jan-1 populations is applied.
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
    decade = 10.
    while decade < limit:
        candidates.extend([-decade, decade])
        decade *= 10
    # Do not label the exact 98th-percentile end caps: these can nearly
    # overlap the adjacent round-number decade ticks (e.g. 142 and 100).
    # The extend='both' colourbar triangles already mark clipped extremes.
    return sorted({float(z) for z in candidates if abs(z) <= limit + 1e-10})


def draw_panel(fig, ax, geom, column: str, limit: float, label: str,
               population_scale: str | None = None,
               linthresh: float = 1.) -> None:
    if population_scale == "symlog":
        norm = SymLogNorm(
            linthresh=linthresh, linscale=1.0, base=10,
            vmin=-limit, vmax=limit, clip=True)
    else:
        norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    cmap = "RdBu_r"  # Growth / improvement red, declines / losses blue.
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
    if population_scale == "symlog":
        ticks = symlog_ticks(limit, linthresh)
        bar.set_ticks(ticks)
        bar.set_ticklabels([f"{v:g}" for v in ticks])
        bar.set_label("Population density change (people/km²) | symmetric log")
    elif population_scale == "linear":
        bar.set_label("Population density change (people/km²) | linear")
    else:
        bar.set_label(label)


def make_map(nuts, data: pd.DataFrame, year: int, kind: str, out: Path,
             linthresh: float) -> dict:
    pop_col = f"change_density_jan1_fixed_land_2014_{year}"
    cols = [pop_col, f"cumulative_delta_spider_prc_{year}",
            f"cumulative_delta_moth_prc_{year}"]
    geo = nuts.merge(data[["NUTS_ID", *cols]], how="left",
                     on="NUTS_ID", validate="one_to_one")
    limits = {
        pop_col: robust_symmetric_limit(data[pop_col], 5.),
        cols[1]: robust_symmetric_limit(data[cols[1]], 1e-10),
        cols[2]: robust_symmetric_limit(data[cols[2]], 1e-10),
    }
    fig, ax = plt.subplots(1, 3, figsize=(21, 8.3),
                           layout="constrained")
    labels = [
        f"Jan 1 population density change {2014}–{year} (people/km²)",
        f"Spider PRC change {2014}–{year}",
        f"Moth PRC change {2014}–{year}",
    ]
    draw_panel(fig, ax[0], geo, pop_col, limits[pop_col],
               labels[0], population_scale=kind, linthresh=linthresh)
    for i in (1, 2):
        draw_panel(fig, ax[i], geo, cols[i], limits[cols[i]], labels[i])
    fig.suptitle(
        f"Cumulative change 2014–{year} | fixed NUTS2016 geography"
        f" | independently reconstructed Jan-1 population\n"
        f"Population colour scale: {kind}; blue=decrease; red=increase; "
        f"grey=outside complete 1,273-region cohort\n"
        f"Spider/moth PRC: unchanged locked pixel-based estimates; "
        f"all three maps have identical geographical coverage",
        fontsize=12,
    )
    basename = f"cumulative_population_spider_moth_2014_{year}_jan1_{kind}"
    fig.savefig(out / (basename + ".pdf"), bbox_inches="tight")
    fig.savefig(out / (basename + ".png"), bbox_inches="tight", dpi=200)
    plt.close(fig)
    return {"year": year, "scale": kind, "regions": int(len(data)),
            "colour_scale_max_abs_98pct": limits,
            "symlog_linear_threshold_people_per_km2": linthresh if kind=="symlog" else None}


def self_test() -> None:
    a = pd.Series([-400., -10., -1., 0., 1., 10., 400.])
    v = robust_symmetric_limit(a, 1.)
    norm = SymLogNorm(linthresh=1., vmin=-v, vmax=v, base=10)
    assert np.isclose(float(norm(0)), .5)
    assert float(norm(-10)) < .5 < float(norm(10))
    assert float(norm(-1)) < float(norm(0)) < float(norm(1))
    assert symlog_ticks(v, 1.).count(0.) == 1
    assert not any(np.isclose(abs(t), v) for t in symlog_ticks(v, 1.))
    print("Symmetric-log signed population-scale self-test: PASS")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--input", type=Path, default=INPUT)
    p.add_argument("--out", type=Path, default=OUTPUT)
    p.add_argument("--linthresh", type=float, default=1.0,
                   help="Linear part of SymLogNorm around zero (people/km2)")
    args = p.parse_args()
    if args.linthresh <= 0:
        p.error("--linthresh must be positive")
    self_test()
    if args.self_test:
        return
    source = args.input if args.input.is_absolute() else ROOT / args.input
    out = args.out if args.out.is_absolute() else ROOT / args.out
    regions = read_matched_endpoints(source)
    nuts = load_nuts(load_config())[["NUTS_ID", "geometry"]]
    if not set(regions.NUTS_ID).issubset(set(nuts.NUTS_ID)):
        raise ValueError("Reconstructed region missing from fixed NUTS2016 geometries")
    out.mkdir(parents=True, exist_ok=True)
    info = []
    for year in ENDPOINTS:
        for kind in ("linear", "symlog"):
            info.append(make_map(nuts, regions, year, kind, out, args.linthresh))
    (out / "map_normalisation_manifest.json").write_text(json.dumps({
        "population_dataset": str(source.relative_to(ROOT)) if source.is_relative_to(ROOT) else str(source),
        "population_definition": "P_Jan1 / frozen regional land area; no smoothing",
        "PRC_source": "Original fixed NUTS2016 pixel-based spider/moth endpoint values",
        "matched_cohort_regions": len(regions),
        "palette": "RdBu_r: blue decline, red increase",
        "display_clip": "symmetric 98th percentile of absolute values, same as earlier maps",
        "map_extent_longitude": MAP_XLIM,
        "map_extent_latitude": MAP_YLIM,
        "maps": info,
    }, indent=2) + "\n")
    print(f"Recreated 2014–2020 linear and symlog PRC triptychs for {len(regions)} regions")
    print("Also created matched 2014–2019 sensitivity maps.")
    print(f"PDF and PNG results: {out}")


if __name__ == "__main__":
    main()
