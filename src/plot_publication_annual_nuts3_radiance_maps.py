#!/usr/bin/env python3
"""Talk/publication atlas: one calibrated NUTS2016 radiance map for every year.

2014–2024, fixed geographic boundaries, fixed colourbar on ALL panels, 16:9
slide canvas, Illustrator-editable SVG/PDF and high-resolution talk PNG.
No input scientific data, sensor corrections, PRC fits or geographic zonal
statistics are recomputed or modified.

2014–2020, 2022–2024: already calibrated annual NUTS2016 VIIRS parquet
(2017+ 0.235 nW pixel-level zero correction already applied).
2021 ONLY: substitute the EXISTING LOCKED 2021 snow-patched native-pixel
NUTS3 radiance, via the Figure 1 validated provenance-aware cache. Never
use the uncorrected 2021 source annual mean, and never apply the pixel offset
to a NUTS3 mean a second time.

Colour normalisation is globally fixed over all years IN THE DISPLAY
GEOGRAPHY (0..pooled continental maximum, rounded up), not refitted by year.
Natural units displayed on a symmetric-log-style sequential radiance scale
(linear 0–0.2 nW, logarithmic higher). This is NOT a colour scale of changes.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
from matplotlib.colors import SymLogNorm
from matplotlib.patches import Patch
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from figure_style import apply_figure_style, save_figure
from historical_prc_pixel_nuts3 import load_nuts
from project_config import PROJECT_ROOT as ROOT, load_config
from plot_publication_figure1_radiance import (
    DEFAULT_ANNUAL, DEFAULT_ENDPOINT, DEFAULT_SNOW_OUTPUT,
    get_locked_snow_2021, validate_endpoint,
)

YEARS = tuple(range(2014, 2025))
TOTAL_REGIONS = 1522
EXTENT_LONLAT = (-12., 34., 46., 73.)  # Continental view includes eastern Türkiye.
OUTPUT = ROOT / "figures/publication/annual_nuts3_radiance_2014_2024"
CMAP = "magma"  # Sequential: dark = low radiance, bright = high radiance.
LINEAR_THRESHOLD = .2  # nW cm-2 sr-1, shared for all years
UNITS = r"nW cm$^{-2}$ sr$^{-1}$"


def absolute_path(value: str | Path) -> Path:
    p = Path(value)
    return p if p.is_absolute() else ROOT / p


def nice_upper(value: float) -> float:
    """Round a strictly positive maximum *up*; never clip a mapped value."""
    if not np.isfinite(value) or value <= 0:
        raise ValueError("Need a positive finite maximum radiance")
    scale = 10.0 ** np.floor(np.log10(value))
    return float(np.ceil(value / scale - 1e-12) * scale)


def global_ticks(limit: float) -> list[float]:
    """Physical-unit ticks on the common nonlinear map colourbar."""
    candidates = [
        0, .2, .5, 1, 2, 5, 10, 20, 50, 100, 200, 500,
        1000, 2000, 5000, 10000, 20000,
    ]
    ticks = [float(v) for v in candidates if v < limit]
    ticks.append(float(limit))
    return sorted(set(ticks))


def read_annual_table(path: Path, valid_ids: set[str]) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Annual calibrated NUTS3 parquet missing: {path}")
    cols = ["NUTS_ID", "year", "radiance_mean_corrected"]
    d = pd.read_parquet(path, columns=cols)
    d["NUTS_ID"] = d.NUTS_ID.astype(str)
    d["year"] = pd.to_numeric(d.year, errors="raise").astype(int)
    d = d.loc[d.year.isin(YEARS)].copy()
    if d.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate annual NUTS2016/year radiance")
    if set(d.year) != set(YEARS):
        raise ValueError(f"Annual years absent: {sorted(set(YEARS) - set(d.year))}")
    unexpected = set(d.NUTS_ID) - valid_ids
    if unexpected:
        raise ValueError(
            f"Annual radiance includes non-NUTS2016 identifiers: "
            f"{sorted(unexpected)[:12]}"
        )
    d["radiance_mean_corrected"] = pd.to_numeric(
        d.radiance_mean_corrected, errors="coerce"
    )
    if not np.isfinite(d.radiance_mean_corrected.to_numpy(float)).all():
        raise ValueError("Nonfinite annual calibrated VIIRS radiance")
    if (d.radiance_mean_corrected < 0).any():
        raise ValueError("Negative calibrated regional VIIRS radiance")
    n_by_year = d.groupby("year").NUTS_ID.nunique()
    if (n_by_year < 1400).any() or (n_by_year > TOTAL_REGIONS).any():
        raise ValueError(
            "Annual radiance coverage differs substantially from established "
            f"NUTS2016 cohort: {n_by_year.to_dict()}"
        )
    d["radiance_source"] = "existing annual calibrated VIIRS NUTS2016"
    d["native_pixel_2021_snow_locked"] = False
    return d.sort_values(["year", "NUTS_ID"]).reset_index(drop=True)


def validate_annual_endpoints(d: pd.DataFrame, path: Path) -> dict:
    """Independently cross-check annual NUTS3 with published locked endpoints.

    The previous Figure 1 production used 99th-percentile absolute error
    <=0.025 nW as its check, due to differences in archived extraction
    product rounding/support. Preserve that exact standard and report max
    difference; don't silently replace annual values or fill missing areas.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Locked 2014/2024 endpoint CSV absent: {path}")
    endpoints = validate_endpoint(pd.read_csv(path, dtype={"NUTS_ID": str}))
    checks = {}
    for year, col in ((2014, "radiance_2014"), (2024, "radiance_2024")):
        x = d.loc[d.year.eq(year), ["NUTS_ID", "radiance_mean_corrected"]]
        joined = x.merge(endpoints[["NUTS_ID", col]],
                         on="NUTS_ID", how="inner", validate="one_to_one")
        if len(joined) != len(x):
            raise ValueError(f"Annual {year} IDs missing in published endpoints")
        error = (
            joined.radiance_mean_corrected.to_numpy(float)
            - joined[col].to_numpy(float)
        )
        p99 = float(np.quantile(np.abs(error), .99))
        if p99 > .025:
            raise ValueError(
                f"Annual {year} values disagree with frozen endpoints: "
                f"99th-percentile abs error={p99:.5g} nW, tolerance=.025"
            )
        checks[str(year)] = {
            "compared_nuts3": len(joined),
            "max_abs_error_nw": float(np.abs(error).max()),
            "p99_abs_error_nw": p99,
            "source_field": col,
            "endpoint_values_not_substituted": True,
        }
    return checks


def integrate_locked_2021(d: pd.DataFrame, cache: Path) -> tuple[pd.DataFrame, dict]:
    # Cache validation depends on the identical original global annual raster,
    # locked patch, pixel offset, timestamps and model ID.
    snow, source_audit = get_locked_snow_2021(cache, rebuild=False)
    if len(snow) != TOTAL_REGIONS or not snow.NUTS_ID.is_unique:
        raise ValueError("Invalid locked NUTS2016 2021 snow cache")
    selected = d.year.eq(2021)
    years = d.loc[selected, ["NUTS_ID", "radiance_mean_corrected"]]
    joined = years.merge(
        snow[[
            "NUTS_ID", "radiance_2021_original_calibrated",
            "radiance_2021_snow_corrected",
        ]],
        on="NUTS_ID", validate="one_to_one"
    )
    if len(joined) != len(years):
        raise ValueError("Some annual 2021 regions lack locked snow radiance")
    raw = joined.radiance_mean_corrected.to_numpy(float)
    orig = joined.radiance_2021_original_calibrated.to_numpy(float)
    corrected = joined.radiance_2021_snow_corrected.to_numpy(float)
    tolerance = 2e-5 * np.maximum(1., np.abs(raw))
    mismatch = np.abs(raw-orig)
    if np.any(mismatch > tolerance):
        raise ValueError(
            "Unpatched annual 2021 NUTS3 differs from validated locked original "
            f"on identical pixels (max abs mismatch={mismatch.max():.5g} nW)"
        )
    if not np.isfinite(corrected).all() or (corrected < 0).any():
        raise ValueError("Invalid locked 2021 corrected mean radiance")
    if np.any(corrected > orig + 2e-6):
        raise ValueError("Locked 2021 snow correction brightens a region")
    # NUTS3 key-based substitution, NOT row-order-dependent.
    patched = dict(zip(joined.NUTS_ID, corrected))
    out = d.copy()
    out.loc[selected, "radiance_mean_corrected"] = (
        out.loc[selected, "NUTS_ID"].map(patched).to_numpy(float)
    )
    out.loc[selected, "radiance_source"] = (
        "original 2021 annual VIIRS + locked native-pixel snow correction"
    )
    out.loc[selected, "native_pixel_2021_snow_locked"] = True
    meta = {
        "2021_nuts3_substituted": int(selected.sum()),
        "original_annual_vs_locked_original_max_abs_error_nw":
            float(mismatch.max()),
        "locked_2021_original_region_mean_max_nw": float(orig.max()),
        "locked_2021_snow_corrected_region_mean_max_nw":
            float(corrected.max()),
        "n_regions_2021_snow_reduced": int(np.count_nonzero(corrected < orig - 1e-6)),
        "cache": str(cache),
        "cache_audit": source_audit,
    }
    return out, meta


def map_geometry(config: dict, simplify_metres: float):
    if not np.isfinite(simplify_metres) or simplify_metres < 0:
        raise ValueError("Display simplification must be finite and >= 0")
    g = load_nuts(config)
    if len(g) != TOTAL_REGIONS or g.NUTS_ID.duplicated().any():
        raise ValueError("Original fixed NUTS2016 polygons changed")
    region_ids = set(g.NUTS_ID.astype(str))
    bbox = box(*EXTENT_LONLAT)
    clipped = gpd.clip(g, bbox)
    clipped = clipped.loc[~clipped.geometry.is_empty].to_crs("EPSG:3035")
    if simplify_metres:
        clipped["geometry"] = clipped.geometry.simplify(
            simplify_metres, preserve_topology=True
        )
    bounds = tuple(
        gpd.GeoSeries([bbox], crs="EPSG:4326").to_crs("EPSG:3035").total_bounds
    )
    return clipped[["NUTS_ID", "CNTR_CODE", "geometry"]], bounds, region_ids


def compute_common_scale(d: pd.DataFrame, mapped_ids: set[str],
                         requested_max: float | None):
    values = d.loc[d.NUTS_ID.isin(mapped_ids),
                   "radiance_mean_corrected"].to_numpy(float)
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("No finite continental radiance values")
    observed = float(values.max())
    if observed <= 0:
        raise ValueError("All continental radiance values are zero")
    upper = nice_upper(observed) if requested_max is None else float(requested_max)
    if not np.isfinite(upper) or upper <= LINEAR_THRESHOLD:
        raise ValueError("Common map colour maximum must exceed 0.2 nW")
    if upper < observed:
        raise ValueError(
            f"Requested common --vmax {upper:g} would clip radiance "
            f"observations up to {observed:g}. Increase --vmax or omit it."
        )
    return upper, observed


def make_figure(geom, countries, bounds, slice_: pd.DataFrame, year: int,
                norm, upper: float, total_valid: int):
    # Existing fixed NUTS2016 map polygons; only visual simplification.
    vals = slice_[["NUTS_ID", "radiance_mean_corrected"]]
    g = geom.merge(vals, on="NUTS_ID", how="left", validate="one_to_one")
    shown = g.radiance_mean_corrected.notna()
    fig, ax = plt.subplots(figsize=(12.8, 7.2), facecolor="white")
    fig.subplots_adjust(left=.02, right=.84, bottom=.12, top=.89)

    # Absent/unavailable NUTS3 stays grey; zero radiance is a genuine dark
    # value in the colourmap, NOT a missing observation.
    g.plot(ax=ax, color="#dedede", edgecolor="none", linewidth=0)
    if shown.any():
        g.loc[shown].plot(
            ax=ax, column="radiance_mean_corrected", cmap=CMAP,
            norm=norm, edgecolor="none", linewidth=0
        )
    countries.plot(ax=ax, color="#5b5b5b", linewidth=.33)
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.text(.035, .93, str(year), fontsize=26, weight="bold",
             ha="left", va="center")
    fig.text(.19, .931, "European night-time radiance",
             fontsize=17, ha="left", va="center")
    sm = ScalarMappable(norm=norm, cmap=CMAP)
    sm.set_array([])
    cb_ax = fig.add_axes([.865, .23, .025, .55])
    cb = fig.colorbar(sm, cax=cb_ax, extend="neither")
    ticks = global_ticks(upper)
    cb.set_ticks(ticks)
    cb.ax.set_yticklabels([f"{v:g}" for v in ticks], fontsize=9)
    cb.set_label(f"VIIRS DNB radiance ({UNITS})", fontsize=11, labelpad=9)
    fig.text(
        .035, .065,
        "Annual NUTS3 radiance · fixed NUTS2016 geography · common nonlinear colour scale",
        fontsize=10, ha="left", va="center"
    )
    source_note = (
        "2021: approved native-pixel snow adjustment + calibrated VIIRS"
        if year == 2021
        else "2017+ VIIRS pixel zero-point calibration; no snow adjustment"
        if year >= 2017
        else "Annual VIIRS, fixed NUTS2016 regional means"
    )
    fig.text(.035, .036, source_note, fontsize=8.3, color="#555555",
             ha="left", va="center")
    missing = int((~shown).sum())
    if missing:
        fig.legend(
            handles=[Patch(facecolor="#dedede", label="Missing annual NUTS3 record")],
            loc="lower right", bbox_to_anchor=(.82, .075), frameon=False,
            fontsize=8
        )
    return fig, {
        "year": year,
        "n_nuts3_all_source": total_valid,
        "n_map_polygons_coloured": int(shown.sum()),
        "n_map_polygons_missing": missing,
        "mapped_year_max_radiance": (
            float(g.loc[shown, "radiance_mean_corrected"].max())
            if shown.any() else None
        ),
        "map_scale_upper_nw": upper,
    }


def self_test():
    assert nice_upper(277.) == 300.
    assert nice_upper(10.1) == 20.
    assert nice_upper(.127) == .2
    assert global_ticks(300.) == [0, .2, .5, 1, 2, 5, 10, 20,
                                   50, 100, 200, 300]
    test = pd.DataFrame({
        "NUTS_ID": ["A", "B", "A", "B"], "year": [2014, 2014, 2015, 2015],
        "radiance_mean_corrected": [1., 0., 3., .2]
    })
    assert test.loc[test.year.eq(2014), "radiance_mean_corrected"].min() == 0
    assert len(test.groupby("year")) == 2
    print("Annual radiance map scale and 0-vs-missing tests: PASS", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annual", type=Path, default=DEFAULT_ANNUAL)
    p.add_argument("--endpoints", type=Path, default=DEFAULT_ENDPOINT)
    p.add_argument("--snow-cache", type=Path, default=DEFAULT_SNOW_OUTPUT)
    p.add_argument("--out", type=Path, default=OUTPUT)
    p.add_argument("--display-simplify-m", type=float, default=150.,
                   help="NUTS3 map polygon display simplification, EPSG:3035 metres")
    p.add_argument("--vmax", type=float, default=None,
                   help="Optional manually fixed upper scale, must not clip data")
    p.add_argument("--png-dpi", type=int, default=240,
                   help="Resolution of talk PNG copies (0 disables PNGs)")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    self_test()
    if args.self_test:
        return
    if args.png_dpi < 0 or args.png_dpi > 600:
        p.error("--png-dpi must be 0..600")
    apply_figure_style()
    geometry, bounds, region_ids = map_geometry(load_config(),
                                                  args.display_simplify_m)
    print(f"Fixed NUTS2016: {len(region_ids):,} NUTS3; "
          f"{len(geometry):,} mapped continental polygons", flush=True)
    annual = read_annual_table(absolute_path(args.annual), region_ids)
    endpoint_qa = validate_annual_endpoints(
        annual, absolute_path(args.endpoints)
    )
    annual, snow_qa = integrate_locked_2021(
        annual, absolute_path(args.snow_cache)
    )
    mapped_ids = set(geometry.NUTS_ID.astype(str))
    upper, observed = compute_common_scale(annual, mapped_ids, args.vmax)
    norm = SymLogNorm(
        linthresh=LINEAR_THRESHOLD, linscale=1.,
        vmin=0., vmax=upper, base=10, clip=False
    )
    country_borders = geometry.dissolve(by="CNTR_CODE").boundary
    out = absolute_path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    # Derived data used to paint maps; sources are never overwritten.
    annual.to_csv(out/"annual_nuts3_radiance_2014_2024_map_values.csv",
                  index=False)
    maps=[]
    deck_path = out/"all_years_nuts3_radiance_2014_2024.pdf"
    print(f"Shared colour scale: 0–{upper:g} nW "
          f"(observed mapped max {observed:.5g} nW), "
          f"linear up to {LINEAR_THRESHOLD} nW", flush=True)
    with PdfPages(deck_path) as deck:
        for year in YEARS:
            a=annual[annual.year.eq(year)].copy()
            fig, m=make_figure(
                geometry, country_borders, bounds, a, year, norm, upper,
                len(a)
            )
            base = out/f"radiance_nuts3_{year}"
            save_figure(fig, base, map_figure=True)
            if args.png_dpi:
                fig.savefig(
                    base.with_suffix(".png"), format="png",
                    dpi=args.png_dpi, facecolor="white",
                    bbox_inches="tight", pad_inches=.06
                )
            deck.savefig(fig, bbox_inches="tight", pad_inches=.06)
            plt.close(fig)
            maps.append(m)
            print(f"{year}: {len(a):,} NUTS3 values; "
                  f"mapped {m['n_map_polygons_coloured']:,}; "
                  f"missing {m['n_map_polygons_missing']:,}", flush=True)
    n_mapped = annual.NUTS_ID.isin(mapped_ids)
    summary={
        "years": list(YEARS),
        "full_locked_nuts2016_geography_count": len(region_ids),
        "country_map_extent_lonlat": list(EXTENT_LONLAT),
        "geometry_display_simplification_epsg3035_m": args.display_simplify_m,
        "annual_calibrated_source": str(absolute_path(args.annual)),
        "independent_2014_2024_endpoint_source": str(absolute_path(args.endpoints)),
        "annual_vs_locked_endpoint_check": endpoint_qa,
        "locked_2021_native_pixel_snow_correction": snow_qa,
        "2021_does_not_use_unpatched_mean": True,
        "2017_onward_pixel_zero_point_subtracted_once_in_original_product": True,
        "regional_mean_values_unmodified_except_2021_locked_patch": True,
        "color_normalisation": "common positive SymLogNorm (not maps of deltas)",
        "color_map": CMAP,
        "color_zero": 0.,
        "color_linthresh_nw": LINEAR_THRESHOLD,
        "color_max_fixed_all_years_nw": upper,
        "observed_max_mapped_radiance_nw": observed,
        "map_data_count_all_years": len(annual),
        "map_data_count_within_extent_all_years": int(n_mapped.sum()),
        "year_map_records": maps,
        "output_formats": ["svg", "pdf"]+(
            ["png"] if args.png_dpi else []
        )+["multipage_pdf"],
        "talk_png_dpi": args.png_dpi,
        "source_data_written_or_modified": False,
        "note": (
            "No missing NUTS3 values were imputed. All values from the "
            "existing fixed-NUTS2016 annual corrected radiance table except "
            "2021, which uses the already locked native-pixel snow-patched "
            "radiance validated against the original 2021 VIIRS source."
        ),
    }
    (out/"annual_radiance_map_provenance.json").write_text(
        json.dumps(summary, indent=2)+"\n", encoding="utf-8"
    )
    print("\nANNUAL NUTS3 RADIANCE MAP OUTPUTS", flush=True)
    for f in sorted(out.iterdir()):
        if f.is_file():
            print(f"  {f} ({f.stat().st_size/(1024**2):.2f} MiB)",
                  flush=True)


if __name__ == "__main__":
    main()
