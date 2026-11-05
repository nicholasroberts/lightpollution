#!/usr/bin/env python3
"""Publication Figure 3A–F: locked spider/moth PRC functions, maps and histograms.

No new calibration, fitting, sensor correction, native-raster processing or
NUTS3 aggregation. Uses the existing locked spider17/moth16 exponential
transfers, FINAL snow-corrected 2014–2024 PRC NUTS3 tables, and fixed
NUTS2016 geometry. The 2014/2024 endpoints do not depend on the 2021 patch.

3A   Two independently locked calibrated-radiance -> PRC curves.
3B   Species-specific ΔPRC = f(R+ΔR)-f(R) against *starting* calibrated
     radiance, for the SAME fixed positive absolute radiance step ΔR.
     This is a transfer-function illustration, not a regional NUTS mean
     or a claim that a given pixel experienced that step.
3C   Spider 2024 minus 2014 native-pixel-derived NUTS3 PRC mean.
3D   Corresponding moth NUTS3 PRC mean change, independently scaled.
3E   Equal-NUTS3 histogram of signed absolute spider PRC change, 2014–2024.
3F   Equivalent independently calculated moth PRC change distribution.
     Both match Figure 1D histogram styling, with PRC rather than % units.

All panels export editable SVG and vector PDF via figure_style. Production
datasets and locked fit parameters are READ ONLY.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box

from figure_style import (
    CHANGE_CMAP, SERIES_COLOURS, apply_figure_style, figure_size, save_figure
)
from historical_prc_pixel_nuts3 import load_inputs, load_nuts, prc_of_radiance
from project_config import PROJECT_ROOT, load_config

ROOT = PROJECT_ROOT
TRANSFER = {
    "spider": ROOT / "config/prc_transfer.yaml",
    "moth": ROOT / "config/prc_transfer_moth.yaml",
}
ANNUAL = {
    "spider": (ROOT / "data/processed/"
               "historical_prc_2014_2024_spider_snow_corrected/"
               "prc_nuts3_annual_cumulative.csv"),
    "moth": (ROOT / "data/processed/"
             "historical_prc_2014_2024_moth_snow_corrected/"
             "prc_nuts3_annual_cumulative.csv"),
}
JOINED_ENDPOINT_QA = (ROOT / "figures/prc_radiance_spider_moth_2014_2024/"
                      "nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv")
OUTPUT = ROOT / "figures/publication/figure3_prc"

YEARS = (2014, 2024)
ALL_YEARS = tuple(range(2014, 2025))
EXPECTED_REGIONS = 1522
EXTENT_4326 = (-12.0, 34.0, 36.0, 73.0)  # Exactly match publication Fig 1.
RAD_UNIT = r"nW cm$^{-2}$ sr$^{-1}$"
LINEAR_X_TICK_COUNT = 6  # 0, 20, 40, 60, 80, 100 at the default R max.
TAXON_LABEL = {"spider": "Spider", "moth": "Moth"}
PRC_HIST_LINTHRESH_DEFAULT = .001  # Signed-log display only; PRC is dimensionless.
PRC_HIST_N_BINS = 55  # Match Figure 1D.


def resolve_path(value: Path | str) -> Path:
    p = Path(value)
    return p if p.is_absolute() else ROOT / p


def response(r_calibrated: np.ndarray, params: tuple[float, float, float]) -> np.ndarray:
    """Evaluate the locked R_SSasymp_exponential on calibrated R>=0 only.

    No 0.235 subtraction here: the x axis is ALREADY calibrated radiance.
    This must agree with historical_prc_pixel_nuts3.prc_of_radiance
    called with offset=0. The fitted parameters are not modified.
    """
    r = np.asarray(r_calibrated, dtype=np.float64)
    if not np.isfinite(r).all() or (r < 0).any():
        raise ValueError("Transfer curves require finite nonnegative calibrated radiance")
    yf, y0, k = params
    return yf + (y0 - yf) * np.exp(-k * r)


def transfer_functions(config) -> tuple[dict, dict]:
    params, details = {}, {}
    for taxon, filepath in TRANSFER.items():
        model, p, offset, from_year = load_inputs(config, filepath)
        if model.get("status") != "locked":
            raise ValueError(f"{taxon} empirical PRC function is not locked")
        if model.get("model_family") != "R_SSasymp_exponential":
            raise ValueError(f"Unexpected {taxon} transfer-family change")
        if not (0 <= p[0] < p[1] < 1 and p[2] > 0):
            raise ValueError(f"Unexpected parameters for {taxon}")
        params[taxon] = p
        details[taxon] = {
            "model_file": str(filepath.relative_to(ROOT)),
            "model_id": model["model_id"],
            "model_status": "locked",
            "model_family": model["model_family"],
            "parameters_yf_y0_k": list(map(float, p)),
            "sensor_pixel_offset": float(offset),
            "sensor_offset_from_year": int(from_year),
        }
    if np.allclose(params["spider"], params["moth"]):
        raise ValueError("Moth and spider must have independent PRC functions")
    return params, details


def read_locked_panel(filepath: Path, taxon: str) -> pd.DataFrame:
    if not filepath.is_file():
        raise FileNotFoundError(
            f"Final snow-corrected {taxon} annual NUTS3 data missing: {filepath}. "
            "This figure uses the already LOCKED published PRC panel, not an "
            "unapproved 2021 candidate or newly refitted source."
        )
    a = pd.read_csv(filepath, dtype={"NUTS_ID": str})
    required = {
        "NUTS_ID", "year", "prc_mean", "valid_area_km2",
        "cumulative_prc_change_absolute", "snow_status",
    }
    if missing := required - set(a.columns):
        raise ValueError(f"{taxon} annual source lacks {sorted(missing)}")
    if len(a) != EXPECTED_REGIONS * len(ALL_YEARS):
        raise ValueError(f"Unexpected {taxon} annual PRC row count: {len(a)}")
    if a.duplicated(["NUTS_ID", "year"]).any() or a.NUTS_ID.isna().any():
        raise ValueError(f"Duplicate or missing {taxon} region/year")
    if set(a.year.unique()) != set(ALL_YEARS):
        raise ValueError(f"Unexpected {taxon} annual years")
    if a.NUTS_ID.nunique() != EXPECTED_REGIONS:
        raise ValueError(f"Unexpected {taxon} NUTS3 cohort")
    if not a.snow_status.eq("locked").all():
        raise ValueError(f"{taxon} historical panel lacks locked-snow status")
    for key in ("prc_mean", "valid_area_km2", "cumulative_prc_change_absolute"):
        a[key] = pd.to_numeric(a[key], errors="coerce")
        if not np.isfinite(a[key].to_numpy(float)).all():
            raise ValueError(f"{taxon}: non-finite {key}")
    if (a.valid_area_km2 <= 0).any() or not a.prc_mean.between(0, 1).all():
        raise ValueError(f"{taxon}: invalid PRC or area in locked annual panel")
    return a


def endpoint_values(annual: dict[str, pd.DataFrame]) -> pd.DataFrame:
    out = None
    for taxon in ("spider", "moth"):
        a = annual[taxon].loc[
            annual[taxon].year.eq(2014),
            ["NUTS_ID", "prc_mean", "valid_area_km2"],
        ].rename(columns={
            "prc_mean": f"{taxon}_prc_2014",
            "valid_area_km2": f"{taxon}_valid_area_km2_2014",
        })
        b = annual[taxon].loc[
            annual[taxon].year.eq(2024),
            ["NUTS_ID", "prc_mean", "valid_area_km2",
             "cumulative_prc_change_absolute"],
        ].rename(columns={
            "prc_mean": f"{taxon}_prc_2024",
            "valid_area_km2": f"{taxon}_valid_area_km2_2024",
            "cumulative_prc_change_absolute": f"{taxon}_source_change",
        })
        d = a.merge(b, on="NUTS_ID", validate="one_to_one")
        if len(d) != EXPECTED_REGIONS:
            raise ValueError(f"{taxon}: missing endpoint NUTS3")
        d[f"{taxon}_prc_change_abs"] = d[f"{taxon}_prc_2024"] - d[f"{taxon}_prc_2014"]
        if not np.allclose(
            d[f"{taxon}_source_change"], d[f"{taxon}_prc_change_abs"],
            rtol=1e-5, atol=2e-7
        ):
            raise ValueError(f"{taxon}: cumulative change differs from locked endpoints")
        out = d if out is None else out.merge(
            d, on="NUTS_ID", how="outer", validate="one_to_one", indicator=True
        )
        if "_merge" in out:
            if not out._merge.eq("both").all():
                raise ValueError("Locked spider/moth regional cohorts disagree")
            out = out.drop(columns="_merge")
    for year in YEARS:
        if not np.allclose(
            out[f"spider_valid_area_km2_{year}"],
            out[f"moth_valid_area_km2_{year}"],
            rtol=1e-6, atol=.05
        ):
            raise ValueError(f"Different valid NUTS3 pixel support for taxa in {year}")
    return out.sort_values("NUTS_ID").reset_index(drop=True)


def check_joined_endpoint_source(d: pd.DataFrame, qa_path: Path) -> dict:
    """Cross-check published joined radiance/spider/moth endpoints if available."""
    if not qa_path.is_file():
        return {"checked": False, "reason": "Previously joined endpoint CSV absent"}
    q = pd.read_csv(qa_path, dtype={"NUTS_ID": str})
    required = {"NUTS_ID", "spider_prc_2014", "spider_prc_2024",
                "spider_prc_change_abs", "moth_prc_2014", "moth_prc_2024",
                "moth_prc_change_abs"}
    if missing := required - set(q):
        raise ValueError(f"Joined endpoint QA is missing {sorted(missing)}")
    if not q.NUTS_ID.is_unique or len(q) != EXPECTED_REGIONS:
        raise ValueError("Joined endpoint QA no longer has 1,522 distinct regions")
    fields = sorted(required - {"NUTS_ID"})
    j = d[["NUTS_ID", *fields]].merge(
        q[["NUTS_ID", *fields]], on="NUTS_ID",
        suffixes=("_locked", "_joined"), validate="one_to_one"
    )
    if len(j) != EXPECTED_REGIONS:
        raise ValueError("Locked and joined cohorts differ")
    discrepancies = {}
    for col in fields:
        delta = (j[f"{col}_locked"] - j[f"{col}_joined"]).abs()
        discrepancies[col] = float(delta.max())
        if not np.allclose(
            j[f"{col}_locked"], j[f"{col}_joined"], rtol=1e-5, atol=2e-6
        ):
            raise ValueError(
                f"Locked historical endpoint and joined CSV disagree: {col}, "
                f"max absolute difference {delta.max():.5g}"
            )
    return {"checked": True, "source": str(qa_path), "max_abs_error": discrepancies}


def geometry_for_maps(config, simplify_m: float):
    if not np.isfinite(simplify_m) or simplify_m < 0:
        raise ValueError("Map display simplification must be a nonnegative number")
    nuts = load_nuts(config)
    if len(nuts) != EXPECTED_REGIONS or not nuts.NUTS_ID.is_unique:
        raise ValueError("Fixed NUTS2016 geometry is not the 1,522-region cohort")
    bbox = box(*EXTENT_4326)
    clipped = gpd.clip(nuts, bbox)
    clipped = clipped.loc[~clipped.geometry.is_empty].to_crs("EPSG:3035")
    if simplify_m:
        clipped["geometry"] = clipped.geometry.simplify(
            simplify_m, preserve_topology=True
        )
    envelope = gpd.GeoSeries([bbox], crs="EPSG:4326").to_crs("EPSG:3035")
    return clipped[["NUTS_ID", "CNTR_CODE", "geometry"]], tuple(envelope.total_bounds)


def radiance_axis(ax, max_radiance: float):
    """Shared LINEAR calibrated-radiance axis for publication Figures 3A/3B.

    Plot-only modification: neither empirical PRC functions nor their
    evaluation domain and fixed-ΔR calculation are transformed or refitted.
    """
    if not np.isfinite(max_radiance) or max_radiance <= 0:
        raise ValueError("Radiance-axis maximum must be positive and finite")
    ax.set_xscale("linear")
    ticks = np.linspace(0., max_radiance, LINEAR_X_TICK_COUNT)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
    ax.set_xlim(0., max_radiance)
    ax.grid(alpha=.16)


def plot_3a(r: np.ndarray, curves: dict[str, np.ndarray],
            output_stem: Path, max_radiance: float) -> None:
    fig, ax = plt.subplots(figsize=figure_size("double", height_ratio=.68))
    for taxon in ("spider", "moth"):
        ax.plot(
            r, curves[taxon], color=SERIES_COLOURS[taxon],
            linewidth=1.9, label=f"{TAXON_LABEL[taxon]} (locked transfer)"
        )
    radiance_axis(ax, max_radiance)
    ax.set_xlabel(f"Calibrated VIIRS DNB radiance ({RAD_UNIT})")
    ax.set_ylabel("Polarization receptor contrast (dimensionless)")
    ax.set_title("3A  |  Empirical polarization receptor contrast functions", loc="left")
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    save_figure(fig, output_stem)
    plt.close(fig)


def plot_3b(r: np.ndarray, changes: dict[str, np.ndarray],
            output_stem: Path, max_radiance: float, delta_r: float) -> None:
    fig, ax = plt.subplots(figsize=figure_size("double", height_ratio=.68))
    for taxon in ("spider", "moth"):
        ax.plot(
            r, changes[taxon], color=SERIES_COLOURS[taxon],
            linewidth=1.9, label=f"{TAXON_LABEL[taxon]}"
        )
    ax.axhline(0, color="#777777", linewidth=.7)
    radiance_axis(ax, max_radiance)
    ax.set_xlabel(f"Starting calibrated VIIRS DNB radiance ({RAD_UNIT})")
    ax.set_ylabel("Change in receptor contrast, ΔPRC (dimensionless)")
    ax.set_title("3B  |  Equal absolute radiance increase at different baselines",
                 loc="left")
    ax.text(
        .97, .08, f"Fixed illustrative Δradiance = +{delta_r:g} {RAD_UNIT}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8
    )
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    save_figure(fig, output_stem)
    plt.close(fig)


def plot_prc_map(geom: gpd.GeoDataFrame, bounds: tuple, d: pd.DataFrame,
                 taxon: str, output_stem: Path) -> dict:
    metric = f"{taxon}_prc_change_abs"
    a = d[metric].to_numpy(float)
    lim = max(float(np.quantile(np.abs(a), .98)), 1e-8)
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
    g = geom.merge(d[["NUTS_ID", metric]], on="NUTS_ID",
                   how="left", validate="one_to_one")
    present = g[metric].notna()
    fig, ax = plt.subplots(figsize=figure_size("double", height_ratio=.91))
    fig.subplots_adjust(left=.01, right=.83, top=.91, bottom=.09)
    g.plot(ax=ax, color="#e3e3e3", edgecolor="none")
    if present.any():
        g.loc[present].plot(ax=ax, column=metric, cmap=CHANGE_CMAP,
                            norm=norm, edgecolor="none", linewidth=0)
    g.dissolve(by="CNTR_CODE").boundary.plot(
        ax=ax, color="#555555", linewidth=.25
    )
    xmin, ymin, xmax, ymax = bounds
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.set_axis_off()
    panel = "3C" if taxon == "spider" else "3D"
    ax.set_title(f"{panel}  |  {TAXON_LABEL[taxon]} PRC change, 2014–2024",
                 loc="left", pad=8)
    sm = ScalarMappable(cmap=CHANGE_CMAP, norm=norm)
    sm.set_array([])
    cax = fig.add_axes([.86, .19, .025, .61])
    cb = fig.colorbar(sm, cax=cax, extend="both")
    cb.set_label(f"2024 minus 2014 {taxon} PRC (dimensionless)")
    fig.text(
        .02, .035,
        "Fixed NUTS2016; 1,522 regions in source. "
        f"Separate ±98th-percentile colour limit: ±{lim:.4g}. "
        "Grey = outside study/mapped extent.",
        fontsize=7, va="bottom"
    )
    save_figure(fig, output_stem, map_figure=True)
    plt.close(fig)
    return {
        "taxon": taxon,
        "map_clip_quantile_of_abs_delta": .98,
        "symmetric_color_limit": float(lim),
        "mapped_mainland_polygons": int(present.sum()),
        "all_regions_in_source": len(d),
        "blue_sign": "negative PRC change",
        "red_sign": "positive PRC change",
        "map_display_only_simplification": True,
    }



def prc_hist_signed_transform(x: np.ndarray, linear: float) -> np.ndarray:
    """Display-only signed-log coordinate, matching Figure 1D's bin method."""
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.log1p(np.abs(x) / linear)


def prc_hist_signed_inverse(x: np.ndarray, linear: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return np.sign(x) * linear * np.expm1(np.abs(x))


def prc_histogram_data(d: pd.DataFrame, taxon: str,
                       linthresh: float) -> tuple[dict, pd.DataFrame]:
    """All fixed-cohort NUTS3 values, no clipping, density weights or % maths.

    'Absolute change' here means a difference of PRC levels in absolute
    dimensionless units, NOT abs(PRC2024 - PRC2014): preserve gain/loss signs.
    """
    if taxon not in TAXON_LABEL:
        raise ValueError(f"Unknown PRC taxon: {taxon}")
    if not np.isfinite(linthresh) or linthresh <= 0:
        raise ValueError("PRC histogram linear threshold must be positive and finite")
    metric = f"{taxon}_prc_change_abs"
    values = d[metric].to_numpy(dtype=float)
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError(f"Invalid or missing locked {taxon} endpoint PRC changes")
    if not d.NUTS_ID.is_unique or d.NUTS_ID.isna().any():
        raise ValueError("Histogram requires one value per unique NUTS3 region")
    x = prc_hist_signed_transform(values, linthresh)
    lo, hi = float(x.min()), float(x.max())
    if lo == hi:
        lo -= .5
        hi += .5
    # Same signed-log bin construction and bin count as Figure 1D.
    # Tiny padding ensures both observed extrema fall inside the edge range.
    t_edges = np.linspace(lo - 1e-9, hi + 1e-9,
                          PRC_HIST_N_BINS + 1)
    edges = prc_hist_signed_inverse(t_edges, linthresh)
    if not np.all(np.diff(edges) > 0):
        raise ValueError("PRC histogram bin edges must be strictly increasing")
    counts, _ = np.histogram(values, bins=edges)
    if int(counts.sum()) != len(values):
        raise ValueError("Histogram bins did not retain every NUTS3 value")
    summary = {
        "taxon": taxon,
        "years": "2014-2024",
        "statistic": "signed absolute PRC change = regional PRC2024 - PRC2014",
        "n_nuts3": int(len(values)),
        "n_positive_prc_change": int(np.count_nonzero(values > 0)),
        "n_negative_prc_change": int(np.count_nonzero(values < 0)),
        "n_exact_zero_prc_change": int(np.count_nonzero(values == 0)),
        "mean_prc_change_equal_nuts3_weight": float(values.mean()),
        "median_prc_change_equal_nuts3_weight": float(np.median(values)),
        "min_prc_change": float(values.min()),
        "max_prc_change": float(values.max()),
        "x_scale": "symlog",
        "symlog_linthresh_prc_units": float(linthresh),
        "histogram_binning": "55 equal steps in signed-log1p-transformed PRC change",
        "n_bins": PRC_HIST_N_BINS,
        "n_counted": int(counts.sum()),
    }
    bins = pd.DataFrame({
        "taxon": taxon,
        "bin_index": np.arange(PRC_HIST_N_BINS) + 1,
        "bin_lower_delta_prc": edges[:-1],
        "bin_upper_delta_prc": edges[1:],
        "nuts3_count": counts,
    })
    return summary, bins


def plot_prc_histogram(d: pd.DataFrame, taxon: str, output_stem: Path,
                       linthresh: float) -> tuple[dict, pd.DataFrame]:
    """Figure 1D-matched histogram: 55 bins, mean/median/zero, vector export."""
    summary, bins = prc_histogram_data(d, taxon, linthresh)
    values = d[f"{taxon}_prc_change_abs"].to_numpy(dtype=float)
    edges = np.r_[bins.bin_lower_delta_prc.to_numpy(dtype=float),
                  bins.bin_upper_delta_prc.iloc[-1]]
    fig, ax = plt.subplots(figsize=figure_size("double", height_ratio=.62))
    ax.hist(values, bins=edges, color="#638cb7",
            edgecolor="white", linewidth=.22)
    mean = summary["mean_prc_change_equal_nuts3_weight"]
    median = summary["median_prc_change_equal_nuts3_weight"]
    ax.axvline(mean, color="#d55e00", linestyle="--", linewidth=1.5,
               label=f"NUTS3 mean = {mean:+.5f}")
    ax.axvline(median, color="#222222", linestyle="-", linewidth=1.5,
               label=f"NUTS3 median = {median:+.5f}")
    ax.axvline(0, color="#555555", linestyle=":", linewidth=.9)
    ax.set_xscale("symlog", linthresh=linthresh)
    ax.set_xlim(float(edges[0]), float(edges[-1]))
    ax.set_ylabel("NUTS3 regions (count)")
    ax.set_xlabel(
        f"{TAXON_LABEL[taxon]} PRC change, 2024 minus 2014 (dimensionless)"
    )
    panel = "3E" if taxon == "spider" else "3F"
    ax.set_title(
        f"{panel}  |  Distribution of European NUTS3 "
        f"{taxon} PRC changes", loc="left"
    )
    ax.legend(loc="upper right", frameon=False)
    ax.grid(axis="y", alpha=.18)
    fig.text(
        .14, .015,
        f"n={len(values):,} of {len(d):,} fixed NUTS2016 regions; "
        "signed PRC2024 - PRC2014 (not a % or unsigned magnitude). "
        "Signed-log x axis; counts per transformed-width bin.",
        fontsize=7
    )
    fig.tight_layout(rect=(0, .06, 1, 1))
    save_figure(fig, output_stem)
    plt.close(fig)
    return summary, bins


def self_test():
    s = (.1453159020633674, .45612835620720266, .1037115607084461)
    m = (.10128334523974113, .6870948746846356, .0782312036944599)
    r = np.array([0., .2, 1., 5., 20., 100.])
    for p in (s, m):
        observed = response(r, p)
        reference = np.ma.asarray(prc_of_radiance(r, 0., *p))
        assert np.allclose(observed, reference, atol=1e-7)
        assert np.isclose(observed[0], p[1])
        assert np.all(np.diff(observed) < 0)
        step = response(r+1., p) - response(r, p)
        assert (step < 0).all()
        assert abs(step[0]) > abs(step[-1])
    assert not np.allclose(response(r, s), response(r, m))
    # Regression test for Figure 3A and 3B's shared display-only x axis.
    fig_test, ax_test = plt.subplots()
    try:
        radiance_axis(ax_test, 100.)
        assert ax_test.get_xscale() == "linear"
        assert np.allclose(ax_test.get_xticks(), [0., 20., 40., 60., 80., 100.])
        assert np.allclose(ax_test.get_xlim(), [0., 100.])
        radiance_axis(ax_test, 50.)
        assert ax_test.get_xscale() == "linear"
        assert np.allclose(ax_test.get_xticks(), [0., 10., 20., 30., 40., 50.])
    finally:
        plt.close(fig_test)
    test = pd.DataFrame({
        "NUTS_ID": ["A", "B", "C"],
        "year": [2024, 2024, 2024],
        "prc_mean": [.2, .3, .4],
        "cumulative_prc_change_absolute": [.01, -.01, .0],
    })
    assert test.NUTS_ID.is_unique
    # Separate data-only histogram tests: signed values (not magnitudes),
    # reproducible central statistics, no bin losses, exact region counts.
    trial = pd.DataFrame({
        "NUTS_ID": ["A", "B", "C", "D", "E"],
        "spider_prc_change_abs": [-.02, -.001, 0., .003, .01],
        "moth_prc_change_abs": [-.06, -.004, .001, .005, .012],
    })
    for taxon in TAXON_LABEL:
        summary, bins = prc_histogram_data(
            trial, taxon, PRC_HIST_LINTHRESH_DEFAULT
        )
        values = trial[f"{taxon}_prc_change_abs"].to_numpy()
        assert summary["n_nuts3"] == 5
        assert summary["n_counted"] == 5
        assert len(bins) == PRC_HIST_N_BINS
        assert summary["n_positive_prc_change"] == 2 if taxon == "spider" else 3
        assert np.isclose(
            summary["mean_prc_change_equal_nuts3_weight"], values.mean()
        )
        assert np.isclose(
            summary["median_prc_change_equal_nuts3_weight"], np.median(values)
        )
        edges = np.r_[bins.bin_lower_delta_prc.to_numpy(),
                      bins.bin_upper_delta_prc.iloc[-1]]
        assert np.all(np.diff(edges) > 0)
        assert int(bins.nuts3_count.sum()) == 5
        assert edges[0] < values.min() and edges[-1] > values.max()
        assert np.allclose(
            prc_hist_signed_inverse(
                prc_hist_signed_transform(values, PRC_HIST_LINTHRESH_DEFAULT),
                PRC_HIST_LINTHRESH_DEFAULT
            ), values
        )
    print("Figure 3 locked functions, endpoint and signed ΔPRC histogram tests: PASS",
          flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spider-panel", type=Path, default=ANNUAL["spider"])
    p.add_argument("--moth-panel", type=Path, default=ANNUAL["moth"])
    p.add_argument("--joined-endpoint-qa", type=Path, default=JOINED_ENDPOINT_QA)
    p.add_argument("--out", type=Path, default=OUTPUT)
    p.add_argument("--max-radiance", type=float, default=100.,
                   help="Upper calibrated R plotted in 3A/3B (nW cm-2 sr-1)")
    p.add_argument("--delta-radiance", type=float, default=1.,
                   help="Common positive absolute R increase shown in 3B")
    p.add_argument("--display-simplify-m", type=float, default=150.,
                   help="Display-only EPSG:3035 polygon simplification in metres")
    p.add_argument("--prc-hist-linthresh", type=float,
                   default=PRC_HIST_LINTHRESH_DEFAULT,
                   help="Signed-log linear region for 3E/3F (dimensionless PRC; display only)")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    self_test()
    if args.self_test:
        return
    if not np.isfinite(args.max_radiance) or args.max_radiance <= 0:
        p.error("--max-radiance must be positive and finite")
    if not np.isfinite(args.delta_radiance) or args.delta_radiance <= 0:
        p.error("--delta-radiance must be positive and finite")
    if not np.isfinite(args.prc_hist_linthresh) or args.prc_hist_linthresh <= 0:
        p.error("--prc-hist-linthresh must be positive and finite")

    config = load_config()
    parameters, details = transfer_functions(config)
    sources = {
        "spider": resolve_path(args.spider_panel),
        "moth": resolve_path(args.moth_panel),
    }
    annual = {taxon: read_locked_panel(path, taxon)
              for taxon, path in sources.items()}
    endpoints = endpoint_values(annual)
    checked = check_joined_endpoint_source(
        endpoints, resolve_path(args.joined_endpoint_qa)
    )
    nuts, bounds = geometry_for_maps(config, args.display_simplify_m)
    out = resolve_path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    apply_figure_style()

    # The curve and derivative-like difference are EVALUATIONS of the locked
    # transfer functions, not fits to recalculated NUTS3 regional mean R.
    r = np.linspace(0., args.max_radiance, 1601, dtype=np.float64)
    curves = {taxon: response(r, params) for taxon, params in parameters.items()}
    deltas = {
        taxon: response(r + args.delta_radiance, params) - curves[taxon]
        for taxon, params in parameters.items()
    }
    curve_table = pd.DataFrame({
        "starting_calibrated_radiance_nw": r,
        "spider_prc": curves["spider"],
        "moth_prc": curves["moth"],
        "spider_delta_prc_fixed_positive_deltaR": deltas["spider"],
        "moth_delta_prc_fixed_positive_deltaR": deltas["moth"],
    })
    curve_table.to_csv(out/"fig3AB_locked_function_values.csv", index=False)
    endpoints.to_csv(out/"fig3CD_locked_prc_endpoint_values.csv", index=False)

    plot_3a(r, curves, out/"fig3A_locked_spider_moth_prc_functions",
            args.max_radiance)
    plot_3b(r, deltas, out/"fig3B_equal_absolute_radiance_change",
            args.max_radiance, args.delta_radiance)
    maps = [
        plot_prc_map(nuts, bounds, endpoints, taxon,
                     out / f"fig3{panel}_{taxon}_absolute_prc_change_2014_2024")
        for taxon, panel in (("spider", "C"), ("moth", "D"))
    ]
    histograms = []
    bin_tables = []
    for taxon, panel in (("spider", "E"), ("moth", "F")):
        summary, bins = plot_prc_histogram(
            endpoints, taxon,
            out / f"fig3{panel}_{taxon}_absolute_prc_change_distribution_2014_2024",
            args.prc_hist_linthresh
        )
        histograms.append(summary)
        bin_tables.append(bins)
    pd.DataFrame(histograms).to_csv(
        out / "fig3EF_absolute_prc_change_summary.csv", index=False
    )
    pd.concat(bin_tables, ignore_index=True).to_csv(
        out / "fig3EF_absolute_prc_histogram_bins.csv", index=False
    )
    manifest = {
        "years": list(YEARS),
        "all_historical_years_in_locked_inputs": list(ALL_YEARS),
        "historical_snow_status": "locked",
        "native_pixel_transform_before_nuts_aggregation": True,
        "n_nuts3": len(endpoints),
        "locked_annual_sources": {k: str(v) for k, v in sources.items()},
        "joined_endpoint_cross_check": checked,
        "transfer_functions": details,
        "curve_domain_calibrated_radiance_nw": [0, args.max_radiance],
        "curve_3a_3b_x_axis": "linear (0 to maximum calibrated radiance), plot only",
        "curve_3a_3b_x_ticks": [
            float(v) for v in np.linspace(
                0., args.max_radiance, LINEAR_X_TICK_COUNT
            )
        ],
        "fixed_delta_radiance_nw": args.delta_radiance,
        "fixed_radiance_change_sign": "positive increase",
        "fixed_change_plot_is_illustrative": True,
        "map_extent_lonlat": list(EXTENT_4326),
        "display_simplification_epsg3035_metres": args.display_simplify_m,
        "map_color_scales": maps,
        "histogram_3e_3f": histograms,
        "histogram_3e_3f_weight": "equal weight per fixed NUTS2016 NUTS3 region",
        "histogram_3e_3f_source": "same 2014 and 2024 locked pixel-derived PRC means as maps 3C/3D",
        "histogram_3e_3f_no_percentages": True,
        "data_modified": False,
        "vector_formats": ["svg", "pdf"],
    }
    (out/"fig3_provenance.json").write_text(
        json.dumps(manifest, indent=2)+"\n", encoding="utf-8"
    )
    print("\nFIGURE 3 OUTPUTS", flush=True)
    for artifact in sorted(out.glob("fig3*")):
        print(f"  {artifact} ({artifact.stat().st_size/1048576:.2f} MiB)",
              flush=True)
    print(f"Frozen NUTS3 cohort: {len(endpoints):,}; 2021-only snow correction "
          "does not modify endpoints.", flush=True)


if __name__ == "__main__":
    main()
