#!/usr/bin/env python3
"""Illustrator-ready Figure 1 A, B, C + NUTS3 percent-change histogram.

READ ONLY: consumes the *existing* final 1,522-region calibrated
radiance+spider+moth endpoint CSV and the previously extracted 2014–2024
fixed-NUTS2016 annual radiance parquet. Does not extract rasters, recalculate
sensor corrections, refit models, or alter any scientific source tables.

A. Absolute corrected VIIRS radiance change, 2014–2024 (NUTS3).
B. Percentage corrected VIIRS radiance change, 2014–2024, recalculated
   from the LOCKED endpoints for EVERY region with positive baseline,
   including those previously excluded at <0.2 nW. Regions with exactly
   zero baseline are explicitly identified as undefined percentage, not grey.
C. European area-weighted annual radiance, now including 2021 via the
   EXISTING LOCKED native VIIRS snow patch and 0.235 pixel offset once.
D. Full mathematically defined NUTS3 percentage-change histogram with
   equal-NUTS3 mean and median; zero-baseline exclusions reported explicitly.

SVG (editable text) and PDF (embedded TrueType) exported with figure_style.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
from matplotlib.colors import SymLogNorm
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box

from project_config import PROJECT_ROOT, load_config
from historical_prc_pixel_nuts3 import load_nuts
from figure1_2021_locked_snow_radiance import extract_locked_2021, _lock_info
from extract_nuts3_radiance import discover_rasters
from figure_style import (
    apply_figure_style, figure_size, save_figure, CHANGE_CMAP, SERIES_COLOURS
)

DEFAULT_ENDPOINT = (
    "figures/prc_radiance_spider_moth_2014_2024/"
    "nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv"
)
DEFAULT_ANNUAL = "data/processed/fixed_nuts2016_viirs_2014_2024.parquet"
DEFAULT_OUT = "figures/publication/figure1_radiance"
YEARS = tuple(range(2014, 2025))
EPS_ZERO = 0.0  # Only a literal zero starting radiance makes % mathematically undefined.
DEFAULT_SNOW_OUTPUT = "data/processed/publication_figures/fig1_2021_locked_snow_radiance.csv"
RADIANCE_UNIT = r"nW cm$^{-2}$ sr$^{-1}$"
EUROPE_EXTENT_LONLAT = (-12.0, 34.0, 36.0, 73.0)


def as_repo_path(value: str | Path) -> Path:
    p = Path(value)
    return p if p.is_absolute() else PROJECT_ROOT / p


def validate_endpoint(data: pd.DataFrame) -> pd.DataFrame:
    """Reuse the published radiance endpoint columns, no new correction."""
    required = {
        "NUTS_ID", "radiance_2014", "radiance_2024",
        "radiance_change_abs", "radiance_change_pct",
        "valid_area_km2_2014",
    }
    missing = required - set(data)
    if missing:
        raise ValueError(f"Final endpoint CSV is missing columns: {sorted(missing)}")
    if data.NUTS_ID.isna().any() or data.NUTS_ID.duplicated().any():
        raise ValueError("Missing or duplicate NUTS3 IDs in endpoint CSV")
    if len(data) != 1522:
        raise ValueError(f"Expected frozen 1,522-region endpoint cohort, got {len(data)}")
    d = data.copy()
    for col in required - {"NUTS_ID"}:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    for col in ("radiance_2014", "radiance_2024", "radiance_change_abs",
                "valid_area_km2_2014"):
        if not np.isfinite(d[col]).all():
            raise ValueError(f"Non-finite endpoint values for {col}")
    if (d.radiance_2014 < 0).any() or (d.radiance_2024 < 0).any():
        raise ValueError("Negative calibrated radiance in endpoint source")
    if (d.valid_area_km2_2014 <= 0).any():
        raise ValueError("Non-positive area weights")
    expected_abs = d.radiance_2024 - d.radiance_2014
    if not np.allclose(expected_abs, d.radiance_change_abs, atol=1e-6, rtol=1e-7):
        raise ValueError("Endpoint absolute differences disagree with source columns")
    # The archival % column intentionally suppressed R2014 < 0.2. That
    # display policy was NOT missing endpoint data, so reconstruct % from
    # unmodified 2014/2024 values for every strictly positive baseline.
    positive = d.radiance_2014 > EPS_ZERO
    d["percent_baseline_zero"] = ~positive
    d["radiance_change_pct_all"] = np.nan
    d.loc[positive, "radiance_change_pct_all"] = (
        100.0 * expected_abs[positive] / d.loc[positive, "radiance_2014"]
    )
    if not np.isfinite(d.loc[positive, "radiance_change_pct_all"]).all():
        raise ValueError("Invalid recalculated percentage for positive starting radiance")
    # Independently check the already-reported values wherever present.
    valid_archived = d.radiance_change_pct.notna()
    if (valid_archived & ~positive).any():
        raise ValueError("Original source reports a percentage at zero radiance")
    if not np.allclose(
        d.loc[valid_archived, "radiance_change_pct_all"],
        d.loc[valid_archived, "radiance_change_pct"],
        atol=1e-5, rtol=1e-6
    ):
        raise ValueError("Recalculated percentages disagree with original published values")
    return d


def endpoint_summary(d: pd.DataFrame) -> dict:
    valid = d.radiance_change_pct_all.dropna().to_numpy(dtype=float)
    weights = d.valid_area_km2_2014.to_numpy(dtype=float)
    # Distinct estimands: equal-region mean/median vs covered-area-weighted
    # percentage change of the *European mean radiance*.
    r0 = float(np.average(d.radiance_2014.to_numpy(float), weights=weights))
    r1 = float(np.average(d.radiance_2024.to_numpy(float), weights=weights))
    return {
        "n_nuts3_all_endpoints": len(d),
        "n_nuts3_valid_percent_change": int(len(valid)),
        "n_undefined_percent_exact_zero_baseline": int(len(d) - len(valid)),
        "n_reinstated_positive_baseline_below_0p2": int(
            ((d.radiance_2014 > 0) & (d.radiance_2014 < .2)).sum()
        ),
        "percentage_definition": "100 * (R2024-R2014) / R2014, only if R2014 > 0",
        "mean_nuts3_percent_change_equal_region_weight": float(valid.mean()),
        "median_nuts3_percent_change_equal_region_weight": float(np.median(valid)),
        "n_regions_radiance_increase": int((d.radiance_change_abs > 0).sum()),
        "n_regions_radiance_decrease": int((d.radiance_change_abs < 0).sum()),
        "2014_europe_mean_radiance_area_weighted_nw_cm2_sr": r0,
        "2024_europe_mean_radiance_area_weighted_nw_cm2_sr": r1,
        "europe_change_of_area_weighted_radiance_percent": 100.0 * (r1 / r0 - 1),
        "histogram_x_transform": "symmetric log; 10 percentage point linear region",
    }


def display_geometry(config: dict, simplify_m: float) -> tuple[gpd.GeoDataFrame, tuple]:
    """Simplification is VISUAL ONLY, never used for calculations/weights."""
    if simplify_m < 0:
        raise ValueError("Simplification tolerance must be >= 0 metres")
    nuts = load_nuts(config)
    if nuts.NUTS_ID.duplicated().any() or len(nuts) != 1522:
        raise ValueError("Unexpected fixed NUTS2016 geometry cohort")
    bbox = box(*EUROPE_EXTENT_LONLAT)
    clip = gpd.clip(nuts, bbox)
    clip = clip.loc[~clip.geometry.is_empty].to_crs("EPSG:3035")
    if simplify_m > 0:
        clip["geometry"] = clip.geometry.simplify(
            tolerance=simplify_m, preserve_topology=True
        )
    envelope = gpd.GeoSeries([bbox], crs="EPSG:4326").to_crs("EPSG:3035")
    return clip[["NUTS_ID", "CNTR_CODE", "geometry"]], tuple(envelope.total_bounds)


def publication_map(d: pd.DataFrame, geom: gpd.GeoDataFrame,
                    bounds: tuple, metric: str, name: str,
                    label: str, output: Path, linthresh: float) -> None:
    vals = d[metric].dropna().to_numpy(dtype=float)
    lim = float(np.quantile(np.abs(vals), 0.98))
    lim = max(lim, linthresh * 1.05)
    norm = SymLogNorm(
        linthresh=linthresh, linscale=1.0,
        vmin=-lim, vmax=lim, base=10, clip=False
    )
    g = geom.merge(
        d[["NUTS_ID", metric, "percent_baseline_zero"]], on="NUTS_ID",
        how="left", validate="one_to_one"
    )
    fig, ax = plt.subplots(figsize=figure_size("double", height_ratio=0.91))
    fig.subplots_adjust(left=0.01, right=0.83, top=0.91, bottom=0.09)
    g.plot(ax=ax, color="#e3e3e3", edgecolor="none")
    valid = g[metric].notna()
    g.loc[valid].plot(
        ax=ax, column=metric, cmap=CHANGE_CMAP, norm=norm,
        edgecolor="none", linewidth=0
    )
    # Zero starting radiance has undefined %, NOT missing geographic data.
    zero_baseline = g["percent_baseline_zero"] if metric == "radiance_change_pct_all" else pd.Series(False, index=g.index)
    if zero_baseline.any():
        g.loc[zero_baseline].plot(
            ax=ax, color="#8C75A2", edgecolor="none", linewidth=0
        )
        ax.legend(
            handles=[Patch(facecolor="#8C75A2",
                           label="2014 radiance = 0; percentage undefined")],
            loc="lower left", frameon=False, fontsize=7
        )
    borders = g.dissolve(by="CNTR_CODE").boundary
    borders.plot(ax=ax, color="#555555", linewidth=0.25)
    x0, y0, x1, y1 = bounds
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title(name, loc="left", pad=8)
    sm = ScalarMappable(norm=norm, cmap=CHANGE_CMAP)
    sm.set_array([])
    cax = fig.add_axes([0.86, 0.19, 0.025, 0.61])
    cb = fig.colorbar(sm, cax=cax, extend="both")
    cb.set_label(label)
    fig.text(
        0.02, 0.035,
        f"Fixed NUTS2016; n={len(d):,} regions. "
        f"Purple = exactly zero 2014 radiance; percentage undefined. "
        f"Display saturation beyond ±98th percentile (±{lim:.3g}).",
        fontsize=6.8, va="bottom"
    )
    save_figure(fig, output, map_figure=True)
    plt.close(fig)
    print(
        f"{output.name}: mapped {int(valid.sum())} numeric polygons; "
        f"{int(zero_baseline.sum())} zero-baseline polygons shown separately; "
        f"colour limits ±{lim:.4g}", flush=True
    )


def get_locked_snow_2021(cache_csv: Path, *, rebuild: bool = False):
    """Read a verified cached 2021 zonal result or extract the locked patch.

    The cached table is reused ONLY if the original VIIRS and patch files,
    their modification times and byte sizes, model ID and pixel offset agree
    with the current locked sources.
    """
    cfg = load_config()
    lock, offset = _lock_info(cfg)
    patch = PROJECT_ROOT / lock["patch"]
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
    audit_path = cache_csv.with_suffix(".json")
    if not rebuild and cache_csv.is_file() and audit_path.is_file():
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if all(audit.get(k) == v for k, v in expected.items()):
            snow = pd.read_csv(cache_csv, dtype={"NUTS_ID": str})
            print("Reusing validated locked 2021 radiance cache:", cache_csv, flush=True)
        else:
            print("Locked 2021 sources changed; regenerating cached extraction", flush=True)
            rebuild = True
    else:
        rebuild = True
    if rebuild:
        print("Extracting locked 2021 snow-patched VIIRS radiance on NUTS2016...",
              flush=True)
        snow, audit = extract_locked_2021()
        cache_csv.parent.mkdir(parents=True, exist_ok=True)
        snow.to_csv(cache_csv, index=False)
        audit_path.write_text(json.dumps(audit, indent=2) + "\n",
                              encoding="utf-8")
    necessary = {
        "NUTS_ID", "radiance_2021_original_calibrated",
        "radiance_2021_snow_corrected",
        "valid_area_km2_2021_original", "valid_area_km2_2021_snow"
    }
    if len(snow) != 1522 or snow.NUTS_ID.duplicated().any() or not necessary.issubset(snow):
        raise ValueError("2021 locked snow radiance cache has invalid cohort/columns")
    for c in necessary - {"NUTS_ID"}:
        if not np.isfinite(pd.to_numeric(snow[c], errors="coerce")).all():
            raise ValueError(f"Non-finite cached 2021 snow input: {c}")
    if (snow.radiance_2021_snow_corrected > snow.radiance_2021_original_calibrated+2e-6).any():
        raise ValueError("Locked snow correction unexpectedly brightens a region")
    if not np.allclose(
        snow.valid_area_km2_2021_original,
        snow.valid_area_km2_2021_snow,
        rtol=1e-7, atol=.05
    ):
        raise ValueError("2021 locked snow extraction changed pixel-area support")
    return snow, audit


def load_annual(path: Path, d: pd.DataFrame,
                snow_2021: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Use previously calibrated annual NUTS means, replacing 2021 with
    the extracted radiance of the EXISTING approved native-pixel snow patch.
    Every other year remains untouched; 2021 original values are checked.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Annual NUTS2016 radiance table missing: {path}."
        )
    cols = ["NUTS_ID", "year", "radiance_mean_corrected"]
    raw = pd.read_parquet(path, columns=cols)
    if raw.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate NUTS3/year in annual radiance source")
    raw["NUTS_ID"] = raw.NUTS_ID.astype(str)
    raw["year"] = pd.to_numeric(raw.year, errors="raise").astype(int)
    raw = raw.loc[raw.year.isin(YEARS)].copy()
    if set(raw.year) != set(YEARS):
        raise ValueError(
            f"Missing one or more annual years: {sorted(set(YEARS)-set(raw.year))}"
        )
    common = set.intersection(*[
        set(raw.loc[raw.year.eq(y), "NUTS_ID"]) for y in YEARS
    ])
    common &= set(d.NUTS_ID)
    if len(common) < 1400:
        raise ValueError(f"Only {len(common)} annual regions have complete histories")
    weights = d.loc[d.NUTS_ID.isin(common), [
        "NUTS_ID", "valid_area_km2_2014", "radiance_2014", "radiance_2024"
    ]]
    a = raw.merge(weights, on="NUTS_ID", how="inner", validate="many_to_one")
    a["radiance_mean_corrected"] = pd.to_numeric(
        a.radiance_mean_corrected, errors="coerce"
    )
    if len(a) != len(common) * len(YEARS):
        raise ValueError("Annual cohort not fully balanced")
    if not np.isfinite(a.radiance_mean_corrected).all() or (
        a.radiance_mean_corrected < 0
    ).any():
        raise ValueError("Invalid annual NUTS3 calibrated radiance")
    # Endpoint QA: compare archived annual table against independently
    # extracted final 2014 and 2024 source (no recalibration).
    for yr, col in ((2014, "radiance_2014"), (2024, "radiance_2024")):
        sl = a.loc[a.year.eq(yr)]
        diff = (sl.radiance_mean_corrected - sl[col]).abs()
        if diff.quantile(.99) > .025:
            raise ValueError(
                f"Annual vs published endpoint mismatch in {yr}: "
                f"99th % |difference|={diff.quantile(.99):.5g} nW; inspect sources."
            )
    # 2021 must use locked native patch and the identical once-only offset.
    z = snow_2021[[
        "NUTS_ID", "radiance_2021_original_calibrated",
        "radiance_2021_snow_corrected"
    ]]
    a = a.merge(z, how="left", on="NUTS_ID", validate="many_to_one")
    idx = a.year.eq(2021)
    if not np.isfinite(a.loc[idx, "radiance_2021_snow_corrected"]).all():
        raise ValueError("Locked 2021 snow radiance is missing for annual cohort")
    # The old annual NUTS2016 extraction and the overlay's unmodified
    # original radiance must agree. If they do not, do not splice curves.
    original_gap = (
        a.loc[idx, "radiance_mean_corrected"]
        - a.loc[idx, "radiance_2021_original_calibrated"]
    ).abs()
    rel_tolerance = 2e-5 * np.maximum(
        1.0, a.loc[idx, "radiance_mean_corrected"].abs()
    )
    if (original_gap > rel_tolerance).any():
        raise ValueError(
            f"Old 2021 table vs identical-grid corrected extraction mismatch "
            f"(max={original_gap.max():.6g} nW); do not silently splice histories."
        )
    a.loc[idx, "radiance_mean_corrected"] = (
        a.loc[idx, "radiance_2021_snow_corrected"].to_numpy()
    )
    a["radiance_times_area"] = (
        a.radiance_mean_corrected * a.valid_area_km2_2014
    )
    out = a.groupby("year", as_index=False).agg(
        weighted_sum=("radiance_times_area", "sum"),
        area_sum=("valid_area_km2_2014", "sum")
    )
    out["area_weighted_radiance_nw"] = out.weighted_sum / out.area_sum
    out = out[["year", "area_weighted_radiance_nw"]]
    out["snow_correction_applied"] = out.year.eq(2021)
    # Export original 2021 Europe mean alongside the corrected one for QA.
    pre = a.loc[idx].copy()
    original_europe = float(np.average(
        pre.radiance_2021_original_calibrated,
        weights=pre.valid_area_km2_2014
    ))
    new_europe = float(out.loc[out.year.eq(2021),
                               "area_weighted_radiance_nw"].iloc[0])
    if new_europe > original_europe + 2e-6:
        raise ValueError("2021 European mean unexpectedly increased after snow correction")
    out["2021_original_area_weighted_radiance_nw"] = np.where(
        out.year.eq(2021), original_europe, np.nan
    )
    meta = {
        "annual_nuts3_balanced_cohort": len(common),
        "endpoint_nuts3_cohort": len(d),
        "area_weight": "fixed valid_area_km2_2014 of matched endpoint NUTS3",
        "2021_annual_radiance_status": "LOCKED native-pixel snow overlay then −0.235 pixel correction",
        "2021_uncorrected_area_weighted_radiance_nw": original_europe,
        "2021_snow_corrected_area_weighted_radiance_nw": new_europe,
        "2021_snow_effect_area_weighted_nw": new_europe-original_europe,
        "2021_annual_unmodified_source_match_max_error_nw": float(original_gap.max()),
    }
    return out, meta


def publication_annual(annual: pd.DataFrame, output: Path, n: int) -> None:
    fig, ax = plt.subplots(figsize=figure_size("double", height_ratio=0.64))
    ax.plot(
        annual.year, annual.area_weighted_radiance_nw,
        color=SERIES_COLOURS["radiance"],
        marker="o", markersize=4, linewidth=1.6
    )
    ax.set_xticks(YEARS)
    ax.set_xlim(2013.7, 2024.3)
    ax.set_ylabel(f"Area-weighted VIIRS DNB radiance ({RADIANCE_UNIT})")
    ax.set_xlabel("Year")
    ax.set_title("European annual corrected VIIRS DNB radiance, 2014–2024", loc="left")
    ax.grid(axis="y", color="#e5e5e5", linewidth=0.65)
    # 2021 value explicitly identifies the locked correction.
    y2021 = float(annual.loc[annual.year.eq(2021),
                             "area_weighted_radiance_nw"].iloc[0])
    ax.annotate(
        "2021 snow-adjusted",
        xy=(2021, y2021), xytext=(2021.25, y2021),
        fontsize=8, ha="left", va="bottom",
        arrowprops={"arrowstyle": "-", "linewidth": .7, "color": "#555555"}
    )
    fig.text(
        .14, .02,
        f"Area-weighted mean of {n:,} fixed NUTS2016 regions with complete "
        "annual records; 2021 uses locked native-pixel snow correction.",
        fontsize=7
    )
    fig.tight_layout(rect=(0, .055, 1, 1))
    save_figure(fig, output)
    plt.close(fig)


def signed_transform(x: np.ndarray, linear: float = 10.0) -> np.ndarray:
    return np.sign(x) * np.log1p(np.abs(x) / linear)


def signed_inverse(x: np.ndarray, linear: float = 10.0) -> np.ndarray:
    return np.sign(x) * linear * np.expm1(np.abs(x))


def publication_hist(d: pd.DataFrame, summary: dict, output: Path) -> None:
    vals = d.radiance_change_pct_all.dropna().to_numpy(dtype=float)
    t = signed_transform(vals)
    # Log-spaced tails ensure that every valid region is represented despite
    # the long positive percentage-change tail. Counts, not densities.
    lo, hi = float(np.min(t)), float(np.max(t))
    edges = signed_inverse(np.linspace(lo - 1e-9, hi + 1e-9, 56))
    fig, ax = plt.subplots(figsize=figure_size("double", height_ratio=0.62))
    ax.hist(vals, bins=edges, color="#638cb7", edgecolor="white", linewidth=.22)
    mean = summary["mean_nuts3_percent_change_equal_region_weight"]
    median = summary["median_nuts3_percent_change_equal_region_weight"]
    ax.axvline(mean, color="#d55e00", linestyle="--", linewidth=1.5,
               label=f"NUTS3 mean = {mean:+.1f}%")
    ax.axvline(median, color="#222222", linestyle="-", linewidth=1.5,
               label=f"NUTS3 median = {median:+.1f}%")
    ax.axvline(0, color="#555555", linestyle=":", linewidth=.9)
    ax.set_xscale("symlog", linthresh=10)
    ax.set_xlim(float(np.min(edges)), float(np.max(edges)))
    ax.set_ylabel("NUTS3 regions (count)")
    ax.set_xlabel("Change in regional mean VIIRS radiance, 2014–2024 (%)")
    ax.set_title("Distribution of European NUTS3 percentage radiance changes", loc="left")
    ax.legend(loc="upper right", frameon=False)
    ax.grid(axis="y", alpha=.18)
    fig.text(
        .14, .015,
        f"n={len(vals):,} of {len(d):,} NUTS3; "
        f"{len(d)-len(vals)} have exactly zero baseline, so % undefined. "
        "Signed-log x axis; counts per transformed-width bin.",
        fontsize=7
    )
    fig.tight_layout(rect=(0, .06, 1, 1))
    save_figure(fig, output)
    plt.close(fig)


def self_test() -> None:
    trial = pd.DataFrame({
        "NUTS_ID": ["X", "Y", "Z"],
        "radiance_2014": [1.0, .1, 3.0],
        "radiance_2024": [.7, .2, 3.6],
        "radiance_change_abs": [-.3, .1, .6],
        "radiance_change_pct": [-30.0, np.nan, 20.0],
        "valid_area_km2_2014": [1., 2., 3.],
    })
    # The original threshold is a source-figure display mask, not % math.
    eligible = trial.radiance_2014.gt(0)
    assert eligible.all()
    expected = 100 * (
        trial.loc[eligible, "radiance_2024"] /
        trial.loc[eligible, "radiance_2014"] - 1
    )
    assert np.allclose(expected, [-30, 100, 20])
    zero = trial.copy()
    zero.loc[0, "radiance_2014"] = 0
    assert not zero.radiance_2014.gt(0).iloc[0]
    t = np.array([-99., -10., 0., 1., 20., 300.])
    assert np.allclose(signed_inverse(signed_transform(t)), t)
    print("Figure 1 endpoints, baseline mask and signed-histogram test: PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint-csv", default=DEFAULT_ENDPOINT)
    parser.add_argument("--annual-parquet", default=DEFAULT_ANNUAL)
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--snow-2021-cache", default=DEFAULT_SNOW_OUTPUT)
    parser.add_argument("--rebuild-snow-2021", action="store_true",
                        help="Repeat approved 2021 overlay extraction even if cache verified")
    parser.add_argument("--simplify-map-m", type=float, default=150.0,
                        help="Vector map display-only simplification, projected metres")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    self_test()
    if args.self_test:
        return

    endpoint = as_repo_path(args.endpoint_csv)
    annual_source = as_repo_path(args.annual_parquet)
    out = as_repo_path(args.out)
    if not endpoint.is_file():
        raise FileNotFoundError(f"Final 1,522-region endpoint CSV is missing: {endpoint}")
    apply_figure_style()
    data = validate_endpoint(pd.read_csv(endpoint, dtype={"NUTS_ID": str}))
    out.mkdir(parents=True, exist_ok=True)
    summary = endpoint_summary(data)
    data[[
        "NUTS_ID", "radiance_2014", "radiance_2024",
        "radiance_change_abs", "radiance_change_pct",
        "radiance_change_pct_all", "percent_baseline_zero"
    ]].to_csv(
        out / "fig1B_and_hist_all_nuts3_percent_changes.csv", index=False
    )
    summary.update({
        "source_endpoint_csv": str(endpoint),
        "source_annual_parquet": str(annual_source),
        "map_display_simplification_metres": args.simplify_map_m,
        "geographic_extent_lonlat": list(EUROPE_EXTENT_LONLAT),
        "vector_output_formats": ["svg", "pdf"],
    })
    geom, bounds = display_geometry(load_config(), args.simplify_map_m)
    publication_map(
        data, geom, bounds, "radiance_change_abs",
        "A  |  Absolute radiance change, 2014–2024",
        f"2024 minus 2014 ({RADIANCE_UNIT})",
        out / "fig1A_absolute_radiance_change_2014_2024",
        linthresh=.1
    )
    publication_map(
        data, geom, bounds, "radiance_change_pct_all",
        "B  |  Percentage radiance change, 2014–2024",
        "Radiance change (%)",
        out / "fig1B_percent_radiance_change_2014_2024",
        linthresh=10.
    )
    publication_hist(
        data, summary, out / "fig1D_nuts3_percentage_change_distribution_2014_2024"
    )
    # Do annual data only after endpoint plots are generated: if its older
    # annual cohort is unavailable/inconsistent the completed A/B/D survive.
    (out / "fig1_provenance_and_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    snow, snow_audit = get_locked_snow_2021(
        as_repo_path(args.snow_2021_cache), rebuild=args.rebuild_snow_2021
    )
    summary["locked_2021_snow_source"] = snow_audit
    annual, annual_meta = load_annual(annual_source, data, snow)
    summary.update(annual_meta)
    annual.to_csv(out / "fig1C_europe_annual_radiance_support.csv", index=False)
    publication_annual(
        annual, out / "fig1C_europe_annual_radiance_2014_2024",
        annual_meta["annual_nuts3_balanced_cohort"]
    )
    (out / "fig1_provenance_and_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print("\nFINAL OUTPUT FILES", flush=True)
    for p in sorted(out.glob("fig1*")):
        print(f"  {p} ({p.stat().st_size / 1048576:.2f} MiB)", flush=True)
    print("\nNUTS3 mean percent change:", round(
        summary["mean_nuts3_percent_change_equal_region_weight"], 3))
    print("NUTS3 median percent change:", round(
        summary["median_nuts3_percent_change_equal_region_weight"], 3))
    print("Europe area-weighted endpoint change:", round(
        summary["europe_change_of_area_weighted_radiance_percent"], 3))


if __name__ == "__main__":
    main()
