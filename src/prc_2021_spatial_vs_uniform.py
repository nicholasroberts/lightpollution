#!/usr/bin/env python3
"""PRC sensitivity: spatial 2021 winter correction versus matched-mean NUTS uniform.

Uses out-of-region snow/albedo model predictions from the existing 18-region
monthly pilot. Reconstructs matched VIIRS pixels, applies each coarse-cell/
brightness-class predicted surface effect per pixel, and compares with a single
NUTS3 multiplicative factor calibrated to remove THE SAME MONTHLY MEAN radiance.

Monthly EOG is not the annual VNL product: this is a sensitivity test, NOT
a validated 2021 annual VIIRS or PRC correction. No original files overwritten.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask, geometry_window

from diagnose_fixed_nuts2016_monthly_product_atlas import discover_inputs
from extract_nuts3_radiance import nuts_path_for_release
from historical_prc_pixel_nuts3 import load_inputs as prc_inputs
from project_config import PROJECT_ROOT, configured_path, load_config
from test_fixed_nuts2016_snow_viirs_2021 import (
    discover_mod10cm, ensure_snow_geotiff,
)
from test_prc_2021_pixel_snow_albedo import (
    YEARS, MONTHS, reproject_to_region, equal_grid, path,
)


def weighted_mean(values, weights):
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    if not len(values) or not np.all(np.isfinite(values)) or not np.all(
        np.isfinite(weights) & (weights > 0)
    ):
        raise ValueError("Invalid weights or values")
    return float(np.average(values, weights=weights))


def response(raw, offset, parameters):
    yf, y0, k = parameters
    return yf + (y0 - yf) * np.exp(-k * np.maximum(np.asarray(raw) - offset, 0))


def spatial_vs_uniform(raw, effect, weights, offset, parameters):
    """Exactly equal corrected mean *raw monthly radiance*, nonlinear PRC."""
    raw = np.asarray(raw, float)
    effect = np.asarray(effect, float)
    corrected = np.maximum(raw - effect, 0.0)
    original_mean = weighted_mean(raw, weights)
    corrected_mean = weighted_mean(corrected, weights)
    factor = corrected_mean / original_mean if original_mean > 0 else 1.0
    uniform = raw * factor
    if not np.isclose(weighted_mean(uniform, weights), corrected_mean, atol=1e-8):
        raise AssertionError("Spatial and uniform radiance adjustments not matched")
    baseline_prc = weighted_mean(response(raw, offset, parameters), weights)
    spatial_prc = weighted_mean(response(corrected, offset, parameters), weights)
    uniform_prc = weighted_mean(response(uniform, offset, parameters), weights)
    return {
        "raw_mean_nw": original_mean,
        "corrected_mean_nw": corrected_mean,
        "uniform_factor_equal_mean": factor,
        "prc_before": baseline_prc,
        "prc_spatial": spatial_prc,
        "prc_uniform_equal_mean": uniform_prc,
        "delta_prc_spatial": spatial_prc - baseline_prc,
        "delta_prc_uniform_equal_mean": uniform_prc - baseline_prc,
        "delta_prc_spatial_minus_uniform": spatial_prc - uniform_prc,
        "fraction_negative_radiance_clipped": float(np.mean(raw - effect < 0)),
    }


def test_math():
    r = np.array([0.0, 1.0, 4.0, 30.0], float)
    w = np.array([1., 2., 2., 1.])
    p = (0.144845, 0.461175, 0.106483)
    z = spatial_vs_uniform(r, np.zeros(4), w, .235, p)
    assert np.isclose(z["delta_prc_spatial"], 0)
    assert np.isclose(z["delta_prc_uniform_equal_mean"], 0)
    c = spatial_vs_uniform(r, np.array([0.0, .1, 1.0, 4.0]), w, .235, p)
    assert np.isclose(c["raw_mean_nw"] * c["uniform_factor_equal_mean"],
                      c["corrected_mean_nw"])
    assert all(p[0] - 1e-7 <= v <= p[1] + 1e-7
               for v in (c["prc_spatial"], c["prc_uniform_equal_mean"]))
    print("PRC same-mean comparison synthetic test: PASS", flush=True)


def make_pixels(geom, rad, cf, snow, snow_native_grid):
    ref = rad[2021]
    window = geometry_window(ref, [geom])
    inside = geometry_mask([geom],
                           out_shape=(int(window.height), int(window.width)),
                           transform=ref.window_transform(window), invert=True)
    d = {}
    valid = inside.copy()
    for year in YEARS:
        r = np.ma.filled(rad[year].read(1, window=window, masked=True)
                         .astype("float64"), np.nan)
        c = np.ma.filled(cf[year].read(1, window=window, masked=True)
                         .astype("float64"), np.nan)
        s = reproject_to_region(snow[year], ref, window)
        d[year] = (r, c, s)
        valid &= (np.isfinite(r) & (r >= 0) & np.isfinite(c) & (c > 0)
                  & np.isfinite(s) & (s >= 0) & (s <= 100))
    if int(valid.sum()) < 20:
        return pd.DataFrame(), int(inside.sum()), int(valid.sum())
    rr, cc = np.indices(inside.shape)
    rr += int(window.row_off)
    cc += int(window.col_off)
    trans = ref.transform
    lon = trans.c + trans.a * (cc + .5) + trans.b * (rr + .5)
    lat = trans.f + trans.d * (cc + .5) + trans.e * (rr + .5)
    snow_trans, snow_width = snow_native_grid
    snow_row = np.floor((lat - snow_trans.f) / snow_trans.e).astype("int64")
    snow_col = np.floor((lon - snow_trans.c) / snow_trans.a).astype("int64")
    snow_cell = snow_row * snow_width + snow_col
    baseline = (d[2020][0] + d[2022][0]) / 2.0
    return pd.DataFrame({
        "snow_cell_id": snow_cell[valid].astype(str),
        "brightness_bin": np.digitize(baseline[valid], [1, 5, 20]),
        "r2021": d[2021][0][valid],
        "latitude": lat[valid],
        "weight": np.cos(np.deg2rad(lat[valid])),
    }), int(inside.sum()), int(valid.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--predictions",
        default="data/processed/historical_prc_2021_pixel_snow_test/"
                "nuts_disjoint_holdout_predictions.csv")
    ap.add_argument("--output-dir",
        default="data/processed/historical_prc_2021_prc_sensitivity")
    ap.add_argument("--model", choices=["snow_plus_albedo", "snow_only"],
                    default="snow_plus_albedo")
    ap.add_argument("--min-group-albedo", type=float, default=0.5,
                    help="Minimum share of QA-valid albedo pixels in each group")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    test_math()
    if args.self_test:
        return
    if not (0 <= args.min_group_albedo <= 1):
        raise ValueError("--min-group-albedo outside 0..1")

    config = load_config()
    _, parameters, offset, _ = prc_inputs(config, path("config/prc_transfer.yaml"))
    pred_path = path(args.predictions)
    if not pred_path.exists():
        raise FileNotFoundError(
            f"{pred_path} — run src/test_prc_2021_pixel_snow_albedo.py first.")
    pred = pd.read_csv(pred_path, dtype={"snow_cell_id": str})
    pred = pred[pred.model.eq(args.model) & pred.month.isin([1,2,3])].copy()
    keys = ["NUTS_ID", "month", "snow_cell_id", "brightness_bin"]
    if pred.duplicated(keys).any() or pred.empty:
        raise ValueError("Held-out model predictions missing or duplicate group IDs")
    pred["brightness_bin"] = pred.brightness_bin.astype(int)
    ids = sorted(pred.NUTS_ID.unique())
    if args.limit:
        ids = ids[:args.limit]
    pred = pred[pred.NUTS_ID.isin(ids)]
    nuts = gpd.read_file(nuts_path_for_release(config, 2016)).to_crs("EPSG:4326")
    nuts = nuts[nuts.NUTS_ID.isin(ids)]
    if len(nuts) != len(ids):
        raise ValueError("Missing NUTS2016 geometry for held-out region")
    pair, _ = discover_inputs(configured_path(config, "viirs_raw"))
    needed = {("vcmslcfg", y, m) for y in YEARS for m in MONTHS}
    if needed - set(pair):
        raise FileNotFoundError(f"Missing VIIRS monthly pairs: {needed-set(pair)}")
    if any("_j01_" in pair[key].radiance.name.lower() for key in needed):
        raise ValueError("Mixed NPP/J01 monthly products; must reprocess consistently")
    snow_root = path(config["analysis"]["fixed_nuts2016_snow_viirs_test"]
                     ["mod10cm_raw_dir"])
    snow_inventory = discover_mod10cm(snow_root)
    snow_cache = configured_path(config, "processed") / "mod10cm_061_geotiff_cache"

    rows, bands = [], []
    for month in MONTHS:
        print(f"PRC test month {month}/3", flush=True)
        with ExitStack() as stack:
            rad = {y: stack.enter_context(rasterio.open(pair["vcmslcfg", y, month].radiance))
                   for y in YEARS}
            cf = {y: stack.enter_context(rasterio.open(pair["vcmslcfg", y, month].cloudfree))
                  for y in YEARS}
            equal_grid([*rad.values(), *cf.values()])
            if rad[2021].crs.to_epsg() != 4326:
                raise ValueError("Native snow cell indexing needs VIIRS EPSG:4326")
            snow, native_grid = {}, None
            for y in YEARS:
                tif = ensure_snow_geotiff(snow_inventory[y,month],snow_cache,"EPSG:4326")
                snow[y] = stack.enter_context(rasterio.open(tif))
                grid = (snow[y].transform, snow[y].width)
                if native_grid is None:
                    native_grid = grid
                elif grid != native_grid:
                    raise ValueError("Native snow cell grids differ")
            for i, region in enumerate(nuts.itertuples(index=False), 1):
                subset = pred[(pred.NUTS_ID == region.NUTS_ID) &
                              (pred.month == month)]
                if subset.empty:
                    continue
                pixels, inside_count, support_count = make_pixels(
                    region.geometry, rad, cf, snow, native_grid)
                if pixels.empty:
                    continue
                joincols = ["snow_cell_id", "brightness_bin",
                            "surface_component_nw","albedo_valid_pixel_fraction"]
                merged = pixels.merge(subset[joincols],
                                      on=["snow_cell_id", "brightness_bin"],
                                      how="inner", validate="many_to_one")
                if args.model == "snow_plus_albedo":
                    merged = merged[
                        merged.albedo_valid_pixel_fraction >= args.min_group_albedo
                    ]
                if merged.empty:
                    continue
                matched_frac = len(merged) / len(pixels)
                if matched_frac < 0.05:
                    print(f"LOW matched fraction {region.NUTS_ID} month {month}: "
                          f"{matched_frac:.3f}", flush=True)
                for assumed_offset in (0.0, offset):
                    z = spatial_vs_uniform(
                        merged.r2021.to_numpy(float),
                        merged.surface_component_nw.to_numpy(float),
                        merged.weight.to_numpy(float),
                        assumed_offset, parameters)
                    rows.append({
                        "model": args.model, "NUTS_ID": region.NUTS_ID,
                        "NUTS_NAME": getattr(region, "NUTS_NAME", ""),
                        "CNTR_CODE": getattr(region, "CNTR_CODE", ""),
                        "month": month, "offset_sensitivity_nw": assumed_offset,
                        "pixels_with_heldout_prediction":len(merged),
                        "pixels_common_3years_snow":len(pixels),
                        "total_polygon_pixels":inside_count,
                        "predicted_pixel_fraction_common": matched_frac,
                        **z,
                    })
                    for b, sub in merged.groupby("brightness_bin"):
                        bz = spatial_vs_uniform(
                            sub.r2021.to_numpy(float),
                            sub.surface_component_nw.to_numpy(float),
                            sub.weight.to_numpy(float), assumed_offset,
                            parameters)
                        bands.append({
                            "model":args.model,"NUTS_ID":region.NUTS_ID,
                            "month":month,"brightness_bin":int(b),
                            "offset_sensitivity_nw":assumed_offset,
                            "n_pixels":len(sub), **bz,
                        })
                print(f" {month:02d} {i:02d}/{len(nuts)} {region.NUTS_ID}: "
                      f"{len(merged)}/{len(pixels)} matched pixels",flush=True)

    if not rows:
        raise RuntimeError("No PRC estimates; inspect matched group IDs and data masks")
    output = path(args.output_dir)
    output.mkdir(parents=True,exist_ok=True)
    tab = pd.DataFrame(rows)
    sub = pd.DataFrame(bands)
    tab.to_csv(output/"prc_spatial_vs_equal_mean_uniform_regions.csv",index=False)
    sub.to_csv(output/"prc_spatial_vs_equal_mean_uniform_brightness.csv",index=False)
    summary = [
        "2021 monthly, matched-satellite-pixel PRC SENSITIVITY ONLY",
        "Spatial surface estimates are out-of-NUTS3-held-out snow/albedo models.",
        "Within each NUTS and month, uniform factor preserves EXACT SAME",
        "area-weighted corrected monthly raw radiance as spatial adjustment.",
        "PRC is applied on every valid 2021 VIIRS pixel BEFORE aggregation.",
        "The 0.235 subtraction from annual VNL is a separate sensitivity",
        "scenario here: monthly radiance may be on a different product scale.",
        "Neither option is a validated annual 2021 snow correction.",
        "This test does NOT compare annual 2021 NUTS factor magnitudes directly",
        "with January/February monthly model magnitudes; temporal supports differ.",
        f"Regions sampled={tab.NUTS_ID.nunique()}, region-month comparisons="
        f"{tab.groupby(['NUTS_ID','month']).ngroups}, rows={len(tab)}",
    ]
    for ofs, group in tab.groupby("offset_sensitivity_nw"):
        w = group.pixels_with_heldout_prediction.to_numpy(float)
        for month, monthly in group.groupby("month"):
            vals = monthly.delta_prc_spatial_minus_uniform.to_numpy()
            summary.append(
                f"offset={ofs:.3f} month={month}: n_region={len(monthly)}, "
                f"median absolute(spatial - uniform) mean PRC="
                f"{np.median(np.abs(vals)):.7f}, "
                f"max absolute difference={np.max(np.abs(vals)):.7f}, "
                f"median matched VIIRS pixel fraction="
                f"{np.median(monthly.predicted_pixel_fraction_common):.3f}"
            )
    (output/"prc_spatial_vs_uniform_summary.txt").write_text(
        "\n".join(summary)+"\n",encoding="utf-8")

    base = tab[np.isclose(tab.offset_sensitivity_nw, offset)]
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.8),layout="constrained")
    markers={1:"o",2:"s",3:"^"}
    for month, part in base.groupby("month"):
        axs[0].scatter(part.delta_prc_uniform_equal_mean,
                       part.delta_prc_spatial,marker=markers[month],
                       label=f"Month {month}",alpha=.8)
        axs[1].scatter(part.raw_mean_nw,
                       part.delta_prc_spatial_minus_uniform,
                       marker=markers[month],label=f"Month {month}",alpha=.8)
    mn = min(base.delta_prc_uniform_equal_mean.min(),base.delta_prc_spatial.min())
    mx = max(base.delta_prc_uniform_equal_mean.max(),base.delta_prc_spatial.max())
    axs[0].plot([mn,mx],[mn,mx],linestyle="--",color="gray")
    axs[0].set(xlabel="Uniform adjustment: change in mean PRC",
               ylabel="Spatial adjustment: change in mean PRC",
               title="Matched radiance decrement, different spatial allocation")
    axs[1].axhline(0,color="gray",linewidth=.8)
    axs[1].set(xscale="symlog",xlabel="Region-month mean raw VIIRS radiance",
               ylabel="Mean PRC: spatial minus uniform",
               title="Biological importance of spatial redistribution")
    for ax in axs:
        ax.grid(alpha=.2)
        ax.legend()
    fig.savefig(output/"prc_spatial_vs_uniform_2021_pilot.pdf",bbox_inches="tight")
    plt.close(fig)
    print("\n".join(summary))
    print("Outputs:",output,flush=True)


if __name__ == "__main__":
    main()
