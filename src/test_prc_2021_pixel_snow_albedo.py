#!/usr/bin/env python3
"""Pixel-matched 2021 VIIRS / MOD10CM snow / VNP43MA3 albedo pilot.

OBSERVATIONAL DIAGNOSTIC ONLY: this does not change annual PRC or produce a
validated corrected annual pixel field.  Uses NUTS-disjoint validation and
spatial MOD10CM-cell aggregation to reduce 500m-pixel pseudoreplication.
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
from rasterio.enums import Resampling
from rasterio.features import geometry_mask, geometry_window
from rasterio.vrt import WarpedVRT
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

from diagnose_fixed_nuts2016_monthly_product_atlas import discover_inputs
from extract_nuts3_radiance import nuts_path_for_release
from project_config import PROJECT_ROOT, configured_path, load_config
from test_fixed_nuts2016_snow_viirs_2021 import (
    discover_mod10cm, ensure_snow_geotiff,
)

YEARS = (2020, 2021, 2022)
MONTHS = (1, 2, 3)
MODEL_FORMULA = (
    "delta_radiance ~ snow_delta_pp + snow_delta_pp:np.log1p(baseline_radiance)"
    " + albedo_delta_pp + albedo_delta_pp:np.log1p(baseline_radiance)"
    " + delta_cloudfree + is_february"
)
SNOW_FORMULA = (
    "delta_radiance ~ snow_delta_pp + snow_delta_pp:np.log1p(baseline_radiance)"
    " + delta_cloudfree + is_february"
)


def path(x):
    p = Path(x)
    return p if p.is_absolute() else PROJECT_ROOT / p


def select_regions(geo, regional, limit):
    d = pd.read_csv(regional)
    d = d.loc[
        d.month.eq(2) &
        d.albedo_type.astype(str).str.upper().eq("WSA") &
        d.quality_mode.astype(str).str.lower().eq("qa01")
    ].copy()
    key = "snow_anomaly_2021_vs_flanks_pp"
    if key not in d:
        raise ValueError(f"Missing regional snow anomaly {key} in {regional}")
    d = d.drop_duplicates("NUTS_ID")
    geo = geo.copy()
    projected = geo.to_crs("EPSG:3035")
    geo["area_km2"] = projected.area / 1e6
    # Centroids calculated in a projected CRS and transformed to lon/lat.
    centres = projected.geometry.centroid.to_crs("EPSG:4326")
    geo = geo.loc[
        geo.area_km2.between(100, 2500) &
        centres.x.between(-12, 32) &
        centres.y.between(36, 72)
    ]
    d = geo.merge(d[["NUTS_ID", key]], on="NUTS_ID", how="inner")
    d = d[np.isfinite(d[key])].copy()
    high = d.sort_values(key, ascending=False)
    near_zero = d.iloc[np.argsort(np.abs(d[key].to_numpy()))]
    low = d.sort_values(key)
    chosen, countries = [], {}
    lists = [list(high.index), list(near_zero.index), list(low.index)]
    while len(chosen) < limit and any(lists):
        progressed = False
        for group in lists:
            while group:
                ix = group.pop(0)
                if ix in chosen:
                    continue
                c = str(d.loc[ix, "CNTR_CODE"])
                if countries.get(c, 0) >= 3:
                    continue
                chosen.append(ix)
                countries[c] = countries.get(c, 0) + 1
                progressed = True
                break
            if len(chosen) >= limit:
                break
        if not progressed:
            break
    if len(chosen) < min(limit, 6):
        raise RuntimeError("Not enough eligible contrasting NUTS regions")
    picked = d.loc[chosen].copy()
    picked["selection_snow_group"] = np.where(
        picked[key] > 2, "snow_increase",
        np.where(picked[key] < -2, "snow_decrease", "near_zero"))
    return picked


def equal_grid(sources):
    ref = sources[0]
    for s in sources[1:]:
        if (ref.crs != s.crs or ref.width != s.width or
                ref.height != s.height or
                not ref.transform.almost_equals(s.transform)):
            raise ValueError(f"Monthly VIIRS pixel grid differs: {ref.name} / {s.name}")


def arrays_for_region(geom, rad, cf, snow, albedo, native_snow_grid):
    ref = rad[2021]
    if ref.crs is None or ref.crs.to_epsg() != 4326:
        raise ValueError("MOD10CM mapping expects EPSG:4326 VIIRS pixels")
    try:
        w = geometry_window(ref, [geom], pad_x=0, pad_y=0)
    except rasterio.errors.WindowError:
        return None
    inside = geometry_mask([geom], out_shape=(int(w.height), int(w.width)),
                           transform=ref.window_transform(w), invert=True)
    if inside.sum() < 20:
        return None
    data = {}
    for year in YEARS:
        a = rad[year].read(1, window=w, masked=True)
        b = cf[year].read(1, window=w, masked=True)
        r = np.ma.filled(a.astype("float64"), np.nan)
        c = np.ma.filled(b.astype("float64"), np.nan)
        data[year] = {"r":r, "c":c,
                      "snow":np.ma.filled(snow[year].read(1, window=w, masked=True)
                                         .astype("float64"), np.nan),
                      "albedo":np.ma.filled(albedo[year].read(1,window=w,masked=True)
                                           .astype("float64"), np.nan)}
    common = inside.copy()
    for year in YEARS:
        common &= (np.isfinite(data[year]["r"]) &
                   (data[year]["r"] >= 0) &
                   np.isfinite(data[year]["c"]) & (data[year]["c"] > 0))
    support_frac = float(common.sum() / max(1, inside.sum()))
    # MOD10CM 0..100 are physical snow percent; categorical codes are invalid.
    for year in YEARS:
        snowval = data[year]["snow"]
        common &= np.isfinite(snowval) & (snowval >= 0) & (snowval <= 100)
    if common.sum() < 20:
        return None
    shape = inside.shape
    rows, cols = np.indices(shape)
    rows = rows + int(w.row_off)
    cols = cols + int(w.col_off)
    t = ref.transform
    lon = t.c + t.a * (cols + .5) + t.b * (rows + .5)
    lat = t.f + t.d * (cols + .5) + t.e * (rows + .5)
    # snow[2021] is a WarpedVRT on the VIIRS grid. Using its transform
    # would spuriously define a separate snow observation per VIIRS pixel.
    st, snow_width = native_snow_grid
    sr = np.floor((lat-st.f) / st.e).astype("int64")
    sc = np.floor((lon-st.c) / st.a).astype("int64")
    cell = sr * snow_width + sc
    flank = .5 * (data[2020]["r"]+data[2022]["r"])
    snow_delta = data[2021]["snow"] - .5*(data[2020]["snow"]+data[2022]["snow"])
    albedo_delta = 100 * (data[2021]["albedo"] -
                           .5*(data[2020]["albedo"]+data[2022]["albedo"]))
    cf_delta = data[2021]["c"] - .5*(data[2020]["c"]+data[2022]["c"])
    outcome = data[2021]["r"]-flank
    brightness_bins = np.digitize(flank, [1, 5, 20])
    # Spherical area is proportional to cos(latitude) for regular lon/lat cells.
    wt = np.cos(np.deg2rad(lat))
    albedo_ok = np.ones(shape, dtype=bool)
    for year in YEARS:
        x = data[year]["albedo"]
        albedo_ok &= np.isfinite(x) & (x >= 0) & (x <= 1)
    sel = common
    pixel = pd.DataFrame({
        "snow_cell_id":cell[sel],
        "brightness_bin":brightness_bins[sel],
        "pixel_weight":wt[sel],
        "delta_radiance":outcome[sel],
        "baseline_radiance":flank[sel],
        "snow_delta_pp":snow_delta[sel],
        "albedo_delta_pp":np.where(albedo_ok, albedo_delta, np.nan)[sel],
        "delta_cloudfree":cf_delta[sel],
        "albedo_valid":albedo_ok[sel],
    })
    # Collapse many VIIRS pixels to coarse snow cells and four fixed brightness
    # classes; each snow cell remains the spatial unit for inference.
    agg = []
    for (cell_id, bb), g in pixel.groupby(["snow_cell_id", "brightness_bin"],
                                        sort=False):
        weights = g.pixel_weight.to_numpy()
        def wm(col, subset=None):
            part = g if subset is None else g.loc[subset]
            if part.empty: return np.nan
            return np.average(part[col], weights=part.pixel_weight)
        valid_albedo = g.albedo_valid & g.albedo_delta_pp.notna()
        agg.append({
            "snow_cell_id":str(int(cell_id)),
            "brightness_bin":int(bb),
            "n_pixels":len(g),
            "area_weight":float(weights.sum()),
            "delta_radiance":wm("delta_radiance"),
            "baseline_radiance":wm("baseline_radiance"),
            "snow_delta_pp":wm("snow_delta_pp"),
            "albedo_delta_pp":wm("albedo_delta_pp", valid_albedo),
            "delta_cloudfree":wm("delta_cloudfree"),
            "albedo_valid_pixel_fraction":float(valid_albedo.mean()),
        })
    return pd.DataFrame(agg), int(inside.sum()), int(common.sum()), support_frac


def fit_and_predict(d, formula, train_mask, test_mask):
    train = d.loc[train_mask].copy()
    test = d.loc[test_mask].copy()
    predictors = (["snow_delta_pp", "baseline_radiance",
                   "delta_cloudfree"] +
                  (["albedo_delta_pp"] if "albedo_delta_pp" in formula else []))
    train = train.dropna(subset=predictors+["delta_radiance"])
    test = test.dropna(subset=predictors+["delta_radiance"])
    if train.NUTS_ID.nunique() < 3 or len(train) < 20:
        raise ValueError("Insufficient independent regions/cells for pilot model")
    model = smf.wls(formula, data=train, weights=train.area_weight).fit()
    out = test.copy()
    out["predicted_total_anomaly"] = model.predict(out)
    p = model.params
    lum = np.log1p(out.baseline_radiance.to_numpy())
    snow_effect = (
        p.get("snow_delta_pp", 0.0) * out.snow_delta_pp.to_numpy() +
        p.get("snow_delta_pp:np.log1p(baseline_radiance)", 0.0) *
        out.snow_delta_pp.to_numpy() * lum
    )
    albedo_effect = (
        p.get("albedo_delta_pp", 0.0) * out.albedo_delta_pp.to_numpy() +
        p.get("albedo_delta_pp:np.log1p(baseline_radiance)", 0.0) *
        out.albedo_delta_pp.to_numpy() * lum
    ) if "albedo_delta_pp" in formula else np.zeros(len(out))
    out["surface_component_nw"] = snow_effect + albedo_effect
    out["outcome_minus_surface"] = out.delta_radiance - out.surface_component_nw
    return model, out


def lobo(d, formula):
    chunks = []
    # Leave entire NUTS3 polygons out, never individual nearby pixels.
    for held in sorted(d.NUTS_ID.unique()):
        train = d.NUTS_ID.ne(held) & d.month.isin([1, 2])
        test = d.NUTS_ID.eq(held)
        try:
            _, pred = fit_and_predict(d, formula, train, test)
        except ValueError:
            continue
        pred["held_out_NUTS_ID"] = held
        chunks.append(pred)
    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


def test_functions():
    v = pd.DataFrame({
        "NUTS_ID":["A"]*4+["B"]*4+["C"]*4+["D"]*4,
        "month":[1,1,2,2]*4,
        "area_weight":[1.0]*16,
        "delta_radiance":[0.0,1.0,1.0,2.0]*4,
        "snow_delta_pp":[0,10,0,10]*4,
        "albedo_delta_pp":[0,1,0,1]*4,
        "baseline_radiance":[1,2,3,4]*4,
        "delta_cloudfree":[1,0,1,0]*4,
    })
    if np.abs((v.delta_radiance-v.snow_delta_pp/10).mean()-.5)>1e-10:
        raise AssertionError("Basic synthetic model sanity check")
    print("Synthetic preprocessing sanity check: PASS")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config/pipeline.yaml")
    ap.add_argument("--n-regions",type=int,default=18)
    ap.add_argument("--limit",type=int,default=0,help="Smoke-test subset")
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--output-dir",default=
                    "data/processed/historical_prc_2021_pixel_snow_test")
    args = ap.parse_args()
    test_functions()
    if args.self_test:
        return
    config = load_config(args.config)
    out = path(args.output_dir)
    out.mkdir(parents=True,exist_ok=True)
    nuts = gpd.read_file(nuts_path_for_release(config,2016)).to_crs("EPSG:4326")
    selected = select_regions(
        nuts, configured_path(config,"processed") /
        "fixed_nuts2016_albedo_snow_viirs_anomaly_regional.csv",
        args.limit if args.limit else args.n_regions)
    selected.drop(columns=["geometry"]).to_csv(
        out / "selected_nuts3_pilot_regions.csv",index=False)
    pairs, missing = discover_inputs(configured_path(config,"viirs_raw"))
    expected = {("vcmslcfg",y,m) for y in YEARS for m in MONTHS}
    if expected-set(pairs):
        raise FileNotFoundError(f"Missing monthly VIIRS: {sorted(expected-set(pairs))}")
    # Explicitly reject mixed satellite platforms for the matched monthly test.
    for key in expected:
        if "_j01_" in pairs[key].radiance.name.lower():
            raise ValueError(f"Mixed-platform monthly product: {pairs[key].radiance}")

    snow_dir = configured_path(config,"processed") / "mod10cm_061_geotiff_cache"
    # Raw-path config is relative to PROJECT_ROOT, never processed/ parent.
    modis_raw_dir = path(
        config["analysis"]["fixed_nuts2016_snow_viirs_test"]["mod10cm_raw_dir"])
    if not modis_raw_dir.is_dir():
        raise FileNotFoundError(f"MOD10CM directory not found: {modis_raw_dir}")
    print(f"MOD10CM inventory: {modis_raw_dir}",flush=True)
    snow_raw = discover_mod10cm(modis_raw_dir)
    albedo_dir = configured_path(config,"processed") / "vnp43ma3_002_monthly_composite_cache"
    albedo_paths = {(y,m):albedo_dir / f"vnp43ma3_002_{y}_{m:02d}_monthly.tif"
                    for y in YEARS for m in MONTHS}
    if any(not p.exists() for p in albedo_paths.values()):
        raise FileNotFoundError(
            "Required QA-filtered VNP43MA3 monthly composites missing. "
            "Run the existing albedo test first. Missing: " +
            ", ".join(str(p) for p in albedo_paths.values() if not p.exists()))

    record, counts = [], []
    for month in MONTHS:
        print(f"Month {month}: common-pixel VIIRS + MOD10CM + VNP43MA3",flush=True)
        with ExitStack() as stack:
            rd = {y:stack.enter_context(rasterio.open(pairs["vcmslcfg",y,month].radiance))
                  for y in YEARS}
            cf = {y:stack.enter_context(rasterio.open(pairs["vcmslcfg",y,month].cloudfree))
                  for y in YEARS}
            equal_grid([*rd.values(),*cf.values()])
            snowds={}
            native_snow_grid = None
            for y in YEARS:
                hdf = snow_raw[(y,month)]
                snowtif = ensure_snow_geotiff(hdf,snow_dir,"EPSG:4326")
                src = stack.enter_context(rasterio.open(snowtif))
                source_grid = (src.transform,src.width)
                if native_snow_grid is None:
                    native_snow_grid = source_grid
                elif native_snow_grid != source_grid:
                    raise ValueError("Native MOD10CM source grids differ by year")
                snowds[y] = stack.enter_context(WarpedVRT(
                    src, crs=rd[2021].crs, transform=rd[2021].transform,
                    width=rd[2021].width,height=rd[2021].height,
                    resampling=Resampling.nearest, src_nodata=255,
                    nodata=255))
            albedods={}
            for y in YEARS:
                src = stack.enter_context(rasterio.open(albedo_paths[y,month]))
                if src.count < 1 or src.descriptions[0] != "wsa_qa01_mean":
                    raise ValueError("Albedo composite bands are not the expected QA0/1 WSA")
                albedods[y] = stack.enter_context(WarpedVRT(
                    src, crs=rd[2021].crs, transform=rd[2021].transform,
                    width=rd[2021].width,height=rd[2021].height,
                    resampling=Resampling.nearest, nodata=np.nan))
            for n,row in enumerate(selected.itertuples(index=False),1):
                try:
                    result = arrays_for_region(
                        row.geometry,rd,cf,snowds,albedods,native_snow_grid)
                except rasterio.errors.WindowError:
                    result=None
                if result is None:
                    counts.append({"NUTS_ID":row.NUTS_ID,"month":month,
                                   "status":"no_snow_viirs_support"})
                    continue
                g,n_inside,n_common,frac = result
                g["NUTS_ID"] = row.NUTS_ID
                g["CNTR_CODE"] = row.CNTR_CODE
                g["NUTS_NAME"] = row.NUTS_NAME
                g["selection_snow_group"] = row.selection_snow_group
                g["month"]=month
                record.append(g)
                counts.append({"NUTS_ID":row.NUTS_ID,"month":month,
                               "n_polygon_pixels":n_inside,
                               "n_viirs_snow_pixels":n_common,
                               "common_viirs_support_fraction":frac,
                               "n_coarse_groups":len(g),"status":"ok"})
                print(f"  {month:02d} {n:02d}/{len(selected)} {row.NUTS_ID}: "
                      f"{n_common} pixels, {len(g)} coarse groups",flush=True)
    coverage = pd.DataFrame(counts)
    coverage.to_csv(out/"spatial_support_quality.csv",index=False)
    if not record:
        raise RuntimeError("No matched VIIRS pixels: see spatial_support_quality.csv")
    groups = pd.concat(record,ignore_index=True)
    # March is held-out: C(month) fitted to Jan/Feb cannot predict unseen
    # March. Use a February indicator, with Jan as March's reference intercept.
    groups["is_february"] = groups.month.eq(2).astype(int)
    groups.to_csv(out/"matched_snowcell_brightness_groups.csv",index=False)

    winter = groups[groups.month.isin([1,2])].copy()
    models, holdouts = [], []
    for label, formula in [("snow_only",SNOW_FORMULA),
                           ("snow_plus_albedo",MODEL_FORMULA)]:
        cols=["delta_radiance","snow_delta_pp","baseline_radiance",
              "delta_cloudfree"] + (["albedo_delta_pp"] if label=="snow_plus_albedo"
                                    else [])
        sample = groups.dropna(subset=cols).copy()
        train = sample.month.isin([1,2])
        if sample.loc[train,"NUTS_ID"].nunique()<3:
            continue
        fit, predictions = fit_and_predict(sample,formula,train,
                                           sample.index.isin(sample.index))
        params=fit.params
        for key,value in params.items():
            models.append({"model":label,"term":key,"coefficient":value,
                           "p_value_descriptive_only":fit.pvalues[key],
                           "training_n_groups":sum(train),
                           "training_n_regions":sample.loc[train,"NUTS_ID"].nunique()})
        hold = lobo(sample,formula)
        if not hold.empty:
            hold["model"]=label
            holdouts.append(hold)
    pd.DataFrame(models).to_csv(out/"snow_albedo_model_coefficients.csv",
                               index=False)
    if not holdouts:
        raise RuntimeError("Insufficient model holdouts; add pilot regions")
    held = pd.concat(holdouts,ignore_index=True)
    held.to_csv(out/"nuts_disjoint_holdout_predictions.csv",index=False)

    diagnostics=[]
    for (model,month),v in held.groupby(["model","month"]):
        w=v.area_weight.to_numpy(float)
        err=v.delta_radiance.to_numpy(float)
        corrected=v.outcome_minus_surface.to_numpy(float)
        diagnostics.append({
            "model":model,"month":int(month),
            "role":"winter_fit_period" if month<=2 else "March_negative_control",
            "n_regions":v.NUTS_ID.nunique(),
            "n_snow_cell_brightness_groups":len(v),
            "rmse_raw_to_zero_nw":np.sqrt(np.average(err**2,weights=w)),
            "rmse_after_surface_subtraction_nw":
                np.sqrt(np.average(corrected**2,weights=w)),
            "mean_predicted_surface_component_nw":np.average(
                v.surface_component_nw,weights=w),
            "mean_observed_radiance_anomaly_nw":np.average(err,weights=w),
            "mean_residual_after_surface_nw":np.average(corrected,weights=w),
        })
    diag=pd.DataFrame(diagnostics)
    diag.to_csv(out/"spatial_holdout_diagnostics.csv",index=False)

    region_rows=[]
    for (model, nuts_id, month), g in held.groupby(["model","NUTS_ID","month"]):
        region_rows.append({
            "model":model,"NUTS_ID":nuts_id,"month":month,
            "predicted_surface_component_nw":np.average(
                g.surface_component_nw,weights=g.area_weight),
            "observed_month_anomaly_nw":np.average(
                g.delta_radiance,weights=g.area_weight),
            "residual_after_surface_nw":np.average(
                g.outcome_minus_surface,weights=g.area_weight),
            "n_groups":len(g),
        })
    region=pd.DataFrame(region_rows)
    region.to_csv(out/"regional_monthly_spatial_contrast.csv",index=False)

    fig,axs=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    for name,sub in held.groupby("model"):
        win=sub[sub.month.eq(2)]
        axs[0].scatter(win.snow_delta_pp,win.delta_radiance,
                       alpha=.2,s=8,label=name)
        axs[1].scatter(win.delta_radiance,win.outcome_minus_surface,
                       alpha=.2,s=8,label=name)
    axs[0].axhline(0,linewidth=.7,color="gray")
    axs[1].axline((0,0),(1,1),color="gray",linestyle="--")
    axs[0].set(xlabel="February MOD10CM snow anomaly (percentage points)",
               ylabel="Same-pixel VIIRS radiance anomaly (nW cm⁻² sr⁻¹)",
               title="February matched-pixel diagnostic")
    axs[1].set(xlabel="Held-out original radiance anomaly",
               ylabel="After predicted surface component removed",
               title="Held-out regions; observational, not corrected annual VNL")
    for ax in axs:
        ax.grid(alpha=.2)
        ax.legend()
    fig.savefig(out/"pixel_snow_albedo_spatial_holdout.pdf")
    plt.close(fig)
    summary=[
        "2021 matched-month VIIRS/snow/albedo PILOT (no annual correction applied)",
        f"NUTS2016 regions chosen from independent February snow anomaly: {len(selected)}",
        "VIIRS VCMSLCFG 2020/21/22 same calendar month, matched cf_cvg>0 pixels.",
        "MOD10CM at 0.05 degree: NOT 500m snow mapping; grouped by original cell.",
        "Albedo: VNP43MA3 WSA visible monthly QA0/1, native sinusoidal warped nearest",
        "to VIIRS grid. Reported apparent 500m values retain their source support.",
        "2022 J01 monthly products explicitly rejected; no monthly 0.235 subtracted.",
        "2017+ -0.235 annual VNL calibration remains part of original PRC histories.",
        "Train Jan/Feb with February indicator; March uses January intercept.",
        "March is a temporal negative control; leave whole NUTS3 regions out.",
        "The surface-effect regression is observational and not a causal correction.",
        "No monthly-to-annual correction or corrected PRC produced here.",
        "Do not select model solely for making 2021 anomaly smaller.",
        "",
        diag.to_string(index=False),
    ]
    (out/"pixel_snow_albedo_summary.txt").write_text(
        "\n".join(summary)+"\n",encoding="utf-8")
    print("\n".join(summary))
    print("Outputs:",out)


if __name__=="__main__":
    main()
