#!/usr/bin/env python3
"""August 2020–2022 native-snow-cell snow-free control for 2021 VIIRS.

This is an observational, pre-lock QC of the 2021 winter snow/albedo effect.
It uses the SAME 18 NUTS2016 pilot regions, matched EOG VCMSLCFG VIIRS pixels,
and native 0.05-degree MOD10CM snow cells in all three Augusts. August 2022
VIIRS is only used if its sensor matches August 2020 and 2021 (NPP vs J01).
It requires no August visible-albedo data and makes no snow/albedo correction.
"""
from __future__ import annotations

import argparse
import json
import re
from contextlib import ExitStack
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask, geometry_window

from diagnose_fixed_nuts2016_monthly_product_atlas import discover_inputs
from extract_nuts3_radiance import nuts_path_for_release
from project_config import PROJECT_ROOT, configured_path, load_config
from test_fixed_nuts2016_snow_viirs_2021 import (
    discover_mod10cm, ensure_snow_geotiff,
)
from test_prc_2021_pixel_snow_albedo import equal_grid, reproject_to_region

YEARS=(2020,2021,2022)
MONTH=8
BINS=(1.0,5.0,20.0)


def platform(p):
    match=re.search(r"(?<![a-z0-9])(npp|j01)(?![a-z0-9])",
                    Path(p).name.lower())
    return match.group(1) if match else "unknown"


def check_control_threshold(threshold):
    if not 0 <= threshold <= 5:
        raise ValueError("Snow-free cover tolerance must be 0..5 percent")


def self_test():
    check_control_threshold(1.)
    snow=np.array([[0,0,0],[0,.5,.8],[0,2,0],[0,np.nan,0],[0,100,0]])
    valid=np.isfinite(snow).all(axis=1) & ((snow>=0)&(snow<=100)).all(axis=1)
    free=valid & (snow <= 1).all(axis=1)
    assert free.tolist()==[True,True,False,False,False]
    assert platform("SVDNB_npp_20210801.vcmslcfg.tif")=="npp"
    assert platform("SVDNB_j01_20220801.vcmslcfg.tif")=="j01"
    print("August control QA self-test: PASS",flush=True)


def weighted_mean(x,w):
    a=np.asarray(x,float)
    wt=np.asarray(w,float)
    ok=np.isfinite(a)&np.isfinite(wt)&(wt>0)
    return float(np.average(a[ok],weights=wt[ok])) if np.any(ok) else float("nan")


def group_region(nuts_id,geom,radiance,cf,snow,grid,same_2022,threshold):
    ref=radiance[2021]
    try:
        win=geometry_window(ref,[geom],pad_x=0,pad_y=0)
    except rasterio.errors.WindowError:
        return None, {"NUTS_ID":nuts_id,"reason":"outside_raster"}
    transform=ref.window_transform(win)
    inside=geometry_mask([geom],out_shape=(int(win.height),int(win.width)),
                         transform=transform,invert=True)
    region_pixels=int(inside.sum())
    if region_pixels<20:
        return None,{"NUTS_ID":nuts_id,"reason":"too_few_polygon_pixels"}
    r={}
    c={}
    sn={}
    for year in YEARS:
        sn[year]=reproject_to_region(snow[year],ref,win)
    active=(2020,2021,2022) if same_2022 else (2020,2021)
    for year in active:
        r[year]=np.ma.filled(
            radiance[year].read(1,window=win,masked=True).astype(float),np.nan)
        c[year]=np.ma.filled(
            cf[year].read(1,window=win,masked=True).astype(float),np.nan)
    support=inside.copy()
    for year in active:
        support &= (np.isfinite(r[year]) & (r[year]>=0) &
                    np.isfinite(c[year]) & (c[year]>0))
    n_viirs=int(support.sum())
    valid_snow=support.copy()
    for year in YEARS:
        valid_snow &= (np.isfinite(sn[year]) &
                       (sn[year]>=0)&(sn[year]<=100))
    n_valid=int(valid_snow.sum())
    if n_valid<20:
        return None,{"NUTS_ID":nuts_id,"reason":"insufficient_snow_matched_pixels",
                     "n_polygon_pixels":region_pixels,"n_viirs_common":n_viirs,
                     "n_snow_valid":n_valid}
    snowfree=valid_snow.copy()
    for year in YEARS:
        snowfree &= sn[year]<=threshold
    rr,cc=np.indices(inside.shape)
    rr=rr+int(win.row_off)
    cc=cc+int(win.col_off)
    t=ref.transform
    lat=t.f+t.d*(cc+.5)+t.e*(rr+.5)
    lon=t.c+t.a*(cc+.5)+t.b*(rr+.5)
    native, width=grid
    cellr=np.floor((lat-native.f)/native.e).astype(np.int64)
    cellc=np.floor((lon-native.c)/native.a).astype(np.int64)
    cell=cellr*width+cellc
    baseline=(r[2020]+r[2022])/2 if same_2022 else r[2020]
    delta=(r[2021]-baseline)
    valid=valid_snow
    data=pd.DataFrame({
        "NUTS_ID":nuts_id,
        "snow_cell_id":cell[valid].astype(str),
        "brightness_bin":np.digitize(baseline[valid],BINS),
        "pixel_weight":np.cos(np.deg2rad(lat[valid])),
        "radiance_2020":r[2020][valid],
        "radiance_2021":r[2021][valid],
        "radiance_reference":baseline[valid],
        "radiance_delta":delta[valid],
        "cf_cvg_2020":c[2020][valid],
        "cf_cvg_2021":c[2021][valid],
        "snow_2020":sn[2020][valid],
        "snow_2021":sn[2021][valid],
        "snow_2022":sn[2022][valid],
        "snowfree_all_three":snowfree[valid],
    })
    groups=[]
    for (cellid,bb),d in data.groupby(["snow_cell_id","brightness_bin"],
                                      sort=False):
        w=d.pixel_weight.to_numpy(float)
        g={"NUTS_ID":nuts_id,"month":8,"snow_cell_id":cellid,
           "brightness_bin":int(bb),"n_pixels":len(d),
           "area_weight":float(w.sum()),
           "snowfree_all_three":bool(d.snowfree_all_three.all())}
        for col in ("radiance_2020","radiance_2021",
                    "radiance_reference","radiance_delta",
                    "cf_cvg_2020","cf_cvg_2021",
                    "snow_2020","snow_2021","snow_2022"):
            g[col]=weighted_mean(d[col],w)
        groups.append(g)
    return pd.DataFrame(groups),{
        "NUTS_ID":nuts_id,"reason":"ok",
        "n_polygon_pixels":region_pixels,
        "n_viirs_common":n_viirs,
        "n_snow_valid":n_valid,
        "n_snowfree_pixels":int(snowfree.sum()),
        "snowfree_fraction_of_valid":float(snowfree.sum()/n_valid),
        "viirs_matched_fraction_of_polygon":float(n_viirs/region_pixels),
    }


def summary_regions(groups):
    rows=[]
    for (region,is_free),g in groups.groupby(["NUTS_ID","snowfree_all_three"]):
        w=g.area_weight
        rows.append({
            "NUTS_ID":region,"snowfree_all_three":bool(is_free),
            "n_native_snowcell_brightness_groups":len(g),
            "n_pixels":int(g.n_pixels.sum()),
            "area_weight":float(w.sum()),
            "radiance_2020":weighted_mean(g.radiance_2020,w),
            "radiance_2021":weighted_mean(g.radiance_2021,w),
            "radiance_reference":weighted_mean(g.radiance_reference,w),
            "radiance_delta":weighted_mean(g.radiance_delta,w),
            "snow_2021":weighted_mean(g.snow_2021,w),
            "snow_change_2021_vs_2020":weighted_mean(g.snow_2021-g.snow_2020,w),
        })
    return pd.DataFrame(rows)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config",default="config/pipeline.yaml")
    ap.add_argument("--selected",default=
        "data/processed/historical_prc_2021_pixel_snow_test/"
        "selected_nuts3_pilot_regions.csv")
    ap.add_argument("--output-dir",default=
        "data/processed/historical_prc_2021_august_control")
    ap.add_argument("--limit",type=int,default=0)
    ap.add_argument("--max-snow-percent",type=float,default=1.)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    self_test()
    check_control_threshold(args.max_snow_percent)
    if args.self_test:
        return
    if args.limit<0:
        raise ValueError("--limit must be nonnegative")
    config=load_config(args.config)
    chosen=Path(args.selected)
    if not chosen.is_absolute():
        chosen=PROJECT_ROOT/chosen
    ids=pd.read_csv(chosen,dtype={"NUTS_ID":str}).NUTS_ID.astype(str).tolist()
    if args.limit:
        ids=ids[:args.limit]
    if not ids:
        raise ValueError("Selected pilot NUTS3 list is empty")
    nuts=gpd.read_file(nuts_path_for_release(config,2016))
    nuts=nuts.to_crs("EPSG:4326").set_index("NUTS_ID")
    missing=set(ids)-set(nuts.index)
    if missing:
        raise ValueError(f"Selected IDs absent from NUTS2016: {sorted(missing)}")
    pairs,_=discover_inputs(configured_path(config,"viirs_raw"))
    need={("vcmslcfg",year,MONTH) for year in YEARS}
    if need-set(pairs):
        raise FileNotFoundError(f"Missing August VIIRS products: {sorted(need-set(pairs))}")
    src={year:pairs[("vcmslcfg",year,MONTH)] for year in YEARS}
    platforms={year:platform(src[year].radiance) for year in YEARS}
    if "unknown" in (platforms[2020],platforms[2021]) or \
            platforms[2020]!=platforms[2021]:
        raise ValueError(f"2020/2021 VIIRS source platform mismatch: {platforms}")
    same_2022=platforms[2022]==platforms[2021]
    snowraw=discover_mod10cm(PROJECT_ROOT/
        config["analysis"]["fixed_nuts2016_snow_viirs_test"]["mod10cm_raw_dir"])
    snowcache=configured_path(config,"processed")/"mod10cm_061_geotiff_cache"
    snowpaths={year:ensure_snow_geotiff(snowraw[(year,8)],snowcache,"EPSG:4326")
               for year in YEARS}
    out=Path(args.output_dir)
    if not out.is_absolute():
        out=PROJECT_ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    manifest={
        "date_month":"August","viirs_product":"VCMSLCFG",
        "platforms":platforms,
        "radiance_reference":"mean(2020,2022)" if same_2022 else "2020 only",
        "include_2022_viirs":same_2022,
        "snow_periods":[f"{year}-08" for year in YEARS],
        "snowfree_max_percent_each_year":args.max_snow_percent,
        "radiance_files":{y:str(src[y].radiance) for y in YEARS},
        "cf_cvg_files":{y:str(src[y].cloudfree) for y in YEARS},
        "snow_files":{y:str(snowpaths[y]) for y in YEARS},
        "note":"No August VNP43MA3 albedo: this is a snow-only control."
    }
    (out/"august_source_platforms.json").write_text(
        json.dumps(manifest,indent=2)+"\n")
    frames,coverage=[],[]
    with ExitStack() as stack:
        active=(2020,2021,2022) if same_2022 else (2020,2021)
        radiance={y:stack.enter_context(rasterio.open(src[y].radiance))
                  for y in active}
        cf={y:stack.enter_context(rasterio.open(src[y].cloudfree))
            for y in active}
        equal_grid(list(radiance.values())+list(cf.values()))
        snow={y:stack.enter_context(rasterio.open(snowpaths[y]))
              for y in YEARS}
        grid=(snow[2021].transform,snow[2021].width)
        for year in YEARS:
            if (snow[year].transform,snow[year].width)!=grid:
                raise ValueError("MOD10CM native grids vary by year")
        for idx,nuts_id in enumerate(ids,1):
            groups,c=group_region(nuts_id,nuts.loc[nuts_id].geometry,
                                  radiance,cf,snow,grid,same_2022,
                                  args.max_snow_percent)
            coverage.append(c)
            if groups is not None:
                frames.append(groups)
            print(f"August {idx}/{len(ids)} {nuts_id}: {c}",flush=True)
    cov=pd.DataFrame(coverage)
    cov.to_csv(out/"august_snowfree_valid_pixel_coverage.csv",index=False)
    if not frames:
        raise RuntimeError("No August matched regions; review coverage CSV")
    g=pd.concat(frames,ignore_index=True)
    g.to_csv(out/"august_native_snowcell_brightness_groups.csv",index=False)
    reg=summary_regions(g)
    reg.to_csv(out/"august_snowfree_region_comparison.csv",index=False)
    lines=[
        "AUGUST SNOW-FREE NEGATIVE CONTROL — 2021 VIIRS",
        f"NUTS3 pilot regions requested={len(ids)}, usable={g.NUTS_ID.nunique()}",
        f"VIIRS sensors 2020={platforms[2020]},2021={platforms[2021]},"
        f"2022={platforms[2022]}",
        f"2022 VIIRS matched platform used={same_2022}",
        f"Radiance reference={manifest['radiance_reference']}",
        f"MOD10CM snow-free: <= {args.max_snow_percent:.2f}% in ALL "
        "August 2020, August 2021 and August 2022",
        f"Snow-free groups {int(g.snowfree_all_three.sum())}/{len(g)}",
        f"Snow-free fraction of matched pixels "
        f"{cov.n_snowfree_pixels.sum()/cov.n_snow_valid.sum():.3f}",
        "No 2021 radiance or PRC correction was applied.",
        "August visible WSA albedo is NOT available in the original winter archive.",
        "August 2021 vs August 2020 is NOT a 2021 vs 2020/22 flank "
        "comparison when 2022 VIIRS is J01 or otherwise sensor mismatched.",
    ]
    for is_free,part in g.groupby("snowfree_all_three"):
        lines.append(
            f"snowfree={bool(is_free)}: cells/brightness groups {len(part)}, "
            f"mean 2021-minus-reference nW="
            f"{weighted_mean(part.radiance_delta,part.area_weight):+.5f}")
    lines.extend(["","Region-level snow-free versus snow-present contrasts:",
                  reg.to_string(index=False)])
    (out/"august_snowfree_control_summary.txt").write_text(
        "\n".join(lines)+"\n",encoding="utf-8")
    print("\n".join(lines),flush=True)
    print("Outputs:",out,flush=True)


if __name__=="__main__":
    main()
