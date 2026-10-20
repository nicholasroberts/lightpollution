#!/usr/bin/env python3
"""One annual/monthly VIIRS source bridge check before 2021 snow correction.

Compare cloud-free weighted 12-month NPP VCMSLCFG radiance (2020/2021) to
official annual EOG VNL V2 average_masked on exactly the same geographic
pixels in the 18 original NUTS2016 pilot regions. Official annual products
include outlier removal and lit-background masking; a weighted monthly mean
is a DIAGNOSTIC PROXY, not an exact re-creation.

This script makes NO snow or PRC correction and never rewrites any raster.
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
from rasterio.enums import Resampling
from rasterio.features import geometry_mask, geometry_window
from rasterio.warp import reproject

from diagnose_fixed_nuts2016_monthly_product_atlas import discover_inputs
from extract_nuts3_radiance import discover_rasters, nuts_path_for_release
from project_config import PROJECT_ROOT, configured_path, load_config

YEARS=(2020,2021)
MONTHS=tuple(range(1,13))


def source_platform(src):
    stem=Path(src).name.lower()
    pattern=re.search(r"(?<![a-z0-9])(npp|j01)(?![a-z0-9])",stem)
    return pattern.group(1) if pattern else "unknown"


def warp_window(ds, reference, win):
    if ds.count != 1:
        raise ValueError(f"Expected one-band monthly source: {ds.name}")
    if ds.crs is None or reference.crs is None:
        raise ValueError(f"Missing source CRS: {ds.name}")
    result=np.full((int(win.height),int(win.width)),np.nan,dtype=np.float32)
    kw=dict(source=rasterio.band(ds,1),destination=result,
            dst_transform=reference.window_transform(win),
            dst_crs=reference.crs,dst_nodata=np.nan,
            init_dest_nodata=True,resampling=Resampling.nearest,num_threads=1)
    if ds.nodata is not None:
        kw["src_nodata"]=ds.nodata
    reproject(**kw)
    return result.astype("float64")


def weighted_mean(a,weights):
    ok=np.isfinite(a)&np.isfinite(weights)&(weights>0)
    return float(np.average(a[ok],weights=weights[ok])) if ok.any() else float("nan")


def inspect_one_region(geom,annual,monthly,year,minimum_lit=0.2):
    try:
        win=geometry_window(annual,[geom],pad_x=0,pad_y=0)
    except rasterio.errors.WindowError:
        return {"year":year,"status":"no_overlap"}
    hh,ww=int(win.height),int(win.width)
    inside=geometry_mask([geom],out_shape=(hh,ww),
                         transform=annual.window_transform(win),invert=True)
    if inside.sum()<20:
        return {"year":year,"status":"too_few_pixels"}
    old=annual.read(1,window=win,masked=True)
    a=np.ma.filled(old.astype("float64"),np.nan)
    raw_good=np.isfinite(a)&(a>=0)&inside
    raw_lit=raw_good&(a>=minimum_lit)
    weighted_sum=np.zeros((hh,ww),dtype=np.float64)
    cf_sum=np.zeros((hh,ww),dtype=np.float64)
    winter_cf=np.zeros((hh,ww),dtype=np.float64)
    winter_sum=np.zeros((hh,ww),dtype=np.float64)
    used_month_count=np.zeros((hh,ww),dtype=np.int16)
    for month in MONTHS:
        pair=monthly[(year,month)]
        rad=warp_window(pair[0],annual,win)
        cf=warp_window(pair[1],annual,win)
        valid=(inside & np.isfinite(rad) & (rad>=0) &
               np.isfinite(cf) & (cf>0))
        weighted_sum[valid]+=rad[valid]*cf[valid]
        cf_sum[valid]+=cf[valid]
        used_month_count[valid]+=1
        if month in (1,2):
            winter_cf[valid]+=cf[valid]
            winter_sum[valid]+=rad[valid]*cf[valid]
    annual_proxy=np.divide(weighted_sum,cf_sum,
                           out=np.full_like(weighted_sum,np.nan),
                           where=cf_sum>0)
    annual_winter_fraction=np.divide(winter_cf,cf_sum,
                           out=np.full_like(winter_cf,np.nan),
                           where=cf_sum>0)
    winter_average=np.divide(winter_sum,winter_cf,
                           out=np.full_like(winter_sum,np.nan),
                           where=winter_cf>0)
    yes=raw_lit & (cf_sum>0) & np.isfinite(annual_proxy)
    if yes.sum()<30:
        return {"year":year,"status":"insufficient_comparable_lit_pixels",
                "total_inside_pixels":int(inside.sum()),
                "comparable_lit_pixels":int(yes.sum())}
    rr,cc=np.indices((hh,ww))
    center=annual.window_transform(win)
    lat=center.f+center.d*(cc+.5)+center.e*(rr+.5)
    area=np.cos(np.deg2rad(lat))
    residual=annual_proxy[yes]-a[yes]
    yarea=area[yes]
    ae=np.abs(residual)
    aa=a[yes]
    rmse=float(np.sqrt(np.average(residual**2,weights=yarea)))
    result={
        "year":year,"status":"ok",
        "total_inside_pixels":int(inside.sum()),
        "valid_annual_pixels":int(raw_good.sum()),
        "official_annual_zero_pixels":int((raw_good&(a==0)).sum()),
        "monthly_annual_common_pixels":int((raw_good&(cf_sum>0)).sum()),
        "comparable_lit_pixels":int(yes.sum()),
        "lit_area_fraction_comparable":float(yes.sum()/raw_lit.sum()) if raw_lit.any() else np.nan,
        "official_annual_lit_mean":weighted_mean(aa,yarea),
        "monthly_weighted_proxy_mean":weighted_mean(annual_proxy[yes],yarea),
        "mean_proxy_minus_annual":weighted_mean(residual,yarea),
        "rmse_proxy_vs_annual":rmse,
        "median_abs_proxy_minus_annual":float(np.median(ae)),
        "median_abs_pct_error_proxy":float(np.median(100*ae/aa)),
        "median_months_with_cloudfree_data":float(np.median(used_month_count[yes])),
        "mean_jan_feb_cf_weight_fraction":weighted_mean(annual_winter_fraction[yes],yarea),
        "median_jan_feb_cf_weight_fraction":float(np.median(annual_winter_fraction[yes])),
        "mean_jan_feb_radiance":weighted_mean(winter_average[yes],yarea),
        "fraction_lit_pixels_all_12months":float(np.mean(used_month_count[yes]==12)),
        "fraction_lit_pixels_gt_20pct_proxy_error":float(np.mean(
            ae > 0.2 * aa)),
        "fraction_winter_cf_missing":float(np.mean(winter_cf[yes]==0)),
    }
    return result


def self_test():
    month_rad=np.array([8.,12.,3.])
    count=np.array([2.,2.,6.])
    observed=(month_rad*count).sum()/count.sum()
    assert observed==5.8
    expected_winter=(count[:2].sum()/count.sum())
    assert abs(expected_winter-.4)<1e-9
    g=2.
    decrement=np.sum(month_rad[:2]*count[:2]*(1.-1./g))/sum(count)
    assert abs(decrement-2.0)<1e-9
    assert decrement>=0
    assert source_platform("SVDNB_npp_20210101_vcmslcfg.tif")=="npp"
    assert source_platform("SVDNB_j01_20220101_vcmslcfg.tif")=="j01"
    print("Annual/monthly bridge arithmetic self-test: PASS",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config",default="config/pipeline.yaml")
    ap.add_argument("--selected",default=
        "data/processed/historical_prc_2021_pixel_snow_test/"
        "selected_nuts3_pilot_regions.csv")
    ap.add_argument("--output-dir",default=
        "data/processed/prc_2021_annual_monthly_bridge")
    ap.add_argument("--limit",type=int,default=6)
    ap.add_argument("--minimum-lit",type=float,default=.2)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    self_test()
    if args.self_test:
        return
    if args.limit<0:
        raise ValueError("Limit must be >=0 (0=all selected NUTS3)")
    if args.minimum_lit<=0:
        raise ValueError("Lit threshold must be >0")
    cfg=load_config(args.config)
    sel=Path(args.selected)
    if not sel.is_absolute():
        sel=PROJECT_ROOT/sel
    ids=pd.read_csv(sel,dtype={"NUTS_ID":str}).NUTS_ID.astype(str).tolist()
    if args.limit:
        ids=ids[:args.limit]
    if not ids:
        raise ValueError("No regions to test")
    nuts=gpd.read_file(nuts_path_for_release(cfg,2016)).to_crs("EPSG:4326")
    nuts=nuts.set_index("NUTS_ID")
    if unknown:=set(ids)-set(nuts.index):
        raise ValueError(f"NUTS2016 IDs missing: {unknown}")
    pairs,_=discover_inputs(configured_path(cfg,"viirs_raw"))
    require={("vcmslcfg",year,month) for year in YEARS for month in MONTHS}
    if absent:=require-set(pairs):
        raise FileNotFoundError(f"Missing monthly data: {sorted(absent)}")
    platforms={}
    for year in YEARS:
        used={source_platform(pairs["vcmslcfg",year,month].radiance)
              for month in MONTHS}
        platforms[year]=sorted(used)
        if used!={"npp"}:
            raise ValueError(f"Cannot bridge NPP annual with monthly platforms {year}: {used}")
    annualpaths=discover_rasters(cfg,list(YEARS))
    out=Path(args.output_dir)
    if not out.is_absolute():
        out=PROJECT_ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    manifest={
        "annual_rasters":{y:str(annualpaths[y]) for y in YEARS},
        "monthly_platforms":platforms,
        "monthly_pairs":{f"{y}-{m:02d}":{
            "radiance":str(pairs["vcmslcfg",y,m].radiance),
            "cloudfree":str(pairs["vcmslcfg",y,m].cloudfree)}
            for y in YEARS for m in MONTHS},
        "source_note":"Annual VNL V2 average_masked includes median-based outlier removal and background mask; weighted monthly annual proxy is not exact.",
        "correction_applied":False,
    }
    (out/"source_audit.json").write_text(json.dumps(manifest,indent=2)+"\n")
    rows=[]
    with ExitStack() as stack:
        ann={y:stack.enter_context(rasterio.open(annualpaths[y])) for y in YEARS}
        sources={}
        for y in YEARS:
            for m in MONTHS:
                pair=pairs["vcmslcfg",y,m]
                sources[y,m]=(
                    stack.enter_context(rasterio.open(pair.radiance)),
                    stack.enter_context(rasterio.open(pair.cloudfree)))
        for year in YEARS:
            if ann[year].crs is None or ann[year].crs.to_epsg()!=4326:
                raise ValueError(f"Unexpected annual raster CRS: {ann[year].crs}")
            for j,nuts_id in enumerate(ids,1):
                local={ (y,m):sources[y,m] for y in [year] for m in MONTHS }
                result=inspect_one_region(nuts.loc[nuts_id].geometry,ann[year],
                                          local,year,args.minimum_lit)
                result["NUTS_ID"]=nuts_id
                rows.append(result)
                print(f"{year} region {j}/{len(ids)} {nuts_id}: {result}",
                      flush=True)
    d=pd.DataFrame(rows)
    d.to_csv(out/"annual_monthly_bridge_by_region.csv",index=False)
    valid=d[d.status.eq("ok")]
    lines=[
        "2020/2021 ANNUAL VNL V2 versus MONTHLY VCMSLCFG BRIDGE (NO CORRECTION)",
        f"Regions requested={len(ids)}; region-years valid={len(valid)}/{len(d)}",
        f"Monthly platform: NPP only in 2020 and 2021",
        "Monthly proxy=sum(monthly_rad * monthly_cloudfree)/sum(monthly_cloudfree)",
        "Official annual average_masked additionally removes monthly outliers",
        "and applies a lit-background mask. Equality is NOT expected.",
        f"Lit-area threshold: official raw annual >= {args.minimum_lit} nW",
        "",
    ]
    for year in YEARS:
        yr=valid[valid.year.eq(year)]
        if yr.empty:continue
        w=yr.comparable_lit_pixels.to_numpy(float)
        lines.extend([
            f"YEAR {year}: regions={len(yr)}; evaluated lit pixels={int(w.sum())}",
            f"  annual official mean (pixel count weighted) = "
            f"{weighted_mean(yr.official_annual_lit_mean,w):.6f}",
            f"  monthly CF-weighted proxy mean = "
            f"{weighted_mean(yr.monthly_weighted_proxy_mean,w):.6f}",
            f"  mean monthly proxy - official annual = "
            f"{weighted_mean(yr.mean_proxy_minus_annual,w):+.6f}",
            f"  RMSE proxy - annual (approx pooled RMS) = "
            f"{np.sqrt(np.average(yr.rmse_proxy_vs_annual**2,weights=w)):.6f}",
            f"  mean Jan+Feb cloudfree fraction of annual = "
            f"{weighted_mean(yr.mean_jan_feb_cf_weight_fraction,w):.4f}",
            f"  fraction pixels with >20% monthly proxy error = "
            f"{weighted_mean(yr.fraction_lit_pixels_gt_20pct_proxy_error,w):.4f}",
        ])
    lines.extend([
        "",
        "Important: The 2021 candidate winter correction would use",
        " sum_{Jan,Feb}(rad_month*cf_month*(1-1/G_month))/sum_{12 months}(cf_month)",
        "only after establishing correspondence to the official annual VNL",
        "and accounting for VNL median-based outlier rejection + lit mask.",
        "This script does NOT apply, fit or approve the snow/albedo gain.",
    ])
    (out/"annual_monthly_bridge_summary.txt").write_text("\n".join(lines)+"\n")
    print("\n".join(lines),flush=True)
    print("Outputs:",out)


if __name__=="__main__":
    main()
