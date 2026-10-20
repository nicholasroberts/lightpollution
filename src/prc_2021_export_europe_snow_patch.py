#!/usr/bin/env python3
"""Export the frozen-gain 2021 annual spatial correction as a native VIIRS patch.

A 4-band Europe-window raster at EXACT native EOG 2021 grid resolution:
  1 raw corrected EOG 2021 annual radiance, 2 raw radiance removed,
  3 annual radiance fraction removed, 4 winter QA months (0..2).
Outside the patch the original global 2021 VNL raster remains authoritative.
DO NOT treat the Europe-only TIFF as an entire annual input to global NUTS.
No 2017+ 0.235 pixel offset is applied here; PRC applies it downstream.
"""
from __future__ import annotations

import argparse
import json
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

from historical_prc_pixel_nuts3 import load_inputs
from project_config import PROJECT_ROOT, load_config
from prc_2021_annual_monthly_bridge import warp_window
from prc_2021_spatial_annual_snow_trial import (
    source_paths, calculate_decrement, YEARS, MONTHS, WINTER,
    BETA_SNOW, BETA_ALBEDO,
)
from test_prc_2021_pixel_snow_albedo import reproject_to_region

DEFAULT_BOUNDS=(-25.,34.,45.,72.5)


def rect_from_bounds(ds,bounds):
    left,bottom,right,top=bounds
    if not left<right or not bottom<top:
        raise ValueError("Bounds require west south east north")
    if ds.crs is None or ds.crs.to_epsg()!=4326:
        raise ValueError("2021 annual VIIRS must be EPSG:4326")
    r0,c0=ds.index(left,top)
    r1,c1=ds.index(right,bottom)
    c0=max(0,min(c0,c1))
    c1=min(ds.width,max(ds.index(left,top)[1],ds.index(right,bottom)[1])+1)
    r0=max(0,min(r0,r1))
    r1=min(ds.height,max(ds.index(left,top)[0],ds.index(right,bottom)[0])+1)
    if c1<=c0 or r1<=r0:
        raise ValueError("Bounds outside annual raster")
    return Window(c0,r0,c1-c0,r1-r0)


def blocks(full,tile):
    for y in range(0,int(full.height),tile):
        for x in range(0,int(full.width),tile):
            w=min(tile,int(full.width)-x)
            h=min(tile,int(full.height)-y)
            yield (Window(int(full.col_off)+x,int(full.row_off)+y,w,h),
                   Window(x,y,w,h))


def calc_tile(ref,w,monthly,snow,albedo):
    raw=np.ma.filled(ref.read(1,window=w,masked=True).astype(float),np.nan)
    valid=np.isfinite(raw)&(raw>=0)
    proxy_numer=np.zeros(raw.shape,float)
    cf_total=np.zeros(raw.shape,float)
    excess_numer=np.zeros(raw.shape,float)
    qa_months=np.zeros(raw.shape,np.uint8)
    for month in MONTHS:
        rm=warp_window(monthly[month][0],ref,w)
        cf=warp_window(monthly[month][1],ref,w)
        use=np.isfinite(rm)&(rm>=0)&np.isfinite(cf)&(cf>0)
        proxy_numer[use]+=rm[use]*cf[use]
        cf_total[use]+=cf[use]
        if month not in WINTER:
            continue
        snowm={y:reproject_to_region(snow[y,month],ref,w) for y in YEARS}
        albm={y:reproject_to_region(albedo[y,month],ref,w) for y in YEARS}
        qa=use.copy()
        for y in YEARS:
            qa&=np.isfinite(snowm[y])&(snowm[y]>=0)&(snowm[y]<=100)
            qa&=np.isfinite(albm[y])&(albm[y]>=0)&(albm[y]<=1)
        qa_months[qa]+=1
        s=np.maximum(snowm[2021]-(snowm[2020]+snowm[2022])/2,0)/100
        a=np.maximum(albm[2021]-(albm[2020]+albm[2022])/2,0)
        gain=1+BETA_SNOW*s+BETA_ALBEDO*a
        good=qa&np.isfinite(gain)&(gain>=1)
        excess_numer[good]+=(rm[good]*cf[good]*(1-1/gain[good]))
    proxy=np.divide(proxy_numer,cf_total,
                    out=np.full(raw.shape,np.nan),where=cf_total>0)
    excess=np.divide(excess_numer,cf_total,
                     out=np.zeros(raw.shape,float),where=cf_total>0)
    corrected,fraction=calculate_decrement(raw,proxy,excess)
    corrected[~valid]=np.nan
    fraction[~valid]=np.nan
    removed=raw-corrected
    qa_float=qa_months.astype("float32")
    qa_float[~valid]=np.nan
    return (corrected,removed,fraction,qa_float)


def test():
    ann=np.array([0,10,20,np.nan],float)
    corrected,frac=calculate_decrement(
        ann,np.array([2,8,25,1]),np.array([0,2,5,2]))
    assert np.allclose(corrected[:3],[0,7.5,16])
    assert np.isnan(corrected[3])
    assert np.allclose(frac[:3],[0,.25,.20])
    from rasterio.io import MemoryFile
    from rasterio.transform import from_origin
    with MemoryFile() as m:
        with m.open(driver="GTiff",width=200,height=200,count=1,
                    dtype="float32",crs="EPSG:4326",
                    transform=from_origin(-30,80,.5,.5)) as ds:
            win=rect_from_bounds(ds,(-25,34,45,72.5))
            area=sum(int(a.width*a.height) for a,b in blocks(win,64))
            assert area==int(win.width*win.height)
    print("Native-grid patch and correction self-test: PASS",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--bounds",type=float,nargs=4,default=DEFAULT_BOUNDS)
    ap.add_argument("--tile",type=int,default=512)
    ap.add_argument("--output-dir",default=
                    "data/processed/prc_2021_spatial_europe_candidate")
    args=ap.parse_args()
    test()
    if args.self_test:return
    if args.tile<128 or args.tile>2048:
        raise ValueError("Tile size must be 128..2048")
    cfg=load_config()
    transfer,_,_,_=load_inputs(cfg,PROJECT_ROOT/"config/prc_transfer.yaml")
    if transfer.get("status")!="locked" or \
            transfer.get("fit",{}).get("included_prc_observations")!=17:
        raise RuntimeError("Use the locked 17-site PRC function")
    pairs,annual,snowpath,albpath=source_paths(cfg)
    out=Path(args.output_dir)
    if not out.is_absolute():out=PROJECT_ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    tif=out/"vnl_2021_snow_albedo_europe_patch.tif"
    if tif.exists():
        raise FileExistsError(f"Do not overwrite a prior correction: {tif}")
    temp=tif.with_suffix(".incomplete.tif")
    if temp.exists():
        raise FileExistsError(f"Remove stale incomplete TIFF after inspecting: {temp}")
    counters={k:0 for k in
        ["n_valid_lit","n_adjusted","n_reduction_gt25","n_reduction_gt50",
         "n_any_winter_QA","n_both_winter_QA"]}
    counters.update(area_lit=0.,area_fraction_removed_sum=0.,
                    annual_before_area_sum=0.,annual_after_area_sum=0.)
    with ExitStack() as stack:
        annual_ds=stack.enter_context(rasterio.open(annual[2021]))
        rect=rect_from_bounds(annual_ds,args.bounds)
        m={mon:(
            stack.enter_context(rasterio.open(pairs["vcmslcfg",2021,mon].radiance)),
            stack.enter_context(rasterio.open(pairs["vcmslcfg",2021,mon].cloudfree)))
           for mon in MONTHS}
        snow={key:stack.enter_context(rasterio.open(p)) for key,p in snowpath.items()}
        alb={key:stack.enter_context(rasterio.open(p)) for key,p in albpath.items()}
        if any(ds.descriptions[0]!="wsa_qa01_mean" for ds in alb.values()):
            raise ValueError("Expected QA 0/1 VNP43MA3 WSA visible mean band")
        profile=annual_ds.profile.copy()
        profile.update(driver="GTiff",dtype="float32",count=4,
            width=int(rect.width),height=int(rect.height),
            transform=annual_ds.window_transform(rect),
            nodata=np.nan,tiled=True,blockxsize=512,blockysize=512,
            compress="deflate",predictor=3,BIGTIFF="IF_SAFER")
        total=((int(rect.width)+args.tile-1)//args.tile)*(
               (int(rect.height)+args.tile-1)//args.tile)
        with rasterio.open(temp,"w",**profile) as dst:
            for n,desc in enumerate([
                "2021 native annual raw corrected for Jan-Feb surface gain",
                "Raw annual 2021 radiance removed",
                "Fraction of annual radiance removed",
                "Winter months of valid all-three-year MOD10CM and VNP43MA3 QA"
            ],1):
                dst.set_band_description(n,desc)
            for i,(source,local) in enumerate(blocks(rect,args.tile),1):
                adjusted,removed,frac,qa=calc_tile(
                    annual_ds,source,m,snow,alb)
                for n,a in enumerate((adjusted,removed,frac,qa),1):
                    dst.write(a.astype(np.float32),n,window=local)
                good=np.isfinite(adjusted)&(adjusted>0)
                original=adjusted+removed
                lit=np.isfinite(original)&(original>0)
                window_transform=annual_ds.window_transform(source)
                yy=np.arange(int(source.height),dtype=float)
                lat=window_transform.f+window_transform.e*(yy+.5)
                weight=np.broadcast_to(np.cos(np.deg2rad(lat))[:,None],lit.shape)
                w=weight[lit]
                counters["n_valid_lit"]+=int(lit.sum())
                counters["n_adjusted"]+=int(np.sum(lit&(frac>0)))
                counters["n_reduction_gt25"]+=int(np.sum(lit&(frac>.25)))
                counters["n_reduction_gt50"]+=int(np.sum(lit&(frac>.50)))
                counters["n_any_winter_QA"]+=int(np.sum(lit&(qa>=1)))
                counters["n_both_winter_QA"]+=int(np.sum(lit&(qa>=2)))
                counters["area_lit"]+=float(w.sum())
                counters["area_fraction_removed_sum"]+=float((frac[lit]*w).sum())
                counters["annual_before_area_sum"]+=float((original[lit]*w).sum())
                counters["annual_after_area_sum"]+=float((adjusted[lit]*w).sum())
                if i%20==0 or i==total:
                    print(f"Europe 2021 snow patch {i}/{total} tiles",flush=True)
    temp.replace(tif)
    valid=counters["n_valid_lit"]
    area=counters["area_lit"]
    counters.update({
        "status":"candidate_not_locked",
        "output_tif":str(tif),
        "annual_vnl_original":str(annual[2021]),
        "global_annual_overlay_required":True,
        "global_window_col":int(rect.col_off),
        "global_window_row":int(rect.row_off),
        "width":int(rect.width),"height":int(rect.height),
        "bounds":list(map(float,args.bounds)),
        "gain_snow":BETA_SNOW,"gain_albedo":BETA_ALBEDO,
        "fraction_lit_adjusted":counters["n_adjusted"]/valid if valid else None,
        "fraction_lit_gt25":counters["n_reduction_gt25"]/valid if valid else None,
        "fraction_lit_gt50":counters["n_reduction_gt50"]/valid if valid else None,
        "fraction_lit_any_QA":counters["n_any_winter_QA"]/valid if valid else None,
        "fraction_lit_both_QA":counters["n_both_winter_QA"]/valid if valid else None,
        "area_mean_fraction_removed":counters["area_fraction_removed_sum"]/area if area else None,
        "area_mean_raw_before":counters["annual_before_area_sum"]/area if area else None,
        "area_mean_raw_after":counters["annual_after_area_sum"]/area if area else None,
        "note":"Raw VIIRS 0.235 background subtraction is NOT applied; annual VNL outlier months not reconstructed",
    })
    path=out/"europe_patch_summary.json"
    path.write_text(json.dumps(counters,indent=2)+"\n")
    print(json.dumps(counters,indent=2),flush=True)
    print("Wrote native-grid Europe-only candidate patch:",tif,flush=True)


if __name__=="__main__":
    main()
