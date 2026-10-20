#!/usr/bin/env python3
"""Trial spatial snow/albedo 2021 annual VIIRS correction, no parameter fitting.

Pilot on original 18 NUTS2016 regions. Derive cloud-free-weighted Jan-Feb
2021 snow/albedo radiance-excess fraction at each annual EOG pixel and apply
that fraction to original 2021 Annual VNL masked radiance. This preserves the
official annual product's radiance scale, lit/background zeros and nodata.
It is a PRODUCT-APPROXIMATION, not an exact undo of VNL outlier rejection.
Report before/after locked-17-site PRC and 2020/22 annual comparisons.
"""
from __future__ import annotations

import argparse
import json
from contextlib import ExitStack
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_window, geometry_mask

from diagnose_fixed_nuts2016_monthly_product_atlas import discover_inputs
from extract_nuts3_radiance import discover_rasters, nuts_path_for_release
from historical_prc_pixel_nuts3 import load_inputs, prc_of_radiance
from project_config import PROJECT_ROOT, configured_path, load_config
from test_fixed_nuts2016_snow_viirs_2021 import (
    discover_mod10cm, ensure_snow_geotiff,
)
from test_prc_2021_pixel_snow_albedo import reproject_to_region
from prc_2021_annual_monthly_bridge import source_platform, warp_window

YEARS=(2020,2021,2022)
MONTHS=tuple(range(1,13))
WINTER=(1,2)
BETA_SNOW=1.626479
BETA_ALBEDO=3.788066


def weighted_mean(data,weights):
    x=np.asarray(data,float)
    w=np.asarray(weights,float)
    good=np.isfinite(x)&np.isfinite(w)&(w>0)
    return float(np.average(x[good],weights=w[good])) if good.any() else float("nan")


def calculate_decrement(annual, proxy, winter_excess):
    """Transfer frozen winter gain fraction from monthlies to official annual."""
    eligible=np.isfinite(annual)&(annual>0)&np.isfinite(proxy)&(proxy>0)
    fraction=np.zeros(annual.shape,float)
    fraction[eligible]=np.clip(
        np.maximum(winter_excess[eligible],0)/proxy[eligible],0,1)
    corrected=np.array(annual,copy=True)
    corrected[eligible]=annual[eligible]*(1-fraction[eligible])
    if np.any(corrected[eligible]>annual[eligible]+1e-9) or \
            np.any(corrected[eligible]<-1e-9):
        raise AssertionError("Proposed radiance correction increased or inverted radiance")
    return corrected, fraction


def self_test():
    annual=np.array([0.,10.,20.,np.nan])
    proxy=np.array([3.,12.,16.,1.])
    excess=np.array([1.,3.,100.,2.])
    new,frac=calculate_decrement(annual,proxy,excess)
    assert np.allclose(new[:3],[0,7.5,0])
    assert np.allclose(frac[:3],[0,.25,1])
    assert np.isnan(new[3])
    assert np.isclose(1.626479*.25+3.788066*.1,
                      .78542635)
    assert source_platform("SVDNB_npp_20210101_x.tif")=="npp"
    # Regression test: a perfectly valid VIIRS-style GeoTIFF without an
    # explicit nodata attribute must still reproject on the annual grid.
    from rasterio.io import MemoryFile
    from rasterio.transform import from_origin
    from rasterio.windows import Window
    src_values=np.array([[1.5,0.0],[4.25,10.0]],dtype="float32")
    with MemoryFile() as memory:
        with memory.open(driver="GTiff",height=2,width=2,count=1,
                         dtype="float32",crs="EPSG:4326",
                         transform=from_origin(-2,2,1,1)) as ds:
            ds.write(src_values,1)
            assert ds.nodata is None
            regridded=warp_window(ds,ds,Window(0,0,2,2))
            assert np.allclose(regridded,src_values,equal_nan=True),regridded
    print("Spatial annual winter-gain arithmetic, no-nodata VIIRS warp and no-brightening: PASS",flush=True)


def source_paths(config):
    pairs,_=discover_inputs(configured_path(config,"viirs_raw"))
    needed={("vcmslcfg",2021,m) for m in MONTHS}
    if needed-set(pairs):
        raise FileNotFoundError(f"Missing monthly 2021 VCMSLCFG: {sorted(needed-set(pairs))}")
    for m in MONTHS:
        name=pairs["vcmslcfg",2021,m].radiance
        if source_platform(name)!="npp":
            raise ValueError(f"2021 monthly must be NPP: {name}")
    annual=discover_rasters(config,list(YEARS))
    snow_raw_dir=PROJECT_ROOT / config["analysis"]["fixed_nuts2016_snow_viirs_test"]["mod10cm_raw_dir"]
    snow_raw=discover_mod10cm(snow_raw_dir)
    snow_dir=configured_path(config,"processed")/"mod10cm_061_geotiff_cache"
    snow={(y,m):ensure_snow_geotiff(snow_raw[y,m],snow_dir,"EPSG:4326")
          for y in YEARS for m in WINTER}
    alb_dir=configured_path(config,"processed")/"vnp43ma3_002_monthly_composite_cache"
    alb={(y,m):alb_dir/f"vnp43ma3_002_{y}_{m:02d}_monthly.tif"
         for y in YEARS for m in WINTER}
    missing=[str(p) for p in alb.values() if not p.exists()]
    if missing:
        raise FileNotFoundError("QA VNP43MA3 winter albedo not available: "+", ".join(missing))
    return pairs,annual,snow,alb


def one_region(nuts_id,geom,ds_ann,monthly,snow,albedo,params,offset,output,
               raster_exports=False):
    ref=ds_ann[2021]
    try:
        win=geometry_window(ref,[geom])
    except rasterio.errors.WindowError:
        return {"NUTS_ID":nuts_id,"status":"outside_annual_raster"}
    transform=ref.window_transform(win)
    inside=geometry_mask([geom],out_shape=(int(win.height),int(win.width)),
                         transform=transform,invert=True)
    old_raw={}
    for y in YEARS:
        a=ds_ann[y].read(1,window=win,masked=True)
        old_raw[y]=np.ma.filled(a.astype(float),np.nan)
    if not any(inside.flat):
        return {"NUTS_ID":nuts_id,"status":"no_pixels"}
    a2021=old_raw[2021]
    base_valid=inside&np.isfinite(a2021)&(a2021>=0)
    lit=base_valid&(a2021>0)
    proxy_sum=np.zeros(inside.shape,float)
    count_sum=np.zeros(inside.shape,float)
    winter_excess_sum=np.zeros(inside.shape,float)
    winter_count=np.zeros(inside.shape,float)
    snow_affected=np.zeros(inside.shape,bool)
    qa_supported=np.zeros(inside.shape,bool)
    # NUTS2016 is always queried on the same fixed annual pixel grid.
    for m in MONTHS:
        p=monthly[m]
        # Monthly VIIRS .avg_rade9h / .cf_cvg GeoTIFFs may have nodata=None.
        # Use the same valid-pixel, CRS-aware warp verified by the completed
        # annual/monthly source bridge; retain the QA-specific helper below
        # for the snow and albedo rasters, which carry explicit nodata.
        rm=warp_window(p[0],ref,win)
        cm=warp_window(p[1],ref,win)
        valid=(inside & np.isfinite(rm)&(rm>=0)&
               np.isfinite(cm)&(cm>0))
        count_sum[valid]+=cm[valid]
        proxy_sum[valid]+=(rm[valid]*cm[valid])
        if m not in WINTER:
            continue
        winter_count[valid]+=cm[valid]
        sval={}
        aval={}
        for y in YEARS:
            sval[y]=reproject_to_region(snow[(y,m)],ref,win)
            aval[y]=reproject_to_region(albedo[(y,m)],ref,win)
        qa=valid.copy()
        for y in YEARS:
            qa &= np.isfinite(sval[y])&(sval[y]>=0)&(sval[y]<=100)
            qa &= np.isfinite(aval[y])&(aval[y]>=0)&(aval[y]<=1)
        qa_supported|=qa
        sn=np.maximum(sval[2021]-.5*(sval[2020]+sval[2022]),0)/100
        al=np.maximum(aval[2021]-.5*(aval[2020]+aval[2022]),0)
        gain=1+BETA_SNOW*sn+BETA_ALBEDO*al
        snow_affected|=qa & ((sn>0)|(al>0))
        # No albedo extrapolation at pixels lacking all-three-year QA.
        valid_gain=qa&(gain>=1)&np.isfinite(gain)
        winter_excess_sum[valid_gain]+=(
            rm[valid_gain]*cm[valid_gain]*
            (1-1/gain[valid_gain]))
    proxy=np.divide(proxy_sum,count_sum,
        out=np.full_like(proxy_sum,np.nan),where=count_sum>0)
    # numerator/denominator are both in units radiance (since proxy_sum=R*CF)
    excess=np.divide(winter_excess_sum,count_sum,
        out=np.zeros_like(proxy_sum),where=count_sum>0)
    corrected,fraction=calculate_decrement(a2021,proxy,excess)
    # Guarantee no values outside original annual lit mask are changed.
    corrected[~lit]=a2021[~lit]
    before_prc=np.asarray(prc_of_radiance(a2021,offset,*params))
    after_prc=np.asarray(prc_of_radiance(corrected,offset,*params))
    flank_prc={}
    for y in (2020,2022):
        flank_prc[y]=np.asarray(prc_of_radiance(old_raw[y],offset,*params))
    # In all years the annual 2017+ offset acts only inside PRC transform.
    rows,cols=np.indices(inside.shape)
    lat=transform.f+transform.d*(cols+.5)+transform.e*(rows+.5)
    area=np.cos(np.deg2rad(lat))
    valid_out=lit&np.isfinite(before_prc)&np.isfinite(after_prc)
    full=base_valid&np.isfinite(before_prc)&np.isfinite(after_prc)
    valid_control=(full&np.isfinite(flank_prc[2020])&
                   np.isfinite(flank_prc[2022]))
    flank=.5*(flank_prc[2020]+flank_prc[2022])
    status="ok" if full.sum()>20 else "too_few_valid_pixels"
    result={
        "NUTS_ID":nuts_id,"status":status,
        "n_inside":int(inside.sum()),"n_valid_annual":int(full.sum()),
        "n_lit_annual":int(valid_out.sum()),
        "fraction_lit_with_both_winter_qa":float(np.mean(
            qa_supported[valid_out])) if valid_out.any() else np.nan,
        "fraction_lit_with_positive_surface_anomaly":float(np.mean(
            snow_affected[valid_out])) if valid_out.any() else np.nan,
        "fraction_lit_actually_reduced":float(np.mean(
            fraction[valid_out]>1e-10)) if valid_out.any() else np.nan,
        "mean_fraction_annual_radiance_removed":weighted_mean(fraction[valid_out],area[valid_out]),
        "max_fraction_annual_radiance_removed":float(np.max(fraction[valid_out])) if valid_out.any() else np.nan,
        "monthly_proxy_radiance_mean":weighted_mean(proxy[valid_out],area[valid_out]),
        "official_annual_radiance_original_mean":weighted_mean(a2021[valid_out],area[valid_out]),
        "official_annual_radiance_corrected_mean":weighted_mean(corrected[valid_out],area[valid_out]),
        "prc_2021_before":weighted_mean(before_prc[full],area[full]),
        "prc_2021_after":weighted_mean(after_prc[full],area[full]),
        "prc_flanks_2020_2022":weighted_mean(flank[valid_control],area[valid_control]),
        "mean_2021_prc_correction":weighted_mean(after_prc[full]-before_prc[full],area[full]),
        "mean_jan_feb_cf_share":weighted_mean(np.divide(winter_count,count_sum,
            out=np.zeros_like(count_sum),where=count_sum>0)[valid_out],area[valid_out]),
        "mean_proxy_annual_difference":weighted_mean((proxy-a2021)[valid_out],area[valid_out])
    }
    if raster_exports:
        def band(x):
            y=np.asarray(x,dtype="float32").copy()
            y[~inside]=np.nan
            return y
        profile=ref.profile.copy()
        profile.update(driver="GTiff",height=int(win.height),width=int(win.width),
            transform=transform,count=5,dtype="float32",nodata=np.nan,
            tiled=False,compress="deflate",predictor=3)
        folder=output/"pilot_region_pixels"
        folder.mkdir(parents=True,exist_ok=True)
        rp=folder/f"{nuts_id}_2021_snow_albedo_trial.tif"
        with rasterio.open(rp,"w",**profile) as dest:
            bands=[a2021,corrected,a2021-corrected,before_prc,after_prc]
            labels=["Original 2021 raw VIIRS","Snow-adjusted 2021 raw VIIRS",
                    "Estimated annual winter-excess removed",
                    "Original locked17 PRC","Trial snow-adjusted locked17 PRC"]
            for i,(val,label) in enumerate(zip(bands,labels),1):
                dest.write(band(val),i)
                dest.set_band_description(i,label)
        result["regional_pixel_tif"]=str(rp.relative_to(PROJECT_ROOT))
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selected",default=
        "data/processed/historical_prc_2021_pixel_snow_test/"
        "selected_nuts3_pilot_regions.csv")
    ap.add_argument("--output-dir",default=
        "data/processed/prc_2021_spatial_annual_snow_trial")
    ap.add_argument("--limit",type=int,default=6,
                    help="Selected NUTS3 subset; 0=all 18")
    ap.add_argument("--export-region-rasters",action="store_true")
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    self_test()
    if args.self_test:return
    config=load_config()
    transfer,params,offset,from_year=load_inputs(
        config,PROJECT_ROOT/"config/prc_transfer.yaml")
    if transfer.get("status")!="locked" or \
            int(transfer["fit"]["included_prc_observations"])!=17:
        raise RuntimeError("Spatial trial requires locked 17-site PRC calibration")
    if from_year>2020:
        raise RuntimeError("Annual baseline offset year mismatched")
    p=Path(args.selected)
    if not p.is_absolute():p=PROJECT_ROOT/p
    ids=pd.read_csv(p,dtype={"NUTS_ID":str}).NUTS_ID.astype(str).tolist()
    if args.limit<0:raise ValueError("limit must be >=0")
    if args.limit:ids=ids[:args.limit]
    nuts=gpd.read_file(nuts_path_for_release(config,2016))
    nuts=nuts.to_crs("EPSG:4326").set_index("NUTS_ID")
    if missing:=set(ids)-set(nuts.index):
        raise ValueError(f"Unknown NUTS2016: {missing}")
    pairs,annual_paths,snow_paths,alb_paths=source_paths(config)
    out=Path(args.output_dir)
    if not out.is_absolute():out=PROJECT_ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    provenance={
        "frozen_gain_coefficients":{"snow":BETA_SNOW,"albedo":BETA_ALBEDO},
        "frozen_prc_parameters":{"yf":params[0],"y0":params[1],"k":params[2]},
        "monthly_2021_product":"NPP VCMSLCFG V1",
        "annual_2021_product":"EOG VNL V2 average_masked",
        "algorithm":"official_annual_raw * (1 - winter_excess_monthly_raw / annual_monthly_proxy_raw)",
        "winter_excess":"sum(Jan/Feb cloudfree * R_m*(1-1/G_m)) / sum(Jan-Dec cloudfree)",
        "source_annual":{y:str(annual_paths[y]) for y in YEARS},
        "source_monthly":{m:{"radiance":str(pairs["vcmslcfg",2021,m].radiance),
                             "cloudfree":str(pairs["vcmslcfg",2021,m].cloudfree)}
                          for m in MONTHS},
        "source_snow":{f"{y}-{m:02d}":str(snow_paths[y,m])
                       for y in YEARS for m in WINTER},
        "source_albedo":{f"{y}-{m:02d}":str(alb_paths[y,m])
                         for y in YEARS for m in WINTER},
        "annual_vnl_outlier_limitation":"Annual V2 rejected outlier months. The ratio-based trial is a source-product approximation, not exact recompositing.",
        "missing_snow_or_albedo_QA":"Do not correct pixel for that month",
        "annual_zero_background":"preserve original pixel's zero state",
        "correction_accepted_as_final":False,
    }
    (out/"spatial_trial_source_audit.json").write_text(
        json.dumps(provenance,indent=2)+"\n")
    results=[]
    with ExitStack() as stack:
        annual={y:stack.enter_context(rasterio.open(annual_paths[y]))
                for y in YEARS}
        if any(annual[y].crs!=annual[2021].crs or
               annual[y].width!=annual[2021].width or
               annual[y].height!=annual[2021].height or
               not annual[y].transform.almost_equals(annual[2021].transform)
               for y in YEARS):
            raise ValueError("2020/21/22 annual rasters do not share a native grid")
        monthly={m:(
            stack.enter_context(rasterio.open(pairs["vcmslcfg",2021,m].radiance)),
            stack.enter_context(rasterio.open(pairs["vcmslcfg",2021,m].cloudfree)))
                 for m in MONTHS}
        snow={(y,m):stack.enter_context(rasterio.open(snow_paths[y,m]))
              for y in YEARS for m in WINTER}
        albedo={(y,m):stack.enter_context(rasterio.open(alb_paths[y,m]))
                for y in YEARS for m in WINTER}
        for ds in albedo.values():
            if ds.descriptions[0]!="wsa_qa01_mean":
                raise ValueError(f"Expected QA0/1 WSA visible albedo first band: {ds.name}")
        for n,nutsid in enumerate(ids,1):
            record=one_region(nutsid,nuts.loc[nutsid].geometry,annual,
                              monthly,snow,albedo,params,offset,out,
                              args.export_region_rasters)
            results.append(record)
            print(f"{n}/{len(ids)} {nutsid}: {record}",flush=True)
    df=pd.DataFrame(results)
    df.to_csv(out/"spatial_trial_region_diagnostics.csv",index=False)
    valid=df[df.status.eq("ok")]
    if valid.empty:
        raise RuntimeError("No valid pilot region results")
    weights=valid.n_valid_annual.to_numpy(float)
    summary=[
        "2021 SPATIAL ANNUAL SNOW/ALBEDO FRACTIONAL CORRECTION (PILOT ONLY)",
        f"Regions requested={len(ids)}, valid={len(valid)}",
        "Frozen snow+albedo G=1+1.626479*positive_snow_fraction"
        "+3.788066*positive_visible_albedo_fraction",
        "Annual correction scales official 2021 raw VIIRS pixel with fraction",
        "derived from observed Jan-Feb raw monthly VIIRS contribution.",
        "Original background zeros, nodata and all other years unchanged.",
        "Annual VNL outlier month rejection NOT explicitly reproduced.",
        "",
        f"Lit pixels with winter QA (at least one month)=% "
        f"{100*weighted_mean(valid.fraction_lit_with_both_winter_qa,weights):.2f}",
        f"Lit pixels with positive snow/albedo anomaly=% "
        f"{100*weighted_mean(valid.fraction_lit_with_positive_surface_anomaly,weights):.2f}",
        f"Lit pixels actually reduced=% "
        f"{100*weighted_mean(valid.fraction_lit_actually_reduced,weights):.2f}",
        f"Mean fraction of annual radiance removed in lit pixels=% "
        f"{100*weighted_mean(valid.mean_fraction_annual_radiance_removed,weights):.3f}",
        f"Mean locked17 PRC before={weighted_mean(valid.prc_2021_before,weights):.8f}",
        f"Mean locked17 PRC after={weighted_mean(valid.prc_2021_after,weights):.8f}",
        f"Mean change in locked17 PRC={weighted_mean(valid.mean_2021_prc_correction,weights):+.8f}",
        f"Mean of 2020/22 reference PRC={weighted_mean(valid.prc_flanks_2020_2022,weights):.8f}",
        "",
        "A smoother PRC year is a downstream check, NEVER the optimization target.",
        "This is not a validated Europe-wide corrected annual raster.",
        "Next: inspect corrected spatial maps/magnitudes, QA support, source bias",
        "and compare to 2020/22. If satisfactory, export one Europe-wide",
        "native-grid 2021 correction surface and rerun locked PRC.",
    ]
    (out/"spatial_trial_summary.txt").write_text("\n".join(summary)+"\n")
    print("\n".join(summary),flush=True)
    print("Outputs:",out)


if __name__=="__main__":
    main()
