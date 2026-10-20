#!/usr/bin/env python3
"""Locked-17 PRC 2014–2024 with a reversible spatial 2021 VIIRS overlay.

No original raster files overwritten. Band 1 of the Europe patch replaces
2021 annual pixels ONLY on their exact native grid; elsewhere return the
original 2021 annual raster. PRC transformation occurs before NUTS3 averaging.
The new panel is a CANDIDATE, not an approved snow-free annual VNL product.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource
from historical_prc_pixel_nuts3 import (
    PRCRasterSource,prc_of_radiance,load_inputs,load_nuts,
    build_panel,europe_table,
)
from extract_nuts3_radiance import discover_rasters,validate_raster
from project_config import PROJECT_ROOT,load_config
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


class OverlaySource(PRCRasterSource):
    def __init__(self,original,patch,xoff,yoff,params,offset):
        super().__init__(original,offset,params)
        self.patch=patch
        self.xoff=int(xoff)
        self.yoff=int(yoff)
    def read_window(self,x0,y0,nx,ny):
        x0,y0,nx,ny=map(int,(x0,y0,nx,ny))
        raw=np.ma.array(RasterioRasterSource.read_window(
            self,x0,y0,nx,ny),copy=True)
        left=max(x0,self.xoff)
        right=min(x0+nx,self.xoff+self.patch.width)
        top=max(y0,self.yoff)
        bottom=min(y0+ny,self.yoff+self.patch.height)
        if right>left and bottom>top:
            crop=self.patch.read(1,window=Window(
                left-self.xoff,top-self.yoff,right-left,bottom-top),
                masked=True)
            target=raw[top-y0:bottom-y0,left-x0:right-x0]
            source_data=np.asarray(crop.data,float)
            old_data=np.asarray(target.data,float)
            apply=(~np.ma.getmaskarray(crop)&
                   ~np.ma.getmaskarray(target)&
                   np.isfinite(source_data)&(source_data>=0)&
                   np.isfinite(old_data)&(old_data>0))
            if np.any(source_data[apply]>old_data[apply]+1e-4):
                raise ValueError("Spatial patch attempts to brighten a VIIRS pixel")
            target.data[apply]=source_data[apply]
        return prc_of_radiance(raw,self.offset,*self.params)


def verify(original,patch,manifest):
    if patch.count!=4 or original.count!=1:
        raise ValueError("Expected annual band1 and 4-band snow patch")
    if original.crs!=patch.crs or original.crs.to_epsg()!=4326:
        raise ValueError("Annual/patch CRS mismatch")
    row,col=original.index(patch.bounds.left,patch.bounds.top)
    if (col,row)!=(manifest["global_window_col"],manifest["global_window_row"]):
        raise ValueError("Patch origin disagrees with annual-grid manifest")
    if (patch.width,patch.height)!=(manifest["width"],manifest["height"]):
        raise ValueError("Patch dimensions disagree with manifest")
    if not original.window_transform(
        Window(col,row,patch.width,patch.height)).almost_equals(patch.transform):
        raise ValueError("Patch resolution/grid misalignment")
    if Path(original.name).resolve()!=Path(manifest["annual_vnl_original"]).resolve():
        raise ValueError("Patch was derived from a different annual source")
    return col,row


def test():
    from rasterio.io import MemoryFile
    from rasterio.transform import from_origin
    p=(.145315902063,.456128356207,.103711560708)
    with MemoryFile() as x,MemoryFile() as y:
        with x.open(driver="GTiff",height=4,width=4,count=1,
                    dtype="float32",crs="EPSG:4326",
                    transform=from_origin(0,4,1,1)) as old:
            a=np.arange(16,dtype="float32").reshape(4,4)
            old.write(a,1)
            with y.open(driver="GTiff",height=2,width=2,count=4,
                        dtype="float32",nodata=np.nan,crs="EPSG:4326",
                        transform=from_origin(1,3,1,1)) as patch:
                patch.write(np.array([[1,2],[5,6]],dtype="float32"),1)
                layer=OverlaySource(old,patch,1,1,p,.235)
                v=layer.read_window(0,0,4,4)
                baseline=prc_of_radiance(a,.235,*p)
                assert np.isclose(v[1,1],prc_of_radiance(
                    np.array([1]),.235,*p)[0])
                assert np.isclose(v[0,0],baseline[0,0])
                assert np.isclose(v[3,3],baseline[3,3])
    print("2021 patch pixel alignment/PRC overlay self-test: PASS")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--baseline",default=
        "data/processed/historical_prc_2014_2024_locked17/prc_nuts3_annual_cumulative.csv")
    ap.add_argument("--patch",default=
        "data/processed/prc_2021_spatial_europe_candidate/vnl_2021_snow_albedo_europe_patch.tif")
    ap.add_argument("--patch-summary",default=
        "data/processed/prc_2021_spatial_europe_candidate/europe_patch_summary.json")
    ap.add_argument("--out",default="data/processed/prc_2021_locked17_combined")
    args=ap.parse_args()
    test()
    if args.self_test:return
    cfg=load_config()
    transfer,params,offset,first=load_inputs(cfg,PROJECT_ROOT/"config/prc_transfer.yaml")
    if transfer.get("status")!="locked" or transfer["fit"]["included_prc_observations"]!=17:
        raise ValueError("Expected locked 17-site PRC function")
    baseline=pd.read_csv(PROJECT_ROOT/args.baseline,dtype={"NUTS_ID":str})
    if sorted(baseline.year.unique())!=list(range(2014,2025)):
        raise ValueError("Baseline must include 2014–2024")
    nuts=load_nuts(cfg)
    raw=discover_rasters(cfg,[2021])[2021]
    validate_raster(raw,cfg["nuts"]["crs"])
    summary=json.loads((PROJECT_ROOT/args.patch_summary).read_text())
    with rasterio.open(raw) as original,rasterio.open(PROJECT_ROOT/args.patch) as patch:
        col,row=verify(original,patch,summary)
        source=OverlaySource(original,patch,col,row,params,
                             offset if 2021>=first else 0)
        cols=[q for q in ("NUTS_ID","CNTR_CODE","NUTS_NAME","NAME_LATN")
              if q in nuts.columns]
        values=exact_extract(source,nuts,
            ["mean(coverage_weight=area_spherical_m2)",
             "count(coverage_weight=area_spherical_km2)"],
            include_cols=cols,output="pandas",progress=True)
    rename={}
    for new,candidates in (("prc_mean",("prc_mean","mean")),
                           ("valid_area_km2",("prc_count","count"))):
        found=[x for x in candidates if x in values]
        if len(found)!=1:raise ValueError(f"Unexpected exactextract fields: {list(values)}")
        rename[found[0]]=new
    values=values.rename(columns=rename)
    prev=baseline[baseline.year.eq(2021)]
    compare=prev.merge(values[["NUTS_ID","prc_mean","valid_area_km2"]],
        on="NUTS_ID",how="outer",validate="one_to_one",
        suffixes=("_original","_candidate"),indicator=True)
    if not compare._merge.eq("both").all():
        raise ValueError("Different 2021 NUTS3 cohort after spatial correction")
    if not np.allclose(compare.valid_area_km2_original,
                       compare.valid_area_km2_candidate,rtol=1e-7,atol=.05):
        raise ValueError("Native patch altered VIIRS spatial-validity support")
    if (compare.prc_mean_candidate < compare.prc_mean_original-5e-6).any():
        raise ValueError("2021 snow correction unexpectedly reduced PRC")
    updated=baseline.copy()
    changed=updated.year.eq(2021)
    changed_values=compare.set_index("NUTS_ID").prc_mean_candidate
    updated.loc[changed,"prc_mean"]=updated.loc[changed,"NUTS_ID"].map(changed_values)
    updated.loc[changed,"viirs_file"]="2021 original annual VNL + spatial snow/albedo overlay"
    updated=build_panel(updated,list(range(2014,2025)),.95)
    output=PROJECT_ROOT/args.out
    output.mkdir(parents=True,exist_ok=True)
    updated.to_csv(output/"prc_nuts3_annual_cumulative_candidate.csv",index=False)
    europe_old=europe_table(baseline)
    europe_new=europe_table(updated)
    europe_new.to_csv(output/"prc_europe_annual_cumulative_candidate.csv",index=False)
    c=europe_old[["year","area_weighted_prc"]].merge(
        europe_new[["year","area_weighted_prc"]],
        on="year",suffixes=("_original","_candidate"),validate="one_to_one")
    c["correction"]=c.area_weighted_prc_candidate-c.area_weighted_prc_original
    c.to_csv(output/"europe_2021_original_vs_candidate.csv",index=False)
    compare["prc_change_snow"]=compare.prc_mean_candidate-compare.prc_mean_original
    compare.to_csv(output/"nuts3_2021_spatial_adjustment.csv",index=False)
    fig,ax=plt.subplots(figsize=(10,5))
    ax.plot(c.year,c.area_weighted_prc_original,marker="o",label="Original 2021")
    ax.plot(c.year,c.area_weighted_prc_candidate,marker="o",label="Spatially adjusted 2021")
    ax.set(title="Europe PRC: locked17 response and spatial 2021 snow/albedo candidate",
           xlabel="Year",ylabel="Area-weighted polarization-receptor contrast")
    ax.grid(alpha=.2);ax.legend()
    fig.savefig(output/"europe_original_vs_candidate_PR C.pdf".replace(" ",""),bbox_inches="tight")
    plt.close(fig)
    t=c.set_index("year")
    old=float(t.loc[2021,"area_weighted_prc_original"]-
        (t.loc[2020,"area_weighted_prc_original"]+
         t.loc[2022,"area_weighted_prc_original"])/2)
    new=float(t.loc[2021,"area_weighted_prc_candidate"]-
        (t.loc[2020,"area_weighted_prc_candidate"]+
         t.loc[2022,"area_weighted_prc_candidate"])/2)
    verdict={"status":"candidate_not_locked","regions":len(compare),
        "2021_europe_prc_original":float(t.loc[2021,"area_weighted_prc_original"]),
        "2021_europe_prc_candidate":float(t.loc[2021,"area_weighted_prc_candidate"]),
        "anomaly_vs_2020_22_original":old,
        "anomaly_vs_2020_22_candidate":new,
        "fraction_anomaly_removed":1-abs(new)/abs(old) if old else None,
        "outside_2021_values_unchanged":True,
        "official_annual_outlier_processing_not_reconstructed":True}
    (output/"prc_2021_spatial_overlay_summary.json").write_text(
        json.dumps(verdict,indent=2)+"\n")
    print(json.dumps(verdict,indent=2))


if __name__=="__main__":
    main()
