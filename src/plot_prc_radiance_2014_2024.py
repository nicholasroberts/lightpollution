#!/usr/bin/env python3
"""Compare NUTS2016 2014->2024 radiance and pixel-transformed locked17 PRC.

Four maps: absolute and percentage endpoint changes, radiance left, PRC right.
One country-faceted PDF scatter (x Δradiance, y ΔPRC); merged CSV.
Radiance per-region is area-weighted raw annual EOG pixel mean with the SAME
0.235 post-2017 calibration as PRC, applied per pixel and clipped to zero.
PRC is pixel-transformed BEFORE NUTS aggregation, not f(mean radiance).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource
from extract_nuts3_radiance import discover_rasters
from historical_prc_pixel_nuts3 import load_inputs,load_nuts
from project_config import PROJECT_ROOT,load_config


class CorrectedRadianceSource(RasterioRasterSource):
    def __init__(self,ds,offset):
        super().__init__(ds,1,name="radiance")
        self.offset=offset
    def read_window(self,x0,y0,nx,ny):
        a=np.ma.asarray(super().read_window(x0,y0,nx,ny))
        values=np.asarray(a.data,float)
        invalid=np.ma.getmaskarray(a)|~np.isfinite(values)|(values<0)
        return np.ma.array(np.maximum(values-self.offset,0),mask=invalid)


def extract_rad(config,nuts,year,offset,from_year):
    src=discover_rasters(config,[year])[year]
    with rasterio.open(src) as ds:
        out=exact_extract(CorrectedRadianceSource(
            ds,offset if year>=from_year else 0),nuts,
            ["mean(coverage_weight=area_spherical_m2)",
             "count(coverage_weight=area_spherical_km2)"],
            include_cols=["NUTS_ID"],output="pandas",progress=True)
    for target,variants in (
        ("radiance",["radiance_mean","mean"]),
        ("area_rad",["radiance_count","count"])):
        found=[col for col in variants if col in out.columns]
        if len(found)!=1:raise ValueError(f"{target}: unexpected fields {list(out)}")
        out=out.rename(columns={found[0]:target})
    out["year"]=year
    return out


def prepare(prc,rad,minimum_baseline):
    if prc.duplicated(["NUTS_ID","year"]).any():raise ValueError("Duplicate PRC site-year")
    if rad.duplicated(["NUTS_ID","year"]).any():raise ValueError("Duplicate radiance site-year")
    p=prc[prc.year.isin([2014,2024])][["NUTS_ID","year","CNTR_CODE",
                  "NUTS_NAME","prc_mean","valid_area_km2","coverage_stable"]]
    d=p.merge(rad,on=["NUTS_ID","year"],validate="one_to_one")
    if d.groupby("NUTS_ID").year.nunique().min()!=2:
        raise ValueError("Missing NUTS endpoint in radiance or PRC")
    a=d[d.year.eq(2014)].set_index("NUTS_ID")
    b=d[d.year.eq(2024)].set_index("NUTS_ID")
    out=a[["CNTR_CODE","NUTS_NAME","coverage_stable",
           "radiance","prc_mean","valid_area_km2","area_rad"]].join(
        b[["radiance","prc_mean","valid_area_km2","area_rad"]],
        lsuffix="_2014",rsuffix="_2024",validate="one_to_one")
    out=out.rename(columns={
        "CNTR_CODE":"country","NUTS_NAME":"name",
        "radiance_2014":"radiance_2014","radiance_2024":"radiance_2024",
        "prc_mean_2014":"prc_2014","prc_mean_2024":"prc_2024"})
    out["coverage_stable"]=out.coverage_stable_2014 & out.coverage_stable_2024
    out["radiance_change_abs"]=out.radiance_2024-out.radiance_2014
    out["prc_change_abs"]=out.prc_2024-out.prc_2014
    out["radiance_change_pct"]=np.where(
        out.radiance_2014>=minimum_baseline,
        100*out.radiance_change_abs/out.radiance_2014,np.nan)
    out["prc_change_pct"]=100*out.prc_change_abs/out.prc_2014
    out["rad_pct_available"]=out.radiance_2014>=minimum_baseline
    out.index.name="NUTS_ID"
    return out.reset_index()


def map_pair(data,nuts,columns,labels,out,title):
    g=nuts.merge(data[["NUTS_ID",*columns]],how="left",on="NUTS_ID")
    fig,axs=plt.subplots(1,2,figsize=(16,9),layout="constrained")
    for ax,key,label in zip(axs,columns,labels):
        vals=np.abs(data[key].to_numpy(float))
        vals=vals[np.isfinite(vals)]
        vmax=float(np.quantile(vals,.98)) if len(vals) else 1.
        vmax=max(vmax,1.e-10)
        g.plot(ax=ax,column=key,cmap="RdBu",norm=TwoSlopeNorm(vcenter=0,
               vmin=-vmax,vmax=vmax),linewidth=0,
               missing_kwds={"color":"#d5d5d5"},
               legend=True,legend_kwds={"label":label,"shrink":.62})
        ax.set(xlim=(-12,36),ylim=(35,73),title=label,
               xlabel="Longitude",ylabel="Latitude")
    fig.suptitle(title+"\nSame NUTS2016 boundaries; independent symmetric 98th-percentile metric limits")
    fig.savefig(out,bbox_inches="tight")
    plt.close(fig)


def facets(d,out):
    allcountries=sorted(d.country.dropna().unique())
    with PdfPages(out) as pdf:
        for start in range(0,len(allcountries),12):
            section=allcountries[start:start+12]
            fig,axes=plt.subplots(3,4,figsize=(15,11),
                                  sharex=False,sharey=False,
                                  layout="constrained")
            for ax,country in zip(axes.flat,section):
                x=d[(d.country==country)&d.coverage_stable].dropna(
                    subset=["radiance_change_abs","prc_change_abs"])
                ax.scatter(x.radiance_change_abs,x.prc_change_abs,
                           alpha=.6,s=9,rasterized=True)
                ax.axhline(0,lw=.5,color="gray")
                ax.axvline(0,lw=.5,color="gray")
                ax.set_title(f"{country} (n={len(x)})")
                ax.set_xlabel("Δ calibrated radiance")
                ax.set_ylabel("Δ pixel-derived PRC")
                ax.grid(alpha=.15)
            for ax in list(axes.flat)[len(section):]:ax.set_visible(False)
            fig.suptitle("2014–2024 NUTS3 ΔPRC versus Δradiance, by country")
            pdf.savefig(fig)
            plt.close(fig)


def self_test():
    x=np.array([.0,.235,1.,5.])
    result=np.maximum(x-.235,0)
    assert np.allclose(result,[0,0,.765,4.765])
    example=pd.DataFrame({"NUTS_ID":["X","X","Y","Y"],
        "year":[2014,2024,2014,2024],"CNTR_CODE":["DE"]*4,
        "NUTS_NAME":["sample"]*4,"prc_mean":[.4,.3,.44,.42],
        "valid_area_km2":[1.]*4,"coverage_stable":[True]*4})
    rr=pd.DataFrame({"NUTS_ID":["X","X","Y","Y"],
        "year":[2014,2024,2014,2024],"radiance":[3.,5.,0.,1.],
        "area_rad":[1.]*4})
    d=prepare(example,rr,.2)
    assert d.loc[d.NUTS_ID.eq("X"),"radiance_change_abs"].iloc[0]==2
    assert np.isnan(d.loc[d.NUTS_ID.eq("Y"),"radiance_change_pct"].iloc[0])
    print("End-point alignment and low-baseline percentage QA: PASS",flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--prc-input",default="data/processed/historical_prc_2014_2024_locked17/prc_nuts3_annual_cumulative.csv")
    p.add_argument("--out",default="figures/prc_radiance_comparison_locked17")
    p.add_argument("--min-baseline-radiance",type=float,default=.2)
    args=p.parse_args()
    self_test()
    if args.self_test:return
    cfg=load_config()
    tr,params,offset,first=load_inputs(cfg,
                        PROJECT_ROOT/"config/prc_transfer.yaml")
    if tr.get("status")!="locked" or tr["fit"]["included_prc_observations"]!=17:
        raise RuntimeError("Expected locked 17-site PRC transfer")
    prc=pd.read_csv(PROJECT_ROOT/args.prc_input,dtype={"NUTS_ID":str})
    nuts=load_nuts(cfg)
    rad=pd.concat([extract_rad(cfg,nuts,y,offset,first) for y in (2014,2024)],
                  ignore_index=True)
    d=prepare(prc,rad,args.min_baseline_radiance)
    out=PROJECT_ROOT/args.out
    out.mkdir(parents=True,exist_ok=True)
    d.to_csv(out/"nuts3_radiance_prc_cumulative_change_2014_2024.csv",index=False)
    map_pair(d,nuts,["radiance_change_abs","prc_change_abs"],
             ["Δ radiance (nW cm⁻² sr⁻¹)","Δ polarization-receptor contrast"],
             out/"compare_cumulative_absolute_radiance_vs_prc_2014_2024.pdf",
             "Absolute change, 2014 → 2024")
    map_pair(d,nuts,["radiance_change_pct","prc_change_pct"],
             ["Δ radiance (%)","Δ polarization-receptor contrast (%)"],
             out/"compare_cumulative_percent_radiance_vs_prc_2014_2024.pdf",
             f"Percentage change, 2014 → 2024; radiance baseline ≥{args.min_baseline_radiance:g}")
    facets(d,out/"scatter_prc_vs_radiance_change_by_country_2014_2024.pdf")
    print("Four comparison outputs written under",out,flush=True)


if __name__=="__main__":
    main()
