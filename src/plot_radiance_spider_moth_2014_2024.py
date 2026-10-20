#!/usr/bin/env python3
"""Three-column 2014–2024 maps: radiance, locked spider PRC, locked moth PRC.

Reuses previously computed calibrated NUTS3 radiance endpoint means (not another
full VIIRS raster run); adds separate native-pixel moth PRC extraction results.
Produces absolute maps, percentage maps, country-faceted scatter (both taxa),
and joined source/derived NUTS3 table, with identical fixed NUTS2016 cohort.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from historical_prc_pixel_nuts3 import load_nuts, load_inputs
from project_config import PROJECT_ROOT, load_config


def path(s):
    p=Path(s)
    return p if p.is_absolute() else PROJECT_ROOT/p


def endpoints(moth):
    if moth.duplicated(["NUTS_ID","year"]).any():
        raise ValueError("Duplicate moth NUTS3/year")
    frames={}
    for year in (2014,2024):
        a=moth[moth.year.eq(year)].copy()
        if a.NUTS_ID.duplicated().any():
            raise ValueError("Duplicate site")
        frames[year]=a.set_index("NUTS_ID")
    a,b=frames[2014],frames[2024]
    if len(a)!=len(b) or len(a)!=1522 or set(a.index)!=set(b.index):
        raise ValueError("Expected identical 1522 NUTS2016 NUTS3 at both endpoints")
    result=a[["prc_mean","valid_area_km2","coverage_stable"]].join(
        b[["prc_mean","valid_area_km2","coverage_stable"]],
        lsuffix="_2014",rsuffix="_2024",validate="one_to_one")
    result=result.rename(columns={"prc_mean_2014":"moth_prc_2014",
                                  "prc_mean_2024":"moth_prc_2024"})
    result["moth_prc_change_abs"]=result.moth_prc_2024-result.moth_prc_2014
    result["moth_prc_change_pct"]=100*result.moth_prc_change_abs/result.moth_prc_2014
    result["moth_coverage_stable"]=result.coverage_stable_2014 & result.coverage_stable_2024
    result["moth_area_2014"]=result.valid_area_km2_2014
    result["moth_area_2024"]=result.valid_area_km2_2024
    return result.reset_index()[[
        "NUTS_ID","moth_prc_2014","moth_prc_2024",
        "moth_prc_change_abs","moth_prc_change_pct",
        "moth_coverage_stable","moth_area_2014","moth_area_2024"]]


def join_data(spider_and_radiance,moth_frame):
    req={"NUTS_ID","country","radiance_2014","radiance_2024",
         "radiance_change_abs","radiance_change_pct","prc_2014","prc_2024",
         "prc_change_abs","prc_change_pct","coverage_stable",
         "area_rad_2014","area_rad_2024","valid_area_km2_2014","valid_area_km2_2024"}
    if missing:=req-set(spider_and_radiance.columns):
        raise ValueError(f"Existing spider+radiance table missing fields: {missing}")
    if spider_and_radiance.NUTS_ID.duplicated().any():
        raise ValueError("Duplicate spider region")
    d=spider_and_radiance.merge(moth_frame,on="NUTS_ID",
                                validate="one_to_one",how="outer",indicator=True)
    if not d._merge.eq("both").all() or len(d)!=1522:
        raise ValueError("Moth and spider NUTS3 cohorts differ")
    for value in ("2014","2024"):
        x=d[f"moth_area_{value}"].to_numpy(float)
        y=d[f"valid_area_km2_{value}"].to_numpy(float)
        if not np.allclose(x,y,rtol=1e-6,atol=.05,equal_nan=True):
            raise ValueError(f"Moth and spider valid-area cohort differs in {value}")
    d["both_taxa_coverage_stable"]=d.coverage_stable & d.moth_coverage_stable
    d=d.rename(columns={
        "prc_2014":"spider_prc_2014","prc_2024":"spider_prc_2024",
        "prc_change_abs":"spider_prc_change_abs",
        "prc_change_pct":"spider_prc_change_pct"})
    return d.drop(columns=["_merge"])


def map_three(d,nuts,columns,labels,path_out,title):
    g=nuts.merge(d[["NUTS_ID",*columns]],on="NUTS_ID",how="left")
    fig,axs=plt.subplots(1,3,figsize=(21,8),layout="constrained")
    for ax,key,label in zip(axs,columns,labels):
        v=np.abs(d[key].to_numpy(float))
        v=v[np.isfinite(v)]
        vmax=max(float(np.quantile(v,.98)) if len(v) else 0,1e-10)
        g.plot(ax=ax,column=key,cmap="RdBu",
               norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax),
               linewidth=0,legend=True,legend_kwds={"shrink":.64,"label":label},
               missing_kwds={"color":"#d4d4d4"})
        ax.set(title=label,xlim=(-12,36),ylim=(35,73),
               xlabel="Longitude",ylabel="Latitude")
    fig.suptitle(title+"\nFixed NUTS2016; symmetric 98th-percentile limits separately per metric")
    fig.savefig(path_out,bbox_inches="tight")
    plt.close(fig)


def country_scatter(d,out):
    countries=sorted(d.country.dropna().unique())
    with PdfPages(out) as pdf:
        for ix in range(0,len(countries),8):
            subset=countries[ix:ix+8]
            fig,axes=plt.subplots(4,2,figsize=(13,15),layout="constrained")
            for ax,country in zip(axes.flat,subset):
                selected=d[d.country.eq(country)&d.both_taxa_coverage_stable]
                ax.scatter(selected.radiance_change_abs,selected.spider_prc_change_abs,
                           s=11,alpha=.65,label="Spider PRC")
                ax.scatter(selected.radiance_change_abs,selected.moth_prc_change_abs,
                           s=11,alpha=.65,marker="x",label="Moth PRC")
                ax.axhline(0,lw=.65,color="gray")
                ax.axvline(0,lw=.65,color="gray")
                ax.set_title(f"{country} (n={len(selected)})")
                ax.set_xlabel("Δ calibrated radiance (nW cm⁻² sr⁻¹)")
                ax.set_ylabel("Δ PRC")
                ax.grid(alpha=.15)
                ax.legend(fontsize=7,loc="best")
            for ax in list(axes.flat)[len(subset):]:
                ax.set_visible(False)
            fig.suptitle("2014–2024 NUTS3 ΔPRC vs Δradiance, separate moth and spider responses")
            pdf.savefig(fig)
            plt.close(fig)


def self_test():
    # Test joins/renaming and conservation of two species' independent values.
    n=1522
    rows=[]
    for yr in (2014,2024):
        for idx in range(n):
            rows.append({"year":yr,"NUTS_ID":f"X{idx:04d}",
                         "prc_mean":.6 if yr==2014 else .5,
                         "valid_area_km2":1.,"coverage_stable":True})
    frame=endpoints(pd.DataFrame(rows))
    assert len(frame)==n
    assert np.allclose(frame.moth_prc_change_abs,-.1)
    assert frame.moth_coverage_stable.all()
    print("Moth endpoint pairing: PASS",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    parser.add_argument("--moth-panel",default=
        "data/processed/historical_prc_2014_2024_locked_moth16/prc_nuts3_annual_cumulative.csv")
    parser.add_argument("--spider-radiance",default=
        "figures/prc_radiance_comparison_locked17/nuts3_radiance_prc_cumulative_change_2014_2024.csv")
    parser.add_argument("--output-dir",default="figures/prc_radiance_spider_moth_2014_2024")
    args=parser.parse_args()
    self_test()
    if args.self_test:return
    cfg=load_config()
    _,s,so,sstart=load_inputs(cfg,PROJECT_ROOT/"config/prc_transfer.yaml")
    _,m,mo,mstart=load_inputs(cfg,PROJECT_ROOT/"config/prc_transfer_moth.yaml")
    assert (so,sstart)==(mo,mstart)
    assert not np.allclose(s,m), "Moth and spider transfer functions must differ"
    d=join_data(pd.read_csv(path(args.spider_radiance),dtype={"NUTS_ID":str}),
                endpoints(pd.read_csv(path(args.moth_panel),dtype={"NUTS_ID":str})))
    nuts=load_nuts(cfg)
    if len(nuts)!=1522:
        raise ValueError("NUTS2016 cohort changed")
    out=path(args.output_dir)
    out.mkdir(parents=True,exist_ok=True)
    d.to_csv(out/"nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv",index=False)
    map_three(d,nuts,["radiance_change_abs","spider_prc_change_abs",
                     "moth_prc_change_abs"],
              ["Δ radiance (nW cm⁻² sr⁻¹)","Δ spider PRC","Δ moth PRC"],
              out/"compare_absolute_radiance_spider_moth_2014_2024.pdf",
              "Absolute cumulative change 2014 → 2024")
    map_three(d,nuts,["radiance_change_pct","spider_prc_change_pct",
                     "moth_prc_change_pct"],
              ["Δ radiance (%)","Δ spider PRC (%)","Δ moth PRC (%)"],
              out/"compare_percent_radiance_spider_moth_2014_2024.pdf",
              "Percentage cumulative change 2014 → 2024 (radiance threshold ≥0.2)")
    country_scatter(d,out/"scatter_spider_moth_prc_vs_radiance_by_country_2014_2024.pdf")
    print(f"Completed {len(d)} NUTS3: triple maps, country scatter and combined CSV in {out}")


if __name__=="__main__":
    main()
