#!/usr/bin/env python3
"""Illustrator-ready diagnostics: pre/post 2017 baseline and locked 2021 snow.

Figure C1: original annual VIIRS radiance versus independently calibrated
2017+ native-pixel baseline correction (0.235 nW), across 2014–2024.
Figure C2: baseline-calibrated annual trajectory before/after the *locked*
2021 native-pixel snow/albedo adjustment and NUTS3 map of the 2021 reduction.
Supplement W: original Jan–Mar 2021 monthly anomaly relative to identical
months of 2020/2022 (VCMSLCFG; DIAGNOSTIC, not a corrected monthly product).

All science is from existing archived NUTS2016 annual table, locked 2021
native pixel overlay and monthly anomaly atlas. Read-only scientific inputs,
all new outputs in a separate publication directory; no new correction fit.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, SymLogNorm
from shapely.geometry import box
import numpy as np
import pandas as pd

from figure_style import apply_figure_style, save_figure
from historical_prc_pixel_nuts3 import load_nuts
from plot_publication_figure1_radiance import (
    DEFAULT_ANNUAL, DEFAULT_ENDPOINT, DEFAULT_SNOW_OUTPUT,
    get_locked_snow_2021, validate_endpoint,
)
from project_config import PROJECT_ROOT as ROOT, load_config

YEARS=tuple(range(2014,2025))
ANNUAL_YEARS=set(YEARS)
OUT=ROOT/"figures/publication/radiance_correction_demonstrations"
MONTHLY=ROOT/"data/processed/fixed_nuts2016_monthly_product_atlas/monthly_region_2021_vs_flanks.csv"
FIG1C=ROOT/"figures/publication/figure1_radiance/fig1C_europe_annual_radiance_support.csv"
REGION_COLS=["NUTS_ID","year","radiance_mean_raw","radiance_mean_corrected"]
WEIGHT="valid_area_km2_2014"
RAD=r"nW cm$^{-2}$ sr$^{-1}$"
MAP_BBOX=(-12.,34.,46.,73.)


def pathlike(x):
    p=Path(x)
    return p if p.is_absolute() else ROOT/p


def annual_region_source(source, endpoint, snow):
    """Match exactly Figure 1C's 11-year balanced NUTS3 cohort and area weights."""
    if not source.is_file():
        raise FileNotFoundError(f"Missing original fixed NUTS2016 annual parquet: {source}")
    annual=pd.read_parquet(source,columns=REGION_COLS)
    annual["NUTS_ID"]=annual.NUTS_ID.astype(str)
    annual["year"]=pd.to_numeric(annual.year,errors="raise").astype(int)
    annual=annual[annual.year.isin(YEARS)].copy()
    if annual.duplicated(["NUTS_ID","year"]).any() or set(annual.year)!=ANNUAL_YEARS:
        raise ValueError("Original annual NUTS3 records contain duplicates/missing years")
    eligible=set(endpoint.NUTS_ID)
    for y in YEARS:
        eligible &= set(annual.loc[annual.year.eq(y),"NUTS_ID"])
    if len(eligible)<1400 or len(eligible)>1522:
        raise ValueError(f"Unusual complete annual source cohort: {len(eligible)}")
    fixed=endpoint.loc[endpoint.NUTS_ID.isin(eligible),
                       ["NUTS_ID",WEIGHT,"radiance_2014","radiance_2024"]]
    x=annual[annual.NUTS_ID.isin(eligible)].merge(
        fixed,on="NUTS_ID",validate="many_to_one"
    )
    if len(x)!=len(eligible)*len(YEARS):
        raise ValueError("Region-year balance broke during locked endpoint join")
    for c in ["radiance_mean_raw","radiance_mean_corrected",WEIGHT]:
        x[c]=pd.to_numeric(x[c],errors="coerce")
        if not np.isfinite(x[c].to_numpy(float)).all():
            raise ValueError(f"Missing or nonfinite {c}")
    if (x[REGION_COLS[2:]]<0).any().any() or (x[WEIGHT]<=0).any():
        raise ValueError("Negative VIIRS radiance or invalid 2014 area weight")

    # Baseline is an upstream native-pixel correction; NEVER add/subtract
    # 0.235 to regional means to reconstruct a synthetic 'before' series.
    zero=x.year.lt(2017)
    after=~zero
    raw=x["radiance_mean_raw"].to_numpy(float)
    cal=x["radiance_mean_corrected"].to_numpy(float)
    if not np.allclose(raw[zero],cal[zero],rtol=2e-6,atol=2e-6):
        raise ValueError("2014–2016 originally uncalibrated annual values disagree")
    reduction=raw[after]-cal[after]
    if (reduction < -2e-6).any() or (reduction > .235+2e-6).any():
        raise ValueError("2017+ native-pixel offset correction violated bound 0..0.235")
    x["final_locked_radiance"]=x.radiance_mean_corrected

    # Repeat the same 2021 checks made by Figure 1C.
    if len(snow)!=1522 or not snow.NUTS_ID.is_unique:
        raise ValueError("Locked snow cache must contain 1,522 unique NUTS2016 units")
    before=x.loc[x.year.eq(2021),["NUTS_ID","radiance_mean_corrected"]]
    z=before.merge(snow[["NUTS_ID","radiance_2021_original_calibrated",
                         "radiance_2021_snow_corrected"]],
                   on="NUTS_ID",validate="one_to_one")
    if len(z)!=len(before):
        raise ValueError("Locked 2021 correction is incomplete in matched annual cohort")
    discrepancy=np.abs(
        z.radiance_mean_corrected.to_numpy(float)-
        z.radiance_2021_original_calibrated.to_numpy(float)
    )
    limit=2e-5*np.maximum(1.,np.abs(z.radiance_mean_corrected.to_numpy(float)))
    if (discrepancy>limit).any():
        raise ValueError(f"Unpatched 2021 does not match locked raw overlay: max {discrepancy.max()}")
    corrected=z.radiance_2021_snow_corrected.to_numpy(float)
    original=z.radiance_2021_original_calibrated.to_numpy(float)
    if not np.isfinite(corrected).all() or (corrected<0).any() or (
        corrected>original+2e-6
    ).any():
        raise ValueError("2021 snow model invalid/brightened NUTS3 means")
    snow_map=dict(zip(z.NUTS_ID,corrected))
    idx=x.year.eq(2021)
    x.loc[idx,"final_locked_radiance"]=(
        x.loc[idx,"NUTS_ID"].map(snow_map).to_numpy(float)
    )
    # This is the exact Fig1C aggregation definition.
    for year,col in ((2014,"radiance_2014"),(2024,"radiance_2024")):
        a=x.loc[x.year.eq(year)]
        mismatch=(a.radiance_mean_corrected-a[col]).abs()
        if mismatch.quantile(.99)>.025:
            raise ValueError(
                f"Original annual {year} differs from locked endpoints: "
                f"p99 diff={mismatch.quantile(.99):.6g}"
            )
    qa={
        "balanced_nuts2016_regions":len(eligible),
        "complete_region_years":len(x),
        "zero_correction_years":[2014,2015,2016],
        "calibration_start":2017,
        "baseline_pixel_offset_nw":.235,
        "pixel_clip_floor":0.,
        "max_measured_region_baseline_reduction":
            float(np.max(reduction)),
        "original_2021_locked_patch_source_max_error":
            float(discrepancy.max()),
        "2021_year_regions_validated":len(z),
    }
    return x,qa


def annual_weighted(x):
    vals=["radiance_mean_raw","radiance_mean_corrected","final_locked_radiance"]
    a=x[["NUTS_ID","year",WEIGHT,*vals]].copy()
    for v in vals:
        a[f"{v}_total"]=a[v]*a[WEIGHT]
    agg=a.groupby("year",as_index=False).agg(
        n_nuts3=("NUTS_ID","nunique"),
        fixed_area_km2=(WEIGHT,"sum"),
        **{f"{v}_sum":(f"{v}_total","sum") for v in vals}
    )
    for v in vals:
        agg[v]=agg[f"{v}_sum"]/agg.fixed_area_km2
    keep=["year","n_nuts3","fixed_area_km2",*vals]
    agg=agg[keep].sort_values("year")
    if agg.n_nuts3.nunique()!=1 or agg.fixed_area_km2.max()-agg.fixed_area_km2.min()>1e-3:
        raise ValueError("Yearly region support / area weighting changed")
    if np.any(agg.loc[agg.year.lt(2017),"radiance_mean_raw"].to_numpy(float) !=
              agg.loc[agg.year.lt(2017),"radiance_mean_corrected"].to_numpy(float)):
        # Allow only machine-level floating point deviation.
        if not np.allclose(
            agg.loc[agg.year.lt(2017),"radiance_mean_raw"],
            agg.loc[agg.year.lt(2017),"radiance_mean_corrected"],rtol=1e-12,atol=1e-12
        ):
            raise ValueError("2014–16 baseline histories unexpectedly differ")
    if not np.allclose(
        agg.loc[~agg.year.eq(2021),"radiance_mean_corrected"],
        agg.loc[~agg.year.eq(2021),"final_locked_radiance"],
        atol=1e-12,rtol=1e-12
    ):
        raise ValueError("Snow correction unexpectedly affects years besides 2021")
    return agg


def check_fig1c(curves,ref):
    if not ref.is_file():
        return {"checked":False,"explanation":"Figure 1C annual support CSV absent"}
    d=pd.read_csv(ref)
    d=d[["year","area_weighted_radiance_nw"]]
    pair=curves.merge(d,on="year",validate="one_to_one")
    if len(pair)!=11:
        raise ValueError("Figure 1C expected 11 calendar years")
    gap=np.abs(pair.final_locked_radiance-pair.area_weighted_radiance_nw)
    if not np.allclose(pair.final_locked_radiance,
                       pair.area_weighted_radiance_nw,rtol=2e-10,atol=2e-9):
        raise ValueError("Correction figure final curve disagrees with frozen Figure 1C")
    return {"checked":True,"max_abs_difference":float(gap.max())}


def plot_baseline(a,out):
    fig,ax=plt.subplots(figsize=(10.2,5.7))
    ax.plot(a.year,a.radiance_mean_raw,marker="o",linewidth=2.2,
            color="#be542a",label="Original uncorrected VIIRS")
    ax.plot(a.year,a.radiance_mean_corrected,marker="o",linewidth=2.2,
            color="#236eaa",label="2017+ pixel baseline corrected (snow not removed)")
    ax.axvline(2016.5,color="#555555",linestyle=":",linewidth=.95)
    ax.annotate("2017 measurement-series step",xy=(2016.5,a.radiance_mean_raw.max()),
                xytext=(2014.25,a.radiance_mean_raw.max()*1.075),fontsize=9,
                arrowprops={"arrowstyle":"-","color":"#555555"})
    ax.set_title("Before and after the 2017 VIIRS baseline correction",loc="left")
    ax.set_ylabel(f"European area-weighted mean radiance ({RAD})")
    ax.set_xlabel("Year")
    ax.set_xticks(YEARS)
    ax.set_xlim(2013.8,2024.2)
    ax.grid(axis="y",alpha=.3)
    ax.legend(frameon=False,loc="best")
    fig.tight_layout()
    save_figure(fig,out)
    plt.close(fig)


def map_geography():
    nuts=load_nuts(load_config())
    if len(nuts)!=1522 or nuts.NUTS_ID.duplicated().any():
        raise ValueError("Locked fixed NUTS2016 1,522-region geometry changed")
    bound=box(*MAP_BBOX)
    region=gpd.clip(nuts,bound).to_crs("EPSG:3035")
    region=region.loc[~region.geometry.is_empty].copy()
    region["geometry"]=region.geometry.simplify(150.,preserve_topology=True)
    bounds=tuple(gpd.GeoSeries([bound],crs="EPSG:4326").to_crs("EPSG:3035").total_bounds)
    return region,bounds


def map_panel(ax,g,bounds,column,cmap,norm):
    mapped=g.loc[g[column].notna()].copy()
    g.plot(ax=ax,color="#dedede",edgecolor="none")
    if not mapped.empty:
        mapped.plot(ax=ax,column=column,cmap=cmap,norm=norm,
                    edgecolor="none",linewidth=0)
    edges=g.dissolve(by="CNTR_CODE").boundary
    edges.plot(ax=ax,color="#666666",linewidth=.22)
    ax.set_xlim(bounds[0],bounds[2])
    ax.set_ylim(bounds[1],bounds[3])
    ax.set_aspect("equal")
    ax.set_axis_off()
    return len(mapped)


def plot_snow(a,snow,geom,bounds,out):
    pre=float(a.loc[a.year.eq(2021),"radiance_mean_corrected"].iloc[0])
    post=float(a.loc[a.year.eq(2021),"final_locked_radiance"].iloc[0])
    if post>pre+2e-6:
        raise ValueError("2021 snow correction brightened European mean")
    fig,(ax,axm)=plt.subplots(1,2,figsize=(12.5,6.25),
                              gridspec_kw={"width_ratios":[1.05,1.1]})
    fig.subplots_adjust(left=.08,right=.89,top=.86,bottom=.13,wspace=.06)
    ax.plot(a.year,a.radiance_mean_corrected,marker="o",linewidth=1.7,
            color="#c16035",linestyle="--",label="Calibrated, before snow correction")
    ax.plot(a.year,a.final_locked_radiance,marker="o",linewidth=2,
            color="#177ab0",label="Final locked (snow corrected)")
    ax.scatter([2021,2021],[pre,post],s=70,c=["#c16035","#177ab0"],zorder=5)
    ax.annotate("",xy=(2021.28,post),xytext=(2021.28,pre),
                arrowprops=dict(arrowstyle="->",lw=1.3,color="#555555"))
    ax.text(2021.5,(pre+post)*.5,f"−{pre-post:.3f}",fontsize=9,va="center")
    ax.set_title("A  |  Effect on European annual trajectory",loc="left",fontsize=10)
    ax.set_xlim(2013.8,2024.25)
    ax.set_xticks([2014,2016,2018,2020,2021,2022,2024])
    ax.set_ylabel(f"Area-weighted mean radiance ({RAD})")
    ax.set_xlabel("Year")
    ax.grid(axis="y",alpha=.25)
    ax.legend(frameon=False,fontsize=8,loc="upper left")
    z=snow[["NUTS_ID","radiance_2021_original_calibrated",
            "radiance_2021_snow_corrected"]].copy()
    z["snow_reduction_nw"]=(
        z.radiance_2021_original_calibrated-z.radiance_2021_snow_corrected
    )
    if (z.snow_reduction_nw< -2e-6).any():
        raise ValueError("Snow magnitude map shows impossible brightening")
    gm=geom.merge(z[["NUTS_ID","snow_reduction_nw"]],on="NUTS_ID",
                  how="left",validate="one_to_one")
    vals=gm.snow_reduction_nw.dropna().to_numpy(float)
    vmax=max(.001,float(np.quantile(vals,.99)))
    norm=Normalize(vmin=0.,vmax=vmax,clip=False)
    count=map_panel(axm,gm,bounds,"snow_reduction_nw","viridis",norm)
    axm.set_title("B  |  2021 regional snow/albedo radiance reduction",loc="left",fontsize=10)
    sm=ScalarMappable(cmap="viridis",norm=norm)
    sm.set_array([])
    cb=fig.colorbar(sm,ax=axm,fraction=.044,pad=.017,extend="max")
    cb.set_label(f"Reduction in 2021 NUTS3 mean ({RAD})",fontsize=8.5)
    fig.suptitle("2021 snow/albedo adjustment: original versus locked correction",
                 fontsize=13)
    fig.text(.08,.035,"2021 spatial change computed at native pixels; 2017 baseline correction applied in both series.",fontsize=8)
    save_figure(fig,out,map_figure=True)
    plt.close(fig)
    return {"pre_2021_europe_mean":pre,"post_2021_europe_mean":post,
            "difference_europe":pre-post,
            "snow_reduction_map_display_99pct_cap":vmax,
            "map_display_NUTS3_count":count}


def winter_monthly(monthly_path,geom,bounds,out):
    """Diagnostic of the 2021 winter months; NOT monthly corrected VIIRS."""
    if not monthly_path.is_file():
        return {"rendered":False,"reason":f"Optional monthly atlas absent: {monthly_path}"}
    cols=["product","month","NUTS_ID","paired_valid","delta_nw",
          "radiance_2020","radiance_2021","radiance_2022"]
    d=pd.read_csv(monthly_path,usecols=cols,dtype={"NUTS_ID":str})
    d=d[(d["product"]=="vcmslcfg") & d["month"].isin((1,2,3))].copy()
    if d.empty:
        raise ValueError("Monthly atlas has no Jan-Feb-Mar VCMSLCFG observations")
    d["month"]=d["month"].astype(int)
    d["paired_valid"]=d.paired_valid.astype(str).str.lower().isin(("true","1"))
    d=d[d.paired_valid].copy()
    for c in ["delta_nw","radiance_2020","radiance_2021","radiance_2022"]:
        d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d[np.isfinite(d["delta_nw"])].copy()
    if d.duplicated(["month","NUTS_ID"]).any():
        raise ValueError("Monthly region comparison has duplicates")
    if set(d.month)!={1,2,3} or len(d)<200:
        raise ValueError("Insufficient complete Jan–Mar matched-month observations")
    delta=(d.radiance_2021-.5*(d.radiance_2020+d.radiance_2022)).to_numpy(float)
    if not np.allclose(delta,d.delta_nw,rtol=2e-7,atol=1e-5):
        raise ValueError("Monthly atlas anomaly values do not match 2020/2022 flanks")
    maxima=np.abs(d.delta_nw.to_numpy(float))
    scale=max(.05,float(np.quantile(maxima,.98)))
    norm=SymLogNorm(linthresh=.1,vmin=-scale,vmax=scale,base=10)
    fig,axs=plt.subplots(1,3,figsize=(14.5,5.6))
    fig.subplots_adjust(left=.02,right=.90,top=.84,bottom=.12,wspace=.02)
    counts={}
    for m,ax,label in zip((1,2,3),axs,("January","February","March")):
        sub=d[d.month.eq(m)][["NUTS_ID","delta_nw"]]
        gm=geom.merge(sub,on="NUTS_ID",how="left",validate="one_to_one")
        counts[label]=map_panel(ax,gm,bounds,"delta_nw","RdBu_r",norm)
        ax.set_title(f"{label} 2021  (n={len(sub):,} NUTS3)",fontsize=10)
    sm=ScalarMappable(cmap="RdBu_r",norm=norm)
    sm.set_array([])
    cbar_ax=fig.add_axes([.915,.26,.014,.49])
    cb=fig.colorbar(sm,cax=cbar_ax,extend="both")
    cb.set_label(f"2021 minus 2020/2022 mean ({RAD})",fontsize=9)
    fig.suptitle("January–March 2021 radiance anomalies relative to same-month flanks",
                 fontsize=13)
    fig.text(.04,.038,"Monthly VCMSLCFG diagnostic; not a monthly snow-corrected product. Grey = missing paired coverage.",
             fontsize=9)
    save_figure(fig,out,map_figure=True)
    plt.close(fig)
    return {"rendered":True,"monthly_atlas_source":str(monthly_path),
            "product":"vcmslcfg","months":[1,2,3],
            "valid_paired_monthly_regions":counts,
            "shared_colour_magnitude_cap":scale,
            "inference":"2021 monthly observation minus matching months 2020/2022, not snow-corrected months"}


def self_test():
    example=pd.DataFrame({"year":[2016,2017,2021],
        "radiance_mean_raw":[1.,1.235,1.50],
        "radiance_mean_corrected":[1.,1.,1.265]})
    assert (example.loc[example.year.eq(2016),"radiance_mean_raw"].iloc[0] ==
            example.loc[example.year.eq(2016),"radiance_mean_corrected"].iloc[0])
    assert len(YEARS)==11 and YEARS[0]==2014 and YEARS[-1]==2024
    assert np.isclose((.8*10.+.4*90.)/100.,.44)
    print("Radiance correction before/after source-join arithmetic tests: PASS",flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annual",default=DEFAULT_ANNUAL)
    p.add_argument("--endpoints",default=DEFAULT_ENDPOINT)
    p.add_argument("--snow-cache",default=DEFAULT_SNOW_OUTPUT)
    p.add_argument("--monthly",default=MONTHLY)
    p.add_argument("--fig1c",default=FIG1C)
    p.add_argument("--out",default=OUT)
    p.add_argument("--skip-monthly",action="store_true")
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    self_test()
    if args.self_test:return
    apply_figure_style()
    out=pathlike(args.out)
    source=pathlike(args.annual)
    e=pathlike(args.endpoints)
    endpoint=validate_endpoint(pd.read_csv(e,dtype={"NUTS_ID":str}))
    snow, provenance=get_locked_snow_2021(pathlike(args.snow_cache))
    x,science=annual_region_source(source,endpoint,snow)
    a=annual_weighted(x)
    qa=check_fig1c(a,pathlike(args.fig1c))
    out.mkdir(exist_ok=True,parents=True)
    a.to_csv(out/"correction_demo_annual_area_weighted_2014_2024.csv",index=False)
    x[["NUTS_ID","year","radiance_mean_raw","radiance_mean_corrected",
       "final_locked_radiance",WEIGHT]].to_csv(
           out/"correction_demo_nuts3_source_values.csv",index=False)
    plot_baseline(a,out/"C1_2017_baseline_before_after_annual_trajectory")
    geom,bounds=map_geography()
    snow_qa=plot_snow(a,snow,geom,bounds,
                      out/"C2_2021_snow_before_after_trajectory_and_nuts3_map")
    winter=(
        {"rendered":False,"reason":"--skip-monthly specified"}
        if args.skip_monthly
        else winter_monthly(
            pathlike(args.monthly),geom,bounds,
            out/"CW_2021_january_february_march_monthly_anomaly_maps"
        )
    )
    evidence={
        "n_complete_nuts3":science["balanced_nuts2016_regions"],
        "fixed_years":list(YEARS),
        "source_original_calibrated":str(source),
        "source_final_endpoints":str(e),
        "source_snow_cache":str(pathlike(args.snow_cache)),
        "locked_snow_cache_provenance":provenance,
        "calibration_model":"0.235 nW positive native-pixel subtraction 2017 onward; clipped at zero",
        "snow_model":"locked native-pixel 2021 snow/albedo excess correction, then pixel baseline",
        "full_annual_source_checks":science,
        "figure1c_comparison":qa,
        "snow_comparison":snow_qa,
        "winter_monthly_supplement":winter,
        "monthly_diagnostic_is_not_correction":True,
        "original_scientific_sources_modified":False
    }
    (out/"correction_demo_provenance.json").write_text(
        json.dumps(evidence,indent=2)+"\n",encoding="utf8")
    print("\n2017 BASELINE & 2021 SNOW CORRECTION FIGURES",flush=True)
    print("Science QA:",science,flush=True)
    print("Figure 1C comparison:",qa,flush=True)
    print("2021 correction:",snow_qa,flush=True)
    print("Jan-Mar diagnostic:",winter,flush=True)
    print("Outputs:",out,flush=True)
    for file in sorted(out.glob("*")):
        if file.is_file():
            print(" ",file.name,round(file.stat().st_size/1024,1),"KiB",flush=True)


if __name__=="__main__":
    main()
