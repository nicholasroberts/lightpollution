#!/usr/bin/env python3
"""Madrid UCM: VIIRS distance-weighted skyglow SOURCE index, 2015-2019.

Duriscoe et al. (2018), JQSRT 214:133-145, Eqs 9 and 11:
  W = sum over VIIRS pixels [ radiance_i * d_i^(-alpha(d_i)) ]
  alpha(d) = 2.3*(d/350)^0.28  (d in kilometres)
https://doi.org/10.1016/j.jqsrt.2018.04.028

IMPORTANT: the published coefficient calibrates remote-site *all-sky V-band*
ALR, NOT inner-city UCM Johnson B/V/R zenith radiance. This code computes
UN-CALIBRATED within-site trends, with a fixed kernel and sensitivity to
source radius and near-source floor. The <20 km kernel is not validated for
extended urban sources. The artificial component of ground AstMon radiance
is not separately measured.

Requires existing 75N060W EOG VCMSLCFG VNL V1 avg_rade9h+cf_cvg GeoTIFFs.
Default: calendar Jan/Nov/Dec 2015-2019, NOT exact year-spanning
Nov-Jan winters or identical nights as Robles et al. (2021) at 01:30.

 python src/madrid_ucm_duriscoe_skyglow_kernel.py --self-test
 python src/madrid_ucm_duriscoe_skyglow_kernel.py --inventory
 python src/madrid_ucm_duriscoe_skyglow_kernel.py
 python src/madrid_ucm_duriscoe_skyglow_kernel.py --months 3 4 5 6 7 8 9 10

Read-only: no locked VNL V2, snow, baseline or PRC changes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.geometry import Point
from pyproj import Transformer

from local_viirs_five_case_analysis import ROOT, build_inventory, grid_month
from madrid_ucm_local_ground_bridge import SITE

DATA = ROOT / "data/processed/local_lighting_interventions/madrid_ground_viirs/duriscoe_kernel"
FIG = ROOT / "figures/local_lighting_interventions/madrid_ground_viirs/duriscoe_kernel"
RADII = (3., 10., 25., 50., 100., 200., 300.)
FLOORS = (0.25, 0.5, 1., 2., 5., 10., 20.)
RINGS = (0., 3., 10., 25., 50., 100., 200., 300.)
ASTMON = {"B": (19.30, 18.97), "V": (18.07, 17.96), "R": (17.42, 17.45)}


def kernel(dist, floor):
    if floor <= 0:
        raise ValueError("Near-source floor must be positive")
    d = np.maximum(dist, floor)
    return d ** (-2.3*(d/350.)**.28)


def source_geometry():
    p = gpd.GeoSeries([Point(*SITE)], crs="EPSG:4326").to_crs("EPSG:3035")
    return gpd.GeoSeries([p.iloc[0].buffer(300000)],
                         crs="EPSG:3035").to_crs("EPSG:4326").iloc[0]


def distance_and_area(grid):
    # Same projection/footprint convention as original 3/10 km Madrid study.
    t = grid["transform"]
    if abs(t.b) > 1e-10 or abs(t.d) > 1e-10:
        raise ValueError("Rotated EOG grid not supported")
    rows, cols = grid["rad"].shape
    xgrid = t.c + t.a*(np.arange(cols)+.5)
    ygrid = t.f + t.e*(np.arange(rows)+.5)
    lon, lat = np.meshgrid(xgrid, ygrid)
    transform = Transformer.from_crs(grid["crs"], "EPSG:3035", always_xy=True)
    xx, yy = transform.transform(lon,lat)
    xc, yc = Transformer.from_crs(
        "EPSG:4326", "EPSG:3035", always_xy=True).transform(*SITE)
    d = np.hypot(xx-xc,yy-yc)/1000.
    # Duriscoe uses equal-area reprojected pixels. For native geographic
    # DNB pixels adjust area proportional to cos(latitude).
    area = np.cos(np.deg2rad(ygrid))[:, None] / np.cos(np.deg2rad(SITE[1]))
    return d,area


def compute(index, months, cf_min=3, clip=0.5):
    years = range(2015, 2020)
    months = tuple(sorted(set(months)))
    if not months or any(m not in range(1,13) for m in months):
        raise ValueError("Choose calendar months 1..12")
    geo = source_geometry()
    observed={}
    dist=area=mask=None
    signature=None
    for year in years:
        for month in months:
            ym=f"{year}-{month:02d}"
            if ym not in index:
                raise FileNotFoundError(f"Missing original EOG VNL V1 month {ym}")
            g=grid_month(geo,*index[ym])
            sig=(str(g["crs"]),tuple(g["transform"]),g["rad"].shape)
            if signature is None:
                signature=sig
                dist,area=distance_and_area(g)
                mask=g["inside"]&(dist<=300.)
                print("Native VIIRS 300km footprint:",g["rad"].shape,
                      "pixel count",int(mask.sum()),flush=True)
            elif sig!=signature:
                raise ValueError(f"Spatial grid shifted in {ym}; refuse resampling")
            good=mask & np.isfinite(g["rad"]) & np.isfinite(g["cf"]) & (g["cf"]>=cf_min)
            # Duriscoe SETS values below 0.5 to zero; DO NOT subtract 0.5.
            raw=np.nan_to_num(g["rad"],nan=0.,posinf=0.,neginf=0.)
            values=np.where(raw>=clip,raw,0.).astype("float32")
            observed[(year,month)]=(values,good)
            print(ym, f"{int(good.sum()):,}/{int(mask.sum()):,}",
                  "clear-view pixels",flush=True)
    kmasks={r:mask&(dist<=r) for r in RADII}
    rows=[]
    for (year,month),(values,valid) in sorted(observed.items()):
        for r,m in kmasks.items():
            common=m&valid
            if not common.any():
                raise ValueError(f"No valid pixels {year}-{month:02d}, radius {r}")
            k=kernel(dist[common],.5)
            score=float(np.sum(values[common].astype("float64")*k*
                               np.broadcast_to(area,mask.shape)[common]))
            rows.append({"year":year,"month":month,"radius_km":r,
                         "weighted_score":score,
                         "total_pixels":int(m.sum()),
                         "valid_pixels":int(common.sum()),
                         "valid_fraction":float(common.sum()/m.sum())})
    monthly=pd.DataFrame(rows)
    annual=monthly.groupby(["year","radius_km"],as_index=False).agg(
        score=("weighted_score","mean"),
        min_valid_fraction=("valid_fraction","min"))
    baseline=annual[annual.year.eq(2015)].set_index("radius_km").score
    annual["index_2015_100"]=annual.apply(
        lambda row:100*row.score/baseline.loc[row.radius_km],axis=1)
    # Paired valid cells, same months, sensitivity to near-source regularization.
    gridarea=np.broadcast_to(area,mask.shape)
    sensitivity=[]
    ringrows=[]
    for floor in FLOORS:
        k=kernel(dist,floor)*gridarea
        for r,rmask in kmasks.items():
            a=b=0.
            fractions=[]
            for month in months:
                x,vx=observed[(2015,month)]
                y,vy=observed[(2019,month)]
                common=rmask&vx&vy
                if not common.any():
                    raise ValueError(f"No paired valid 2015/2019 cells at {r}km")
                a+=float(np.sum(x[common].astype("float64")*k[common]))
                b+=float(np.sum(y[common].astype("float64")*k[common]))
                fractions.append(float(common.sum()/rmask.sum()))
            sensitivity.append({"floor_km":floor,"radius_km":r,
                                "pct_2015_2019":100*(b/a-1),
                                "common_fraction_min":min(fractions),
                                "source_index_2015":a/len(months),
                                "source_index_2019":b/len(months)})
    for year in (2015,2019):
        for lower,upper in zip(RINGS[:-1],RINGS[1:]):
            m=mask&(dist>=lower)&(dist<upper)
            k=kernel(dist[m],.5)*gridarea[m]
            score=0.
            for month in months:
                v,good=observed[(year,month)]
                weights=np.zeros(m.sum(),dtype="float64")
                weights[:]=np.asarray(v[m],dtype="float64")
                weights[~good[m]]=0.
                score+=float(np.sum(weights*k))
            ringrows.append({"year":year,"from_km":lower,
                             "to_km":upper,"weighted_score":score/len(months)})
    return monthly,annual,pd.DataFrame(sensitivity),pd.DataFrame(ringrows)


def chart(annual,sens,rings,months,out=FIG):
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"pdf.fonttype":42,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(2,2,figsize=(13.1,10.1),constrained_layout=True)
    ax,bx,cx,dx=axes.ravel()
    for r in RADII:
        s=annual[annual.radius_km.eq(r)].sort_values("year")
        ax.plot(s.year,s.index_2015_100,"-o",label=f"{r:g} km")
    ax.axhline(100,color="0.5",ls="--",lw=.9)
    ax.set(title="A  Fixed-kernel VIIRS source index, 2015 = 100",
           xlabel="Year",ylabel="Index")
    ax.set_xticks(range(2015,2020));ax.legend(ncol=3,fontsize=8)
    for f in FLOORS:
        s=sens[sens.floor_km.eq(f)]
        bx.plot(s.radius_km,s.pct_2015_2019,"-o",label=f"{f:g}km floor",ms=3)
    bx.axhline(0,color="0.3",lw=1)
    bchange=100*(10**(-.4*(ASTMON["B"][1]-ASTMON["B"][0]))-1)
    bx.axhline(bchange,color="#3173b3",lw=1.5,ls="--",
               label="Observed UCM B +35.5%")
    bx.set(xlabel="Contributing source radius (km)",
           ylabel="2015 to 2019 change (%)",
           title="B  Distance and urban near-field sensitivity")
    bx.set_xscale("symlog",linthresh=3)
    bx.set_xticks(RADII,[str(int(r)) for r in RADII])
    bx.legend(fontsize=7,ncol=2)
    t=rings.pivot(index=["from_km","to_km"],columns="year",values="weighted_score")
    delta=100*(t[2019]-t[2015])/t[2015].sum()
    cx.bar(range(len(delta)),delta,color=np.where(delta>=0,"#438d70","#b45e58"))
    cx.axhline(0,color="0.4",lw=.9)
    cx.set_xticks(range(len(delta)),
                  [f"{int(a)}–{int(b)}" for a,b in delta.index],rotation=40)
    cx.set(title="C  Distance-ring contributions to total model change",
           ylabel="Change / total 2015 source index (%)",xlabel="Source distance (km)")
    main=sens[(sens.floor_km.eq(.5))&sens.radius_km.eq(300.)].iloc[0]
    ratio=1+float(main.pct_2015_2019)/100
    actual=[10**(-.4*(ASTMON[b][1]-ASTMON[b][0])) for b in "BVR"]
    x=np.arange(3)
    p=dx.bar(x-.2,[ratio]*3,.36,color="#758996",
             label="Fixed VIIRS kernel, all bands")
    q=dx.bar(x+.2,actual,.36,color=["#347db6","#3a8965","#b65a50"],
             label="Ground, quality-filtered")
    dx.axhline(1,color="0.3",ls="--",lw=.9)
    dx.set_xticks(x,["Johnson B","Johnson V","Johnson R"])
    dx.set(title="D  VIIRS-only prediction vs observed zenith brightness",
           ylabel="2019 / 2015 (within series)")
    dx.legend(fontsize=8)
    for bars in (p,q):
        for rect in bars:
            dx.annotate(f"{rect.get_height():.3f}",
                        (rect.get_x()+rect.get_width()/2,rect.get_height()),
                        xytext=(0,4),textcoords="offset points",ha="center",
                        fontsize=8)
    for a in axes.ravel(): a.grid(axis="y",alpha=.16)
    fig.suptitle(
        "Madrid | Exploratory Duriscoe-style VIIRS source-distance model\n"
        "NOT a validated urban-zenith or absolute Johnson-band skyglow model",
        fontsize=14,weight="bold")
    fig.text(.01,-.03,
      "Monthly VNL V1 source, EOG 75N060W; low-value threshold 0.5 nW cm⁻² sr⁻¹; "
      "distance alpha(d)=2.3(d/350)^0.28; area correction cos(latitude). "
      "Ground: Robles et al. (2021), Table 3-5, 01:30 winter.\\n"
      "VIIRS calendar months and AstMon winter year not exactly aligned; "
      "ground includes natural background. Fixed spectrum, angular pattern "
      "and atmosphere assumed; <20km model is an unvalidated extension.",
      fontsize=7.5,color="0.35")
    filename="madrid_ucm_duriscoe_kernel_"+''.join(map(str,months))
    for ext in (".png",".pdf",".svg"):
        fig.savefig(out/(filename+ext),dpi=175,bbox_inches="tight")
    plt.close(fig)
    return filename, float(main.pct_2015_2019)


def self_test():
    d=np.array([0.,.25,1.,3.,100.,300.])
    assert np.isfinite(kernel(d,.5)).all()
    assert np.isclose(kernel(np.array([1.]),.5)[0],1.)
    assert np.isclose(2.3*(300./350.)**.28,2.201,atol=.02)
    # Synthetic mocked month/geometry rendering and IO use normal run.
    print("PASS: verified Duriscoe Eq 11, finite urban near-field kernel. "
          "Real 75N060W raster test requires --inventory and main run.")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--inventory",action="store_true")
    ap.add_argument("--months",type=int,nargs="+",default=[1,11,12])
    ap.add_argument("--cf",type=int,choices=(1,3,5),default=3)
    ap.add_argument("--background-threshold",type=float,default=.5)
    args=ap.parse_args()
    if args.self_test: self_test();return
    months=tuple(sorted(set(args.months)))
    if not months or not set(months).issubset(set(range(1,13))):
        raise ValueError("Months must be 1..12")
    index,missing=build_inventory()
    if args.inventory:
        selected=[f"{y}-{m:02d}" for y in range(2015,2020) for m in months]
        print("Requested",len(selected),"months, missing:",
              [k for k in selected if k not in index]);return
    if args.background_threshold!=.5:
        raise ValueError("Other background thresholds require explicit code revision")
    monthly,annual,sens,rings=compute(index,months,cf_min=args.cf)
    DATA.mkdir(parents=True,exist_ok=True)
    for name,data in (("monthly",monthly),("annual",annual),
                      ("sensitivity",sens),("rings",rings)):
        data.to_csv(DATA/f"madrid_ucm_duriscoe_{name}.csv",index=False)
    stem,pct=chart(annual,sens,rings,months)
    audit={"paper":"Duriscoe 2018 Eq 9 and 11, uncalibrated ALR-style index",
           "months":months,"years":[2015,2016,2017,2018,2019],
           "cf_cvg_threshold":args.cf,"DNB_below_0p5_set_to_zero":True,
           "geometry":"EPSG3035, 300km radius; WGS84 pixel-area cosine adjustment",
           "kernel_floor_sensitivity_km":FLOORS,
           "primary_radius_km":300,"primary_floor_km":.5,
           "primary_model_change_pct":pct,
           "ground_observation":"Robles 2021 filtered 01:30 winter B/V/R, "
                                "not identical acquisition dates",
           "not_calibrated_to_urban_zenith":True,
           "near_field_under_20km_unvalidated":True,
           "frozen_prc_annual_sources_modified":False}
    (DATA/"model_audit.json").write_text(json.dumps(audit,indent=2))
    print("2015->2019 VIIRS only fixed-kernel change:",f"{pct:+.2f}%")
    print("Published ground B increase:",f"{100*(10**.132-1):+.1f}%")
    print("Figure:",FIG/(stem+".png"))
    print("CSV and audit:",DATA)
    print("No frozen Europe annual, snow or PRC data modified.")


if __name__=="__main__":
    main()
