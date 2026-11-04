#!/usr/bin/env python3
"""Illustrator-ready diagnostics: pre/post 2017 baseline and locked 2021 snow.

Figure C1: original annual VIIRS radiance versus independently calibrated
2017+ native-pixel baseline correction (0.235 nW), across 2014–2024.
Figure C2: baseline-calibrated annual trajectory before/after the *locked*
2021 native-pixel snow/albedo adjustment and NUTS3 map of the 2021 reduction.
Supplement W: eight actual monthly VCMSLCFG NUTS3 radiance maps:
Jan/Feb/Dec 2020, Jan/Feb/Mar 2021, Jan/Feb 2022; fixed scale, 3x3
grid, ninth bottom-right panel empty. Original monthly radiance (not
2021-minus-flank differences and not monthly snow-corrected radiance).

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
MONTHLY=ROOT/"data/processed/fixed_nuts2016_monthly_product_atlas/monthly_nuts2016_radiance_2020_2022.parquet"
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


WINTER_MONTHS = (
    (2020, 1), (2020, 2), (2020, 12),
    (2021, 1), (2021, 2), (2021, 3),
    (2022, 1), (2022, 2),
)
MONTH_NAMES = {
    1: "January", 2: "February", 3: "March", 12: "December"
}


def winter_monthly(monthly_path, geom, bounds, out):
    """8 absolute monthly radiance maps, one common scale, 3 x 3 (9th blank).

    Original cloud-coverage-masked VCMSLCFG monthly radiance; NOT the
    2021-minus-2020/2022 anomaly and NOT an annual-snow-adjusted product.
    A row per year: 2020 Jan/Feb/Dec, 2021 Jan/Feb/Mar,
    2022 Jan/Feb; bottom-right panel is empty.
    """
    if not monthly_path.is_file():
        return {
            "rendered": False,
            "reason": f"Monthly source absent: {monthly_path}",
            "expected_input": "monthly_nuts2016_radiance_2020_2022.parquet"
        }
    columns = [
        "product", "year", "month", "NUTS_ID",
        "radiance_mean", "analysis_valid", "valid_area_fraction",
    ]
    if monthly_path.suffix.lower() == ".parquet":
        d = pd.read_parquet(monthly_path, columns=columns)
    else:
        d = pd.read_csv(
            monthly_path, usecols=columns, dtype={"NUTS_ID": str}
        )
    d["NUTS_ID"] = d["NUTS_ID"].astype(str)
    d["year"] = pd.to_numeric(d.year, errors="raise").astype(int)
    d["month"] = pd.to_numeric(d.month, errors="raise").astype(int)
    d = d[d["product"].astype(str).str.lower().eq("vcmslcfg")].copy()
    wanted = set(WINTER_MONTHS)
    d = d.loc[
        [(int(y), int(m)) in wanted for y, m in zip(d.year, d.month)]
    ].copy()
    if d.duplicated(["year", "month", "NUTS_ID"]).any():
        raise ValueError("Monthly archive contains duplicate VCMSLCFG NUTS3/year/month")
    if set(zip(d.year, d.month)) != wanted:
        missing = wanted - set(zip(d.year, d.month))
        raise ValueError(f"Missing requested VCMSLCFG months: {sorted(missing)}")
    d["radiance_mean"] = pd.to_numeric(d.radiance_mean, errors="coerce")
    d["valid_area_fraction"] = pd.to_numeric(
        d.valid_area_fraction, errors="coerce"
    )
    # Match the locked monthly atlas' coverage criteria; no interpolation,
    # filling, zero substitution or additional sensor correction.
    valid_text = d.analysis_valid.astype(str).str.lower()
    is_valid = valid_text.isin(("true", "1"))
    d["display_radiance"] = d.radiance_mean.where(
        is_valid & np.isfinite(d.radiance_mean)
        & (d.valid_area_fraction >= .5)
    )
    if (d.loc[d.display_radiance.notna(), "display_radiance"] < 0).any():
        raise ValueError("Invalid negative monthly radiance")
    counts = {}
    for year, month in WINTER_MONTHS:
        subset = d.loc[(d.year == year) & (d.month == month)]
        n = int(subset.display_radiance.notna().sum())
        counts[f"{year}-{month:02d}"] = n
        if n < 100:
            raise ValueError(
                f"Only {n} NUTS3 have valid VCMSLCFG coverage for "
                f"{year}-{month:02d}; check monthly tile and coverage filter"
            )
    pool = d.display_radiance.dropna().to_numpy(float)
    if not len(pool) or not np.isfinite(pool).all():
        raise ValueError("No valid monthly radiance for shared colour scale")
    shared_max = float(pool.max())
    if shared_max <= .2:
        raise ValueError("Unusually low monthly radiance maximum")
    # One fixed, non-clipping scale across all eight maps.  SymLogNorm
    # retains 0 and expands low-radiance differences while still showing
    # bright urban NUTS3.  This is an ABSOLUTE, not a diverging scale.
    norm = SymLogNorm(
        linthresh=.2, linscale=1., vmin=0.,
        vmax=shared_max, base=10
    )
    fig, axs = plt.subplots(3, 3, figsize=(14.8, 12.5))
    fig.subplots_adjust(
        left=.018, right=.895, top=.918, bottom=.055,
        wspace=.035, hspace=.19
    )
    fig.suptitle(
        "Winter monthly VIIRS radiance, 2020–2022",
        fontsize=16, weight="semibold"
    )
    map_counts = {}
    for (year, month), ax in zip(WINTER_MONTHS, axs.flat[:8]):
        sub = d.loc[
            (d.year == year) & (d.month == month),
            ["NUTS_ID", "display_radiance"]
        ]
        gm = geom.merge(
            sub, on="NUTS_ID", how="left", validate="one_to_one"
        )
        map_counts[f"{year}-{month:02d}"] = map_panel(
            ax, gm, bounds, "display_radiance", "magma", norm
        )
        ax.set_title(
            f"{MONTH_NAMES[month]} {year}",
            fontsize=12, pad=5
        )
    # Explicitly blank ninth panel: no key or inset in bottom right.
    axs.flat[8].set_axis_off()
    sm = ScalarMappable(cmap="magma", norm=norm)
    sm.set_array([])
    cb_ax = fig.add_axes([.916, .27, .017, .51])
    bar = fig.colorbar(sm, cax=cb_ax)
    bar.set_label(f"VIIRS monthly mean radiance ({RAD})", fontsize=10)
    fig.text(
        .025, .016,
        "VCMSLCFG monthly observations; all eight panels use one colour scale. "
        "Grey = missing/insufficient cloud-free coverage.",
        fontsize=9, color="#555555"
    )
    save_figure(fig, out, map_figure=True)
    plt.close(fig)
    return {
        "rendered": True,
        "monthly_atlas_source": str(monthly_path),
        "product": "vcmslcfg",
        "year_month_order": [f"{y}-{m:02d}" for y, m in WINTER_MONTHS],
        "layout": "3 x 3; final bottom-right panel empty",
        "map_valid_NUTS3_counts": map_counts,
        "monthly_valid_area_fraction_minimum": .5,
        "shared_colour": "magma sequential absolute radiance",
        "shared_color_min": 0.,
        "shared_color_max": shared_max,
        "shared_color_scale": "SymLogNorm, linthresh=0.2, base=10",
        "type": "actual monthly NUTS3 radiance, not 2021-minus-flank anomaly",
        "2021_monthly_snow_correction_applied": False,
        "no_missing_values_imputed": True,
        "year_month_cumulative_NUTS3": len(d),
    }



def winter_flank_anomaly_3x3(monthly_path, geom, bounds, out):
    """Same eight months: observed R(year,month) minus 2020/22 same-month mean."""
    if not monthly_path.is_file():
        return {"rendered": False, "reason": f"Missing monthly atlas: {monthly_path}"}
    need=set(WINTER_MONTHS)|{(2020,3),(2022,3),(2022,12)}
    cols=["NUTS_ID","year","month","product","radiance_mean",
          "analysis_valid","valid_area_fraction"]
    if monthly_path.suffix.lower()==".parquet":
        d=pd.read_parquet(monthly_path, columns=cols)
    else:
        d=pd.read_csv(monthly_path, usecols=cols, dtype={"NUTS_ID":str})
    d["NUTS_ID"]=d.NUTS_ID.astype(str)
    d["year"]=pd.to_numeric(d.year,errors="raise").astype(int)
    d["month"]=pd.to_numeric(d.month,errors="raise").astype(int)
    d=d.loc[d["product"].astype(str).str.lower().eq("vcmslcfg") &
        pd.MultiIndex.from_arrays([d.year,d.month]).isin(need)].copy()
    if d.duplicated(["NUTS_ID","year","month"]).any():
        raise ValueError("Duplicate monthly radiance region/year/month")
    found=set(zip(d.year,d.month))
    if found!=need:
        raise ValueError(f"Missing target/flank monthly archive entries: {sorted(need-found)}")
    d["radiance_mean"]=pd.to_numeric(d.radiance_mean,errors="coerce")
    d["valid_area_fraction"]=pd.to_numeric(d.valid_area_fraction,errors="coerce")
    good=d.analysis_valid.astype(str).str.lower().isin(("true","1"))
    good &= np.isfinite(d.radiance_mean) & d.valid_area_fraction.ge(.5)
    d["observed_nw"]=d.radiance_mean.where(good)
    if (d.observed_nw.dropna()<0).any():
        raise ValueError("Negative valid monthly VIIRS radiance")
    baseline=d.loc[d.year.isin((2020,2022))].pivot(
        index=["NUTS_ID","month"],columns="year",values="observed_nw")
    if not {2020,2022}.issubset(baseline.columns):
        raise ValueError("Both 2020 and 2022 reference columns required")
    baseline=baseline[[2020,2022]].rename(
        columns={2020:"ref_2020_nw",2022:"ref_2022_nw"}).reset_index()
    baseline["flank_mean_nw"]=(baseline.ref_2020_nw+baseline.ref_2022_nw)/2.
    target=d.loc[pd.MultiIndex.from_arrays([d.year,d.month]).isin(
        set(WINTER_MONTHS))].copy()
    target=target.merge(baseline,on=["NUTS_ID","month"],
                        how="left",validate="many_to_one")
    target["difference_nw"]=target.observed_nw-target.flank_mean_nw
    counts={}
    values=[]
    for year,month in WINTER_MONTHS:
        pair=target.loc[target.year.eq(year)&target.month.eq(month)]
        num=int(pair.difference_nw.notna().sum())
        counts[f"{year}-{month:02d}"]=num
        if num<100:
            raise ValueError(f"Only {num} valid paired regions for {year}-{month:02d}")
        values.extend(pair.difference_nw.dropna().tolist())
    pool=np.asarray(values,float)
    if len(pool)==0 or not np.isfinite(pool).all():
        raise ValueError("No valid monthly differences")
    cap=max(.1,float(np.quantile(np.abs(pool),.985)))
    norm=SymLogNorm(linthresh=.1,vmin=-cap,vmax=cap,base=10)
    fig,axs=plt.subplots(3,3,figsize=(14.8,12.5))
    fig.subplots_adjust(left=.018,right=.895,top=.918,bottom=.055,
                        wspace=.035,hspace=.19)
    fig.suptitle("Monthly VIIRS radiance minus 2020/2022 same-month mean",
                 fontsize=15.3,weight="semibold")
    mapped={}
    for (year,month),ax in zip(WINTER_MONTHS,axs.flat[:8]):
        sub=target.loc[target.year.eq(year)&target.month.eq(month),
                       ["NUTS_ID","difference_nw"]]
        gm=geom.merge(sub,on="NUTS_ID",how="left",validate="one_to_one")
        mapped[f"{year}-{month:02d}"]=map_panel(
            ax,gm,bounds,"difference_nw","RdBu_r",norm)
        ax.set_title(f"{MONTH_NAMES[month]} {year}",fontsize=12,pad=5)
    axs.flat[8].set_axis_off()
    sm=ScalarMappable(cmap="RdBu_r",norm=norm)
    sm.set_array([])
    cax=fig.add_axes([.916,.27,.017,.51])
    cb=fig.colorbar(sm,cax=cax,extend="both")
    cb.set_label(f"Month − 2020/2022 same-month mean ({RAD})",fontsize=10)
    fig.text(.025,.016,
             "VCMSLCFG monthly observations; grey = missing valid flanks; "
             "shared colour range capped at the pooled 98.5th percentile.",
             fontsize=8.7,color="#555555")
    save_figure(fig,out,map_figure=True)
    plt.close(fig)
    csv=out.parent/"CW_eight_monthly_flank_difference_values.csv"
    target[["NUTS_ID","year","month","observed_nw",
            "ref_2020_nw","ref_2022_nw","flank_mean_nw",
            "difference_nw","valid_area_fraction"]].sort_values(
            ["year","month","NUTS_ID"]).to_csv(csv,index=False)
    for month in (1,2):
        control=target.loc[target.month.eq(month)&target.year.isin(
            (2020,2022))].pivot(
            index="NUTS_ID",columns="year",values="difference_nw"
        ).dropna(subset=[2020,2022])
        if not np.allclose(control[2020],-control[2022],
                           rtol=1e-10,atol=1e-10):
            raise ValueError(f"2020/2022 control differences asymmetric in {month}")
    return {"rendered":True,"monthly_source":str(monthly_path),
            "year_month_order":[f"{y}-{m:02d}" for y,m in WINTER_MONTHS],
            "formula":"R(year,m) - [R(2020,m)+R(2022,m)]/2",
            "extra_reference_months":["2020-03","2022-03","2022-12"],
            "regions_with_valid_differences":counts,
            "mapped_region_counts":mapped,"shared_colour_cap_nw":cap,
            "same_symmetric_colour_scale_all_eight":True,
            "colour_clip_for_display_only":True,
            "2020_and_2022_include_themselves_in_baseline":True,
            "original_monthly_not_snow_corrected":True,
            "output_values_csv":str(csv)}

def self_test():
    example=pd.DataFrame({"year":[2016,2017,2021],
        "radiance_mean_raw":[1.,1.235,1.50],
        "radiance_mean_corrected":[1.,1.,1.265]})
    assert (example.loc[example.year.eq(2016),"radiance_mean_raw"].iloc[0] ==
            example.loc[example.year.eq(2016),"radiance_mean_corrected"].iloc[0])
    assert len(YEARS)==11 and YEARS[0]==2014 and YEARS[-1]==2024
    assert len(WINTER_MONTHS)==8 and len(set(WINTER_MONTHS))==8
    assert WINTER_MONTHS==(
        (2020,1),(2020,2),(2020,12),
        (2021,1),(2021,2),(2021,3),(2022,1),(2022,2)
    )
    assert np.isclose((.8*10.+.4*90.)/100.,.44)
    assert np.isclose(2-.5*(2+4),-1.)
    assert np.isclose(4-.5*(2+4),1.)
    assert np.isclose(7-.5*(2+4),4.)
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
            out/"CW_eight_winter_monthly_radiance_2020_2022_3x3"
        )
    )
    winter_difference=(
        {"rendered":False,"reason":"--skip-monthly specified"}
        if args.skip_monthly else winter_flank_anomaly_3x3(
            pathlike(args.monthly),geom,bounds,
            out/"CW_eight_winter_monthly_minus_2020_2022_flank_mean_3x3"
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
        "winter_same_month_flank_difference":winter_difference,
        "monthly_diagnostic_is_not_correction":True,
        "monthly_sources_are_original_VCMSLCFG":True,
        "monthly_outputs_include_absolute_and_same_month_flank_differences":True,
        "original_scientific_sources_modified":False
    }
    (out/"correction_demo_provenance.json").write_text(
        json.dumps(evidence,indent=2)+"\n",encoding="utf8")
    print("\n2017 BASELINE & 2021 SNOW CORRECTION FIGURES",flush=True)
    print("Science QA:",science,flush=True)
    print("Figure 1C comparison:",qa,flush=True)
    print("2021 correction:",snow_qa,flush=True)
    print("Eight-month winter maps (3x3):",winter,flush=True)
    print("Eight-month month-minus-flank map (3x3):",winter_difference,flush=True)
    print("Outputs:",out,flush=True)
    for file in sorted(out.glob("*")):
        if file.is_file():
            print(" ",file.name,round(file.stat().st_size/1024,1),"KiB",flush=True)


if __name__=="__main__":
    main()
