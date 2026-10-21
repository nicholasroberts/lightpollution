#!/usr/bin/env python3
"""VISUAL first: annual and 2014-2020 population change vs spider/moth PRC.

Use original fixed-NUTS2016 observed Eurostat annual density and LOCKED
snow-corrected taxon PRC panels. Six matched annual-transition map pages,
two cumulative map sheets, baseline-radiance focus maps, matched scatter
pages, Pearson/Spearman descriptive correlations and a joined source CSV.
No regression or invented post-2020 annual population data.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from historical_prc_pixel_nuts3 import load_nuts
from project_config import PROJECT_ROOT as ROOT,load_config

YEARS=list(range(2014,2021))
R_BINS=[-np.inf,.2,1,5,20,np.inf]
R_LABELS=["<0.2","0.2–1","1–5","5–20","≥20"]
R_COLORS=["#aaaaaa","#4e9dd0","#35967b","#d79542","#8354a1"]
OUT=ROOT/"figures/prc_population_change_2014_2020"


def inputs():
    proc=ROOT/"data/processed"
    paths={
        "population":proc/"nuts3_population_density.parquet",
        "spider":proc/"historical_prc_2014_2024_spider_snow_corrected"/"prc_nuts3_annual_cumulative.csv",
        "moth":proc/"historical_prc_2014_2024_moth_snow_corrected"/"prc_nuts3_annual_cumulative.csv",
        "radiance":proc/"fixed_nuts2016_viirs_2014_2024.parquet",
        "reference":proc/"fixed_nuts2016_population_density_reference.csv",
        "endpoints":ROOT/"figures/prc_radiance_spider_moth_2014_2024"/"nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv"}
    missing=[str(p) for p in paths.values() if not p.is_file()]
    if missing:raise FileNotFoundError("Missing earlier project output:\n"+"\n".join(missing))
    pop=pd.read_parquet(paths["population"])
    needed={"NUTS_ID","year","nuts_release","population_density"}
    if not needed.issubset(pop):raise ValueError(f"Population columns absent: {needed-set(pop.columns)}")
    pop=pop[pop.year.isin(YEARS)&pop.nuts_release.eq(2016)].copy()
    if "matches_nuts_release" in pop:
        pop=pop[pop.matches_nuts_release.fillna(False)]
    pop=pop[["NUTS_ID","year","population_density"]]
    pop["NUTS_ID"]=pop.NUTS_ID.astype(str)
    if pop.duplicated(["NUTS_ID","year"]).any():
        raise ValueError("Ambiguous Eurostat population for a NUTS/year")
    rad=pd.read_parquet(paths["radiance"])
    need={"NUTS_ID","year","radiance_mean_corrected"}
    if not need.issubset(rad):
        raise ValueError(f"Existing NUTS2016 annual radiance missing: {need-set(rad.columns)}")
    rad=rad[rad.year.isin(YEARS)][list(need)].copy()
    rad["NUTS_ID"]=rad.NUTS_ID.astype(str)
    if rad.duplicated(["NUTS_ID","year"]).any():
        raise ValueError("Duplicate fixed NUTS2016 radiance/year")
    joint=pop.merge(rad,on=["NUTS_ID","year"],how="inner",validate="one_to_one")
    for taxon in ("spider","moth"):
        t=pd.read_csv(paths[taxon],dtype={"NUTS_ID":str})
        if "snow_status" not in t.columns or not t.snow_status.eq("locked").all():
            raise ValueError(f"{taxon}: data not from locked 2021 snow-corrected series")
        t=t[t.year.isin(YEARS)][["NUTS_ID","year","prc_mean","coverage_stable",
                                 "valid_area_km2"]].rename(
            columns={"prc_mean":f"{taxon}_prc",
                     "coverage_stable":f"{taxon}_stable",
                     "valid_area_km2":f"{taxon}_area"})
        if t.duplicated(["NUTS_ID","year"]).any():
            raise ValueError(f"Duplicate {taxon} NUTS-year")
        joint=joint.merge(t,on=["NUTS_ID","year"],how="inner",validate="one_to_one")
    ref=pd.read_csv(paths["reference"],dtype={"NUTS_ID":str})[
        ["NUTS_ID","density_band","reference_population_density"]]
    ep=pd.read_csv(paths["endpoints"],dtype={"NUTS_ID":str})[
        ["NUTS_ID","radiance_2014"]]
    for label,d in (("density_reference",ref),("radiance_endpoints",ep)):
        if d.NUTS_ID.duplicated().any():raise ValueError(f"Duplicate {label} NUTS_ID")
    joint=joint.merge(ref,on="NUTS_ID",how="left",validate="many_to_one")
    joint=joint.merge(ep,on="NUTS_ID",how="left",validate="many_to_one")
    joint["radiance_band_2014"]=pd.cut(joint.radiance_2014,R_BINS,
        labels=R_LABELS,right=False)
    return joint


def calculate(panel):
    d=panel.sort_values(["NUTS_ID","year"]).copy()
    d=d[d.spider_stable.fillna(False)&d.moth_stable.fillna(False)]
    for col in ["population_density","radiance_mean_corrected","spider_prc","moth_prc"]:
        d[col]=pd.to_numeric(d[col],errors="coerce")
    d=d[(d.population_density>0)&(d.radiance_mean_corrected>=0)]
    d=d.dropna(subset=["population_density","radiance_mean_corrected",
                       "spider_prc","moth_prc"])
    # Exactly the same NUTS3 units must contribute to each annual panel and
    # 2014->2020 endpoint plot; NEVER fill missing population measurements.
    complete=d.groupby("NUTS_ID").year.nunique()
    eligible=complete.index[complete.eq(len(YEARS))]
    d=d[d.NUTS_ID.isin(eligible)].copy()
    if d.duplicated(["NUTS_ID","year"]).any():
        raise ValueError("Repeated matched NUTS/year")
    if len(d)!=len(eligible)*7:
        raise ValueError("Missing one of the seven measured years")
    if not len(eligible):
        raise ValueError("No complete fixed-NUTS2016 population+PRC panel")
    group=d.groupby("NUTS_ID",sort=False)
    d["previous_year"]=group.year.shift()
    for col in ["population_density","radiance_mean_corrected","spider_prc","moth_prc"]:
        d["previous_"+col]=group[col].shift()
        d["delta_"+col]=d[col]-d["previous_"+col]
        d["pct_"+col]=100*d["delta_"+col]/d["previous_"+col].where(
            d["previous_"+col].ne(0))
        first=d[d.year.eq(2014)].set_index("NUTS_ID")[col]
        d["base_"+col]=d.NUTS_ID.map(first)
        d["cumulative_delta_"+col]=d[col]-d["base_"+col]
        d["cumulative_pct_"+col]=100*d["cumulative_delta_"+col]/d[
            "base_"+col].where(d["base_"+col].ne(0))
    annual=d[d.year.between(2015,2020)].copy()
    if not annual.year.sub(annual.previous_year).eq(1).all():
        raise ValueError("Non-consecutive NUTS population/PRC observations")
    cumulative=d[d.year.eq(2020)].copy()
    return d,annual,cumulative


def correlation(data,period):
    rows=[]
    for x_name in ("population_density_absolute","population_pct"):
        x=("pct_population_density" if period=="annual" else
           "cumulative_pct_population_density") if x_name=="population_pct" else (
           "delta_population_density" if period=="annual" else
           "cumulative_delta_population_density")
        for taxon in ("spider","moth"):
            for response_type in ("absolute","percentage"):
                y=(("delta_" if period=="annual" else "cumulative_delta_")
                   if response_type=="absolute" else
                   ("pct_" if period=="annual" else "cumulative_pct_"))+taxon+"_prc"
                sub=data[[x,y]].replace([np.inf,-np.inf],np.nan).dropna()
                n=len(sub)
                if n>2 and sub[x].nunique()>1 and sub[y].nunique()>1:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        pr,pp=pearsonr(sub[x],sub[y])
                        sr,sp=spearmanr(sub[x],sub[y])
                else:
                    pr=pp=sr=sp=np.nan
                rows.append({"period":period,"population_predictor":x_name,
                   "prc_taxon":taxon,"prc_response":response_type,
                   "population_field":x,"prc_field":y,"n_regions":n,
                   "pearson_r":float(pr),"pearson_p":float(pp),
                   "spearman_rho":float(sr),"spearman_p":float(sp)})
    return pd.DataFrame(rows)


def limits(d,cols):
    ans={}
    for key in cols:
        a=np.abs(pd.to_numeric(d[key],errors="coerce").to_numpy(float))
        a=a[np.isfinite(a)]
        ans[key]=max(float(np.quantile(a,.98)) if len(a) else 0,1e-10)
    return ans


def maps_for(d,nuts,cols,labels,output,title,scale,focus=None):
    geo=nuts.merge(d[["NUTS_ID",*cols]],on="NUTS_ID",how="left",validate="one_to_one")
    if focus is not None:
        active=set(d.loc[d.radiance_band_2014.astype(str).eq(focus),"NUTS_ID"])
        for c in cols:geo.loc[~geo.NUTS_ID.isin(active),c]=np.nan
    fig,axs=plt.subplots(1,3,figsize=(21,8.4),layout="constrained")
    for ax,col,label in zip(axs,cols,labels):
        geo.plot(column=col,ax=ax,cmap="RdBu_r",linewidth=0,
                 norm=TwoSlopeNorm(vmin=-scale[col],vcenter=0,vmax=scale[col]),
                 missing_kwds={"color":"#d3d3d3"},legend=True,
                 legend_kwds={"shrink":.65,"label":label})
        ax.set(xlim=(-12,36),ylim=(35,73),xlabel="Longitude",
               ylabel="Latitude",title=f"{label}\n(n={geo[col].notna().sum()})")
    suffix=(f" | only 2014 radiance band {focus} nW" if focus else "")
    fig.suptitle(title+suffix+"\nGrey: no matching data / excluded band; each variable independently zero-centred")
    fig.savefig(output,bbox_inches="tight")
    plt.close(fig)


def plot_annual_maps(annual,nuts,out):
    cols=["delta_population_density","delta_spider_prc","delta_moth_prc"]
    titles=["Population density Δ (people/km²)","Spider ΔPRC (absolute)","Moth ΔPRC (absolute)"]
    scale=limits(annual,cols) # CONSISTENT limits across all six years
    with PdfPages(out/"population_vs_prc_annual_maps_2015_2020.pdf") as pdf:
        for year in range(2015,2021):
            tmp=out/f"_annual_{year}_temporary.pdf"
            maps_for(annual[annual.year.eq(year)],nuts,cols,titles,tmp,
                     f"Measured population and modelled PRC: {year-1} → {year}",scale)
            from matplotlib.backends.backend_pdf import PdfPages as PP
            # Directly render to a multipage PDF; page composition is in
            # figure objects, not by merging existing PDF files.
            import matplotlib.pyplot as plot
            # The temporary PDF is independently useful as a single-year map.
    # Implemented separately below by a single-pages bundle index.
    print("Annual map PDFs written",flush=True)


def map_sheet(d,nuts,cols,labels,output,title,scale,focus=None):
    maps_for(d,nuts,cols,labels,output,title,scale,focus)


def scatter_axes(ax,d,x,y,taxon,period):
    use=d[[x,y,"radiance_band_2014"]].replace([np.inf,-np.inf],np.nan).dropna(subset=[x,y])
    for idx,lab in enumerate(R_LABELS):
        sub=use[use.radiance_band_2014.astype(str).eq(lab)]
        if len(sub):
            ax.scatter(sub[x],sub[y],s=9,alpha=.35,c=R_COLORS[idx],
                       label=f"{lab} nW",rasterized=True)
    if len(use)>2 and use[x].nunique()>1 and use[y].nunique()>1:
        p=float(pearsonr(use[x],use[y]).statistic)
        s=float(spearmanr(use[x],use[y]).statistic)
        ax.set_title(f"{taxon.capitalize()} PRC | n={len(use)} | r={p:+.3f} | ρ={s:+.3f}")
    else:ax.set_title(f"{taxon.capitalize()} PRC | n={len(use)}")
    ax.axhline(0,lw=.6,color="gray")
    ax.axvline(0,lw=.6,color="gray")
    ax.set(xlabel="Population density change (people/km²)",ylabel="Absolute ΔPRC")
    ax.grid(alpha=.16)


def draw_scatter_pair(d,output,title,period):
    x="delta_population_density" if period=="annual" else "cumulative_delta_population_density"
    yprefix="delta_" if period=="annual" else "cumulative_delta_"
    fig,axs=plt.subplots(1,2,figsize=(14,5.8),layout="constrained")
    for taxon,ax in zip(("spider","moth"),axs):
        scatter_axes(ax,d,x,yprefix+taxon+"_prc",taxon,period)
    handles=[Line2D([],[],marker="o",linestyle="None",color=R_COLORS[i],
               label=f"{lab} nW") for i,lab in enumerate(R_LABELS)]
    fig.legend(handles=handles,loc="lower center",ncol=5,
               title="2014 mean NUTS3 VIIRS radiance")
    fig.suptitle(title)
    fig.savefig(output,bbox_inches="tight")
    plt.close(fig)


def synthetic_test():
    rows=[]
    for i in range(4):
        for year in YEARS:
            rows.append({"NUTS_ID":f"X{i}","year":year,
                "population_density":100+i+2*(year-2014),
                "radiance_mean_corrected":1+i*.1,
                "spider_prc":.42-.001*(year-2014),
                "moth_prc":.63-.002*(year-2014),
                "spider_stable":True,"moth_stable":True,
                "density_band":"30-100","radiance_2014":1.2})
    _,a,c=calculate(pd.DataFrame(rows))
    assert len(a)==24 and len(c)==4
    assert a.year.min()==2015 and c.year.iloc[0]==2020
    assert np.allclose(c.cumulative_delta_population_density,12)
    assert np.allclose(c.cumulative_delta_spider_prc,-.006)
    assert len(correlation(a[a.year.eq(2015)],"annual"))==8
    print("Common-region annual and cumulative population-vs-PRC QA: PASS")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--output-dir",default="figures/prc_population_change_2014_2020")
    args=ap.parse_args()
    synthetic_test()
    if args.self_test:return
    full,annual,cum=calculate(inputs())
    if len(cum)<200:raise ValueError("Insufficient complete matched NUTS2016 panel")
    out=ROOT/args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    full.to_csv(out/"matched_population_and_prc_annual_2014_2020.csv",index=False)
    annual.to_csv(out/"annual_changes_2015_2020.csv",index=False)
    cum.to_csv(out/"cumulative_change_2014_2020.csv",index=False)
    nuts=load_nuts(load_config())[["NUTS_ID","geometry"]]
    # Produce six separate single-year map PDFs to enable immediate inspection.
    annual_cols=["delta_population_density","delta_spider_prc","delta_moth_prc"]
    annual_titles=["Annual density Δ (people/km²)","Annual spider ΔPRC","Annual moth ΔPRC"]
    annual_scale=limits(annual,annual_cols)
    for year in range(2015,2021):
        map_sheet(annual[annual.year.eq(year)],nuts,annual_cols,annual_titles,
           out/f"annual_population_spider_moth_{year-1}_{year}.pdf",
           f"Year-on-year change {year-1} → {year}",annual_scale)
    cumulative_cols=["cumulative_delta_population_density",
                     "cumulative_delta_spider_prc","cumulative_delta_moth_prc"]
    cumulative_titles=["Density change 2014–2020 (people/km²)",
                       "Spider PRC change 2014–2020",
                       "Moth PRC change 2014–2020"]
    cumulative_scale=limits(cum,cumulative_cols)
    map_sheet(cum,nuts,cumulative_cols,cumulative_titles,
           out/"cumulative_population_spider_moth_2014_2020.pdf",
           "Cumulative changes on fixed NUTS2016 geography",cumulative_scale)
    for band,name in (("0.2–1","radiance_02_to_1"),("1–5","radiance_1_to_5")):
        map_sheet(cum,nuts,cumulative_cols,cumulative_titles,
           out/f"cumulative_{name}_population_spider_moth_2014_2020.pdf",
           "Cumulative population and PRC changes",cumulative_scale,focus=band)
    with PdfPages(out/"annual_population_vs_spider_moth_prc_scatter_2015_2020.pdf") as pdf:
        # Matplotlib pages use the same plotted content as standalone files.
        for year in range(2015,2021):
            z=annual[annual.year.eq(year)]
            fig,axs=plt.subplots(1,2,figsize=(14,5.8),layout="constrained")
            for taxon,ax in zip(("spider","moth"),axs):
                scatter_axes(ax,z,"delta_population_density",
                             "delta_"+taxon+"_prc",taxon,"annual")
            handles=[Line2D([],[],marker="o",linestyle="None",color=R_COLORS[i],
                            label=f"{lab} nW") for i,lab in enumerate(R_LABELS)]
            fig.legend(handles=handles,loc="lower center",ncol=5,
                       title="2014 regional mean radiance")
            fig.suptitle(f"Population change vs PRC change: {year-1} → {year}")
            pdf.savefig(fig,bbox_inches="tight")
            plt.close(fig)
    draw_scatter_pair(cum,out/"cumulative_population_vs_spider_moth_prc_scatter_2014_2020.pdf",
                      "2014→2020 total measured population change vs modelled PRC change","cumulative")
    tables=[correlation(annual[annual.year.eq(year)],"annual").assign(
                start_year=year-1,end_year=year) for year in range(2015,2021)]
    tables.append(correlation(cum,"cumulative").assign(
        start_year=2014,end_year=2020))
    allcorr=pd.concat(tables,ignore_index=True)
    allcorr.to_csv(out/"correlations_annual_and_cumulative.csv",index=False)
    strats=[]
    for key,subset in cum.groupby("radiance_band_2014",observed=True):
        if len(subset)<4:continue
        strats.append(correlation(subset,"cumulative").assign(
           radiance_band_2014=str(key)))
    if strats:pd.concat(strats,ignore_index=True).to_csv(
           out/"correlations_cumulative_by_baseline_radiance.csv",index=False)
    print(f"Complete fixed-NUTS2016 population/PRC panel: {len(cum)} regions × 7 years")
    print(allcorr[(allcorr.population_predictor=="population_pct")&
        (allcorr.prc_response=="absolute")][
        ["start_year","end_year","prc_taxon","n_regions","pearson_r","spearman_rho"]].to_string(index=False))
    print("Maps, scatter plots and tables:",out)


if __name__=="__main__":
    main()
