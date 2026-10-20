#!/usr/bin/env python3
"""Map locked spider/moth PRC change in fixed NUTS2016 radiance-density classes.

Reuses historical NUTS3 population-density assignment (median 2014-2020),
NOT year-specific NUTS 2021/2024 classifications. No PRC refitting.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from historical_prc_pixel_nuts3 import load_nuts
from project_config import PROJECT_ROOT,load_config

ROOT=PROJECT_ROOT
R_LABELS=["<0.2","0.2–1","1–5","5–20","≥20"]
D_LABELS=["<30","30-100","100-300","300-1000",">=1000"]
COLORS=["#d9ecf1","#91c5c2","#429c9c","#277276","#22464c"]

def inputs():
    processed=ROOT/"data/processed"
    paired=processed/"paired_nuts3_radiance_changes_2014_2024.csv"
    reference=processed/"fixed_nuts2016_population_density_reference.csv"
    raw=(ROOT/"figures/prc_radiance_spider_moth_2014_2024"/
         "nuts3_radiance_spider_moth_cumulative_change_2014_2024.csv")
    if not raw.is_file():raise FileNotFoundError(raw)
    if paired.is_file():
        density=pd.read_csv(paired,usecols=["NUTS_ID","density_band",
                                             "reference_population_density"])
    elif reference.is_file():
        density=pd.read_csv(reference)[["NUTS_ID","density_band",
                                      "reference_population_density"]]
    else:raise FileNotFoundError("Existing NUTS2016 density reference missing")
    source=pd.read_csv(raw,dtype={"NUTS_ID":str})
    density.NUTS_ID=density.NUTS_ID.astype(str)
    density=density.drop_duplicates("NUTS_ID")
    if not density.density_band.dropna().isin(D_LABELS).all():
        raise ValueError("Density classes disagree with historical model")
    data=source.merge(density,on="NUTS_ID",how="left",validate="one_to_one")
    data["radiance_band"]=pd.cut(
        data.radiance_2014,[-np.inf,.2,1,5,20,np.inf],
        labels=R_LABELS,ordered=True,right=False)
    if len(data)!=1522 or not data.NUTS_ID.is_unique:
        raise ValueError("Expected same 1522 region cohort as two-taxon PRC")
    return data

def counts_and_effects(data,out):
    ok=data.dropna(subset=["density_band","radiance_band"]).copy()
    records=[]
    for rb in R_LABELS:
        for db in D_LABELS:
            x=ok[ok.radiance_band.eq(rb)&ok.density_band.eq(db)]
            a=x.valid_area_km2_2014
            records.append({"radiance_band":rb,"density_band":db,"n_nuts3":len(x),
                "area_km2":a.sum(),"median_density":x.reference_population_density.median(),
                "radiance_abs_area_sum":(a*x.radiance_change_abs.abs()).sum(),
                "spider_prc_abs_area_sum":(a*x.spider_prc_change_abs.abs()).sum(),
                "moth_prc_abs_area_sum":(a*x.moth_prc_change_abs.abs()).sum(),
                "median_radiance_change_nw":x.radiance_change_abs.median(),
                "median_spider_prc_change":x.spider_prc_change_abs.median(),
                "median_moth_prc_change":x.moth_prc_change_abs.median()})
    table=pd.DataFrame(records)
    for key in ["area_km2","radiance_abs_area_sum","spider_prc_abs_area_sum",
                "moth_prc_abs_area_sum"]:
        table[key+"_percent"]=100*table[key]/table[key].sum()
    table.to_csv(out/"radiance_density_joint_contributions.csv",index=False)
    print(table.sort_values("spider_prc_abs_area_sum",ascending=False).head(8).to_string(index=False))
    return table

def draw_base(ax,g,col,title,labels,colors):
    for idx,label in enumerate(labels):
        z=g[g[col].astype("string").eq(label)]
        if len(z):z.plot(ax=ax,color=colors[idx],linewidth=0,edgecolor="none")
    ax.set(xlim=(-12,36),ylim=(35,73),title=title,
           xlabel="Longitude",ylabel="Latitude")
    ax.set_aspect("auto")
    handles=[Patch(facecolor=colors[i],label=t) for i,t in enumerate(labels)]
    ax.legend(handles=handles,loc="lower left",fontsize=8,framealpha=.9)

def draw_focus(ax,g,band):
    x=g[g.radiance_band.astype("string").eq(band)]
    g.plot(ax=ax,color="#ebebeb",linewidth=0)
    draw_base(ax,x,"density_band",f"2014 radiance {band} nW: population class",
              D_LABELS,COLORS)

def maps(data,out):
    nuts=load_nuts(load_config())[["NUTS_ID","geometry"]]
    g=nuts.merge(data[["NUTS_ID","radiance_band","density_band"]],on="NUTS_ID",
                 how="left",validate="one_to_one")
    fig,axs=plt.subplots(2,2,figsize=(15,15),layout="constrained")
    draw_base(axs[0,0],g,"radiance_band","2014 NUTS3 mean VIIRS radiance (nW)",
              R_LABELS,["#e3eef4","#b7d8ec","#519dc1","#216a92","#163853"])
    draw_base(axs[0,1],g,"density_band",
              "Historical NUTS2016 median density, 2014–2020 (people/km²)",
              D_LABELS,COLORS)
    draw_focus(axs[1,0],g,"0.2–1")
    draw_focus(axs[1,1],g,"1–5")
    fig.suptitle("European radiance/settlement geography: fixed NUTS2016 NUTS3",
                 fontsize=16)
    fig.savefig(out/"baseline_radiance_density_four_maps.pdf",bbox_inches="tight")
    plt.close(fig)

def observed_population_change(data,out):
    file=ROOT/"data/processed/nuts3_population_density.parquet"
    if not file.exists():
        print("No original Eurostat population parquet: skipping measured 2014–2020 change")
        return
    pop=pd.read_parquet(file)
    for column in ["NUTS_ID","year","population_density","nuts_release"]:
        if column not in pop:raise ValueError(f"Eurostat source missing {column}")
    pop=pop[pop.year.isin([2014,2020]) & pop.nuts_release.eq(2016)].copy()
    if "matches_nuts_release" in pop:
        pop=pop[pop.matches_nuts_release]
    pop=pop.pivot_table(index="NUTS_ID",columns="year",
                        values="population_density",aggfunc="first")
    if not {2014,2020}.issubset(pop):raise ValueError("Need matched 2014/2020 density")
    pop=pop[[2014,2020]].rename(columns={2014:"density_2014",2020:"density_2020"})
    pop["pop_density_change_pct_2014_2020"]=100*(pop.density_2020/pop.density_2014-1)
    merged=data.merge(pop.reset_index(),on="NUTS_ID",how="left",validate="one_to_one")
    merged.to_csv(out/"joined_observed_2014_2020_demography_with_2014_2024_prc.csv",
                  index=False)
    print("Observed fixed-boundary population 2014–2020 matches:",
        merged.pop_density_change_pct_2014_2020.notna().sum())
    print("This is NOT a 2014–2024 population-change inference.")

def main():
    data=inputs()
    out=ROOT/"figures/prc_population_density_2014_2024"
    out.mkdir(parents=True,exist_ok=True)
    data.to_csv(out/"nuts3_radiance_density_spider_moth_joined.csv",index=False)
    counts_and_effects(data,out)
    maps(data,out)
    observed_population_change(data,out)
    print("Files written under:",out)
if __name__=="__main__":main()
