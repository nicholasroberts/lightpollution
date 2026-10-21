#!/usr/bin/env python3
"""Country-by-country ABSOLUTE density change vs PRC analysis (2014–2020).

Visualize whether population GROWTH accompanies PRC LOSS, whereas population
DECLINE does not reliably accompany PRC GAIN. Preserve country baselines, both
taxa, fixed NUTS2016 cohort, all annual transitions and 2014–20 cumulative.
Not a causal attribution model; LED and lighting policies are not controlled.
"""
from __future__ import annotations
import argparse
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, fisher_exact

from visualize_annual_population_prc_changes import inputs, calculate
from historical_prc_pixel_nuts3 import load_nuts
from project_config import PROJECT_ROOT as ROOT, load_config

OUT = ROOT/"figures/prc_population_country_absolute_density_2014_2020"
TAXA=("spider","moth")
# Fixed visual sign thresholds, with an unthresholded sensitivity in tables.
CUTS={"cumulative":(5.0,.001),"annual":(1.0,.0005)}
# Main thresholds: people/km² and absolute PRC. Zero-only signs are also reported.
QUAD={
 "population_up__prc_down": ("#ad3b35","Population ↑ / PRC ↓"),
 "population_down__prc_up": ("#2879a3","Population ↓ / PRC ↑"),
 "population_up__prc_up": ("#d49b36","Population ↑ / PRC ↑"),
 "population_down__prc_down": ("#8059a0","Population ↓ / PRC ↓"),
 "near_zero_or_mixed": ("#e1e1e1","Near-zero change / mixed")
}
MIN_CORR_N=8


def exact_corr(x,y):
    a=np.asarray(x,float)
    b=np.asarray(y,float)
    ok=np.isfinite(a)&np.isfinite(b)
    a,b=a[ok],b[ok]
    if len(a)<MIN_CORR_N or len(np.unique(a))<3 or len(np.unique(b))<3:
        return float("nan"),float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(pearsonr(a,b).statistic),float(spearmanr(a,b).statistic)


def classify(pop,prc,pop_min=0,prc_min=0):
    """Sign classifications use absolute density change (people/km²) and absolute PRC."""
    pop=np.asarray(pop,float)
    prc=np.asarray(prc,float)
    category=np.full(len(pop),"near_zero_or_mixed",dtype=object)
    category[(pop>pop_min)&(prc< -prc_min)]="population_up__prc_down"
    category[(pop< -pop_min)&(prc>prc_min)]="population_down__prc_up"
    category[(pop>pop_min)&(prc>prc_min)]="population_up__prc_up"
    category[(pop< -pop_min)&(prc< -prc_min)]="population_down__prc_down"
    return category


def magnitude_subgroups(x,y,prc_min):
    """Separate magnitude associations without constructing 'correlation' from signs alone."""
    x=np.asarray(x,float); y=np.asarray(y,float)
    loss=y< -prc_min
    gain=y>prc_min
    _,loss_rho=exact_corr(x[loss],-y[loss])
    _,gain_rho=exact_corr(-x[gain],y[gain])
    return int(loss.sum()),loss_rho,int(gain.sum()),gain_rho


def directional_summary(data,country,taxon,period,pop_min,prc_min):
    x=("delta_population_density" if period=="annual"
       else "cumulative_delta_population_density")
    y=(("delta_" if period=="annual" else "cumulative_delta_")
       +taxon+"_prc")
    z=data[["NUTS_ID",x,y]].dropna().copy()
    z=z[np.isfinite(z[x])&np.isfinite(z[y])]
    cats=classify(z[x],z[y],pop_min,prc_min)
    counts={k:int(np.sum(cats==k)) for k in QUAD}
    a=counts["population_up__prc_down"]
    b=counts["population_up__prc_up"]
    c=counts["population_down__prc_down"]
    d=counts["population_down__prc_up"]
    both=a+b+c+d
    raw_pearson,raw_spearman=exact_corr(z[x],z[y])
    nloss,loss_rho,ngain,gain_rho=magnitude_subgroups(z[x],z[y],prc_min)
    return {
        "country":country,"taxon":taxon,"period":period,
        "n_regions":len(z),"pop_threshold_people_per_km2":pop_min,"prc_threshold_abs":prc_min,
        "growth_count":int((z[x]>pop_min).sum()),
        "decline_count":int((z[x]<-pop_min).sum()),
        "prc_loss_count":nloss,"prc_gain_count":ngain,
        "n_growth_prc_loss":a,"n_growth_prc_gain":b,
        "n_decline_prc_loss":c,"n_decline_prc_gain":d,
        "n_signed_both":both,"n_neutral_or_mixed":len(z)-both,
        "fraction_growth_given_prc_loss":a/nloss if nloss else np.nan,
        "fraction_decline_given_prc_gain":d/ngain if ngain else np.nan,
        "fraction_prc_loss_given_growth":a/(a+b) if a+b else np.nan,
        "fraction_prc_gain_given_decline":d/(c+d) if c+d else np.nan,
        "growth_baseline_fraction":float((z[x]>pop_min).sum()/len(z)) if len(z) else np.nan,
        "prc_loss_baseline_fraction":float(nloss/len(z)) if len(z) else np.nan,
        "expected_direction_fraction":(a+d)/both if both else np.nan,
        "fisher_odds_growth_vs_loss":float(fisher_exact([[a,b],[c,d]])[0])
            if min(a+b,c+d,a+c,b+d)>0 else np.nan,
        "fisher_p":float(fisher_exact([[a,b],[c,d]])[1])
            if min(a+b,c+d,a+c,b+d)>0 else np.nan,
        "pearson_r_pop_vs_delta_prc":raw_pearson,
        "spearman_rho_pop_vs_delta_prc":raw_spearman,
        "n_prc_loss_for_magnitude":nloss,
        "spearman_rho_growth_vs_prc_loss_magnitude":loss_rho,
        "n_prc_gain_for_magnitude":ngain,
        "spearman_rho_decline_vs_prc_gain_magnitude":gain_rho
    }


def countries_by_year(annual,cumulative):
    nuts=load_nuts(load_config())
    if "CNTR_CODE" not in nuts.columns:
        raise ValueError("Fixed NUTS2016 polygons lack CNTR_CODE")
    if nuts.NUTS_ID.duplicated().any():
        raise ValueError("Non-unique NUTS2016 polygon NUTS_ID")
    meta=nuts[["NUTS_ID","CNTR_CODE","geometry"]].copy()
    meta.NUTS_ID=meta.NUTS_ID.astype(str)
    for d in (annual,cumulative):
        if d.NUTS_ID.duplicated().any() and d is cumulative:
            raise ValueError("Cumulative NUTS3 duplicates")
    a=annual.merge(meta[["NUTS_ID","CNTR_CODE"]],on="NUTS_ID",validate="many_to_one")
    c=cumulative.merge(meta[["NUTS_ID","CNTR_CODE"]],on="NUTS_ID",validate="one_to_one")
    return a.rename(columns={"CNTR_CODE":"country"}),c.rename(columns={"CNTR_CODE":"country"}),meta


def make_tables(annual,cumulative):
    rows=[]
    for period,df in (("cumulative",cumulative),("annual",annual)):
        years=[(2014,2020,cumulative)] if period=="cumulative" else [
            (y-1,y,annual[annual.year.eq(y)]) for y in range(2015,2021)]
        for start,end,frame in years:
            for country,g in [("ALL",frame),*list(frame.groupby("country"))]:
                for taxon in TAXA:
                    # Explicit nonzero noise floor and unthresholded sign check.
                    for threshold in ("main","zero"):
                        pcut,rcut=CUTS[period] if threshold=="main" else (0.,0.)
                        row=directional_summary(g,country,taxon,period,pcut,rcut)
                        row.update(start_year=start,end_year=end,threshold_case=threshold)
                        rows.append(row)
    return pd.DataFrame(rows)


def plot_country_scatter(cumulative,output):
    selected=["FR","TR","IS"]
    fig,axs=plt.subplots(len(selected),2,figsize=(15,15),layout="constrained")
    x="cumulative_delta_population_density"
    for ridx,country in enumerate(selected):
        region=cumulative[cumulative.country.eq(country)]
        for cidx,taxon in enumerate(TAXA):
            ax=axs[ridx,cidx]
            y="cumulative_delta_"+taxon+"_prc"
            category=classify(region[x],region[y],*CUTS["cumulative"])
            for key,(color,label) in QUAD.items():
                sel=category==key
                if sel.any():
                    ax.scatter(region.loc[sel,x],region.loc[sel,y],
                        s=26,alpha=.75,c=color,label=label,rasterized=True)
            r,rho=exact_corr(region[x],region[y])
            suffix=(f"r={r:+.3f}; ρ={rho:+.3f}"
                    if np.isfinite(r) else "correlation withheld: too few regions")
            ax.axhline(0,lw=.7,color="gray")
            ax.axvline(0,lw=.7,color="gray")
            ax.set(xlabel="Population density change 2014–20 (people/km²)",
                   ylabel="Absolute Δ"+taxon+" PRC",
                   title=f"{country} | {taxon} | n={len(region)} | {suffix}")
            ax.grid(alpha=.15)
    patches=[Patch(facecolor=v[0],label=v[1]) for v in QUAD.values()]
    fig.legend(handles=patches,loc="lower center",ncol=3,fontsize=8)
    fig.suptitle("Within-country direction of observed population and predicted PRC change")
    fig.savefig(output,bbox_inches="tight")
    plt.close(fig)


def country_metric_plot(t,out):
    d=t[(t.period=="cumulative")&(t.threshold_case=="main")&
        (t.country.ne("ALL"))&(t.n_regions.ge(8))].copy()
    for metric,name,caption in [
        ("spearman_rho_pop_vs_delta_prc","country_spearman",
         "2014–20 within-country Spearman ρ: Δpopulation density (people/km²) vs ΔPRC"),
        ("fraction_growth_given_prc_loss","loss_with_growth",
         "Fraction of PRC-loss regions with population growth"),
        ("fraction_decline_given_prc_gain","gain_with_decline",
         "Fraction of PRC-gain regions with population decline")]:
        wide=d.pivot(index="country",columns="taxon",values=metric)
        wide=wide.sort_values("spider",ascending=True)
        if wide.empty:continue
        fig,ax=plt.subplots(figsize=(11,max(5,0.32*len(wide)+1)),layout="constrained")
        y=np.arange(len(wide))
        for i,taxon in enumerate(TAXA):
            ax.scatter(wide[taxon],y+(i-.5)*.20,s=33,
                       label=taxon,marker="o" if taxon=="spider" else "x")
        if metric=="spearman_rho_pop_vs_delta_prc":
            ax.set_xlim(-1,1)
        else:ax.set_xlim(0,1)
        ax.axvline(0,color="gray",lw=.65)
        ax.set_yticks(y,wide.index)
        ax.set_xlabel(metric.replace("_"," "))
        ax.set_title(caption)
        ax.grid(axis="x",alpha=.18)
        ax.legend()
        fig.savefig(out/(name+"_by_country_2014_2020.pdf"),bbox_inches="tight")
        plt.close(fig)


def annual_heatmap(t,out):
    d=t[(t.period=="annual")&(t.threshold_case=="main")&
        (t.country.ne("ALL"))&(t.n_regions.ge(8))]
    for taxon in TAXA:
        q=d[d.taxon.eq(taxon)]
        wide=q.pivot(index="country",columns="end_year",
                     values="spearman_rho_pop_vs_delta_prc")
        if wide.empty:continue
        # Rank countries by number of six-year region matches rather than observed trend.
        sizes=q.groupby("country").n_regions.max()
        wide=wide.reindex(sizes.sort_values(ascending=False).index)
        fig,ax=plt.subplots(figsize=(10,max(6,len(wide)*.29)),layout="constrained")
        arr=wide.to_numpy(float)
        img=ax.imshow(np.ma.masked_invalid(arr),cmap="RdBu",vmin=-.6,vmax=.6,aspect="auto")
        ax.set_xticks(np.arange(len(wide.columns)),
                      [f"{int(y)-1}–{str(y)[-2:]}" for y in wide.columns])
        ax.set_yticks(np.arange(len(wide)),wide.index)
        ax.set_title(f"{taxon.capitalize()} PRC: within-country annual Spearman ρ vs population change\n"
                     "Countries with ≥8 matched NUTS3 regions; 2014–2020")
        fig.colorbar(img,ax=ax,label="Spearman ρ (annual Δdensity (people/km²) vs ΔPRC)")
        fig.savefig(out/(f"country_annual_correlations_{taxon}.pdf"),bbox_inches="tight")
        plt.close(fig)


def agreement_map(cumulative,meta,taxon,out):
    import geopandas as gpd
    key="cumulative_delta_"+taxon+"_prc"
    d=cumulative[["NUTS_ID","cumulative_delta_population_density",key]].copy()
    d["category"]=classify(d.cumulative_delta_population_density,d[key],
                           *CUTS["cumulative"])
    geo=meta.merge(d,on="NUTS_ID",how="left",validate="one_to_one")
    fig,ax=plt.subplots(figsize=(11,8),layout="constrained")
    geo.plot(ax=ax,color="#bdbdbd",linewidth=0)
    for label,(color,desc) in QUAD.items():
        part=geo[geo.category.eq(label)]
        if len(part):part.plot(ax=ax,color=color,linewidth=0)
    ax.set(xlim=(-12,36),ylim=(35,73),xlabel="Longitude",ylabel="Latitude",
           title=f"Population vs {taxon} PRC 2014–2020: matched sign geography")
    patches=[Patch(facecolor=color,label=desc) for color,desc in QUAD.values()]
    patches.append(Patch(facecolor="#bdbdbd",label="Missing NUTS2016 population / PRC"))
    ax.legend(handles=patches,fontsize=8,loc="lower left",framealpha=.95)
    fig.savefig(out/(f"map_population_prc_direction_{taxon}_2014_2020.pdf"),
                bbox_inches="tight")
    plt.close(fig)


def audit(annual,cumulative,out):
    """Outliers must be investigated, not silently deleted or winsorised."""
    allrows=[]
    for period,frame,col,threshold in (
        ("annual",annual,"delta_population_density",25.),
        ("cumulative",cumulative,"cumulative_delta_population_density",100.)):
        ext=frame.loc[frame[col].abs().gt(threshold),
            ["NUTS_ID","country","year",col,"radiance_band_2014"]
            +["delta_spider_prc","delta_moth_prc"] if period=="annual" else
            ["NUTS_ID","country","year",col,"radiance_band_2014",
             "cumulative_delta_spider_prc","cumulative_delta_moth_prc"]].copy()
        ext["period"]=period
        ext["extreme_cutoff_people_per_km2"]=threshold
        allrows.append(ext)
    pd.concat(allrows,ignore_index=True).to_csv(
        out/"extreme_population_change_review.csv",index=False)


def synthetic_test():
    pop=np.array([2,-3,1,-1,0.0])
    prc=np.array([-.003,.004,.001,-.002,.001])
    cats=classify(pop,prc,.5,.0005)
    assert cats.tolist()==[
      "population_up__prc_down","population_down__prc_up",
      "population_up__prc_up","population_down__prc_down","near_zero_or_mixed"]
    fake=pd.DataFrame({"NUTS_ID":list("ABCDE"),
        "delta_population_density":pop,"delta_spider_prc":prc,
        "cumulative_delta_population_density":pop,
        "cumulative_delta_spider_prc":prc})
    result=directional_summary(fake,"TR","spider","cumulative",.5,.0005)
    assert result["n_growth_prc_loss"]==1
    assert result["n_decline_prc_gain"]==1
    assert result["n_neutral_or_mixed"]==1
    assert np.isnan(result["spearman_rho_pop_vs_delta_prc"])
    print("Country directional quadrant classification and denominator QA: PASS")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--out",type=Path,default=OUT)
    args=ap.parse_args()
    synthetic_test()
    if args.self_test:return
    full,annual,cumulative=calculate(inputs())
    annual,cumulative,meta=countries_by_year(annual,cumulative)
    if cumulative.NUTS_ID.nunique()<1000:
        raise ValueError("Unexpectedly few complete fixed NUTS2016 regions")
    out=args.out if args.out.is_absolute() else ROOT/args.out
    out.mkdir(parents=True,exist_ok=True)
    result=make_tables(annual,cumulative)
    # Side-by-side sensitivity: original population % vs absolute density change.
    # Percentages are retained for audit but are no longer the main analysis.
    comparison=[]
    for period,frame in (("cumulative",cumulative),("annual",annual)):
        intervals=[(2014,2020,frame)] if period=="cumulative" else [
            (yr-1,yr,frame[frame.year.eq(yr)]) for yr in range(2015,2021)]
        for start,end,part in intervals:
            for country,g in [("ALL",part),*list(part.groupby("country"))]:
                for taxon in TAXA:
                    y=("delta_" if period=="annual" else "cumulative_delta_")+taxon+"_prc"
                    for scale,field in (
                        ("absolute_people_per_km2", "delta_population_density" if period=="annual" else "cumulative_delta_population_density"),
                        ("percentage", "pct_population_density" if period=="annual" else "cumulative_pct_population_density")):
                        p,rho=exact_corr(g[field],g[y])
                        comparison.append({"start_year":start,"end_year":end,"country":country,
                          "taxon":taxon,"population_measure":scale,"n_regions":len(g),
                          "pearson_r":p,"spearman_rho":rho})
    result.to_csv(out/"country_directional_summary_annual_and_cumulative.csv",index=False)
    pd.DataFrame(comparison).to_csv(out/"absolute_vs_percentage_correlation_sensitivity.csv",index=False)
    result[(result.period=="cumulative")&(result.threshold_case=="main")].to_csv(
        out/"country_directional_summary_2014_2020.csv",index=False)
    audit(annual,cumulative,out)
    for taxon in TAXA:agreement_map(cumulative,meta,taxon,out)
    plot_country_scatter(cumulative,out/"france_turkiye_iceland_scatter_2014_2020.pdf")
    country_metric_plot(result,out)
    annual_heatmap(result,out)
    print("\nCountry directional 2014–2020 summary (spider):",flush=True)
    cols=["country","n_regions","pearson_r_pop_vs_delta_prc",
          "spearman_rho_pop_vs_delta_prc","n_growth_prc_loss",
          "n_decline_prc_gain","fraction_growth_given_prc_loss",
          "fraction_decline_given_prc_gain"]
    relevant=result[(result.period=="cumulative")&
      (result.threshold_case=="main")&(result.taxon=="spider")]
    print(relevant[relevant.country.isin(["ALL","FR","TR","IS"])][cols].to_string(index=False))
    print("\nResults:",out,flush=True)


if __name__=="__main__":
    main()
