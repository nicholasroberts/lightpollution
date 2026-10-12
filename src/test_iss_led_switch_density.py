#!/usr/bin/env python3
"""Test whether dense-area VIIRS decline is consistent with LED spectral switching.

Uses paired ISS colour-ratio change and the current production-calibrated
2013/2020 NUTS3 VIIRS endpoints. The primary mechanistic conversion is from
Sanchez de Miguel et al. (2022):

    G_ISS / VIIRS = 0.21 * 1.5 ** (G_ISS / R_ISS)

Thus, at equal green-band output:

    predicted spectral delta log10(VIIRS)
        = -log10(1.5) * delta(G/R)

Positive delta G/R therefore predicts a negative VIIRS response component.
The ISS post mosaic is a 2014-2020 composite, so this is an endpoint mechanism
analysis rather than a year-resolved LED-adoption model.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from project_config import configured_path, load_config

LOG10_1P5=float(np.log10(1.5))

def load_data(config):
    cfg=config["analysis"]["iss_led_switch_density_test"]
    processed=configured_path(config,"processed")
    spectral=pd.read_csv(processed/"iss_nuts3_spectral_shift.csv")
    merged=pd.read_parquet(processed/config["outputs"]["merged"]).copy()

    need_s=["NUTS_ID","CNTR_CODE","n_common_valid_pixels","delta_gr_median","delta_bg_median"]
    need_m=["year","NUTS_ID","CNTR_CODE","nuts_release","radiance_mean_corrected","population_density"]
    for label,frame,cols in [("spectral",spectral,need_s),("merged",merged,need_m)]:
        miss=[c for c in cols if c not in frame.columns]
        if miss: raise KeyError(f"{label} missing columns: {miss}")

    pre=int(cfg["endpoint_pre_year"]); post=int(cfg["endpoint_post_year"]); rel=int(cfg["nuts_release"])
    e=merged[(merged.year.isin([pre,post]))&(merged.nuts_release==rel)][need_m].copy()
    a=e[e.year==pre].drop(columns=["year","nuts_release"]).rename(columns={
        "radiance_mean_corrected":"radiance_pre","population_density":"population_density_pre"})
    b=e[e.year==post].drop(columns=["year","nuts_release"]).rename(columns={
        "radiance_mean_corrected":"radiance_post","population_density":"population_density_post"})
    d=a.merge(b,on=["NUTS_ID","CNTR_CODE"],validate="one_to_one").merge(
        spectral,on=["NUTS_ID","CNTR_CODE"],validate="one_to_one")

    ok=(d.radiance_pre>0)&(d.radiance_post>0)&(d.population_density_pre>0)&(d.population_density_post>0)
    ok &= d.n_common_valid_pixels>=int(cfg["min_valid_iss_pixels"])
    for c in ["delta_gr_median","delta_bg_median"]:
        ok &= np.isfinite(pd.to_numeric(d[c],errors="coerce"))
    d=d.loc[ok].copy()

    d["log_rad_pre"]=np.log10(d.radiance_pre)
    d["log_rad_post"]=np.log10(d.radiance_post)
    d["delta_log_rad_observed"]=d.log_rad_post-d.log_rad_pre
    d["log_pop_pre"]=np.log10(d.population_density_pre)
    d["delta_log_pop"]=np.log10(d.population_density_post)-d.log_pop_pre
    d["predicted_spectral_delta_log_rad"]=-LOG10_1P5*d.delta_gr_median
    d["delta_log_rad_spectral_adjusted"]=d.delta_log_rad_observed-d.predicted_spectral_delta_log_rad
    d["log_pop_within_country"]=d.log_pop_pre-d.groupby("CNTR_CODE").log_pop_pre.transform("mean")
    d["spectral_delta_within_country"]=d.predicted_spectral_delta_log_rad-d.groupby("CNTR_CODE").predicted_spectral_delta_log_rad.transform("mean")
    d["delta_gr_within_country"]=d.delta_gr_median-d.groupby("CNTR_CODE").delta_gr_median.transform("mean")
    d["delta_bg_within_country"]=d.delta_bg_median-d.groupby("CNTR_CODE").delta_bg_median.transform("mean")

    edges=[-np.inf]+list(cfg["density_band_breaks"])+[np.inf]
    d["density_band"]=pd.cut(d.population_density_pre,bins=edges,labels=cfg["density_band_labels"],right=False,ordered=True)
    cutoff=float(cfg["high_density_cutoff"])
    d["dense_vs_sparse"]=np.where(d.population_density_pre>=cutoff,1,np.where(d.population_density_pre<100,0,np.nan))

    counts=d.groupby("CNTR_CODE").NUTS_ID.nunique()
    d=d[d.CNTR_CODE.isin(counts[counts>=int(cfg["min_country_regions"])].index)].copy()
    d["CNTR_CODE"]=d.CNTR_CODE.astype("category")
    return d

def density_summary(d,config):
    cfg=config["analysis"]["iss_led_switch_density_test"]
    years=int(cfg["endpoint_post_year"])-int(cfg["endpoint_pre_year"])
    cb=d.groupby(["CNTR_CODE","density_band"],observed=True,as_index=False).agg(
        n_nuts3=("NUTS_ID","nunique"),delta_gr=("delta_gr_median","mean"),delta_bg=("delta_bg_median","mean"),
        observed=("delta_log_rad_observed","mean"),spectral=("predicted_spectral_delta_log_rad","mean"),
        adjusted=("delta_log_rad_spectral_adjusted","mean"))
    rows=[]
    for band in cfg["density_band_labels"]:
        g=cb[cb.density_band.astype(str)==band]
        if g.empty: continue
        r={"density_band":band,"n_countries":g.CNTR_CODE.nunique(),"n_nuts3":int(g.n_nuts3.sum()),
           "mean_country_delta_gr":g.delta_gr.mean(),"mean_country_delta_bg":g.delta_bg.mean()}
        for label,col in [("observed","observed"),("spectral_predicted","spectral"),("spectral_adjusted","adjusted")]:
            x=float(g[col].mean())
            r[label+"_delta_log10"]=x
            r[label+"_total_pct"]=100*(10**x-1)
            r[label+"_annualized_pct"]=100*(10**(x/years)-1)
        rows.append(r)
    return pd.DataFrame(rows)

def fit_models(d,config):
    rows=[]
    for metric in ["delta_gr_median","delta_bg_median"]:
        f=smf.ols(f"{metric} ~ C(CNTR_CODE) + log_pop_within_country",d).fit(cov_type="HC3")
        rows.append({"model":metric+"_density_gradient","subset":"all","term":"log_pop_within_country",
                     "estimate":f.params["log_pop_within_country"],"se_hc3":f.bse["log_pop_within_country"],
                     "p_hc3":f.pvalues["log_pop_within_country"],"n":int(f.nobs)})

    base="delta_log_rad_observed ~ C(CNTR_CODE) + delta_log_pop + log_pop_within_country"
    full=base+" + spectral_delta_within_country"
    f=smf.ols(full,d).fit(cov_type="HC3")
    plain0=smf.ols(base,d).fit(); plain1=smf.ols(full,d).fit(); a=anova_lm(plain0,plain1).iloc[-1]
    rows.append({"model":"spectral_component","subset":"all","term":"spectral_delta_within_country",
                 "estimate":f.params["spectral_delta_within_country"],"se_hc3":f.bse["spectral_delta_within_country"],
                 "p_hc3":f.pvalues["spectral_delta_within_country"],"n":int(f.nobs),
                 "block_F":a["F"],"block_p":a["Pr(>F)"],"delta_r2":plain1.rsquared-plain0.rsquared})

    cutoff=float(config["analysis"]["iss_led_switch_density_test"]["high_density_cutoff"])
    dense=d[d.population_density_pre>=cutoff].copy()
    if len(dense)>=20:
        b0="delta_log_rad_observed ~ C(CNTR_CODE) + delta_log_pop"
        b1=b0+" + spectral_delta_within_country"
        f=smf.ols(b1,dense).fit(cov_type="HC3"); p0=smf.ols(b0,dense).fit(); p1=smf.ols(b1,dense).fit(); a=anova_lm(p0,p1).iloc[-1]
        rows.append({"model":"spectral_component","subset":f"density>={cutoff:g}","term":"spectral_delta_within_country",
                     "estimate":f.params["spectral_delta_within_country"],"se_hc3":f.bse["spectral_delta_within_country"],
                     "p_hc3":f.pvalues["spectral_delta_within_country"],"n":int(f.nobs),
                     "block_F":a["F"],"block_p":a["Pr(>F)"],"delta_r2":p1.rsquared-p0.rsquared})

    c=d[np.isfinite(d.dense_vs_sparse)].copy()
    for outcome,label in [("delta_log_rad_observed","observed"),("delta_log_rad_spectral_adjusted","spectral_adjusted")]:
        form=f"{outcome} ~ C(CNTR_CODE) + delta_log_pop + dense_vs_sparse"
        f=smf.ols(form,c).fit(cov_type="HC3")
        rows.append({"model":"dense_vs_sparse_contrast","subset":label,"term":"dense_vs_sparse",
                     "estimate":f.params["dense_vs_sparse"],"se_hc3":f.bse["dense_vs_sparse"],
                     "p_hc3":f.pvalues["dense_vs_sparse"],"n":int(f.nobs)})
    return pd.DataFrame(rows)

def make_plot(s, path, dpi):
    x=np.arange(len(s)); fig,ax=plt.subplots(figsize=(10,6))
    ax.plot(x,s.observed_annualized_pct,marker="o",label="Observed calibrated VIIRS")
    ax.plot(x,s.spectral_predicted_annualized_pct,marker="s",linestyle="--",label="ISS G/R spectral component")
    ax.plot(x,s.spectral_adjusted_annualized_pct,marker="D",linestyle="-.",label="After removing spectral component")
    ax.axhline(0,color="0.4",linewidth=1); ax.set_xticks(x); ax.set_xticklabels(s.density_band.str.replace(">=","≥",regex=False))
    ax.set_xlabel("2013 population-density band (people km$^{-2}$)"); ax.set_ylabel("2013→2020 annualized change (%)")
    ax.set_title("ISS LED-switch decomposition of calibrated VIIRS endpoint change"); ax.grid(axis="y",alpha=.2); ax.legend()
    fig.tight_layout(); fig.savefig(path,dpi=dpi,bbox_inches="tight"); plt.close(fig)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",default="config/pipeline.yaml"); args=ap.parse_args()
    config=load_config(args.config); cfg=config["analysis"]["iss_led_switch_density_test"]
    d=load_data(config); s=density_summary(d,config); m=fit_models(d,config)
    processed=configured_path(config,"processed"); figs=configured_path(config,"figures")/"population_relationship"/"spectral_validation"/"led_switch_density"
    processed.mkdir(parents=True,exist_ok=True); figs.mkdir(parents=True,exist_ok=True)
    d.to_csv(processed/"iss_led_switch_density_nuts3.csv",index=False)
    s.to_csv(processed/"iss_led_switch_density_band_summary.csv",index=False)
    m.to_csv(processed/"iss_led_switch_density_models.csv",index=False)
    make_plot(s,figs/"density_band_led_switch_decomposition.png",int(cfg["figure_dpi"]))
    out=processed/"iss_led_switch_density_summary.txt"
    with out.open("w") as h:
        h.write("ISS COLOUR RATIOS AND THE DENSE-AREA LED-SWITCH HYPOTHESIS\n"+"="*78+"\n\n")
        h.write(f"NUTS3={len(d)}; countries={d.CNTR_CODE.nunique()}\n")
        h.write("ISS post mosaic is a 2014-2020 composite; 2020 VIIRS is an endpoint anchor.\n")
        h.write("predicted spectral delta log10(VIIRS) = -log10(1.5) * paired delta(G/R), at equal green output.\n\n")
        h.write("COUNTRY-BALANCED DENSITY-BAND DECOMPOSITION\n")
        h.write(s.to_string(index=False,float_format=lambda x:f"{x:.6g}"))
        h.write("\n\nMODELS\n")
        h.write(m.to_string(index=False,float_format=lambda x:f"{x:.6g}"))
        h.write("\n\nINTERPRETATION\n")
        h.write("- Positive delta G/R is whitening and predicts a negative VIIRS spectral component.\n")
        h.write("- If whitening increases with density and the dense-vs-sparse contrast shrinks after spectral adjustment, LED switching plausibly contributes to dense-area VIIRS decline.\n")
        h.write("- This is not proof that LED switching explains all decline; equal-green output is a spectral normalization, not measured constant lumens or power.\n")
    print(s.to_string(index=False,float_format=lambda x:f"{x:.6f}"))
    print("\nMODELS\n"+m.to_string(index=False,float_format=lambda x:f"{x:.6g}"))
    print("\nWrote:",out)

if __name__=="__main__": main()
