#!/usr/bin/env python3
"""Measured annual demography vs annual radiance and both taxa's PRC.

Use fixed NUTS2016 measured Eurostat density 2014–2020: six 1-year changes.
Compare nested models (country x transition-year FE and prior-year radiance /
density, with vs without actual annual population growth). Inference uses
NUTS3-clustered SE. Results are associations, NOT proven causal effects.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from project_config import PROJECT_ROOT as ROOT

BASE = "C(country_year) + log_lag_r + I(log_lag_r ** 2) + log_lag_pop"


def load():
    processed = ROOT / "data/processed"
    pop = pd.read_parquet(processed / "nuts3_population_density.parquet")
    rad = pd.read_parquet(processed / "fixed_nuts2016_viirs_2014_2024.parquet")
    ref = pd.read_csv(processed / "fixed_nuts2016_population_density_reference.csv")
    pop = pop[pop.year.between(2014, 2020) & pop.nuts_release.eq(2016)].copy()
    if "matches_nuts_release" in pop:
        pop = pop[pop.matches_nuts_release.fillna(False)]
    pop = pop[["NUTS_ID", "year", "population_density"]]
    rad = rad[rad.year.between(2014, 2020)][
        ["NUTS_ID", "year", "CNTR_CODE", "radiance_mean_corrected"]].copy()
    panel = rad.merge(pop, on=["NUTS_ID","year"], validate="one_to_one")
    for taxon in ("spider", "moth"):
        f = (processed / f"historical_prc_2014_2024_{taxon}_snow_corrected"
             / "prc_nuts3_annual_cumulative.csv")
        t = pd.read_csv(f, dtype={"NUTS_ID":str})
        if "snow_status" not in t or not t.snow_status.eq("locked").all():
            raise ValueError(f"{taxon}: needs locked snow-corrected PRC panel")
        t = t[t.year.between(2014,2020)][
            ["NUTS_ID","year","prc_mean","coverage_stable"]]
        t = t.rename(columns={"prc_mean":f"{taxon}_prc",
                              "coverage_stable":f"{taxon}_stable"})
        panel = panel.merge(t, on=["NUTS_ID","year"], validate="one_to_one")
    panel = panel.merge(ref[["NUTS_ID","density_band"]],
                        on="NUTS_ID",validate="many_to_one")
    return panel


def prepare(panel):
    d = panel.sort_values(["NUTS_ID","year"]).copy()
    if d.duplicated(["NUTS_ID","year"]).any():
        raise ValueError("Duplicate NUTS3/year")
    for key in ("population_density","radiance_mean_corrected",
                "spider_prc","moth_prc"):
        d["lag_"+key] = d.groupby("NUTS_ID")[key].shift()
        d["d_"+key] = d[key] - d["lag_"+key]
    d["previous_year"] = d.groupby("NUTS_ID").year.shift()
    d = d[d.year.sub(d.previous_year).eq(1) & d.year.between(2015,2020)]
    d = d[(d.population_density>0)&(d.lag_population_density>0)&
          (d.radiance_mean_corrected>=0)&(d.lag_radiance_mean_corrected>=0)&
          d.spider_stable & d.moth_stable].copy()
    d["pop_change_pct"] = 100*(d.population_density/d.lag_population_density-1)
    d["log_lag_r"] = np.log1p(d.lag_radiance_mean_corrected)
    d["log_lag_pop"] = np.log(d.lag_population_density)
    d["country_year"] = d.CNTR_CODE.astype(str)+"-"+d.year.astype(str)
    keep = d.groupby("CNTR_CODE").NUTS_ID.nunique()
    d = d[d.CNTR_CODE.isin(keep[keep>=5].index)]
    needed = ["pop_change_pct","log_lag_r","log_lag_pop",
              "d_radiance_mean_corrected","d_spider_prc","d_moth_prc"]
    d = d.replace([np.inf,-np.inf],np.nan).dropna(subset=needed)
    return d.reset_index(drop=True)


def compare(d,response):
    baseline = smf.ols(f"{response} ~ {BASE}",data=d).fit()
    full = smf.ols(f"{response} ~ {BASE} + pop_change_pct",
                   data=d).fit(cov_type="cluster",
                               cov_kwds={"groups":d.NUTS_ID})
    b = float(full.params["pop_change_pct"])
    se = float(full.bse["pop_change_pct"])
    residual_base = np.square(d[response]-baseline.predict(d)).sum()
    residual_full = np.square(d[response]-full.predict(d)).sum()
    row = {"response":response,"n_region_year_transitions":len(d),
           "n_distinct_regions":int(d.NUTS_ID.nunique()),
           "beta_per_population_growth_percent":b,
           "beta_ci_low":b-1.96*se,"beta_ci_high":b+1.96*se,
           "cluster_robust_p":float(full.pvalues["pop_change_pct"]),
           "baseline_R2":float(baseline.rsquared),
           "population_R2":float(full.rsquared),
           "incremental_total_R2":float(full.rsquared-baseline.rsquared),
           "fraction_of_residual_variance_explained":
             float((residual_base-residual_full)/residual_base)}
    # This component is a conditional statistical association, not a
    # validated 'number of animals harmed by population growth'.
    d[response+"_population_association"] = b*d.pop_change_pct
    d[response+"_fitted_unchanged_population"] = (
        full.predict(d)-b*d.pop_change_pct)
    return row


def synthetic_test():
    records=[]
    for i in range(10):
        for year in range(2014,2021):
            records.append({"NUTS_ID":f"A{i}","year":year,
                "CNTR_CODE":"DE","population_density":100+i+year-2014,
                "radiance_mean_corrected":1+i/100,
                "spider_prc":.42,"moth_prc":.66,
                "spider_stable":True,"moth_stable":True})
    result=prepare(pd.DataFrame(records))
    assert len(result)==60 and result.year.min()==2015
    assert np.isfinite(result.pop_change_pct).all()
    print("Annual population and PRC pairing self-test: PASS",flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    args = ap.parse_args()
    synthetic_test()
    if args.self_test:return
    d=prepare(load())
    if d.NUTS_ID.nunique()<200:
        raise ValueError("Too few matched fixed-boundary regions")
    out=ROOT/"data/processed/prc_population_annual_change"
    out.mkdir(parents=True,exist_ok=True)
    results=[compare(d,col) for col in (
        "d_radiance_mean_corrected","d_spider_prc","d_moth_prc")]
    pd.DataFrame(results).to_csv(out/"population_explanatory_power.csv",index=False)
    d.to_csv(out/"annual_matched_population_radiance_spider_moth.csv",index=False)
    (out/"model_summary.json").write_text(json.dumps({
        "years":"2015-2020 annual transitions",
        "population_basis":"NUTS2016 Eurostat observed 2014-2020",
        "models":results,
        "warning":"Population component is observational association, not causal attribution."
    },indent=2)+"\n")
    fig,ax=plt.subplots(figsize=(8,5))
    ax.bar(["Radiance","Spider PRC","Moth PRC"],
           [x["fraction_of_residual_variance_explained"]*100 for x in results])
    ax.set_ylabel("Additional unexplained variance accounted for (%)")
    ax.set_title("Annual population change: incremental association (2015-2020)")
    fig.savefig(out/"annual_population_incremental_explanation.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print("Done:",out,flush=True)


if __name__=="__main__":
    main()
