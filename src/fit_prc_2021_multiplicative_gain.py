#!/usr/bin/env python3
"""Direct multiplicative 2021 snow/albedo–VIIRS calibration, NUTS3 held out.

Fits R21 = B * T(background, month, cloud-free) *
               [1 + beta_s*positive_delta_snow + beta_a*positive_delta_albedo]
       + radiance_scale_errors.
B is same-month mean 2020/2022 VIIRS radiance. T is a separate nuisance factor.
Both beta coefficients constrained nonnegative. Model trained only on Jan/Feb.
The operational counterfactual R21/G therefore uses THE SAME G that was fitted.
No 2021 annual data or PRC history is changed. Observational model only.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from project_config import PROJECT_ROOT

FAMILIES = {
    "background_only": (),
    "snow_only": ("snow",),
    "albedo_only": ("albedo",),
    "snow_plus_albedo": ("snow", "albedo"),
}
BASE_NAMES = ("log_background", "february", "cloudfree_10",
              "log_brightness")
MIN_BASELINE = 0.2      # empirical gain unidentified at effectively zero B
HUBER_SCALE_NW = 2.0    # predeclared robust residual scale; sensitivity follows


def prepare(data, min_albedo_valid=0.5):
    required = ["NUTS_ID","month","area_weight","baseline_radiance",
                "delta_radiance","snow_delta_pp","albedo_delta_pp",
                "delta_cloudfree","albedo_valid_pixel_fraction"]
    missing = set(required).difference(data.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    d = data.copy()
    d["NUTS_ID"] = d.NUTS_ID.astype(str)
    for col in required:
        if col != "NUTS_ID":
            d[col] = pd.to_numeric(d[col], errors="coerce")
    good = np.isfinite(d[required[1:]]).all(axis=1)
    d = d.loc[good & d.month.isin([1,2,3]) &
              (d.area_weight > 0) &
              (d.baseline_radiance >= MIN_BASELINE) &
              (d.albedo_valid_pixel_fraction >= min_albedo_valid)].copy()
    d["radiance_2021"] = d.baseline_radiance + d.delta_radiance
    d = d.loc[d.radiance_2021 >= 0].copy()
    d["snow_pos"] = np.maximum(d.snow_delta_pp, 0) / 100.0
    d["albedo_pos"] = np.maximum(d.albedo_delta_pp, 0) / 100.0
    d["brightness_stratum"] = pd.cut(
        d.baseline_radiance, bins=[MIN_BASELINE-1e-6,1,5,20,np.inf],
        labels=["0.2–1","1–5","5–20",">=20"], include_lowest=True)
    return d


def region_month_weights(d):
    area = d.area_weight.to_numpy(float)
    sums = d.groupby(["NUTS_ID","month"]).area_weight.transform("sum").to_numpy(float)
    return (area / sums) * (len(d) / len(d.groupby(["NUTS_ID","month"])))


def model_data(d):
    return {
        "B":d.baseline_radiance.to_numpy(float),
        "February":d.month.eq(2).to_numpy(float),
        "CF":d.delta_cloudfree.to_numpy(float) / 10.0,
        "L":np.log1p(d.baseline_radiance.to_numpy(float)) / np.log(11.0),
        "snow":d.snow_pos.to_numpy(float),
        "albedo":d.albedo_pos.to_numpy(float),
        "target":d.radiance_2021.to_numpy(float),
    }


def fit_predict(theta, arrays, family, jacobian=False):
    """Direct product (not a linearized anomaly regression)."""
    n_background = len(BASE_NAMES)
    logt = (theta[0] + theta[1]*arrays["February"] +
            theta[2]*arrays["CF"] + theta[3]*arrays["L"])
    # Keep finite numerical representation even if optimizer visits a bound.
    T = np.exp(np.clip(logt, -20, 20))
    gain = np.ones_like(T)
    for pos, term in enumerate(FAMILIES[family], n_background):
        gain = gain + theta[pos]*arrays[term]
    pred = arrays["B"]*T*gain
    if not jacobian:
        return pred, T, gain
    J = np.column_stack([pred, pred*arrays["February"],
                         pred*arrays["CF"], pred*arrays["L"]] +
                        [arrays["B"]*T*arrays[term]
                         for term in FAMILIES[family]])
    return J


def fit_multiplicative(d, family, f_scale=HUBER_SCALE_NW):
    if d.NUTS_ID.nunique() < 6 or len(d) < 40:
        raise ValueError("Too few independent regions/groups")
    arr = model_data(d)
    rootw = np.sqrt(region_month_weights(d))
    n = len(BASE_NAMES)+len(FAMILIES[family])
    lower = np.array([-2.0,-1.2,-1.2,-1.2] + [0.0]*(n-4))
    upper = np.array([2.0, 1.2, 1.2, 1.2] + [15.0]*(n-4))
    initial = np.array([0.,0.,0.,0.] + [1.0]*(n-4))
    def fun(theta):
        pred,_,_ = fit_predict(theta,arr,family)
        return rootw*(pred-arr["target"])
    def jac(theta):
        return rootw[:,None]*fit_predict(theta,arr,family,jacobian=True)
    solved = least_squares(fun, initial, jac=jac, bounds=(lower,upper),
                           loss="soft_l1", f_scale=f_scale,
                           max_nfev=350, x_scale="jac", ftol=1e-9, xtol=1e-9)
    if not solved.success:
        raise RuntimeError(f"{family}: optimizer did not converge: {solved.message}")
    params = dict(zip([*BASE_NAMES,
                      *["beta_"+s for s in FAMILIES[family]]],solved.x))
    for term in FAMILIES[family]:
        if params["beta_"+term]<0:
            raise AssertionError("Gain coefficient violates nonnegative bound")
    return params, float(solved.cost), int(solved.nfev)


def evaluate(d, params, family):
    out = d.copy()
    theta = np.array([params[k] for k in BASE_NAMES] +
                     [params["beta_"+t] for t in FAMILIES[family]])
    arr = model_data(out)
    pred,T,G = fit_predict(theta,arr,family)
    corrected = arr["target"] / G
    if np.any(G < 1-1e-12) or np.any(corrected > arr["target"]+1e-9):
        raise AssertionError("Snow gain cannot increase radiance on correction")
    out["model"] = family
    out["background_factor"] = T
    out["surface_gain"] = G
    out["predicted_2021"] = pred
    out["radiance_surface_corrected"] = corrected
    out["removed_radiance"] = arr["target"]-corrected
    out["error_observed_prediction"] = arr["target"]-pred
    out["error_flank_before"] = arr["target"]-arr["B"]
    out["error_flank_after"] = corrected-arr["B"]
    out["is_surface_control"] = (out.snow_pos <= 0) & (out.albedo_pos <= 0)
    return out


def diagnostics(d):
    rows=[]
    for (family,month), g in d.groupby(["model","month"]):
        w = region_month_weights(g)
        def wrmse(col):
            return float(np.sqrt(np.average(g[col].to_numpy(float)**2,weights=w)))
        reg=[]
        for _,x in g.groupby("NUTS_ID"):
            before=np.average(x.error_flank_before,weights=x.area_weight)
            after=np.average(x.error_flank_after,weights=x.area_weight)
            reg.append(abs(after)<abs(before))
        rows.append({
            "model":family,"month":int(month),
            "n_regions":g.NUTS_ID.nunique(),"n_groups":len(g),
            "rmse_2021_prediction":wrmse("error_observed_prediction"),
            "rmse_to_flanks_before":wrmse("error_flank_before"),
            "rmse_to_flanks_after":wrmse("error_flank_after"),
            "region_fraction_improved":float(np.mean(reg)),
            "fraction_gain_gt2":float(np.mean(g.surface_gain>2)),
            "max_gain":float(g.surface_gain.max()),
            "median_gain":float(g.surface_gain.median()),
            "fraction_surface_control":float(g.is_surface_control.mean()),
            "median_radiance_removed":float(g.removed_radiance.median()),
        })
    return pd.DataFrame(rows)


def holdout_parameter_ranges(coefs):
    d = coefs[coefs.scope.eq("held_out_fold")]
    rows=[]
    for (family,param), g in d.groupby(["model","parameter"]):
        v=g.value.to_numpy(float)
        rows.append({"model":family,"parameter":param,"n_folds":len(v),
                     "min":v.min(),"p10":np.quantile(v,.1),
                     "median":np.median(v),"p90":np.quantile(v,.9),
                     "max":v.max()})
    return pd.DataFrame(rows)


def stratum_diagnostics(d):
    rows=[]
    for (family,month,bright),g in d.groupby(
            ["model","month","brightness_stratum"],observed=True):
        w=region_month_weights(g)
        rows.append({"model":family,"month":int(month),"brightness_stratum":str(bright),
                     "n_regions":g.NUTS_ID.nunique(),"n_groups":len(g),
                     "mean_gain":np.average(g.surface_gain,weights=w),
                     "mean_radiance_removed":np.average(g.removed_radiance,weights=w),
                     "rmse_flank_before":np.sqrt(np.average(g.error_flank_before**2,weights=w)),
                     "rmse_flank_after":np.sqrt(np.average(g.error_flank_after**2,weights=w))})
    return pd.DataFrame(rows)


def unit_tests():
    rng=np.random.default_rng(39)
    n=720
    B=rng.uniform(.4,30,n)
    month=np.tile(np.repeat([1,2],20),18)
    snow=rng.uniform(-.2,.6,n)
    alb=rng.uniform(-.08,.3,n)
    cf=rng.normal(0,3,n)
    df=pd.DataFrame({"NUTS_ID":np.repeat([f"T{i:02d}" for i in range(18)],40),
        "month":month,"area_weight":np.ones(n),
        "baseline_radiance":B,"snow_delta_pp":snow*100,
        "albedo_delta_pp":alb*100,"delta_cloudfree":cf,
        "albedo_valid_pixel_fraction":np.ones(n)})
    T=np.exp(.04+.03*(month==2)+.02*cf/10+.01*np.log1p(B)/np.log(11))
    G=1+.8*np.maximum(snow,0)+2.4*np.maximum(alb,0)
    df["delta_radiance"]=B*T*G-B
    pre=prepare(df)
    params,_,_=fit_multiplicative(pre,"snow_plus_albedo",f_scale=1.0)
    assert abs(params["beta_snow"]-.8)<.025,params
    assert abs(params["beta_albedo"]-2.4)<.025,params
    val=evaluate(pre,params,"snow_plus_albedo")
    assert np.all(val.radiance_surface_corrected<=val.radiance_2021+1e-10)
    z=val[val.is_surface_control]
    if len(z):
        assert np.allclose(z.removed_radiance,0)
    print("Direct-multiplicative synthetic recovery and sign test: PASS",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups",default=
        "data/processed/historical_prc_2021_pixel_snow_test/"
        "matched_snowcell_brightness_groups.csv")
    parser.add_argument("--output-dir",default=
        "data/processed/historical_prc_2021_multiplicative_calibration")
    parser.add_argument("--min-albedo-valid-fraction",type=float,default=.5)
    parser.add_argument("--robust-scale-nw",type=float,default=HUBER_SCALE_NW)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    unit_tests()
    if args.self_test:
        return
    if args.robust_scale_nw<=0:
        raise ValueError("Robust scale must be positive")
    p=Path(args.groups)
    if not p.is_absolute():p=PROJECT_ROOT/p
    if not p.exists():raise FileNotFoundError(p)
    d=prepare(pd.read_csv(p,dtype={"snow_cell_id":str}),
              args.min_albedo_valid_fraction)
    if d.NUTS_ID.nunique()<8:
        raise ValueError("Need >=8 independent regions")
    print(f"Matched QA+brightness cohort: {len(d)} groups in {d.NUTS_ID.nunique()} regions",
          flush=True)
    winter=d[d.month.isin([1,2])]
    coefficients, predictions=[],[]
    for family in FAMILIES:
        full,cost,nfev=fit_multiplicative(winter,family,args.robust_scale_nw)
        coefficients.extend(
            {"model":family,"scope":"all_Jan_Feb","held_out_region":"",
             "parameter":key,"value":v,"cost":cost,"nfev":nfev}
            for key,v in full.items())
        for held in sorted(d.NUTS_ID.unique()):
            subset=d[d.NUTS_ID.ne(held)&d.month.isin([1,2])]
            estimates,fit_cost,fit_steps=fit_multiplicative(
                subset,family,args.robust_scale_nw)
            coefficients.extend(
                {"model":family,"scope":"held_out_fold","held_out_region":held,
                 "parameter":key,"value":v,"cost":fit_cost,"nfev":fit_steps}
                for key,v in estimates.items())
            pred=evaluate(d[d.NUTS_ID.eq(held)],estimates,family)
            predictions.append(pred)
        print(f"{family}: "+", ".join(f"{k}={v:.5f}" for k,v in full.items())+
              f"; optimization calls {nfev}",flush=True)
    cv=pd.concat(predictions,ignore_index=True)
    coefficients=pd.DataFrame(coefficients)
    diag=diagnostics(cv)
    ranges=holdout_parameter_ranges(coefficients)
    strata=stratum_diagnostics(cv)
    reg=[]
    for (model,month,nuts),g in cv.groupby(["model","month","NUTS_ID"]):
        w=g.area_weight.to_numpy(float)
        reg.append({"model":model,"month":int(month),"NUTS_ID":nuts,
                    "n_groups":len(g),"area_sum":sum(w),
                    "radiance_observed":np.average(g.radiance_2021,weights=w),
                    "radiance_snow_adjusted":np.average(g.radiance_surface_corrected,weights=w),
                    "radiance_2020_22_mean":np.average(g.baseline_radiance,weights=w),
                    "model_predicted_2021":np.average(g.predicted_2021,weights=w),
                    "mean_gain":np.average(g.surface_gain,weights=w),
                    "fraction_groups_gain_gt2":float(np.mean(g.surface_gain>2))})
    out=Path(args.output_dir)
    if not out.is_absolute():out=PROJECT_ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    coefficients.to_csv(out/"multiplicative_coefficients_and_folds.csv",index=False)
    diag.to_csv(out/"multiplicative_holdout_diagnostics.csv",index=False)
    ranges.to_csv(out/"multiplicative_coefficient_stability.csv",index=False)
    strata.to_csv(out/"multiplicative_brightness_strata.csv",index=False)
    pd.DataFrame(reg).to_csv(out/"multiplicative_region_month.csv",index=False)
    cv.to_csv(out/"multiplicative_heldout_groups.csv",index=False)
    fig,axes=plt.subplots(1,2,figsize=(11,5),layout="constrained")
    for family in ("snow_only","albedo_only","snow_plus_albedo"):
        parts=diag[diag.model.eq(family)].sort_values("month")
        axes[0].plot(parts.month,parts.rmse_to_flanks_after,marker="o",label=family)
        axes[1].plot(parts.month,parts.rmse_2021_prediction,marker="o",label=family)
    axes[0].set(title="Surface-adjusted anomaly, held-out regions",
                xlabel="Month",ylabel="RMSE vs flanks (nW cm⁻² sr⁻¹)")
    axes[1].set(title="Direct multiplicative model prediction",
                xlabel="Month",ylabel="Held-out 2021 prediction RMSE")
    for ax in axes:
        ax.set_xticks([1,2,3],["January","February","March"])
        ax.grid(alpha=.2)
        ax.legend()
    fig.savefig(out/"multiplicative_holdout_diagnostics.pdf")
    plt.close(fig)
    text=[
        "DIRECT MULTIPLICATIVE 2021 SURFACE GAIN CALIBRATION (NO ANNUAL REWRITE)",
        "Direct training: observed_R21 = B * exp(gamma0+gammaFeb*I_Feb+",
        " gammaCF*(delta_cf/10)+gammaBrightness*log1p(B)/log(11)) * G + error",
        "G=1+betaSnow*max(deltaSnow,0)+betaAlbedo*max(deltaAlbedo,0), beta>=0",
        "B = matched-calendar-month mean of 2020 and 2022 EOG monthly VIIRS.",
        "Snow and albedo anomalies in fractional units, not 0-100 percentages.",
        "Apply the *same fitted G*: R21_standardized=R21/G, no radiance increase.",
        f"Model calibration requires baseline >= {MIN_BASELINE} nW cm^-2 sr^-1.",
        "Nuisance factor T is NOT removed by the surface-only correction.",
        "Raw 2021 brightness vs flanks is a validation, not a fit selection target.",
        "Train Jan/Feb only; leave whole NUTS3 out; March is temporal control.",
        "Region-month weights each normalized to equal contribution; within area.",
        "NUTS-disjoint holdouts still reuse 2020/22 endpoints as observational B.",
        "No change to annual 0.235 baseline correction or existing 2021 PRC.",
        f"Input groups={len(d)}, regions={d.NUTS_ID.nunique()}, Jan+Feb={len(winter)}",
        "",
        "Full Jan/Feb coefficients:",
        coefficients[coefficients.scope.eq("all_Jan_Feb")].to_string(index=False),
        "",
        "Held-out coefficients by parameter:",
        ranges.to_string(index=False),
        "",
        "Held-out diagnostics:",
        diag.to_string(index=False),
        "",
        "IMPORTANT: strong G, sharp coefficient shifts, winter March placebo",
        "errors or high-brightness leverage prevent acceptance as causal bias.",
    ]
    (out/"multiplicative_summary.txt").write_text("\n".join(text)+"\n",encoding="utf-8")
    print("\n".join(text))
    print("Outputs:",out)


if __name__=="__main__":
    main()
