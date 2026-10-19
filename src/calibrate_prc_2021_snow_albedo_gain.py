#!/usr/bin/env python3
"""Estimate a monotone surface-state radiance amplification from matched VIIRS.

Input is the independently assembled 2020/21/22 NUTS3 × month × 0.05-degree
snow cell × brightness-class pilot. Fits January/February with NUTS3-disjoint
validation and evaluates March as a negative temporal control.

The gain is G=1+beta_s*max(dSnow,0)+beta_a*max(dAlbedo,0),
where both differences are FRACTIONS and both betas are >=0.
The candidate surface-only correction R2021/G never increases VIIRS radiance.
No annual PRC or VIIRS files are modified. Not a causal calibration.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

from project_config import PROJECT_ROOT

MODELS = {
    "snow_only": ("snow",),
    "albedo_only": ("albedo",),
    "snow_plus_albedo": ("snow", "albedo"),
}
NUISANCE = ("intercept", "february", "cloudfree_change", "brightness_trend")


def response_matrix(d, terms):
    """Linear in fitted coefficients; multiplicative surface hypothesis."""
    b = d.baseline_radiance.to_numpy(dtype=float)
    x = [np.ones(len(d)), d.month.eq(2).to_numpy(dtype=float),
         d.delta_cloudfree.to_numpy(dtype=float), b]
    names = list(NUISANCE)
    if "snow" in terms:
        x.append(b * d.snow_increase_fraction.to_numpy(dtype=float))
        names.append("beta_snow")
    if "albedo" in terms:
        x.append(b * d.albedo_increase_fraction.to_numpy(dtype=float))
        names.append("beta_albedo")
    return np.column_stack(x), names


def weights_for_region_month(d):
    """Each training region-month has equal total weight, cells area-weighted."""
    area = d.area_weight.to_numpy(dtype=float)
    totals = d.groupby(["NUTS_ID", "month"]).area_weight.transform("sum").to_numpy(float)
    w = area / totals
    return w / np.mean(w)


def constrained_fit(d, terms, robust_iterations=5):
    if len(d) < 20 or d.NUTS_ID.nunique() < 6:
        raise ValueError("Need >=6 training regions and >=20 groups")
    x, names = response_matrix(d, terms)
    y = d.delta_radiance.to_numpy(dtype=float)
    wbase = weights_for_region_month(d)
    bounds_low = np.full(len(names), -np.inf, dtype=float)
    for name in ("beta_snow", "beta_albedo"):
        if name in names:
            bounds_low[names.index(name)] = 0.0
    bounds_high = np.full(len(names), np.inf, dtype=float)
    w = wbase.copy()
    for _ in range(robust_iterations):
        sw = np.sqrt(w)
        fit = lsq_linear(x * sw[:, None], y * sw,
                         bounds=(bounds_low, bounds_high),
                         max_iter=400, tol=1e-11)
        if not fit.success:
            raise RuntimeError("Constrained surface calibration did not converge: "+fit.message)
        r = y - x @ fit.x
        centre = np.median(r)
        sigma = max(1e-6, 1.4826 * np.median(np.abs(r - centre)))
        # Fixed Huber-like IRLS weights, same predeclared robust policy in each fold.
        robust = np.minimum(1.0, 1.345 * sigma / np.maximum(np.abs(r-centre), 1e-12))
        w = wbase * robust
    return dict(zip(names, fit.x)), float(sigma)


def evaluate(d, coefficients, terms):
    out = d.copy()
    beta_s = max(0.0, coefficients.get("beta_snow", 0.0))
    beta_a = max(0.0, coefficients.get("beta_albedo", 0.0))
    snow_term = beta_s * out.snow_increase_fraction.to_numpy(float)
    alb_term = beta_a * out.albedo_increase_fraction.to_numpy(float)
    gain = 1.0 + snow_term + alb_term
    if not np.all(np.isfinite(gain) & (gain >= 1)):
        raise ValueError("Gain outside physical sign constraint")
    r_obs = np.maximum(
        0.0, out.baseline_radiance.to_numpy(float) +
        out.delta_radiance.to_numpy(float)
    )
    corrected = r_obs / gain
    if np.any(corrected > r_obs + 1e-8):
        raise AssertionError("Surface correction increased radiance")
    out["gain_snow_component"] = snow_term
    out["gain_albedo_component"] = alb_term
    out["surface_gain"] = gain
    out["radiance_observed_2021"] = r_obs
    out["radiance_surface_counterfactual"] = corrected
    out["surface_radiance_removed"] = r_obs - corrected
    out["residual_to_flanks_before"] = out.delta_radiance
    out["residual_to_flanks_after"] = (
        corrected - out.baseline_radiance.to_numpy(float))
    out["gain_gt_2"] = gain > 2
    out["snow_albedo_no_increase"] = (
        (out.snow_increase_fraction == 0) &
        (out.albedo_increase_fraction == 0))
    return out


def region_summary(holdouts):
    rows=[]
    for (model,month,nuts), g in holdouts.groupby(["model","month","NUTS_ID"]):
        w = g.area_weight.to_numpy(float)
        rows.append({
            "model":model, "month":int(month), "NUTS_ID":nuts,
            "n_snow_cell_brightness_groups":len(g),
            "area_sum":sum(w),
            "radiance_before":np.average(g.radiance_observed_2021,weights=w),
            "radiance_after":np.average(g.radiance_surface_counterfactual,weights=w),
            "reference_flanks":np.average(g.baseline_radiance,weights=w),
            "mean_removed":np.average(g.surface_radiance_removed,weights=w),
            "mean_gain":np.average(g.surface_gain,weights=w),
            "mean_residual_before":np.average(g.residual_to_flanks_before,weights=w),
            "mean_residual_after":np.average(g.residual_to_flanks_after,weights=w),
        })
    return pd.DataFrame(rows)


def cohort(x, min_albedo):
    required=["NUTS_ID","month","area_weight","snow_delta_pp",
              "albedo_delta_pp","baseline_radiance","delta_cloudfree",
              "delta_radiance","albedo_valid_pixel_fraction"]
    missing=set(required)-set(x.columns)
    if missing:
        raise ValueError(f"Missing columns {sorted(missing)}")
    d=x.copy()
    for n in required:
        if n!="NUTS_ID":
            d[n]=pd.to_numeric(d[n],errors="coerce")
    d=d[d.month.isin([1,2,3])]
    d=d[np.isfinite(d[required[1:]]).all(axis=1) &
        (d.area_weight>0)&(d.baseline_radiance>=0)&
        (d.albedo_valid_pixel_fraction>=min_albedo)].copy()
    d["snow_increase_fraction"] = np.maximum(d.snow_delta_pp,0)/100.0
    d["albedo_increase_fraction"] = np.maximum(d.albedo_delta_pp,0)/100.0
    d["NUTS_ID"]=d.NUTS_ID.astype(str)
    return d


def diagnose(d):
    records=[]
    for (model,month), v in d.groupby(["model","month"]):
        w=weights_for_region_month(v)
        before=v.residual_to_flanks_before.to_numpy(float)
        after=v.residual_to_flanks_after.to_numpy(float)
        records.append({
            "model":model,"month":int(month),
            "n_regions":v.NUTS_ID.nunique(),"n_groups":len(v),
            "rmse_before":float(np.sqrt(np.average(before**2,weights=w))),
            "rmse_after":float(np.sqrt(np.average(after**2,weights=w))),
            "regional_improved_fraction":float(np.mean([
                abs(g.residual_to_flanks_after.mean()) <
                abs(g.residual_to_flanks_before.mean())
                for _,g in v.groupby("NUTS_ID")])),
            "fraction_groups_gain_above_2":float(np.mean(v.gain_gt_2)),
            "max_gain":float(v.surface_gain.max()),
            "median_removed_radiance":float(v.surface_radiance_removed.median()),
            "fraction_groups_surface_free":float(v.snow_albedo_no_increase.mean()),
        })
    return pd.DataFrame(records)


def unit_tests():
    rng=np.random.default_rng(7)
    n=200
    d=pd.DataFrame({
        "NUTS_ID":np.repeat(list("ABCDEFGHIJ"),20),
        "month":np.tile(np.repeat([1,2],10),10),
        "area_weight":np.ones(n),
        "baseline_radiance":np.repeat([1,3,7,14,40],40).astype(float),
        "snow_delta_pp":rng.uniform(-10,60,n),
        "albedo_delta_pp":rng.uniform(-5,25,n),
        "delta_cloudfree":rng.normal(size=n),
        "albedo_valid_pixel_fraction":np.ones(n),
    })
    d["delta_radiance"]=(.2 + .03*d.delta_cloudfree +
                          d.baseline_radiance*(.30*np.maximum(d.snow_delta_pp,0)/100.0
                          +1.1*np.maximum(d.albedo_delta_pp,0)/100.0))
    d=cohort(d,.5)
    fit,_=constrained_fit(d[d.month.isin([1,2])],MODELS["snow_plus_albedo"])
    assert abs(fit["beta_snow"]-.30)<1e-3 and abs(fit["beta_albedo"]-1.1)<1e-3,fit
    correct=evaluate(d,fit,MODELS["snow_plus_albedo"])
    assert (correct.radiance_surface_counterfactual <=
            correct.radiance_observed_2021 +1e-10).all()
    assert (correct.loc[correct.snow_albedo_no_increase,
                        "surface_radiance_removed"].abs()<1e-10).all()
    print("Constrained gain calibration synthetic test: PASS")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--groups",default=
        "data/processed/historical_prc_2021_pixel_snow_test/"
        "matched_snowcell_brightness_groups.csv")
    ap.add_argument("--output-dir",default=
        "data/processed/historical_prc_2021_surface_gain_calibration")
    ap.add_argument("--min-albedo-valid-fraction",type=float,default=.5)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    unit_tests()
    if args.self_test:return
    p=Path(args.groups)
    if not p.is_absolute():p=PROJECT_ROOT/p
    if not p.exists():raise FileNotFoundError(p)
    d=cohort(pd.read_csv(p,dtype={"snow_cell_id":str}),
             args.min_albedo_valid_fraction)
    if d.NUTS_ID.nunique()<8:raise ValueError("Too few regions for spatial holdout")
    print(f"Eligible {len(d)} grouped observations, {d.NUTS_ID.nunique()} regions",
          flush=True)
    records,holdouts=[],[]
    for label,terms in MODELS.items():
        train=d[d.month.isin([1,2])]
        full,robust_sigma=constrained_fit(train,terms)
        for term,value in full.items():
            records.append({
                "model":label,"fit_scope":"all_Jan_Feb",
                "term":term,"coefficient":value,
                "robust_residual_scale_nw":robust_sigma,
                "n_training_regions":train.NUTS_ID.nunique(),
                "n_training_groups":len(train),
            })
        for idx,held_id in enumerate(sorted(d.NUTS_ID.unique()),1):
            training=d[(d.NUTS_ID!=held_id)&d.month.isin([1,2])]
            testing=d[d.NUTS_ID==held_id]
            params,_=constrained_fit(training,terms)
            h=evaluate(testing,params,terms)
            h["model"]=label
            h["held_out_region"]=held_id
            holdouts.append(h)
        print(f"Model {label}: beta snow={full.get('beta_snow',0):.5f}, "
              f"beta albedo={full.get('beta_albedo',0):.5f}",flush=True)
    held=pd.concat(holdouts,ignore_index=True)
    results=diagnose(held)
    by_region=region_summary(held)
    out=Path(args.output_dir)
    if not out.is_absolute():out=PROJECT_ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(records).to_csv(out/"surface_gain_calibration_coefficients.csv",index=False)
    held.to_csv(out/"surface_gain_heldout_groups.csv",index=False)
    results.to_csv(out/"surface_gain_holdout_diagnostics.csv",index=False)
    by_region.to_csv(out/"surface_gain_region_month_diagnostics.csv",index=False)
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.8),layout="constrained")
    s=np.linspace(0,1,101)
    ca=pd.DataFrame(records)
    for label,terms in MODELS.items():
        f=ca[ca.model==label].set_index("term").coefficient.to_dict()
        axes[0].plot(s,1+f.get("beta_snow",0)*s,label=label)
        axes[1].plot(s,1+f.get("beta_albedo",0)*s,label=label)
    axes[0].set(xlabel="Additional snow cover fraction",ylabel="Predicted amplification G",
                title="Snow gain at unchanged albedo")
    axes[1].set(xlabel="Additional visible albedo fraction",ylabel="Predicted amplification G",
                title="Albedo gain at unchanged snow")
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend()
    fig.savefig(out/"surface_gain_fitted_curves.pdf")
    plt.close(fig)
    statements=[
        "2021 VIIRS surface-state empirical amplification CALIBRATION TEST",
        "Correction uses R_corrected=R_observed / G, G>=1. Never increases R.",
        "G=1+beta_s*max(snow anomaly fraction,0)+beta_a*max(albedo anomaly fraction,0)",
        "Beta values fitted NONNEGATIVE with independent nuisance intercept,",
        "February term, cloud-free change and brightness-trend term.",
        "All models fitted to matched 2020/21/22 same-month data, Jan-Feb only.",
        "Validation: leave entire NUTS3 region out; March is temporal placebo.",
        "Each region/month contributes equal model weight; within region uses area.",
        "All 3 model families use the same QA-complete albedo subset.",
        "This is a conditional observational curve, not independent",
        "laboratory validation or a confirmed causal snow-reflectance law.",
        "Snow and albedo can be collinear: evaluate albedo_only vs snow_only",
        "against combined. Do not adopt a fit because 2021 flanks look smoother.",
        "The established annual VNL -0.235 correction is unchanged. This",
        "stage contains only monthly EOG measurements and no annual correction.",
        f"Eligible {len(d)} cell/brightness groups in {d.NUTS_ID.nunique()} regions.",
        "",
        "Full Jan-Feb fitted coefficients:",
        pd.DataFrame(records).to_string(index=False),
        "",
        "Out-of-NUTS3 holdout diagnostics:",
        results.to_string(index=False),
    ]
    (out/"surface_gain_summary.txt").write_text("\n".join(statements)+"\n")
    print("\n".join(statements))
    print("Output folder:",out)


if __name__=="__main__":
    main()
