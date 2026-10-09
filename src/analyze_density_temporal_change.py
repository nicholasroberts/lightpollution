#!/usr/bin/env python3
"""Test whether temporal VIIRS change depends on population density.

Hypothesis: countries differ partly because they sample different portions of
a common rural-to-urban relationship. If low-density NUTS3 regions have gained
radiance faster than dense regions, the European log-log relationship should
rotate through time even without intrinsically different national light demand.
"""

from __future__ import annotations

import argparse
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
import statsmodels.formula.api as smf

from project_config import configured_path, load_config


def fit_mixedlm(formula: str, data: pd.DataFrame):
    model = smf.mixedlm(formula, data=data, groups=data["NUTS_ID"], re_formula="1")
    last = None
    for method in ("lbfgs", "bfgs", "cg"):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = model.fit(reml=False, method=method, maxiter=3000, disp=False)
            if result.converged:
                return result, method
            last = RuntimeError(f"non-converged with {method}")
        except Exception as exc:
            last = exc
    raise RuntimeError(f"MixedLM failed for {formula}: {last}")


def lrt(reduced, full):
    lr = max(0.0, 2.0 * (full.llf - reduced.llf))
    df = int(len(full.fe_params) - len(reduced.fe_params))
    return lr, df, float(chi2.sf(lr, df)) if df > 0 else np.nan


def contrast(result, weights):
    names = list(result.fe_params.index)
    L = np.zeros(len(names))
    for name, value in weights.items():
        L[names.index(name)] = value
    beta = result.fe_params.to_numpy(float)
    cov = result.cov_params().loc[names, names].to_numpy(float)
    est = float(L @ beta)
    se = float(np.sqrt(max(0.0, L @ cov @ L)))
    z = est / se if se > 0 else np.nan
    p = float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
    return est, se, p


def prepare_data(config):
    cfg = config["analysis"]["density_temporal_change"]
    path = configured_path(config, "processed") / config["outputs"]["merged"]
    d = pd.read_parquet(path).copy()
    d = d[d["year"].between(cfg["first_year"], cfg["last_year"])].copy()

    ok = (
        np.isfinite(pd.to_numeric(d["radiance_mean_corrected"], errors="coerce"))
        & np.isfinite(pd.to_numeric(d["population_density"], errors="coerce"))
        & (d["radiance_mean_corrected"] > 0)
        & (d["population_density"] > 0)
        & d["NUTS_ID"].notna()
        & d["CNTR_CODE"].notna()
    )
    d = d.loc[ok].copy()

    repeated = d.groupby("NUTS_ID").size()
    d = d[d["NUTS_ID"].isin(repeated[repeated >= cfg["min_repeated_years"]].index)].copy()
    counts = d.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    d = d[d["CNTR_CODE"].isin(counts[counts >= cfg["min_country_regions"]].index)].copy()

    d["log_rad"] = np.log10(d["radiance_mean_corrected"].astype(float))
    d["log_pop"] = np.log10(d["population_density"].astype(float))
    midpoint = (cfg["first_year"] + cfg["last_year"]) / 2
    d["year_c"] = d["year"] - midpoint
    d["country"] = d["CNTR_CODE"].astype("category")

    ref = (
        d.groupby("NUTS_ID", as_index=False)["population_density"]
        .median()
        .rename(columns={"population_density": "reference_population_density"})
    )
    d = d.merge(ref, on="NUTS_ID", how="left", validate="many_to_one")
    edges = [-np.inf] + list(cfg["density_band_breaks"]) + [np.inf]
    d["density_band"] = pd.cut(
        d["reference_population_density"],
        bins=edges,
        labels=cfg["density_band_labels"],
        right=False,
        ordered=True,
    )

    center = float(d["log_pop"].median())
    d["log_pop_c"] = d["log_pop"] - center
    d.attrs["log_pop_center"] = center
    d.attrs["year_midpoint"] = midpoint
    return d


def annual_band_summary(d):
    rows = []
    for (year, band), g in d.groupby(["year", "density_band"], observed=True):
        rows.append({
            "year": int(year),
            "density_band": str(band),
            "n_nuts3": int(g["NUTS_ID"].nunique()),
            "n_countries": int(g["CNTR_CODE"].nunique()),
            "median_log10_radiance": float(g["log_rad"].median()),
            "mean_log10_radiance": float(g["log_rad"].mean()),
            "median_radiance": float(g["radiance_mean_corrected"].median()),
            "median_population_density": float(g["population_density"].median()),
        })
    return pd.DataFrame(rows)


def country_band_composition(d):
    r = d[["NUTS_ID","CNTR_CODE","density_band","reference_population_density"]].drop_duplicates("NUTS_ID")
    out = (
        r.groupby(["CNTR_CODE","density_band"], observed=True)
        .size().rename("n_nuts3").reset_index()
    )
    totals = r.groupby("CNTR_CODE").size().rename("country_nuts3_total").reset_index()
    out = out.merge(totals, on="CNTR_CODE", how="left")
    out["share_of_country_nuts3"] = out["n_nuts3"] / out["country_nuts3_total"]
    return out.sort_values(["CNTR_CODE","density_band"])


def fit_band_model(d, labels):
    reduced, opt0 = fit_mixedlm("log_rad ~ year_c + C(density_band) + C(country)", d)
    full, opt1 = fit_mixedlm("log_rad ~ year_c * C(density_band) + C(country)", d)
    lr, df, p = lrt(reduced, full)

    rows = []
    ref = labels[0]
    for band in labels:
        w = {"year_c": 1.0}
        if band != ref:
            candidates = [
                f"year_c:C(density_band)[T.{band}]",
                f"C(density_band)[T.{band}]:year_c",
            ]
            name = next(x for x in candidates if x in full.fe_params.index)
            w[name] = 1.0
        est, se, pp = contrast(full, w)
        rows.append({
            "density_band": band,
            "log10_radiance_trend_per_year": est,
            "trend_se": se,
            "trend_p": pp,
            "pct_radiance_change_per_year": 100 * (10**est - 1),
        })
    info = {
        "lr": lr, "df": df, "p": p, "reduced_optimizer": opt0,
        "full_optimizer": opt1, "n_obs": len(d),
        "n_nuts3": d["NUTS_ID"].nunique(),
        "n_countries": d["CNTR_CODE"].nunique(),
    }
    return pd.DataFrame(rows), info


def fit_continuous_model(d, config):
    base = "log_rad ~ C(country) + year_c + log_pop_c + I(log_pop_c ** 2) + I(log_pop_c ** 3)"
    common = "log_rad ~ C(country) + year_c * (log_pop_c + I(log_pop_c ** 2) + I(log_pop_c ** 3))"
    cy = common + " + year_c:C(country)"
    cd = common + " + log_pop_c:C(country)"

    m0, o0 = fit_mixedlm(base, d)
    m1, o1 = fit_mixedlm(common, d)
    m2, o2 = fit_mixedlm(cy, d)
    m3, o3 = fit_mixedlm(cd, d)

    tests = []
    for name, a, b, oa, ob in [
        ("density-dependent temporal trend", m0, m1, o0, o1),
        ("residual country-specific temporal trends", m1, m2, o1, o2),
        ("residual country-specific density slopes", m1, m3, o1, o3),
    ]:
        lr, df, p = lrt(a, b)
        tests.append({"test": name, "lr": lr, "df": df, "p": p,
                      "reduced_optimizer": oa, "full_optimizer": ob})

    cfg = config["analysis"]["density_temporal_change"]
    qlo, qhi = cfg["continuous_density_quantiles"]
    lo = float(d["population_density"].quantile(qlo))
    hi = float(d["population_density"].quantile(qhi))
    grid = np.geomspace(lo, hi, cfg["continuous_grid_points"])
    c = d.attrs["log_pop_center"]

    rows = []
    for density in grid:
        x = np.log10(density) - c
        w = {
            "year_c": 1.0,
            "year_c:log_pop_c": x,
            "year_c:I(log_pop_c ** 2)": x**2,
            "year_c:I(log_pop_c ** 3)": x**3,
        }
        est, se, p = contrast(m1, w)
        rows.append({
            "population_density": density,
            "log10_radiance_trend_per_year": est,
            "trend_se": se,
            "trend_ci_low": est - 1.96*se,
            "trend_ci_high": est + 1.96*se,
            "trend_p": p,
            "pct_radiance_change_per_year": 100*(10**est - 1),
        })
    return pd.DataFrame(rows), pd.DataFrame(tests)


def make_plots(annual, band_trends, continuous, labels, outdir):
    outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.5,6))
    for band in labels:
        g = annual[annual["density_band"] == band].sort_values("year")
        ax.plot(g["year"], g["median_log10_radiance"], marker="o", label=band)
    ax.set_xlabel("Year")
    ax.set_ylabel("Median log10 corrected VIIRS radiance")
    ax.set_title("Radiance through time by fixed NUTS3 population-density band")
    ax.grid(alpha=.2)
    ax.legend(title="People km⁻²")
    fig.tight_layout()
    p = outdir / "density_band_radiance_by_year.pdf"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5,5.6))
    x = np.arange(len(band_trends))
    est = band_trends["pct_radiance_change_per_year"].to_numpy()
    se_log = band_trends["trend_se"].to_numpy()
    lo = 100*(10**(band_trends["log10_radiance_trend_per_year"].to_numpy()-1.96*se_log)-1)
    hi = 100*(10**(band_trends["log10_radiance_trend_per_year"].to_numpy()+1.96*se_log)-1)
    ax.errorbar(x, est, yerr=[est-lo, hi-est], fmt="o", capsize=4)
    ax.axhline(0, linewidth=.9)
    ax.set_xticks(x); ax.set_xticklabels(band_trends["density_band"])
    ax.set_xlabel("Fixed reference population-density band (people km⁻²)")
    ax.set_ylabel("Estimated radiance change per year (%)")
    ax.set_title("Country-adjusted temporal trend by population-density band")
    ax.grid(alpha=.2, axis="y")
    fig.tight_layout()
    p = outdir / "density_band_temporal_slopes.pdf"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5,5.6))
    ax.plot(continuous["population_density"], continuous["pct_radiance_change_per_year"])
    low = 100*(10**continuous["trend_ci_low"]-1)
    high = 100*(10**continuous["trend_ci_high"]-1)
    ax.fill_between(continuous["population_density"], low, high, alpha=.2)
    ax.axhline(0, linewidth=.9)
    ax.set_xscale("log")
    ax.set_xlabel("Population density (people km⁻²)")
    ax.set_ylabel("Estimated corrected-VIIRS radiance change per year (%)")
    ax.set_title("Continuous density dependence of European temporal radiance trend")
    ax.grid(alpha=.2)
    fig.tight_layout()
    p = outdir / "temporal_trend_vs_population_density.pdf"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()
    config = load_config(args.config)
    cfg = config["analysis"]["density_temporal_change"]
    d = prepare_data(config)

    annual = annual_band_summary(d)
    composition = country_band_composition(d)
    bands, band_info = fit_band_model(d, cfg["density_band_labels"])
    continuous, tests = fit_continuous_model(d, config)

    processed = configured_path(config, "processed")
    outdir = configured_path(config, "figures") / "population_relationship" / "density_temporal_change"

    outputs = {
        "annual": processed/"density_band_year_summary.csv",
        "bands": processed/"density_band_temporal_trends.csv",
        "composition": processed/"country_density_band_composition.csv",
        "continuous": processed/"continuous_density_temporal_trend.csv",
        "tests": processed/"density_temporal_model_tests.csv",
        "summary": processed/"density_temporal_change_summary.txt",
    }
    annual.to_csv(outputs["annual"], index=False)
    bands.to_csv(outputs["bands"], index=False)
    composition.to_csv(outputs["composition"], index=False)
    continuous.to_csv(outputs["continuous"], index=False)
    tests.to_csv(outputs["tests"], index=False)
    make_plots(annual, bands, continuous, cfg["density_band_labels"], outdir)

    with outputs["summary"].open("w") as h:
        h.write("POPULATION-DENSITY DEPENDENCE OF TEMPORAL VIIRS CHANGE\n")
        h.write("="*72 + "\n\n")
        h.write(f"Interval: {cfg['first_year']}-{cfg['last_year']}\n")
        h.write(f"Observations: {band_info['n_obs']}; NUTS3: {band_info['n_nuts3']}; countries: {band_info['n_countries']}\n\n")
        h.write("DENSITY-BAND TEMPORAL TRENDS\n")
        h.write(bands.to_string(index=False))
        h.write(f"\n\nBand x year LRT: LR={band_info['lr']:.6g}, df={band_info['df']}, p={band_info['p']:.6g}\n\n")
        h.write("CONTINUOUS MODEL BLOCK TESTS\n")
        h.write(tests.to_string(index=False))
        h.write("\n")

    print("\nDENSITY-BAND TEMPORAL TRENDS")
    print(bands.to_string(index=False, float_format=lambda x:f"{x:.6g}"))
    print(f"\nBand x year LRT: LR={band_info['lr']:.6g}, df={band_info['df']}, p={band_info['p']:.6g}")
    print("\nCONTINUOUS MODEL BLOCK TESTS")
    print(tests[["test","lr","df","p","full_optimizer"]].to_string(index=False, float_format=lambda x:f"{x:.6g}"))
    print("\nWROTE")
    for p in outputs.values():
        print(p)


if __name__ == "__main__":
    main()
