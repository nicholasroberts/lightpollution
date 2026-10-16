#!/usr/bin/env python3
"""Test whether anomalous winter snow amplifies the Jan/Feb 2021 VIIRS pulse.

Primary response:
    2021 monthly radiance - mean(same month 2020, 2022)

Snow exposure:
    MOD10CM 2021 snow cover - mean(same month 2020, 2022)

Primary model:
    absolute anomaly ~ log1p(flank radiance) + density band + country
                       + snow anomaly
                       + snow anomaly x log1p(flank radiance)

VCMSLCFG February is primary, January supportive, March a negative control.
VCMCFG repeats the same models as a product-robustness check. This diagnostic
does not alter the annual production time series.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr
import statsmodels.formula.api as smf

from diagnose_fixed_nuts2016_2021_reversal import add_latitude
from project_config import PROJECT_ROOT, configured_path, load_config
from test_fixed_nuts2016_snow_viirs_2021 import (
    ensure_snow_geotiff,
    extract_snow_month,
    fixed_nuts,
    month_from_year_doy,
    normalize_snow_extract,
)

YEARS = (2020, 2021, 2022)
DENSITY_ORDER = ("<30", "30-100", "100-300", "300-1000", ">=1000")
MONTH_NAME = {1: "Jan", 2: "Feb", 3: "Mar"}
HDF_RE = re.compile(r"^MOD10CM\.A(?P<year>\d{4})(?P<doy>\d{3})\.061\..*\.hdf$")


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def spearman_clean(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    d = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 3 or d["x"].nunique() < 2 or d["y"].nunique() < 2:
        return np.nan, np.nan, len(d)
    r = spearmanr(d["x"], d["y"])
    return float(r.statistic), float(r.pvalue), int(len(d))


def discover_snow(raw_dir: Path, months: tuple[int, ...]) -> dict[tuple[int, int], Path]:
    found: dict[tuple[int, int], Path] = {}
    for p in sorted(raw_dir.glob("MOD10CM*.hdf")):
        m = HDF_RE.match(p.name)
        if not m:
            continue
        year = int(m.group("year"))
        if year not in YEARS:
            continue
        month = month_from_year_doy(year, int(m.group("doy")))
        if month not in months:
            continue
        key = (year, month)
        if key in found:
            raise RuntimeError(f"Duplicate MOD10CM month {key}: {found[key]} and {p}")
        found[key] = p

    missing = sorted(
        {(y, m) for y in YEARS for m in months} - set(found)
    )
    if missing:
        raise FileNotFoundError(
            "Missing MOD10CM months: "
            + ", ".join(f"{y}-{m:02d}" for y, m in missing)
            + ". Run 'python src/run_pipeline.py future-radiance-download-snow-history'."
        )
    return found


def load_atlas(cfg: dict) -> pd.DataFrame:
    path = PROJECT_ROOT / str(cfg["monthly_region_file"])
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run "
            "'python src/run_pipeline.py future-radiance-monthly-product-atlas'."
        )
    d = pd.read_csv(path)
    required = {
        "product", "month", "NUTS_ID", "CNTR_CODE", "density_band",
        "flank_mean", "delta_nw", "pct_change", "paired_valid",
    }
    missing = required - set(d.columns)
    if missing:
        raise RuntimeError(f"Monthly atlas file missing {sorted(missing)}")
    d["NUTS_ID"] = d["NUTS_ID"].astype(str)
    d["CNTR_CODE"] = d["CNTR_CODE"].astype(str)
    d["product"] = d["product"].astype(str).str.lower()
    d["month"] = pd.to_numeric(d["month"], errors="raise").astype(int)
    d["paired_valid"] = as_bool(d["paired_valid"])
    d["density_band"] = pd.Categorical(
        d["density_band"].astype(str),
        categories=list(DENSITY_ORDER),
        ordered=True,
    )
    return d


def build_snow_anomaly(
    config: dict,
    cfg: dict,
    atlas_ids: set[str],
) -> pd.DataFrame:
    raw_dir = PROJECT_ROOT / str(cfg["mod10cm_raw_dir"])
    if not raw_dir.exists():
        raise FileNotFoundError(raw_dir)

    months = tuple(int(x) for x in cfg["months"])
    files = discover_snow(raw_dir, months)
    tif_cache = configured_path(config, "processed") / "mod10cm_061_geotiff_cache"
    zonal_cache = configured_path(config, "processed") / "mod10cm_061_nuts2016_cache"

    nuts = fixed_nuts(config)
    id_field = config["nuts"]["id_field"]
    nuts[id_field] = nuts[id_field].astype(str)

    frames = []
    for (year, month), hdf in sorted(files.items()):
        tif = ensure_snow_geotiff(hdf, tif_cache, config["nuts"]["crs"])
        print(f"Snow {year}-{month:02d}: {hdf.name}", flush=True)
        frames.append(
            extract_snow_month(tif, year, month, nuts, config, zonal_cache)
        )

    snow = normalize_snow_extract(pd.concat(frames, ignore_index=True))
    snow["NUTS_ID"] = snow[id_field].astype(str)
    snow = snow[snow["NUTS_ID"].isin(atlas_ids)].copy()

    rows = []
    for month in months:
        g = snow[snow["month"] == month]
        v = g.pivot(index="NUTS_ID", columns="year", values="snow_mean_pct")
        q = g.pivot(
            index="NUTS_ID",
            columns="year",
            values="snow_valid_area_fraction",
        )
        if not set(YEARS).issubset(v.columns) or not set(YEARS).issubset(q.columns):
            raise RuntimeError(f"Incomplete snow pivot for month {month}.")
        out = v.index.to_frame(index=False, name="NUTS_ID")
        for year in YEARS:
            out[f"snow_{year}_pct"] = v[year].to_numpy()
            out[f"snow_valid_{year}"] = q[year].to_numpy()
        out["snow_flank_mean_pct"] = (
            out["snow_2020_pct"] + out["snow_2022_pct"]
        ) / 2
        out["snow_anomaly_pp"] = (
            out["snow_2021_pct"] - out["snow_flank_mean_pct"]
        )
        out["min_snow_valid_fraction"] = out[
            [f"snow_valid_{y}" for y in YEARS]
        ].min(axis=1)
        out["month"] = month
        out["month_name"] = MONTH_NAME[month]
        rows.append(out)
    return pd.concat(rows, ignore_index=True)


def prepare(
    atlas: pd.DataFrame,
    snow: pd.DataFrame,
    config: dict,
    cfg: dict,
) -> pd.DataFrame:
    months = tuple(int(x) for x in cfg["months"])
    products = tuple(str(x).lower() for x in cfg["products"])
    d = atlas[
        atlas["month"].isin(months) & atlas["product"].isin(products)
    ].copy()
    d = d.merge(
        snow,
        on=["NUTS_ID", "month"],
        how="left",
        validate="many_to_one",
    )
    d = add_latitude(d, config)

    for col in [
        "flank_mean", "delta_nw", "pct_change",
        "snow_anomaly_pp", "min_snow_valid_fraction",
    ]:
        d[col] = pd.to_numeric(d[col], errors="coerce")

    d["baseline_nonnegative_nw"] = d["flank_mean"].clip(lower=0)
    d["baseline_log1p"] = np.log1p(d["baseline_nonnegative_nw"])
    d["eligible"] = (
        d["paired_valid"]
        & d["delta_nw"].notna()
        & d["baseline_log1p"].notna()
        & d["snow_anomaly_pp"].notna()
        & (d["min_snow_valid_fraction"] >= float(cfg["min_snow_valid_area_fraction"]))
        & d["density_band"].notna()
        & d["CNTR_CODE"].notna()
        & d["centroid_latitude"].notna()
    )
    return d


def partial_r2(reduced, full) -> float:
    if reduced.ssr <= 0 or not np.isfinite(reduced.ssr):
        return np.nan
    return float(max(0, (reduced.ssr - full.ssr) / reduced.ssr))


def fit_models(
    d: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    model_rows = []
    effect_rows = []
    coef_rows = []

    for (product, month), g0 in d.groupby(["product", "month"], observed=True):
        g = g0[g0["eligible"]].copy()
        if len(g) < 100:
            raise RuntimeError(f"Too few rows for {product} {month}: {len(g)}")

        mu = float(g["baseline_log1p"].mean())
        sd = float(g["baseline_log1p"].std(ddof=0))
        raw_mu = float(g["baseline_nonnegative_nw"].mean())
        raw_sd = float(g["baseline_nonnegative_nw"].std(ddof=0))
        if sd <= 0 or raw_sd <= 0:
            raise RuntimeError(f"Degenerate baseline for {product} {month}.")
        g["baseline_z"] = (g["baseline_log1p"] - mu) / sd
        g["baseline_raw_z"] = (
            g["baseline_nonnegative_nw"] - raw_mu
        ) / raw_sd

        f0 = "delta_nw ~ baseline_z + C(density_band) + C(CNTR_CODE)"
        f1 = f0 + " + snow_anomaly_pp"
        f2 = (
            "delta_nw ~ snow_anomaly_pp * baseline_z "
            "+ C(density_band) + C(CNTR_CODE)"
        )
        f_raw = (
            "delta_nw ~ snow_anomaly_pp * baseline_raw_z "
            "+ C(density_band) + C(CNTR_CODE)"
        )
        f_lat = (
            "delta_nw ~ snow_anomaly_pp * baseline_z "
            "+ centroid_latitude + C(density_band) + C(CNTR_CODE)"
        )

        m0 = smf.ols(f0, data=g).fit(cov_type="HC3")
        m1 = smf.ols(f1, data=g).fit(cov_type="HC3")
        m2 = smf.ols(f2, data=g).fit(cov_type="HC3")
        m_raw = smf.ols(f_raw, data=g).fit(cov_type="HC3")
        m_lat = smf.ols(f_lat, data=g).fit(cov_type="HC3")

        snow_term = "snow_anomaly_pp"
        int_term = "snow_anomaly_pp:baseline_z"
        raw_int = "snow_anomaly_pp:baseline_raw_z"
        joint = m2.wald_test(
            f"{snow_term} = 0, {int_term} = 0",
            scalar=True,
        )

        role = (
            "primary" if (product == "vcmslcfg" and month == 2)
            else "supportive" if (product == "vcmslcfg" and month == 1)
            else "negative_control" if (product == "vcmslcfg" and month == 3)
            else "product_robustness"
        )
        rb, pb, _ = spearman_clean(g["flank_mean"], g["delta_nw"])
        rs, ps, _ = spearman_clean(g["snow_anomaly_pp"], g["delta_nw"])

        model_rows.append({
            "product": product,
            "month": int(month),
            "month_name": MONTH_NAME[int(month)],
            "role": role,
            "n": int(m2.nobs),
            "n_countries": int(g["CNTR_CODE"].nunique()),
            "r2_base": float(m0.rsquared),
            "r2_plus_snow": float(m1.rsquared),
            "r2_full": float(m2.rsquared),
            "partial_r2_snow_main": partial_r2(m0, m1),
            "partial_r2_interaction_after_snow": partial_r2(m1, m2),
            "partial_r2_full_snow_block": partial_r2(m0, m2),
            "snow_coef_at_mean_brightness": float(m2.params[snow_term]),
            "snow_p": float(m2.pvalues[snow_term]),
            "snow_x_logbrightness_coef": float(m2.params[int_term]),
            "snow_x_logbrightness_p": float(m2.pvalues[int_term]),
            "joint_snow_block_wald": float(np.asarray(joint.statistic).squeeze()),
            "joint_snow_block_p": float(np.asarray(joint.pvalue).squeeze()),
            "raw_baseline_interaction_coef": float(m_raw.params[raw_int]),
            "raw_baseline_interaction_p": float(m_raw.pvalues[raw_int]),
            "latitude_sensitivity_snow_coef": float(m_lat.params[snow_term]),
            "latitude_sensitivity_snow_p": float(m_lat.pvalues[snow_term]),
            "latitude_sensitivity_interaction_coef": float(m_lat.params[int_term]),
            "latitude_sensitivity_interaction_p": float(m_lat.pvalues[int_term]),
            "spearman_baseline_vs_delta": rb,
            "spearman_baseline_vs_delta_p": pb,
            "spearman_snow_vs_delta": rs,
            "spearman_snow_vs_delta_p": ps,
        })

        ci = m2.conf_int()
        for term in m2.params.index:
            coef_rows.append({
                "product": product,
                "month": int(month),
                "month_name": MONTH_NAME[int(month)],
                "role": role,
                "term": str(term),
                "coef": float(m2.params[term]),
                "std_error_hc3": float(m2.bse[term]),
                "ci_low": float(ci.loc[term, 0]),
                "ci_high": float(ci.loc[term, 1]),
                "p_value": float(m2.pvalues[term]),
            })

        cov = m2.cov_params()
        b0 = float(m2.params[snow_term])
        bi = float(m2.params[int_term])
        for quantile in (0.10, 0.25, 0.50, 0.75, 0.90):
            base = float(g["baseline_nonnegative_nw"].quantile(quantile))
            z = (math.log1p(base) - mu) / sd
            eff = b0 + bi * z
            var = (
                float(cov.loc[snow_term, snow_term])
                + z * z * float(cov.loc[int_term, int_term])
                + 2 * z * float(cov.loc[snow_term, int_term])
            )
            se = math.sqrt(max(var, 0))
            zstat = eff / se if se > 0 else np.nan
            p = 2 * norm.sf(abs(zstat)) if np.isfinite(zstat) else np.nan
            effect_rows.append({
                "product": product,
                "month": int(month),
                "month_name": MONTH_NAME[int(month)],
                "baseline_quantile": quantile,
                "baseline_radiance_nw": base,
                "snow_effect_nw_per_pp": eff,
                "ci_low": eff - 1.96 * se,
                "ci_high": eff + 1.96 * se,
                "p_value": p,
            })

    return pd.DataFrame(model_rows), pd.DataFrame(effect_rows), pd.DataFrame(coef_rows)


def summaries(d: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    e = d[d["eligible"]].copy()
    density = (
        e.groupby(["product", "month", "density_band"], observed=True)
        .agg(
            n=("NUTS_ID", "nunique"),
            n_countries=("CNTR_CODE", "nunique"),
            median_baseline_nw=("flank_mean", "median"),
            mean_delta_nw=("delta_nw", "mean"),
            median_delta_nw=("delta_nw", "median"),
            median_pct_change=("pct_change", "median"),
            median_snow_anomaly_pp=("snow_anomaly_pp", "median"),
        )
        .reset_index()
    )
    density.insert(2, "month_name", density["month"].map(MONTH_NAME))

    country = (
        e.groupby(["product", "month", "CNTR_CODE"], observed=True)
        .agg(
            n=("NUTS_ID", "nunique"),
            median_baseline_nw=("flank_mean", "median"),
            median_delta_nw=("delta_nw", "median"),
            median_pct_change=("pct_change", "median"),
            median_snow_anomaly_pp=("snow_anomaly_pp", "median"),
        )
        .reset_index()
    )
    country.insert(2, "month_name", country["month"].map(MONTH_NAME))

    frames = []
    for (product, month), g in e.groupby(["product", "month"], observed=True):
        x = g.copy()
        x["baseline_decile"] = (
            pd.qcut(
                x["baseline_nonnegative_nw"],
                q=10,
                labels=False,
                duplicates="drop",
            )
            + 1
        )
        z = (
            x.groupby("baseline_decile", observed=True)
            .agg(
                n=("NUTS_ID", "nunique"),
                median_baseline_nw=("baseline_nonnegative_nw", "median"),
                mean_delta_nw=("delta_nw", "mean"),
                median_delta_nw=("delta_nw", "median"),
                median_snow_anomaly_pp=("snow_anomaly_pp", "median"),
            )
            .reset_index()
        )
        z["product"] = product
        z["month"] = int(month)
        z["month_name"] = MONTH_NAME[int(month)]
        frames.append(z)
    return density, country, pd.concat(frames, ignore_index=True)


def make_plots(
    d: pd.DataFrame,
    models: pd.DataFrame,
    effects: pd.DataFrame,
    deciles: pd.DataFrame,
    out_dir: Path,
    dpi: int,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    q = deciles[deciles["product"] == "vcmslcfg"]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    for month, g in q.groupby("month", observed=True):
        ax.plot(
            g["baseline_decile"],
            g["mean_delta_nw"],
            marker="o",
            label=MONTH_NAME[int(month)],
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Baseline-radiance decile")
    ax.set_ylabel("Mean absolute anomaly (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("VCMSLCFG 2021 anomaly versus normal regional brightness")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_dir / "baseline_decile_absolute_anomaly.pdf", dpi=dpi)
    plt.close(fig)

    feb = d[
        (d["product"] == "vcmslcfg")
        & (d["month"] == 2)
        & d["eligible"]
    ]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    sc = ax.scatter(
        feb["snow_anomaly_pp"],
        feb["delta_nw"],
        c=feb["baseline_log1p"],
        s=14,
        alpha=0.55,
    )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel("February snow-cover anomaly (percentage points)")
    ax.set_ylabel("February absolute radiance anomaly (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Primary snow-amplification test")
    fig.colorbar(sc, ax=ax, label="log1p normal February radiance")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_dir / "february_snow_vs_radiance_anomaly.pdf", dpi=dpi)
    plt.close(fig)

    e = effects[effects["product"] == "vcmslcfg"]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    for month, g in e.groupby("month", observed=True):
        ax.plot(
            g["baseline_radiance_nw"],
            g["snow_effect_nw_per_pp"],
            marker="o",
            label=MONTH_NAME[int(month)],
        )
        ax.fill_between(
            g["baseline_radiance_nw"],
            g["ci_low"],
            g["ci_high"],
            alpha=0.15,
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xscale("symlog", linthresh=0.05)
    ax.set_xlabel("Normal monthly radiance (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_ylabel("Estimated anomaly per +1 snow-cover percentage point")
    ax.set_title("Does snow have a larger effect in brighter regions?")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_dir / "snow_effect_by_baseline_brightness.pdf", dpi=dpi)
    plt.close(fig)

    m = models[models["product"] == "vcmslcfg"].sort_values("month")
    fig, ax = plt.subplots(figsize=(7.5, 5.4))
    ax.bar(m["month_name"], 100 * m["partial_r2_full_snow_block"])
    ax.set_xlabel("Month")
    ax.set_ylabel("Snow block partial R² (%)")
    ax.set_title("Incremental explanatory power of snow + snow×brightness")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(out_dir / "snow_block_partial_r2.pdf", dpi=dpi)
    plt.close(fig)


def write_summary(
    models: pd.DataFrame,
    effects: pd.DataFrame,
    density: pd.DataFrame,
    cfg: dict,
    out: Path,
) -> None:
    feb = models[(models["product"] == "vcmslcfg") & (models["month"] == 2)].iloc[0]
    jan = models[(models["product"] == "vcmslcfg") & (models["month"] == 1)].iloc[0]
    mar = models[(models["product"] == "vcmslcfg") & (models["month"] == 3)].iloc[0]

    with out.open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 MONTHLY SNOW-AMPLIFICATION TEST\n")
        h.write("=" * 78 + "\n\n")
        h.write(
            "Question: does the Jan/Feb 2021 absolute VIIRS pulse track "
            "same-month anomalous snow cover, and is the snow effect larger "
            "where normal artificial-light radiance is higher?\n\n"
        )
        h.write("DESIGN\n")
        h.write("- response: 2021 - mean(2020, 2022), same calendar month\n")
        h.write("- snow: MOD10CM 2021 - mean(2020, 2022), same month\n")
        h.write("- primary: VCMSLCFG February; supportive: January; control: March\n")
        h.write("- VCMCFG repeats the analysis as product robustness\n")
        h.write(
            "- model: log1p baseline radiance + fixed density band + country "
            "+ snow + snow x baseline radiance\n"
        )
        h.write("- HC3 robust covariance; latitude-adjusted sensitivity retained\n")
        h.write(
            f"- minimum MOD10CM valid area: "
            f"{100*float(cfg['min_snow_valid_area_fraction']):.0f}%\n\n"
        )
        h.write("MODEL RESULTS\n")
        h.write(models.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nMARGINAL SNOW EFFECTS\n")
        h.write(effects.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nDENSITY SUMMARY\n")
        h.write(density.to_string(index=False, float_format=lambda x: f"{x:.6g}"))

        h.write("\n\nPRIMARY DECISION CHECK\n")
        h.write(
            f"- February snow-block partial R2={feb['partial_r2_full_snow_block']:.4f}; "
            f"joint p={feb['joint_snow_block_p']:.4g}; "
            f"snow x brightness p={feb['snow_x_logbrightness_p']:.4g}.\n"
        )
        h.write(
            f"- January snow-block partial R2={jan['partial_r2_full_snow_block']:.4f}; "
            f"joint p={jan['joint_snow_block_p']:.4g}.\n"
        )
        h.write(
            f"- March control snow-block partial R2={mar['partial_r2_full_snow_block']:.4f}; "
            f"joint p={mar['joint_snow_block_p']:.4g}.\n"
        )
        h.write(
            "Support for snow amplification requires a meaningful positive "
            "February snow block after brightness, density and country controls, "
            "a positive snow x brightness interaction, and a materially weaker "
            "March control. Agreement in VCMCFG strengthens robustness. This "
            "diagnostic does not by itself establish causality.\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_monthly_snow_amplification"]
    atlas = load_atlas(cfg)

    products = tuple(str(x).lower() for x in cfg["products"])
    months = tuple(int(x) for x in cfg["months"])
    use = atlas[
        atlas["product"].isin(products) & atlas["month"].isin(months)
    ].copy()

    ids = set(use["NUTS_ID"])
    expected_n = int(cfg["expected_monthly_tile_nuts3"])
    if len(ids) != expected_n:
        raise RuntimeError(
            f"Monthly atlas has {len(ids):,} NUTS3; expected {expected_n:,}."
        )
    for product in products:
        for month in months:
            n = use[
                (use["product"] == product) & (use["month"] == month)
            ]["NUTS_ID"].nunique()
            if n != expected_n:
                raise RuntimeError(
                    f"{product} {month}: {n} NUTS3; expected {expected_n}."
                )

    snow = build_snow_anomaly(config, cfg, ids)
    analysis = prepare(use, snow, config, cfg)
    models, effects, coefficients = fit_models(analysis)
    density, country, deciles = summaries(analysis)

    processed = (
        configured_path(config, "processed")
        / "fixed_nuts2016_monthly_snow_amplification"
    )
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_monthly_snow_amplification"
    )
    processed.mkdir(parents=True, exist_ok=True)

    outputs = {
        "regional": processed / "monthly_snow_radiance_regional.csv",
        "models": processed / "monthly_snow_amplification_models.csv",
        "effects": processed / "monthly_snow_marginal_effects.csv",
        "coefficients": processed / "monthly_snow_amplification_coefficients.csv",
        "density": processed / "monthly_snow_density_summary.csv",
        "country": processed / "monthly_snow_country_summary.csv",
        "deciles": processed / "monthly_snow_baseline_deciles.csv",
        "summary": processed / "monthly_snow_amplification_summary.txt",
    }
    analysis.to_csv(outputs["regional"], index=False)
    models.to_csv(outputs["models"], index=False)
    effects.to_csv(outputs["effects"], index=False)
    coefficients.to_csv(outputs["coefficients"], index=False)
    density.to_csv(outputs["density"], index=False)
    country.to_csv(outputs["country"], index=False)
    deciles.to_csv(outputs["deciles"], index=False)

    make_plots(analysis, models, effects, deciles, figures, int(cfg["figure_dpi"]))
    write_summary(models, effects, density, cfg, outputs["summary"])

    print("\nMONTHLY SNOW-AMPLIFICATION TEST COMPLETE")
    print(models.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for path in [*outputs.values(), figures]:
        print(path)


if __name__ == "__main__":
    main()
