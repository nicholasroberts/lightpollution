#!/usr/bin/env python3
"""Test a general historical winter surface-state control on fixed NUTS2016 trends.

Scientific question
-------------------
Does controlling transient winter snow/surface state materially change the
2014-2024 density-dependent trend in calibrated annual VIIRS DNB radiance?

This is deliberately NOT a special correction for 2021.  The same January-
February MOD10CM-derived covariate is constructed for every available region-
year.  The primary surface-state index is the mean of month-specific snow
anomalies for the available January/February months.  Month-specific anomalies
are centred against the same region/month's other years, so the genuine
February-2016 archive gap can be represented by January alone without treating
January snow as if it were a January-February raw mean.

Two matched-sample model pairs are fitted:
  1. primary: at least one valid Jan/Feb month;
  2. sensitivity: both Jan and Feb valid (therefore excluding 2016).

For each sample, the same model is fitted before and after adding:
    winter snow anomaly
    winter snow anomaly x fixed baseline brightness

The locked fixed-NUTS2016 cohort and population-density classes are unchanged.
Both proportional (log10 radiance -> % yr-1) and absolute (nW cm-2 sr-1 yr-1)
trend estimates are retained.

No annual production data are modified.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress

from analyze_density_temporal_change import contrast, fit_mixedlm, lrt
from diagnose_fixed_nuts2016_2021_reversal import load_locked_panel
from diagnose_fixed_nuts2016_mask_occupancy import attach_fixed_brightness
from project_config import PROJECT_ROOT, configured_path, load_config
from test_fixed_nuts2016_snow_viirs_2021 import (
    ensure_snow_geotiff,
    extract_snow_month,
    fixed_nuts,
    month_from_year_doy,
    normalize_snow_extract,
)


HDF_RE = re.compile(r"^MOD10CM\.A(?P<year>\d{4})(?P<doy>\d{3})\.061\..*\.hdf$")
MONTH_NAMES = {1: "Jan", 2: "Feb"}


def discover_winter_hdfs(
    raw_dir: Path,
    first_year: int,
    last_year: int,
    months: tuple[int, ...],
    allowed_missing: set[tuple[int, int]],
) -> dict[tuple[int, int], Path]:
    found: dict[tuple[int, int], Path] = {}
    for p in sorted(raw_dir.glob("MOD10CM*.hdf")):
        m = HDF_RE.match(p.name)
        if not m:
            continue
        year = int(m.group("year"))
        if not (first_year <= year <= last_year):
            continue
        month = month_from_year_doy(year, int(m.group("doy")))
        if month not in months:
            continue
        key = (year, month)
        if key in found:
            raise RuntimeError(
                f"Duplicate MOD10CM winter month {key}: {found[key]} and {p}"
            )
        found[key] = p

    expected = {
        (year, month)
        for year in range(first_year, last_year + 1)
        for month in months
    }
    missing = expected - set(found)
    unexpected = sorted(missing - allowed_missing)
    if unexpected:
        raise RuntimeError(
            "Historical winter MOD10CM inventory has unexpected gaps: "
            + repr(unexpected)
        )
    extra_allowed = sorted(allowed_missing - missing)
    if extra_allowed:
        print(
            "Previously allowed missing winter month(s) are now present: "
            + ", ".join(f"{y}-{m:02d}" for y, m in extra_allowed),
            flush=True,
        )
    print(
        "Winter MOD10CM inventory: "
        f"{len(found)} files; missing allowed="
        + (
            ", ".join(f"{y}-{m:02d}" for y, m in sorted(missing))
            if missing
            else "none"
        ),
        flush=True,
    )
    return found


def extract_historical_snow(
    config: dict,
    panel: pd.DataFrame,
    cfg: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw_dir = PROJECT_ROOT / str(cfg["mod10cm_raw_dir"])
    if not raw_dir.exists():
        raise FileNotFoundError(raw_dir)

    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    months = tuple(int(x) for x in cfg["winter_months"])
    allowed_missing = {
        tuple(int(v) for v in str(x).split("-"))
        for x in cfg.get("allowed_missing_year_month", [])
    }
    hdfs = discover_winter_hdfs(
        raw_dir, first_year, last_year, months, allowed_missing
    )

    cache_tif = configured_path(config, "processed") / "mod10cm_061_geotiff_cache"
    cache_zonal = configured_path(config, "processed") / "mod10cm_061_nuts2016_cache"

    nuts = fixed_nuts(config)
    id_field = config["nuts"]["id_field"]
    nuts[id_field] = nuts[id_field].astype(str)
    keep_ids = set(panel["NUTS_ID"].astype(str).unique())
    nuts = nuts[nuts[id_field].isin(keep_ids)].copy()
    if set(nuts[id_field]) != keep_ids:
        missing_ids = sorted(keep_ids - set(nuts[id_field]))
        raise RuntimeError(
            "Locked panel contains NUTS IDs absent from fixed geometry: "
            + repr(missing_ids[:20])
        )

    frames = []
    manifest_rows = []
    for (year, month), hdf in sorted(hdfs.items()):
        tif = ensure_snow_geotiff(hdf, cache_tif, config["nuts"]["crs"])
        print(f"Extracting historical MOD10CM {year}-{month:02d}", flush=True)
        z = extract_snow_month(
            tif, year, month, nuts, config, cache_zonal
        )
        frames.append(z)
        manifest_rows.append(
            {
                "year": year,
                "month": month,
                "month_name": MONTH_NAMES[month],
                "hdf_file": hdf.name,
                "hdf_size_bytes": int(hdf.stat().st_size),
                "cached_geotiff": tif.name,
            }
        )

    snow = normalize_snow_extract(pd.concat(frames, ignore_index=True))
    snow["NUTS_ID"] = snow["NUTS_ID"].astype(str)
    return snow, pd.DataFrame(manifest_rows)


def add_leave_one_out_month_anomaly(
    snow: pd.DataFrame,
    min_valid_area: float,
    min_climatology_years: int,
) -> pd.DataFrame:
    d = snow.copy()
    d["eligible_month"] = (
        d["snow_mean_pct"].notna()
        & d["snow_valid_area_fraction"].notna()
        & (d["snow_valid_area_fraction"] >= min_valid_area)
    )
    d["snow_for_climatology"] = d["snow_mean_pct"].where(d["eligible_month"])

    stats = (
        d.groupby(["NUTS_ID", "month"], observed=True)["snow_for_climatology"]
        .agg(["sum", "count"])
        .rename(columns={"sum": "snow_month_sum", "count": "snow_month_n"})
        .reset_index()
    )
    d = d.merge(
        stats,
        on=["NUTS_ID", "month"],
        how="left",
        validate="many_to_one",
    )

    d["snow_climatology_loo_pct"] = np.nan
    ok = (
        d["eligible_month"]
        & (d["snow_month_n"] - 1 >= int(min_climatology_years))
    )
    d.loc[ok, "snow_climatology_loo_pct"] = (
        d.loc[ok, "snow_month_sum"] - d.loc[ok, "snow_mean_pct"]
    ) / (d.loc[ok, "snow_month_n"] - 1)

    d["snow_month_anomaly_pp"] = (
        d["snow_mean_pct"] - d["snow_climatology_loo_pct"]
    )
    d.loc[~ok, "snow_month_anomaly_pp"] = np.nan

    finite = d["snow_month_anomaly_pp"].dropna()
    if finite.empty:
        raise RuntimeError("No valid historical month-specific snow anomalies.")
    return d


def build_winter_index(
    snow_monthly: pd.DataFrame,
    panel: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    d = add_leave_one_out_month_anomaly(
        snow_monthly,
        float(cfg["min_snow_valid_area_fraction"]),
        int(cfg["min_climatology_other_years"]),
    )

    valid = d[d["snow_month_anomaly_pp"].notna()].copy()
    winter = (
        valid.groupby(["NUTS_ID", "year"], observed=True)
        .agg(
            winter_snow_anomaly_pp=("snow_month_anomaly_pp", "mean"),
            winter_snow_mean_pct=("snow_mean_pct", "mean"),
            winter_months_available=("month", "nunique"),
            winter_min_valid_area_fraction=("snow_valid_area_fraction", "min"),
            winter_mean_valid_area_fraction=("snow_valid_area_fraction", "mean"),
        )
        .reset_index()
    )

    # Preserve explicit Jan/Feb values and anomaly values for audit.
    for value, prefix in [
        ("snow_mean_pct", "snow_pct"),
        ("snow_month_anomaly_pp", "snow_anomaly_pp"),
        ("snow_valid_area_fraction", "snow_valid_fraction"),
    ]:
        p = d.pivot(index=["NUTS_ID", "year"], columns="month", values=value)
        p = p.rename(
            columns={1: f"{prefix}_jan", 2: f"{prefix}_feb"}
        ).reset_index()
        winter = winter.merge(
            p, on=["NUTS_ID", "year"], how="outer", validate="one_to_one"
        )

    grid = panel[["NUTS_ID", "year"]].drop_duplicates().copy()
    grid["NUTS_ID"] = grid["NUTS_ID"].astype(str)
    winter = grid.merge(
        winter, on=["NUTS_ID", "year"], how="left", validate="one_to_one"
    )
    winter["winter_months_available"] = (
        winter["winter_months_available"].fillna(0).astype(int)
    )
    winter["eligible_primary"] = (
        winter["winter_months_available"]
        >= int(cfg["primary_min_winter_months"])
    )
    winter["eligible_complete_janfeb"] = (
        winter["winter_months_available"]
        >= int(cfg["sensitivity_min_winter_months"])
    )
    return winter.sort_values(["NUTS_ID", "year"]).reset_index(drop=True)


def prepare_model_panel(
    panel: pd.DataFrame,
    winter: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    d = attach_fixed_brightness(panel).merge(
        winter,
        on=["NUTS_ID", "year"],
        how="left",
        validate="one_to_one",
    )
    d["log_rad"] = np.log10(d["radiance_mean_corrected"].astype(float))
    midpoint = (int(cfg["first_year"]) + int(cfg["last_year"])) / 2.0
    d["year_c"] = d["year"].astype(float) - midpoint
    d["country"] = d["CNTR_CODE"].astype("category")

    log_base = np.log10(d["baseline_raw_nw"].astype(float))
    centre = float(
        d[["NUTS_ID"]]
        .assign(log_base=log_base)
        .drop_duplicates("NUTS_ID")["log_base"]
        .median()
    )
    d["log10_baseline_c"] = log_base - centre
    d.attrs["log10_baseline_center"] = centre

    if d["density_band"].isna().any():
        raise RuntimeError("Density band missing after winter-index merge.")
    return d


def _trend_contrast(result, band: str, labels: list[str]):
    weights = {"year_c": 1.0}
    if band != labels[0]:
        candidates = [
            f"year_c:C(density_band)[T.{band}]",
            f"C(density_band)[T.{band}]:year_c",
        ]
        name = next((x for x in candidates if x in result.fe_params.index), None)
        if name is None:
            raise RuntimeError(
                f"Could not find density-band trend term for {band}; "
                f"parameters={list(result.fe_params.index)}"
            )
        weights[name] = 1.0
    return contrast(result, weights)


def _snow_term_name(result) -> tuple[str, str]:
    snow = "winter_snow_anomaly_pp"
    if snow not in result.fe_params.index:
        raise RuntimeError("Controlled model is missing the snow main effect.")
    candidates = [
        "winter_snow_anomaly_pp:log10_baseline_c",
        "log10_baseline_c:winter_snow_anomaly_pp",
    ]
    inter = next((x for x in candidates if x in result.fe_params.index), None)
    if inter is None:
        raise RuntimeError("Controlled model is missing snow x brightness.")
    return snow, inter


def fit_model_pair(
    data: pd.DataFrame,
    outcome: str,
    outcome_scale: str,
    sample_name: str,
    labels: list[str],
) -> tuple[pd.DataFrame, dict, object, object, pd.DataFrame]:
    cols = [
        outcome,
        "NUTS_ID",
        "CNTR_CODE",
        "country",
        "density_band",
        "year",
        "year_c",
        "winter_snow_anomaly_pp",
        "log10_baseline_c",
        "baseline_raw_nw",
    ]
    d = (
        data[cols]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .copy()
    )
    if d.empty:
        raise RuntimeError(f"No rows for {sample_name}/{outcome_scale}.")
    if d["NUTS_ID"].nunique() < 100:
        raise RuntimeError(
            f"Unexpectedly small historical snow-control cohort: "
            f"{d['NUTS_ID'].nunique()} regions."
        )

    base_formula = (
        f"{outcome} ~ year_c * C(density_band) + C(country)"
    )
    controlled_formula = (
        base_formula
        + " + winter_snow_anomaly_pp * log10_baseline_c"
    )

    base, opt0 = fit_mixedlm(base_formula, d)
    controlled, opt1 = fit_mixedlm(controlled_formula, d)
    lr, df, p = lrt(base, controlled)

    rows = []
    for band in labels:
        b_est, b_se, b_p = _trend_contrast(base, band, labels)
        c_est, c_se, c_p = _trend_contrast(controlled, band, labels)
        region = d[d["density_band"].astype(str) == band][
            ["NUTS_ID", "baseline_raw_nw"]
        ].drop_duplicates("NUTS_ID")
        baseline_median = float(region["baseline_raw_nw"].median())

        row = {
            "sample": sample_name,
            "outcome_scale": outcome_scale,
            "density_band": band,
            "n_observations": int(
                (d["density_band"].astype(str) == band).sum()
            ),
            "n_nuts3": int(
                d.loc[d["density_band"].astype(str) == band, "NUTS_ID"].nunique()
            ),
            "median_baseline_raw_nw": baseline_median,
            "trend_unadjusted": b_est,
            "trend_unadjusted_se": b_se,
            "trend_unadjusted_p": b_p,
            "trend_controlled": c_est,
            "trend_controlled_se": c_se,
            "trend_controlled_p": c_p,
            "trend_shift_controlled_minus_unadjusted": c_est - b_est,
        }
        if outcome_scale == "proportional":
            row.update(
                {
                    "pct_per_year_unadjusted": 100.0 * (10.0 ** b_est - 1.0),
                    "pct_per_year_controlled": 100.0 * (10.0 ** c_est - 1.0),
                    "pct_per_year_unadjusted_ci_low": 100.0
                    * (10.0 ** (b_est - 1.96 * b_se) - 1.0),
                    "pct_per_year_unadjusted_ci_high": 100.0
                    * (10.0 ** (b_est + 1.96 * b_se) - 1.0),
                    "pct_per_year_controlled_ci_low": 100.0
                    * (10.0 ** (c_est - 1.96 * c_se) - 1.0),
                    "pct_per_year_controlled_ci_high": 100.0
                    * (10.0 ** (c_est + 1.96 * c_se) - 1.0),
                }
            )
            row["trend_shift_display"] = (
                row["pct_per_year_controlled"]
                - row["pct_per_year_unadjusted"]
            )
            row["trend_shift_display_units"] = "percentage_points_per_year"
        else:
            row.update(
                {
                    "nw_per_year_unadjusted": b_est,
                    "nw_per_year_controlled": c_est,
                    "nw_per_year_unadjusted_ci_low": b_est - 1.96 * b_se,
                    "nw_per_year_unadjusted_ci_high": b_est + 1.96 * b_se,
                    "nw_per_year_controlled_ci_low": c_est - 1.96 * c_se,
                    "nw_per_year_controlled_ci_high": c_est + 1.96 * c_se,
                    "unadjusted_nw_trend_pct_of_baseline_per_year": (
                        100.0 * b_est / baseline_median
                        if baseline_median > 0 else np.nan
                    ),
                    "controlled_nw_trend_pct_of_baseline_per_year": (
                        100.0 * c_est / baseline_median
                        if baseline_median > 0 else np.nan
                    ),
                }
            )
            row["trend_shift_display"] = c_est - b_est
            row["trend_shift_display_units"] = "nW_per_year"
            row["normalized_absolute_shift_pct_of_baseline_per_year"] = (
                100.0 * (c_est - b_est) / baseline_median
                if baseline_median > 0 else np.nan
            )
        rows.append(row)

    snow_name, int_name = _snow_term_name(controlled)
    info = {
        "sample": sample_name,
        "outcome_scale": outcome_scale,
        "n_observations": int(len(d)),
        "n_nuts3": int(d["NUTS_ID"].nunique()),
        "n_countries": int(d["CNTR_CODE"].nunique()),
        "n_single_month_rows": int(
            (data.loc[d.index, "winter_months_available"] == 1).sum()
        ),
        "snow_block_lr": lr,
        "snow_block_df": df,
        "snow_block_p": p,
        "aic_unadjusted": float(base.aic),
        "aic_controlled": float(controlled.aic),
        "delta_aic_controlled_minus_unadjusted": float(
            controlled.aic - base.aic
        ),
        "snow_main_coef": float(controlled.fe_params[snow_name]),
        "snow_main_p": float(controlled.pvalues[snow_name]),
        "snow_x_brightness_coef": float(controlled.fe_params[int_name]),
        "snow_x_brightness_p": float(controlled.pvalues[int_name]),
        "unadjusted_optimizer": opt0,
        "controlled_optimizer": opt1,
        "unadjusted_formula": base_formula,
        "controlled_formula": controlled_formula,
    }

    residuals = d[
        ["NUTS_ID", "CNTR_CODE", "density_band", "year"]
    ].copy()
    residuals["sample"] = sample_name
    residuals["outcome_scale"] = outcome_scale
    residuals["residual_unadjusted"] = np.asarray(base.resid, dtype=float)
    residuals["residual_controlled"] = np.asarray(controlled.resid, dtype=float)

    return pd.DataFrame(rows), info, base, controlled, d, residuals


def summarise_year_residuals(residual_rows: list[pd.DataFrame]) -> pd.DataFrame:
    allr = pd.concat(residual_rows, ignore_index=True)
    rows = []
    for (sample, scale, year), g in allr.groupby(
        ["sample", "outcome_scale", "year"], observed=True
    ):
        countries = (
            g.groupby("CNTR_CODE", observed=True)[
                ["residual_unadjusted", "residual_controlled"]
            ]
            .median()
            .reset_index()
        )
        row = {
            "sample": sample,
            "outcome_scale": scale,
            "year": int(year),
            "n_nuts3": int(g["NUTS_ID"].nunique()),
            "n_countries": int(g["CNTR_CODE"].nunique()),
            "median_region_residual_unadjusted": float(
                g["residual_unadjusted"].median()
            ),
            "median_region_residual_controlled": float(
                g["residual_controlled"].median()
            ),
            "country_balanced_median_residual_unadjusted": float(
                countries["residual_unadjusted"].median()
            ),
            "country_balanced_median_residual_controlled": float(
                countries["residual_controlled"].median()
            ),
        }
        if scale == "proportional":
            for key in [
                "median_region_residual_unadjusted",
                "median_region_residual_controlled",
                "country_balanced_median_residual_unadjusted",
                "country_balanced_median_residual_controlled",
            ]:
                row[key.replace("residual", "residual_pct")] = (
                    100.0 * (10.0 ** row[key] - 1.0)
                )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["sample", "outcome_scale", "year"]
    )


def adjusted_region_slopes(
    controlled_models: dict[tuple[str, str], tuple[object, pd.DataFrame]],
    labels: list[str],
) -> pd.DataFrame:
    # Region-level paired slopes are reported for the primary sample only.
    rows_by_scale: dict[str, pd.DataFrame] = {}
    for scale in ("proportional", "absolute_nw"):
        result, d = controlled_models[("available_winter", scale)]
        snow_name, int_name = _snow_term_name(result)
        beta_snow = float(result.fe_params[snow_name])
        beta_int = float(result.fe_params[int_name])

        x = d.copy()
        x["snow_component"] = (
            beta_snow * x["winter_snow_anomaly_pp"]
            + beta_int
            * x["winter_snow_anomaly_pp"]
            * x["log10_baseline_c"]
        )
        outcome = "log_rad" if scale == "proportional" else "radiance_mean_corrected"
        x["surface_adjusted_outcome"] = x[outcome] - x["snow_component"]

        rr = []
        for nuts_id, g in x.groupby("NUTS_ID", observed=True):
            g = g.sort_values("year")
            if len(g) < 5:
                continue
            obs = linregress(g["year"], g[outcome])
            adj = linregress(g["year"], g["surface_adjusted_outcome"])
            one = {
                "NUTS_ID": nuts_id,
                "CNTR_CODE": str(g["CNTR_CODE"].iloc[0]),
                "density_band": str(g["density_band"].iloc[0]),
                "baseline_raw_nw": float(g["baseline_raw_nw"].iloc[0]),
                "n_years": int(len(g)),
                "observed_slope": float(obs.slope),
                "surface_adjusted_slope": float(adj.slope),
            }
            if scale == "proportional":
                one["observed_pct_per_year"] = 100.0 * (10.0 ** obs.slope - 1.0)
                one["surface_adjusted_pct_per_year"] = (
                    100.0 * (10.0 ** adj.slope - 1.0)
                )
            else:
                one["observed_nw_per_year"] = float(obs.slope)
                one["surface_adjusted_nw_per_year"] = float(adj.slope)
            rr.append(one)
        rows_by_scale[scale] = pd.DataFrame(rr)

    a = rows_by_scale["proportional"].rename(
        columns={
            "observed_slope": "observed_log10_slope_per_year",
            "surface_adjusted_slope": "surface_adjusted_log10_slope_per_year",
            "n_years": "n_years_proportional",
        }
    )
    b = rows_by_scale["absolute_nw"][
        [
            "NUTS_ID",
            "n_years",
            "observed_nw_per_year",
            "surface_adjusted_nw_per_year",
        ]
    ].rename(columns={"n_years": "n_years_absolute"})
    return a.merge(b, on="NUTS_ID", how="inner", validate="one_to_one")


def snow_by_year_density(data: pd.DataFrame) -> pd.DataFrame:
    g = data[data["eligible_primary"]].copy()
    rows = []
    for (year, band), x in g.groupby(
        ["year", "density_band"], observed=True
    ):
        rows.append(
            {
                "year": int(year),
                "density_band": str(band),
                "n_nuts3": int(x["NUTS_ID"].nunique()),
                "median_winter_snow_anomaly_pp": float(
                    x["winter_snow_anomaly_pp"].median()
                ),
                "median_winter_snow_mean_pct": float(
                    x["winter_snow_mean_pct"].median()
                ),
                "pct_single_month": float(
                    100.0 * (x["winter_months_available"] == 1).mean()
                ),
            }
        )
    return pd.DataFrame(rows)


def apply_materiality(
    trends: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    d = trends.copy()
    pct_threshold = float(cfg["material_pct_shift_pp_per_year"])
    relative_threshold = float(cfg["material_relative_shift_fraction"])
    baseline_threshold = float(
        cfg["material_absolute_shift_pct_of_baseline_per_year"]
    )

    d["relative_shift_fraction"] = np.nan
    d["material_shift"] = False

    prop = d["outcome_scale"] == "proportional"
    if prop.any():
        un = d.loc[prop, "pct_per_year_unadjusted"].astype(float)
        co = d.loc[prop, "pct_per_year_controlled"].astype(float)
        shift = co - un
        rel = np.where(
            np.abs(un) >= pct_threshold,
            np.abs(shift) / np.abs(un),
            np.nan,
        )
        sign_change = (
            (np.sign(un) != np.sign(co))
            & (np.maximum(np.abs(un), np.abs(co)) >= pct_threshold)
        )
        material = (
            (np.abs(shift) >= pct_threshold)
            | (rel >= relative_threshold)
            | sign_change
        )
        d.loc[prop, "relative_shift_fraction"] = rel
        d.loc[prop, "material_shift"] = material

    absolute = d["outcome_scale"] == "absolute_nw"
    if absolute.any():
        un = d.loc[absolute, "nw_per_year_unadjusted"].astype(float)
        co = d.loc[absolute, "nw_per_year_controlled"].astype(float)
        shift = co - un
        rel = np.where(np.abs(un) > 1e-12, np.abs(shift) / np.abs(un), np.nan)
        normalized = d.loc[
            absolute, "normalized_absolute_shift_pct_of_baseline_per_year"
        ].astype(float)
        sign_change = (
            (np.sign(un) != np.sign(co))
            & (np.abs(normalized) >= baseline_threshold)
        )
        material = (
            (np.abs(normalized) >= baseline_threshold)
            | (
                (rel >= relative_threshold)
                & (
                    np.abs(
                        d.loc[
                            absolute,
                            "unadjusted_nw_trend_pct_of_baseline_per_year",
                        ].astype(float)
                    )
                    >= baseline_threshold
                )
            )
            | sign_change
        )
        d.loc[absolute, "relative_shift_fraction"] = rel
        d.loc[absolute, "material_shift"] = material

    return d


def decision_summary(
    trends: pd.DataFrame,
    tests: pd.DataFrame,
    cfg: dict,
) -> dict:
    primary = trends[trends["sample"] == "available_winter"]
    test_primary = tests[tests["sample"] == "available_winter"]
    alpha = float(cfg["selection_alpha"])
    min_classes = int(cfg["material_min_density_classes"])

    n_prop = int(
        primary.loc[
            primary["outcome_scale"] == "proportional", "material_shift"
        ].sum()
    )
    n_abs = int(
        primary.loc[
            primary["outcome_scale"] == "absolute_nw", "material_shift"
        ].sum()
    )
    significant_block = bool(
        (test_primary["snow_block_p"].astype(float) < alpha).any()
    )
    proceed = significant_block and max(n_prop, n_abs) >= min_classes

    return {
        "snow_block_significant_any_primary_scale": significant_block,
        "material_classes_proportional": n_prop,
        "material_classes_absolute": n_abs,
        "material_min_density_classes": min_classes,
        "selection_alpha": alpha,
        "historical_albedo_warranted_by_preregistered_rule": proceed,
    }


def make_plots(
    trends: pd.DataFrame,
    residuals: pd.DataFrame,
    snow_year: pd.DataFrame,
    labels: list[str],
    outdir: Path,
    dpi: int,
) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = []

    primary = trends[trends["sample"] == "available_winter"].copy()
    order = {b: i for i, b in enumerate(labels)}
    primary["_order"] = primary["density_band"].map(order)
    primary = primary.sort_values("_order")

    p = primary[primary["outcome_scale"] == "proportional"]
    x = np.arange(len(labels), dtype=float)
    fig, ax = plt.subplots(figsize=(9.3, 5.9))
    for shift, suffix, marker, label in [
        (-0.12, "unadjusted", "o", "Unadjusted"),
        (0.12, "controlled", "s", "Winter snow controlled"),
    ]:
        y = p[f"pct_per_year_{suffix}"].to_numpy(float)
        lo = p[f"pct_per_year_{suffix}_ci_low"].to_numpy(float)
        hi = p[f"pct_per_year_{suffix}_ci_high"].to_numpy(float)
        ax.errorbar(
            x + shift, y, yerr=[y - lo, hi - y],
            fmt=marker, capsize=4, label=label
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Fixed population-density class (people km$^{-2}$)")
    ax.set_ylabel("Calibrated VIIRS radiance trend (% yr$^{-1}$)")
    ax.set_title("Historical density trends before and after winter snow control")
    ax.grid(alpha=0.2, axis="y")
    ax.legend()
    fig.tight_layout()
    path = outdir / "proportional_trends_before_after_snow.pdf"
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    outputs.append(path)

    a = primary[primary["outcome_scale"] == "absolute_nw"]
    fig, ax = plt.subplots(figsize=(9.3, 5.9))
    for shift, suffix, marker, label in [
        (-0.12, "unadjusted", "o", "Unadjusted"),
        (0.12, "controlled", "s", "Winter snow controlled"),
    ]:
        y = a[f"nw_per_year_{suffix}"].to_numpy(float)
        lo = a[f"nw_per_year_{suffix}_ci_low"].to_numpy(float)
        hi = a[f"nw_per_year_{suffix}_ci_high"].to_numpy(float)
        ax.errorbar(
            x + shift, y, yerr=[y - lo, hi - y],
            fmt=marker, capsize=4, label=label
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Fixed population-density class (people km$^{-2}$)")
    ax.set_ylabel("Calibrated VIIRS radiance trend (nW cm$^{-2}$ sr$^{-1}$ yr$^{-1}$)")
    ax.set_title("Absolute historical trends before and after winter snow control")
    ax.grid(alpha=0.2, axis="y")
    ax.legend()
    fig.tight_layout()
    path = outdir / "absolute_trends_before_after_snow.pdf"
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    outputs.append(path)

    r = residuals[
        (residuals["sample"] == "available_winter")
        & (residuals["outcome_scale"] == "proportional")
    ].sort_values("year")
    fig, ax = plt.subplots(figsize=(9.3, 5.8))
    ax.plot(
        r["year"],
        r["country_balanced_median_residual_pct_unadjusted"],
        marker="o",
        label="Unadjusted",
    )
    ax.plot(
        r["year"],
        r["country_balanced_median_residual_pct_controlled"],
        marker="s",
        label="Winter snow controlled",
    )
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Country-balanced residual from fitted trend (%)")
    ax.set_title("Annual departures from density-class trend")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    path = outdir / "annual_proportional_residuals_before_after_snow.pdf"
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    outputs.append(path)

    sy = snow_year.copy()
    fig, ax = plt.subplots(figsize=(9.5, 5.9))
    for band in labels:
        g = sy[sy["density_band"] == band].sort_values("year")
        ax.plot(
            g["year"],
            g["median_winter_snow_anomaly_pp"],
            marker="o",
            label=band,
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Median Jan-Feb snow anomaly (percentage points)")
    ax.set_title("Historical winter surface-state anomaly by density class")
    ax.grid(alpha=0.2)
    ax.legend(title="People km$^{-2}$")
    fig.tight_layout()
    path = outdir / "winter_snow_anomaly_by_year_density.pdf"
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    outputs.append(path)

    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_historical_snow_control"]
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )

    panel = load_locked_panel(config)
    expected_years = list(
        range(int(cfg["first_year"]), int(cfg["last_year"]) + 1)
    )
    observed_years = sorted(panel["year"].astype(int).unique())
    if observed_years != expected_years:
        raise RuntimeError(
            f"Historical snow-control stage expects years {expected_years}; "
            f"locked panel has {observed_years}."
        )

    snow_monthly, manifest = extract_historical_snow(
        config, panel, cfg
    )
    winter = build_winter_index(snow_monthly, panel, cfg)
    model_panel = prepare_model_panel(panel, winter, cfg)

    # Explicit availability audit before modelling.
    availability = (
        winter.groupby("year", observed=True)
        .agg(
            n_nuts3=("NUTS_ID", "nunique"),
            n_primary=("eligible_primary", "sum"),
            n_complete_janfeb=("eligible_complete_janfeb", "sum"),
            n_zero_months=("winter_months_available", lambda s: int((s == 0).sum())),
            n_one_month=("winter_months_available", lambda s: int((s == 1).sum())),
            n_two_months=("winter_months_available", lambda s: int((s >= 2).sum())),
        )
        .reset_index()
    )

    trend_rows = []
    test_rows = []
    residual_rows = []
    controlled_models = {}

    samples = [
        ("available_winter", "eligible_primary"),
        ("complete_janfeb", "eligible_complete_janfeb"),
    ]
    outcomes = [
        ("log_rad", "proportional"),
        ("radiance_mean_corrected", "absolute_nw"),
    ]

    for sample_name, eligibility_col in samples:
        sample = model_panel[model_panel[eligibility_col]].copy()
        print(
            f"\nFitting sample={sample_name}: {len(sample)} rows, "
            f"{sample['NUTS_ID'].nunique()} NUTS3",
            flush=True,
        )
        for outcome, scale in outcomes:
            trends, info, base, controlled, d, residuals = fit_model_pair(
                sample, outcome, scale, sample_name, labels
            )
            trend_rows.append(trends)
            test_rows.append(info)
            residual_rows.append(residuals)
            controlled_models[(sample_name, scale)] = (controlled, d)

    trends = apply_materiality(pd.concat(trend_rows, ignore_index=True), cfg)
    tests = pd.DataFrame(test_rows)
    year_residuals = summarise_year_residuals(residual_rows)
    region_slopes = adjusted_region_slopes(controlled_models, labels)
    snow_year = snow_by_year_density(model_panel)
    decision = decision_summary(trends, tests, cfg)

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_historical_snow_control"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    outputs = {
        "manifest": processed / "fixed_nuts2016_historical_snow_manifest.csv",
        "monthly": processed / "fixed_nuts2016_mod10cm_janfeb_2014_2024.csv",
        "winter": processed / "fixed_nuts2016_winter_surface_state_2014_2024.csv",
        "availability": processed / "fixed_nuts2016_winter_surface_state_availability.csv",
        "trends": processed / "fixed_nuts2016_snow_control_trend_comparison.csv",
        "tests": processed / "fixed_nuts2016_snow_control_model_tests.csv",
        "region_slopes": processed / "fixed_nuts2016_snow_control_region_slopes.csv",
        "year_residuals": processed / "fixed_nuts2016_snow_control_year_residuals.csv",
        "snow_year": processed / "fixed_nuts2016_snow_anomaly_by_year_density.csv",
        "summary": processed / "fixed_nuts2016_historical_snow_control_summary.txt",
    }

    manifest.to_csv(outputs["manifest"], index=False)
    snow_monthly.to_csv(outputs["monthly"], index=False)
    winter.to_csv(outputs["winter"], index=False)
    availability.to_csv(outputs["availability"], index=False)
    trends.to_csv(outputs["trends"], index=False)
    tests.to_csv(outputs["tests"], index=False)
    region_slopes.to_csv(outputs["region_slopes"], index=False)
    year_residuals.to_csv(outputs["year_residuals"], index=False)
    snow_year.to_csv(outputs["snow_year"], index=False)

    figure_paths = make_plots(
        trends,
        year_residuals,
        snow_year,
        labels,
        figures,
        int(cfg["figure_dpi"]),
    )

    primary_trends = trends[trends["sample"] == "available_winter"].copy()
    primary_tests = tests[tests["sample"] == "available_winter"].copy()
    r2021 = year_residuals[
        (year_residuals["sample"] == "available_winter")
        & (year_residuals["year"] == 2021)
    ]

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 HISTORICAL WINTER SURFACE-STATE CONTROL\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            "Question: does a general Jan-Feb snow/surface-state covariate "
            "materially change the 2014-2024 density-dependent annual VIIRS "
            "trend? 2021 receives no special correction.\n"
        )
        h.write(
            "Primary snow index: mean available Jan/Feb month-specific anomaly, "
            "where each month is centred against that region/month's other "
            "years. This makes the genuine 2016-02 archive gap auditable without "
            "substituting a raw January value for a raw two-month mean.\n"
        )
        h.write(
            "Primary sample: >=1 eligible winter month. Sensitivity: both Jan "
            "and Feb eligible.\n"
        )
        h.write(
            "Both proportional log-radiance and absolute nW trends are primary "
            "reporting spaces.\n\n"
        )

        h.write("WINTER DATA AVAILABILITY\n")
        h.write(availability.to_string(index=False))
        h.write("\n\nPRIMARY TREND COMPARISON\n")
        h.write(
            primary_trends.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\nPRIMARY MODEL TESTS\n")
        h.write(
            primary_tests.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\n2021 RESIDUAL CHECK\n")
        h.write(
            r2021.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
            if len(r2021)
            else "(no 2021 residual row)"
        )

        h.write("\n\nPRE-REGISTERED MATERIALITY RULE\n")
        h.write(
            f"- Proportional class shift threshold: "
            f"{cfg['material_pct_shift_pp_per_year']} percentage points yr-1, "
            f"or {100*float(cfg['material_relative_shift_fraction']):.0f}% "
            "relative change in a non-trivial trend, or a non-trivial sign change.\n"
        )
        h.write(
            f"- Absolute class shift threshold: "
            f"{cfg['material_absolute_shift_pct_of_baseline_per_year']}% of "
            "class median baseline radiance per year, with the same relative/sign "
            "safeguards.\n"
        )
        h.write(
            f"- Historical albedo is warranted only if the primary snow block "
            f"is significant at alpha={cfg['selection_alpha']} in at least one "
            f"reporting space AND at least {cfg['material_min_density_classes']} "
            "density classes shift materially in proportional or absolute space.\n"
        )

        h.write("\nDECISION\n")
        for key, value in decision.items():
            h.write(f"- {key}: {value}\n")
        h.write(
            "\nThis stage is diagnostic. It does not overwrite annual VIIRS, "
            "remove 2021, or create a corrected production time series.\n"
        )

    print("\nHISTORICAL WINTER SURFACE-STATE CONTROL")
    print("\nPRIMARY MODEL TESTS")
    print(primary_tests.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nPRIMARY DENSITY-CLASS TREND COMPARISON")
    cols = [
        "outcome_scale",
        "density_band",
        "pct_per_year_unadjusted",
        "pct_per_year_controlled",
        "nw_per_year_unadjusted",
        "nw_per_year_controlled",
        "trend_shift_display",
        "material_shift",
    ]
    cols = [c for c in cols if c in primary_trends.columns]
    print(
        primary_trends[cols].to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print("\nPRE-REGISTERED DECISION")
    for key, value in decision.items():
        print(f"{key}: {value}")
    print("\nWROTE")
    for p in [*outputs.values(), *figure_paths]:
        print(p)


if __name__ == "__main__":
    main()
