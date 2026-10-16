#!/usr/bin/env python3
"""Build a physically informed 2021 surface-state correction for annual VIIRS.

Purpose
-------
This stage converts the validated Jan/Feb 2021 snow/albedo disturbance into a
counterfactual annual 2021 NUTS3 radiance and immediately re-fits the locked
fixed-NUTS2016 historical trends.

The design is deliberately correction-oriented rather than diagnostic:

1. Refit the validated Jan/Feb WSA QA0/1 absolute-radiance model.
2. For each NUTS3-month, predict the surface-state contribution by setting the
   2021 snow and albedo anomalies to zero while holding all controls fixed.
3. Bridge the common-support monthly outcome used by the albedo model to the
   complete monthly VCMSLCFG atlas.
4. Derive annual month weights from the actual monthly VCMSLCFG cf_cvg counts
   (not an unweighted mean of monthly radiance).
5. Validate the reconstructed monthly-to-annual mapping against the observed
   2020 and 2022 Annual VNL fixed-NUTS2016 products.
6. Propagate the Jan/Feb surface-state component into annual raw and calibrated
   radiance using the validated control-year mapping.
7. Change only 2021 in a copy of the fixed 2014-2024 panel and re-fit the same
   density-dependent historical model.

Primary correction:
    combined snow + WSA visible albedo model where albedo is eligible,
    snow-only fallback where albedo is unavailable.

Two sensitivity corrections are also retained (snow-only and albedo-primary)
so the final historical trend can be checked without launching another
mechanism-search workflow.

No existing production file is overwritten.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from exactextract.raster import RasterSource
from rasterio.windows import Window
from scipy.stats import pearsonr
import statsmodels.formula.api as smf

from diagnose_fixed_nuts2016_monthly_product_atlas import discover_inputs
from project_config import configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import (
    fit_window,
    load_fixed_nuts,
    observed_trajectories,
    prepare_analysis_panel,
)


MONTHS_TO_CORRECT = (1, 2)
MONTH_NAMES = {1: "Jan", 2: "Feb"}
YEARS_RECON = (2020, 2021, 2022)


class RadianceCoverageSource(RasterSource):
    """Area-extraction source for cf_cvg and avg_rade9h * cf_cvg."""

    def __init__(self, radiance: Path, cloudfree: Path, mode: str, name: str):
        super().__init__()
        if mode not in {"cf", "radcf"}:
            raise ValueError(mode)
        self.rad_path = Path(radiance)
        self.cf_path = Path(cloudfree)
        self.rad = rasterio.open(self.rad_path)
        self.cf = rasterio.open(self.cf_path)
        self.mode = mode
        if self.rad.crs != self.cf.crs or self.rad.transform != self.cf.transform:
            raise RuntimeError(
                f"Monthly VIIRS grid mismatch: {self.rad_path} / {self.cf_path}"
            )
        if (self.rad.width, self.rad.height) != (self.cf.width, self.cf.height):
            raise RuntimeError("Monthly VIIRS shape mismatch.")
        self.set_name(name)

    def res(self):
        return (abs(self.rad.transform.a), abs(self.rad.transform.e))

    def extent(self):
        b = self.rad.bounds
        return (b.left, b.bottom, b.right, b.top)

    def srs_wkt(self):
        return self.rad.crs.wkt

    def nodata_value(self):
        return None

    def read_window(self, x0, y0, nx, ny):
        if nx == 0 or ny == 0:
            return np.empty((0, 0), dtype=np.float32)
        w = Window(int(x0), int(y0), int(nx), int(ny))
        rad = self.rad.read(1, window=w, masked=True)
        cf = self.cf.read(1, window=w, masked=True)
        r = np.asarray(rad.data, dtype=np.float64)
        c = np.asarray(cf.data, dtype=np.float64)
        valid = (
            ~np.ma.getmaskarray(rad)
            & ~np.ma.getmaskarray(cf)
            & np.isfinite(r)
            & np.isfinite(c)
            & (c > 0)
        )
        out = np.zeros(r.shape, dtype=np.float32)
        if self.mode == "cf":
            out[valid] = c[valid].astype(np.float32)
        else:
            out[valid] = (r[valid] * c[valid]).astype(np.float32)
        return out

    def close(self):
        self.rad.close()
        self.cf.close()


def load_inputs(config: dict, cfg: dict):
    processed = configured_path(config, "processed")
    annual_path = processed / str(cfg["annual_fixed_file"])
    albedo_path = processed / str(cfg["albedo_analysis_file"])
    atlas_region_path = (
        processed
        / "fixed_nuts2016_monthly_product_atlas"
        / "monthly_region_2021_vs_flanks.csv"
    )
    atlas_monthly_path = (
        processed
        / "fixed_nuts2016_monthly_product_atlas"
        / "monthly_nuts2016_radiance_2020_2022.parquet"
    )
    reference_path = processed / "fixed_nuts2016_population_density_reference.csv"

    for p in (
        annual_path, albedo_path, atlas_region_path, atlas_monthly_path, reference_path
    ):
        if not p.exists():
            raise FileNotFoundError(p)

    annual = pd.read_parquet(annual_path)
    albedo = pd.read_csv(albedo_path)
    atlas_region = pd.read_csv(atlas_region_path)
    atlas_monthly = pd.read_parquet(atlas_monthly_path)
    reference = pd.read_csv(reference_path)

    for d in (annual, albedo, atlas_region, atlas_monthly, reference):
        if "NUTS_ID" in d.columns:
            d["NUTS_ID"] = d["NUTS_ID"].astype(str)

    return annual, albedo, atlas_region, atlas_monthly, reference


def fit_surface_components(
    albedo: pd.DataFrame,
    cfg: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Estimate signed Jan/Feb surface-state components on common support."""
    use = albedo[
        (albedo["albedo_type"].astype(str).str.upper() == "WSA")
        & (albedo["quality_mode"].astype(str) == "qa01")
        & albedo["month"].isin(MONTHS_TO_CORRECT)
    ].copy()

    rows = []
    model_rows = []

    for month in MONTHS_TO_CORRECT:
        g0 = use[use["month"] == month].copy()

        common_cols = [
            "viirs_anomaly_abs_nw",
            "snow_anomaly_2021_vs_flanks_pp",
            "baseline_raw_nw",
            "centroid_latitude",
            "cf_cvg_anomaly_2021_vs_flanks",
            "CNTR_CODE",
        ]
        snow_train = (
            g0[common_cols]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .copy()
        )
        snow_mu = float(
            np.log10(snow_train["baseline_raw_nw"].clip(lower=1e-6)).mean()
        )
        snow_train["log10_baseline_c"] = (
            np.log10(snow_train["baseline_raw_nw"].clip(lower=1e-6)) - snow_mu
        )
        snow_formula = (
            "viirs_anomaly_abs_nw ~ "
            "snow_anomaly_2021_vs_flanks_pp * log10_baseline_c + "
            "centroid_latitude + cf_cvg_anomaly_2021_vs_flanks + C(CNTR_CODE)"
        )
        snow_model = smf.ols(snow_formula, data=snow_train).fit(cov_type="HC3")

        comb_cols = common_cols + ["albedo_anomaly_pp"]
        comb_train = (
            g0[g0["eligible_albedo"].fillna(False).astype(bool)][comb_cols]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .copy()
        )
        comb_mu = float(
            np.log10(comb_train["baseline_raw_nw"].clip(lower=1e-6)).mean()
        )
        comb_train["log10_baseline_c"] = (
            np.log10(comb_train["baseline_raw_nw"].clip(lower=1e-6)) - comb_mu
        )
        comb_formula = (
            "viirs_anomaly_abs_nw ~ "
            "snow_anomaly_2021_vs_flanks_pp * log10_baseline_c + "
            "albedo_anomaly_pp * log10_baseline_c + "
            "centroid_latitude + cf_cvg_anomaly_2021_vs_flanks + C(CNTR_CODE)"
        )
        combined = smf.ols(comb_formula, data=comb_train).fit(cov_type="HC3")

        alb_formula = (
            "viirs_anomaly_abs_nw ~ "
            "albedo_anomaly_pp * log10_baseline_c + "
            "centroid_latitude + cf_cvg_anomaly_2021_vs_flanks + C(CNTR_CODE)"
        )
        albedo_only = smf.ols(alb_formula, data=comb_train).fit(cov_type="HC3")

        # Snow-only component for every row with the required fields.
        snow_pred = g0.copy()
        snow_pred["log10_baseline_c"] = (
            np.log10(pd.to_numeric(snow_pred["baseline_raw_nw"], errors="coerce").clip(lower=1e-6))
            - snow_mu
        )
        snow_ok = snow_pred[common_cols].replace(
            [np.inf, -np.inf], np.nan
        ).notna().all(axis=1)
        snow_actual = snow_pred.loc[snow_ok].copy()
        snow_zero = snow_actual.copy()
        snow_zero["snow_anomaly_2021_vs_flanks_pp"] = 0.0
        snow_component = pd.Series(np.nan, index=g0.index, dtype=float)
        snow_component.loc[snow_actual.index] = (
            snow_model.predict(snow_actual) - snow_model.predict(snow_zero)
        )

        # Combined and albedo-only components only where primary albedo is valid.
        comb_pred = g0.copy()
        comb_pred["log10_baseline_c"] = (
            np.log10(pd.to_numeric(comb_pred["baseline_raw_nw"], errors="coerce").clip(lower=1e-6))
            - comb_mu
        )
        comb_ok = (
            comb_pred["eligible_albedo"].fillna(False).astype(bool)
            & comb_pred[comb_cols].replace(
                [np.inf, -np.inf], np.nan
            ).notna().all(axis=1)
        )
        comb_actual = comb_pred.loc[comb_ok].copy()
        comb_zero = comb_actual.copy()
        comb_zero["snow_anomaly_2021_vs_flanks_pp"] = 0.0
        comb_zero["albedo_anomaly_pp"] = 0.0
        combined_component = pd.Series(np.nan, index=g0.index, dtype=float)
        combined_component.loc[comb_actual.index] = (
            combined.predict(comb_actual) - combined.predict(comb_zero)
        )

        alb_zero = comb_actual.copy()
        alb_zero["albedo_anomaly_pp"] = 0.0
        albedo_component = pd.Series(np.nan, index=g0.index, dtype=float)
        albedo_component.loc[comb_actual.index] = (
            albedo_only.predict(comb_actual) - albedo_only.predict(alb_zero)
        )

        out = g0[
            [
                "NUTS_ID", "CNTR_CODE", "density_band", "month",
                "viirs_anomaly_abs_nw", "snow_anomaly_2021_vs_flanks_pp",
                "albedo_anomaly_pp", "eligible_albedo",
            ]
        ].copy()
        out["month_name"] = MONTH_NAMES[month]
        out["surface_component_common_snow_nw"] = snow_component.to_numpy()
        out["surface_component_common_combined_nw"] = combined_component.to_numpy()
        out["surface_component_common_albedo_nw"] = albedo_component.to_numpy()
        out["surface_component_common_primary_nw"] = np.where(
            np.isfinite(out["surface_component_common_combined_nw"]),
            out["surface_component_common_combined_nw"],
            out["surface_component_common_snow_nw"],
        )
        out["surface_component_primary_method"] = np.where(
            np.isfinite(out["surface_component_common_combined_nw"]),
            "snow_plus_albedo",
            np.where(
                np.isfinite(out["surface_component_common_snow_nw"]),
                "snow_fallback",
                "no_surface_estimate",
            ),
        )
        rows.append(out)

        model_rows.append(
            {
                "month": month,
                "month_name": MONTH_NAMES[month],
                "n_snow_model": int(snow_model.nobs),
                "n_combined_model": int(combined.nobs),
                "r2_snow_model": float(snow_model.rsquared),
                "r2_combined_model": float(combined.rsquared),
                "snow_model_baseline_log10_mean": snow_mu,
                "combined_model_baseline_log10_mean": comb_mu,
            }
        )

    return pd.concat(rows, ignore_index=True), pd.DataFrame(model_rows)


def bridge_common_to_full_atlas(
    surface: pd.DataFrame,
    atlas_region: pd.DataFrame,
    cfg: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    full = atlas_region[
        (atlas_region["product"].astype(str).str.lower() == "vcmslcfg")
        & atlas_region["month"].isin(MONTHS_TO_CORRECT)
        & atlas_region["paired_valid"].fillna(False).astype(bool)
    ][["NUTS_ID", "month", "delta_nw", "flank_mean"]].copy()
    merged = surface.merge(
        full,
        on=["NUTS_ID", "month"],
        how="left",
        validate="one_to_one",
    )

    bridge_rows = []
    for month in MONTHS_TO_CORRECT:
        g = merged[merged["month"] == month][
            ["viirs_anomaly_abs_nw", "delta_nw"]
        ].replace([np.inf, -np.inf], np.nan).dropna()
        model = smf.ols(
            "delta_nw ~ viirs_anomaly_abs_nw",
            data=g,
        ).fit(cov_type="HC3")
        slope = float(model.params["viirs_anomaly_abs_nw"])
        bridge_rows.append(
            {
                "month": month,
                "month_name": MONTH_NAMES[month],
                "n": int(model.nobs),
                "slope_full_per_common": slope,
                "intercept_nw": float(model.params["Intercept"]),
                "r2": float(model.rsquared),
                "rmse_nw": float(np.sqrt(np.mean(model.resid.to_numpy(float) ** 2))),
            }
        )
        idx = merged["month"].eq(month)
        for stem in (
            "primary", "snow", "combined", "albedo",
        ):
            src = f"surface_component_common_{stem}_nw"
            if src in merged.columns:
                merged.loc[idx, f"surface_component_full_{stem}_nw"] = (
                    slope * merged.loc[idx, src]
                )

    bridge = pd.DataFrame(bridge_rows)
    min_r2 = float(cfg["min_common_to_full_r2"])
    if (bridge["r2"] < min_r2).any():
        raise RuntimeError(
            "Common-support to full-month bridge failed validation:\n"
            + bridge.to_string(index=False)
        )
    return merged, bridge


def load_tile_nuts(
    config: dict,
    atlas_monthly: pd.DataFrame,
) -> gpd.GeoDataFrame:
    ids = sorted(set(atlas_monthly["NUTS_ID"].astype(str)))
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release).copy()
    id_field = config["nuts"]["id_field"]
    if id_field != "NUTS_ID":
        nuts = nuts.rename(columns={id_field: "NUTS_ID"})
    nuts["NUTS_ID"] = nuts["NUTS_ID"].astype(str)
    nuts = nuts[nuts["NUTS_ID"].isin(ids)][["NUTS_ID", "geometry"]].copy()
    if len(nuts) != len(ids):
        raise RuntimeError(
            f"Only {len(nuts)} of {len(ids)} monthly-atlas NUTS3 geometries found."
        )
    return nuts


def extract_monthly_cf_weights(
    config: dict,
    atlas_monthly: pd.DataFrame,
    cache_dir: Path,
) -> pd.DataFrame:
    """Extract area-weighted mean(cf) and mean(rad*cf) for all 2020-2022 months."""
    pairs, missing = discover_inputs(configured_path(config, "viirs_raw"))
    need = {
        ("vcmslcfg", year, month)
        for year in YEARS_RECON
        for month in range(1, 13)
    }
    absent = sorted(need - set(pairs))
    if absent:
        raise FileNotFoundError(f"Missing VCMSLCFG monthly inputs: {absent}")

    nuts = load_tile_nuts(config, atlas_monthly)
    cache_dir.mkdir(parents=True, exist_ok=True)
    frames = []

    for year in YEARS_RECON:
        for month in range(1, 13):
            cache = cache_dir / f"vcmslcfg_{year}{month:02d}_cf_weight.parquet"
            pair = pairs[("vcmslcfg", year, month)]
            sig = "|".join(
                [
                    pair.radiance.name,
                    str(pair.radiance.stat().st_size),
                    pair.cloudfree.name,
                    str(pair.cloudfree.stat().st_size),
                ]
            )

            if cache.exists():
                d = pd.read_parquet(cache)
                if (
                    "source_signature" in d.columns
                    and d["source_signature"].dropna().astype(str).unique().tolist() == [sig]
                    and set(d["NUTS_ID"].astype(str)) == set(nuts["NUTS_ID"])
                ):
                    frames.append(d)
                    continue

            print(f"Extracting cf weights {year}-{month:02d}", flush=True)
            cf_src = RadianceCoverageSource(
                pair.radiance, pair.cloudfree, "cf", "cf"
            )
            rc_src = RadianceCoverageSource(
                pair.radiance, pair.cloudfree, "radcf", "radcf"
            )
            try:
                d = exact_extract(
                    [cf_src, rc_src],
                    nuts,
                    [
                        "mean(coverage_weight=area_spherical_m2)",
                    ],
                    include_cols=["NUTS_ID"],
                    output="pandas",
                    progress=True,
                )
            finally:
                cf_src.close()
                rc_src.close()

            d = d.rename(
                columns={
                    "cf_mean": "cf_area_mean",
                    "radcf_mean": "radcf_area_mean",
                }
            )
            required = {"NUTS_ID", "cf_area_mean", "radcf_area_mean"}
            if required - set(d.columns):
                raise RuntimeError(
                    f"Unexpected cf-weight extraction columns: {list(d.columns)}"
                )
            d["year"] = year
            d["month"] = month
            d["source_signature"] = sig
            d.to_parquet(cache, index=False)
            frames.append(d)

    out = pd.concat(frames, ignore_index=True)
    out["NUTS_ID"] = out["NUTS_ID"].astype(str)
    return out


def annual_reconstruction(weights: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (nuts_id, year), g in weights.groupby(["NUTS_ID", "year"], observed=True):
        den = float(g["cf_area_mean"].sum())
        num = float(g["radcf_area_mean"].sum())
        rows.append(
            {
                "NUTS_ID": nuts_id,
                "year": int(year),
                "sum_cf_area_mean": den,
                "sum_radcf_area_mean": num,
                "reconstructed_monthly_weighted_raw_nw": (
                    num / den if den > 0 else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def validate_annual_mapping(
    recon: pd.DataFrame,
    annual: pd.DataFrame,
    cfg: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    a = annual[
        annual["year"].isin([2020, 2022])
    ][
        [
            "NUTS_ID", "year",
            "radiance_mean_raw", "radiance_mean_corrected",
        ]
    ].copy()
    d = recon.merge(a, on=["NUTS_ID", "year"], how="inner", validate="one_to_one")
    d = d.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[
            "reconstructed_monthly_weighted_raw_nw",
            "radiance_mean_raw",
            "radiance_mean_corrected",
        ]
    )

    rows = []
    models = {}
    for outcome in ("radiance_mean_raw", "radiance_mean_corrected"):
        model = smf.ols(
            f"{outcome} ~ reconstructed_monthly_weighted_raw_nw + C(year)",
            data=d,
        ).fit(cov_type="HC3")
        models[outcome] = model
        rows.append(
            {
                "outcome": outcome,
                "n": int(model.nobs),
                "slope": float(
                    model.params["reconstructed_monthly_weighted_raw_nw"]
                ),
                "intercept": float(model.params["Intercept"]),
                "r2": float(model.rsquared),
                "rmse_nw": float(np.sqrt(np.mean(model.resid.to_numpy(float) ** 2))),
            }
        )

    metrics = pd.DataFrame(rows)

    wide = d.pivot(
        index="NUTS_ID",
        columns="year",
        values=[
            "reconstructed_monthly_weighted_raw_nw",
            "radiance_mean_corrected",
        ],
    ).dropna()
    recon_delta = (
        wide["reconstructed_monthly_weighted_raw_nw"][2022]
        - wide["reconstructed_monthly_weighted_raw_nw"][2020]
    )
    annual_delta = (
        wide["radiance_mean_corrected"][2022]
        - wide["radiance_mean_corrected"][2020]
    )
    corr = float(pearsonr(recon_delta, annual_delta).statistic)
    beta = float(
        models["radiance_mean_corrected"].params[
            "reconstructed_monthly_weighted_raw_nw"
        ]
    )
    predicted_delta = beta * recon_delta
    delta_rmse = float(
        np.sqrt(np.mean((annual_delta.to_numpy() - predicted_delta.to_numpy()) ** 2))
    )
    metrics["control_delta_pearson_r"] = corr
    metrics["control_delta_rmse_nw"] = delta_rmse

    if (
        metrics.loc[
            metrics["outcome"].eq("radiance_mean_corrected"), "r2"
        ].iloc[0]
        < float(cfg["min_annual_reconstruction_r2"])
    ):
        raise RuntimeError(
            "Monthly-to-annual corrected-radiance reconstruction is too weak:\n"
            + metrics.to_string(index=False)
        )
    if corr < float(cfg["min_control_delta_correlation"]):
        raise RuntimeError(
            "Monthly-to-annual change reconstruction is too weak:\n"
            + metrics.to_string(index=False)
        )

    return metrics, d


def build_annual_corrections(
    surface: pd.DataFrame,
    weights: pd.DataFrame,
    mapping: pd.DataFrame,
) -> pd.DataFrame:
    w21 = weights[weights["year"].eq(2021)].copy()
    totals = (
        w21.groupby("NUTS_ID", observed=True)["cf_area_mean"]
        .sum()
        .rename("annual_cf_area_mean")
        .reset_index()
    )
    w21 = w21.merge(totals, on="NUTS_ID", validate="many_to_one")
    w21["annual_month_weight"] = np.where(
        w21["annual_cf_area_mean"] > 0,
        w21["cf_area_mean"] / w21["annual_cf_area_mean"],
        np.nan,
    )

    s = surface.merge(
        w21[["NUTS_ID", "month", "annual_month_weight"]],
        on=["NUTS_ID", "month"],
        how="left",
        validate="one_to_one",
    )

    stems = ("primary", "snow", "albedo")
    for stem in stems:
        col = f"surface_component_full_{stem}_nw"
        s[f"weighted_{stem}_component_nw"] = (
            s[col] * s["annual_month_weight"]
        )

    agg_dict = {
        f"weighted_{stem}_component_nw": "sum"
        for stem in stems
    }
    agg_dict.update(
        {
            "month": "count",
        }
    )
    annual = (
        s.groupby("NUTS_ID", observed=True)
        .agg(agg_dict)
        .rename(columns={"month": "n_surface_month_rows"})
        .reset_index()
    )

    beta_raw = float(
        mapping.loc[mapping["outcome"].eq("radiance_mean_raw"), "slope"].iloc[0]
    )
    beta_corr = float(
        mapping.loc[
            mapping["outcome"].eq("radiance_mean_corrected"), "slope"
        ].iloc[0]
    )

    for stem in stems:
        recon_col = f"weighted_{stem}_component_nw"
        annual[f"annual_surface_correction_raw_{stem}_nw"] = beta_raw * annual[recon_col]
        annual[f"annual_surface_correction_corrected_{stem}_nw"] = (
            beta_corr * annual[recon_col]
        )

    method = (
        s.groupby("NUTS_ID", observed=True)["surface_component_primary_method"]
        .agg(lambda x: "+".join(sorted(set(x.dropna().astype(str)))))
        .rename("surface_component_methods")
        .reset_index()
    )
    annual = annual.merge(method, on="NUTS_ID", how="left")
    return annual


def apply_to_history(
    annual: pd.DataFrame,
    corrections: pd.DataFrame,
) -> pd.DataFrame:
    d = annual.merge(corrections, on="NUTS_ID", how="left", validate="many_to_one")

    scenarios = ("primary", "snow", "albedo")
    for scenario in scenarios:
        raw_out = f"radiance_mean_raw_surface_{scenario}"
        corr_out = f"radiance_mean_corrected_surface_{scenario}"
        d[raw_out] = d["radiance_mean_raw"]
        d[corr_out] = d["radiance_mean_corrected"]

        idx = d["year"].eq(2021)
        raw_c = d[f"annual_surface_correction_raw_{scenario}_nw"].fillna(0.0)
        corr_c = d[
            f"annual_surface_correction_corrected_{scenario}_nw"
        ].fillna(0.0)
        d.loc[idx, raw_out] = (
            d.loc[idx, "radiance_mean_raw"] - raw_c.loc[idx]
        ).clip(lower=0.0)
        d.loc[idx, corr_out] = (
            d.loc[idx, "radiance_mean_corrected"] - corr_c.loc[idx]
        ).clip(lower=1e-8)

    d["surface_correction_applied"] = (
        d["year"].eq(2021)
        & d["annual_surface_correction_corrected_primary_nw"].notna()
    )
    return d


def fit_corrected_trends(
    history: pd.DataFrame,
    reference: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    cfg = config["analysis"]["fixed_nuts2016_2014_2024"]
    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    pre_last = int(cfg["comparison_last_year"])

    tables = []
    trajectories = {}

    base = prepare_analysis_panel(
        history[
            [
                c for c in history.columns
                if not c.startswith("radiance_mean_corrected_surface_")
            ]
        ].copy(),
        reference,
        config,
    )
    t_pre, _ = fit_window(base, first_year, pre_last, config)
    t_orig, _ = fit_window(base, first_year, last_year, config)
    for label, table in (("2014_2020_reference", t_pre), ("original_2014_2024", t_orig)):
        x = table.copy()
        x["scenario"] = label
        tables.append(x)
    trajectories["original"] = observed_trajectories(base, baseline_year=first_year)

    for scenario in ("primary", "snow", "albedo"):
        x = history.copy()
        x["radiance_mean_corrected"] = x[
            f"radiance_mean_corrected_surface_{scenario}"
        ]
        panel = prepare_analysis_panel(x, reference, config)
        trends, _ = fit_window(panel, first_year, last_year, config)
        trends["scenario"] = f"corrected_{scenario}"
        tables.append(trends)
        trajectories[scenario] = observed_trajectories(
            panel, baseline_year=first_year
        )

    return pd.concat(tables, ignore_index=True), trajectories


def correction_summary_by_density(
    history: pd.DataFrame,
    reference: pd.DataFrame,
) -> pd.DataFrame:
    y = history[history["year"].eq(2021)].merge(
        reference[["NUTS_ID", "density_band"]],
        on="NUTS_ID",
        how="left",
        validate="one_to_one",
    )
    return (
        y.groupby("density_band", observed=True)
        .agg(
            n_nuts3=("NUTS_ID", "nunique"),
            mean_original_corrected_nw=("radiance_mean_corrected", "mean"),
            median_original_corrected_nw=("radiance_mean_corrected", "median"),
            mean_primary_correction_nw=(
                "annual_surface_correction_corrected_primary_nw", "mean"
            ),
            median_primary_correction_nw=(
                "annual_surface_correction_corrected_primary_nw", "median"
            ),
            mean_primary_corrected_2021_nw=(
                "radiance_mean_corrected_surface_primary", "mean"
            ),
            median_primary_corrected_2021_nw=(
                "radiance_mean_corrected_surface_primary", "median"
            ),
        )
        .reset_index()
    )


def plot_trajectories(
    trajectories: dict[str, pd.DataFrame],
    config: dict,
    output: Path,
):
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0), sharex=True)
    axes = axes.ravel()
    for i, band in enumerate(labels):
        ax = axes[i]
        for key, label, marker in (
            ("original", "Original", "o"),
            ("primary", "2021 surface-corrected", "s"),
        ):
            g = trajectories[key]
            g = g[g["density_band"].astype(str).eq(str(band))].sort_values("year")
            ax.plot(
                g["year"],
                g["observed_index_2014_100"],
                marker=marker,
                linewidth=1.4,
                label=label,
            )
        ax.axhline(100, linewidth=0.7)
        ax.set_title(str(band))
        ax.grid(alpha=0.2)
    axes[5].axis("off")
    handles, labs = axes[0].get_legend_handles_labels()
    axes[5].legend(handles, labs, loc="center", frameon=False)
    fig.supylabel("Calibrated VIIRS radiance index (2014 = 100)")
    fig.supxlabel("Year")
    fig.suptitle("Fixed NUTS2016 trajectories before and after 2021 surface correction")
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.95))
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_2021_surface_correction"]
    (
        annual, albedo, atlas_region, atlas_monthly, reference
    ) = load_inputs(config, cfg)

    surface_common, surface_models = fit_surface_components(albedo, cfg)
    surface_full, bridge = bridge_common_to_full_atlas(
        surface_common, atlas_region, cfg
    )

    processed = configured_path(config, "processed")
    out_dir = processed / "fixed_nuts2016_2021_surface_correction"
    cache_dir = out_dir / "monthly_cf_weight_cache"
    figure_dir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_2021_surface_correction"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    weights = extract_monthly_cf_weights(
        config, atlas_monthly, cache_dir
    )
    recon = annual_reconstruction(weights)
    mapping, mapping_rows = validate_annual_mapping(recon, annual, cfg)
    corrections = build_annual_corrections(surface_full, weights, mapping)
    history = apply_to_history(annual, corrections)
    trends, trajectories = fit_corrected_trends(
        history, reference, config
    )
    density = correction_summary_by_density(history, reference)

    outputs = {
        "surface_monthly": out_dir / "surface_state_monthly_components.csv",
        "surface_models": out_dir / "surface_state_model_summary.csv",
        "bridge": out_dir / "common_to_full_month_bridge.csv",
        "weights": out_dir / "monthly_cf_weights_2020_2022.csv",
        "reconstruction": out_dir / "monthly_weighted_annual_reconstruction.csv",
        "mapping": out_dir / "annual_reconstruction_validation.csv",
        "mapping_rows": out_dir / "annual_reconstruction_control_rows.csv",
        "corrections": out_dir / "annual_2021_surface_corrections.csv",
        "history": out_dir / "fixed_nuts2016_viirs_2014_2024_surface_corrected.parquet",
        "history_csv": out_dir / "fixed_nuts2016_viirs_2014_2024_surface_corrected.csv",
        "trends": out_dir / "density_trends_original_vs_surface_corrected.csv",
        "density": out_dir / "2021_correction_by_density.csv",
        "summary": out_dir / "surface_correction_summary.txt",
    }

    surface_full.to_csv(outputs["surface_monthly"], index=False)
    surface_models.to_csv(outputs["surface_models"], index=False)
    bridge.to_csv(outputs["bridge"], index=False)
    weights.to_csv(outputs["weights"], index=False)
    recon.to_csv(outputs["reconstruction"], index=False)
    mapping.to_csv(outputs["mapping"], index=False)
    mapping_rows.to_csv(outputs["mapping_rows"], index=False)
    corrections.to_csv(outputs["corrections"], index=False)
    history.to_parquet(outputs["history"], index=False)
    history.to_csv(outputs["history_csv"], index=False)
    trends.to_csv(outputs["trends"], index=False)
    density.to_csv(outputs["density"], index=False)

    fig = figure_dir / "trajectories_original_vs_surface_corrected.pdf"
    plot_trajectories(trajectories, config, fig)

    primary_corr = corrections[
        "annual_surface_correction_corrected_primary_nw"
    ].dropna()
    primary_methods = (
        surface_full["surface_component_primary_method"]
        .value_counts(dropna=False)
        .to_dict()
    )

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 2021 SURFACE-STATE CORRECTION\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            "Primary correction: Jan/Feb 2021 snow + WSA visible-albedo "
            "counterfactual (surface anomalies set to zero), with snow-only "
            "fallback when primary albedo is unavailable. March and all other "
            "months are left unchanged.\n"
        )
        h.write(
            "Annual propagation uses actual monthly VCMSLCFG cf_cvg observation "
            "weights. It does not use an unweighted mean of monthly radiance.\n"
        )
        h.write(
            "Only a copy of the fixed 2014-2024 panel is changed; original "
            "production outputs remain untouched.\n\n"
        )

        h.write("COMMON-SUPPORT -> FULL-MONTH BRIDGE\n")
        h.write(bridge.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nMONTHLY -> ANNUAL CONTROL-YEAR VALIDATION\n")
        h.write(mapping.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nSURFACE MODEL COVERAGE\n")
        h.write(repr(primary_methods) + "\n")
        h.write(
            f"Regions with annual primary correction: {len(primary_corr)}; "
            f"median={primary_corr.median():+.6g} nW; "
            f"mean={primary_corr.mean():+.6g} nW.\n"
        )
        h.write("\n2021 CORRECTION BY FIXED DENSITY BAND\n")
        h.write(density.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nHISTORICAL TREND COMPARISON\n")
        h.write(trends.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write(
            "\n\nInterpretation: this stage is the candidate 2021 correction "
            "used to decide the final historical-response model. It should be "
            "accepted only if the control-year annual reconstruction passes the "
            "pre-specified validation thresholds and the corrected 2021 "
            "trajectory removes the surface-state discontinuity without creating "
            "a new discontinuity in neighbouring years.\n"
        )

    print("\n2021 SURFACE-STATE CORRECTION COMPLETE")
    print("\nANNUAL RECONSTRUCTION VALIDATION")
    print(mapping.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nCORRECTION BY DENSITY")
    print(density.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nTREND COMPARISON")
    print(trends.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for p in [*outputs.values(), fig]:
        print(p)


if __name__ == "__main__":
    main()
