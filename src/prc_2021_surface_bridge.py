#!/usr/bin/env python3
"""Sensitivity: translate the existing 2021 NUTS3 snow/albedo correction into PRC.

This is a regional COUNTERFACTUAL test, not a calibrated pixel snow map.
For each fixed NUTS2016 region, scale its 2021 calibrated VIIRS pixels so that
their corrected-radiance mean agrees with the independently modelled 2021
NUTS3 counterfactual radiance mean, then average PRC(pixel) over the polygon.
The nonlinear transformation is therefore evaluated at pixel level.

Do not interpret within-NUTS3 scaling as a spatially resolved snow correction.
No VIIRS GeoTIFFs or existing historical PRC outputs are overwritten.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource

from apply_radiance_calibration import CalibrationRasterSource
from extract_nuts3_radiance import discover_rasters, nuts_path_for_release
from historical_prc_pixel_nuts3 import prc_of_radiance, load_inputs
from project_config import PROJECT_ROOT, load_config


class ScaledPRCRaster(RasterioRasterSource):
    def __init__(self, dataset, offset, params, factor, name):
        super().__init__(dataset, 1, name=name)
        self.offset = float(offset)
        self.params = params
        self.factor = float(factor)

    def read_window(self, x0, y0, nx, ny):
        arr = super().read_window(x0, y0, nx, ny)
        masked = np.ma.asarray(arr)
        data = np.asarray(masked.data, dtype=np.float32)
        invalid = (np.ma.getmaskarray(masked) |
                   ~np.isfinite(data) | (data < 0))
        calibrated = np.maximum(data - self.offset, 0.0)
        # Params from the same pinned transfer config as the historical atlas.
        yf, y0, k = self.params
        response = yf + (y0 - yf) * np.exp(-k * self.factor * calibrated)
        return np.ma.array(response.astype(np.float32), mask=invalid)


def locate(path):
    p = Path(path)
    return p if p.is_absolute() else PROJECT_ROOT / p


def reference_table(path):
    if not path.is_file():
        raise FileNotFoundError(f"Run prior historical surface stage first: {path}")
    d = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    if d.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate region-year rows in original correction history")
    d = d.loc[d.year.eq(2021)].copy()
    required = [
        "NUTS_ID", "radiance_mean_corrected",
        "radiance_mean_corrected_surface_logratio_primary",
        "radiance_mean_corrected_surface_logratio_snow",
    ]
    missing = set(required) - set(d.columns)
    if missing:
        raise ValueError(f"Historical correction lacks columns: {sorted(missing)}")
    return d[required]


def get_inputs(args):
    config = load_config(args.config)
    transfer, params, offset, from_year = load_inputs(
        config, locate(args.transfer_config))
    if from_year > 2021:
        raise ValueError("Unexpected 2021 calibration start year")
    raster = discover_rasters(config, [2021])[2021]
    nuts = gpd.read_file(nuts_path_for_release(config, 2016))
    if nuts.crs is None:
        raise ValueError("NUTS2016 CRS is undefined")
    nuts = nuts.to_crs("EPSG:4326")
    if nuts.NUTS_ID.duplicated().any():
        raise ValueError("Duplicate NUTS ids")
    frame = pd.read_csv(locate(args.prc_history))
    if frame.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate PRC region-years")
    needed = frame.loc[frame.year.isin([2020, 2021, 2022]),
                       ["NUTS_ID", "year", "prc_mean", "valid_area_km2"]]
    w = needed.pivot(index="NUTS_ID", columns="year", values="prc_mean")
    for y in (2020, 2021, 2022):
        if y not in w:
            raise ValueError(f"Historical PRC series missing {y}")
    inputs = nuts.merge(
        reference_table(locate(args.surface_history)), on="NUTS_ID",
        how="left", validate="one_to_one")
    inputs = inputs.merge(
        w.rename(columns={2020:"prc_2020", 2021:"prc_2021",
                          2022:"prc_2022"}).reset_index(),
        on="NUTS_ID", how="left", validate="one_to_one")
    a21 = needed.loc[needed.year.eq(2021),
                     ["NUTS_ID", "valid_area_km2"]]
    inputs = inputs.merge(a21, on="NUTS_ID", how="left", validate="one_to_one")
    return transfer, params, offset, raster, inputs


def factors_from_model(d, factor_min, factor_max):
    obs = pd.to_numeric(d.radiance_mean_corrected, errors="coerce")
    denom = obs.replace(0.0, np.nan)
    for scenario, field in [
        ("primary", "radiance_mean_corrected_surface_logratio_primary"),
        ("snow_only", "radiance_mean_corrected_surface_logratio_snow"),
    ]:
        adjusted = pd.to_numeric(d[field], errors="coerce")
        factor = adjusted / denom
        valid = (np.isfinite(factor) & (factor >= factor_min) &
                 (factor <= factor_max) & np.isfinite(adjusted) &
                 (adjusted >= 0))
        d[f"factor_{scenario}"] = factor
        d[f"eligible_{scenario}"] = valid
    return d


def numeric_value(frame, expected):
    if expected in frame.columns:
        return float(frame.iloc[0][expected])
    # exactextract may omit named RasterSource prefixes for one source only.
    raise ValueError(f"Missing {expected} in columns {list(frame.columns)}")


def compute_one(dataset, geo, offset, params, q1, q2):
    sources = [
        CalibrationRasterSource(
            dataset, offset=offset, mode="corrected", name="radiance_anchor"),
        ScaledPRCRaster(dataset, offset, params, 1.0, "prc_reference"),
        ScaledPRCRaster(dataset, offset, params, q1, "prc_primary"),
        ScaledPRCRaster(dataset, offset, params, q2, "prc_snow_only"),
    ]
    frame = exact_extract(
        sources, geo, ["mean(coverage_weight=area_spherical_m2)"],
        include_cols=["NUTS_ID"], output="pandas", progress=False)
    if len(frame) != 1:
        raise ValueError("Expected exactly one NUTS polygon extraction")
    return {n: numeric_value(frame, n + "_mean") for n in
            ["radiance_anchor", "prc_reference", "prc_primary",
             "prc_snow_only"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config/pipeline.yaml")
    ap.add_argument("--transfer-config", default="config/prc_transfer.yaml")
    ap.add_argument("--surface-history",
        default="data/processed/fixed_nuts2016_2021_surface_correction/"
                "fixed_nuts2016_viirs_2014_2024_surface_corrected.parquet")
    ap.add_argument("--prc-history",
        default="data/processed/historical_prc_2014_2024/"
                "prc_nuts3_annual_cumulative.csv")
    ap.add_argument("--output-dir",
        default="data/processed/historical_prc_2021_surface_bridge")
    ap.add_argument("--factor-min", type=float, default=.25)
    ap.add_argument("--factor-max", type=float, default=4.0)
    ap.add_argument("--limit", type=int, default=0,
                    help="For initial smoke test only; first NUTS regions.")
    args = ap.parse_args()
    if not (0 < args.factor_min < 1 < args.factor_max):
        raise ValueError("Factor safeguards must encompass 1")
    _, params, offset, raster, inputs = get_inputs(args)
    inputs = factors_from_model(inputs, args.factor_min, args.factor_max)
    if args.limit:
        inputs = inputs.head(args.limit).copy()
    path = locate(args.output_dir)
    path.mkdir(parents=True, exist_ok=True)
    rows = []
    print(f"Regional test for EOG 2021: {raster}", flush=True)
    print("Existing post-2017 radiance offset:", offset, "(not reapplied twice)")
    print(f"NUTS regions to evaluate: {len(inputs)}", flush=True)
    with rasterio.open(raster) as ds:
        if ds.crs is None or ds.crs.to_epsg() != 4326:
            raise ValueError("VIIRS must be EPSG:4326")
        for i, (_, r) in enumerate(inputs.iterrows(), 1):
            q1 = (float(r.factor_primary)
                  if r.eligible_primary else 1.0)
            q2 = (float(r.factor_snow_only)
                  if r.eligible_snow_only else 1.0)
            geo = inputs.loc[[r.name], ["NUTS_ID", "geometry"]]
            sample = compute_one(ds, geo, offset, params, q1, q2)
            original = float(r.prc_2021)
            orig_rad = float(r.radiance_mean_corrected)
            if np.isfinite(original) and abs(sample["prc_reference"] - original) > 2e-5:
                raise RuntimeError(
                    f"{r.NUTS_ID}: original 2021 PRC polygon mean does not "
                    f"reproduce the historical PRC atlas: extracted "
                    f"{sample['prc_reference']:.7f} vs saved {original:.7f}")
            if np.isfinite(orig_rad) and abs(
                sample["radiance_anchor"] - orig_rad) > max(.01, .01 * abs(orig_rad)):
                raise RuntimeError(
                    f"{r.NUTS_ID}: pixel-calibrated annual radiance disagrees "
                    f"with existing surface-correction input: extracted "
                    f"{sample['radiance_anchor']:.7f} vs saved {orig_rad:.7f}")
            rows.append({
                "NUTS_ID": r.NUTS_ID,
                "CNTR_CODE": r.get("CNTR_CODE", None),
                "NUTS_NAME": r.get("NUTS_NAME", None),
                "valid_area_km2": r.valid_area_km2,
                "prc_2020": r.prc_2020,
                "prc_2021_original": original,
                "prc_2022": r.prc_2022,
                "radiance_2021_observed_corrected": orig_rad,
                "radiance_2021_surface_model_primary":
                    r.radiance_mean_corrected_surface_logratio_primary,
                "factor_primary": r.factor_primary,
                "factor_snow_only": r.factor_snow_only,
                "eligible_primary": bool(r.eligible_primary),
                "eligible_snow_only": bool(r.eligible_snow_only),
                "prc_2021_primary": sample["prc_primary"] if
                    r.eligible_primary else np.nan,
                "prc_2021_snow_only": sample["prc_snow_only"] if
                    r.eligible_snow_only else np.nan,
            })
            if i % 100 == 0 or i == len(inputs):
                print(f"  processed {i}/{len(inputs)} regions", flush=True)
    d = pd.DataFrame(rows)
    d["prc_2021_flanks_mean"] = .5 * (d.prc_2020 + d.prc_2022)
    for mode in ("original", "primary", "snow_only"):
        prc_col = f"prc_2021_{mode}"
        d[f"prc_2021_anomaly_{mode}"] = d[prc_col] - d.prc_2021_flanks_mean
    d["prc_correction_primary"] = d.prc_2021_primary - d.prc_2021_original
    d["prc_correction_snow_only"] = d.prc_2021_snow_only - d.prc_2021_original
    out = path / ("prc_2021_region_bridge_pilot.csv"
                  if args.limit else "prc_2021_region_bridge.csv")
    d.to_csv(out, index=False)
    print("Saved", out, flush=True)
    if args.limit:
        return
    valid = d.loc[d.eligible_primary &
                  np.isfinite(d.prc_2021_flanks_mean) &
                  np.isfinite(d.valid_area_km2) &
                  (d.valid_area_km2 > 0)].copy()
    if len(valid) < 1000:
        raise RuntimeError(
            f"Only {len(valid)} eligible NUTS3; review model coverage, "
            "factor limits and original correction input")
    def wmean(x, col):
        return float(np.average(x[col], weights=x.valid_area_km2))
    a0 = wmean(valid, "prc_2021_anomaly_original")
    a1 = wmean(valid, "prc_2021_anomaly_primary")
    a2 = wmean(valid, "prc_2021_anomaly_snow_only")
    ae0 = float(np.median(np.abs(valid.prc_2021_anomaly_original)))
    ae1 = float(np.median(np.abs(valid.prc_2021_anomaly_primary)))
    comparisons = pd.DataFrame([
        {"scenario":"original", "n_regions":len(valid),
         "area_weighted_2021_prc":wmean(valid, "prc_2021_original"),
         "area_weighted_2021_vs_flanks":a0,
         "median_region_abs_deviation_from_flanks":ae0},
        {"scenario":"NUTS3_surface_primary", "n_regions":len(valid),
         "area_weighted_2021_prc":wmean(valid, "prc_2021_primary"),
         "area_weighted_2021_vs_flanks":a1,
         "median_region_abs_deviation_from_flanks":ae1},
        {"scenario":"NUTS3_surface_snow_only", "n_regions":len(valid),
         "area_weighted_2021_prc":wmean(valid, "prc_2021_snow_only"),
         "area_weighted_2021_vs_flanks":a2,
         "median_region_abs_deviation_from_flanks":float(
             np.median(np.abs(valid.prc_2021_anomaly_snow_only)))},
    ])
    comparisons.to_csv(path / "prc_2021_bridge_comparison.csv", index=False)
    c = d.groupby("CNTR_CODE", dropna=False).agg(
        n_regions=("NUTS_ID", "size"),
        n_eligible=("eligible_primary", "sum"),
        median_prc_correction=("prc_correction_primary", "median"),
        median_original_2021_anomaly=("prc_2021_anomaly_original", "median"),
        median_corrected_2021_anomaly=("prc_2021_anomaly_primary", "median")
    ).reset_index()
    c.to_csv(path / "prc_2021_bridge_by_country.csv", index=False)

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    axs[0].scatter(valid.prc_2021_anomaly_original,
                   valid.prc_2021_anomaly_primary, s=7, alpha=.35)
    # axline((0, 0), (1, 1)) expands axes to 0..1, hiding the ~0.001
    # PRC-scale cloud; draw a bounded identity segment instead.
    both = np.concatenate([
        valid.prc_2021_anomaly_original.to_numpy(float),
        valid.prc_2021_anomaly_primary.to_numpy(float),
    ])
    finite = both[np.isfinite(both)]
    lim = (float(finite.min() - .001), float(finite.max() + .001))
    axs[0].plot(lim, lim, color="gray", linestyle="--")
    axs[0].axhline(0, color="gray", linewidth=.6)
    axs[0].axvline(0, color="gray", linewidth=.6)
    axs[0].set(xlim=lim, ylim=lim,
               xlabel="Original 2021 PRC anomaly vs 2020/22 flanks",
               ylabel="2021 anomaly after regional correction")
    axs[1].scatter(valid.radiance_2021_observed_corrected,
                   valid.prc_correction_primary, s=7, alpha=.4)
    axs[1].set(xscale="symlog", xlabel="Corrected 2021 mean radiance (nW cm⁻² sr⁻¹)",
               ylabel="Change to NUTS mean PRC from regional correction")
    for ax in axs:ax.grid(alpha=.2)
    fig.savefig(path / "prc_2021_bridge_diagnostics.pdf")
    plt.close(fig)
    summary=[
        "2021 regional surface correction applied to calibrated *pixel* radiance",
        "before NUTS3 mean PRC extraction. SENSITIVITY, not adopted pixel snow map.",
        "NUTS 2021 correction from prior regional snow/albedo model; no refitting.",
        "A single factor per region is assumed; no within-region snow localization.",
        f"Total regions={len(d)}, primary eligible={int(d.eligible_primary.sum())}, "
        f"snow-only eligible={int(d.eligible_snow_only.sum())}",
        "2021 PRC anomaly = PRC(2021) - 0.5*(PRC(2020)+PRC(2022)).",
        "2022 version changes can invalidate the flanking-year diagnostic.",
        comparisons.to_string(index=False),
        "Do not choose a correction purely to minimize flanking-year discrepancy.",
        "Use monthly snow/albedo pixel diagnostics, independent controls, and",
        "cross-product tests before elevating this sensitivity to main analysis.",
    ]
    (path / "prc_2021_bridge_summary.txt").write_text(
        "\n".join(summary)+"\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
