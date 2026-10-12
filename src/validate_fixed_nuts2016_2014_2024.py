#!/usr/bin/env python3
"""Fixed-geography 2014-2024 validation of density-dependent VIIRS change.

This is the decisive historical validation for the future-radiance response.

Why this stage exists
---------------------
The year-specific production pipeline changes NUTS release after 2020. The
2013 endpoint also has strong leverage on the historical trend. This stage
removes both issues by:

1. excluding 2013;
2. re-extracting every Annual VNL V2.2 raster from 2014 through 2024 onto the
   SAME NUTS2016 polygons;
3. applying the production 0.235 nW cm-2 sr-1 pixel correction from 2017
   onward before polygon aggregation;
4. assigning each NUTS2016 region to a fixed population-density class using
   its median 2014-2020 Eurostat population density; and
5. comparing country-adjusted trends over 2014-2020 versus 2014-2024.

The response is calibrated VIIRS DNB radiance at NUTS3 scale. This analysis
does not claim that city-centre pixels themselves change at the same rate as
the NUTS3 mean.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import linregress
from shapely.geometry import box

from analyze_density_temporal_change import fit_band_model
from apply_radiance_calibration import CalibrationRasterSource
from extract_nuts3_radiance import (
    discover_rasters,
    nuts_path_for_release,
    validate_raster,
)
from production_validation import validate_production_merged
from project_config import configured_path, load_config


def load_fixed_nuts(config: dict, release: int) -> gpd.GeoDataFrame:
    path = nuts_path_for_release(config, release)
    if not path.exists():
        raise FileNotFoundError(
            f"NUTS {release} Level-3 file not found: {path}"
        )
    nuts = gpd.read_file(path)
    if nuts.crs is None:
        raise ValueError(f"NUTS file has no CRS: {path}")
    target_crs = config["nuts"]["crs"]
    if nuts.crs.to_string() != target_crs:
        nuts = nuts.to_crs(target_crs)

    id_field = config["nuts"]["id_field"]
    if id_field not in nuts.columns:
        raise ValueError(f"Missing NUTS identifier field {id_field!r}")
    if nuts[id_field].duplicated().any():
        raise ValueError("Duplicate NUTS identifiers in fixed geometry.")
    return nuts


def extract_fixed_year(
    raster_path: Path,
    year: int,
    nuts: gpd.GeoDataFrame,
    release: int,
    config: dict,
) -> pd.DataFrame:
    """Extract raw and production-corrected radiance on fixed polygons."""
    id_field = config["nuts"]["id_field"]
    cal = config["calibration"]
    apply_from = int(cal["apply_from_year"])
    offset = float(cal["pixel_offset_nW_cm2_sr"])
    expected_crs = config["nuts"]["crs"]

    validate_raster(raster_path, expected_crs)

    include_cols = [
        c for c in [
            id_field,
            "CNTR_CODE",
            "NUTS_NAME",
            "NAME_LATN",
            "EU_STAT",
            "EFTA_STAT",
            "CC_STAT",
        ]
        if c in nuts.columns
    ]

    effective_offset = offset if year >= apply_from else 0.0
    print(
        f"Fixed NUTS {release}: extracting {year} "
        f"({raster_path.name}); pixel offset={effective_offset:.3f}",
        flush=True,
    )

    with rasterio.open(raster_path) as ds:
        sources = [
            RasterioRasterSource(ds, 1, name="raw"),
            CalibrationRasterSource(
                ds,
                offset=effective_offset,
                mode="corrected",
                name="corrected",
            ),
        ]
        result = exact_extract(
            sources,
            nuts,
            ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=include_cols,
            output="pandas",
            progress=True,
        )

    expected = {id_field, "raw_mean", "corrected_mean"}
    missing = expected - set(result.columns)
    if missing:
        raise RuntimeError(
            "Unexpected exactextract output. "
            f"Missing {sorted(missing)}; found {list(result.columns)}"
        )

    result = result.rename(
        columns={
            "raw_mean": "radiance_mean_raw",
            "corrected_mean": "radiance_mean_corrected",
        }
    )
    result.insert(0, "year", int(year))
    result["nuts_release"] = int(release)
    result["calibration_correction_applied"] = year >= apply_from
    result["calibration_pixel_offset"] = offset
    result["calibration_effective_pixel_offset"] = effective_offset
    result["calibration_method"] = cal["method"]
    result["calibration_version"] = cal["version"]
    result["radiance_units"] = config["viirs"]["units"]
    result["viirs_file"] = raster_path.name
    result["zonal_method"] = config["analysis"]["area_weighting"]["method"]

    tolerance = 2e-6
    raw = result["radiance_mean_raw"].to_numpy(dtype=float)
    corrected = result["radiance_mean_corrected"].to_numpy(dtype=float)
    finite = np.isfinite(raw) & np.isfinite(corrected)
    if year < apply_from:
        limit = tolerance * np.maximum(1.0, np.abs(raw))
        if (finite & (np.abs(corrected - raw) > limit)).any():
            raise RuntimeError(
                f"Pre-calibration fixed extraction differs from raw in {year}."
            )
    else:
        if (finite & (corrected > raw + tolerance)).any():
            raise RuntimeError(
                f"Pixel calibration increased at least one mean in {year}."
            )
        if (finite & (corrected < -tolerance)).any():
            raise RuntimeError(
                f"Negative corrected regional mean detected in {year}."
            )

    return result


def validate_against_production(
    fixed: pd.DataFrame,
    config: dict,
    first_year: int,
    comparison_last_year: int,
    release: int,
) -> pd.DataFrame:
    """Check that 2014-2020 fixed extraction reproduces production exactly."""
    processed = configured_path(config, "processed")
    production = pd.read_parquet(
        processed / config["outputs"]["merged"]
    ).copy()
    validate_production_merged(production, config)

    prod = production[
        production["year"].between(first_year, comparison_last_year)
        & (production["nuts_release"].astype(int) == release)
    ][
        ["year", "NUTS_ID", "radiance_mean_corrected"]
    ].rename(
        columns={
            "radiance_mean_corrected": "production_radiance_mean_corrected"
        }
    )

    check = fixed[
        fixed["year"].between(first_year, comparison_last_year)
    ][
        ["year", "NUTS_ID", "radiance_mean_corrected"]
    ].merge(
        prod,
        on=["year", "NUTS_ID"],
        how="inner",
        validate="one_to_one",
    )

    check["difference"] = (
        check["radiance_mean_corrected"]
        - check["production_radiance_mean_corrected"]
    )
    check["abs_difference"] = check["difference"].abs()

    if check.empty:
        raise RuntimeError(
            "No overlap available for fixed-extraction production validation."
        )

    max_abs = float(check["abs_difference"].max())
    raw_scale = np.maximum(
        1.0,
        np.abs(check["production_radiance_mean_corrected"].to_numpy(float)),
    )
    tolerance = 2e-6 * raw_scale
    if (
        check["abs_difference"].to_numpy(float) > tolerance
    ).any():
        bad = check.loc[
            check["abs_difference"].to_numpy(float) > tolerance
        ].head(10)
        raise RuntimeError(
            "Fixed NUTS2016 extraction does not reproduce the existing "
            "production 2014-2020 extraction. Examples:\n"
            + bad.to_string(index=False)
        )

    return check


def build_population_reference(
    config: dict,
    release: int,
    first_year: int,
    last_year: int,
) -> pd.DataFrame:
    """Fixed population-density classes from 2014-2020 NUTS2016 data."""
    processed = configured_path(config, "processed")
    production = pd.read_parquet(
        processed / config["outputs"]["merged"]
    ).copy()
    validate_production_merged(production, config)

    pop = production[
        production["year"].between(first_year, last_year)
        & (production["nuts_release"].astype(int) == release)
    ][
        ["year", "NUTS_ID", "CNTR_CODE", "population_density"]
    ].copy()
    pop["population_density"] = pd.to_numeric(
        pop["population_density"], errors="coerce"
    )
    pop = pop[
        np.isfinite(pop["population_density"])
        & (pop["population_density"] > 0)
    ].copy()

    ref = (
        pop.groupby(["NUTS_ID", "CNTR_CODE"], as_index=False)
        .agg(
            reference_population_density=(
                "population_density", "median"
            ),
            n_population_years=("year", "nunique"),
        )
    )

    cfg = config["analysis"]["density_temporal_change"]
    edges = [-np.inf] + list(cfg["density_band_breaks"]) + [np.inf]
    ref["density_band"] = pd.cut(
        ref["reference_population_density"],
        bins=edges,
        labels=cfg["density_band_labels"],
        right=False,
        ordered=True,
    )
    return ref


def prepare_analysis_panel(
    fixed: pd.DataFrame,
    reference: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    d = fixed.merge(
        reference,
        on=["NUTS_ID", "CNTR_CODE"],
        how="inner",
        validate="many_to_one",
    )
    d["radiance_mean_corrected"] = pd.to_numeric(
        d["radiance_mean_corrected"], errors="coerce"
    )
    d = d[
        np.isfinite(d["radiance_mean_corrected"])
        & (d["radiance_mean_corrected"] > 0)
        & d["density_band"].notna()
    ].copy()

    cfg = config["analysis"]["fixed_nuts2016_2014_2024"]
    expected_years = int(cfg["last_year"]) - int(cfg["first_year"]) + 1
    counts = d.groupby("NUTS_ID")["year"].nunique()
    if bool(cfg["require_complete_panel"]):
        keep = counts[counts == expected_years].index
    else:
        keep = counts[
            counts >= int(cfg["minimum_panel_years"])
        ].index
    d = d[d["NUTS_ID"].isin(keep)].copy()

    country_counts = d.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    keep_countries = country_counts[
        country_counts >= int(cfg["min_country_regions"])
    ].index
    d = d[d["CNTR_CODE"].isin(keep_countries)].copy()

    d["log_rad"] = np.log10(d["radiance_mean_corrected"])
    d["country"] = d["CNTR_CODE"].astype("category")
    return d


def fit_window(panel, first_year, last_year, config):
    d = panel[panel["year"].between(first_year, last_year)].copy()
    d["year_c"] = d["year"] - ((first_year + last_year) / 2.0)
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    trends, info = fit_band_model(d, labels)
    trends["first_year"] = first_year
    trends["last_year"] = last_year
    trends["interval"] = f"{first_year}-{last_year}"
    trends["pct_change_per_year_ci_low"] = 100.0 * (
        10.0 ** (
            trends["log10_radiance_trend_per_year"]
            - 1.96 * trends["trend_se"]
        ) - 1.0
    )
    trends["pct_change_per_year_ci_high"] = 100.0 * (
        10.0 ** (
            trends["log10_radiance_trend_per_year"]
            + 1.96 * trends["trend_se"]
        ) - 1.0
    )
    return trends, info


def observed_trajectories(panel, baseline_year=2014):
    baseline = panel[panel["year"] == baseline_year][
        ["NUTS_ID", "log_rad"]
    ].rename(columns={"log_rad": "baseline_log_rad"})
    d = panel.merge(
        baseline, on="NUTS_ID", how="inner", validate="many_to_one"
    )
    d["delta_log10"] = d["log_rad"] - d["baseline_log_rad"]

    rows = []
    for (band, year), g in d.groupby(
        ["density_band", "year"], observed=True
    ):
        vals = g["delta_log10"].to_numpy(float)
        mean_log = float(vals.mean())
        se = (
            float(vals.std(ddof=1) / np.sqrt(len(vals)))
            if len(vals) > 1 else np.nan
        )
        rows.append({
            "density_band": str(band),
            "year": int(year),
            "n_nuts3": int(g["NUTS_ID"].nunique()),
            "mean_delta_log10": mean_log,
            "se_delta_log10": se,
            "observed_index_2014_100": 100.0 * 10.0 ** mean_log,
            "observed_index_ci_low": (
                100.0 * 10.0 ** (mean_log - 1.96 * se)
                if np.isfinite(se) else np.nan
            ),
            "observed_index_ci_high": (
                100.0 * 10.0 ** (mean_log + 1.96 * se)
                if np.isfinite(se) else np.nan
            ),
        })
    return pd.DataFrame(rows)


def region_slopes(panel, first_year, last_year):
    d = panel[panel["year"].between(first_year, last_year)].copy()
    rows = []
    for nuts_id, g in d.groupby("NUTS_ID"):
        g = g.sort_values("year")
        fit = linregress(g["year"], g["log_rad"])
        rows.append({
            "NUTS_ID": nuts_id,
            "CNTR_CODE": str(g["CNTR_CODE"].iloc[0]),
            "density_band": str(g["density_band"].iloc[0]),
            "reference_population_density": float(
                g["reference_population_density"].iloc[0]
            ),
            "n_years": int(len(g)),
            "observed_log10_slope_per_year": float(fit.slope),
            "observed_pct_change_per_year": float(
                100.0 * (10.0 ** fit.slope - 1.0)
            ),
        })
    return pd.DataFrame(rows)


def comparison_table(t20, t24):
    a = t20[[
        "density_band",
        "pct_radiance_change_per_year",
        "pct_change_per_year_ci_low",
        "pct_change_per_year_ci_high",
        "trend_p",
    ]].rename(columns={
        "pct_radiance_change_per_year": "pct_per_year_2014_2020",
        "pct_change_per_year_ci_low": "ci_low_2014_2020",
        "pct_change_per_year_ci_high": "ci_high_2014_2020",
        "trend_p": "p_2014_2020",
    })
    b = t24[[
        "density_band",
        "pct_radiance_change_per_year",
        "pct_change_per_year_ci_low",
        "pct_change_per_year_ci_high",
        "trend_p",
    ]].rename(columns={
        "pct_radiance_change_per_year": "pct_per_year_2014_2024",
        "pct_change_per_year_ci_low": "ci_low_2014_2024",
        "pct_change_per_year_ci_high": "ci_high_2014_2024",
        "trend_p": "p_2014_2024",
    })
    out = a.merge(b, on="density_band", validate="one_to_one")
    out["change_2024_vs_2020_pp_per_year"] = (
        out["pct_per_year_2014_2024"]
        - out["pct_per_year_2014_2020"]
    )
    return out


def plot_trend_comparison(t20, t24, config, output):
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    x = np.arange(len(labels), dtype=float)
    fig, ax = plt.subplots(figsize=(9.5, 5.8))

    for trends, xpos, label, marker in [
        (t20, x - 0.13, "2014–2020", "o"),
        (t24, x + 0.13, "2014–2024", "s"),
    ]:
        g = trends.set_index("density_band").loc[labels].reset_index()
        est = g["pct_radiance_change_per_year"].to_numpy(float)
        lo = g["pct_change_per_year_ci_low"].to_numpy(float)
        hi = g["pct_change_per_year_ci_high"].to_numpy(float)
        ax.errorbar(
            xpos, est, yerr=[est-lo, hi-est], fmt=marker,
            capsize=4, linewidth=1.2, label=label
        )

    ax.axhline(0, linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Fixed population-density class (people km$^{-2}$)")
    ax.set_ylabel(
        "Country-adjusted change in calibrated VIIRS DNB radiance (% yr$^{-1}$)"
    )
    ax.set_title(
        "Fixed-NUTS2016 density trends: 2014–2020 versus 2014–2024"
    )
    ax.grid(axis="y", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_trajectories(trajectories, t24, config, output):
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    lookup = t24.set_index("density_band")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0), sharex=True)
    axes = axes.ravel()

    for i, band in enumerate(labels):
        ax = axes[i]
        g = trajectories[
            trajectories["density_band"] == str(band)
        ].sort_values("year")
        years = g["year"].to_numpy(float)
        rate = float(
            lookup.loc[str(band), "pct_radiance_change_per_year"]
        )
        ax.plot(
            years, g["observed_index_2014_100"],
            marker="o", linewidth=1.5, label="Observed fixed NUTS2016 data"
        )
        ax.fill_between(
            years, g["observed_index_ci_low"], g["observed_index_ci_high"],
            alpha=0.18
        )
        fitted = 100.0 * (1.0 + rate / 100.0) ** (years - 2014)
        ax.plot(
            years, fitted, linestyle="--", linewidth=1.8,
            label="2014–2024 fitted class trend"
        )
        ax.axhline(100, linewidth=0.7)
        ax.set_title(
            f"{band} people km$^{{-2}}$\n{rate:+.2f}% yr$^{{-1}}$"
        )
        ax.grid(alpha=0.2)

    axes[5].axis("off")
    handles, labs = axes[0].get_legend_handles_labels()
    axes[5].legend(handles, labs, loc="center", frameon=False)
    fig.supylabel("Calibrated VIIRS DNB radiance index (2014 = 100)")
    fig.supxlabel("Year")
    fig.suptitle(
        "Observed 2014–2024 radiance trajectories on fixed NUTS2016 geography",
        fontsize=14,
    )
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.95))
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_maps(s20, s24, config, output):
    cfg = config["analysis"]["fixed_nuts2016_2014_2024"]
    release = int(cfg["nuts_release"])
    nuts = load_fixed_nuts(config, release)

    west, south, east, north = map(
        float, config["analysis"]["europe_maps"]["extent_4326"]
    )
    europe_box = gpd.GeoDataFrame(
        geometry=[box(west, south, east, north)],
        crs="EPSG:4326",
    )
    nuts = gpd.clip(nuts, europe_box)

    a = s20[["NUTS_ID", "observed_pct_change_per_year"]].rename(
        columns={"observed_pct_change_per_year": "trend_2014_2020"}
    )
    b = s24[["NUTS_ID", "observed_pct_change_per_year"]].rename(
        columns={"observed_pct_change_per_year": "trend_2014_2024"}
    )
    vals = a.merge(b, on="NUTS_ID", validate="one_to_one")
    vals["change_after_2020"] = (
        vals["trend_2014_2024"] - vals["trend_2014_2020"]
    )

    mapped = nuts.merge(
        vals, on="NUTS_ID", how="left", validate="one_to_one"
    ).to_crs("EPSG:3035")
    box3035 = europe_box.to_crs("EPSG:3035")
    xmin, ymin, xmax, ymax = box3035.total_bounds

    q = float(cfg["map_absolute_quantile"])
    common = np.concatenate([
        vals["trend_2014_2020"].to_numpy(float),
        vals["trend_2014_2024"].to_numpy(float),
    ])
    common = common[np.isfinite(common)]
    lim = max(float(np.quantile(np.abs(common), q)), 0.25)

    diff = vals["change_after_2020"].to_numpy(float)
    diff = diff[np.isfinite(diff)]
    diff_lim = max(float(np.quantile(np.abs(diff), q)), 0.25)

    cmap = "coolwarm"
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
    diff_norm = TwoSlopeNorm(
        vmin=-diff_lim, vcenter=0, vmax=diff_lim
    )

    panels = [
        ("trend_2014_2020", "Observed trend\n2014–2020", norm),
        ("trend_2014_2024", "Observed trend\n2014–2024", norm),
        (
            "change_after_2020",
            "Change when 2021–2024 are added\n(2014–2024 minus 2014–2020)",
            diff_norm,
        ),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(17, 7))
    for ax, (col, title, this_norm) in zip(axes, panels):
        mapped.plot(
            ax=ax, color="0.92", edgecolor="0.78", linewidth=0.15
        )
        good = mapped[np.isfinite(mapped[col])].copy()
        good.plot(
            column=col, ax=ax, cmap=cmap, norm=this_norm,
            edgecolor="white", linewidth=0.08
        )
        if "CNTR_CODE" in mapped.columns:
            mapped.dissolve(by="CNTR_CODE").boundary.plot(
                ax=ax, color="0.18", linewidth=0.35
            )
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_axis_off()
        ax.set_title(title)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(
        sm, ax=axes[:2].tolist(), orientation="horizontal",
        fraction=0.045, pad=0.035
    )
    cb.set_label(
        "Observed annual change in calibrated VIIRS DNB radiance (% yr$^{-1}$)"
    )

    sm2 = ScalarMappable(norm=diff_norm, cmap=cmap)
    sm2.set_array([])
    cb2 = fig.colorbar(
        sm2, ax=[axes[2]], orientation="horizontal",
        fraction=0.045, pad=0.035
    )
    cb2.set_label(
        "Change in fitted trend (percentage points yr$^{-1}$)"
    )

    fig.suptitle(
        "Effect of extending the fixed-geography record from 2020 to 2024",
        fontsize=15, y=0.97
    )
    fig.savefig(
        output, dpi=int(cfg["figure_dpi"]), bbox_inches="tight"
    )
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()
    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_2014_2024"]

    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    compare_last = int(cfg["comparison_last_year"])
    release = int(cfg["nuts_release"])

    years = list(range(first_year, last_year + 1))
    rasters = discover_rasters(config, years)
    nuts = load_fixed_nuts(config, release)

    frames = [
        extract_fixed_year(
            rasters[year], year, nuts, release, config
        )
        for year in years
    ]
    fixed = pd.concat(frames, ignore_index=True)
    fixed = fixed.sort_values(["year", "NUTS_ID"]).reset_index(drop=True)

    validation = validate_against_production(
        fixed, config, first_year, compare_last, release
    )

    reference = build_population_reference(
        config,
        release,
        int(cfg["population_reference_first_year"]),
        int(cfg["population_reference_last_year"]),
    )
    panel = prepare_analysis_panel(fixed, reference, config)

    t20, info20 = fit_window(
        panel, first_year, compare_last, config
    )
    t24, info24 = fit_window(
        panel, first_year, last_year, config
    )
    comparison = comparison_table(t20, t24)
    trajectories = observed_trajectories(panel, baseline_year=first_year)
    s20 = region_slopes(panel, first_year, compare_last)
    s24 = region_slopes(panel, first_year, last_year)

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_2014_2024"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    paths = {
        "fixed": processed / "fixed_nuts2016_viirs_2014_2024.parquet",
        "fixed_csv": processed / "fixed_nuts2016_viirs_2014_2024.csv",
        "validation": processed
        / "fixed_nuts2016_2014_2020_production_validation.csv",
        "reference": processed
        / "fixed_nuts2016_population_density_reference.csv",
        "comparison": processed
        / "fixed_nuts2016_density_trends_2014_2020_vs_2014_2024.csv",
        "trajectories": processed
        / "fixed_nuts2016_density_trajectories_2014_2024.csv",
        "region_slopes": processed
        / "fixed_nuts2016_region_slopes_2014_2024.csv",
        "summary": processed
        / "fixed_nuts2016_2014_2024_summary.txt",
    }

    fixed.to_parquet(paths["fixed"], index=False)
    fixed.to_csv(paths["fixed_csv"], index=False)
    validation.to_csv(paths["validation"], index=False)
    reference.to_csv(paths["reference"], index=False)
    comparison.to_csv(paths["comparison"], index=False)
    trajectories.to_csv(paths["trajectories"], index=False)
    s24.to_csv(paths["region_slopes"], index=False)

    fig_compare = figures / "density_trends_2014_2020_vs_2014_2024.pdf"
    fig_traj = figures / "real_trajectories_2014_2024.pdf"
    fig_maps = figures / "maps_2014_2020_vs_2014_2024.pdf"

    plot_trend_comparison(t20, t24, config, fig_compare)
    plot_trajectories(trajectories, t24, config, fig_traj)
    plot_maps(s20, s24, config, fig_maps)

    with paths["summary"].open("w", encoding="utf-8") as h:
        h.write(
            "FIXED NUTS2016 VALIDATION OF 2014-2024 CALIBRATED VIIRS TRENDS\n"
        )
        h.write("=" * 80 + "\n\n")
        h.write(
            f"Fixed geography: NUTS {release}; years {first_year}-{last_year}.\n"
        )
        h.write(
            f"Population-density classes: median Eurostat density over "
            f"{cfg['population_reference_first_year']}-"
            f"{cfg['population_reference_last_year']}.\n"
        )
        h.write(
            f"Complete analysis panel: {panel['NUTS_ID'].nunique()} NUTS3; "
            f"{panel['CNTR_CODE'].nunique()} countries; {len(panel)} annual rows.\n"
        )
        h.write(
            f"Production cross-check 2014-2020: {len(validation)} matched rows; "
            f"maximum absolute corrected-radiance difference="
            f"{validation['abs_difference'].max():.12g}.\n\n"
        )
        h.write("COUNTRY-ADJUSTED DENSITY-CLASS TREND COMPARISON\n")
        h.write(
            comparison.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\nDENSITY-BAND × YEAR TESTS\n")
        h.write(
            f"2014-2020: LR={info20['lr']:.6g}, "
            f"df={info20['df']}, p={info20['p']:.6g}\n"
        )
        h.write(
            f"2014-2024: LR={info24['lr']:.6g}, "
            f"df={info24['df']}, p={info24['p']:.6g}\n"
        )
        h.write(
            "\nINTERPRETATION\n"
            "- This comparison removes both the 2013 endpoint and changing "
            "post-2020 NUTS geography.\n"
            "- The 2014-2020 and 2014-2024 estimates use the same fixed "
            "NUTS2016 regions and the same fixed population-density classes.\n"
            "- Any change in the class-specific trends after extending to 2024 "
            "therefore reflects the additional VIIRS observations rather than "
            "a change in administrative geography or class membership.\n"
        )

    print("\nFIXED NUTS2016 2014-2024 TREND COMPARISON")
    print(
        comparison.to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print(
        f"\n2014-2020 production cross-check max abs difference: "
        f"{validation['abs_difference'].max():.12g}"
    )
    print("\nWROTE")
    for p in [*paths.values(), fig_compare, fig_traj, fig_maps]:
        print(p)


if __name__ == "__main__":
    main()
