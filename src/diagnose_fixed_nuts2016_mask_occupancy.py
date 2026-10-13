#!/usr/bin/env python3
"""Audit annual masked-lit area occupancy on fixed NUTS2016 geography.

Scientific purpose
------------------
The 2021 annual VNL pulse is unusually large in absolute radiance units, but it
is not perfectly additive across the full brightness range. Because
average_masked sets background pixels to zero, a direct way to test whether
the 2021 event is primarily a MASK-OCCUPANCY problem is to measure, for every
fixed NUTS2016 region and every year, the physical-area-weighted fraction of
the masked annual raster that is positive.

If 2020->2021 shows an exceptional increase in positive masked area that
reverses in 2022 and covaries with the radiance pulse, that supports a
mask/background-classification contribution. If positive masked area remains
stable while radiance jumps, the main disturbance lies within already-positive
pixels rather than in the binary mask occupancy.

This stage:
- reuses the same locked 1,514-region NUTS2016 cohort;
- uses the existing raw EOG average_masked rasters only;
- does not modify any raster or production output;
- caches each expensive annual zonal extraction immediately.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from scipy.stats import spearmanr

from apply_radiance_calibration import CalibrationRasterSource
from diagnose_fixed_nuts2016_2021_reversal import (
    add_latitude,
    load_locked_panel,
)
from extract_nuts3_radiance import discover_rasters, validate_raster
from project_config import configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import load_fixed_nuts


YEARS = list(range(2014, 2025))
Q_LABELS = ["Q1 darkest", "Q2", "Q3", "Q4", "Q5 brightest"]


def attach_fixed_brightness(panel: pd.DataFrame) -> pd.DataFrame:
    base = (
        panel[panel["year"].isin([2014, 2015, 2016])]
        .groupby("NUTS_ID", as_index=False)
        .agg(baseline_raw_nw=("radiance_mean_raw", "median"))
    )
    if len(base) != panel["NUTS_ID"].nunique():
        raise RuntimeError("Incomplete 2014-2016 baseline brightness reference.")
    if (~np.isfinite(base["baseline_raw_nw"]) | (base["baseline_raw_nw"] <= 0)).any():
        raise RuntimeError("Invalid baseline raw radiance.")

    base["brightness_quintile"] = pd.qcut(
        base["baseline_raw_nw"],
        q=5,
        labels=Q_LABELS,
        duplicates="raise",
    )
    out = panel.merge(base, on="NUTS_ID", how="left", validate="many_to_one")
    if out["brightness_quintile"].isna().any():
        raise RuntimeError("Brightness quintile merge failed.")
    return out


def normalize_exactextract_mean(
    result: pd.DataFrame,
    id_field: str,
    out_col: str,
) -> pd.DataFrame:
    """Normalize single-source exactextract output and validate uniqueness."""
    if out_col not in result.columns:
        if "mean" in result.columns:
            result = result.rename(columns={"mean": out_col})
        else:
            raise RuntimeError(
                f"Expected 'mean' or {out_col}; found {list(result.columns)}"
            )
    result[id_field] = result[id_field].astype(str)
    result = result[[id_field, out_col]].copy()
    if result[id_field].duplicated().any():
        raise RuntimeError("Duplicate NUTS_ID from mask-occupancy extraction.")
    return result


def validate_cached(
    cached: pd.DataFrame,
    year: int,
    keep_ids: set[str],
    source_name: str,
) -> pd.DataFrame:
    required = {
        "year",
        "NUTS_ID",
        "positive_masked_area_fraction",
        "viirs_file",
    }
    if not required.issubset(cached.columns):
        raise RuntimeError(
            f"Cached {year} mask occupancy has invalid schema: "
            f"{list(cached.columns)}"
        )
    out = cached[list(required)].copy()
    out["NUTS_ID"] = out["NUTS_ID"].astype(str)
    if set(out["year"].astype(int)) != {year}:
        raise RuntimeError(f"Cached mask occupancy has wrong year for {year}.")
    if set(out["NUTS_ID"]) != keep_ids:
        raise RuntimeError(f"Cached {year} mask occupancy has wrong cohort.")
    if out["NUTS_ID"].duplicated().any():
        raise RuntimeError(f"Cached {year} mask occupancy has duplicate regions.")
    if set(out["viirs_file"].astype(str)) != {source_name}:
        raise RuntimeError(
            f"Cached {year} mask occupancy comes from a different raster."
        )
    vals = pd.to_numeric(out["positive_masked_area_fraction"], errors="coerce")
    if (~np.isfinite(vals) | (vals < -1e-8) | (vals > 1 + 1e-8)).any():
        raise RuntimeError(f"Invalid cached positive-area fraction in {year}.")
    out["positive_masked_area_fraction"] = vals
    return out


def extract_positive_fraction(
    path: Path,
    year: int,
    nuts,
    keep_ids: set[str],
    config: dict,
    cache_dir: Path,
) -> pd.DataFrame:
    cache = cache_dir / f"nuts2016_{year}_positive_masked_area.parquet"
    if cache.exists():
        print(f"Using cached {year} mask occupancy: {cache}", flush=True)
        return validate_cached(pd.read_parquet(cache), year, keep_ids, path.name)

    print(f"Extracting {year} positive masked area: {path.name}", flush=True)
    validate_raster(path, config["nuts"]["crs"])
    id_field = config["nuts"]["id_field"]

    with rasterio.open(path) as ds:
        source = CalibrationRasterSource(
            ds,
            offset=0.0,
            mode="raw_positive_indicator",
            name="positive_masked",
        )
        result = exact_extract(
            source,
            nuts,
            ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=[id_field],
            output="pandas",
            progress=True,
        )

    result = normalize_exactextract_mean(
        result,
        id_field,
        "positive_masked_area_fraction",
    )
    result = result[result[id_field].isin(keep_ids)].copy()
    if set(result[id_field]) != keep_ids:
        raise RuntimeError(
            f"{year} positive-area extraction does not match locked cohort."
        )

    vals = pd.to_numeric(
        result["positive_masked_area_fraction"], errors="coerce"
    )
    if (~np.isfinite(vals) | (vals < -1e-8) | (vals > 1 + 1e-8)).any():
        raise RuntimeError(f"Invalid positive-area fraction extracted in {year}.")

    result.insert(0, "year", year)
    result["viirs_file"] = path.name
    result.to_parquet(cache, index=False)
    print(f"Cached {year} mask occupancy: {cache}", flush=True)
    return result


def extract_all(panel: pd.DataFrame, config: dict) -> pd.DataFrame:
    rasters = discover_rasters(config, YEARS)
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release)
    keep_ids = set(panel["NUTS_ID"].astype(str))

    cache_dir = configured_path(config, "processed") / "viirs_mask_occupancy_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    pieces = [
        extract_positive_fraction(
            rasters[year],
            year,
            nuts,
            keep_ids,
            config,
            cache_dir,
        )
        for year in YEARS
    ]
    out = pd.concat(pieces, ignore_index=True)
    expected = len(keep_ids) * len(YEARS)
    if len(out) != expected:
        raise RuntimeError(
            f"Mask occupancy table has {len(out)} rows; expected {expected}."
        )
    return out


def make_transitions(panel: pd.DataFrame, occupancy: pd.DataFrame) -> pd.DataFrame:
    d = panel.merge(
        occupancy[
            ["year", "NUTS_ID", "positive_masked_area_fraction"]
        ],
        on=["year", "NUTS_ID"],
        how="left",
        validate="one_to_one",
    ).sort_values(["NUTS_ID", "year"])

    d["previous_positive_fraction"] = d.groupby("NUTS_ID")[
        "positive_masked_area_fraction"
    ].shift(1)
    d["delta_positive_fraction"] = (
        d["positive_masked_area_fraction"] - d["previous_positive_fraction"]
    )
    d["delta_positive_percentage_points"] = 100.0 * d["delta_positive_fraction"]

    d["previous_raw_nw"] = d.groupby("NUTS_ID")["radiance_mean_raw"].shift(1)
    d["delta_abs_raw_nw"] = d["radiance_mean_raw"] - d["previous_raw_nw"]
    d["previous_year"] = d.groupby("NUTS_ID")["year"].shift(1)

    t = d[d["previous_year"].notna()].copy()
    if ((t["year"].astype(int) - t["previous_year"].astype(int)) != 1).any():
        raise RuntimeError("Non-consecutive annual transition.")
    return t


def _spearman(a: pd.Series, b: pd.Series) -> tuple[float, float]:
    ok = np.isfinite(a.to_numpy(float)) & np.isfinite(b.to_numpy(float))
    r = spearmanr(a[ok], b[ok])
    return float(r.statistic), float(r.pvalue)


def country_summary(trans: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (year, country), g in trans.groupby(["year", "CNTR_CODE"], observed=True):
        rows.append(
            {
                "end_year": int(year),
                "CNTR_CODE": str(country),
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "median_delta_positive_pp": float(
                    g["delta_positive_percentage_points"].median()
                ),
                "median_delta_raw_nw": float(g["delta_abs_raw_nw"].median()),
            }
        )
    return pd.DataFrame(rows)


def annual_summary(
    trans: pd.DataFrame,
    countries: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for year, g in trans.groupby("year", observed=True):
        c = countries[countries["end_year"] == int(year)]
        rho, p = _spearman(
            g["delta_positive_percentage_points"],
            g["delta_abs_raw_nw"],
        )
        rows.append(
            {
                "end_year": int(year),
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "median_delta_positive_pp": float(
                    g["delta_positive_percentage_points"].median()
                ),
                "mean_delta_positive_pp": float(
                    g["delta_positive_percentage_points"].mean()
                ),
                "pct_regions_mask_area_increasing": float(
                    100.0 * (g["delta_positive_percentage_points"] > 0).mean()
                ),
                "country_balanced_mean_median_delta_positive_pp": float(
                    c["median_delta_positive_pp"].mean()
                ),
                "country_median_delta_positive_pp": float(
                    c["median_delta_positive_pp"].median()
                ),
                "median_delta_raw_nw": float(g["delta_abs_raw_nw"].median()),
                "country_balanced_mean_median_delta_raw_nw": float(
                    c["median_delta_raw_nw"].mean()
                ),
                "rho_mask_change_vs_abs_radiance_change": rho,
                "p_mask_change_vs_abs_radiance_change": p,
            }
        )

    out = pd.DataFrame(rows).sort_values("end_year").reset_index(drop=True)
    key = "country_balanced_mean_median_delta_positive_pp"
    out["abs_mask_shift_rank"] = (
        out[key].abs().rank(method="min", ascending=False).astype(int)
    )
    return out


def grouped_summary(trans: pd.DataFrame, group: str) -> pd.DataFrame:
    rows = []
    for (year, label), g in trans.groupby(["year", group], observed=True):
        rows.append(
            {
                "end_year": int(year),
                group: str(label),
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "median_delta_positive_pp": float(
                    g["delta_positive_percentage_points"].median()
                ),
                "median_delta_raw_nw": float(g["delta_abs_raw_nw"].median()),
            }
        )
    return pd.DataFrame(rows)


def synthetic_self_test() -> None:
    rows = []
    occ = []
    for i in range(20):
        nuts = f"R{i:02d}"
        country = "AA" if i < 10 else "BB"
        baseline = 0.2 + i * 0.2
        for year in YEARS:
            frac = 0.1 + i * 0.002
            raw = baseline
            if year == 2021:
                frac += 0.03
                raw += 0.05
            if year == 2022:
                frac += 0.005
                raw += 0.01
            rows.append(
                {
                    "NUTS_ID": nuts,
                    "CNTR_CODE": country,
                    "density_band": "<30" if i < 10 else "30-100",
                    "year": year,
                    "radiance_mean_raw": raw,
                    "baseline_raw_nw": baseline,
                    "brightness_quintile": Q_LABELS[min(i // 4, 4)],
                    "centroid_latitude": 45 + i * 0.2,
                }
            )
            occ.append(
                {
                    "year": year,
                    "NUTS_ID": nuts,
                    "positive_masked_area_fraction": frac,
                }
            )

    panel = pd.DataFrame(rows)
    occupancy = pd.DataFrame(occ)
    trans = make_transitions(panel, occupancy)
    countries = country_summary(trans)
    annual = annual_summary(trans, countries)
    by_brightness = grouped_summary(trans, "brightness_quintile")
    by_density = grouped_summary(trans, "density_band")

    a21 = annual.loc[annual["end_year"] == 2021].iloc[0]
    if a21["median_delta_positive_pp"] <= 2.5:
        raise RuntimeError("Self-test failed to recover synthetic mask jump.")
    if a21["median_delta_raw_nw"] <= 0.04:
        raise RuntimeError("Self-test failed to recover synthetic radiance jump.")
    if by_brightness.empty or by_density.empty:
        raise RuntimeError("Self-test grouped summaries are empty.")
    print("Mask-occupancy downstream synthetic self-test: PASS", flush=True)


def make_plots(
    trans: pd.DataFrame,
    annual: pd.DataFrame,
    by_brightness: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    ax.plot(
        annual["end_year"],
        annual["country_balanced_mean_median_delta_positive_pp"],
        marker="o",
    )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(2021, linestyle="--", linewidth=0.9)
    ax.set_xlabel("Transition ending in year")
    ax.set_ylabel("Country-balanced median change in positive masked area (percentage points)")
    ax.set_title("Annual change in masked-positive area on fixed NUTS2016 geography")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "annual_mask_occupancy_change_2015_2024.pdf", dpi=dpi)
    plt.close(fig)

    g = trans[trans["year"] == 2021]
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    ax.scatter(
        g["delta_positive_percentage_points"],
        g["delta_abs_raw_nw"],
        s=12,
        alpha=0.4,
    )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel("Change in positive masked area, 2020→2021 (percentage points)")
    ax.set_ylabel("Raw absolute radiance change, 2020→2021 (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Does the 2021 radiance pulse track mask occupancy?")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "radiance_vs_mask_change_2021.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    for label in Q_LABELS:
        g = by_brightness[
            by_brightness["brightness_quintile"].astype(str) == label
        ]
        ax.plot(
            g["end_year"],
            g["median_delta_positive_pp"],
            marker="o",
            label=label,
        )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(2021, linestyle="--", linewidth=0.9)
    ax.set_xlabel("Transition ending in year")
    ax.set_ylabel("Median change in positive masked area (percentage points)")
    ax.set_title("Mask-occupancy change by fixed baseline brightness quintile")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "mask_change_by_baseline_brightness.pdf", dpi=dpi)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    synthetic_self_test()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_2021_diagnostic"]
    panel = load_locked_panel(config)
    panel = add_latitude(panel, config)
    panel = attach_fixed_brightness(panel)

    occupancy = extract_all(panel, config)
    trans = make_transitions(panel, occupancy)
    countries = country_summary(trans)
    annual = annual_summary(trans, countries)
    by_brightness = grouped_summary(trans, "brightness_quintile")
    by_density = grouped_summary(trans, "density_band")

    expected_n = int(cfg["expected_n_nuts3"])
    if trans["NUTS_ID"].nunique() != expected_n:
        raise RuntimeError("Unexpected NUTS3 cohort after mask audit.")
    if len(trans) != expected_n * 10:
        raise RuntimeError(
            f"Expected {expected_n * 10} transition rows; found {len(trans)}."
        )

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_mask_occupancy_audit"
    )
    processed.mkdir(parents=True, exist_ok=True)

    outputs = {
        "annual_occupancy": processed
        / "fixed_nuts2016_mask_positive_area_2014_2024.csv",
        "transitions": processed
        / "fixed_nuts2016_mask_occupancy_transitions_2015_2024.csv",
        "annual": processed
        / "fixed_nuts2016_mask_occupancy_annual_summary_2015_2024.csv",
        "country": processed
        / "fixed_nuts2016_mask_occupancy_country_medians_2015_2024.csv",
        "brightness": processed
        / "fixed_nuts2016_mask_occupancy_by_brightness_2015_2024.csv",
        "density": processed
        / "fixed_nuts2016_mask_occupancy_by_density_2015_2024.csv",
        "summary": processed / "fixed_nuts2016_mask_occupancy_audit_summary.txt",
    }

    occupancy.to_csv(outputs["annual_occupancy"], index=False)
    trans.to_csv(outputs["transitions"], index=False)
    annual.to_csv(outputs["annual"], index=False)
    countries.to_csv(outputs["country"], index=False)
    by_brightness.to_csv(outputs["brightness"], index=False)
    by_density.to_csv(outputs["density"], index=False)

    make_plots(
        trans,
        annual,
        by_brightness,
        figures,
        int(cfg["figure_dpi"]),
    )

    a21 = annual.loc[annual["end_year"] == 2021].iloc[0]
    a22 = annual.loc[annual["end_year"] == 2022].iloc[0]

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 MASK-OCCUPANCY AUDIT: 2014-2024\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            f"Locked panel: {trans['NUTS_ID'].nunique()} NUTS3 regions; "
            f"{trans['CNTR_CODE'].nunique()} countries.\n"
        )
        h.write(
            "Mask diagnostic: physical-area-weighted fraction of each NUTS3 "
            "with raw average_masked radiance > 0.\n"
        )
        h.write(
            "Interpretation: a jump in this fraction means more polygon area "
            "is classified/retained as positive in the masked annual raster.\n\n"
        )

        h.write("ANNUAL MASK-OCCUPANCY SUMMARY\n")
        h.write(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nMASK CHANGE BY FIXED BASELINE BRIGHTNESS\n")
        h.write(
            by_brightness.to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )
        h.write("\n\nMASK CHANGE BY FIXED POPULATION-DENSITY BAND\n")
        h.write(
            by_density.to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )

        h.write("\n\nPRE-SPECIFIED 2021/2022 CHECKS\n")
        h.write(
            f"- 2021 country-balanced median mask-area shift: "
            f"{a21['country_balanced_mean_median_delta_positive_pp']:+.6f} pp.\n"
        )
        h.write(
            f"- 2021 absolute mask-shift rank among 2015-2024 transitions: "
            f"{int(a21['abs_mask_shift_rank'])} (1 = largest).\n"
        )
        h.write(
            f"- 2021 mask-change versus absolute-radiance-change Spearman rho: "
            f"{a21['rho_mask_change_vs_abs_radiance_change']:+.4f}.\n"
        )
        h.write(
            f"- 2022 country-balanced median mask-area shift: "
            f"{a22['country_balanced_mean_median_delta_positive_pp']:+.6f} pp.\n"
        )
        h.write(
            f"- 2022 mask-change versus absolute-radiance-change Spearman rho: "
            f"{a22['rho_mask_change_vs_abs_radiance_change']:+.4f}.\n"
        )

        h.write(
            "\nDecision rule:\n"
            "- If 2021 has an exceptional positive mask-area shift, followed "
            "by a substantial reversal in 2022, and region-level mask changes "
            "track the radiance pulse, masking/background occupancy is a likely "
            "major contributor.\n"
            "- If mask occupancy is stable or weakly associated with the "
            "radiance pulse, the 2021 disturbance primarily reflects radiance "
            "values within already-positive pixels; the next check should then "
            "compare masked with unmasked annual layers or monthly inputs.\n"
            "- No 2021 correction is justified by this diagnostic alone.\n"
        )

    print("\nFIXED NUTS2016 MASK-OCCUPANCY AUDIT")
    print(
        f"Locked panel: {trans['NUTS_ID'].nunique()} NUTS3; "
        f"{trans['CNTR_CODE'].nunique()} countries"
    )
    print("\nANNUAL MASK-OCCUPANCY SUMMARY")
    print(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\n2021/2022 CHECK")
    print(
        annual[annual["end_year"].isin([2021, 2022])].to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
