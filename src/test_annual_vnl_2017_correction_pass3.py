#!/usr/bin/env python3
"""Third-pass 2017 VIIRS calibration diagnostic using segmented time series.

Pass 2 showed a remarkably flat retained-pixel 2017 excess of about
0.23-0.24 nW cm-2 sr-1 for independent baseline radiances below ~5, followed
by a breakdown of the simple one-step pattern at brighter radiances.

This stage asks whether that pattern survives when the calibration step is
estimated from the full pre/post time series rather than only the three
adjacent transitions.

Primary interval
----------------
2013-2019 on fixed NUTS2016 geography. 2020 is deliberately excluded from the
primary calibration estimate to avoid the COVID-era year influencing a sensor
step estimate.

For each NUTS3 and independent 2013-2014 brightness bin, pixels must be
positive in every year 2013-2019. The annual median radiance of that fixed
pixel set is then fitted with:

    R_t = a + b * (year - 2016) + delta * I(year >= 2017)

where delta is the additive step in nW cm-2 sr-1.

A log10 model is fitted to the same annual medians:

    log10(R_t) = a + b * (year - 2016) + delta_log * I(year >= 2017)

which gives a multiplicative step:

    100 * (10**delta_log - 1)

A sensitivity model allows the post-2017 slope to differ:

    R_t = a + b*t + delta*post + gamma*max(year-2017, 0)

The linear and log models let us distinguish an approximately additive sensor
step (constant delta across brightness) from a multiplicative step (constant
percentage across brightness).

Candidate correction shapes
----------------------------
Using the segmented estimates, the script compares:

- no correction;
- 0.15 nW constant;
- a pass-2 0.23478 nW constant;
- a pass-3 empirical constant estimated from baseline <=5 nW;
- the same empirical value only for baseline <=5 nW;
- the same value flat to 5 nW, tapering linearly to zero by 10 nW.

These are diagnostics, not production corrections. The next production step
should only be taken if one simple candidate removes the segmented step across
countries without creating systematic brightness-dependent residuals.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask, geometry_window
from scipy.stats import spearmanr
from shapely.geometry import mapping

from extract_nuts3_radiance import discover_rasters, nuts_path_for_release
from project_config import configured_path, load_config


def validate_alignment(sources: dict[int, rasterio.io.DatasetReader]) -> None:
    years = sorted(sources)
    ref = sources[years[0]]
    for year in years[1:]:
        src = sources[year]
        problems: list[str] = []
        if src.crs != ref.crs:
            problems.append(f"CRS {src.crs} != {ref.crs}")
        if (src.width, src.height) != (ref.width, ref.height):
            problems.append(
                f"shape {src.width}x{src.height} != {ref.width}x{ref.height}"
            )
        if not src.transform.almost_equals(ref.transform):
            problems.append("affine transform differs")
        if problems:
            raise ValueError(
                f"VIIRS raster grid for {year} is not aligned with "
                f"{years[0]}: " + "; ".join(problems)
            )


def read_region(
    geom,
    sources: dict[int, rasterio.io.DatasetReader],
) -> tuple[dict[int, np.ndarray], np.ndarray] | None:
    years = sorted(sources)
    ref = sources[years[0]]
    try:
        window = geometry_window(ref, [mapping(geom)])
    except Exception:
        return None

    if window.width <= 0 or window.height <= 0:
        return None

    shape = (int(window.height), int(window.width))
    transform = ref.window_transform(window)
    inside = geometry_mask(
        [mapping(geom)],
        out_shape=shape,
        transform=transform,
        invert=True,
        all_touched=False,
    )

    arrays: dict[int, np.ndarray] = {}
    for year, src in sources.items():
        arr = src.read(1, window=window, masked=True)
        arrays[year] = arr.filled(np.nan).astype(np.float32, copy=False)

    return arrays, inside


def fit_step_model(
    years: np.ndarray,
    y: np.ndarray,
    slope_break: bool = False,
) -> dict[str, float]:
    t = years.astype(float) - 2016.0
    post = (years >= 2017).astype(float)

    if slope_break:
        post_time = np.maximum(years.astype(float) - 2017.0, 0.0)
        X = np.column_stack([np.ones_like(t), t, post, post_time])
        step_idx = 2
    else:
        X = np.column_stack([np.ones_like(t), t, post])
        step_idx = 2

    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    fitted = X @ beta
    residual = y - fitted
    rmse = float(np.sqrt(np.mean(residual**2)))

    return {
        "step": float(beta[step_idx]),
        "slope": float(beta[1]),
        "rmse": rmse,
        "n_years": int(len(years)),
    }


def region_bin_models(
    nuts_id: str,
    country: str,
    geom,
    sources: dict[int, rasterio.io.DatasetReader],
    baseline_years: tuple[int, int],
    bin_edges: np.ndarray,
    bin_labels: list[str],
    min_bin_pixels: int,
) -> list[dict]:
    read = read_region(geom, sources)
    if read is None:
        return []

    arrays, inside = read
    years = np.array(sorted(sources), dtype=int)

    support = inside.copy()
    for year in years:
        support &= np.isfinite(arrays[int(year)])
        support &= arrays[int(year)] > 0

    if not support.any():
        return []

    b0, b1 = baseline_years
    baseline = 0.5 * (
        arrays[b0][support].astype(np.float64, copy=False)
        + arrays[b1][support].astype(np.float64, copy=False)
    )

    annual = {
        int(year): arrays[int(year)][support].astype(np.float64, copy=False)
        for year in years
    }

    bin_index = np.digitize(
        baseline,
        bin_edges[1:-1],
        right=False,
    )

    rows: list[dict] = []

    for i, label in enumerate(bin_labels):
        take = bin_index == i
        n_pixels = int(take.sum())
        if n_pixels < min_bin_pixels:
            continue

        medians = np.array(
            [np.median(annual[int(year)][take]) for year in years],
            dtype=float,
        )

        linear = fit_step_model(years, medians, slope_break=False)
        linear_break = fit_step_model(years, medians, slope_break=True)
        log = fit_step_model(years, np.log10(medians), slope_break=False)

        rows.append(
            {
                "NUTS_ID": nuts_id,
                "CNTR_CODE": country,
                "brightness_bin": label,
                "n_pixels": n_pixels,
                "median_baseline_2013_2014": float(np.median(baseline[take])),
                "step_additive_nW": linear["step"],
                "step_additive_slopebreak_nW": linear_break["step"],
                "prepost_common_slope_nW_per_year": linear["slope"],
                "linear_rmse": linear["rmse"],
                "linear_slopebreak_rmse": linear_break["rmse"],
                "step_log10": log["step"],
                "step_multiplicative_pct": 100.0 * (10.0 ** log["step"] - 1.0),
                "log_rmse": log["rmse"],
                **{
                    f"median_radiance_{int(year)}": float(value)
                    for year, value in zip(years, medians)
                },
            }
        )

    return rows


def aggregate_country_bins(
    region_bins: pd.DataFrame,
    min_regions: int,
) -> pd.DataFrame:
    rows = []

    for (country, brightness_bin), g in region_bins.groupby(
        ["CNTR_CODE", "brightness_bin"],
        observed=True,
    ):
        if len(g) < min_regions:
            continue

        rows.append(
            {
                "CNTR_CODE": country,
                "brightness_bin": brightness_bin,
                "n_regions": int(len(g)),
                "n_pixels": int(g["n_pixels"].sum()),
                "median_baseline_2013_2014": float(
                    g["median_baseline_2013_2014"].median()
                ),
                "median_additive_step_nW": float(
                    g["step_additive_nW"].median()
                ),
                "median_additive_step_slopebreak_nW": float(
                    g["step_additive_slopebreak_nW"].median()
                ),
                "median_multiplicative_step_pct": float(
                    g["step_multiplicative_pct"].median()
                ),
                "median_linear_rmse": float(g["linear_rmse"].median()),
                "median_linear_slopebreak_rmse": float(
                    g["linear_slopebreak_rmse"].median()
                ),
            }
        )

    return pd.DataFrame(rows)


def aggregate_europe_bins(
    country_bins: pd.DataFrame,
    bin_labels: list[str],
    min_countries: int,
) -> pd.DataFrame:
    rows = []

    for label in bin_labels:
        g = country_bins[country_bins["brightness_bin"] == label]
        if len(g) < min_countries:
            continue

        rows.append(
            {
                "brightness_bin": label,
                "n_countries": int(len(g)),
                "n_regions": int(g["n_regions"].sum()),
                "n_pixels": int(g["n_pixels"].sum()),
                "median_country_baseline": float(
                    g["median_baseline_2013_2014"].median()
                ),
                "median_country_additive_step_nW": float(
                    g["median_additive_step_nW"].median()
                ),
                "additive_step_q25": float(
                    g["median_additive_step_nW"].quantile(0.25)
                ),
                "additive_step_q75": float(
                    g["median_additive_step_nW"].quantile(0.75)
                ),
                "median_country_additive_step_slopebreak_nW": float(
                    g["median_additive_step_slopebreak_nW"].median()
                ),
                "median_country_multiplicative_step_pct": float(
                    g["median_multiplicative_step_pct"].median()
                ),
                "multiplicative_step_q25": float(
                    g["median_multiplicative_step_pct"].quantile(0.25)
                ),
                "multiplicative_step_q75": float(
                    g["median_multiplicative_step_pct"].quantile(0.75)
                ),
            }
        )

    return pd.DataFrame(rows)


def dependence(
    country_bins: pd.DataFrame,
    value_col: str,
    max_baseline: float | None = None,
) -> dict[str, float]:
    x = country_bins.dropna(
        subset=["median_baseline_2013_2014", value_col]
    ).copy()
    x = x[x["median_baseline_2013_2014"] > 0]

    if max_baseline is not None:
        x = x[x["median_baseline_2013_2014"] <= max_baseline]

    if len(x) < 3:
        return {
            "n_country_bin_cells": int(len(x)),
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
        }

    rho, p = spearmanr(
        np.log10(x["median_baseline_2013_2014"]),
        x[value_col],
    )

    return {
        "n_country_bin_cells": int(len(x)),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
    }


def empirical_step(
    country_bins: pd.DataFrame,
    max_baseline: float,
) -> tuple[float, pd.DataFrame]:
    x = country_bins[
        country_bins["median_baseline_2013_2014"] <= max_baseline
    ].copy()

    by_country = (
        x.groupby("CNTR_CODE", as_index=False)["median_additive_step_nW"]
        .median()
        .rename(columns={"median_additive_step_nW": "country_empirical_step"})
    )

    return float(by_country["country_empirical_step"].median()), by_country


def correction_value(
    scenario: str,
    baseline: np.ndarray,
    empirical: float,
    pass2: float,
    taper_start: float,
    taper_end: float,
) -> np.ndarray:
    b = np.asarray(baseline, dtype=float)

    if scenario == "none":
        return np.zeros_like(b)
    if scenario == "constant_0p15":
        return np.full_like(b, 0.15)
    if scenario == "constant_pass2_0p23478":
        return np.full_like(b, pass2)
    if scenario == "constant_pass3_empirical":
        return np.full_like(b, empirical)
    if scenario == "pass3_le5_only":
        return np.where(b <= taper_start, empirical, 0.0)
    if scenario == "pass3_taper_5_to_10":
        out = np.zeros_like(b)
        out[b <= taper_start] = empirical
        middle = (b > taper_start) & (b < taper_end)
        out[middle] = empirical * (
            (taper_end - b[middle]) / (taper_end - taper_start)
        )
        return out

    raise KeyError(f"Unknown scenario: {scenario}")


def evaluate_scenarios(
    country_bins: pd.DataFrame,
    empirical: float,
    pass2: float,
    taper_start: float,
    taper_end: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scenarios = [
        "none",
        "constant_0p15",
        "constant_pass2_0p23478",
        "constant_pass3_empirical",
        "pass3_le5_only",
        "pass3_taper_5_to_10",
    ]

    detail_rows = []

    for scenario in scenarios:
        x = country_bins.copy()
        corr = correction_value(
            scenario,
            x["median_baseline_2013_2014"].to_numpy(),
            empirical,
            pass2,
            taper_start,
            taper_end,
        )
        residual = x["median_additive_step_nW"].to_numpy() - corr

        for (_, row), c, r in zip(x.iterrows(), corr, residual):
            detail_rows.append(
                {
                    "scenario": scenario,
                    "CNTR_CODE": row["CNTR_CODE"],
                    "brightness_bin": row["brightness_bin"],
                    "median_baseline_2013_2014": row[
                        "median_baseline_2013_2014"
                    ],
                    "observed_segmented_step_nW": row[
                        "median_additive_step_nW"
                    ],
                    "correction_nW": float(c),
                    "residual_step_nW": float(r),
                    "abs_residual_step_nW": float(abs(r)),
                }
            )

    detail = pd.DataFrame(detail_rows)

    summary_rows = []

    for scenario, g in detail.groupby("scenario", sort=False):
        row = {
            "scenario": scenario,
            "median_residual_all": float(g["residual_step_nW"].median()),
            "median_abs_residual_all": float(
                g["abs_residual_step_nW"].median()
            ),
        }

        for label, mask in [
            ("le5", g["median_baseline_2013_2014"] <= taper_start),
            ("le10", g["median_baseline_2013_2014"] <= taper_end),
            ("gt10", g["median_baseline_2013_2014"] > taper_end),
        ]:
            z = g[mask]
            row[f"n_cells_{label}"] = int(len(z))
            row[f"median_residual_{label}"] = (
                float(z["residual_step_nW"].median())
                if len(z)
                else np.nan
            )
            row[f"median_abs_residual_{label}"] = (
                float(z["abs_residual_step_nW"].median())
                if len(z)
                else np.nan
            )

        # Primary ranking favours the brightness range in which pass 2
        # suggested an additive calibration step is identifiable.
        row["primary_score"] = (
            row["median_abs_residual_le10"]
            + 0.5 * row["median_abs_residual_all"]
        )
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows).sort_values("primary_score")
    summary["rank"] = np.arange(1, len(summary) + 1)

    return detail, summary


def plot_steps(
    europe_bins: pd.DataFrame,
    output: Path,
    dpi: int,
) -> None:
    b = europe_bins.sort_values("median_country_baseline").copy()

    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    ax.plot(
        b["median_country_baseline"],
        b["median_country_additive_step_nW"],
        marker="o",
        linewidth=1.8,
        label="Common-slope segmented step",
    )
    ax.plot(
        b["median_country_baseline"],
        b["median_country_additive_step_slopebreak_nW"],
        marker="s",
        linestyle="--",
        linewidth=1.5,
        label="Post-slope-break sensitivity",
    )
    ax.fill_between(
        b["median_country_baseline"],
        b["additive_step_q25"],
        b["additive_step_q75"],
        alpha=0.15,
        linewidth=0,
    )
    ax.axhline(0, color="0.4", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel(
        "Independent baseline radiance: mean 2013-2014 "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_ylabel(
        "Estimated 2017 additive step (nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title(
        "Segmented 2013-2019 estimate of the 2017 Annual VNL step"
    )
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_additive_vs_multiplicative(
    europe_bins: pd.DataFrame,
    output: Path,
    dpi: int,
) -> None:
    b = europe_bins.sort_values("median_country_baseline").copy()

    fig, ax1 = plt.subplots(figsize=(9.5, 6.2))
    ax1.plot(
        b["median_country_baseline"],
        b["median_country_additive_step_nW"],
        marker="o",
        linewidth=1.8,
        label="Additive step (left axis)",
    )
    ax1.axhline(0, color="0.4", linewidth=1)
    ax1.set_xscale("log")
    ax1.set_xlabel(
        "Independent baseline radiance: mean 2013-2014 "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax1.set_ylabel("Additive step (nW cm$^{-2}$ sr$^{-1}$)")

    ax2 = ax1.twinx()
    ax2.plot(
        b["median_country_baseline"],
        b["median_country_multiplicative_step_pct"],
        marker="s",
        linestyle="--",
        linewidth=1.5,
        label="Multiplicative step (right axis)",
    )
    ax2.set_ylabel("Multiplicative step (%)")

    lines = ax1.get_lines() + ax2.get_lines()
    labels = [line.get_label() for line in lines if not line.get_label().startswith("_")]
    lines = [line for line in lines if not line.get_label().startswith("_")]
    ax1.legend(lines, labels, loc="best")
    ax1.set_title("Additive versus multiplicative form of the 2017 step")
    ax1.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_candidate_residuals(
    detail: pd.DataFrame,
    output: Path,
    dpi: int,
) -> None:
    fig, ax = plt.subplots(figsize=(10.2, 6.5))

    grouped = (
        detail.groupby(
            ["scenario", "brightness_bin"],
            observed=True,
            as_index=False,
        )
        .agg(
            median_baseline=("median_baseline_2013_2014", "median"),
            median_residual=("residual_step_nW", "median"),
        )
    )

    for scenario, g in grouped.groupby("scenario", sort=False):
        g = g.sort_values("median_baseline")
        ax.plot(
            g["median_baseline"],
            g["median_residual"],
            marker="o",
            linewidth=1.4,
            label=scenario,
        )

    ax.axhline(0, color="0.4", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel(
        "Independent baseline radiance: mean 2013-2014 "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_ylabel("Residual segmented 2017 step (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Candidate correction shapes: residual 2017 step")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    path: Path,
    europe_bins: pd.DataFrame,
    additive_all: dict[str, float],
    additive_le5: dict[str, float],
    multiplicative_all: dict[str, float],
    multiplicative_le5: dict[str, float],
    empirical: float,
    scenario_summary: pd.DataFrame,
    cfg: dict,
) -> None:
    with path.open("w", encoding="utf-8") as h:
        h.write("2017 VIIRS CALIBRATION PASS 3: SEGMENTED 2013-2019 TIME SERIES\n")
        h.write("=" * 78 + "\n\n")
        h.write(
            "Diagnostic only: no production raster or radiance table was modified.\n\n"
        )

        h.write("EUROPE-WIDE BRIGHTNESS SUMMARY\n")
        h.write(
            europe_bins.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\n")

        h.write("ADDITIVE STEP BRIGHTNESS DEPENDENCE\n")
        h.write("All bins:\n")
        for key, value in additive_all.items():
            h.write(f"{key}: {value}\n")
        h.write(f"\nBaseline <= {cfg['empirical_max_baseline']}:\n")
        for key, value in additive_le5.items():
            h.write(f"{key}: {value}\n")

        h.write("\nMULTIPLICATIVE STEP BRIGHTNESS DEPENDENCE\n")
        h.write("All bins:\n")
        for key, value in multiplicative_all.items():
            h.write(f"{key}: {value}\n")
        h.write(f"\nBaseline <= {cfg['empirical_max_baseline']}:\n")
        for key, value in multiplicative_le5.items():
            h.write(f"{key}: {value}\n")

        h.write("\nEMPIRICAL ADDITIVE STEP\n")
        h.write(
            f"country-balanced median for baseline <= "
            f"{cfg['empirical_max_baseline']} = {empirical:.9g} "
            "nW cm-2 sr-1\n"
        )

        h.write("\nCANDIDATE CORRECTION PERFORMANCE\n")
        h.write(
            scenario_summary.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nINTERPRETATION GUIDE\n")
        h.write(
            "- If the additive step is approximately flat with independent "
            "brightness while the multiplicative percentage is not, the "
            "calibration discontinuity is better represented as additive.\n"
        )
        h.write(
            "- Agreement between the simple common-slope segmented estimate and "
            "the post-slope-break sensitivity argues that the inferred step is "
            "not just a change in secular lighting trend.\n"
        )
        h.write(
            "- The candidate ranking prioritises residual performance through "
            "baseline <=10 nW because pass 2 showed that this is where a coherent "
            "positive calibration step is identifiable.\n"
        )
        h.write(
            "- A tapering correction should be preferred over a universal "
            "constant only if it clearly reduces residual brightness structure "
            "without worsening the low-to-moderate brightness range.\n"
        )
        h.write(
            "- Do not promote any correction directly from this table. The "
            "winning candidate should next be applied to 2017-2024 rasters in a "
            "diagnostic extraction, followed by the NUTS3 rotation test.\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["annual_vnl_2017_correction_pass3"]

    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    years = list(range(first_year, last_year + 1))
    baseline_years = tuple(int(x) for x in cfg["brightness_reference_years"])
    if len(baseline_years) != 2:
        raise ValueError("brightness_reference_years must contain two years")

    rasters = discover_rasters(config, years)
    nuts_release = int(cfg["nuts_release"])
    nuts_path = nuts_path_for_release(config, nuts_release)
    nuts = gpd.read_file(nuts_path)

    id_field = config["nuts"]["id_field"]
    for col in [id_field, "CNTR_CODE", "geometry"]:
        if col not in nuts.columns:
            raise KeyError(f"NUTS file missing required column {col!r}")

    finite_edges = [float(x) for x in cfg["brightness_bin_edges"]]
    bin_edges = np.array([-np.inf] + finite_edges + [np.inf], dtype=float)
    bin_labels = list(cfg["brightness_bin_labels"])
    if len(bin_labels) != len(bin_edges) - 1:
        raise ValueError("brightness_bin_labels does not match bin edges")

    min_bin_pixels = int(cfg["min_region_bin_pixels"])
    min_country_regions = int(cfg["min_country_bin_regions"])
    min_europe_countries = int(cfg["min_europe_bin_countries"])
    empirical_max = float(cfg["empirical_max_baseline"])
    taper_start = float(cfg["taper_start"])
    taper_end = float(cfg["taper_end"])
    pass2_reference = float(cfg["pass2_reference_offset"])
    dpi = int(cfg["figure_dpi"])

    print("2017 VIIRS CALIBRATION PASS 3: SEGMENTED 2013-2019 TIME SERIES")
    print("Diagnostic only; production data will not be modified.")
    print(f"Years: {first_year}-{last_year}")
    print(f"Brightness baseline: mean {baseline_years[0]}-{baseline_years[1]}")

    rows: list[dict] = []

    with ExitStack() as stack:
        sources = {
            year: stack.enter_context(rasterio.open(rasters[year]))
            for year in years
        }
        validate_alignment(sources)

        ref_crs = sources[years[0]].crs
        if nuts.crs != ref_crs:
            nuts = nuts.to_crs(ref_crs)

        total = len(nuts)
        for idx, row in enumerate(nuts.itertuples(index=False), start=1):
            geom = getattr(row, "geometry")
            if geom is None or geom.is_empty:
                continue

            rows.extend(
                region_bin_models(
                    str(getattr(row, id_field)),
                    str(getattr(row, "CNTR_CODE")),
                    geom,
                    sources,
                    baseline_years,
                    bin_edges,
                    bin_labels,
                    min_bin_pixels,
                )
            )

            if idx % 50 == 0 or idx == total:
                print(
                    f"  processed {idx:,}/{total:,} NUTS3 regions",
                    flush=True,
                )

    region_bins = pd.DataFrame(rows)
    country_bins = aggregate_country_bins(
        region_bins,
        min_regions=min_country_regions,
    )
    europe_bins = aggregate_europe_bins(
        country_bins,
        bin_labels=bin_labels,
        min_countries=min_europe_countries,
    )

    additive_all = dependence(
        country_bins,
        "median_additive_step_nW",
    )
    additive_le5 = dependence(
        country_bins,
        "median_additive_step_nW",
        max_baseline=empirical_max,
    )
    multiplicative_all = dependence(
        country_bins,
        "median_multiplicative_step_pct",
    )
    multiplicative_le5 = dependence(
        country_bins,
        "median_multiplicative_step_pct",
        max_baseline=empirical_max,
    )

    empirical, empirical_by_country = empirical_step(
        country_bins,
        max_baseline=empirical_max,
    )

    detail, scenario_summary = evaluate_scenarios(
        country_bins,
        empirical=empirical,
        pass2=pass2_reference,
        taper_start=taper_start,
        taper_end=taper_end,
    )

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "annual_change"
        / "annual_vnl_2017_correction_pass3"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    paths = {
        "region_bins": processed / "annual_vnl_2017_pass3_nuts3_brightness.csv",
        "country_bins": processed / "annual_vnl_2017_pass3_country_brightness.csv",
        "europe_bins": processed / "annual_vnl_2017_pass3_brightness_summary.csv",
        "empirical_by_country": processed / "annual_vnl_2017_pass3_empirical_country.csv",
        "scenario_detail": processed / "annual_vnl_2017_pass3_scenario_detail.csv",
        "scenario_summary": processed / "annual_vnl_2017_pass3_scenario_summary.csv",
        "summary": processed / "annual_vnl_2017_correction_pass3_summary.txt",
    }

    region_bins.to_csv(paths["region_bins"], index=False)
    country_bins.to_csv(paths["country_bins"], index=False)
    europe_bins.to_csv(paths["europe_bins"], index=False)
    empirical_by_country.to_csv(paths["empirical_by_country"], index=False)
    detail.to_csv(paths["scenario_detail"], index=False)
    scenario_summary.to_csv(paths["scenario_summary"], index=False)

    fig_steps = figures / "segmented_additive_step_by_brightness.png"
    fig_form = figures / "additive_vs_multiplicative_step.png"
    fig_residuals = figures / "candidate_residuals_by_brightness.png"

    plot_steps(europe_bins, fig_steps, dpi)
    plot_additive_vs_multiplicative(europe_bins, fig_form, dpi)
    plot_candidate_residuals(detail, fig_residuals, dpi)

    write_summary(
        paths["summary"],
        europe_bins,
        additive_all,
        additive_le5,
        multiplicative_all,
        multiplicative_le5,
        empirical,
        scenario_summary,
        cfg,
    )

    print("\nEUROPE-WIDE SEGMENTED STEPS")
    print(
        europe_bins.to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )

    print("\nEMPIRICAL COUNTRY-BALANCED ADDITIVE STEP")
    print(
        f"  baseline <= {empirical_max}: {empirical:.6f} "
        "nW cm-2 sr-1"
    )

    print("\nCANDIDATE CORRECTION PERFORMANCE")
    print(
        scenario_summary.to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )

    print("\nWROTE")
    for path in [*paths.values(), fig_steps, fig_form, fig_residuals]:
        print(f"  {path}")


if __name__ == "__main__":
    main()
