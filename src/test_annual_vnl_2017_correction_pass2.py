#!/usr/bin/env python3
"""Second-pass diagnostic for the 2017 Annual VNL V2.2 calibration step.

This stage fixes the main weakness in the first retained-pixel diagnostic:
brightness strata are now defined from an INDEPENDENT pre-event reference,
the mean of 2013 and 2014 radiance, rather than from 2016 itself.

The primary analysis keeps only pixels that are positive in every year from
2013 through 2018. For each such pixel:

    baseline = mean(R2013, R2014)

    d15_16 = R2016 - R2015
    d16_17 = R2017 - R2016
    d17_18 = R2018 - R2017

    excess_2017 = d16_17 - 0.5 * (d15_16 + d17_18)

The independent baseline avoids mathematical coupling / regression-to-the-mean
created when 2016 radiance is used both to define a bin and in the adjacent
transition differences.

The stage then makes a second pass at the *form and magnitude* of a candidate
production correction. It compares four classes of pixel-level additive
correction, all applied only to positive post-2017 pixels and clipped at zero:

- no correction;
- literature sensitivity values configured in pipeline.yaml;
- an empirical Europe-wide value: median of country-level retained-pixel
  excesses;
- an empirical dim-to-moderate value: for each country take the median excess
  across independent brightness bins up to a configured maximum baseline
  radiance, then take the median across countries.

Each candidate is evaluated by applying it to 2017 and 2018 pixel radiances and
recomputing the residual 2017 excess. The best constant correction is the one
that leaves country residuals centred on zero with the smallest country-level
absolute residual while not creating a strong brightness-dependent residual.

IMPORTANT
---------
This is diagnostic only. It does not overwrite any production raster or
radiance table. A production correction should be implemented only after the
independent-brightness and candidate-residual results have been inspected.
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
from scipy.stats import binomtest, spearmanr, wilcoxon
from shapely.geometry import mapping

from extract_nuts3_radiance import discover_rasters, nuts_path_for_release
from project_config import configured_path, load_config


ANALYSIS_YEARS = (2013, 2014, 2015, 2016, 2017, 2018)


def median(values: np.ndarray) -> float:
    return float(np.nanmedian(values)) if values.size else np.nan


def q(values: np.ndarray, prob: float) -> float:
    return float(np.nanquantile(values, prob)) if values.size else np.nan


def validate_alignment(sources: dict[int, rasterio.io.DatasetReader]) -> None:
    ref = sources[ANALYSIS_YEARS[0]]
    for year in ANALYSIS_YEARS[1:]:
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
                f"{ANALYSIS_YEARS[0]}: " + "; ".join(problems)
            )


def read_region(
    geom,
    sources: dict[int, rasterio.io.DatasetReader],
) -> tuple[dict[int, np.ndarray], np.ndarray] | None:
    ref = sources[ANALYSIS_YEARS[0]]
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


def persistent_support(
    arrays: dict[int, np.ndarray],
    inside: np.ndarray,
) -> np.ndarray:
    support = inside.copy()
    for year in ANALYSIS_YEARS:
        support &= np.isfinite(arrays[year])
        support &= arrays[year] > 0
    return support


def pixel_components(
    arrays: dict[int, np.ndarray],
    support: np.ndarray,
) -> dict[str, np.ndarray]:
    r13 = arrays[2013][support].astype(np.float64, copy=False)
    r14 = arrays[2014][support].astype(np.float64, copy=False)
    r15 = arrays[2015][support].astype(np.float64, copy=False)
    r16 = arrays[2016][support].astype(np.float64, copy=False)
    r17 = arrays[2017][support].astype(np.float64, copy=False)
    r18 = arrays[2018][support].astype(np.float64, copy=False)

    baseline = 0.5 * (r13 + r14)
    d15 = r16 - r15
    d16 = r17 - r16
    d17 = r18 - r17
    excess = d16 - 0.5 * (d15 + d17)

    return {
        "baseline": baseline,
        "r15": r15,
        "r16": r16,
        "r17": r17,
        "r18": r18,
        "d15": d15,
        "d16": d16,
        "d17": d17,
        "excess": excess,
    }


def first_pass_region(
    nuts_id: str,
    country: str,
    geom,
    sources: dict[int, rasterio.io.DatasetReader],
    bin_edges: np.ndarray,
    bin_labels: list[str],
    min_bin_pixels: int,
) -> tuple[dict, list[dict]]:
    read = read_region(geom, sources)
    if read is None:
        return {"NUTS_ID": nuts_id, "CNTR_CODE": country, "read_ok": False}, []

    arrays, inside = read
    support = persistent_support(arrays, inside)
    n = int(support.sum())

    row: dict[str, float | int | str | bool] = {
        "NUTS_ID": nuts_id,
        "CNTR_CODE": country,
        "read_ok": True,
        "n_persistent_2013_2018_pixels": n,
    }

    if n == 0:
        return row, []

    z = pixel_components(arrays, support)
    row.update(
        {
            "median_baseline_2013_2014": median(z["baseline"]),
            "median_d15_16": median(z["d15"]),
            "median_d16_17": median(z["d16"]),
            "median_d17_18": median(z["d17"]),
            "candidate_offset": (
                median(z["d16"])
                - 0.5 * (median(z["d15"]) + median(z["d17"]))
            ),
            "median_pixel_excess_2017": median(z["excess"]),
            "pixel_excess_q10": q(z["excess"], 0.10),
            "pixel_excess_q90": q(z["excess"], 0.90),
            "fraction_positive_excess": float(np.mean(z["excess"] > 0)),
        }
    )

    bin_index = np.digitize(
        z["baseline"],
        bin_edges[1:-1],
        right=False,
    )

    bin_rows: list[dict] = []
    for i, label in enumerate(bin_labels):
        take = bin_index == i
        n_bin = int(take.sum())
        if n_bin < min_bin_pixels:
            continue

        d15 = z["d15"][take]
        d16 = z["d16"][take]
        d17 = z["d17"][take]
        ex = z["excess"][take]

        med15 = median(d15)
        med16 = median(d16)
        med17 = median(d17)

        bin_rows.append(
            {
                "NUTS_ID": nuts_id,
                "CNTR_CODE": country,
                "brightness_bin": label,
                "n_pixels": n_bin,
                "median_baseline_2013_2014": median(z["baseline"][take]),
                "median_d15_16": med15,
                "median_d16_17": med16,
                "median_d17_18": med17,
                "candidate_offset": med16 - 0.5 * (med15 + med17),
                "median_pixel_excess_2017": median(ex),
                "fraction_positive_excess": float(np.mean(ex > 0)),
            }
        )

    return row, bin_rows


def aggregate_country_region_offsets(
    regions: pd.DataFrame,
    min_region_pixels: int,
    min_country_regions: int,
) -> pd.DataFrame:
    x = regions[
        regions["read_ok"].fillna(False)
        & (
            regions["n_persistent_2013_2018_pixels"].fillna(0)
            >= min_region_pixels
        )
    ].copy()

    rows = []
    for country, g in x.groupby("CNTR_CODE"):
        if len(g) < min_country_regions:
            continue
        rows.append(
            {
                "CNTR_CODE": country,
                "n_regions": int(len(g)),
                "n_pixels": int(g["n_persistent_2013_2018_pixels"].sum()),
                "median_region_baseline": float(
                    g["median_baseline_2013_2014"].median()
                ),
                "median_region_d15_16": float(g["median_d15_16"].median()),
                "median_region_d16_17": float(g["median_d16_17"].median()),
                "median_region_d17_18": float(g["median_d17_18"].median()),
                "median_region_candidate_offset": float(
                    g["candidate_offset"].median()
                ),
                "median_region_pixel_excess": float(
                    g["median_pixel_excess_2017"].median()
                ),
            }
        )

    return pd.DataFrame(rows).sort_values("median_region_candidate_offset")


def aggregate_country_bins(
    bin_regions: pd.DataFrame,
    min_country_bin_regions: int,
) -> pd.DataFrame:
    rows = []
    for (country, brightness_bin), g in bin_regions.groupby(
        ["CNTR_CODE", "brightness_bin"],
        observed=True,
    ):
        if len(g) < min_country_bin_regions:
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
                "median_d15_16": float(g["median_d15_16"].median()),
                "median_d16_17": float(g["median_d16_17"].median()),
                "median_d17_18": float(g["median_d17_18"].median()),
                "candidate_offset": float(g["candidate_offset"].median()),
                "median_pixel_excess_2017": float(
                    g["median_pixel_excess_2017"].median()
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
                "median_country_candidate_offset": float(
                    g["candidate_offset"].median()
                ),
                "candidate_offset_q25": float(
                    g["candidate_offset"].quantile(0.25)
                ),
                "candidate_offset_q75": float(
                    g["candidate_offset"].quantile(0.75)
                ),
                "median_country_pixel_excess": float(
                    g["median_pixel_excess_2017"].median()
                ),
            }
        )

    return pd.DataFrame(rows)


def brightness_dependence(
    country_bins: pd.DataFrame,
    max_baseline: float | None = None,
) -> dict[str, float]:
    x = country_bins.dropna(
        subset=["median_baseline_2013_2014", "candidate_offset"]
    ).copy()

    if max_baseline is not None:
        x = x[x["median_baseline_2013_2014"] <= max_baseline].copy()

    x = x[x["median_baseline_2013_2014"] > 0].copy()

    if len(x) < 3:
        return {
            "n_country_bin_cells": int(len(x)),
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
        }

    rho, p = spearmanr(
        np.log10(x["median_baseline_2013_2014"]),
        x["candidate_offset"],
    )
    return {
        "n_country_bin_cells": int(len(x)),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
    }


def empirical_candidates(
    country: pd.DataFrame,
    country_bins: pd.DataFrame,
    max_baseline: float,
) -> dict[str, float]:
    global_country_median = float(
        country["median_region_candidate_offset"].median()
    )

    dim = country_bins[
        country_bins["median_baseline_2013_2014"] <= max_baseline
    ].copy()

    per_country_dim = (
        dim.groupby("CNTR_CODE", as_index=False)["candidate_offset"]
        .median()
        .rename(columns={"candidate_offset": "country_dim_offset"})
    )
    dim_country_median = float(
        per_country_dim["country_dim_offset"].median()
    )

    return {
        "empirical_country_median": global_country_median,
        "empirical_dim_country_median": dim_country_median,
    }


def build_scenarios(
    cfg: dict,
    empirical: dict[str, float],
) -> pd.DataFrame:
    rows = [
        {
            "scenario": "none",
            "offset_nW_cm2_sr": 0.0,
            "source": "no correction",
        }
    ]

    for label, value in cfg["literature_pixel_offsets"].items():
        rows.append(
            {
                "scenario": str(label),
                "offset_nW_cm2_sr": float(value),
                "source": "literature sensitivity",
            }
        )

    for label, value in empirical.items():
        rows.append(
            {
                "scenario": label,
                "offset_nW_cm2_sr": float(value),
                "source": "empirical independent-brightness diagnostic",
            }
        )

    return pd.DataFrame(rows)


def corrected_excess(
    z: dict[str, np.ndarray],
    offset: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply post-2017 positive-pixel subtraction with clipping at zero."""
    r17c = np.maximum(z["r17"] - offset, 0.0)
    r18c = np.maximum(z["r18"] - offset, 0.0)

    d16c = r17c - z["r16"]
    d17c = r18c - r17c
    residual = d16c - 0.5 * (z["d15"] + d17c)

    clipped17 = z["r17"] <= offset
    clipped18 = z["r18"] <= offset
    return residual, clipped17, clipped18


def second_pass_region(
    nuts_id: str,
    country: str,
    geom,
    sources: dict[int, rasterio.io.DatasetReader],
    scenarios: pd.DataFrame,
    bin_edges: np.ndarray,
    bin_labels: list[str],
    min_bin_pixels: int,
) -> tuple[list[dict], list[dict]]:
    read = read_region(geom, sources)
    if read is None:
        return [], []

    arrays, inside = read
    support = persistent_support(arrays, inside)
    n = int(support.sum())
    if n == 0:
        return [], []

    z = pixel_components(arrays, support)
    bin_index = np.digitize(
        z["baseline"],
        bin_edges[1:-1],
        right=False,
    )

    region_rows: list[dict] = []
    bin_rows: list[dict] = []

    for scenario in scenarios.itertuples(index=False):
        offset = float(scenario.offset_nW_cm2_sr)
        residual, clipped17, clipped18 = corrected_excess(z, offset)

        region_rows.append(
            {
                "scenario": scenario.scenario,
                "offset_nW_cm2_sr": offset,
                "NUTS_ID": nuts_id,
                "CNTR_CODE": country,
                "n_pixels": n,
                "median_residual_excess": median(residual),
                "residual_q10": q(residual, 0.10),
                "residual_q90": q(residual, 0.90),
                "fraction_positive_residual": float(
                    np.mean(residual > 0)
                ),
                "fraction_clipped_2017": float(np.mean(clipped17)),
                "fraction_clipped_2018": float(np.mean(clipped18)),
            }
        )

        for i, label in enumerate(bin_labels):
            take = bin_index == i
            n_bin = int(take.sum())
            if n_bin < min_bin_pixels:
                continue

            bin_rows.append(
                {
                    "scenario": scenario.scenario,
                    "offset_nW_cm2_sr": offset,
                    "NUTS_ID": nuts_id,
                    "CNTR_CODE": country,
                    "brightness_bin": label,
                    "n_pixels": n_bin,
                    "median_baseline_2013_2014": median(
                        z["baseline"][take]
                    ),
                    "median_residual_excess": median(residual[take]),
                    "fraction_positive_residual": float(
                        np.mean(residual[take] > 0)
                    ),
                    "fraction_clipped_2017": float(
                        np.mean(clipped17[take])
                    ),
                    "fraction_clipped_2018": float(
                        np.mean(clipped18[take])
                    ),
                }
            )

    return region_rows, bin_rows


def aggregate_scenario_countries(
    region_scenarios: pd.DataFrame,
    min_region_pixels: int,
    min_country_regions: int,
) -> pd.DataFrame:
    x = region_scenarios[
        region_scenarios["n_pixels"] >= min_region_pixels
    ].copy()

    rows = []
    for (scenario, offset, country), g in x.groupby(
        ["scenario", "offset_nW_cm2_sr", "CNTR_CODE"],
        observed=True,
    ):
        if len(g) < min_country_regions:
            continue

        rows.append(
            {
                "scenario": scenario,
                "offset_nW_cm2_sr": float(offset),
                "CNTR_CODE": country,
                "n_regions": int(len(g)),
                "n_pixels": int(g["n_pixels"].sum()),
                "median_region_residual_excess": float(
                    g["median_residual_excess"].median()
                ),
                "median_region_abs_residual": float(
                    g["median_residual_excess"].abs().median()
                ),
                "median_region_fraction_clipped_2017": float(
                    g["fraction_clipped_2017"].median()
                ),
                "median_region_fraction_clipped_2018": float(
                    g["fraction_clipped_2018"].median()
                ),
            }
        )

    return pd.DataFrame(rows)


def aggregate_scenario_summary(
    country_scenarios: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for (scenario, offset), g in country_scenarios.groupby(
        ["scenario", "offset_nW_cm2_sr"],
        observed=True,
    ):
        values = g["median_region_residual_excess"].to_numpy(dtype=float)
        n = len(values)
        n_positive = int(np.sum(values > 0))

        if n:
            sign_p = float(
                binomtest(
                    n_positive,
                    n,
                    p=0.5,
                    alternative="two-sided",
                ).pvalue
            )
            try:
                wilcox_p = float(
                    wilcoxon(
                        values,
                        alternative="two-sided",
                        zero_method="wilcox",
                    ).pvalue
                )
            except ValueError:
                wilcox_p = np.nan
        else:
            sign_p = np.nan
            wilcox_p = np.nan

        rows.append(
            {
                "scenario": scenario,
                "offset_nW_cm2_sr": float(offset),
                "n_countries": int(n),
                "median_country_residual": (
                    float(np.median(values)) if n else np.nan
                ),
                "median_abs_country_residual": (
                    float(np.median(np.abs(values))) if n else np.nan
                ),
                "country_residual_q25": (
                    float(np.quantile(values, 0.25)) if n else np.nan
                ),
                "country_residual_q75": (
                    float(np.quantile(values, 0.75)) if n else np.nan
                ),
                "fraction_positive_country_residual": (
                    n_positive / n if n else np.nan
                ),
                "country_sign_two_sided_p": sign_p,
                "country_wilcoxon_two_sided_p": wilcox_p,
                "median_country_fraction_clipped_2017": float(
                    g["median_region_fraction_clipped_2017"].median()
                ),
                "median_country_fraction_clipped_2018": float(
                    g["median_region_fraction_clipped_2018"].median()
                ),
            }
        )

    out = pd.DataFrame(rows)
    out["rank_abs_country_residual"] = out[
        "median_abs_country_residual"
    ].rank(method="min")
    return out.sort_values(
        ["rank_abs_country_residual", "median_abs_country_residual"]
    )


def aggregate_scenario_brightness(
    bin_scenarios: pd.DataFrame,
    min_country_bin_regions: int,
    min_countries: int,
    bin_labels: list[str],
) -> pd.DataFrame:
    country_rows = []

    for (
        scenario,
        offset,
        country,
        brightness_bin,
    ), g in bin_scenarios.groupby(
        [
            "scenario",
            "offset_nW_cm2_sr",
            "CNTR_CODE",
            "brightness_bin",
        ],
        observed=True,
    ):
        if len(g) < min_country_bin_regions:
            continue

        country_rows.append(
            {
                "scenario": scenario,
                "offset_nW_cm2_sr": float(offset),
                "CNTR_CODE": country,
                "brightness_bin": brightness_bin,
                "n_regions": int(len(g)),
                "median_baseline_2013_2014": float(
                    g["median_baseline_2013_2014"].median()
                ),
                "median_region_residual_excess": float(
                    g["median_residual_excess"].median()
                ),
                "median_fraction_clipped_2017": float(
                    g["fraction_clipped_2017"].median()
                ),
            }
        )

    country_bins = pd.DataFrame(country_rows)
    europe_rows = []

    for scenario in country_bins["scenario"].drop_duplicates():
        s = country_bins[country_bins["scenario"] == scenario]
        offset = float(s["offset_nW_cm2_sr"].iloc[0])

        for label in bin_labels:
            g = s[s["brightness_bin"] == label]
            if len(g) < min_countries:
                continue

            europe_rows.append(
                {
                    "scenario": scenario,
                    "offset_nW_cm2_sr": offset,
                    "brightness_bin": label,
                    "n_countries": int(len(g)),
                    "median_country_baseline": float(
                        g["median_baseline_2013_2014"].median()
                    ),
                    "median_country_residual_excess": float(
                        g["median_region_residual_excess"].median()
                    ),
                    "residual_q25": float(
                        g["median_region_residual_excess"].quantile(0.25)
                    ),
                    "residual_q75": float(
                        g["median_region_residual_excess"].quantile(0.75)
                    ),
                    "median_country_fraction_clipped_2017": float(
                        g["median_fraction_clipped_2017"].median()
                    ),
                }
            )

    return country_bins, pd.DataFrame(europe_rows)


def plot_independent_brightness(
    europe_bins: pd.DataFrame,
    output: Path,
    dpi: int,
) -> None:
    if europe_bins.empty:
        return

    b = europe_bins.sort_values("median_country_baseline").copy()

    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    ax.plot(
        b["median_country_baseline"],
        b["median_country_candidate_offset"],
        marker="o",
        linewidth=1.8,
    )
    ax.fill_between(
        b["median_country_baseline"],
        b["candidate_offset_q25"],
        b["candidate_offset_q75"],
        alpha=0.18,
        linewidth=0,
    )
    ax.axhline(0, color="0.4", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel(
        "Independent baseline radiance: mean 2013-2014 "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_ylabel(
        "2017 excess contrast (nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title(
        "2017 retained-pixel excess by independent pre-event brightness"
    )
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_scenario_summary(
    summary: pd.DataFrame,
    output: Path,
    dpi: int,
) -> None:
    s = summary.sort_values("offset_nW_cm2_sr").copy()
    x = np.arange(len(s))

    fig, ax = plt.subplots(figsize=(10.5, 6.3))
    ax.errorbar(
        x,
        s["median_country_residual"],
        yerr=np.vstack(
            [
                s["median_country_residual"]
                - s["country_residual_q25"],
                s["country_residual_q75"]
                - s["median_country_residual"],
            ]
        ),
        fmt="o",
        capsize=4,
    )
    ax.axhline(0, color="0.4", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [
            f"{name}\n({offset:.3f})"
            for name, offset in zip(
                s["scenario"],
                s["offset_nW_cm2_sr"],
            )
        ],
        rotation=25,
        ha="right",
    )
    ax.set_ylabel(
        "Residual country 2017 excess (nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title(
        "Candidate post-2017 pixel corrections: residual country discontinuity"
    )
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_scenario_brightness(
    scenario_bins: pd.DataFrame,
    output: Path,
    dpi: int,
) -> None:
    if scenario_bins.empty:
        return

    fig, ax = plt.subplots(figsize=(10.2, 6.5))
    for scenario, g in scenario_bins.groupby("scenario", sort=False):
        g = g.sort_values("median_country_baseline")
        ax.plot(
            g["median_country_baseline"],
            g["median_country_residual_excess"],
            marker="o",
            linewidth=1.5,
            label=scenario,
        )

    ax.axhline(0, color="0.4", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel(
        "Independent baseline radiance: mean 2013-2014 "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_ylabel(
        "Residual 2017 excess after correction "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title(
        "Does a constant pixel correction remove the step across brightness?"
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    path: Path,
    country: pd.DataFrame,
    europe_bins: pd.DataFrame,
    dependence_all: dict[str, float],
    dependence_dim: dict[str, float],
    scenarios: pd.DataFrame,
    scenario_summary: pd.DataFrame,
    scenario_bins: pd.DataFrame,
    cfg: dict,
) -> None:
    with path.open("w", encoding="utf-8") as h:
        h.write(
            "ANNUAL VNL V2.2 2017 CALIBRATION: INDEPENDENT-BRIGHTNESS PASS 2\n"
        )
        h.write("=" * 78 + "\n\n")
        h.write(
            "Diagnostic only: no production raster or radiance table was modified.\n\n"
        )

        h.write("WHY THIS PASS\n")
        h.write(
            "Brightness is defined from mean 2013-2014 radiance, independent of "
            "the 2015-2018 transition contrast. Pixels must be positive in all "
            "six years 2013-2018. This removes the mathematical coupling present "
            "when 2016 radiance was used to define brightness strata.\n\n"
        )

        h.write("COUNTRY OFFSETS ON SIX-YEAR PERSISTENT-LIT SUPPORT\n")
        h.write(
            country.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\n")

        h.write("INDEPENDENT BRIGHTNESS DEPENDENCE\n")
        h.write(
            europe_bins.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nAll brightness bins:\n")
        for key, value in dependence_all.items():
            h.write(f"{key}: {value}\n")
        h.write(
            f"\nBrightness <= {cfg['candidate_brightness_max']} "
            "nW cm-2 sr-1:\n"
        )
        for key, value in dependence_dim.items():
            h.write(f"{key}: {value}\n")

        h.write("\nCANDIDATE PIXEL CORRECTIONS\n")
        h.write(
            scenarios.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nCANDIDATE PERFORMANCE: COUNTRY RESIDUALS\n")
        h.write(
            scenario_summary.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nCANDIDATE PERFORMANCE: BRIGHTNESS RESIDUALS\n")
        h.write(
            scenario_bins.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nDECISION GUIDE\n")
        h.write(
            "- First ask whether the independent brightness curve is roughly "
            "flat through the dim-to-moderate range. If yes, one additive pixel "
            "offset is plausible.\n"
        )
        h.write(
            "- Do not choose an offset merely because its median residual is "
            "zero by construction. Compare median absolute country residual, "
            "country IQR, clipping fraction and residual brightness structure.\n"
        )
        h.write(
            "- Literature 0.125 and 0.15 values are sensitivity candidates, not "
            "assumed truths for Annual VNL V2.2 average-masked.\n"
        )
        h.write(
            "- If one constant pixel correction performs well, the next step is "
            "to apply it diagnostically to 2017+ rasters and re-extract NUTS3 "
            "means with the production exact-overlap area weighting before "
            "promoting it to production.\n"
        )
        h.write(
            "- If residuals retain strong brightness dependence, do not force a "
            "single global additive correction; investigate a brightness-aware "
            "or upstream-product correction instead.\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["annual_vnl_2017_correction_pass2"]

    years = tuple(int(x) for x in cfg["years"])
    if years != ANALYSIS_YEARS:
        raise ValueError(
            f"This diagnostic expects years {ANALYSIS_YEARS}; found {years}"
        )

    rasters = discover_rasters(config, list(ANALYSIS_YEARS))
    nuts_release = int(cfg["nuts_release"])
    nuts_path = nuts_path_for_release(config, nuts_release)

    nuts = gpd.read_file(nuts_path)
    if nuts.crs is None:
        raise ValueError(f"NUTS file has no CRS: {nuts_path}")

    id_field = config["nuts"]["id_field"]
    for col in [id_field, "CNTR_CODE", "geometry"]:
        if col not in nuts.columns:
            raise KeyError(f"NUTS file missing required column {col!r}")

    finite_edges = [float(x) for x in cfg["brightness_bin_edges"]]
    bin_edges = np.array([-np.inf] + finite_edges + [np.inf], dtype=float)
    bin_labels = list(cfg["brightness_bin_labels"])
    if len(bin_labels) != len(bin_edges) - 1:
        raise ValueError("brightness_bin_labels does not match bin edges")

    min_region_pixels = int(cfg["min_region_persistent_pixels"])
    min_country_regions = int(cfg["min_country_regions"])
    min_bin_pixels = int(cfg["min_region_bin_pixels"])
    min_country_bin_regions = int(cfg["min_country_bin_regions"])
    min_europe_bin_countries = int(cfg["min_europe_bin_countries"])
    candidate_max = float(cfg["candidate_brightness_max"])
    dpi = int(cfg["figure_dpi"])

    print("ANNUAL VNL V2.2 2017 CALIBRATION: INDEPENDENT-BRIGHTNESS PASS 2")
    print("Diagnostic only; production data will not be modified.")
    print("Brightness reference: mean 2013-2014")
    for year in ANALYSIS_YEARS:
        print(f"  {year}: {rasters[year]}")

    first_region_rows: list[dict] = []
    first_bin_rows: list[dict] = []

    with ExitStack() as stack:
        sources = {
            year: stack.enter_context(rasterio.open(rasters[year]))
            for year in ANALYSIS_YEARS
        }
        validate_alignment(sources)

        raster_crs = sources[ANALYSIS_YEARS[0]].crs
        if nuts.crs != raster_crs:
            nuts = nuts.to_crs(raster_crs)

        total = len(nuts)
        for idx, row in enumerate(nuts.itertuples(index=False), start=1):
            geom = getattr(row, "geometry")
            if geom is None or geom.is_empty:
                continue

            region, bins = first_pass_region(
                str(getattr(row, id_field)),
                str(getattr(row, "CNTR_CODE")),
                geom,
                sources,
                bin_edges,
                bin_labels,
                min_bin_pixels,
            )
            first_region_rows.append(region)
            first_bin_rows.extend(bins)

            if idx % 50 == 0 or idx == total:
                print(
                    f"  first pass: {idx:,}/{total:,} NUTS3 regions",
                    flush=True,
                )

    regions = pd.DataFrame(first_region_rows)
    bin_regions = pd.DataFrame(first_bin_rows)

    country = aggregate_country_region_offsets(
        regions,
        min_region_pixels=min_region_pixels,
        min_country_regions=min_country_regions,
    )
    country_bins = aggregate_country_bins(
        bin_regions,
        min_country_bin_regions=min_country_bin_regions,
    )
    europe_bins = aggregate_europe_bins(
        country_bins,
        bin_labels=bin_labels,
        min_countries=min_europe_bin_countries,
    )

    dependence_all = brightness_dependence(country_bins)
    dependence_dim = brightness_dependence(
        country_bins,
        max_baseline=candidate_max,
    )

    empirical = empirical_candidates(
        country,
        country_bins,
        max_baseline=candidate_max,
    )
    scenarios = build_scenarios(cfg, empirical)

    print("\nCANDIDATE OFFSETS")
    print(
        scenarios.to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )

    second_region_rows: list[dict] = []
    second_bin_rows: list[dict] = []

    with ExitStack() as stack:
        sources = {
            year: stack.enter_context(rasterio.open(rasters[year]))
            for year in ANALYSIS_YEARS
        }
        validate_alignment(sources)

        total = len(nuts)
        for idx, row in enumerate(nuts.itertuples(index=False), start=1):
            geom = getattr(row, "geometry")
            if geom is None or geom.is_empty:
                continue

            region_rows, bin_rows = second_pass_region(
                str(getattr(row, id_field)),
                str(getattr(row, "CNTR_CODE")),
                geom,
                sources,
                scenarios,
                bin_edges,
                bin_labels,
                min_bin_pixels,
            )
            second_region_rows.extend(region_rows)
            second_bin_rows.extend(bin_rows)

            if idx % 50 == 0 or idx == total:
                print(
                    f"  correction pass: {idx:,}/{total:,} NUTS3 regions",
                    flush=True,
                )

    region_scenarios = pd.DataFrame(second_region_rows)
    bin_scenarios = pd.DataFrame(second_bin_rows)

    country_scenarios = aggregate_scenario_countries(
        region_scenarios,
        min_region_pixels=min_region_pixels,
        min_country_regions=min_country_regions,
    )
    scenario_summary = aggregate_scenario_summary(country_scenarios)
    scenario_country_bins, scenario_bins = aggregate_scenario_brightness(
        bin_scenarios,
        min_country_bin_regions=min_country_bin_regions,
        min_countries=min_europe_bin_countries,
        bin_labels=bin_labels,
    )

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "annual_change"
        / "annual_vnl_2017_correction_pass2"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    paths = {
        "regions": processed / "annual_vnl_2017_pass2_nuts3.csv",
        "region_bins": processed / "annual_vnl_2017_pass2_nuts3_brightness.csv",
        "country": processed / "annual_vnl_2017_pass2_country.csv",
        "country_bins": processed / "annual_vnl_2017_pass2_country_brightness.csv",
        "europe_bins": processed / "annual_vnl_2017_pass2_brightness_summary.csv",
        "scenarios": processed / "annual_vnl_2017_pass2_candidate_offsets.csv",
        "region_scenarios": processed / "annual_vnl_2017_pass2_scenario_nuts3.csv",
        "country_scenarios": processed / "annual_vnl_2017_pass2_scenario_country.csv",
        "scenario_summary": processed / "annual_vnl_2017_pass2_scenario_summary.csv",
        "scenario_country_bins": processed / "annual_vnl_2017_pass2_scenario_country_brightness.csv",
        "scenario_bins": processed / "annual_vnl_2017_pass2_scenario_brightness.csv",
        "summary": processed / "annual_vnl_2017_correction_pass2_summary.txt",
    }

    regions.to_csv(paths["regions"], index=False)
    bin_regions.to_csv(paths["region_bins"], index=False)
    country.to_csv(paths["country"], index=False)
    country_bins.to_csv(paths["country_bins"], index=False)
    europe_bins.to_csv(paths["europe_bins"], index=False)
    scenarios.to_csv(paths["scenarios"], index=False)
    region_scenarios.to_csv(paths["region_scenarios"], index=False)
    country_scenarios.to_csv(paths["country_scenarios"], index=False)
    scenario_summary.to_csv(paths["scenario_summary"], index=False)
    scenario_country_bins.to_csv(paths["scenario_country_bins"], index=False)
    scenario_bins.to_csv(paths["scenario_bins"], index=False)

    fig_independent = figures / "independent_brightness_dependence.png"
    fig_scenarios = figures / "candidate_country_residuals.png"
    fig_scenario_bins = figures / "candidate_brightness_residuals.png"

    plot_independent_brightness(europe_bins, fig_independent, dpi)
    plot_scenario_summary(scenario_summary, fig_scenarios, dpi)
    plot_scenario_brightness(scenario_bins, fig_scenario_bins, dpi)

    write_summary(
        paths["summary"],
        country,
        europe_bins,
        dependence_all,
        dependence_dim,
        scenarios,
        scenario_summary,
        scenario_bins,
        cfg,
    )

    print("\nINDEPENDENT BRIGHTNESS DEPENDENCE")
    print("All bins:")
    for key, value in dependence_all.items():
        print(f"  {key}: {value}")
    print(f"Baseline <= {candidate_max}:")
    for key, value in dependence_dim.items():
        print(f"  {key}: {value}")

    print("\nCANDIDATE PERFORMANCE")
    print(
        scenario_summary.to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )

    print("\nWROTE")
    for path in [*paths.values(), fig_independent, fig_scenarios, fig_scenario_bins]:
        print(f"  {path}")


if __name__ == "__main__":
    main()
