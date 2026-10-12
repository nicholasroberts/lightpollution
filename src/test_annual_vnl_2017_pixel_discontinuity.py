#!/usr/bin/env python3
"""Test whether Annual VNL V2.2 retains a 2017 pixel-radiance discontinuity.

This is a DIAGNOSTIC ONLY. It does not alter production radiance.

Why this test exists
--------------------
EOG Annual VNL V2 average-masked rasters set background pixels to zero. That
raises a specific question: if the background mask already removes the changed
2017 DNB noise floor from dark pixels, does a calibration step remain in the
radiances of pixels that are persistently classified as lit?

The analysis therefore separates two possible manifestations of the 2017
change:

1. MASK/SUPPORT CHANGE
   Did the fraction of positive (retained/lit) pixels change unusually from
   2016 to 2017 compared with the flanking 2015->2016 and 2017->2018
   transitions?

2. RETAINED-PIXEL RADIANCE STEP
   Among exactly the same pixels that are positive in ALL FOUR years
   2015-2018, is the 2016->2017 radiance change unusually positive?

The second test is the key one. It cannot be explained by pixels merely moving
from zero/background into the lit mask because every tested pixel is positive
in every year.

The script works on fixed NUTS2016 polygons. Each NUTS3 is processed directly
from the four annual rasters, then country summaries are formed from medians of
NUTS3 statistics so that very large countries/regions do not dominate the
diagnostic.

For a persistent-lit pixel p, define

    d15_16 = R2016 - R2015
    d16_17 = R2017 - R2016
    d17_18 = R2018 - R2017

and the local 2017 excess contrast

    excess_2017 =
        d16_17 - 0.5 * (d15_16 + d17_18)

A simple calibration step predicts a positive excess that is broadly
consistent across countries and, in absolute linear-radiance units, reasonably
stable across the 2016 brightness range.

The script also reconnects this pixel test to the previously observed
country-to-country jump in NUTS3 mean radiance by calculating the same
2015-2018 excess contrast from the existing RAW NUTS3 zonal means.

No correction is applied here. The purpose is to determine whether the Annual
VNL V2.2 average-masked product itself still contains a material 2017 step and,
if so, how that step behaves before choosing a production correction.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from pathlib import Path
from typing import Iterable

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from rasterio.windows import Window
from rasterio.features import geometry_window
from scipy.stats import binomtest, spearmanr, wilcoxon
from shapely.geometry import mapping

from extract_nuts3_radiance import discover_rasters, nuts_path_for_release
from project_config import configured_path, load_config


YEARS = (2015, 2016, 2017, 2018)
TRANSITIONS = (
    ("2015-2016", 2015, 2016),
    ("2016-2017", 2016, 2017),
    ("2017-2018", 2017, 2018),
)


def _float_or_nan(value) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return np.nan
    return value if np.isfinite(value) else np.nan


def _median(values: np.ndarray) -> float:
    if values.size == 0:
        return np.nan
    return float(np.nanmedian(values))


def _quantile(values: np.ndarray, q: float) -> float:
    if values.size == 0:
        return np.nan
    return float(np.nanquantile(values, q))


def _window_shape(window: Window) -> tuple[int, int]:
    return int(window.height), int(window.width)


def validate_alignment(sources: dict[int, rasterio.io.DatasetReader]) -> None:
    ref = sources[YEARS[0]]
    for year in YEARS[1:]:
        src = sources[year]
        problems = []
        if src.crs != ref.crs:
            problems.append(f"CRS {src.crs} != {ref.crs}")
        if src.width != ref.width or src.height != ref.height:
            problems.append(
                f"shape {src.width}x{src.height} != "
                f"{ref.width}x{ref.height}"
            )
        if not src.transform.almost_equals(ref.transform):
            problems.append("affine transform differs")
        if problems:
            raise ValueError(
                f"VIIRS raster grid for {year} is not aligned with "
                f"{YEARS[0]}: " + "; ".join(problems)
            )


def read_region_arrays(
    geom,
    sources: dict[int, rasterio.io.DatasetReader],
) -> tuple[dict[int, np.ndarray], np.ndarray] | None:
    """Read one NUTS3 bounding window and return four arrays + polygon mask."""
    ref = sources[YEARS[0]]
    try:
        window = geometry_window(ref, [mapping(geom)])
    except Exception:
        return None

    if window.width <= 0 or window.height <= 0:
        return None

    shape = _window_shape(window)
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
        arr = src.read(
            1,
            window=window,
            masked=True,
        )
        arrays[year] = arr.filled(np.nan).astype(np.float32, copy=False)

    return arrays, inside


def region_pixel_statistics(
    nuts_id: str,
    country: str,
    geom,
    sources: dict[int, rasterio.io.DatasetReader],
    brightness_edges: np.ndarray,
    brightness_labels: list[str],
    min_bin_pixels: int,
) -> tuple[dict, list[dict]]:
    read = read_region_arrays(geom, sources)
    if read is None:
        return (
            {
                "NUTS_ID": nuts_id,
                "CNTR_CODE": country,
                "read_ok": False,
            },
            [],
        )

    arrays, inside = read
    common_valid = inside.copy()
    for year in YEARS:
        common_valid &= np.isfinite(arrays[year])

    n_valid = int(common_valid.sum())
    if n_valid == 0:
        return (
            {
                "NUTS_ID": nuts_id,
                "CNTR_CODE": country,
                "read_ok": True,
                "n_common_valid_pixels": 0,
            },
            [],
        )

    lit = {
        year: common_valid & (arrays[year] > 0)
        for year in YEARS
    }
    persistent = common_valid.copy()
    for year in YEARS:
        persistent &= arrays[year] > 0

    n_persistent = int(persistent.sum())

    row: dict[str, float | int | str | bool] = {
        "NUTS_ID": nuts_id,
        "CNTR_CODE": country,
        "read_ok": True,
        "n_common_valid_pixels": n_valid,
        "n_persistent_lit_pixels": n_persistent,
    }

    for year in YEARS:
        row[f"lit_pixel_fraction_{year}"] = float(
            lit[year].sum() / n_valid
        )

    for label, y0, y1 in TRANSITIONS:
        new = common_valid & ~(arrays[y0] > 0) & (arrays[y1] > 0)
        lost = common_valid & (arrays[y0] > 0) & ~(arrays[y1] > 0)
        row[f"new_positive_fraction_{label}"] = float(new.sum() / n_valid)
        row[f"lost_positive_fraction_{label}"] = float(lost.sum() / n_valid)
        row[f"net_positive_fraction_change_{label}"] = float(
            lit[y1].sum() / n_valid - lit[y0].sum() / n_valid
        )

    if n_persistent == 0:
        return row, []

    r15 = arrays[2015][persistent].astype(np.float64, copy=False)
    r16 = arrays[2016][persistent].astype(np.float64, copy=False)
    r17 = arrays[2017][persistent].astype(np.float64, copy=False)
    r18 = arrays[2018][persistent].astype(np.float64, copy=False)

    d15 = r16 - r15
    d16 = r17 - r16
    d17 = r18 - r17
    excess = d16 - 0.5 * (d15 + d17)

    med15 = _median(d15)
    med16 = _median(d16)
    med17 = _median(d17)
    candidate_from_medians = med16 - 0.5 * (med15 + med17)

    row.update(
        {
            "median_persistent_radiance_2016": _median(r16),
            "median_d15_16_pixel": med15,
            "median_d16_17_pixel": med16,
            "median_d17_18_pixel": med17,
            "candidate_offset_from_transition_medians": (
                candidate_from_medians
            ),
            "median_pixel_excess_2017": _median(excess),
            "pixel_excess_q10": _quantile(excess, 0.10),
            "pixel_excess_q90": _quantile(excess, 0.90),
            "fraction_persistent_pixels_positive_excess": float(
                np.mean(excess > 0)
            ),
        }
    )

    bin_rows: list[dict] = []
    bin_index = np.digitize(
        r16,
        brightness_edges[1:-1],
        right=False,
    )

    for i, bin_label in enumerate(brightness_labels):
        take = bin_index == i
        n = int(take.sum())
        if n < min_bin_pixels:
            continue

        bd15 = d15[take]
        bd16 = d16[take]
        bd17 = d17[take]
        bexcess = excess[take]

        bmed15 = _median(bd15)
        bmed16 = _median(bd16)
        bmed17 = _median(bd17)

        bin_rows.append(
            {
                "NUTS_ID": nuts_id,
                "CNTR_CODE": country,
                "brightness_bin": bin_label,
                "n_pixels": n,
                "median_radiance_2016": _median(r16[take]),
                "median_d15_16_pixel": bmed15,
                "median_d16_17_pixel": bmed16,
                "median_d17_18_pixel": bmed17,
                "candidate_offset_from_transition_medians": (
                    bmed16 - 0.5 * (bmed15 + bmed17)
                ),
                "median_pixel_excess_2017": _median(bexcess),
                "fraction_pixels_positive_excess": float(
                    np.mean(bexcess > 0)
                ),
            }
        )

    return row, bin_rows


def raw_zonal_excess(config: dict) -> pd.DataFrame:
    """Calculate the same contrast from existing raw NUTS3 zonal means."""
    raw_path = (
        configured_path(config, "processed")
        / config["outputs"]["radiance_by_nuts3"]
    )
    raw = pd.read_parquet(raw_path)
    x = raw[
        raw["year"].isin(YEARS)
        & (raw["nuts_release"] == 2016)
    ][
        [
            "year",
            "NUTS_ID",
            "CNTR_CODE",
            "radiance_mean",
        ]
    ].copy()

    wide = (
        x.pivot_table(
            index=["NUTS_ID", "CNTR_CODE"],
            columns="year",
            values="radiance_mean",
            aggfunc="first",
        )
        .dropna(subset=list(YEARS))
        .reset_index()
    )

    wide["zonal_d15_16"] = wide[2016] - wide[2015]
    wide["zonal_d16_17"] = wide[2017] - wide[2016]
    wide["zonal_d17_18"] = wide[2018] - wide[2017]
    wide["zonal_candidate_offset"] = (
        wide["zonal_d16_17"]
        - 0.5 * (wide["zonal_d15_16"] + wide["zonal_d17_18"])
    )

    return wide[
        [
            "NUTS_ID",
            "CNTR_CODE",
            "zonal_d15_16",
            "zonal_d16_17",
            "zonal_d17_18",
            "zonal_candidate_offset",
        ]
    ]


def aggregate_country_regions(
    regions: pd.DataFrame,
    zonal: pd.DataFrame,
    min_persistent_pixels: int,
    min_regions: int,
) -> pd.DataFrame:
    x = regions[
        regions["read_ok"].fillna(False)
        & (regions["n_persistent_lit_pixels"].fillna(0) >= min_persistent_pixels)
    ].copy()

    x = x.merge(
        zonal,
        on=["NUTS_ID", "CNTR_CODE"],
        how="left",
        validate="one_to_one",
    )

    metric_cols = [
        "median_persistent_radiance_2016",
        "median_d15_16_pixel",
        "median_d16_17_pixel",
        "median_d17_18_pixel",
        "candidate_offset_from_transition_medians",
        "median_pixel_excess_2017",
        "fraction_persistent_pixels_positive_excess",
        "net_positive_fraction_change_2015-2016",
        "net_positive_fraction_change_2016-2017",
        "net_positive_fraction_change_2017-2018",
        "new_positive_fraction_2015-2016",
        "new_positive_fraction_2016-2017",
        "new_positive_fraction_2017-2018",
        "lost_positive_fraction_2015-2016",
        "lost_positive_fraction_2016-2017",
        "lost_positive_fraction_2017-2018",
        "zonal_d15_16",
        "zonal_d16_17",
        "zonal_d17_18",
        "zonal_candidate_offset",
    ]

    rows = []
    for country, g in x.groupby("CNTR_CODE"):
        if len(g) < min_regions:
            continue
        row: dict[str, float | int | str] = {
            "CNTR_CODE": country,
            "n_regions": int(len(g)),
            "persistent_pixels_total": int(
                g["n_persistent_lit_pixels"].sum()
            ),
        }
        for col in metric_cols:
            row[f"median_region_{col}"] = float(g[col].median())
        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        "median_region_candidate_offset_from_transition_medians"
    )


def aggregate_brightness(
    bin_regions: pd.DataFrame,
    min_regions_per_country_bin: int,
    min_countries_per_bin: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    country_rows = []

    for (country, brightness_bin), g in bin_regions.groupby(
        ["CNTR_CODE", "brightness_bin"],
        observed=True,
    ):
        if len(g) < min_regions_per_country_bin:
            continue
        country_rows.append(
            {
                "CNTR_CODE": country,
                "brightness_bin": brightness_bin,
                "n_regions": int(len(g)),
                "n_pixels": int(g["n_pixels"].sum()),
                "median_radiance_2016": float(
                    g["median_radiance_2016"].median()
                ),
                "median_d15_16_pixel": float(
                    g["median_d15_16_pixel"].median()
                ),
                "median_d16_17_pixel": float(
                    g["median_d16_17_pixel"].median()
                ),
                "median_d17_18_pixel": float(
                    g["median_d17_18_pixel"].median()
                ),
                "candidate_offset_from_transition_medians": float(
                    g["candidate_offset_from_transition_medians"].median()
                ),
                "median_pixel_excess_2017": float(
                    g["median_pixel_excess_2017"].median()
                ),
                "fraction_pixels_positive_excess": float(
                    g["fraction_pixels_positive_excess"].median()
                ),
            }
        )

    country_bins = pd.DataFrame(country_rows)

    europe_rows = []
    for brightness_bin, g in country_bins.groupby(
        "brightness_bin",
        observed=True,
    ):
        if len(g) < min_countries_per_bin:
            continue
        europe_rows.append(
            {
                "brightness_bin": brightness_bin,
                "n_countries": int(len(g)),
                "n_regions": int(g["n_regions"].sum()),
                "n_pixels": int(g["n_pixels"].sum()),
                "median_country_radiance_2016": float(
                    g["median_radiance_2016"].median()
                ),
                "median_country_d15_16": float(
                    g["median_d15_16_pixel"].median()
                ),
                "median_country_d16_17": float(
                    g["median_d16_17_pixel"].median()
                ),
                "median_country_d17_18": float(
                    g["median_d17_18_pixel"].median()
                ),
                "median_country_candidate_offset": float(
                    g["candidate_offset_from_transition_medians"].median()
                ),
                "candidate_offset_q25": float(
                    g["candidate_offset_from_transition_medians"].quantile(0.25)
                ),
                "candidate_offset_q75": float(
                    g["candidate_offset_from_transition_medians"].quantile(0.75)
                ),
                "median_country_pixel_excess": float(
                    g["median_pixel_excess_2017"].median()
                ),
            }
        )

    return country_bins, pd.DataFrame(europe_rows)


def country_significance(country: pd.DataFrame) -> dict[str, float]:
    values = country[
        "median_region_candidate_offset_from_transition_medians"
    ].dropna().to_numpy(dtype=float)

    n = len(values)
    n_positive = int(np.sum(values > 0))
    sign_p = (
        float(binomtest(n_positive, n, p=0.5, alternative="greater").pvalue)
        if n
        else np.nan
    )

    try:
        wilcox_p = float(
            wilcoxon(values, alternative="greater", zero_method="wilcox").pvalue
        )
    except ValueError:
        wilcox_p = np.nan

    return {
        "n_countries": n,
        "n_positive_country_offsets": n_positive,
        "fraction_positive_country_offsets": (
            n_positive / n if n else np.nan
        ),
        "country_sign_test_p": sign_p,
        "country_wilcoxon_greater_p": wilcox_p,
        "median_country_offset": (
            float(np.median(values)) if n else np.nan
        ),
        "country_offset_q25": (
            float(np.quantile(values, 0.25)) if n else np.nan
        ),
        "country_offset_q75": (
            float(np.quantile(values, 0.75)) if n else np.nan
        ),
    }


def region_zonal_association(
    regions: pd.DataFrame,
    zonal: pd.DataFrame,
    min_persistent_pixels: int,
) -> dict[str, float]:
    x = regions[
        regions["read_ok"].fillna(False)
        & (regions["n_persistent_lit_pixels"].fillna(0) >= min_persistent_pixels)
    ][
        [
            "NUTS_ID",
            "CNTR_CODE",
            "candidate_offset_from_transition_medians",
        ]
    ].merge(
        zonal[
            [
                "NUTS_ID",
                "CNTR_CODE",
                "zonal_candidate_offset",
            ]
        ],
        on=["NUTS_ID", "CNTR_CODE"],
        how="inner",
    ).dropna()

    if len(x) < 3:
        return {
            "n_regions": len(x),
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
        }

    rho, p = spearmanr(
        x["candidate_offset_from_transition_medians"],
        x["zonal_candidate_offset"],
    )
    return {
        "n_regions": int(len(x)),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
    }


def plot_country_transitions(country: pd.DataFrame, output: Path, dpi: int) -> None:
    c = country.sort_values(
        "median_region_candidate_offset_from_transition_medians"
    ).copy()

    x = np.arange(len(c))
    fig, ax = plt.subplots(figsize=(max(11, len(c) * 0.34), 6.4))

    series = [
        (
            "median_region_median_d15_16_pixel",
            "2015→2016",
            "o",
            "-",
        ),
        (
            "median_region_median_d16_17_pixel",
            "2016→2017",
            "s",
            "-",
        ),
        (
            "median_region_median_d17_18_pixel",
            "2017→2018",
            "^",
            "-",
        ),
    ]

    for col, label, marker, linestyle in series:
        ax.plot(
            x,
            c[col],
            marker=marker,
            linestyle=linestyle,
            linewidth=1.2,
            markersize=4.5,
            label=label,
        )

    ax.axhline(0, color="0.4", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(c["CNTR_CODE"], rotation=90)
    ax.set_ylabel(
        "Country median of NUTS3 persistent-lit pixel change\n"
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title(
        "Annual VNL V2.2: country consistency of the 2016→2017 radiance jump"
    )
    ax.grid(axis="y", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_country_excess(country: pd.DataFrame, output: Path, dpi: int) -> None:
    c = country.sort_values(
        "median_region_candidate_offset_from_transition_medians"
    ).copy()

    fig, ax = plt.subplots(figsize=(9.2, max(6, len(c) * 0.25)))
    y = np.arange(len(c))

    ax.scatter(
        c["median_region_candidate_offset_from_transition_medians"],
        y,
        s=28,
    )
    ax.axvline(0, color="0.4", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(c["CNTR_CODE"])
    ax.set_xlabel(
        "2017 excess contrast in persistent-lit pixels\n"
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_ylabel("Country")
    ax.set_title(
        "Country-level 2017 excess after removing flanking-year change"
    )
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_brightness(
    europe_bins: pd.DataFrame,
    output: Path,
    dpi: int,
) -> None:
    if europe_bins.empty:
        return

    b = europe_bins.sort_values("median_country_radiance_2016").copy()

    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    ax.plot(
        b["median_country_radiance_2016"],
        b["median_country_candidate_offset"],
        marker="o",
        linewidth=1.8,
        label="Median country excess",
    )
    ax.fill_between(
        b["median_country_radiance_2016"],
        b["candidate_offset_q25"],
        b["candidate_offset_q75"],
        alpha=0.18,
        linewidth=0,
        label="Country IQR",
    )
    ax.axhline(0, color="0.4", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel(
        "Median 2016 radiance of persistent-lit pixels "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_ylabel(
        "2017 excess contrast (nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title(
        "Does the retained-pixel 2017 step depend on pixel brightness?"
    )
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_zonal_vs_pixel(
    regions: pd.DataFrame,
    zonal: pd.DataFrame,
    min_persistent_pixels: int,
    output: Path,
    dpi: int,
) -> None:
    x = regions[
        regions["read_ok"].fillna(False)
        & (regions["n_persistent_lit_pixels"].fillna(0) >= min_persistent_pixels)
    ][
        [
            "NUTS_ID",
            "CNTR_CODE",
            "candidate_offset_from_transition_medians",
        ]
    ].merge(
        zonal[
            [
                "NUTS_ID",
                "CNTR_CODE",
                "zonal_candidate_offset",
            ]
        ],
        on=["NUTS_ID", "CNTR_CODE"],
        how="inner",
    ).dropna()

    if x.empty:
        return

    fig, ax = plt.subplots(figsize=(7.2, 7.0))
    ax.scatter(
        x["candidate_offset_from_transition_medians"],
        x["zonal_candidate_offset"],
        s=9,
        alpha=0.35,
    )
    limits = np.nanquantile(
        np.r_[
            x["candidate_offset_from_transition_medians"].to_numpy(),
            x["zonal_candidate_offset"].to_numpy(),
        ],
        [0.01, 0.99],
    )
    lo, hi = float(limits[0]), float(limits[1])
    ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(
        "Persistent-lit pixel excess contrast\n"
        "(median within NUTS3)"
    )
    ax.set_ylabel(
        "Existing raw NUTS3-mean excess contrast"
    )
    ax.set_title(
        "Does the pixel-level 2017 step track the observed NUTS3-mean jump?"
    )
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    path: Path,
    country: pd.DataFrame,
    europe_bins: pd.DataFrame,
    region_stats: pd.DataFrame,
    sign_stats: dict[str, float],
    assoc: dict[str, float],
) -> None:
    with path.open("w", encoding="utf-8") as h:
        h.write(
            "ANNUAL VNL V2.2 2017 PIXEL-DISCONTINUITY DIAGNOSTIC\n"
        )
        h.write("=" * 76 + "\n\n")
        h.write(
            "Diagnostic only: no production radiance values were modified.\n\n"
        )

        h.write("PRIMARY QUESTION\n")
        h.write(
            "Do pixels that are positive in ALL of 2015, 2016, 2017 and 2018 "
            "still show an anomalous 2016->2017 radiance jump?\n\n"
        )

        h.write("COUNTRY CONSISTENCY\n")
        for key, value in sign_stats.items():
            h.write(f"{key}: {value}\n")

        h.write("\nCOUNTRY TABLE\n")
        cols = [
            "CNTR_CODE",
            "n_regions",
            "median_region_median_d15_16_pixel",
            "median_region_median_d16_17_pixel",
            "median_region_median_d17_18_pixel",
            "median_region_candidate_offset_from_transition_medians",
            "median_region_median_pixel_excess_2017",
            "median_region_zonal_candidate_offset",
            "median_region_net_positive_fraction_change_2015-2016",
            "median_region_net_positive_fraction_change_2016-2017",
            "median_region_net_positive_fraction_change_2017-2018",
        ]
        h.write(
            country[cols].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\n")

        h.write("BRIGHTNESS DEPENDENCE\n")
        if europe_bins.empty:
            h.write("No brightness bins passed minimum support.\n")
        else:
            h.write(
                europe_bins.to_string(
                    index=False,
                    float_format=lambda x: f"{x:.6g}",
                )
            )
            h.write("\n")

        h.write("\nPIXEL-LEVEL VS EXISTING NUTS3-MEAN JUMP\n")
        for key, value in assoc.items():
            h.write(f"{key}: {value}\n")

        h.write("\nMASK/SUPPORT CHECK\n")
        valid_regions = region_stats[
            region_stats["read_ok"].fillna(False)
            & (region_stats["n_common_valid_pixels"].fillna(0) > 0)
        ]
        for label, _, _ in TRANSITIONS:
            net = valid_regions[
                f"net_positive_fraction_change_{label}"
            ].dropna()
            h.write(
                f"{label}: median NUTS3 net positive-pixel fraction change "
                f"= {net.median():.6g}; "
                f"Q25={net.quantile(0.25):.6g}; "
                f"Q75={net.quantile(0.75):.6g}\n"
            )

        h.write("\nINTERPRETATION GUIDE\n")
        h.write(
            "- A large 2016->2017 change in positive-pixel support would "
            "suggest that mask/support behaviour contributes to the jump.\n"
        )
        h.write(
            "- A positive 2017 excess among pixels positive in all four years "
            "demonstrates a radiance step in retained lit pixels; background "
            "zeroing alone cannot explain it.\n"
        )
        h.write(
            "- Similar positive excess across countries would agree with the "
            "previous observation that the 2016->2017 jump was geographically "
            "consistent rather than country-specific.\n"
        )
        h.write(
            "- Brightness dependence determines whether a simple additive "
            "offset is plausible for retained pixels. A strong dependence on "
            "2016 radiance would argue against one universal additive value.\n"
        )
        h.write(
            "- The existing raw NUTS3-mean excess is included only to connect "
            "this direct pixel diagnostic to the earlier country/NUTS3 result.\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["annual_vnl_2017_pixel_test"]

    configured_years = tuple(int(x) for x in cfg["years"])
    if configured_years != YEARS:
        raise ValueError(
            f"This diagnostic expects years {YEARS}; found {configured_years}"
        )

    rasters = discover_rasters(config, list(YEARS))
    nuts_release = int(cfg["nuts_release"])
    nuts_path = nuts_path_for_release(config, nuts_release)

    if not nuts_path.exists():
        raise FileNotFoundError(
            f"NUTS {nuts_release} file not found: {nuts_path}"
        )

    nuts = gpd.read_file(nuts_path)
    if nuts.crs is None:
        raise ValueError(f"NUTS file has no CRS: {nuts_path}")

    id_field = config["nuts"]["id_field"]
    required = [id_field, "CNTR_CODE", "geometry"]
    missing = [c for c in required if c not in nuts.columns]
    if missing:
        raise KeyError(f"NUTS file missing columns: {missing}")

    min_region_persistent = int(cfg["min_region_persistent_pixels"])
    min_bin_pixels = int(cfg["min_region_bin_pixels"])
    min_country_regions = int(cfg["min_country_regions"])
    min_country_bin_regions = int(cfg["min_country_bin_regions"])
    min_europe_bin_countries = int(cfg["min_europe_bin_countries"])
    dpi = int(cfg["figure_dpi"])

    finite_edges = [float(x) for x in cfg["brightness_bin_edges"]]
    brightness_edges = np.array(
        [-np.inf] + finite_edges + [np.inf],
        dtype=float,
    )
    brightness_labels = list(cfg["brightness_bin_labels"])
    if len(brightness_labels) != len(brightness_edges) - 1:
        raise ValueError(
            "brightness_bin_labels must have one fewer entry than the "
            "full bin-edge vector including -inf/+inf"
        )

    print("ANNUAL VNL V2.2 2017 PIXEL-DISCONTINUITY DIAGNOSTIC")
    print("Diagnostic only; production radiance will not be changed.")
    print("Rasters:")
    for year in YEARS:
        print(f"  {year}: {rasters[year]}")

    region_rows: list[dict] = []
    bin_rows: list[dict] = []

    with ExitStack() as stack:
        sources = {
            year: stack.enter_context(rasterio.open(rasters[year]))
            for year in YEARS
        }
        validate_alignment(sources)

        raster_crs = sources[YEARS[0]].crs
        if nuts.crs != raster_crs:
            nuts = nuts.to_crs(raster_crs)

        total = len(nuts)
        for idx, row in enumerate(nuts.itertuples(index=False), start=1):
            nuts_id = getattr(row, id_field)
            country = getattr(row, "CNTR_CODE")
            geom = getattr(row, "geometry")

            if geom is None or geom.is_empty:
                continue

            region, bins = region_pixel_statistics(
                str(nuts_id),
                str(country),
                geom,
                sources,
                brightness_edges,
                brightness_labels,
                min_bin_pixels,
            )
            region_rows.append(region)
            bin_rows.extend(bins)

            if idx % 50 == 0 or idx == total:
                print(
                    f"  processed {idx:,}/{total:,} NUTS3 regions",
                    flush=True,
                )

    region_stats = pd.DataFrame(region_rows)
    bin_region_stats = pd.DataFrame(bin_rows)
    zonal = raw_zonal_excess(config)

    country = aggregate_country_regions(
        region_stats,
        zonal,
        min_persistent_pixels=min_region_persistent,
        min_regions=min_country_regions,
    )

    country_bins, europe_bins = aggregate_brightness(
        bin_region_stats,
        min_regions_per_country_bin=min_country_bin_regions,
        min_countries_per_bin=min_europe_bin_countries,
    )

    sign_stats = country_significance(country)
    assoc = region_zonal_association(
        region_stats,
        zonal,
        min_persistent_pixels=min_region_persistent,
    )

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "annual_change"
        / "annual_vnl_2017_pixel_test"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    region_path = processed / "annual_vnl_2017_pixel_nuts3.csv"
    country_path = processed / "annual_vnl_2017_pixel_country.csv"
    region_bin_path = processed / "annual_vnl_2017_pixel_nuts3_brightness.csv"
    country_bin_path = processed / "annual_vnl_2017_pixel_country_brightness.csv"
    europe_bin_path = processed / "annual_vnl_2017_pixel_brightness_summary.csv"
    summary_path = processed / "annual_vnl_2017_pixel_summary.txt"

    region_stats.to_csv(region_path, index=False)
    country.to_csv(country_path, index=False)
    bin_region_stats.to_csv(region_bin_path, index=False)
    country_bins.to_csv(country_bin_path, index=False)
    europe_bins.to_csv(europe_bin_path, index=False)

    country_transitions_fig = (
        figures / "country_transition_medians.png"
    )
    country_excess_fig = figures / "country_2017_excess.png"
    brightness_fig = figures / "brightness_dependence.png"
    zonal_pixel_fig = figures / "nuts3_zonal_vs_persistent_pixel_excess.png"

    plot_country_transitions(country, country_transitions_fig, dpi)
    plot_country_excess(country, country_excess_fig, dpi)
    plot_brightness(europe_bins, brightness_fig, dpi)
    plot_zonal_vs_pixel(
        region_stats,
        zonal,
        min_region_persistent,
        zonal_pixel_fig,
        dpi,
    )

    write_summary(
        summary_path,
        country,
        europe_bins,
        region_stats,
        sign_stats,
        assoc,
    )

    print("\nCOUNTRY CONSISTENCY")
    for key, value in sign_stats.items():
        print(f"  {key}: {value}")

    print("\nCOUNTRY MEDIAN TRANSITIONS")
    print(
        country[
            [
                "CNTR_CODE",
                "n_regions",
                "median_region_median_d15_16_pixel",
                "median_region_median_d16_17_pixel",
                "median_region_median_d17_18_pixel",
                "median_region_candidate_offset_from_transition_medians",
                "median_region_zonal_candidate_offset",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )

    print("\nBRIGHTNESS DEPENDENCE")
    if europe_bins.empty:
        print("  No bins passed minimum support.")
    else:
        print(
            europe_bins.to_string(
                index=False,
                float_format=lambda x: f"{x:.6f}",
            )
        )

    print("\nPIXEL VS EXISTING NUTS3-MEAN EXCESS")
    for key, value in assoc.items():
        print(f"  {key}: {value}")

    print("\nWROTE")
    for path in [
        region_path,
        country_path,
        region_bin_path,
        country_bin_path,
        europe_bin_path,
        summary_path,
        country_transitions_fig,
        country_excess_fig,
        brightness_fig,
        zonal_pixel_fig,
    ]:
        print(f"  {path}")


if __name__ == "__main__":
    main()
