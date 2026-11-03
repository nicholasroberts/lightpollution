#!/usr/bin/env python3
"""Figure 1E: seven country trajectories of corrected NUTS3 radiance change.

2014–2024: same locked radiance history, fixed 2014 valid-area weights
and 2021 approved native-pixel snow patch as Figure 1C. Countries match
Figure 2D: TR, DE, FR, IE, IT, NL, UK. NOT the smaller Figure 2D population
matched cohort: these country trajectories use all NUTS2016 NUTS3 regions
with complete calibrated radiance for 2014–2024 (same as Fig 1C).

Plot ΔR(t) = area-weighted country mean R(t) - area-weighted R(2014),
so all seven series start at zero. Absolute radiance rather than percent
is chosen to avoid unstable ratios in low-radiance regions; the annual
CSV includes both absolute and percentage changes for optional review.

Scientific input files, calibration, native-pixel 2021 patch and original
Figure 1A–D artwork are read-only. SVG and PDF editable; 16:9 PNG for talk.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figure_style import apply_figure_style, save_figure
from historical_prc_pixel_nuts3 import load_nuts
from project_config import PROJECT_ROOT as ROOT, load_config
from plot_publication_figure1_radiance import (
    DEFAULT_ENDPOINT, DEFAULT_ANNUAL, DEFAULT_SNOW_OUTPUT,
    get_locked_snow_2021, validate_endpoint
)

YEARS = tuple(range(2014, 2025))
COUNTRIES = ("TR", "DE", "FR", "IE", "IT", "NL", "UK")
NAMES = {
    "TR": "Türkiye", "DE": "Germany", "FR": "France",
    "IE": "Ireland", "IT": "Italy", "NL": "Netherlands",
    "UK": "United Kingdom",
}
COLOURS = {
    "TR": "#0072B2",  # Figure 2 country colours retained
    "DE": "#E69F00",
    "FR": "#009E73",
    "IE": "#CC79A7",
    "IT": "#D55E00",
    "NL": "#7B61A8",
    "UK": "#333333",
}
RADIANCE = "radiance_mean_corrected"
WEIGHT = "valid_area_km2_2014"
INPUT_FIG1C_QA = (
    ROOT / "figures/publication/figure1_radiance/"
    "fig1C_europe_annual_radiance_support.csv"
)
OUTPUT = ROOT / "figures/publication/figure1_radiance"
UNITS = r"nW cm$^{-2}$ sr$^{-1}$"


def resolve(path: str | Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def verified_regional_annual(annual_file: Path, endpoints: pd.DataFrame,
                             snow: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Replicate Fig1C regional cohort and 2021 QA, with region rows retained."""
    raw = pd.read_parquet(
        annual_file, columns=["NUTS_ID", "year", RADIANCE]
    )
    raw["NUTS_ID"] = raw.NUTS_ID.astype(str)
    raw["year"] = pd.to_numeric(raw.year, errors="raise").astype(int)
    raw = raw.loc[raw.year.isin(YEARS)].copy()
    if raw.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate annual NUTS2016 regional radiance records")
    if set(raw.year) != set(YEARS):
        raise ValueError(f"Missing years: {sorted(set(YEARS)-set(raw.year))}")
    cohort = set(endpoints.NUTS_ID)
    for year in YEARS:
        cohort &= set(raw.loc[raw.year.eq(year), "NUTS_ID"])
    if len(cohort) < 1400:
        raise ValueError(f"Only {len(cohort)} regions have complete annual history")
    # Exact same source endpoint area weights and all-years cohort as Fig1C.
    reference = endpoints.loc[
        endpoints.NUTS_ID.isin(cohort),
        ["NUTS_ID", WEIGHT, "radiance_2014", "radiance_2024"]
    ]
    regional = raw.loc[raw.NUTS_ID.isin(cohort)].merge(
        reference, on="NUTS_ID", validate="many_to_one"
    )
    if len(regional) != len(cohort) * len(YEARS):
        raise ValueError("Incomplete balanced NUTS2016 annual region panel")
    regional[RADIANCE] = pd.to_numeric(regional[RADIANCE], errors="coerce")
    if not np.isfinite(regional[[RADIANCE, WEIGHT]].to_numpy(float)).all():
        raise ValueError("Missing calibrated radiance or valid 2014 area")
    if (regional[RADIANCE] < 0).any() or (regional[WEIGHT] <= 0).any():
        raise ValueError("Negative radiance or non-positive fixed 2014 area")
    endpoint_checks = {}
    for year, column in ((2014, "radiance_2014"), (2024, "radiance_2024")):
        subset = regional.loc[regional.year.eq(year)]
        errors = (subset[RADIANCE] - subset[column]).abs().to_numpy(float)
        p99 = float(np.quantile(errors, .99))
        if p99 > .025:
            raise ValueError(
                f"{year} mean radiance does not match locked endpoints "
                f"(p99 absolute difference {p99:.5g} > .025 nW)"
            )
        endpoint_checks[str(year)] = {
            "n_regions": len(subset),
            "max_abs_error_nw": float(errors.max()),
            "p99_abs_error_nw": p99
        }
    if len(snow) != 1522 or snow.NUTS_ID.duplicated().any():
        raise ValueError("Locked 2021 NUTS3 source does not contain 1,522 regions")
    winter = regional.loc[regional.year.eq(2021), ["NUTS_ID", RADIANCE]]
    joined = winter.merge(
        snow[["NUTS_ID", "radiance_2021_original_calibrated",
              "radiance_2021_snow_corrected"]],
        on="NUTS_ID", validate="one_to_one"
    )
    if len(joined) != len(winter):
        raise ValueError("Locked 2021 input missing annual-cohort region(s)")
    current = joined[RADIANCE].to_numpy(float)
    original = joined.radiance_2021_original_calibrated.to_numpy(float)
    corrected = joined.radiance_2021_snow_corrected.to_numpy(float)
    difference = np.abs(current - original)
    tolerance = 2e-5 * np.maximum(1., np.abs(current))
    if (difference > tolerance).any():
        raise ValueError(
            f"Annual 2021 original vs locked source disagree (max={difference.max():.7g})"
        )
    if not np.isfinite(corrected).all() or (corrected < 0).any():
        raise ValueError("Non-finite or negative locked 2021 radiance")
    if (corrected > original + 2e-6).any():
        raise ValueError("Snow correction must not brighten NUTS3 radiance")
    patched = pd.Series(corrected, index=joined.NUTS_ID)
    use = regional.year.eq(2021)
    regional.loc[use, RADIANCE] = (
        regional.loc[use, "NUTS_ID"].map(patched).to_numpy(float)
    )
    regional["locked_snow_2021_used"] = regional.year.eq(2021)
    qa = {
        "all_years_balanced_cohort_regions": len(cohort),
        "2014_2024_locked_endpoint_comparison": endpoint_checks,
        "2021_locked_patch_nuts3_matched": len(joined),
        "2021_unpatched_vs_locked_original_max_abs_error_nw":
            float(difference.max()),
        "2021_snow_reduced_regions": int((corrected < original - 1e-6).sum()),
        "valid_area_weights": WEIGHT,
        "region_areas_constant_throughout": True,
        "pixel_offset_applied_once_upstream": True,
    }
    return regional, qa


def country_and_europe_weighted(regional: pd.DataFrame,
                                nuts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Compute weighted regional means independently within country and Europe."""
    name_map = nuts[["NUTS_ID", "CNTR_CODE"]]
    if nuts.NUTS_ID.duplicated().any() or len(nuts) != 1522:
        raise ValueError("Changed fixed NUTS2016 geometry")
    area = regional.merge(name_map, on="NUTS_ID",
                          validate="many_to_one")
    if len(area) != len(regional) or area.CNTR_CODE.isna().any():
        raise ValueError("Country join dropped one or more NUTS3 regions")
    area["weighted_total"] = area[RADIANCE] * area[WEIGHT]
    result = area.loc[area.CNTR_CODE.isin(COUNTRIES)].groupby(
        ["CNTR_CODE", "year"], as_index=False
    ).agg(
        weighted_total=("weighted_total", "sum"),
        fixed_area_km2=(WEIGHT, "sum"),
        n_nuts3=("NUTS_ID", "nunique"),
    )
    result["mean_radiance_nw"] = (
        result.weighted_total / result.fixed_area_km2
    )
    expected = pd.MultiIndex.from_product(
        [COUNTRIES, YEARS], names=["CNTR_CODE", "year"]
    )
    got = pd.MultiIndex.from_frame(result[["CNTR_CODE","year"]])
    if not got.is_unique or set(got) != set(expected):
        raise ValueError("Country/year series not complete for all seven countries")
    n = result.groupby("CNTR_CODE").n_nuts3.nunique()
    if not n.eq(1).all():
        raise ValueError("A country changed NUTS3 support between years")
    sizes = result.groupby("CNTR_CODE").fixed_area_km2.agg(["min","max"])
    if not np.allclose(sizes["min"], sizes["max"], rtol=0, atol=1e-5):
        raise ValueError("Fixed area weights varied across year in a country")
    baseline = result.loc[result.year.eq(2014)].set_index("CNTR_CODE")
    result["country"] = result.CNTR_CODE.map(NAMES)
    result["baseline_2014_mean_nw"] = result.CNTR_CODE.map(
        baseline.mean_radiance_nw
    )
    result["change_from_2014_nw"] = (
        result.mean_radiance_nw - result.baseline_2014_mean_nw
    )
    result["pct_change_from_2014"] = (
        100 * result.change_from_2014_nw / result.baseline_2014_mean_nw
    )
    if (result.baseline_2014_mean_nw <= 0).any():
        raise ValueError("One country has zero 2014 radiance baseline")
    if not np.allclose(
        result.loc[result.year.eq(2014), "change_from_2014_nw"], 0,
        atol=1e-12
    ):
        raise ValueError("A country's 2014 delta does not equal zero")
    # Independent entire-cohort aggregate, expected to reproduce Figure 1C.
    europe = area.groupby("year", as_index=False).agg(
        weighted_total=("weighted_total", "sum"),
        fixed_area_km2=(WEIGHT, "sum"),
        n_nuts3=("NUTS_ID", "nunique"),
    )
    europe["area_weighted_radiance_nw"] = (
        europe.weighted_total / europe.fixed_area_km2
    )
    baseline_europe = europe.loc[europe.year.eq(2014),
                                 "area_weighted_radiance_nw"].iloc[0]
    europe["change_from_2014_nw"] = (
        europe.area_weighted_radiance_nw - baseline_europe
    )
    country_counts = (
        result.loc[result.year.eq(2014), ["CNTR_CODE","n_nuts3","fixed_area_km2"]]
        .set_index("CNTR_CODE").to_dict(orient="index")
    )
    return result, europe, country_counts


def compare_figure1c(europe: pd.DataFrame, reference: Path) -> dict:
    if not reference.is_file():
        return {
            "original_fig1c_csv_present": False,
            "comparison_passed": None,
            "note": "Figure 1C CSV not found; independent EU values exported for comparison."
        }
    old = pd.read_csv(reference)
    required = {"year","area_weighted_radiance_nw"}
    if not required.issubset(old):
        raise ValueError("Existing Fig1C CSV lacks its annual weighted radiance")
    paired = europe.merge(
        old[["year","area_weighted_radiance_nw"]].rename(
            columns={"area_weighted_radiance_nw":"fig1c_mean_nw"}
        ), on="year", validate="one_to_one"
    )
    if len(paired) != len(YEARS):
        raise ValueError("Existing Fig1C not a complete 2014-2024 reference")
    error = np.abs(
        paired.area_weighted_radiance_nw.to_numpy(float)
        - paired.fig1c_mean_nw.to_numpy(float)
    )
    if not np.allclose(paired.area_weighted_radiance_nw,
                       paired.fig1c_mean_nw,
                       rtol=2e-10, atol=2e-9):
        raise ValueError(
            f"Reproduced Europe mean deviates from Figure 1C: "
            f"max abs difference={error.max():.8g} nW"
        )
    return {
        "original_fig1c_csv_present": True,
        "comparison_passed": True,
        "source": str(reference),
        "max_abs_difference_nw": float(error.max()),
        "n_years_compared": len(paired),
    }


def plot_figure(country: pd.DataFrame, path: Path, png_dpi: int):
    fig, ax = plt.subplots(figsize=(12.8, 7.2), facecolor="white")
    fig.subplots_adjust(left=.115, right=.965, top=.875, bottom=.255)
    for code in COUNTRIES:
        sub = country.loc[country.CNTR_CODE.eq(code)].sort_values("year")
        ax.plot(sub.year, sub.change_from_2014_nw,
                color=COLOURS[code], linewidth=2.25,
                marker="o", markersize=4, label=NAMES[code])
    ax.axhline(0, color="#777777", linewidth=.9, zorder=0)
    ax.axvline(2021, color="#999999", linewidth=.7,
               linestyle=":", zorder=0)
    ax.set_xlim(2013.8, 2024.2)
    ax.set_xticks(YEARS)
    ax.set_xlabel("Year", labelpad=8)
    ax.set_ylabel(f"Change in country mean radiance since 2014 ({UNITS})")
    ax.set_title(
        "1E  |  Country trajectories in calibrated night-time radiance, 2014–2024",
        loc="left", fontsize=16, pad=16
    )
    ax.grid(axis="y", linewidth=.6, color="#e7e7e7")
    ax.legend(
        loc="upper center", bbox_to_anchor=(.5, -.145),
        ncol=4, frameon=False, fontsize=11,
        handlelength=2.4, columnspacing=1.35
    )
    fig.text(.115, .035,
             "Fixed NUTS2016 regions; 2014 area-weighted country means = zero baseline. "
             "2021 uses locked native-pixel snow correction.",
             fontsize=8.8, color="#555555")
    saved = save_figure(fig, path)
    if png_dpi > 0:
        # Standard 16:9 slide canvas; do not trim.
        png = path.with_suffix(".png")
        fig.savefig(png, format="png", dpi=png_dpi, facecolor="white")
        saved["png"] = png
    plt.close(fig)
    return saved


def self_test():
    # A real land-area mean, not mean of NUTS units:
    # 1 km² at 10 nW and 9 km² at 0 nW => 1 nW, not 5.
    means = np.array([10., 0.])
    areas = np.array([1., 9.])
    assert np.isclose(np.average(means, weights=areas), 1.)
    assert not np.isclose(np.average(means), 1.)
    assert len(YEARS) == 11 and len(COUNTRIES) == 7
    assert set(COUNTRIES) == {"TR","DE","FR","IE","IT","NL","UK"}
    print("Figure 1E seven-country area weighting and zero-baseline tests: PASS",
          flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--endpoints", default=DEFAULT_ENDPOINT)
    p.add_argument("--annual", default=DEFAULT_ANNUAL)
    p.add_argument("--snow-cache", default=DEFAULT_SNOW_OUTPUT)
    p.add_argument("--fig1c-support", default=INPUT_FIG1C_QA)
    p.add_argument("--out", default=OUTPUT)
    p.add_argument("--png-dpi", type=int, default=240)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    self_test()
    if args.self_test:
        return
    if not 0 <= args.png_dpi <= 600:
        p.error("PNG dpi must be 0..600")
    apply_figure_style()
    endpoints = resolve(args.endpoints)
    annual_file = resolve(args.annual)
    out = resolve(args.out)
    data = validate_endpoint(
        pd.read_csv(endpoints, dtype={"NUTS_ID": str})
    )
    # Reuse Figure 1C's validated native-2021 locked patch provenance.
    snow, snow_meta = get_locked_snow_2021(resolve(args.snow_cache))
    regional, qa = verified_regional_annual(annual_file, data, snow)
    nuts = load_nuts(load_config())
    country, europe, counts = country_and_europe_weighted(regional, nuts)
    europe_check = compare_figure1c(europe, resolve(args.fig1c_support))

    # No change to any Figure 1A–D output, source dataset or prior cache.
    out.mkdir(parents=True, exist_ok=True)
    country.to_csv(out/"fig1E_seven_country_annual_radiance_changes.csv",
                   index=False)
    europe.to_csv(out/"fig1E_europe_reference_area_weighted_qa.csv",
                  index=False)
    stem=out/"fig1E_seven_country_radiance_change_2014_2024"
    plots=plot_figure(country, stem, args.png_dpi)
    europe2014=float(
        europe.loc[europe.year.eq(2014), "area_weighted_radiance_nw"].iloc[0]
    )
    status={
        "publication_subfigure":"1E",
        "definition": "R_country(year) - R_country(2014)",
        "country_codes": list(COUNTRIES),
        "country_names": NAMES,
        "country_nuts3_count_area": counts,
        "years": list(YEARS),
        "reference_dataset_2014_2024":"same full radiance cohort as Figure 1C",
        "NOT_population_matched_Fig2D_cohort": True,
        "radiance_units":"nW cm^-2 sr^-1",
        "area_weight":"valid_area_km2_2014, fixed across all years",
        "2014_start_absolute_europe_mean_nw": europe2014,
        "2014_baseline_change_for_all_countries":0.,
        "source_endpoint_csv":str(endpoints),
        "source_annual_parquet":str(annual_file),
        "2021_locked_snow_patch_audit":snow_meta,
        "regional_scientific_qa":qa,
        "cross_check_against_existing_Fig1C":europe_check,
        "vector_format":["svg","pdf"],
        "talk_png_dpi":args.png_dpi,
        "source_tables_or_other_figures_modified":False,
    }
    (out/"fig1E_seven_country_trajectories_provenance.json").write_text(
        json.dumps(status, indent=2)+"\n", encoding="utf-8"
    )
    print("\nFIGURE 1E RESULTS", flush=True)
    print(f"2014 Europe area-weighted mean, same source as Fig1C: "
          f"{europe2014:.6f} nW cm^-2 sr^-1", flush=True)
    print("Independent Fig1C cross-check:", europe_check, flush=True)
    show=country.loc[country.year.isin((2014, 2020, 2024)),
        ["CNTR_CODE","year","n_nuts3","mean_radiance_nw",
         "change_from_2014_nw"]]
    print(show.sort_values(["CNTR_CODE","year"]).to_string(index=False),
          flush=True)
    print("\nOUTPUTS", flush=True)
    for filename in [
        *plots.values(),
        out/"fig1E_seven_country_annual_radiance_changes.csv",
        out/"fig1E_europe_reference_area_weighted_qa.csv",
        out/"fig1E_seven_country_trajectories_provenance.json",
    ]:
        print(f"  {filename} ({filename.stat().st_size/1048576:.2f} MiB)",
              flush=True)


if __name__ == "__main__":
    main()
