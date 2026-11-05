#!/usr/bin/env python3
"""Country 2014–2024 area shares with NUTS3 mean spider PRC <= DoP=0.23 scenario.

Standalone extension of the existing PRC threshold atlas; writes additional
COUNTRY time-series figures/tables to exactly the same publication folder.
Reads the same locked 1,522 × 11 spider PRC panel. Country IDs are taken
from the original fixed NUTS2016 geometry, NOT from population cohorts.

DO NOT interpret as the fraction of individual VIIRS pixels (or animals) below
a behavioural threshold. It is the fixed-valid-area fraction of whole NUTS3
regions whose annual native-pixel-derived MEAN PRC falls below the scenario.
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
from project_config import load_config
from historical_prc_pixel_nuts3 import load_nuts
from plot_publication_spider_prc_dop023_threshold import (
    ANNUAL, JOINED_ENDPOINT_QA, OUTPUT, DOP_THRESHOLD,
    POLARIZATION_SENSITIVITY, YEARS, EXPECTED_REGIONS, dop_to_prc,
    annotate_regions, read_spider_history, independent_endpoint_qa, resolve,
)

# Country exemplars fixed in the existing Fig. 2/Fig. 4 analyses,
# with UK INCLUDED because the PRC cohort does not require population data.
COUNTRIES = ("NL", "UK", "TR", "DE", "FR", "IE", "IT")
COUNTRY_NAMES = {
    "NL": "Netherlands",
    "UK": "United Kingdom",
    "TR": "Türkiye",
    "DE": "Germany",
    "FR": "France",
    "IE": "Ireland",
    "IT": "Italy",
    "EL": "Greece",
    "ES": "Spain",
    "PT": "Portugal",
    "RO": "Romania",
}
COLOURS = {
    "NL": "#875bb4", "UK": "#0072B2", "TR": "#D55E00",
    "DE": "#E69F00", "FR": "#009E73", "IE": "#CC79A7",
    "IT": "#7b5c49",
}
METRIC = "pct_fixed_2014_area_at_or_below"


def country_geometry_lookup(geometry: pd.DataFrame, expected: int) -> pd.DataFrame:
    """Every locked NUTS2016 ID mapped to its original country code once."""
    needed = {"NUTS_ID", "CNTR_CODE"}
    if needed - set(geometry):
        raise ValueError("NUTS2016 geometry lacks NUTS_ID / CNTR_CODE")
    d = geometry[["NUTS_ID", "CNTR_CODE"]].copy()
    if len(d) != expected or not d.NUTS_ID.is_unique or d.NUTS_ID.isna().any():
        raise ValueError("Country lookup must contain each frozen NUTS2016 region exactly once")
    if d.CNTR_CODE.isna().any() or d.CNTR_CODE.astype(str).str.len().ne(2).any():
        raise ValueError("Missing or invalid NUTS2016 country identifiers")
    if not d.NUTS_ID.astype(str).str[:2].eq(d.CNTR_CODE.astype(str)).all():
        raise ValueError("NUTS2016 code and assigned country mismatch")
    return d


def country_area_series(classified: pd.DataFrame, geometry: pd.DataFrame,
                        expected_regions: int = EXPECTED_REGIONS) -> pd.DataFrame:
    """Country-level copy of the full-Europe fixed 2014 valid-area estimator."""
    lookup = country_geometry_lookup(geometry, expected_regions)
    x = classified.merge(
        lookup, on="NUTS_ID", how="left", validate="many_to_one",
        indicator=True
    )
    if not x._merge.eq("both").all() or len(x) != expected_regions * len(YEARS):
        raise ValueError("Country assignment lost annual PRC regions")
    x.drop(columns="_merge", inplace=True)
    if x.groupby("year").NUTS_ID.nunique().ne(expected_regions).any():
        raise ValueError("Unexpected annual NUTS2016 country counts")
    # Classify each NUTS3 before weighting its fixed VALID AREA; never compare
    # national mean PRC against the threshold.
    result = x.groupby(["CNTR_CODE", "year"], as_index=False).agg(
        n_nuts3=("NUTS_ID", "size"),
        n_nuts3_at_or_below=("prc_at_or_below_threshold", "sum"),
        fixed_2014_valid_area_all_km2=("fixed_2014_valid_area_km2", "sum"),
        fixed_2014_valid_area_at_or_below_km2=("fixed_area_at_or_below_km2", "sum"),
        year_valid_area_all_km2=("valid_area_km2", "sum"),
        year_valid_area_at_or_below_km2=("year_area_at_or_below_km2", "sum"),
    )
    result["country"] = result.CNTR_CODE.map(
        lambda c: COUNTRY_NAMES.get(c, c)
    )
    result["pct_fixed_2014_area_at_or_below"] = (
        100 * result.fixed_2014_valid_area_at_or_below_km2
        / result.fixed_2014_valid_area_all_km2
    )
    result["pct_year_valid_area_at_or_below_sensitivity"] = (
        100 * result.year_valid_area_at_or_below_km2
        / result.year_valid_area_all_km2
    )
    result["pct_nuts3_at_or_below_equal_region"] = (
        100 * result.n_nuts3_at_or_below / result.n_nuts3
    )
    result["area_coverage_ratio_to_2014"] = (
        result.year_valid_area_all_km2 / result.fixed_2014_valid_area_all_km2
    )
    if (result.fixed_2014_valid_area_all_km2 <= 0).any():
        raise ValueError("Country has non-positive 2014 fixed valid area")
    if not result[METRIC].between(0., 100.).all():
        raise ValueError("Country percentage exceeds 0–100")
    if not result.pct_year_valid_area_at_or_below_sensitivity.between(0., 100.).all():
        raise ValueError("Invalid annual valid-area sensitivity")
    if result.groupby("CNTR_CODE").year.nunique().ne(len(YEARS)).any():
        raise ValueError("Missing country-year in locked PRC history")
    if result.groupby("CNTR_CODE").n_nuts3.nunique().ne(1).any():
        raise ValueError("Country NUTS3 cohort changed between years")
    if result.groupby("CNTR_CODE").fixed_2014_valid_area_all_km2.agg(
        lambda y: float(y.max() - y.min()) > 1e-5
    ).any():
        raise ValueError("Fixed national denominator changed over time")
    # Critical identity: disjoint country partitions reproduce the exact
    # full-Europe estimator numerators AND denominators in every year.
    country_sums = result.groupby("year")[
        ["fixed_2014_valid_area_all_km2",
         "fixed_2014_valid_area_at_or_below_km2", "n_nuts3"]
    ].sum().sort_index()
    continent_sums = x.groupby("year").agg(
        fixed_2014_valid_area_all_km2=("fixed_2014_valid_area_km2", "sum"),
        fixed_2014_valid_area_at_or_below_km2=("fixed_area_at_or_below_km2", "sum"),
        n_nuts3=("NUTS_ID", "size")
    ).sort_index()
    if not np.allclose(
        country_sums.to_numpy(float), continent_sums.to_numpy(float),
        rtol=1e-10, atol=1e-5
    ):
        raise ValueError("Country weighted totals do not reconcile to locked European study")
    if country_sums.n_nuts3.ne(expected_regions).any():
        raise ValueError("Country partition lost NUTS3 regions")
    return result.sort_values(["CNTR_CODE", "year"]).reset_index(drop=True)


def full_europe_reference(country_series: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct precisely the EXISTING European area series from countries."""
    sums = country_series.groupby("year", as_index=False).agg(
        fixed_2014_valid_area_all_km2=("fixed_2014_valid_area_all_km2", "sum"),
        fixed_2014_valid_area_at_or_below_km2=(
            "fixed_2014_valid_area_at_or_below_km2", "sum"
        ),
        n_nuts3=("n_nuts3", "sum"),
        n_nuts3_at_or_below=("n_nuts3_at_or_below", "sum"),
    )
    sums["pct_fixed_2014_area_at_or_below"] = (
        100 * sums.fixed_2014_valid_area_at_or_below_km2
        / sums.fixed_2014_valid_area_all_km2
    )
    return sums


def country_change_summary(country_series: pd.DataFrame) -> pd.DataFrame:
    first = country_series.loc[country_series.year.eq(2014)].set_index("CNTR_CODE")
    last = country_series.loc[country_series.year.eq(2024)].set_index("CNTR_CODE")
    if set(first.index) != set(last.index):
        raise ValueError("Country endpoints differ")
    out = first[[
        "country", "n_nuts3", "fixed_2014_valid_area_all_km2",
        METRIC
    ]].rename(columns={METRIC: "pct_area_2014"})
    out["pct_area_2024"] = last[METRIC]
    out["change_percentage_points_2014_2024"] = (
        out.pct_area_2024 - out.pct_area_2014
    )
    return out.reset_index().sort_values(
        "change_percentage_points_2014_2024", ascending=False
    )


def panel_axis(ax, frame: pd.DataFrame, label: str, colour: str,
               *, zoom: bool) -> None:
    ordered = frame.sort_values("year")
    x = ordered.year.to_numpy(int)
    values = ordered[METRIC].to_numpy(float)
    ax.plot(x, values, marker="o", markersize=3.7,
            color=colour, linewidth=1.7)
    ax.axhline(values[0], color="#777777", linestyle=":", linewidth=.9)
    ax.set_xlim(2013.6, 2024.4)
    ax.set_xticks((2014, 2016, 2018, 2020, 2022, 2024))
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.set_ylim(*zoomed_limits(values) if zoom else (0, 100))
    ax.grid(axis="y", alpha=.2)
    ax.set_title(
        f"{label}  ·  n={int(ordered.n_nuts3.iloc[0])}  "
        f"({values[-1]-values[0]:+.2f} pp)", loc="left", fontsize=10
    )


def zoomed_limits(values: np.ndarray) -> tuple[float, float]:
    """Country-specific bounds; only use with explicit zoomed-axis labeling."""
    data = np.asarray(values, dtype=float)
    if len(data) != len(YEARS) or not np.isfinite(data).all():
        raise ValueError("Invalid annual threshold series")
    width = max(float(np.ptp(data)), .15)
    pad = max(.1, .18 * width)
    bottom = max(0., float(data.min()) - pad)
    top = min(100., float(data.max()) + pad)
    if top - bottom < .2:
        bottom = max(0., bottom-.1)
        top = min(100., top+.1)
    if top <= bottom:
        bottom, top = (0., .25) if np.isclose(bottom, 0.) else (99.75, 100.)
    return bottom, top


def plot_country_grid(series: pd.DataFrame, europe: pd.DataFrame,
                      codes: list[str], output: Path, *,
                      zoom: bool, threshold: float, dpi: int) -> None:
    """Seven selected countries plus the Europe reference, 2×4 panels."""
    display = [*codes, "EUROPE"]
    cols = 4
    rows = int(np.ceil(len(display)/cols))
    fig, axes = plt.subplots(rows, cols, figsize=(16., rows*4.0),
                             squeeze=False, facecolor="white")
    for i, code in enumerate(display):
        ax = axes.flat[i]
        if code == "EUROPE":
            # Mimic selected-country frame fields required by panel_axis.
            region = europe.copy()
            colour, name = "#333333", "Europe · locked NUTS2016"
        else:
            region = series.loc[series.CNTR_CODE.eq(code)]
            colour, name = COLOURS.get(code, "#396999"), COUNTRY_NAMES.get(code, code)
        panel_axis(ax, region, name, colour, zoom=zoom)
        if i % cols == 0:
            ax.set_ylabel("Area at/below threshold (%)")
        if i >= (rows-1)*cols:
            ax.set_xlabel("Year")
    for i in range(len(display), rows*cols):
        axes.flat[i].axis("off")
    suffix = "zoomed" if zoom else "common_0_100_scale"
    fig.suptitle(
        f"Share of fixed 2014 valid area in NUTS3 with spider PRC ≤ {threshold:.3f}\n"
        "2014–2024 · illustrative DoP 0.23 scenario",
        fontsize=15, y=.996
    )
    fig.text(
        .5, .018,
        ("Each panel has its OWN ZOOMED y-axis. Do not compare absolute "
         "curve heights across countries."
         if zoom else "All panels use the identical 0–100% area scale.")
        + " NUTS3 means classified; not native-pixel fraction.",
        ha="center", fontsize=9
    )
    fig.tight_layout(rect=(.01, .045, 1., .94))
    stem = output / f"spider_prc_dop023_threshold_selected_countries_{suffix}"
    save_figure(fig, stem)
    if dpi:
        fig.savefig(stem.with_suffix(".png"), dpi=dpi, facecolor="white")
    plt.close(fig)


def plot_country_individual(series: pd.DataFrame, codes: list[str],
                            output: Path, *, threshold: float, dpi: int) -> None:
    for code in codes:
        frame = series.loc[series.CNTR_CODE.eq(code)].sort_values("year")
        values = frame[METRIC].to_numpy(float)
        fig, ax = plt.subplots(figsize=(9.0, 5.2), facecolor="white")
        fig.subplots_adjust(left=.12, right=.96, top=.82, bottom=.19)
        ax.plot(YEARS, values, color=COLOURS.get(code, "#396999"),
                linewidth=2., marker="o", markersize=5)
        ax.axhline(values[0], linestyle=":", linewidth=1, color="#666666",
                   label="2014 baseline")
        ax.set_xlim(2013.6, 2024.4)
        ax.set_xticks(YEARS)
        ax.set_ylim(*zoomed_limits(values))
        ax.grid(axis="y", alpha=.2)
        ax.set_xlabel("Year")
        ax.set_ylabel("Fixed 2014 study area at/below threshold (%)")
        ax.set_title(
            f"{COUNTRY_NAMES.get(code, code)}: spider PRC ≤ {threshold:.3f}"
            f" · n={int(frame.n_nuts3.iloc[0])} NUTS3",
            loc="left"
        )
        ax.legend(loc="best", frameon=False)
        fig.text(
            .12, .045, f"Indicative DoP=0.23; y-axis zoomed. "
            f"2014: {values[0]:.3f}% · 2024: {values[-1]:.3f}% · "
            f"Δ: {values[-1]-values[0]:+.3f} percentage points. "
            "NUTS3 mean PRC, fixed valid area.",
            fontsize=9
        )
        stem = output / f"spider_prc_dop023_threshold_country_{code}_2014_2024"
        save_figure(fig, stem)
        if dpi:
            fig.savefig(stem.with_suffix(".png"), dpi=dpi, facecolor="white")
        plt.close(fig)


def self_test() -> None:
    from plot_publication_spider_prc_dop023_threshold import (
        DOP_THRESHOLD, POLARIZATION_SENSITIVITY
    )
    threshold = dop_to_prc(DOP_THRESHOLD, POLARIZATION_SENSITIVITY)
    # A and B are Dutch; C is a UK region. Fixed areas 2+3+5 km².
    toy = pd.DataFrame({
        "NUTS_ID": ["NLA", "NLB", "UKC"] * len(YEARS),
        "year": np.repeat(YEARS, 3),
        "prc_mean": ([.2, .4, threshold] + [.4, .2, .4]*(len(YEARS)-1)),
        "valid_area_km2": ([2.,3.,5.] * len(YEARS)),
    })
    geometry = pd.DataFrame({
        "NUTS_ID": ["NLA", "NLB", "UKC"],
        "CNTR_CODE": ["NL", "NL", "UK"],
    })
    classified = annotate_regions(toy, threshold, expected_regions=3)
    by_country = country_area_series(classified, geometry, expected_regions=3)
    nl = by_country.loc[by_country.CNTR_CODE.eq("NL")].set_index("year")
    uk = by_country.loc[by_country.CNTR_CODE.eq("UK")].set_index("year")
    assert np.isclose(nl.loc[2014, METRIC], 40.)
    assert np.isclose(nl.loc[2015, METRIC], 60.)
    assert np.isclose(uk.loc[2014, METRIC], 100.)
    assert np.isclose(uk.loc[2015, METRIC], 0.)
    assert nl.n_nuts3.eq(2).all() and uk.n_nuts3.eq(1).all()
    europe = full_europe_reference(by_country).set_index("year")
    assert np.isclose(europe.loc[2014, METRIC], 70.)
    assert np.isclose(europe.loc[2015, METRIC], 30.)
    assert np.allclose(europe.fixed_2014_valid_area_all_km2, 10.)
    summary = country_change_summary(by_country)
    assert len(summary) == 2
    assert np.isclose(summary.set_index("CNTR_CODE").loc["NL",
        "change_percentage_points_2014_2024"], 20.)
    assert 0 <= zoomed_limits(np.array([0.]*len(YEARS)))[0] < zoomed_limits(
        np.array([0.]*len(YEARS)))[1]
    print("Country PRC DoP=0.23 threshold area and reconciliation test: PASS")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annual", type=Path, default=ANNUAL["spider"])
    p.add_argument("--endpoint-qa", type=Path, default=JOINED_ENDPOINT_QA)
    p.add_argument("--out", type=Path, default=OUTPUT)
    p.add_argument("--countries", nargs="+", default=list(COUNTRIES),
                   help="Country codes to plot (all available countries always in CSV)")
    p.add_argument("--dop-threshold", type=float, default=DOP_THRESHOLD)
    p.add_argument("--polarization-sensitivity", type=float,
                   default=POLARIZATION_SENSITIVITY)
    p.add_argument("--png-dpi", type=int, default=220)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    self_test()
    if args.self_test:
        return
    if not 0 <= args.png_dpi <= 600:
        p.error("PNG DPI must be between 0 and 600")
    if not 0 < args.dop_threshold < 1:
        p.error("DoP scenario must be strictly between 0 and 1")
    if not args.polarization_sensitivity > 1:
        p.error("Polarization sensitivity must exceed 1")
    country_codes = [str(c).upper() for c in args.countries]
    if not country_codes or len(country_codes) != len(set(country_codes)):
        p.error("Provide distinct country codes")
    threshold = dop_to_prc(args.dop_threshold, args.polarization_sensitivity)
    apply_figure_style()
    source, source_qa = read_spider_history(resolve(args.annual))
    endpoint_qa = independent_endpoint_qa(source, resolve(args.endpoint_qa))
    classified = annotate_regions(source, threshold)
    geometry = load_nuts(load_config())
    series = country_area_series(classified, geometry)
    available = set(series.CNTR_CODE)
    missing = set(country_codes) - available
    if missing:
        raise ValueError(
            f"Requested country codes absent from locked NUTS2016: {sorted(missing)}. "
            f"Available: {sorted(available)}"
        )
    europe = full_europe_reference(series)
    # For the full locked study area, this is identical to the previously
    # published European series. It cannot silently omit UK or other regions.
    if not europe.n_nuts3.eq(EXPECTED_REGIONS).all():
        raise ValueError("Europe reference changed its locked cohort")
    summary = country_change_summary(series)
    output = resolve(args.out)
    output.mkdir(parents=True, exist_ok=True)
    series.to_csv(
        output/"spider_prc_dop023_threshold_country_area_percent_2014_2024.csv",
        index=False
    )
    summary.to_csv(
        output/"spider_prc_dop023_threshold_country_2014_vs_2024.csv",
        index=False
    )
    summary.loc[summary.CNTR_CODE.isin(country_codes)].to_csv(
        output/"spider_prc_dop023_threshold_selected_country_summary.csv",
        index=False
    )
    plot_country_grid(series, europe, country_codes, output,
                      zoom=True, threshold=threshold, dpi=args.png_dpi)
    plot_country_grid(series, europe, country_codes, output,
                      zoom=False, threshold=threshold, dpi=args.png_dpi)
    plot_country_individual(
        series, country_codes, output, threshold=threshold, dpi=args.png_dpi
    )
    provenance = {
        "years": list(YEARS),
        "locked_spider_source": source_qa,
        "independent_endpoint_QA": endpoint_qa,
        "scenario_threshold_DOP": float(args.dop_threshold),
        "opponent_polarization_sensitivity": float(args.polarization_sensitivity),
        "threshold_PRC": float(threshold),
        "original_Europe_threshold_analysis_unchanged": True,
        "all_NUTS2016_country_codes": sorted(available),
        "selected_country_codes": country_codes,
        "selected_country_nuts3_count": {
            c: int(series.loc[series.CNTR_CODE.eq(c), "n_nuts3"].iloc[0])
            for c in country_codes
        },
        "country_boundaries_source": "original fixed NUTS2016 NUTS3 CNTR_CODE",
        "not_population_matched_cohort": True,
        "UK_included_from_full_1522_region_PRC_cohort": "UK" in available,
        "area_estimand": (
            "Share of each country's fixed 2014 valid pixel-covered area in "
            "entire NUTS2016 NUTS3 regions with annual region-mean spider PRC "
            "at or below the scenario threshold. NOT the proportion of "
            "native VIIRS pixels below threshold."
        ),
        "fixed_country_denominator_over_time": True,
        "country_sums_reconcile_to_full_Europe": True,
        "no_spatial_or_temporal_interpolation": True,
        "locked_2021_snow_product_reused": True,
        "input_data_or_calibration_modified": False,
        "individual_plots_y_axis": "zoomed to each country's observed range, labelled",
        "grid_comparison_scales": ["independently zoomed", "shared 0–100 percent"],
        "output_folder": str(output),
    }
    (output/"spider_prc_dop023_threshold_country_provenance.json").write_text(
        json.dumps(provenance, indent=2)+"\n", encoding="utf-8"
    )
    print("\nCountry PRC threshold: ", f"{threshold:.8f}")
    for c in country_codes:
        row = summary.loc[summary.CNTR_CODE.eq(c)].iloc[0]
        print(
            f"{c} {row.country}: n={int(row.n_nuts3)}; "
            f"2014 {row.pct_area_2014:.5f}%; "
            f"2024 {row.pct_area_2024:.5f}%; "
            f"change {row.change_percentage_points_2014_2024:+.5f} pp.",
            flush=True
        )
    print(f"Country plots, tables and provenance saved: {output}", flush=True)


if __name__ == "__main__":
    main()
