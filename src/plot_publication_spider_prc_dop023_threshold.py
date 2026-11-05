#!/usr/bin/env python3
"""Exploratory DoP=0.23 equivalent threshold for locked spider PRC, 2014–2024.

Produces 11 threshold-centred annual NUTS2016 PRC LEVEL maps and fixed-area
fraction of NUTS3 regions with mean PRC at/below the illustrative threshold.
The fraction is NOT area of native VIIRS pixels below threshold: the source
contains region means, and thresholding a mean is not pixel classification.

All science inputs are read only, including locked 2021 snow-corrected PRC.
Same clean 16:9 atlas geometry as the annual radiance/spider PRC maps.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd

from figure_style import apply_figure_style, save_figure
from project_config import PROJECT_ROOT as ROOT, load_config
from plot_publication_annual_nuts3_spider_prc_maps import (
    YEARS, EXPECTED_REGIONS, read_spider_history, independent_endpoint_qa,
    resolve, shared_linear_limits, map_geometry, EXTENT_LONLAT
)
from plot_publication_figure3_prc import ANNUAL, JOINED_ENDPOINT_QA

OUTPUT = ROOT / "figures/publication/annual_nuts3_spider_prc_dop023_threshold_2014_2024"
DOP_THRESHOLD = 0.23
POLARIZATION_SENSITIVITY = 10.0
PALETTE = "below_red_amber_above_green_teal"
# Below the centre: red through orange. Just above: green/teal.
THRESHOLD_CMAP = LinearSegmentedColormap.from_list(
    PALETTE,
    [(0., "#941d37"), (.28, "#d75d43"), (.499, "#efbd67"),
     (.5, "#f3d06b"), (.501, "#d9e799"), (.74, "#70bd89"),
     (1., "#196d76")],
    N=1024,
)


def dop_to_prc(dop: float, polarization_sensitivity: float = 10.) -> float:
    """Ideal opposing-channel log contrast at *optimal e-vector alignment*.

    PRC = ln(R_parallel/R_perp) =
          2 atanh(DoP * (Sp-1)/(Sp+1)).
    A DoP navigation-performance datum from OTHER TAXA is NOT a spider
    detectability threshold. This formula is a scenario, not new calibration.
    """
    d, sp = float(dop), float(polarization_sensitivity)
    if not np.isfinite(d) or not 0 <= d <= 1:
        raise ValueError("DoP must be finite and between 0 and 1")
    if not np.isfinite(sp) or sp <= 1:
        raise ValueError("Receptor polarization sensitivity must exceed 1")
    return float(2 * np.arctanh(d * (sp - 1) / (sp + 1)))


def threshold_palette(annual: pd.DataFrame, mapped_ids: set[str],
                      threshold: float) -> tuple[TwoSlopeNorm, dict]:
    """Use one threshold-centred scale for EVERY year; never clip observations."""
    low, high, original = shared_linear_limits(
        annual, mapped_ids, None, None
    )
    # Centre is always visible, even if all observed regional means end up
    # on one side. The colourbar is constant through time.
    low = min(low, max(0., threshold - .01))
    high = max(high, min(1., threshold + .01))
    if not low < threshold < high:
        raise ValueError("Threshold must be interior to the global PRC scale")
    norm = TwoSlopeNorm(vmin=low, vcenter=threshold, vmax=high)
    vals = annual.loc[annual.NUTS_ID.isin(mapped_ids), "prc_mean"]
    if vals.min() < low - 1e-12 or vals.max() > high + 1e-12:
        raise ValueError("Common threshold palette clips mapped data")
    return norm, {
        **original,
        "threshold_prc": threshold,
        "scale_vmin": low,
        "scale_vmax": high,
        "normalization": "TwoSlopeNorm, separately linear below/above threshold",
        "palette": PALETTE,
        "below_at_threshold": "red/orange/yellow at or below threshold",
        "above_threshold": "yellow-green through dark teal",
        "shared_norm_every_year": True,
        "yearly_rescaling": False,
        "data_clipped": False,
    }


def annotate_regions(annual: pd.DataFrame, threshold: float,
                     expected_regions: int = EXPECTED_REGIONS) -> pd.DataFrame:
    """Classify each *regional mean* PRC; freeze 2014 valid-area denominator."""
    need = {"NUTS_ID", "year", "prc_mean", "valid_area_km2"}
    if missing := need - set(annual):
        raise ValueError(f"Missing fields from locked annual PRC: {missing}")
    d = annual.copy()
    if d.duplicated(["year", "NUTS_ID"]).any():
        raise ValueError("Duplicated NUTS3/year")
    if sorted(d.year.unique().tolist()) != list(YEARS):
        raise ValueError("Expected 2014–2024 locked years")
    if d.groupby("year").NUTS_ID.nunique().ne(expected_regions).any():
        raise ValueError("Lost locked NUTS3 cohort in threshold classification")
    baseline = d.loc[d.year.eq(2014),
                     ["NUTS_ID", "valid_area_km2"]].rename(
        columns={"valid_area_km2": "fixed_2014_valid_area_km2"}
    )
    if len(baseline) != expected_regions:
        raise ValueError("Wrong 2014 reference area cohort")
    d = d.merge(baseline, on="NUTS_ID", validate="many_to_one")
    if not np.isfinite(d.prc_mean.to_numpy(float)).all():
        raise ValueError("PRC contains nonfinite values")
    if not np.isfinite(d.fixed_2014_valid_area_km2.to_numpy(float)).all():
        raise ValueError("Nonfinite reference area")
    if (d.fixed_2014_valid_area_km2 <= 0).any():
        raise ValueError("Invalid fixed 2014 mapped-pixel area")
    if (d.valid_area_km2 <= 0).any():
        raise ValueError("Invalid yearly valid area")
    d["prc_at_or_below_threshold"] = d.prc_mean.le(threshold)
    d["fixed_area_at_or_below_km2"] = np.where(
        d.prc_at_or_below_threshold, d.fixed_2014_valid_area_km2, 0.
    )
    d["year_area_at_or_below_km2"] = np.where(
        d.prc_at_or_below_threshold, d.valid_area_km2, 0.
    )
    d["prc_minus_threshold"] = d.prc_mean - threshold
    return d


def area_time_series(data: pd.DataFrame,
                     visible_ids: set[str],
                     expected_regions: int = EXPECTED_REGIONS) -> pd.DataFrame:
    """Primary: full 1,522 locked regions, constant 2014 valid-area weights.

    Also record the visible continental subset (membership by NUTS3 polygon
    intersection with the fixed display extent, not clipped pixel-level area)
    and a year-specific coverage sensitivity. No averaging of PRC before
    threshold classification: classify each NUTS3 mean individually.
    """
    if not visible_ids:
        raise ValueError("Empty continental NUTS3 display cohort")
    groups = [
        ("full_locked_nuts2016", data),
        ("continental_display_regions", data[data.NUTS_ID.isin(visible_ids)]),
    ]
    rows = []
    for label, group in groups:
        for year in YEARS:
            sample = group.loc[group.year.eq(year)]
            fixed_total = float(sample.fixed_2014_valid_area_km2.sum())
            fixed_low = float(sample.fixed_area_at_or_below_km2.sum())
            year_total = float(sample.valid_area_km2.sum())
            year_low = float(sample.year_area_at_or_below_km2.sum())
            if not fixed_total > 0 or not year_total > 0:
                raise ValueError(f"No valid denominator: {label} {year}")
            rows.append({
                "geography": label,
                "year": year,
                "n_nuts3": int(len(sample)),
                "n_nuts3_at_or_below": int(
                    sample.prc_at_or_below_threshold.sum()
                ),
                "pct_nuts3_at_or_below_equal_region": (
                    100 * sample.prc_at_or_below_threshold.mean()
                ),
                "fixed_2014_valid_area_all_km2": fixed_total,
                "fixed_2014_valid_area_at_or_below_km2": fixed_low,
                "pct_fixed_2014_area_at_or_below": 100 * fixed_low / fixed_total,
                "year_valid_area_all_km2": year_total,
                "year_valid_area_at_or_below_km2": year_low,
                "pct_year_valid_area_at_or_below_sensitivity":
                    100 * year_low / year_total,
                "area_coverage_ratio_to_2014": year_total / fixed_total,
            })
    result = pd.DataFrame(rows)
    full = result[result.geography.eq("full_locked_nuts2016")]
    if (full.n_nuts3 != expected_regions).any():
        raise ValueError("Full-study area time series lost expected NUTS3 cohort")
    if full.fixed_2014_valid_area_all_km2.max() - (
        full.fixed_2014_valid_area_all_km2.min()
    ) > 1e-5:
        raise ValueError("Area denominator changes despite frozen 2014 support")
    for metric in ("pct_fixed_2014_area_at_or_below",
                   "pct_year_valid_area_at_or_below_sensitivity"):
        if not result[metric].between(0., 100.).all():
            raise ValueError(f"Invalid percentage from area series: {metric}")
    return result


def make_year_map(geom, borders, bounds, values: pd.DataFrame,
                  year: int, norm: TwoSlopeNorm, threshold: float,
                  dop: float):
    g = geom.merge(
        values[["NUTS_ID", "prc_mean", "prc_at_or_below_threshold"]],
        on="NUTS_ID", how="left", validate="one_to_one"
    )
    present = g.prc_mean.notna()
    fig, ax = plt.subplots(figsize=(12.8, 7.2), facecolor="white")
    fig.subplots_adjust(left=.02, right=.84, bottom=.035, top=.89)
    g.plot(ax=ax, color="#dedede", edgecolor="none", linewidth=0)
    if present.any():
        g.loc[present].plot(
            ax=ax, column="prc_mean", norm=norm, cmap=THRESHOLD_CMAP,
            edgecolor="none", linewidth=0
        )
    borders.plot(ax=ax, color="#555555", linewidth=.33)
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.text(.035, .93, str(year), fontsize=26, weight="bold",
             ha="left", va="center")
    fig.text(.19, .931, "Spider PRC · illustrative DoP threshold",
             fontsize=17, ha="left", va="center")
    sm = ScalarMappable(norm=norm, cmap=THRESHOLD_CMAP)
    sm.set_array([])
    bar_ax = fig.add_axes([.865, .23, .025, .55])
    cb = fig.colorbar(sm, cax=bar_ax, extend="neither")
    raw_ticks = [
        norm.vmin, (norm.vmin + threshold)/2, threshold,
        (threshold + norm.vmax)/2, norm.vmax
    ]
    cb.set_ticks(raw_ticks)
    cb.ax.set_yticklabels(
        [f"{v:.3f}" for v in raw_ticks], fontsize=9
    )
    cb.ax.axhline(threshold, color="#272727", linewidth=1.8, zorder=10)
    cb.set_label("Spider PRC (dimensionless)", fontsize=11, labelpad=9)
    fig.text(.855, .145,
             f"DoP {dop:.2f} ≈ PRC {threshold:.3f}\n"
             "Warm: at/below\nGreen–teal: above",
             fontsize=9, ha="left", va="top")
    return fig, {
        "year": int(year),
        "source_nuts3": int(len(values)),
        "display_nuts3": int(present.sum()),
        "display_nuts3_at_or_below": int(
            g.loc[present, "prc_at_or_below_threshold"].sum()
        ),
        "display_nuts3_missing": int((~present).sum()),
        "shared_fixed_threshold_prc": float(threshold),
    }


def plot_time_series(data: pd.DataFrame, out: Path,
                     threshold: float, dop: float):
    whole = data[data.geography.eq("full_locked_nuts2016")]
    cont = data[data.geography.eq("continental_display_regions")]
    if len(whole) != len(YEARS) or len(cont) != len(YEARS):
        raise ValueError("Both geographies must contain 11 years")
    x = whole.year.to_numpy(int)
    y = whole.pct_fixed_2014_area_at_or_below.to_numpy(float)
    yc = cont.pct_fixed_2014_area_at_or_below.to_numpy(float)
    fig, ax = plt.subplots(figsize=(12.8, 6.4), facecolor="white")
    fig.subplots_adjust(left=.11, right=.96, top=.86, bottom=.20)
    ax.plot(x, y, color="#ba4046", marker="o", linewidth=2.2,
            markersize=6, label="Full locked NUTS2016 area (primary)")
    ax.plot(x, yc, color="#126b6e", marker="s", linewidth=1.4,
            markersize=4, linestyle="--",
            label="Continental display-region cohort (sensitivity)")
    ax.axhline(y[0], color="#777777", linewidth=.9,
               linestyle=":", label="2014 full-study baseline")
    ax.set_xlim(2013.7, 2024.3)
    ax.set_xticks(list(YEARS))
    bounds = np.r_[y, yc]
    span = max(float(np.ptp(bounds)), .15)
    padding = max(.1, span * .14)
    lo = max(0., float(bounds.min()) - padding)
    hi = min(100., float(bounds.max()) + padding)
    if hi <= lo:
        hi = min(100., lo + .5)
    ax.set_ylim(lo, hi)
    ax.set_ylabel("Study area at/below threshold (%)")
    ax.set_xlabel("Year")
    ax.set_title(
        f"NUTS3 mean spider PRC ≤ {threshold:.3f} "
        f"(indicative DoP {dop:.2f}) · 2014–2024",
        loc="left"
    )
    ax.grid(axis="y", alpha=.2)
    ax.legend(loc="best", frameon=False)
    fig.text(
        .11, .045,
        "Fixed 2014 valid-pixel area weights; entire NUTS3 classified by its "
        "mean PRC. Not pixel-level area. Y-axis zoomed; inspect CSV for values.",
        fontsize=9
    )
    save_figure(fig, out / "spider_prc_dop023_threshold_area_percent_time_series")
    fig.savefig(
        out / "spider_prc_dop023_threshold_area_percent_time_series.png",
        dpi=220, facecolor="white"
    )
    plt.close(fig)


def self_test() -> None:
    threshold = dop_to_prc(.23, 10.)
    assert np.isclose(threshold, .3809031327477229)
    assert np.isclose(
        (11/9) * np.tanh(threshold/2), .23
    )
    assert dop_to_prc(.11, 10.) < threshold < dop_to_prc(.32, 10.)
    toy = pd.DataFrame({
        "NUTS_ID": ["A", "B", "C"] * 11,
        "year": np.repeat(YEARS, 3),
        "prc_mean": ([.20, .40, threshold] +
                     [.40, .20, .40] * 10),
        "valid_area_km2": ([2., 3., 5.] + [2., 3., 5.] * 10),
    })
    # Test the ACTUAL production classifier and area calculations, using
    # a small synthetic 11-year dataset rather than only copying its formula.
    classified = annotate_regions(toy, threshold, expected_regions=3)
    series = area_time_series(classified, {"A", "B"}, expected_regions=3)
    full = series.loc[series.geography.eq("full_locked_nuts2016")].set_index("year")
    cont = series.loc[series.geography.eq("continental_display_regions")].set_index("year")
    assert np.isclose(full.loc[2014, "pct_fixed_2014_area_at_or_below"], 70.)
    assert np.isclose(full.loc[2015, "pct_fixed_2014_area_at_or_below"], 30.)
    assert np.isclose(cont.loc[2014, "pct_fixed_2014_area_at_or_below"], 40.)
    assert np.isclose(cont.loc[2015, "pct_fixed_2014_area_at_or_below"], 60.)
    assert (full.fixed_2014_valid_area_all_km2 == 10.).all()
    assert classified.loc[
        classified.year.eq(2014) & classified.NUTS_ID.eq("C"),
        "prc_at_or_below_threshold"
    ].iloc[0]
    small = pd.DataFrame({
        "NUTS_ID": ["A", "B", "C", "A", "B", "C"],
        "prc_mean": [.16, .42, .39, .18, .41, .40],
    })
    norm, meta = threshold_palette(small, {"A","B","C"}, threshold)
    assert np.isclose(norm(threshold), .5)
    assert meta["shared_norm_every_year"] and not meta["data_clipped"]
    print("DoP-to-PRC scenario and fixed-area classification QA: PASS",
          flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annual", type=Path, default=ANNUAL["spider"])
    parser.add_argument("--endpoint-qa", type=Path, default=JOINED_ENDPOINT_QA)
    parser.add_argument("--out", type=Path, default=OUTPUT)
    parser.add_argument("--dop-threshold", type=float, default=DOP_THRESHOLD)
    parser.add_argument("--polarization-sensitivity", type=float,
                        default=POLARIZATION_SENSITIVITY)
    parser.add_argument("--display-simplify-m", type=float, default=150.)
    parser.add_argument("--png-dpi", type=int, default=240)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    self_test()
    if args.self_test:
        return
    if not 0 < args.dop_threshold < 1:
        parser.error("--dop-threshold must be between 0 and 1")
    if not args.polarization_sensitivity > 1:
        parser.error("--polarization-sensitivity must exceed 1")
    if not 0 <= args.png_dpi <= 600:
        parser.error("--png-dpi must be between 0 and 600")
    threshold = dop_to_prc(
        args.dop_threshold, args.polarization_sensitivity
    )
    apply_figure_style()
    geom, bounds, full_region_ids = map_geometry(
        load_config(), args.display_simplify_m
    )
    annual, source_qa = read_spider_history(resolve(args.annual))
    if set(annual.NUTS_ID) != full_region_ids:
        raise ValueError("Locked spider PRC IDs differ from fixed NUTS2016 geometry")
    endpoint_qa = independent_endpoint_qa(
        annual, resolve(args.endpoint_qa)
    )
    annotated = annotate_regions(annual, threshold)
    visible_ids = set(geom.NUTS_ID)
    series = area_time_series(annotated, visible_ids)
    norm, palette_meta = threshold_palette(
        annotated, visible_ids, threshold
    )
    output = resolve(args.out)
    output.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(
        output / "annual_nuts3_spider_prc_threshold_classification.csv",
        index=False
    )
    series.to_csv(
        output / "spider_prc_dop023_threshold_area_percent_2014_2024.csv",
        index=False
    )
    plot_time_series(series, output, threshold, args.dop_threshold)
    borders = geom.dissolve(by="CNTR_CODE").boundary
    figure_info = []
    deck_path = output / "all_years_spider_prc_dop023_threshold_2014_2024.pdf"
    with PdfPages(deck_path) as deck:
        for year in YEARS:
            region = annotated[annotated.year.eq(year)]
            fig, info = make_year_map(
                geom, borders, bounds, region, year,
                norm, threshold, args.dop_threshold
            )
            stem = output / f"spider_prc_dop023_threshold_{year}"
            save_figure(fig, stem, map_figure=True)
            if args.png_dpi:
                fig.savefig(
                    stem.with_suffix(".png"), format="png",
                    dpi=args.png_dpi, facecolor="white"
                )
            deck.savefig(fig, facecolor="white")
            plt.close(fig)
            figure_info.append(info)
            primary = series.loc[
                (series.geography.eq("full_locked_nuts2016"))
                & series.year.eq(year),
                "pct_fixed_2014_area_at_or_below"
            ].iloc[0]
            print(f"{year}: {primary:.5f}% of fixed study area "
                  f"in NUTS3 mean PRC <= {threshold:.5f}", flush=True)
    full_series = series[series.geography.eq("full_locked_nuts2016")]
    manifest = {
        "years": list(YEARS),
        "locked_spider_source": source_qa,
        "endpoint_QA": endpoint_qa,
        "dop_scenario": float(args.dop_threshold),
        "assumed_spider_receptor_polarization_sensitivity":
            float(args.polarization_sensitivity),
        "receptor_alignment_assumption":
            "optimally aligned e-vector with orthogonal opponent receptor pair",
        "scenario_equation":
            "threshold_PRC=2*atanh(DoP*(Sp-1)/(Sp+1))",
        "threshold_prc": float(threshold),
        "note_dop_023":
            "Dung beetle weaker orientation below 0.23 DoP is not a "
            "validated Drassodes spider sensory loss or behaviour threshold",
        "threshold_quantity": "regional mean spider PRC, not native-pixel PRC",
        "area_estimator":
            "Classify each full fixed NUTS3 region by its annual pixel-derived "
            "mean PRC <= threshold; weight by its fixed 2014 valid-pixel km2. "
            "Not exact fraction of pixels below threshold.",
        "primary_geography": "all 1,522 locked NUTS2016 regions",
        "secondary_geography":
            "NUTS2016 regions intersecting the fixed continental display "
            "extent; whole NUTS3 reference area weights, not clipped polygon areas",
        "region_fraction_not_primary": True,
        "annual_area_weight_2014_reference_year": True,
        "sensitivity_year_specific_valid_area":
            "present as CSV column only",
        "fixed_study_area_km2":
            float(full_series.fixed_2014_valid_area_all_km2.iloc[0]),
        "n_regions_full": EXPECTED_REGIONS,
        "n_regions_display": len(visible_ids),
        "display_extent_lonlat": list(EXTENT_LONLAT),
        "map_16_9_layout_same_as_radiance_atlas": True,
        "palette": palette_meta,
        "map_info": figure_info,
        "annual_values_modified": False,
        "snow_correction_2021_already_in_locked_source": True,
        "no_behavioural_threshold_validation": True,
        "output": str(output),
    }
    (output / "spider_prc_dop023_threshold_provenance.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print("\nThreshold PRC = ", f"{threshold:.8f}")
    print("Study area 2014: ",
          f"{full_series.pct_fixed_2014_area_at_or_below.iloc[0]:.6f}%")
    print("Study area 2024: ",
          f"{full_series.pct_fixed_2014_area_at_or_below.iloc[-1]:.6f}%")
    print(f"Maps, area series and metadata: {output}", flush=True)


if __name__ == "__main__":
    main()
