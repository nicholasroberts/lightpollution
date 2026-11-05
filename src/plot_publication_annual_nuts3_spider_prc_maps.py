#!/usr/bin/env python3
"""2014–2024 spider polarization receptor contrast NUTS3 talk/publication atlas.

Match the existing annual VIIRS radiance maps: fixed NUTS2016 boundaries and
continental extent, identical 16:9 slide framing and ONE common colour scale
for all 11 years. Show annual PRC LEVELS, not PRC change from 2014.

READ ONLY: final locked native-pixel-derived spider PRC NUTS2016 annual panel.
This already includes the approved 2021 snow correction. No radiance/PRC
recalibration, native-pixel transformation, area aggregation or interpolation.
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
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

from figure_style import apply_figure_style, save_figure
from project_config import PROJECT_ROOT as ROOT, load_config
from plot_publication_annual_nuts3_radiance_maps import (
    map_geometry, EXTENT_LONLAT
)
from plot_publication_figure3_prc import (
    ANNUAL, JOINED_ENDPOINT_QA, EXPECTED_REGIONS, read_locked_panel
)

YEARS = tuple(range(2014, 2025))
OUTPUT = ROOT / "figures/publication/annual_nuts3_spider_prc_2014_2024"
CMAP = "viridis"  # Dark = lower available PRC, yellow = higher PRC.
SCALE_STEP = 0.01  # Presentation-only rounding of pooled PRC extrema.


def resolve(path: Path | str) -> Path:
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def read_spider_history(path: Path) -> tuple[pd.DataFrame, dict]:
    """All annual values originate from the existing locked spider PRC panel."""
    d = read_locked_panel(path, "spider")
    if len(d) != EXPECTED_REGIONS * len(YEARS):
        raise ValueError("Not the complete 1,522 × 11 spider PRC panel")
    if sorted(d.year.unique()) != list(YEARS):
        raise ValueError("Unexpected spider PRC years")
    if not d.NUTS_ID.nunique() == EXPECTED_REGIONS:
        raise ValueError("Wrong frozen spider NUTS3 cohort")
    ids = [set(d.loc[d.year.eq(y), "NUTS_ID"]) for y in YEARS]
    if any(regions != ids[0] for regions in ids):
        raise ValueError("Unequal NUTS3 support among annual spider PRC maps")
    d = d.sort_values(["year", "NUTS_ID"]).reset_index(drop=True)

    # Cross-check the locked historical source against its own cumulative
    # difference. No new PRC response modelling is performed.
    base = d.loc[d.year.eq(2014)].set_index("NUTS_ID").prc_mean
    d["source_endpoint_difference_from_2014"] = (
        d.prc_mean - d.NUTS_ID.map(base)
    )
    for year in (2014, 2024):
        sample = d.loc[d.year.eq(year)]
        if not np.allclose(
            sample.source_endpoint_difference_from_2014,
            sample.cumulative_prc_change_absolute,
            atol=2e-7, rtol=1e-5
        ):
            raise ValueError(
                f"{year} stored cumulative spider PRC disagrees with endpoints"
            )
    d["locked_2021_snow_year"] = d.year.eq(2021)
    return d, {
        "source": str(path),
        "n_regions_each_year": int(EXPECTED_REGIONS),
        "n_years": len(YEARS),
        "source_snow_status_values": sorted(d.snow_status.unique().tolist()),
        "2021_from_final_snow_corrected_spider_panel": True,
        "2014_and_2024_stored_cumulative_change_consistent": True,
    }


def independent_endpoint_qa(annual: pd.DataFrame,
                            endpoint_file: Path) -> dict:
    """Check the two plotted endpoints against the existing Figure 3 source."""
    if not endpoint_file.exists():
        return {
            "checked": False,
            "reason": "Optional previously joined spider/radiance endpoint CSV absent"
        }
    other = pd.read_csv(endpoint_file, dtype={"NUTS_ID": str})
    required = {
        "NUTS_ID", "spider_prc_2014", "spider_prc_2024",
        "spider_prc_change_abs"
    }
    if required - set(other):
        raise ValueError("Joined source missing spider endpoint QA columns")
    if len(other) != EXPECTED_REGIONS or not other.NUTS_ID.is_unique:
        raise ValueError("Joined spider endpoint QA has changed region cohort")
    results = {}
    for year in (2014, 2024):
        col = f"spider_prc_{year}"
        current = annual.loc[annual.year.eq(year),
                             ["NUTS_ID", "prc_mean"]]
        matched = current.merge(
            other[["NUTS_ID", col]], on="NUTS_ID",
            validate="one_to_one"
        )
        if len(matched) != EXPECTED_REGIONS:
            raise ValueError(f"Unmatched {year} spider PRC endpoints")
        err = np.abs(matched.prc_mean.to_numpy(float) -
                     matched[col].to_numpy(float))
        if not np.allclose(
            matched.prc_mean, matched[col], atol=2e-6, rtol=1e-5
        ):
            raise ValueError(f"Locked spider {year} values differ from joined source")
        results[str(year)] = {"n": len(matched), "max_abs_diff": float(err.max())}
    return {"checked": True, "source": str(endpoint_file), "years": results}


def shared_linear_limits(d: pd.DataFrame, mapped_ids: set[str],
                         override_low: float | None,
                         override_high: float | None) -> tuple[float, float, dict]:
    """A single finite, non-clipping colour scale for all continental years."""
    values = d.loc[d.NUTS_ID.isin(mapped_ids), "prc_mean"].to_numpy(float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError("Missing/nonfinite mapped annual PRC values")
    observed_low, observed_high = float(values.min()), float(values.max())
    if not (0 <= observed_low < observed_high <= 1):
        raise ValueError("Invalid pooled biological PRC range")
    if (override_low is None) != (override_high is None):
        raise ValueError("Provide both --vmin and --vmax, or neither")
    if override_low is None:
        low = max(0., np.floor(observed_low / SCALE_STEP) * SCALE_STEP)
        high = min(1., np.ceil(observed_high / SCALE_STEP) * SCALE_STEP)
    else:
        low, high = float(override_low), float(override_high)
    if not (np.isfinite(low) and np.isfinite(high) and
            0 <= low < high <= 1):
        raise ValueError("Colour limits must obey 0 <= vmin < vmax <= 1")
    tolerance = 1e-12
    if low > observed_low + tolerance or high < observed_high - tolerance:
        raise ValueError(
            f"Common scale {low:g}–{high:g} would clip observed "
            f"continental spider PRC range {observed_low:.6g}–{observed_high:.6g}"
        )
    return float(low), float(high), {
        "observed_min_continental": observed_low,
        "observed_max_continental": observed_high,
        "fixed_vmin_all_years": float(low),
        "fixed_vmax_all_years": float(high),
        "colour_normalisation": "linear, absolute annual spider PRC",
        "scale_rounding_step": SCALE_STEP if override_low is None else None,
        "data_clipped": False,
    }


def make_map(geometry, borders, bounds: tuple, values: pd.DataFrame,
             year: int, norm: Normalize) -> tuple[plt.Figure, dict]:
    g = geometry.merge(
        values[["NUTS_ID", "prc_mean"]], on="NUTS_ID",
        how="left", validate="one_to_one"
    )
    valid = g.prc_mean.notna()
    fig, ax = plt.subplots(figsize=(12.8, 7.2), facecolor="white")
    fig.subplots_adjust(left=.02, right=.84, bottom=.035, top=.89)
    g.plot(ax=ax, color="#dedede", edgecolor="none", linewidth=0)
    if valid.any():
        g.loc[valid].plot(
            ax=ax, column="prc_mean", cmap=CMAP, norm=norm,
            edgecolor="none", linewidth=0
        )
    borders.plot(ax=ax, color="#5b5b5b", linewidth=.33)
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.text(.035, .93, str(year), fontsize=26, weight="bold",
             ha="left", va="center")
    fig.text(.19, .931, "Spider polarization receptor contrast",
             fontsize=17, ha="left", va="center")
    sm = ScalarMappable(norm=norm, cmap=CMAP)
    sm.set_array([])
    cb_ax = fig.add_axes([.865, .23, .025, .55])
    cb = fig.colorbar(sm, cax=cb_ax, extend="neither")
    ticks = np.linspace(norm.vmin, norm.vmax, 7)
    cb.set_ticks(ticks)
    cb.ax.set_yticklabels([f"{t:.3f}" for t in ticks], fontsize=9)
    cb.set_label("Spider PRC (dimensionless)", fontsize=11, labelpad=9)
    return fig, {
        "year": year,
        "n_source_nuts3": int(len(values)),
        "n_mapped_polygons_with_prc": int(valid.sum()),
        "n_mapped_polygons_without_prc": int((~valid).sum()),
        "mapped_min_prc": (
            float(g.loc[valid, "prc_mean"].min()) if valid.any() else None
        ),
        "mapped_max_prc": (
            float(g.loc[valid, "prc_mean"].max()) if valid.any() else None
        ),
        "locked_2021_snow_year": year == 2021,
    }


def self_test() -> None:
    assert YEARS == tuple(range(2014, 2025))
    trial = pd.DataFrame({
        "NUTS_ID": ["A", "B", "A", "B"],
        "year": [2014, 2014, 2015, 2015],
        "prc_mean": [.14, .46, .16, .43],
    })
    low, high, meta = shared_linear_limits(
        trial, {"A", "B"}, None, None
    )
    assert np.isclose(low, .14)
    assert np.isclose(high, .46)
    assert not meta["data_clipped"]
    assert np.isclose(shared_linear_limits(
        trial, {"A", "B"}, .1, .5
    )[0], .1)
    for lo, hi in ((.15, .5), (.1, .45)):
        try:
            shared_linear_limits(trial, {"A", "B"}, lo, hi)
        except ValueError:
            pass
        else:
            raise AssertionError("Clipped pooled PRC value accepted")
    norm = Normalize(vmin=low, vmax=high, clip=False)
    assert np.isclose(norm(low), 0.) and np.isclose(norm(high), 1.)
    assert len(trial.groupby("year")) == 2
    print("Locked annual spider PRC shared linear scale self-test: PASS",
          flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annual", type=Path, default=ANNUAL["spider"])
    parser.add_argument("--endpoint-qa", type=Path, default=JOINED_ENDPOINT_QA)
    parser.add_argument("--out", type=Path, default=OUTPUT)
    parser.add_argument("--display-simplify-m", type=float, default=150.)
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--png-dpi", type=int, default=240)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    self_test()
    if args.self_test:
        return
    if not 0 <= args.png_dpi <= 600:
        parser.error("--png-dpi must be 0..600")
    apply_figure_style()
    geometry, bounds, region_ids = map_geometry(
        load_config(), args.display_simplify_m
    )
    annual, source_qa = read_spider_history(resolve(args.annual))
    if set(annual.NUTS_ID) != region_ids:
        raise ValueError("Spider PRC and mapped NUTS2016 geography differ")
    endpoints = independent_endpoint_qa(annual, resolve(args.endpoint_qa))
    mapped_ids = set(geometry.NUTS_ID)
    low, high, scale_qa = shared_linear_limits(
        annual, mapped_ids, args.vmin, args.vmax
    )
    norm = Normalize(vmin=low, vmax=high, clip=False)
    borders = geometry.dissolve(by="CNTR_CODE").boundary
    out = resolve(args.out)
    out.mkdir(parents=True, exist_ok=True)
    # Map-ready copy of the untouched source values, including year and
    # locked 2021 status; original scientific table remains unchanged.
    annual.to_csv(
        out / "annual_nuts3_spider_prc_2014_2024_map_values.csv",
        index=False
    )
    maps = []
    deck = out / "all_years_nuts3_spider_prc_2014_2024.pdf"
    print(
        f"Rendering 2014–2024 spider PRC on {EXPECTED_REGIONS:,} NUTS3. "
        f"Fixed linear PRC scale {low:.3f}–{high:.3f}; "
        f"observed {scale_qa['observed_min_continental']:.6f}–"
        f"{scale_qa['observed_max_continental']:.6f}.",
        flush=True
    )
    with PdfPages(deck) as document:
        for year in YEARS:
            values = annual.loc[annual.year.eq(year)].copy()
            if len(values) != EXPECTED_REGIONS:
                raise ValueError(f"Spider PRC {year} cohort unexpectedly changed")
            fig, meta = make_map(
                geometry, borders, bounds, values, year, norm
            )
            stem = out / f"spider_prc_nuts3_{year}"
            save_figure(fig, stem, map_figure=True)
            if args.png_dpi:
                fig.savefig(
                    stem.with_suffix(".png"), format="png",
                    dpi=args.png_dpi, facecolor="white"
                )
            document.savefig(fig, facecolor="white")
            plt.close(fig)
            maps.append(meta)
            print(
                f"  {year}: {meta['n_mapped_polygons_with_prc']} coloured, "
                f"{meta['n_mapped_polygons_without_prc']} grey "
                f"(locked-snow year: {year == 2021})",
                flush=True
            )
    summary = {
        "years": list(YEARS),
        "n_nuts3_locked_cohort": EXPECTED_REGIONS,
        "n_rows_annual_source": len(annual),
        "source_qa": source_qa,
        "independent_endpoint_qa": endpoints,
        "geography_extent_lonlat": list(EXTENT_LONLAT),
        "map_display_simplification_epsg3035_m": args.display_simplify_m,
        "map_presentation": "same 16:9 year/title, geography and colourbar positions as annual radiance maps; no footer",
        "colour_map": CMAP,
        "shared_linear_scale": scale_qa,
        "yearly_maps": maps,
        "map_quantity": "annual level of native-pixel-derived spider NUTS3 mean PRC; NOT annual delta or percentage",
        "locked_2021_snow_corrected_source_not_reprocessed": True,
        "fixed_nuts2016_no_region_imputation": True,
        "original_scientific_data_or_calibration_modified": False,
        "outputs": ["svg", "pdf", "png"] if args.png_dpi else ["svg", "pdf"],
        "combined_11_page_pdf": str(deck),
    }
    (out / "annual_spider_prc_map_provenance.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"\n11 annual spider PRC maps, shared scale {low:.3f}–{high:.3f}, "
          f"saved in: {out}", flush=True)
    print(f"11-page PDF deck: {deck}", flush=True)


if __name__ == "__main__":
    main()
