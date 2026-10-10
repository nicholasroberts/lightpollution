#!/usr/bin/env python3
"""Plot paired annual radiance changes by fixed density band and habitability.

This stage is designed as a cleaner descriptive companion to the mixed-effects
density-temporal model.

The earlier exploratory plots calculated year-on-year change in a group mean
whose membership could change between years.  That is undesirable because a
large apparent annual change can then reflect changing group composition rather
than radiance change within the same NUTS3 regions.

This implementation therefore:

1. restricts the first analysis to fixed NUTS2016 geography (2013-2020);
2. optionally requires a complete eight-year region panel;
3. assigns each NUTS3 region ONE fixed population-density band from its median
   population density across the panel;
4. assigns each NUTS3 region ONE fixed H category from static H_core;
5. calculates radiance change WITHIN REGION between consecutive years first;
6. only then summarizes those paired changes within density x H groups.

The primary group estimate is calculated on the same log10 scale used by the
main temporal models:

    delta_log10_R_it = log10(R_it) - log10(R_i,t-1)

    group_pct_change = 100 * (10 ** mean(delta_log10_R_it) - 1)

This is the geometric-mean proportional change across paired regions.  Arithmetic
mean and median individual percentage changes are also retained for diagnostics.

Habitability groups are, by default, quintiles WITHIN each population-density
band.  This produces balanced H strata for comparing habitability effects within
a density class and avoids tiny high-H/low-density cells caused by applying one
set of Europe-wide H quintiles.  Exact H ranges are written to a cut-point table.

IMPORTANT CALIBRATION NOTE
--------------------------
The current branch still contains a provisional five-country post-2017
calibration correction.  This stage reads the configured radiance column and
will therefore need to be rerun after the production global VIIRS calibration
correction is fixed.  No raw VIIRS extraction is required.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from project_config import configured_path, load_config


def require_columns(data: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [c for c in columns if c not in data.columns]
    if missing:
        raise KeyError(f"{label} is missing required columns: {missing}")


def load_and_prepare(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = config["analysis"]["paired_density_h_trends"]
    processed = configured_path(config, "processed")

    merged_path = processed / config["outputs"]["merged"]
    habit_path = processed / config["outputs"]["nuts3_habitability"]

    merged = pd.read_parquet(merged_path).copy()
    habit = pd.read_parquet(habit_path).copy()

    radiance_col = str(cfg["radiance_column"])
    require_columns(
        merged,
        [
            "year",
            "NUTS_ID",
            "CNTR_CODE",
            "nuts_release",
            "population_density",
            radiance_col,
        ],
        "merged radiance/population table",
    )
    require_columns(
        habit,
        [
            "nuts_release",
            "NUTS_ID",
            "habitable_fraction_core",
        ],
        "habitability table",
    )

    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    nuts_release = int(cfg["nuts_release"])

    data = merged[
        merged["year"].between(first_year, last_year, inclusive="both")
        & (merged["nuts_release"] == nuts_release)
    ].copy()

    h = (
        habit[
            habit["nuts_release"] == nuts_release
        ][
            [
                "nuts_release",
                "NUTS_ID",
                "habitable_fraction_core",
            ]
        ]
        .drop_duplicates(["nuts_release", "NUTS_ID"])
        .rename(columns={"habitable_fraction_core": "H_core"})
    )

    data = data.merge(
        h,
        on=["nuts_release", "NUTS_ID"],
        how="inner",
        validate="many_to_one",
    )

    data[radiance_col] = pd.to_numeric(data[radiance_col], errors="coerce")
    data["population_density"] = pd.to_numeric(
        data["population_density"],
        errors="coerce",
    )
    data["H_core"] = pd.to_numeric(data["H_core"], errors="coerce")

    valid = (
        np.isfinite(data[radiance_col])
        & np.isfinite(data["population_density"])
        & np.isfinite(data["H_core"])
        & (data[radiance_col] > 0)
        & (data["population_density"] > 0)
        & (data["H_core"] >= 0)
        & (data["H_core"] <= 1)
        & data["NUTS_ID"].notna()
        & data["CNTR_CODE"].notna()
    )
    data = data.loc[valid].copy()

    expected_years = last_year - first_year + 1
    panel_counts = data.groupby("NUTS_ID")["year"].nunique()

    if bool(cfg["require_complete_panel"]):
        keep_ids = panel_counts[panel_counts == expected_years].index
    else:
        keep_ids = panel_counts[
            panel_counts >= int(cfg["minimum_panel_years"])
        ].index

    data = data[data["NUTS_ID"].isin(keep_ids)].copy()

    # Fixed region-level reference values.
    region = (
        data.groupby(["NUTS_ID", "CNTR_CODE"], as_index=False)
        .agg(
            reference_population_density=("population_density", "median"),
            H_core=("H_core", "first"),
            n_years=("year", "nunique"),
        )
    )

    edges = (
        [-np.inf]
        + [float(x) for x in cfg["density_band_breaks"]]
        + [np.inf]
    )
    labels = list(cfg["density_band_labels"])

    region["density_band"] = pd.cut(
        region["reference_population_density"],
        bins=edges,
        labels=labels,
        right=False,
        ordered=True,
    )

    region = region[region["density_band"].notna()].copy()

    return data, region


def assign_h_groups(
    region: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assign fixed H categories and return region assignments + cutpoints."""
    cfg = config["analysis"]["paired_density_h_trends"]
    n_groups = int(cfg["h_quantile_groups"])
    grouping = str(cfg["h_grouping"])

    if n_groups < 2:
        raise ValueError("h_quantile_groups must be >= 2")

    region = region.copy()
    region["H_group"] = pd.NA

    cut_rows: list[dict] = []

    if grouping == "within_density_band_quantiles":
        grouped = region.groupby("density_band", observed=True)
    elif grouping == "global_quantiles":
        grouped = [("all_density_bands", region)]
    else:
        raise ValueError(
            "analysis.paired_density_h_trends.h_grouping must be "
            "'within_density_band_quantiles' or 'global_quantiles'"
        )

    global_assignment = None

    for group_name, g in grouped:
        g = g.copy()

        # qcut on rank rather than H itself guarantees equal-count groups even
        # if several regions have nearly identical H values close to 1.
        rank = g["H_core"].rank(method="first")
        q = pd.qcut(
            rank,
            q=n_groups,
            labels=False,
            duplicates="raise",
        )
        q = q.astype(int) + 1

        labels = [
            (
                f"H Q{k} (lowest)"
                if k == 1
                else f"H Q{k} (highest)"
                if k == n_groups
                else f"H Q{k}"
            )
            for k in range(1, n_groups + 1)
        ]
        mapping = {k: labels[k - 1] for k in range(1, n_groups + 1)}

        assignment = pd.Series(
            q.map(mapping).to_numpy(),
            index=g.index,
            dtype="object",
        )

        if grouping == "within_density_band_quantiles":
            region.loc[g.index, "H_group"] = assignment
            density_name = str(group_name)
        else:
            global_assignment = assignment
            density_name = "all_density_bands"

        for k in range(1, n_groups + 1):
            members = g[q.to_numpy() == k]
            cut_rows.append(
                {
                    "density_band_for_H_grouping": density_name,
                    "H_group_number": k,
                    "H_group": mapping[k],
                    "n_regions": int(len(members)),
                    "H_min": float(members["H_core"].min()),
                    "H_median": float(members["H_core"].median()),
                    "H_max": float(members["H_core"].max()),
                }
            )

    if grouping == "global_quantiles":
        region["H_group"] = global_assignment

    labels = [
        (
            f"H Q{k} (lowest)"
            if k == 1
            else f"H Q{k} (highest)"
            if k == n_groups
            else f"H Q{k}"
        )
        for k in range(1, n_groups + 1)
    ]
    region["H_group"] = pd.Categorical(
        region["H_group"],
        categories=labels,
        ordered=True,
    )

    return region, pd.DataFrame(cut_rows)


def make_paired_changes(
    data: pd.DataFrame,
    region: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    cfg = config["analysis"]["paired_density_h_trends"]
    radiance_col = str(cfg["radiance_column"])

    d = data.merge(
        region[
            [
                "NUTS_ID",
                "reference_population_density",
                "density_band",
                "H_core",
                "H_group",
            ]
        ],
        on="NUTS_ID",
        how="inner",
        validate="many_to_one",
        suffixes=("", "_region"),
    )

    # Use region-level H from the assignment table.
    if "H_core_region" in d.columns:
        d["H_core"] = d["H_core_region"]
        d = d.drop(columns=["H_core_region"])

    d = d.sort_values(["NUTS_ID", "year"]).copy()
    d["log10_radiance"] = np.log10(d[radiance_col].astype(float))

    d["previous_year"] = d.groupby("NUTS_ID")["year"].shift(1)
    d["previous_radiance"] = d.groupby("NUTS_ID")[radiance_col].shift(1)
    d["previous_log10_radiance"] = (
        d.groupby("NUTS_ID")["log10_radiance"].shift(1)
    )

    d["delta_log10_radiance"] = (
        d["log10_radiance"] - d["previous_log10_radiance"]
    )
    d["individual_pct_change"] = 100.0 * (
        d[radiance_col] / d["previous_radiance"] - 1.0
    )

    # Retain only genuine consecutive-year comparisons of the same NUTS3.
    d = d[
        (d["year"] - d["previous_year"]) == 1
    ].copy()

    d = d.rename(columns={"year": "transition_end_year"})
    d["transition_start_year"] = d["transition_end_year"] - 1

    columns = [
        "transition_start_year",
        "transition_end_year",
        "NUTS_ID",
        "CNTR_CODE",
        "nuts_release",
        "reference_population_density",
        "density_band",
        "H_core",
        "H_group",
        "previous_radiance",
        radiance_col,
        "delta_log10_radiance",
        "individual_pct_change",
    ]
    return d[columns].copy()


def summarize_pairs(
    pairs: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    cfg = config["analysis"]["paired_density_h_trends"]

    rows = []
    group_cols = [
        "transition_end_year",
        "density_band",
        "H_group",
    ]

    for keys, g in pairs.groupby(group_cols, observed=True):
        year, density_band, h_group = keys
        x = g["delta_log10_radiance"].to_numpy(dtype=float)
        pct = g["individual_pct_change"].to_numpy(dtype=float)

        n = len(g)
        mean_log = float(np.mean(x))
        sd_log = float(np.std(x, ddof=1)) if n > 1 else np.nan
        se_log = sd_log / np.sqrt(n) if n > 1 else np.nan

        if n > 1:
            log_lo = mean_log - 1.96 * se_log
            log_hi = mean_log + 1.96 * se_log
        else:
            log_lo = np.nan
            log_hi = np.nan

        rows.append(
            {
                "transition_end_year": int(year),
                "transition_label": f"{int(year)-1}-{int(year)}",
                "density_band": str(density_band),
                "H_group": str(h_group),
                "n_regions": int(n),
                "mean_delta_log10_radiance": mean_log,
                "sd_delta_log10_radiance": sd_log,
                "se_delta_log10_radiance": se_log,
                "geometric_mean_pct_change": (
                    100.0 * (10.0**mean_log - 1.0)
                ),
                "geometric_mean_pct_ci_low": (
                    100.0 * (10.0**log_lo - 1.0)
                    if np.isfinite(log_lo)
                    else np.nan
                ),
                "geometric_mean_pct_ci_high": (
                    100.0 * (10.0**log_hi - 1.0)
                    if np.isfinite(log_hi)
                    else np.nan
                ),
                "mean_individual_pct_change": float(np.mean(pct)),
                "median_individual_pct_change": float(np.median(pct)),
                "individual_pct_q10": float(np.quantile(pct, 0.10)),
                "individual_pct_q90": float(np.quantile(pct, 0.90)),
            }
        )

    out = pd.DataFrame(rows)

    min_n = int(cfg["diagnostic_min_regions"])
    threshold = float(cfg["diagnostic_abs_pct_threshold"])
    out["diagnostic_small_n"] = out["n_regions"] < min_n
    out["diagnostic_large_change"] = (
        out["geometric_mean_pct_change"].abs() >= threshold
    )
    out["diagnostic_flag"] = (
        out["diagnostic_small_n"] | out["diagnostic_large_change"]
    )

    return out


def plot_summary(
    summary: pd.DataFrame,
    config: dict,
    output_png: Path,
    output_pdf: Path,
) -> None:
    cfg = config["analysis"]["paired_density_h_trends"]
    density_labels = list(cfg["density_band_labels"])

    h_groups = list(
        dict.fromkeys(
            summary["H_group"].dropna().astype(str).tolist()
        )
    )

    # Deliberately use both colour and shape/line type so the H categories
    # remain distinguishable in print and for colour-vision deficiencies.
    colours = [
        "#1b4f72",
        "#2e86ab",
        "#48a868",
        "#d98e04",
        "#b03a2e",
        "#6c3483",
        "#566573",
    ]
    markers = ["o", "s", "D", "^", "X", "v", "P"]
    linestyles = [
        "-",
        "--",
        ":",
        "-.",
        (0, (5, 1, 1, 1)),
        (0, (2, 1)),
        (0, (7, 2)),
    ]

    style = {
        h: (
            colours[i % len(colours)],
            markers[i % len(markers)],
            linestyles[i % len(linestyles)],
        )
        for i, h in enumerate(h_groups)
    }

    fig, axes = plt.subplots(
        3,
        2,
        figsize=(15.5, 13.0),
        sharex=True,
    )
    axes = axes.ravel()
    for ax in axes:
        ax.set_visible(False)

    for i, density_band in enumerate(density_labels):
        ax = axes[i]
        ax.set_visible(True)
        panel = summary[
            summary["density_band"] == density_band
        ].copy()

        ax.axhline(0.0, color="0.45", linewidth=1.0)

        ymin = np.inf
        ymax = -np.inf

        for h_group in h_groups:
            g = panel[
                panel["H_group"] == h_group
            ].sort_values("transition_end_year")

            if g.empty:
                continue

            colour, marker, linestyle = style[h_group]

            x = g["transition_end_year"].to_numpy(dtype=float)
            y = g["geometric_mean_pct_change"].to_numpy(dtype=float)
            lo = g["geometric_mean_pct_ci_low"].to_numpy(dtype=float)
            hi = g["geometric_mean_pct_ci_high"].to_numpy(dtype=float)

            ax.plot(
                x,
                y,
                color=colour,
                marker=marker,
                linestyle=linestyle,
                linewidth=2.0,
                markersize=5.5,
            )
            ax.fill_between(
                x,
                lo,
                hi,
                color=colour,
                alpha=0.08,
                linewidth=0,
            )

            finite = np.isfinite(np.r_[lo, hi])
            if finite.any():
                vals = np.r_[lo, hi][finite]
                ymin = min(ymin, float(vals.min()))
                ymax = max(ymax, float(vals.max()))

            flagged = g[g["diagnostic_flag"]]
            for _, row in flagged.iterrows():
                ax.annotate(
                    f"n={int(row['n_regions'])}",
                    (
                        row["transition_end_year"],
                        row["geometric_mean_pct_change"],
                    ),
                    xytext=(4, 5),
                    textcoords="offset points",
                    fontsize=7.5,
                    color=colour,
                )

        title_label = density_band.replace(">=", "≥").replace("-", "–")
        ax.set_title(
            f"{title_label} people km$^{{-2}}$",
            fontsize=12,
        )
        ax.grid(True, alpha=0.25)
        ax.set_xlim(
            int(cfg["first_year"]) + 1,
            int(cfg["last_year"]),
        )
        ax.set_xticks(
            range(
                int(cfg["first_year"]) + 1,
                int(cfg["last_year"]) + 1,
            )
        )

        if np.isfinite(ymin) and np.isfinite(ymax):
            pad = max(
                1.0,
                0.10 * (ymax - ymin if ymax > ymin else 10.0),
            )
            ax.set_ylim(ymin - pad, ymax + pad)

        if i % 2 == 0:
            ax.set_ylabel(
                "Paired annual change in\ncorrected radiance (%)",
                fontsize=10.5,
            )

    axes[5].axis("off")

    handles = []
    for h_group in h_groups:
        colour, marker, linestyle = style[h_group]
        handles.append(
            Line2D(
                [0],
                [0],
                color=colour,
                marker=marker,
                linestyle=linestyle,
                linewidth=2.2,
                markersize=6,
                label=h_group,
            )
        )

    fig.suptitle(
        "Paired annual change in corrected VIIRS radiance",
        fontsize=18,
        y=0.975,
    )
    fig.text(
        0.5,
        0.947,
        (
            "Fixed NUTS2016 panel; fixed population-density bands; "
            "fixed within-band H quantiles"
        ),
        ha="center",
        fontsize=11.5,
    )
    fig.text(
        0.5,
        0.928,
        (
            "Each point is the geometric-mean change of the same NUTS3 "
            "regions observed in consecutive years; ribbons are approximate 95% CIs"
        ),
        ha="center",
        fontsize=10,
    )
    fig.legend(
        handles=handles,
        title="Habitability category",
        loc="lower center",
        ncol=min(len(handles), 5),
        frameon=True,
        bbox_to_anchor=(0.5, 0.025),
    )
    fig.supxlabel("Transition end year", y=0.075, fontsize=12)

    fig.subplots_adjust(
        top=0.89,
        bottom=0.14,
        left=0.08,
        right=0.98,
        hspace=0.30,
        wspace=0.10,
    )

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_png,
        dpi=int(cfg["figure_dpi"]),
        bbox_inches="tight",
    )
    fig.savefig(
        output_pdf,
        bbox_inches="tight",
    )
    plt.close(fig)


def write_summary(
    path: Path,
    region: pd.DataFrame,
    cutpoints: pd.DataFrame,
    pairs: pd.DataFrame,
    summary: pd.DataFrame,
    config: dict,
) -> None:
    cfg = config["analysis"]["paired_density_h_trends"]

    with path.open("w", encoding="utf-8") as h:
        h.write("PAIRED FIXED-MEMBERSHIP DENSITY x HABITABILITY TRENDS\n")
        h.write("=" * 76 + "\n\n")

        h.write(
            f"Interval: {cfg['first_year']}-{cfg['last_year']}; "
            f"NUTS release: {cfg['nuts_release']}.\n"
        )
        h.write(
            f"Radiance column: {cfg['radiance_column']}.\n"
        )
        h.write(
            "WARNING: the current production corrected-radiance column is "
            "still provisional until the dataset-wide 2017 calibration "
            "correction is fixed. Rerun this stage afterwards.\n\n"
        )

        h.write(
            f"Regions assigned: {region['NUTS_ID'].nunique()}.\n"
        )
        h.write(
            f"Paired consecutive-year observations: {len(pairs)}.\n"
        )
        h.write(
            "Population-density membership is fixed from each region's "
            "median density across the complete analysis interval.\n"
        )
        h.write(
            f"H grouping: {cfg['h_grouping']}; "
            f"{cfg['h_quantile_groups']} groups.\n\n"
        )

        h.write("FIXED REGION COUNTS BY DENSITY x H GROUP\n")
        counts = (
            region.groupby(
                ["density_band", "H_group"],
                observed=True,
            )
            .size()
            .rename("n_regions")
            .reset_index()
        )
        h.write(counts.to_string(index=False))
        h.write("\n\nH GROUP RANGES\n")
        h.write(
            cutpoints.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        h.write("\n\nANNUAL PAIRED GROUP SUMMARY\n")
        cols = [
            "transition_label",
            "density_band",
            "H_group",
            "n_regions",
            "geometric_mean_pct_change",
            "geometric_mean_pct_ci_low",
            "geometric_mean_pct_ci_high",
            "mean_individual_pct_change",
            "median_individual_pct_change",
            "diagnostic_flag",
        ]
        h.write(
            summary[cols].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

        flagged = summary[summary["diagnostic_flag"]].copy()
        h.write("\n\nFLAGGED CELLS\n")
        if flagged.empty:
            h.write("None.\n")
        else:
            h.write(
                flagged[cols].to_string(
                    index=False,
                    float_format=lambda x: f"{x:.6g}",
                )
            )
            h.write("\n")

        h.write("\nINTERPRETATION\n")
        h.write(
            "- These plots avoid the earlier changing-membership problem: "
            "annual changes are calculated within each NUTS3 first.\n"
        )
        h.write(
            "- Fixed density bands prevent a region crossing a threshold from "
            "changing group between adjacent years.\n"
        )
        h.write(
            "- H categories are fixed for the entire interval because H is a "
            "static physical-geography metric.\n"
        )
        h.write(
            "- The primary percentage estimate is the back-transformed mean "
            "log10 change, matching the scale used in the main temporal models.\n"
        )
        h.write(
            "- This stage is descriptive; it does not replace the mixed model "
            "that estimates the long-term density-dependent temporal trend.\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["paired_density_h_trends"]

    data, region = load_and_prepare(config)
    region, cutpoints = assign_h_groups(region, config)
    pairs = make_paired_changes(data, region, config)
    summary = summarize_pairs(pairs, config)

    processed = configured_path(config, "processed")
    figdir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "paired_density_h_trends"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    region_path = processed / "paired_density_h_region_assignments.csv"
    cutpoints_path = processed / "paired_density_h_cutpoints.csv"
    pairs_path = processed / "paired_density_h_annual_changes.csv"
    group_path = processed / "paired_density_h_group_summary.csv"
    diagnostics_path = processed / "paired_density_h_diagnostics.csv"
    summary_path = processed / "paired_density_h_trends_summary.txt"
    png_path = figdir / "paired_density_h_trends_2013_2020.png"
    pdf_path = figdir / "paired_density_h_trends_2013_2020.pdf"

    region.to_csv(region_path, index=False)
    cutpoints.to_csv(cutpoints_path, index=False)
    pairs.to_csv(pairs_path, index=False)
    summary.to_csv(group_path, index=False)
    summary[summary["diagnostic_flag"]].to_csv(
        diagnostics_path,
        index=False,
    )

    plot_summary(
        summary,
        config,
        png_path,
        pdf_path,
    )
    write_summary(
        summary_path,
        region,
        cutpoints,
        pairs,
        summary,
        config,
    )

    print("\nPAIRED FIXED-MEMBERSHIP DENSITY x H TRENDS")
    print(
        f"Interval: {cfg['first_year']}-{cfg['last_year']} "
        f"(NUTS {cfg['nuts_release']})"
    )
    print(
        f"Regions: {region['NUTS_ID'].nunique():,}; "
        f"paired observations: {len(pairs):,}"
    )
    print(
        "\nWARNING: radiance_mean_corrected remains provisional until the "
        "global 2017 VIIRS calibration correction is fixed. Rerun this stage "
        "after that production update."
    )

    print("\nFIXED REGION COUNTS")
    print(
        region.groupby(
            ["density_band", "H_group"],
            observed=True,
        )
        .size()
        .rename("n_regions")
        .reset_index()
        .to_string(index=False)
    )

    flagged = summary[summary["diagnostic_flag"]]
    print(f"\nFlagged group-year cells: {len(flagged)}")
    if not flagged.empty:
        print(
            flagged[
                [
                    "transition_label",
                    "density_band",
                    "H_group",
                    "n_regions",
                    "geometric_mean_pct_change",
                ]
            ].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

    print("\nWROTE")
    for path in [
        region_path,
        cutpoints_path,
        pairs_path,
        group_path,
        diagnostics_path,
        summary_path,
        png_path,
        pdf_path,
    ]:
        print(path)


if __name__ == "__main__":
    main()
