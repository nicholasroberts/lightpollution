#!/usr/bin/env python3
"""Compare current country fits against the archived 2016-2019 benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_config import PROJECT_ROOT, configured_path, load_config


LEGACY_FITS = (
    PROJECT_ROOT
    / "archive"
    / "legacy_2016_2019"
    / "outputs"
    / "6.rad_density_relationships.csv"
)


def load_legacy() -> pd.DataFrame:
    legacy = pd.read_csv(LEGACY_FITS)

    # Drop an R-written row-number column if present.
    unnamed = [c for c in legacy.columns if c.startswith("Unnamed")]
    if unnamed:
        legacy = legacy.drop(columns=unnamed)

    slope_cols = ["slope2016", "slope2017", "slope2018", "slope2019"]
    intercept_cols = [
        "intercept2016",
        "intercept2017",
        "intercept2018",
        "intercept2019",
    ]

    legacy["legacy_slope_mean"] = legacy[slope_cols].mean(axis=1, skipna=True)
    legacy["legacy_slope_min"] = legacy[slope_cols].min(axis=1, skipna=True)
    legacy["legacy_slope_max"] = legacy[slope_cols].max(axis=1, skipna=True)

    legacy["legacy_intercept_mean"] = legacy[intercept_cols].mean(axis=1, skipna=True)
    legacy["legacy_intercept_min"] = legacy[intercept_cols].min(axis=1, skipna=True)
    legacy["legacy_intercept_max"] = legacy[intercept_cols].max(axis=1, skipna=True)

    return legacy.rename(columns={"country": "CNTR_CODE"})


def make_figure(frame: pd.DataFrame, output: Path) -> None:
    frame = frame.sort_values("CNTR_CODE").reset_index(drop=True)
    x = np.arange(len(frame))

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13, 5),
        constrained_layout=True,
    )

    # Slope comparison.
    axes[0].errorbar(
        x,
        frame["legacy_slope_mean"],
        yerr=[
            frame["legacy_slope_mean"] - frame["legacy_slope_min"],
            frame["legacy_slope_max"] - frame["legacy_slope_mean"],
        ],
        fmt="o",
        capsize=4,
        label="Legacy 2016-2019 mean ± range",
    )
    axes[0].scatter(
        x,
        frame["slope_2024"],
        marker="D",
        s=45,
        label="New 2024",
    )
    axes[0].set_xticks(x, frame["CNTR_CODE"])
    axes[0].set_ylabel("Log-log slope")
    axes[0].set_title("Country slope comparison")
    axes[0].grid(alpha=0.2)
    axes[0].legend()

    # Intercept comparison.
    axes[1].errorbar(
        x,
        frame["legacy_intercept_mean"],
        yerr=[
            frame["legacy_intercept_mean"] - frame["legacy_intercept_min"],
            frame["legacy_intercept_max"] - frame["legacy_intercept_mean"],
        ],
        fmt="o",
        capsize=4,
        label="Legacy 2016-2019 mean ± range",
    )
    axes[1].scatter(
        x,
        frame["intercept_2024"],
        marker="D",
        s=45,
        label="New 2024",
    )
    axes[1].set_xticks(x, frame["CNTR_CODE"])
    axes[1].set_ylabel("Log-log intercept")
    axes[1].set_title("Country intercept comparison")
    axes[1].grid(alpha=0.2)
    axes[1].legend()

    fig.suptitle(
        "Legacy 2016-2019 benchmark vs rebuilt 2024 pipeline",
        fontsize=14,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--year", type=int, default=2024)
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")
    figures = configured_path(config, "figures")

    current_path = processed / "country_loglog_fits.csv"
    current = pd.read_csv(current_path)
    current = current[current["year"] == args.year].copy()

    current = current.rename(
        columns={
            "slope": f"slope_{args.year}",
            "intercept": f"intercept_{args.year}",
            "r2": f"r2_{args.year}",
            "n_fit": f"n_fit_{args.year}",
        }
    )

    legacy = load_legacy()

    keep_current = [
        "CNTR_CODE",
        f"slope_{args.year}",
        f"intercept_{args.year}",
        f"r2_{args.year}",
        f"n_fit_{args.year}",
    ]
    comparison = legacy.merge(
        current[keep_current],
        on="CNTR_CODE",
        how="inner",
        validate="one_to_one",
    )

    comparison["slope_difference_from_legacy_mean"] = (
        comparison[f"slope_{args.year}"] - comparison["legacy_slope_mean"]
    )
    comparison["intercept_difference_from_legacy_mean"] = (
        comparison[f"intercept_{args.year}"] - comparison["legacy_intercept_mean"]
    )
    comparison["slope_2024_within_legacy_range"] = (
        comparison[f"slope_{args.year}"].between(
            comparison["legacy_slope_min"],
            comparison["legacy_slope_max"],
            inclusive="both",
        )
    )
    comparison["intercept_2024_within_legacy_range"] = (
        comparison[f"intercept_{args.year}"].between(
            comparison["legacy_intercept_min"],
            comparison["legacy_intercept_max"],
            inclusive="both",
        )
    )

    output_csv = processed / f"legacy_vs_{args.year}_country_fits.csv"
    comparison.to_csv(output_csv, index=False)

    benchmark_codes = ["DE", "IT", "NL", "FR", "ES"]
    benchmark = comparison[
        comparison["CNTR_CODE"].isin(benchmark_codes)
    ].copy()

    output_pdf = (
        figures
        / "population_relationship"
        / f"legacy_vs_{args.year}_country_fits.pdf"
    )
    make_figure(benchmark, output_pdf)

    display_cols = [
        "CNTR_CODE",
        "legacy_slope_mean",
        "legacy_slope_min",
        "legacy_slope_max",
        f"slope_{args.year}",
        "slope_difference_from_legacy_mean",
        "slope_2024_within_legacy_range",
        "legacy_intercept_mean",
        f"intercept_{args.year}",
        "intercept_difference_from_legacy_mean",
        "intercept_2024_within_legacy_range",
    ]

    print(
        comparison[display_cols]
        .sort_values("CNTR_CODE")
        .to_string(
            index=False,
            float_format=lambda value: f"{value:.4f}",
        )
    )
    print()
    print(f"Wrote: {output_csv}")
    print(f"Wrote: {output_pdf}")


if __name__ == "__main__":
    main()
