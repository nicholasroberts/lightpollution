#!/usr/bin/env python3
"""Sensitivity analysis for country-specific additive VIIRS calibration offsets.

Tests dark-region fractions 10, 15, 20, 25 and 30%.

For each threshold:
- estimate country-specific additive offsets from dark matched 2015-2018 regions;
- apply the offset to 2017+ radiance for diagnostic fits;
- quantify the 2016->2017 slope/intercept discontinuity;
- compare 2017-2019 temporal changes relative to 2016 against the archived
  legacy analysis, which had an upstream zero-point correction.

Production radiance is not modified.
"""

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

FRACTIONS = [0.10, 0.15, 0.20, 0.25, 0.30]


def fit_line(g: pd.DataFrame, radiance_col: str):
    x = g["population_density"].to_numpy(dtype=float)
    y = g[radiance_col].to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = np.log10(x[ok])
    y = np.log10(y[ok])

    if len(x) < 2:
        return len(x), np.nan, np.nan, np.nan

    slope, intercept = np.polyfit(x, y, 1)
    pred = intercept + slope * x
    denom = ((y - y.mean()) ** 2).sum()
    r2 = np.nan if denom == 0 else 1 - ((y - pred) ** 2).sum() / denom
    return len(x), float(slope), float(intercept), float(r2)


def estimate_offset(g: pd.DataFrame, fraction: float):
    wide = (
        g[g["year"].isin([2015, 2016, 2017, 2018])]
        .pivot_table(
            index="NUTS_ID",
            columns="year",
            values="radiance_mean",
            aggfunc="first",
        )
        .dropna(subset=[2015, 2016, 2017, 2018])
        .copy()
    )

    if len(wide) < 10:
        return np.nan, len(wide), 0, np.nan

    threshold = wide[2016].quantile(fraction)
    dark = wide[wide[2016] <= threshold].copy()

    d15 = (dark[2016] - dark[2015]).median()
    d16 = (dark[2017] - dark[2016]).median()
    d17 = (dark[2018] - dark[2017]).median()
    offset = d16 - np.nanmean([d15, d17])

    return float(offset), len(wide), len(dark), float(threshold)


def load_legacy_long(countries):
    legacy = pd.read_csv(LEGACY_FITS)
    unnamed = [c for c in legacy.columns if c.startswith("Unnamed")]
    if unnamed:
        legacy = legacy.drop(columns=unnamed)

    legacy = legacy[legacy["country"].isin(countries)].copy()

    rows = []
    for row in legacy.itertuples(index=False):
        for year in [2016, 2017, 2018, 2019]:
            rows.append(
                {
                    "CNTR_CODE": row.country,
                    "year": year,
                    "legacy_slope": getattr(row, f"slope{year}"),
                    "legacy_intercept": getattr(row, f"intercept{year}"),
                }
            )
    out = pd.DataFrame(rows).sort_values(["CNTR_CODE", "year"])
    out["legacy_delta_slope_from_2016"] = (
        out["legacy_slope"]
        - out.groupby("CNTR_CODE")["legacy_slope"].transform("first")
    )
    out["legacy_delta_intercept_from_2016"] = (
        out["legacy_intercept"]
        - out.groupby("CNTR_CODE")["legacy_intercept"].transform("first")
    )
    return out


def fit_threshold(merged, countries, fraction):
    offsets = []
    for country in countries:
        g = merged[merged["CNTR_CODE"] == country]
        offset, n_matched, n_dark, threshold = estimate_offset(g, fraction)
        offsets.append(
            {
                "dark_fraction": fraction,
                "CNTR_CODE": country,
                "offset": offset,
                "n_matched": n_matched,
                "n_dark": n_dark,
                "dark_threshold_2016": threshold,
            }
        )

    odf = pd.DataFrame(offsets)
    offset_map = dict(zip(odf["CNTR_CODE"], odf["offset"]))

    work = merged[merged["CNTR_CODE"].isin(countries)].copy()
    work["radiance_corrected_test"] = work["radiance_mean"]

    mapped = work["CNTR_CODE"].map(offset_map)
    apply = (work["year"] >= 2017) & mapped.notna()
    work.loc[apply, "radiance_corrected_test"] = (
        work.loc[apply, "radiance_mean"] - mapped[apply]
    )

    rows = []
    for country in countries:
        for year in range(2013, 2025):
            g = work[
                (work["CNTR_CODE"] == country)
                & (work["year"] == year)
            ]
            n0, s0, i0, r0 = fit_line(g, "radiance_mean")
            n1, s1, i1, r1 = fit_line(g, "radiance_corrected_test")
            rows.append(
                {
                    "dark_fraction": fraction,
                    "CNTR_CODE": country,
                    "year": year,
                    "n_uncorrected": n0,
                    "slope_uncorrected": s0,
                    "intercept_uncorrected": i0,
                    "r2_uncorrected": r0,
                    "n_corrected": n1,
                    "slope_corrected": s1,
                    "intercept_corrected": i1,
                    "r2_corrected": r1,
                }
            )

    fits = pd.DataFrame(rows).sort_values(["CNTR_CODE", "year"])

    for prefix in ["uncorrected", "corrected"]:
        fits[f"delta_slope_{prefix}"] = (
            fits.groupby("CNTR_CODE")[f"slope_{prefix}"].diff()
        )
        fits[f"delta_intercept_{prefix}"] = (
            fits.groupby("CNTR_CODE")[f"intercept_{prefix}"].diff()
        )

        baseline_slope = fits.groupby("CNTR_CODE")[f"slope_{prefix}"].transform(
            lambda s: s.loc[fits.loc[s.index, "year"].eq(2016)].iloc[0]
            if fits.loc[s.index, "year"].eq(2016).any()
            else np.nan
        )
        baseline_intercept = fits.groupby("CNTR_CODE")[
            f"intercept_{prefix}"
        ].transform(
            lambda s: s.loc[fits.loc[s.index, "year"].eq(2016)].iloc[0]
            if fits.loc[s.index, "year"].eq(2016).any()
            else np.nan
        )
        fits[f"delta_slope_from_2016_{prefix}"] = (
            fits[f"slope_{prefix}"] - baseline_slope
        )
        fits[f"delta_intercept_from_2016_{prefix}"] = (
            fits[f"intercept_{prefix}"] - baseline_intercept
        )

    return odf, fits


def rmse(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return np.sqrt(np.mean(x**2)) if len(x) else np.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument(
        "--countries",
        nargs="+",
        default=None,
        help="Defaults to analysis.comparison_countries.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    countries = (
        args.countries
        if args.countries
        else config["analysis"]["comparison_countries"]
    )

    processed = configured_path(config, "processed")
    figures = configured_path(config, "figures") / "annual_change"
    figures.mkdir(parents=True, exist_ok=True)

    merged = pd.read_parquet(
        processed / config["outputs"]["merged"]
    )
    legacy = load_legacy_long(countries)

    all_offsets = []
    all_fits = []
    summary_rows = []

    for fraction in FRACTIONS:
        offsets, fits = fit_threshold(merged, countries, fraction)
        all_offsets.append(offsets)
        all_fits.append(fits)

        jump = fits[fits["year"] == 2017].copy()

        legacy_compare = fits[
            fits["year"].isin([2017, 2018, 2019])
        ].merge(
            legacy,
            on=["CNTR_CODE", "year"],
            how="inner",
            validate="one_to_one",
        )

        legacy_compare["slope_change_error_uncorrected"] = (
            legacy_compare["delta_slope_from_2016_uncorrected"]
            - legacy_compare["legacy_delta_slope_from_2016"]
        )
        legacy_compare["slope_change_error_corrected"] = (
            legacy_compare["delta_slope_from_2016_corrected"]
            - legacy_compare["legacy_delta_slope_from_2016"]
        )
        legacy_compare["intercept_change_error_uncorrected"] = (
            legacy_compare["delta_intercept_from_2016_uncorrected"]
            - legacy_compare["legacy_delta_intercept_from_2016"]
        )
        legacy_compare["intercept_change_error_corrected"] = (
            legacy_compare["delta_intercept_from_2016_corrected"]
            - legacy_compare["legacy_delta_intercept_from_2016"]
        )

        summary_rows.append(
            {
                "dark_fraction": fraction,
                "mean_abs_2017_slope_jump_uncorrected": jump[
                    "delta_slope_uncorrected"
                ].abs().mean(),
                "mean_abs_2017_slope_jump_corrected": jump[
                    "delta_slope_corrected"
                ].abs().mean(),
                "mean_abs_2017_intercept_jump_uncorrected": jump[
                    "delta_intercept_uncorrected"
                ].abs().mean(),
                "mean_abs_2017_intercept_jump_corrected": jump[
                    "delta_intercept_corrected"
                ].abs().mean(),
                "legacy_slope_change_rmse_uncorrected": rmse(
                    legacy_compare["slope_change_error_uncorrected"]
                ),
                "legacy_slope_change_rmse_corrected": rmse(
                    legacy_compare["slope_change_error_corrected"]
                ),
                "legacy_intercept_change_rmse_uncorrected": rmse(
                    legacy_compare["intercept_change_error_uncorrected"]
                ),
                "legacy_intercept_change_rmse_corrected": rmse(
                    legacy_compare["intercept_change_error_corrected"]
                ),
            }
        )

    offsets = pd.concat(all_offsets, ignore_index=True)
    fits = pd.concat(all_fits, ignore_index=True)
    summary = pd.DataFrame(summary_rows)

    # A simple combined rank: lower is better across the four corrected metrics.
    rank_cols = [
        "mean_abs_2017_slope_jump_corrected",
        "mean_abs_2017_intercept_jump_corrected",
        "legacy_slope_change_rmse_corrected",
        "legacy_intercept_change_rmse_corrected",
    ]
    for col in rank_cols:
        summary[f"rank_{col}"] = summary[col].rank(method="average")
    summary["combined_rank"] = summary[
        [f"rank_{c}" for c in rank_cols]
    ].mean(axis=1)

    summary = summary.sort_values(
        ["combined_rank", "dark_fraction"]
    ).reset_index(drop=True)

    offsets_path = processed / "calibration_dark_fraction_offsets.csv"
    fits_path = processed / "calibration_dark_fraction_fits.csv"
    summary_path = processed / "calibration_dark_fraction_sensitivity.csv"

    offsets.to_csv(offsets_path, index=False)
    fits.to_csv(fits_path, index=False)
    summary.to_csv(summary_path, index=False)

    print("DARK-FRACTION SENSITIVITY")
    print("=" * 100)
    print(
        summary[
            [
                "dark_fraction",
                "mean_abs_2017_slope_jump_corrected",
                "mean_abs_2017_intercept_jump_corrected",
                "legacy_slope_change_rmse_corrected",
                "legacy_intercept_change_rmse_corrected",
                "combined_rank",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.6f}")
    )

    print("\nOFFSETS BY FRACTION AND COUNTRY")
    print("=" * 100)
    print(
        offsets.pivot(
            index="CNTR_CODE",
            columns="dark_fraction",
            values="offset",
        ).to_string(float_format=lambda x: f"{x:.6f}")
    )

    best = summary.iloc[0]
    print(
        f"\nBest combined threshold: {best['dark_fraction']:.0%} "
        f"(combined rank {best['combined_rank']:.3f})"
    )

    fig, ax = plt.subplots(figsize=(9, 6))
    for country in countries:
        g = offsets[offsets["CNTR_CODE"] == country].sort_values("dark_fraction")
        ax.plot(
            g["dark_fraction"] * 100,
            g["offset"],
            marker="o",
            label=country,
        )
    ax.set_xlabel("Darkest fraction used to estimate offset (%)")
    ax.set_ylabel("Estimated additive offset (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Sensitivity of country calibration offset to dark-region threshold")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    offset_fig = figures / "calibration_dark_fraction_offsets.pdf"
    fig.savefig(offset_fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6))
    x = summary["dark_fraction"] * 100
    ax.plot(
        x,
        summary["mean_abs_2017_slope_jump_corrected"],
        marker="o",
        label="2017 slope jump",
    )
    ax.plot(
        x,
        summary["legacy_slope_change_rmse_corrected"],
        marker="o",
        label="Legacy slope-change RMSE",
    )
    ax.set_xlabel("Darkest fraction used to estimate offset (%)")
    ax.set_ylabel("Slope error metric")
    ax.set_title("Calibration threshold sensitivity: slope metrics")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    slope_fig = figures / "calibration_dark_fraction_slope_metrics.pdf"
    fig.savefig(slope_fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print()
    print(f"Wrote: {offsets_path}")
    print(f"Wrote: {fits_path}")
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {offset_fig}")
    print(f"Wrote: {slope_fig}")
    print("\nDiagnostic only: production radiance was NOT modified.")


if __name__ == "__main__":
    main()
