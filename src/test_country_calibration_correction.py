#!/usr/bin/env python3
"""Test country-specific additive VIIRS offsets estimated from dark regions.

Diagnostic only: production radiance is not modified.

For each country with sufficient matched 2015-2018 NUTS3 regions:
1. select the darkest fraction of regions by 2016 radiance;
2. estimate the excess 2016->2017 change relative to neighbouring transitions;
3. subtract that additive offset from 2017+ radiance;
4. re-fit annual log10(radiance) ~ log10(population density);
5. compare 2016->2017 slope/intercept discontinuities before and after correction.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_config import configured_path, load_config


def fit_line(g: pd.DataFrame, radiance_col: str) -> tuple[int, float, float, float]:
    x = g["population_density"].to_numpy(dtype=float)
    y = g[radiance_col].to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = np.log10(x[ok])
    y = np.log10(y[ok])

    n = len(x)
    if n < 2:
        return n, np.nan, np.nan, np.nan

    slope, intercept = np.polyfit(x, y, 1)
    pred = intercept + slope * x
    denom = ((y - y.mean()) ** 2).sum()
    r2 = np.nan if denom == 0 else 1 - ((y - pred) ** 2).sum() / denom
    return n, float(slope), float(intercept), float(r2)


def estimate_offset(g: pd.DataFrame, dark_fraction: float) -> dict:
    wide = (
        g[g["year"].isin([2015, 2016, 2017, 2018])]
        .pivot_table(index="NUTS_ID", columns="year", values="radiance_mean", aggfunc="first")
        .dropna(subset=[2015, 2016, 2017, 2018])
        .copy()
    )

    if len(wide) < 10:
        return {
            "n_matched": len(wide),
            "n_dark": 0,
            "dark_threshold_2016": np.nan,
            "offset": np.nan,
            "median_d15_16": np.nan,
            "median_d16_17": np.nan,
            "median_d17_18": np.nan,
        }

    threshold = wide[2016].quantile(dark_fraction)
    d = wide[wide[2016] <= threshold].copy()

    d15 = (d[2016] - d[2015]).median()
    d16 = (d[2017] - d[2016]).median()
    d17 = (d[2018] - d[2017]).median()
    expected = np.nanmean([d15, d17])
    offset = d16 - expected

    return {
        "n_matched": len(wide),
        "n_dark": len(d),
        "dark_threshold_2016": float(threshold),
        "offset": float(offset),
        "median_d15_16": float(d15),
        "median_d16_17": float(d16),
        "median_d17_18": float(d17),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--dark-fraction", type=float, default=0.20)
    parser.add_argument(
        "--countries",
        nargs="+",
        default=None,
        help="Defaults to analysis.comparison_countries.",
    )
    args = parser.parse_args()

    if not (0 < args.dark_fraction < 1):
        raise ValueError("--dark-fraction must be between 0 and 1")

    config = load_config(args.config)
    processed = configured_path(config, "processed")
    figures = configured_path(config, "figures") / "annual_change"
    figures.mkdir(parents=True, exist_ok=True)

    merged = pd.read_parquet(
        processed / config["outputs"]["merged"]
    ).copy()

    countries = (
        args.countries
        if args.countries
        else config["analysis"]["comparison_countries"]
    )

    offsets = []
    for country in countries:
        g = merged[merged["CNTR_CODE"] == country].copy()
        est = estimate_offset(g, args.dark_fraction)
        offsets.append({"CNTR_CODE": country, **est})

    odf = pd.DataFrame(offsets)
    offset_map = dict(zip(odf["CNTR_CODE"], odf["offset"]))

    merged["radiance_corrected_test"] = merged["radiance_mean"]
    post = merged["year"] >= 2017
    mapped = merged["CNTR_CODE"].map(offset_map)
    apply = post & mapped.notna()
    merged.loc[apply, "radiance_corrected_test"] = (
        merged.loc[apply, "radiance_mean"] - mapped[apply]
    )

    rows = []
    for country in countries:
        for year in sorted(merged["year"].unique()):
            g = merged[
                (merged["CNTR_CODE"] == country)
                & (merged["year"] == year)
            ]
            n0, s0, i0, r0 = fit_line(g, "radiance_mean")
            n1, s1, i1, r1 = fit_line(g, "radiance_corrected_test")
            rows.append(
                {
                    "year": year,
                    "CNTR_CODE": country,
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
    fits["delta_slope_uncorrected"] = fits.groupby("CNTR_CODE")[
        "slope_uncorrected"
    ].diff()
    fits["delta_intercept_uncorrected"] = fits.groupby("CNTR_CODE")[
        "intercept_uncorrected"
    ].diff()
    fits["delta_slope_corrected"] = fits.groupby("CNTR_CODE")[
        "slope_corrected"
    ].diff()
    fits["delta_intercept_corrected"] = fits.groupby("CNTR_CODE")[
        "intercept_corrected"
    ].diff()

    suffix = f"dark{int(round(args.dark_fraction * 100))}"
    offset_path = processed / f"calibration_country_offsets_{suffix}.csv"
    fit_path = processed / f"calibration_corrected_fit_test_{suffix}.csv"
    odf.to_csv(offset_path, index=False)
    fits.to_csv(fit_path, index=False)

    print(
        f"Country-specific additive offsets from darkest "
        f"{args.dark_fraction:.0%} of 2016 regions:"
    )
    print(
        odf.to_string(index=False, float_format=lambda x: f"{x:.6f}")
    )

    focus = fits[fits["year"].isin([2016, 2017, 2018])].copy()
    print("\n2016-2018 fits before and after test correction:")
    print(
        focus[
            [
                "year",
                "CNTR_CODE",
                "slope_uncorrected",
                "slope_corrected",
                "delta_slope_uncorrected",
                "delta_slope_corrected",
                "intercept_uncorrected",
                "intercept_corrected",
                "delta_intercept_uncorrected",
                "delta_intercept_corrected",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.4f}")
    )

    jump = fits[fits["year"] == 2017].copy()
    print("\n2016->2017 discontinuity summary:")
    print(
        jump[
            [
                "CNTR_CODE",
                "delta_slope_uncorrected",
                "delta_slope_corrected",
                "delta_intercept_uncorrected",
                "delta_intercept_corrected",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.4f}")
    )
    print("\nMean absolute 2016->2017 jump:")
    for col in [
        "delta_slope_uncorrected",
        "delta_slope_corrected",
        "delta_intercept_uncorrected",
        "delta_intercept_corrected",
    ]:
        print(f"  {col}: {jump[col].abs().mean():.6f}")

    fig, ax = plt.subplots(figsize=(9, 6))
    for country in countries:
        g = fits[fits["CNTR_CODE"] == country]
        ax.plot(g["year"], g["slope_uncorrected"], alpha=0.35)
        ax.plot(
            g["year"],
            g["slope_corrected"],
            marker="o",
            label=f"{country} corrected",
        )
    ax.axvline(2016.5, linestyle="--", linewidth=1.2)
    ax.set_xlabel("Year")
    ax.set_ylabel("Log-log slope")
    ax.set_title(
        f"Test correction using country offsets from darkest "
        f"{args.dark_fraction:.0%}"
    )
    ax.grid(alpha=0.2)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig_path = figures / f"calibration_corrected_slope_test_{suffix}.pdf"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print()
    print(f"Wrote: {offset_path}")
    print(f"Wrote: {fit_path}")
    print(f"Wrote: {fig_path}")
    print("\nDiagnostic only: production radiance was NOT modified.")


if __name__ == "__main__":
    main()
