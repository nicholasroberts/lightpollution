#!/usr/bin/env python3
"""Test whether the 2017 VIIRS discontinuity is a simple additive radiance offset.

Uses matched NUTS3 regions on common NUTS 2016 geometry for 2015-2018.
A zero-point shift predicts that the excess 2016->2017 change should be
approximately constant in LINEAR radiance units across the brightness range.

No production radiance values are modified.
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_config import configured_path, load_config


YEARS = [2015, 2016, 2017, 2018]


def candidate_offset(frame: pd.DataFrame) -> dict[str, float]:
    m15 = frame["d15_16"].median()
    m16 = frame["d16_17"].median()
    m17 = frame["d17_18"].median()
    expected = np.nanmean([m15, m17])
    return {
        "median_d15_16": float(m15),
        "median_d16_17": float(m16),
        "median_d17_18": float(m17),
        "expected_background_change": float(expected),
        "candidate_offset": float(m16 - expected),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--quantiles", type=int, default=10)
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")
    figures = configured_path(config, "figures") / "annual_change"
    figures.mkdir(parents=True, exist_ok=True)

    radiance = pd.read_parquet(
        processed / config["outputs"]["radiance_by_nuts3"]
    )

    x = radiance[radiance["year"].isin(YEARS)][
        ["year", "NUTS_ID", "CNTR_CODE", "radiance_mean"]
    ].copy()

    wide = x.pivot_table(
        index=["NUTS_ID", "CNTR_CODE"],
        columns="year",
        values="radiance_mean",
        aggfunc="first",
    ).reset_index()

    for year in YEARS:
        if year not in wide.columns:
            wide[year] = np.nan

    wide = wide.dropna(subset=YEARS).copy()
    wide["d15_16"] = wide[2016] - wide[2015]
    wide["d16_17"] = wide[2017] - wide[2016]
    wide["d17_18"] = wide[2018] - wide[2017]

    overall = candidate_offset(wide)

    print("Matched NUTS3 regions in all four years:", len(wide))
    print("\nOverall transition medians (nW cm-2 sr-1):")
    for key, value in overall.items():
        print(f"  {key}: {value:.6f}")

    # Test constancy across the brightness range.
    q = wide.copy()
    q["radiance_quantile"] = pd.qcut(
        q[2016], q=args.quantiles, duplicates="drop"
    )

    q_rows = []
    for i, (_, g) in enumerate(
        q.groupby("radiance_quantile", observed=True), start=1
    ):
        q_rows.append(
            {
                "quantile": i,
                "radiance_2016_min": float(g[2016].min()),
                "radiance_2016_median": float(g[2016].median()),
                "radiance_2016_max": float(g[2016].max()),
                "n": len(g),
                **candidate_offset(g),
            }
        )
    qdf = pd.DataFrame(q_rows)

    # Country-level robustness.
    country_rows = []
    for country, g in wide.groupby("CNTR_CODE"):
        if len(g) < 10:
            continue
        country_rows.append(
            {"CNTR_CODE": country, "n": len(g), **candidate_offset(g)}
        )
    cdf = pd.DataFrame(country_rows).sort_values("candidate_offset")

    q_path = processed / "calibration_offset_by_radiance_quantile.csv"
    c_path = processed / "calibration_offset_by_country.csv"
    matched_path = processed / "calibration_offset_matched_regions.csv"
    qdf.to_csv(q_path, index=False)
    cdf.to_csv(c_path, index=False)
    wide.to_csv(matched_path, index=False)

    print("\nCandidate offset by 2016 radiance quantile:")
    print(
        qdf[
            [
                "quantile",
                "radiance_2016_median",
                "n",
                "median_d15_16",
                "median_d16_17",
                "median_d17_18",
                "candidate_offset",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.6f}")
    )

    print("\nCandidate offset by country (n >= 10):")
    print(
        cdf[
            [
                "CNTR_CODE",
                "n",
                "median_d15_16",
                "median_d16_17",
                "median_d17_18",
                "candidate_offset",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.6f}")
    )

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(
        qdf["radiance_2016_median"],
        qdf["candidate_offset"],
        marker="o",
    )
    ax.axhline(
        overall["candidate_offset"],
        linestyle="--",
        linewidth=1.2,
        label=f"overall = {overall['candidate_offset']:.3f}",
    )
    ax.set_xscale("log")
    ax.set_xlabel("Median 2016 NUTS3 radiance (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_ylabel(
        "Estimated excess 2016->2017 radiance "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title("Test of simple additive 2017 VIIRS radiance offset")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()

    fig_path = figures / "calibration_offset_by_radiance_quantile.pdf"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print()
    print(f"Wrote: {q_path}")
    print(f"Wrote: {c_path}")
    print(f"Wrote: {matched_path}")
    print(f"Wrote: {fig_path}")


if __name__ == "__main__":
    main()
