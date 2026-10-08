#!/usr/bin/env python3
"""Apply the production VIIRS additive zero-point correction.

The raw zonal radiance extraction is preserved unchanged. Country-specific
additive offsets are estimated from the darkest configured fraction of matched
NUTS3 regions present in 2015-2018, using:

    offset = median(2017-2016)
             - mean[median(2016-2015), median(2018-2017)]

The offset is subtracted from 2017 onward for configured calibration countries.
No clipping or pseudocount is applied; corrected values <= 0 are retained and
are excluded only when a logarithmic fit is performed.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from project_config import configured_path, load_config


CAL_YEARS = [2015, 2016, 2017, 2018]


def estimate_country_offset(
    frame: pd.DataFrame,
    country: str,
    dark_fraction: float,
) -> dict:
    g = frame[
        (frame["CNTR_CODE"] == country)
        & frame["year"].isin(CAL_YEARS)
    ].copy()

    wide = (
        g.pivot_table(
            index="NUTS_ID",
            columns="year",
            values="radiance_mean",
            aggfunc="first",
        )
        .dropna(subset=CAL_YEARS)
        .copy()
    )

    if len(wide) < 10:
        raise ValueError(
            f"{country}: only {len(wide)} matched regions in 2015-2018; "
            "at least 10 are required for production calibration."
        )

    threshold = wide[2016].quantile(dark_fraction)
    dark = wide[wide[2016] <= threshold].copy()

    d15 = (dark[2016] - dark[2015]).median()
    d16 = (dark[2017] - dark[2016]).median()
    d17 = (dark[2018] - dark[2017]).median()
    expected = np.nanmean([d15, d17])
    offset = d16 - expected

    return {
        "CNTR_CODE": country,
        "n_matched": len(wide),
        "n_dark": len(dark),
        "dark_threshold_2016": float(threshold),
        "median_d15_16": float(d15),
        "median_d16_17": float(d16),
        "median_d17_18": float(d17),
        "calibration_offset": float(offset),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")

    cal = config["calibration"]
    if not cal["enabled"]:
        raise ValueError("Production calibration is disabled in configuration.")

    dark_fraction = float(cal["dark_fraction"])
    start_year = int(cal["apply_from_year"])
    countries = list(cal["countries"])

    raw_path = processed / config["outputs"]["radiance_by_nuts3"]
    radiance = pd.read_parquet(raw_path).copy()

    required = {"year", "NUTS_ID", "CNTR_CODE", "radiance_mean"}
    missing = required - set(radiance.columns)
    if missing:
        raise ValueError(
            f"Raw radiance table is missing required columns: {sorted(missing)}"
        )

    offsets = pd.DataFrame(
        [
            estimate_country_offset(radiance, country, dark_fraction)
            for country in countries
        ]
    ).sort_values("CNTR_CODE")

    offset_map = dict(
        zip(offsets["CNTR_CODE"], offsets["calibration_offset"])
    )

    out = radiance.copy()
    out["radiance_mean_raw"] = out["radiance_mean"]
    out["calibration_offset"] = out["CNTR_CODE"].map(offset_map).fillna(0.0)

    eligible_country = out["CNTR_CODE"].isin(countries)
    post_change = out["year"] >= start_year
    out["calibration_correction_applied"] = eligible_country & post_change

    out["radiance_mean_corrected"] = out["radiance_mean_raw"]
    mask = out["calibration_correction_applied"]
    out.loc[mask, "radiance_mean_corrected"] = (
        out.loc[mask, "radiance_mean_raw"]
        - out.loc[mask, "calibration_offset"]
    )

    out["calibration_method"] = cal["method"]
    out["calibration_dark_fraction"] = dark_fraction
    out["calibration_apply_from_year"] = start_year

    corrected_path = (
        processed / config["outputs"]["radiance_by_nuts3_calibrated"]
    )
    offsets_path = processed / config["outputs"]["calibration_offsets"]

    out.to_parquet(corrected_path, index=False)
    out.to_csv(corrected_path.with_suffix(".csv"), index=False)
    offsets.to_csv(offsets_path, index=False)

    print(
        f"Production calibration: darkest {dark_fraction:.0%}, "
        f"applied from {start_year}"
    )
    print()
    print(offsets.to_string(index=False, float_format=lambda x: f"{x:.6f}"))

    affected = out[out["calibration_correction_applied"]]
    nonpositive = (
        affected["radiance_mean_corrected"] <= 0
    ).sum()

    print()
    print(f"Corrected rows: {len(affected):,}")
    print(
        "Corrected rows <= 0 (retained, excluded only from log fits): "
        f"{nonpositive:,}"
    )
    print(f"Wrote: {corrected_path}")
    print(f"Wrote: {corrected_path.with_suffix('.csv')}")
    print(f"Wrote: {offsets_path}")


if __name__ == "__main__":
    main()
