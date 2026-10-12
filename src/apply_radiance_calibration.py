#!/usr/bin/env python3
"""Apply the production Europe-wide 2017 VIIRS pixel calibration correction.

Production rule
---------------
The raw Annual VNL V2.2 average-masked rasters are NEVER modified.

For every raster from 2017 onward, the calibrated value supplied to exactextract
is generated on-the-fly for every raster window:

    if raw <= 0:
        corrected = raw
    else:
        corrected = max(raw - 0.235, 0)

The 0.235 nW cm-2 sr-1 offset is the rounded production value selected after
three diagnostics of the 2017 DNB discontinuity:

- a country/NUTS3 discontinuity diagnostic;
- a direct retained-lit-pixel test;
- an independent-brightness and segmented 2013-2019 time-series test.

The important distinction from the old implementation is that the correction is
applied to PIXELS before the NUTS3 mean is calculated. It is not subtracted from
the NUTS3 mean. EOG background-zero pixels therefore remain zero, while dim
positive pixels can be driven to zero.

The correction is applied to every NUTS3 country/territory in the analysis from
2017 onward. It is not country-specific.

Implementation
--------------
A custom exactextract RasterSource transforms raster windows in memory. No
corrected global GeoTIFF is written and no duplicate multi-gigabyte raster is
stored. exactextract therefore retains the same exact polygon overlap and
spherical physical-area weighting used by the raw production extraction.

The calibrated table retains:
- the original raw NUTS3 mean;
- the corrected NUTS3 mean;
- the effective NUTS3-level reduction caused by the pixel correction;
- the raw-positive area fraction;
- the fraction of area whose positive pixels were clipped to zero;
- complete correction provenance.

The raw zonal table remains unchanged.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource

from extract_nuts3_radiance import (
    discover_rasters,
    load_nuts,
    validate_raster,
)
from project_config import configured_path, load_config


class CalibrationRasterSource(RasterioRasterSource):
    """RasterSource that transforms raw VIIRS values when windows are read."""

    VALID_MODES = {
        "corrected",
        "raw_positive_indicator",
        "clipped_indicator",
    }

    def __init__(
        self,
        ds,
        *,
        offset: float,
        mode: str,
        name: str,
    ):
        if mode not in self.VALID_MODES:
            raise ValueError(f"Unknown calibration raster mode: {mode!r}")
        super().__init__(ds, 1, name=name)
        self.calibration_offset = float(offset)
        self.mode = mode

    def read_window(self, x0, y0, nx, ny):
        arr = super().read_window(x0, y0, nx, ny)

        if np.ma.isMaskedArray(arr):
            mask = np.ma.getmaskarray(arr)
            data = np.asarray(arr.data, dtype=np.float32).copy()
        else:
            mask = None
            data = np.asarray(arr, dtype=np.float32).copy()

        valid = np.isfinite(data)
        if mask is not None:
            valid &= ~mask

        positive = valid & (data > 0)

        if self.mode == "corrected":
            out = data
            out[positive] = np.maximum(
                out[positive] - self.calibration_offset,
                0.0,
            )
        elif self.mode == "raw_positive_indicator":
            out = np.zeros(data.shape, dtype=np.float32)
            out[positive] = 1.0
        elif self.mode == "clipped_indicator":
            out = np.zeros(data.shape, dtype=np.float32)
            clipped = positive & (data <= self.calibration_offset)
            out[clipped] = 1.0
        else:  # pragma: no cover - guarded in __init__
            raise RuntimeError(self.mode)

        if mask is not None:
            return np.ma.masked_array(out, mask=mask, copy=False)
        return out


def extract_calibrated_year(
    raster_path: Path,
    year: int,
    config: dict,
    pixel_offset: float,
) -> pd.DataFrame:
    """Extract corrected mean and calibration diagnostics for one year."""
    nuts, nuts_release = load_nuts(config, year)
    id_field = config["nuts"]["id_field"]
    expected_crs = config["nuts"]["crs"]

    validate_raster(raster_path, expected_crs)

    print(
        f"Calibrating {year} on NUTS {nuts_release}: {raster_path.name}",
        flush=True,
    )

    with rasterio.open(raster_path) as ds:
        sources = [
            RasterioRasterSource(ds, 1, name="raw_check"),
            CalibrationRasterSource(
                ds,
                offset=pixel_offset,
                mode="corrected",
                name="corrected",
            ),
            CalibrationRasterSource(
                ds,
                offset=pixel_offset,
                mode="raw_positive_indicator",
                name="raw_positive",
            ),
            CalibrationRasterSource(
                ds,
                offset=pixel_offset,
                mode="clipped_indicator",
                name="clipped",
            ),
        ]

        result = exact_extract(
            sources,
            nuts,
            ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=[id_field],
            output="pandas",
            progress=True,
        )

    expected_cols = {
        "raw_check_mean",
        "corrected_mean",
        "raw_positive_mean",
        "clipped_mean",
        id_field,
    }
    missing = expected_cols - set(result.columns)
    if missing:
        raise RuntimeError(
            "Unexpected exactextract output while applying pixel calibration. "
            f"Missing {sorted(missing)}; found {list(result.columns)}"
        )

    result = result.rename(
        columns={
            "raw_check_mean": "radiance_mean_raw_check",
            "corrected_mean": "radiance_mean_corrected",
            "raw_positive_mean": "calibration_raw_positive_area_fraction",
            "clipped_mean": "calibration_clipped_area_fraction",
        }
    )
    result.insert(0, "year", int(year))
    result["nuts_release"] = int(nuts_release)

    return result


def validate_calibrated_rows(
    rows: pd.DataFrame,
    pixel_offset: float,
) -> None:
    """Fail loudly if the transformed extraction violates correction bounds."""
    post = rows[rows["calibration_correction_applied"]].copy()
    if post.empty:
        return

    tolerance = 2e-6

    raw_check_error = np.abs(
        post["radiance_mean_raw_check"] - post["radiance_mean_raw"]
    )
    raw_check_limit = tolerance * np.maximum(
        1.0,
        np.abs(post["radiance_mean_raw"]),
    )
    if (raw_check_error > raw_check_limit).any():
        bad = post.loc[
            raw_check_error > raw_check_limit,
            [
                "year",
                "NUTS_ID",
                "radiance_mean_raw",
                "radiance_mean_raw_check",
            ],
        ].head(10)
        raise RuntimeError(
            "Raw re-extraction does not match the stored production raw mean. "
            "This must be resolved before accepting calibration.\n"
            + bad.to_string(index=False)
        )

    if (
        post["radiance_mean_corrected"]
        > post["radiance_mean_raw"] + tolerance
    ).any():
        raise RuntimeError(
            "Pixel calibration increased at least one NUTS3 mean."
        )

    if (post["calibration_effective_offset"] < -tolerance).any():
        raise RuntimeError(
            "Negative effective NUTS3 calibration offset detected."
        )

    theoretical_max = (
        pixel_offset * post["calibration_raw_positive_area_fraction"]
    )
    if (
        post["calibration_effective_offset"]
        > theoretical_max + tolerance
    ).any():
        raise RuntimeError(
            "Effective NUTS3 correction exceeds pixel_offset × positive-area "
            "fraction. Check the transformed RasterSource."
        )

    if (
        post["calibration_clipped_area_fraction"]
        > post["calibration_raw_positive_area_fraction"] + tolerance
    ).any():
        raise RuntimeError(
            "Clipped area fraction exceeds raw positive area fraction."
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")
    cal = config["calibration"]

    if not cal["enabled"]:
        raise ValueError("Production calibration is disabled in configuration.")

    if cal["method"] != "global_pixel_additive_offset_clipped_zero":
        raise ValueError(
            "apply_radiance_calibration.py requires "
            "calibration.method=global_pixel_additive_offset_clipped_zero; "
            f"found {cal['method']!r}"
        )

    start_year = int(cal["apply_from_year"])
    pixel_offset = float(cal["pixel_offset_nW_cm2_sr"])

    if pixel_offset <= 0:
        raise ValueError("calibration.pixel_offset_nW_cm2_sr must be > 0")

    raw_path = processed / config["outputs"]["radiance_by_nuts3"]
    radiance = pd.read_parquet(raw_path).copy()

    required = {
        "year",
        "NUTS_ID",
        "CNTR_CODE",
        "nuts_release",
        "radiance_mean",
    }
    missing = required - set(radiance.columns)
    if missing:
        raise ValueError(
            f"Raw radiance table is missing required columns: {sorted(missing)}"
        )

    first_year = int(radiance["year"].min())
    last_year = int(radiance["year"].max())
    post_years = list(range(max(start_year, first_year), last_year + 1))

    rasters = discover_rasters(config, post_years)

    corrected_frames = []
    for year in post_years:
        corrected_frames.append(
            extract_calibrated_year(
                rasters[year],
                year,
                config,
                pixel_offset,
            )
        )

    corrected = pd.concat(corrected_frames, ignore_index=True)

    join_keys = ["year", "NUTS_ID", "nuts_release"]
    if corrected.duplicated(join_keys).any():
        raise RuntimeError(
            "Calibrated extraction produced duplicate year × NUTS_ID rows."
        )

    out = radiance.copy()
    out["radiance_mean_raw"] = out["radiance_mean"].astype(float)

    out = out.merge(
        corrected,
        on=join_keys,
        how="left",
        validate="one_to_one",
    )

    post = out["year"] >= start_year
    missing_post = post & out["radiance_mean_corrected"].isna()
    if missing_post.any():
        bad = out.loc[
            missing_post,
            ["year", "NUTS_ID", "CNTR_CODE", "nuts_release"],
        ].head(20)
        raise RuntimeError(
            "Post-2017 rows are missing calibrated raster extraction values.\n"
            + bad.to_string(index=False)
        )

    # Before the calibration change, corrected == raw by definition.
    out.loc[~post, "radiance_mean_corrected"] = out.loc[
        ~post,
        "radiance_mean_raw",
    ]
    out.loc[~post, "radiance_mean_raw_check"] = out.loc[
        ~post,
        "radiance_mean_raw",
    ]
    out.loc[~post, "calibration_raw_positive_area_fraction"] = np.nan
    out.loc[~post, "calibration_clipped_area_fraction"] = 0.0

    out["calibration_correction_applied"] = post
    out["calibration_pixel_offset"] = pixel_offset
    out["calibration_effective_offset"] = (
        out["radiance_mean_raw"] - out["radiance_mean_corrected"]
    )

    denom = (
        pixel_offset * out["calibration_raw_positive_area_fraction"]
    )
    out["calibration_effective_fraction_of_theoretical_max"] = np.where(
        post & np.isfinite(denom) & (denom > 0),
        out["calibration_effective_offset"] / denom,
        np.nan,
    )
    out["calibration_effective_pct_of_raw_mean"] = np.where(
        post & (out["radiance_mean_raw"] > 0),
        100.0
        * out["calibration_effective_offset"]
        / out["radiance_mean_raw"],
        np.nan,
    )

    out["calibration_method"] = cal["method"]
    out["calibration_scope"] = cal["scope"]
    out["calibration_pixel_rule"] = cal["pixel_rule"]
    out["calibration_apply_from_year"] = start_year
    out["calibration_version"] = cal["version"]
    out["calibration_evidence"] = cal["evidence"]

    validate_calibrated_rows(out, pixel_offset)

    corrected_path = (
        processed / config["outputs"]["radiance_by_nuts3_calibrated"]
    )
    metadata_path = processed / config["outputs"]["calibration_offsets"]

    out.to_parquet(corrected_path, index=False)
    out.to_csv(corrected_path.with_suffix(".csv"), index=False)

    metadata = pd.DataFrame(
        [
            {
                "calibration_version": cal["version"],
                "method": cal["method"],
                "scope": cal["scope"],
                "apply_from_year": start_year,
                "pixel_offset_nW_cm2_sr": pixel_offset,
                "pixel_rule": cal["pixel_rule"],
                "selection_note": cal["selection_note"],
                "evidence": cal["evidence"],
            }
        ]
    )
    metadata.to_csv(metadata_path, index=False)

    affected = out[out["calibration_correction_applied"]].copy()

    print()
    print("PRODUCTION VIIRS CALIBRATION")
    print(f"Method: {cal['method']}")
    print(f"Scope: {cal['scope']}")
    print(f"Apply from: {start_year}")
    print(
        f"Pixel offset: {pixel_offset:.6f} nW cm-2 sr-1 "
        "(positive pixels only; clipped at zero)"
    )
    print(f"Corrected NUTS3 rows: {len(affected):,}")

    if not affected.empty:
        print()
        print("Effective NUTS3 correction:")
        print(
            affected["calibration_effective_offset"]
            .describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9])
            .to_string()
        )
        print()
        print("Effective correction as % of raw NUTS3 mean:")
        print(
            affected["calibration_effective_pct_of_raw_mean"]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9])
            .to_string()
        )
        print()
        print("Raw-positive area fraction:")
        print(
            affected["calibration_raw_positive_area_fraction"]
            .describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9])
            .to_string()
        )
        print()
        print("Area fraction clipped to zero:")
        print(
            affected["calibration_clipped_area_fraction"]
            .describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9])
            .to_string()
        )

    print()
    print(f"Wrote: {corrected_path}")
    print(f"Wrote: {corrected_path.with_suffix('.csv')}")
    print(f"Wrote: {metadata_path}")


if __name__ == "__main__":
    main()
