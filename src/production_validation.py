#!/usr/bin/env python3
"""Shared validation for the production radiance/population analysis table.

This module centralises the checks required by historical-response and future
projection stages so that they fail early if they are run against stale or
internally inconsistent production data.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


REQUIRED_PRODUCTION_COLUMNS = {
    "year",
    "NUTS_ID",
    "CNTR_CODE",
    "nuts_release",
    "radiance_mean_raw",
    "radiance_mean_corrected",
    "radiance_mean",
    "population_density",
    "calibration_correction_applied",
    "calibration_pixel_offset",
    "calibration_method",
    "calibration_scope",
    "calibration_version",
}


def validate_production_merged(
    data: pd.DataFrame,
    config: dict,
    *,
    required_years: list[int] | tuple[int, ...] | None = None,
) -> dict:
    """Validate calibration provenance, geography and row uniqueness.

    The function raises ValueError if the merged production table is
    inconsistent with the currently configured calibration or NUTS release
    mapping. It returns a compact provenance summary when all checks pass.
    """
    missing = REQUIRED_PRODUCTION_COLUMNS - set(data.columns)
    if missing:
        raise ValueError(
            "Production merged table is missing required columns: "
            f"{sorted(missing)}"
        )

    d = data.copy()
    d["year"] = pd.to_numeric(d["year"], errors="raise").astype(int)
    d["nuts_release"] = pd.to_numeric(
        d["nuts_release"], errors="raise"
    ).astype(int)

    duplicate = d.duplicated(
        ["year", "NUTS_ID", "nuts_release"], keep=False
    )
    if duplicate.any():
        examples = d.loc[
            duplicate,
            ["year", "NUTS_ID", "CNTR_CODE", "nuts_release"],
        ].head(10)
        raise ValueError(
            "Duplicate year x NUTS_ID x NUTS-release rows in production "
            "table. Examples:\n"
            + examples.to_string(index=False)
        )

    if required_years is not None:
        required = {int(y) for y in required_years}
        observed = set(d["year"].unique())
        missing_years = sorted(required - observed)
        if missing_years:
            raise ValueError(
                f"Production table is missing required years: {missing_years}"
            )

    release_by_year = {
        int(k): int(v)
        for k, v in config["nuts"]["release_by_year"].items()
    }
    expected_release = d["year"].map(release_by_year)

    if expected_release.isna().any():
        years = sorted(
            d.loc[expected_release.isna(), "year"].unique().tolist()
        )
        raise ValueError(
            f"No configured NUTS release for production years: {years}"
        )

    bad_release = d["nuts_release"] != expected_release.astype(int)
    if bad_release.any():
        examples = d.loc[
            bad_release,
            ["year", "NUTS_ID", "nuts_release"],
        ].head(10)
        raise ValueError(
            "Production table contains rows on the wrong configured NUTS "
            "release. Examples:\n"
            + examples.to_string(index=False)
        )

    cal = config["calibration"]
    start_year = int(cal["apply_from_year"])
    expected_offset = float(cal["pixel_offset_nW_cm2_sr"])
    expected_method = str(cal["method"])
    expected_scope = str(cal["scope"])
    expected_version = str(cal["version"])

    raw = pd.to_numeric(d["radiance_mean_raw"], errors="coerce")
    corrected = pd.to_numeric(
        d["radiance_mean_corrected"], errors="coerce"
    )
    alias = pd.to_numeric(d["radiance_mean"], errors="coerce")
    offset = pd.to_numeric(
        d["calibration_pixel_offset"], errors="coerce"
    )

    finite_rad = np.isfinite(raw) & np.isfinite(corrected)
    tolerance = 2e-6 * np.maximum(
        1.0, np.abs(raw.fillna(0.0).to_numpy(dtype=float))
    )

    alias_error = np.abs(
        alias.to_numpy(dtype=float)
        - corrected.to_numpy(dtype=float)
    )
    bad_alias = finite_rad.to_numpy() & (
        ~np.isfinite(alias.to_numpy(dtype=float))
        | (alias_error > tolerance)
    )
    if bad_alias.any():
        raise ValueError(
            "radiance_mean is not an exact production alias of "
            "radiance_mean_corrected."
        )

    post = d["year"] >= start_year
    pre = ~post

    if not d.loc[
        post, "calibration_correction_applied"
    ].astype(bool).all():
        raise ValueError(
            "At least one post-calibration row is not flagged as corrected."
        )

    if d.loc[
        pre, "calibration_correction_applied"
    ].astype(bool).any():
        raise ValueError(
            "At least one pre-calibration row is incorrectly flagged as "
            "corrected."
        )

    for column, expected in [
        ("calibration_method", expected_method),
        ("calibration_scope", expected_scope),
        ("calibration_version", expected_version),
    ]:
        bad = d[column].astype(str) != expected
        if bad.any():
            values = sorted(
                d.loc[bad, column].astype(str).unique().tolist()
            )
            raise ValueError(
                f"{column} does not match current configuration "
                f"({expected!r}); observed incompatible values: "
                f"{values[:10]}"
            )

    bad_offset = (
        ~np.isfinite(offset.to_numpy(dtype=float))
        | ~np.isclose(
            offset.to_numpy(dtype=float),
            expected_offset,
            rtol=0.0,
            atol=1e-12,
        )
    )
    if bad_offset.any():
        raise ValueError(
            "calibration_pixel_offset does not match the configured "
            f"production value {expected_offset}."
        )

    pre_error = np.abs(
        corrected.to_numpy(dtype=float)
        - raw.to_numpy(dtype=float)
    )
    if (
        pre.to_numpy()
        & finite_rad.to_numpy()
        & (pre_error > tolerance)
    ).any():
        raise ValueError(
            "Pre-2017 corrected radiance differs from raw radiance."
        )

    corrected_arr = corrected.to_numpy(dtype=float)
    raw_arr = raw.to_numpy(dtype=float)
    post_arr = post.to_numpy()
    finite_arr = finite_rad.to_numpy()

    if (
        post_arr
        & finite_arr
        & (corrected_arr < -tolerance)
    ).any():
        raise ValueError(
            "Post-2017 pixel-calibrated NUTS3 mean is negative; the "
            "production pixel rule clips positive pixels at zero."
        )

    if (
        post_arr
        & finite_arr
        & (corrected_arr > raw_arr + tolerance)
    ).any():
        raise ValueError(
            "Post-2017 calibrated radiance exceeds raw radiance."
        )

    if "calibration_effective_offset" in d.columns:
        effective = pd.to_numeric(
            d["calibration_effective_offset"], errors="coerce"
        ).to_numpy(dtype=float)
        derived = raw_arr - corrected_arr
        bad = (
            finite_arr
            & np.isfinite(effective)
            & (np.abs(effective - derived) > tolerance)
        )
        if bad.any():
            raise ValueError(
                "Stored calibration_effective_offset is inconsistent with "
                "raw - corrected radiance."
            )

    pop = pd.to_numeric(
        d["population_density"], errors="coerce"
    ).to_numpy(dtype=float)
    if (np.isfinite(pop) & (pop < 0)).any():
        raise ValueError("Negative population density detected.")

    return {
        "rows": int(len(d)),
        "first_year": int(d["year"].min()),
        "last_year": int(d["year"].max()),
        "n_nuts3_release_keys": int(
            d[["NUTS_ID", "nuts_release"]]
            .drop_duplicates()
            .shape[0]
        ),
        "n_countries": int(d["CNTR_CODE"].nunique()),
        "calibration_version": expected_version,
        "calibration_method": expected_method,
        "calibration_scope": expected_scope,
        "calibration_apply_from_year": start_year,
        "calibration_pixel_offset_nW_cm2_sr": expected_offset,
    }
