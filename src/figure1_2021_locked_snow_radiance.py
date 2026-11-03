#!/usr/bin/env python3
"""Extract the LOCKED native-pixel snow-corrected 2021 VIIRS radiance once.

This is the radiance counterpart of the already-approved 2021 PRC overlay,
not a new correction. Uses:
  - config/viirs_2021_snow_correction.yaml (status must be locked)
  - original EOG annual 2021 raster
  - existing native-grid Europe patch, band 1 (RAW radiance, pre 0.235)
  - the existing calibration −0.235 from 2017 onward exactly once
  - fixed NUTS2016 and exactextract spherical-area weighted means.
The original VIIRS image and locked patch are NEVER modified.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import yaml
from rasterio.windows import Window
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource

from apply_radiance_calibration import CalibrationRasterSource
from extract_nuts3_radiance import discover_rasters, validate_raster
from historical_prc_pixel_nuts3 import load_nuts
from prc_2021_locked17_overlay_analysis import verify
from project_config import PROJECT_ROOT, load_config


class SnowCorrectedRadianceSource(RasterioRasterSource):
    """Original 2021 raw -> locked raw patch overlay -> existing pixel calibration."""

    def __init__(self, ds, patch, xoff: int, yoff: int, offset: float):
        super().__init__(ds, 1, name="snow")
        self.patch = patch
        self.xoff = int(xoff)
        self.yoff = int(yoff)
        self.offset = float(offset)

    def read_window(self, x0, y0, nx, ny):
        x0, y0, nx, ny = map(int, (x0, y0, nx, ny))
        raw = np.ma.array(
            RasterioRasterSource.read_window(self, x0, y0, nx, ny), copy=True
        )
        left = max(x0, self.xoff)
        right = min(x0 + nx, self.xoff + self.patch.width)
        top = max(y0, self.yoff)
        bottom = min(y0 + ny, self.yoff + self.patch.height)
        if right > left and bottom > top:
            crop = self.patch.read(
                1,
                window=Window(left - self.xoff, top - self.yoff,
                              right - left, bottom - top),
                masked=True,
            )
            target = raw[top - y0:bottom - y0, left - x0:right - x0]
            corrected_raw = np.asarray(crop.data, dtype=np.float64)
            original_raw = np.asarray(target.data, dtype=np.float64)
            replace = (
                ~np.ma.getmaskarray(crop)
                & ~np.ma.getmaskarray(target)
                & np.isfinite(corrected_raw)
                & (corrected_raw >= 0)
                & np.isfinite(original_raw)
                & (original_raw > 0)
            )
            if np.any(corrected_raw[replace] > original_raw[replace] + 1e-4):
                raise ValueError("Locked snow overlay unexpectedly increases raw VIIRS")
            target.data[replace] = corrected_raw[replace]
        # EXACT same transform as the established CalibrationRasterSource:
        # retain zero/background pixels; subtract offset only where >0.
        data = np.asarray(raw.data, dtype=np.float32).copy()
        mask = np.ma.getmaskarray(raw)
        positive = (~mask) & np.isfinite(data) & (data > 0)
        data[positive] = np.maximum(data[positive] - self.offset, 0.0)
        return np.ma.masked_array(data, mask=mask, copy=False)


def _lock_info(config: dict):
    lock_path = PROJECT_ROOT / "config/viirs_2021_snow_correction.yaml"
    lock = yaml.safe_load(lock_path.read_text(encoding="utf-8"))
    if lock.get("status") != "locked" or int(lock.get("applies_to_year", 0)) != 2021:
        raise ValueError("Native 2021 snow correction must be LOCKED")
    offset = float(config["calibration"]["pixel_offset_nW_cm2_sr"])
    year_from = int(config["calibration"]["apply_from_year"])
    if year_from > 2021 or abs(offset - float(lock["downstream_prc_2017_onward_offset"])) > 1e-9:
        raise ValueError("2021 calibration offset disagrees with locked snow config")
    if lock.get("raw_annual_radiance_offset_applied_in_patch") is not False:
        raise ValueError("Locked patch is not labelled as pre-calibration raw radiance")
    return lock, offset


def extract_locked_2021(*, self_test: bool = False) -> tuple[pd.DataFrame, dict] | None:
    if self_test:
        # Checks that raw patch is applied before baseline calibration.
        from rasterio.io import MemoryFile
        from rasterio.transform import from_origin
        with MemoryFile() as m, MemoryFile() as p:
            a = np.array([[0., .4], [2., 5.]], dtype="float32")
            b = np.array([[.2, 1.]], dtype="float32")
            with m.open(driver="GTiff", height=2, width=2, count=1,
                        dtype="float32", crs="EPSG:4326",
                        transform=from_origin(0, 2, 1, 1)) as ds:
                ds.write(a, 1)
            with p.open(driver="GTiff", height=1, width=2, count=1,
                        dtype="float32", crs="EPSG:4326", nodata=np.nan,
                        transform=from_origin(0, 1, 1, 1)) as ds:
                ds.write(np.array([[1., 2.]], dtype="float32"), 1)
            with m.open() as ds, p.open() as patch:
                got = SnowCorrectedRadianceSource(
                    ds, patch, 0, 1, .235
                ).read_window(0, 0, 2, 2)
                assert np.allclose(got, [[0., .165], [.765, 1.765]])
        print("Locked 2021 raw overlay before single calibration: PASS", flush=True)
        return None

    cfg = load_config()
    lock, offset = _lock_info(cfg)
    patch_path = PROJECT_ROOT / lock["patch"]
    manifest_path = PROJECT_ROOT / lock["patch_summary"]
    if not patch_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError(
            f"Locked patch or manifest missing: {patch_path} / {manifest_path}"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    nuts = load_nuts(cfg)
    if len(nuts) != 1522:
        raise ValueError("Expected frozen 1,522-region NUTS2016 geometry")
    original_path = discover_rasters(cfg, [2021])[2021]
    validate_raster(original_path, cfg["nuts"]["crs"])
    with rasterio.open(original_path) as ds, rasterio.open(patch_path) as patch:
        col, row = verify(ds, patch, manifest)
        sources = [
            CalibrationRasterSource(ds, offset=offset, mode="corrected", name="original"),
            SnowCorrectedRadianceSource(ds, patch, col, row, offset),
        ]
        result = exact_extract(
            sources, nuts,
            ["mean(coverage_weight=area_spherical_m2)",
             "count(coverage_weight=area_spherical_km2)"],
            include_cols=["NUTS_ID"], output="pandas", progress=True
        )
    expected = {
        "NUTS_ID", "original_mean", "snow_mean", "original_count", "snow_count"
    }
    if not expected.issubset(result):
        raise ValueError(f"Unexpected 2021 extraction columns: {list(result)}")
    if len(result) != 1522 or result.NUTS_ID.duplicated().any():
        raise ValueError("2021 corrected radiance did not cover 1,522 distinct NUTS3")
    if not np.all(np.isfinite(result[["original_mean", "snow_mean"]].to_numpy(float))):
        raise ValueError("Nonfinite 2021 radiance")
    if (result.snow_mean > result.original_mean + 2e-6).any():
        raise ValueError("2021 correction increased radiance")
    if not np.allclose(result.original_count, result.snow_count, atol=.05, rtol=1e-7):
        raise ValueError("Snow correction changed valid-pixel area support")
    out = result.rename(columns={
        "original_mean": "radiance_2021_original_calibrated",
        "snow_mean": "radiance_2021_snow_corrected",
        "original_count": "valid_area_km2_2021_original",
        "snow_count": "valid_area_km2_2021_snow",
    })
    out["snow_correction_nw"] = (
        out.radiance_2021_snow_corrected - out.radiance_2021_original_calibrated
    )
    audit = {
        "status": "locked",
        "year": 2021,
        "correction_model": lock["model_id"],
        "method": "locked native-pixel raw patch, then 2017+ pixel offset once",
        "original_annual_vnl": str(original_path.resolve()),
        "patch": str(patch_path.resolve()),
        "patch_manifest": str(manifest_path.resolve()),
        "patch_bytes": patch_path.stat().st_size,
        "patch_mtime_ns": patch_path.stat().st_mtime_ns,
        "original_annual_bytes": original_path.stat().st_size,
        "original_annual_mtime_ns": original_path.stat().st_mtime_ns,
        "pixel_offset_nw_cm2_sr": offset,
        "n_nuts3": len(out),
        "all_corrected_means_no_greater_than_original": True,
        "valid_coverage_preserved": True,
        "n_regions_decreased": int((out.snow_correction_nw < -1e-6).sum()),
    }
    return out, audit


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", default=
        "data/processed/publication_figures/fig1_2021_locked_snow_radiance.csv")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    extract_locked_2021(self_test=True)
    if args.self_test:
        return
    rows, audit = extract_locked_2021()
    target = Path(args.output)
    if not target.is_absolute():
        target = PROJECT_ROOT / target
    target.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(target, index=False)
    target.with_suffix(".json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    print("2021 locked radiance:", target, flush=True)
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
