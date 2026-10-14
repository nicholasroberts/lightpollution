#!/usr/bin/env python3
"""Test whether snow explains the January/February 2021 VIIRS anomaly.

This stage is deliberately literature-led. It adapts the approach of
Levin & Zhang (Remote Sensing of Environment 190, 2017, 366-382):
- mean VIIRS radiance within fixed spatial units;
- MOD10CM mean snow cover;
- VIIRS cloud-free coverage;
- percent lit area above published radiance thresholds
  (2, 5, 10, 25, 50, 100, 250 nW cm-2 sr-1);
- Spearman rank associations; and
- winter-versus-summer contrasts.

The primary causal diagnostic is tailored to the anomaly identified in this
project. For January, February and March, it compares:

    2021 - mean(2020, 2022)

for matched-support VIIRS radiance and MOD10CM snow cover. February is the
pre-specified primary month, January is supportive, and March is the negative
control. Country fixed effects, latitude, baseline brightness and cloud-free
coverage anomaly are included in robust sensitivity models. A snow x baseline
brightness interaction tests the physical prediction that snow reflection
should have a larger absolute effect where more artificial light is available.

MOD10CM Collection 6.1 values 0-100 are valid monthly mean percent snow cover;
211 (night), 250 (cloud), 253 (no decision), 254 (water) and 255 (fill) are
treated as invalid. HDF4 is converted once with system gdal_translate to a
cached GeoTIFF because the project Rasterio wheel does not include HDF4.

No annual VIIRS production data are modified and no correction is applied.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from scipy.stats import spearmanr
import statsmodels.formula.api as smf
from shapely.geometry import box

from diagnose_fixed_nuts2016_2021_reversal import add_latitude, load_locked_panel
from diagnose_fixed_nuts2016_mask_occupancy import attach_fixed_brightness
from project_config import PROJECT_ROOT, configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import load_fixed_nuts


YEARS = (2020, 2021, 2022)
PRIMARY_MONTHS = (1, 2, 3)
LEVIN_MONTHS = (1, 7)
AUX_VIIRS_MONTHS = (1, 2, 3, 7)
MONTH_NAMES = ("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
LIT_THRESHOLDS = (2, 5, 10, 25, 50, 100, 250)
SNOW_CLEANING_VERSION = "mod10cm_mask_special_gt100_v2"

HDF_RE = re.compile(r"^MOD10CM\.A(?P<year>\d{4})(?P<doy>\d{3})\.061\..*\.hdf$")
VIIRS_RE = re.compile(
    r"^SVDNB_(?P<satellite>npp|j01)_(?P<start>\d{8})-(?P<end>\d{8})_"
    r"(?P<tile>[^_]+)_vcmslcfg_v10_.*\.avg_rade9h\.tif$",
    re.IGNORECASE,
)


def month_from_year_doy(year: int, doy: int) -> int:
    return int((pd.Timestamp(f"{year}-01-01") + pd.Timedelta(days=doy - 1)).month)


def discover_mod10cm(raw_dir: Path) -> dict[tuple[int, int], Path]:
    found: dict[tuple[int, int], Path] = {}
    for p in sorted(raw_dir.glob("MOD10CM*.hdf")):
        m = HDF_RE.match(p.name)
        if not m:
            continue
        year = int(m.group("year"))
        if year not in YEARS:
            continue
        month = month_from_year_doy(year, int(m.group("doy")))
        key = (year, month)
        if key in found:
            raise RuntimeError(f"Duplicate MOD10CM month {key}: {found[key]} and {p}")
        found[key] = p

    expected = {(y, m) for y in YEARS for m in range(1, 13)}
    missing = sorted(expected - set(found))
    allowed_missing = [(2022, 10)]
    if missing != allowed_missing:
        raise RuntimeError(
            "Unexpected MOD10CM monthly inventory. Expected only 2022-10 to be "
            f"missing; found missing={missing}."
        )
    if len(found) != 35:
        raise RuntimeError(f"Expected 35 MOD10CM files; found {len(found)}.")
    return found


def ensure_snow_geotiff(hdf: Path, cache_dir: Path, expected_crs: str) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{hdf.stem}.Snow_Cover_Monthly_CMG.tif"
    meta_path = out.with_suffix(".json")
    signature = {
        "source": hdf.name,
        "size_bytes": int(hdf.stat().st_size),
        "cleaning_version": SNOW_CLEANING_VERSION,
    }

    if out.exists() and meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            meta = {}
        if meta == signature:
            with rasterio.open(out) as ds:
                if (
                    ds.count == 1
                    and ds.crs is not None
                    and ds.crs.to_string() == expected_crs
                    and ds.width == 7200
                    and ds.height == 3600
                ):
                    return out

    gdal_translate = shutil.which("gdal_translate")
    if not gdal_translate:
        raise RuntimeError(
            "gdal_translate is required for MOD10CM HDF4. Install with: "
            "sudo apt install gdal-bin gdal-plugins"
        )

    subdataset = (
        f'HDF4_EOS:EOS_GRID:"{hdf.resolve()}":'
        "MOD_CMG_Snow_5km:Snow_Cover_Monthly_CMG"
    )
    tmp = out.with_suffix(".tif.part")
    if tmp.exists():
        tmp.unlink()

    print(f"Converting MOD10CM snow layer: {hdf.name}", flush=True)
    # MOD10CM is a regular 0.05-degree longitude/latitude CMG.  The HDF-EOS
    # driver exposes its legacy GCTP geographic grid as an "unknown datum based
    # on Clarke 1866" CRS even though the product coordinates themselves are
    # the global lon/lat grid (-180..180, -90..90).  Assign EPSG:4326 to the
    # cached GeoTIFF without resampling so that the grid can be overlaid with
    # the project's NUTS geometries.  Bounds/shape are validated immediately
    # below, so this cannot silently accept a non-CMG raster.
    subprocess.run(
        [
            gdal_translate,
            "-q",
            "-of", "GTiff",
            "-a_srs", expected_crs,
            "-co", "COMPRESS=DEFLATE",
            "-co", "TILED=YES",
            "-co", "BIGTIFF=IF_SAFER",
            subdataset,
            str(tmp),
        ],
        check=True,
    )
    # IMPORTANT: MOD10CM is uint8 but values above 100 are categorical
    # codes (211 night, 250 cloud, 253 no decision, 254 water, 255 fill),
    # not snow percentages.  Burn every >100 code to a real raster nodata
    # value before exactextract.  This is deliberately physical rather than
    # relying on a Python masked-array subclass, because exactextract/GDAL
    # must see the nodata mask at the raster-source level.
    with rasterio.open(tmp, "r+") as ds:
        ds.nodata = 255
        for _, window in ds.block_windows(1):
            arr = ds.read(1, window=window)
            bad = arr > 100
            if np.any(bad):
                arr[bad] = 255
                ds.write(arr, 1, window=window)

    tmp.replace(out)

    with rasterio.open(out) as ds:
        if ds.count != 1:
            raise RuntimeError(f"Unexpected band count in converted MOD10CM: {out}")
        if ds.crs is None or ds.crs.to_string() != expected_crs:
            raise RuntimeError(
                f"Unexpected MOD10CM CRS {ds.crs}; expected {expected_crs}"
            )
        if (ds.width, ds.height) != (7200, 3600):
            raise RuntimeError(
                f"Unexpected MOD10CM grid {(ds.width, ds.height)}; "
                "expected (7200, 3600)."
            )
        if ds.nodata != 255:
            raise RuntimeError(f"MOD10CM cleaned GeoTIFF has nodata={ds.nodata}, expected 255.")
        b = ds.bounds
        if not (
            abs(b.left + 180) < 0.1 and abs(b.right - 180) < 0.1
            and abs(b.bottom + 90) < 0.1 and abs(b.top - 90) < 0.1
        ):
            raise RuntimeError(f"Unexpected MOD10CM bounds: {b}")

        # Hard invariant: every non-nodata snow value must be a percentage.
        arr = ds.read(1, masked=True)
        if arr.count() and (float(arr.min()) < 0 or float(arr.max()) > 100):
            raise RuntimeError(
                f"MOD10CM cleaning failed for {out.name}: valid range "
                f"{float(arr.min())}..{float(arr.max())}"
            )

    meta_path.write_text(json.dumps(signature, indent=2), encoding="utf-8")
    return out


def _as_data_valid(arr):
    if np.ma.isMaskedArray(arr):
        data = np.asarray(arr.data)
        valid = ~np.ma.getmaskarray(arr)
    else:
        data = np.asarray(arr)
        valid = np.ones(data.shape, dtype=bool)
    valid &= np.isfinite(data)
    return data, valid


class SnowValueSource(RasterioRasterSource):
    """MOD10CM valid snow percent only (0-100); special codes are masked."""

    def __init__(self, ds, *, name="snow"):
        super().__init__(ds, 1, name=name)
        self.src = RasterioRasterSource(ds, 1, name=f"{name}_raw")

    def read_window(self, x0, y0, nx, ny):
        arr = self.src.read_window(x0, y0, nx, ny)
        data, valid = _as_data_valid(arr)
        valid &= (data >= 0) & (data <= 100)
        return np.ma.masked_array(
            data.astype(np.float32), mask=~valid, copy=False
        )


class SnowValidSource(RasterioRasterSource):
    """1 where MOD10CM has a valid 0-100 snow value; 0 for special codes."""

    def __init__(self, ds, *, name="snow_valid"):
        super().__init__(ds, 1, name=name)
        self.src = RasterioRasterSource(ds, 1, name=f"{name}_raw")

    def read_window(self, x0, y0, nx, ny):
        arr = self.src.read_window(x0, y0, nx, ny)
        data, valid = _as_data_valid(arr)
        out = (valid & (data >= 0) & (data <= 100)).astype(np.float32)
        return out


class CommonSupportThresholdSource(RasterioRasterSource):
    """Binary VIIRS threshold source masked outside month-specific common support."""

    def __init__(self, avg_ds, support_ds, threshold: float, *, name: str):
        super().__init__(avg_ds, 1, name=name)
        self.avg = RasterioRasterSource(avg_ds, 1, name=f"{name}_avg")
        self.support = RasterioRasterSource(support_ds, 1, name=f"{name}_support")
        self.threshold = float(threshold)

    def read_window(self, x0, y0, nx, ny):
        avg = self.avg.read_window(x0, y0, nx, ny)
        sup = self.support.read_window(x0, y0, nx, ny)
        avg_data, avg_valid = _as_data_valid(avg)
        sup_data, sup_valid = _as_data_valid(sup)
        valid = avg_valid & sup_valid & (sup_data == 1)
        out = (avg_data > self.threshold).astype(np.float32)
        return np.ma.masked_array(out, mask=~valid, copy=False)


class CommonSupportValueSource(RasterioRasterSource):
    """Numeric source masked outside month-specific common support."""

    def __init__(self, value_ds, support_ds, *, name: str):
        super().__init__(value_ds, 1, name=name)
        self.value = RasterioRasterSource(value_ds, 1, name=f"{name}_value")
        self.support = RasterioRasterSource(support_ds, 1, name=f"{name}_support")

    def read_window(self, x0, y0, nx, ny):
        value = self.value.read_window(x0, y0, nx, ny)
        sup = self.support.read_window(x0, y0, nx, ny)
        value_data, value_valid = _as_data_valid(value)
        sup_data, sup_valid = _as_data_valid(sup)
        valid = value_valid & sup_valid & (sup_data == 1)
        return np.ma.masked_array(
            value_data.astype(np.float32), mask=~valid, copy=False
        )


def integration_self_test() -> None:
    # Snow coding test: one row gives equal-area cells at essentially one latitude.
    snow_mem = MemoryFile()
    avg_mem = MemoryFile()
    support_mem = MemoryFile()
    cf_mem = MemoryFile()
    datasets = []
    try:
        transform = from_origin(0, 1, 1, 1)
        snow_ds = snow_mem.open(
            driver="GTiff", height=1, width=4, count=1, dtype="uint8",
            crs="EPSG:4326", transform=transform
        )
        snow_ds.write(np.array([[0, 50, 100, 211]], dtype="uint8"), 1)
        avg_ds = avg_mem.open(
            driver="GTiff", height=1, width=4, count=1, dtype="float32",
            crs="EPSG:4326", transform=transform
        )
        avg_ds.write(np.array([[1, 3, 6, 10]], dtype="float32"), 1)
        support_ds = support_mem.open(
            driver="GTiff", height=1, width=4, count=1, dtype="uint8",
            crs="EPSG:4326", transform=transform
        )
        support_ds.write(np.array([[1, 1, 1, 0]], dtype="uint8"), 1)
        cf_ds = cf_mem.open(
            driver="GTiff", height=1, width=4, count=1, dtype="float32",
            crs="EPSG:4326", transform=transform
        )
        cf_ds.write(np.array([[2, 4, 6, 8]], dtype="float32"), 1)
        datasets = [snow_ds, avg_ds, support_ds, cf_ds]

        gdf = gpd.GeoDataFrame(
            {"NUTS_ID": ["TEST"]},
            geometry=[box(0, 0, 4, 1)],
            crs="EPSG:4326",
        )
        sources = [
            SnowValueSource(snow_ds, name="snow"),
            SnowValidSource(snow_ds, name="snow_valid"),
            CommonSupportThresholdSource(
                avg_ds, support_ds, 2, name="lit_gt_2"
            ),
            CommonSupportValueSource(
                cf_ds, support_ds, name="cf_cvg"
            ),
        ]
        out = exact_extract(
            sources,
            gdf,
            ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=["NUTS_ID"],
            output="pandas",
            progress=False,
        ).iloc[0]
        if not np.isclose(out["snow_mean"], 50.0, atol=1e-4):
            raise RuntimeError("Self-test snow mean failed.")
        if not np.isclose(out["snow_valid_mean"], 0.75, atol=1e-4):
            raise RuntimeError("Self-test snow valid fraction failed.")
        if not np.isclose(out["lit_gt_2_mean"], 2 / 3, atol=1e-4):
            raise RuntimeError("Self-test percent-lit source failed.")
        if not np.isclose(out["cf_cvg_mean"], 4.0, atol=1e-4):
            raise RuntimeError("Self-test common-support cf_cvg mean failed.")
    finally:
        for ds in datasets:
            ds.close()
        snow_mem.close()
        avg_mem.close()
        support_mem.close()
        cf_mem.close()

    # Guard against pandas label-alignment errors in the year-pivot pattern
    # used by the anomaly construction.
    toy = pd.DataFrame({
        "NUTS_ID": ["A", "A", "A", "B", "B", "B"],
        "year": [2020, 2021, 2022, 2020, 2021, 2022],
        "snow": [10.0, 30.0, 20.0, 40.0, 50.0, 60.0],
    })
    p = toy.pivot(index="NUTS_ID", columns="year", values="snow")
    wide = p.index.to_frame(index=False, name="NUTS_ID")
    for year in YEARS:
        wide[f"snow_{year}"] = p[year].to_numpy()
    expected = np.array([10.0, 40.0])
    if not np.allclose(wide["snow_2020"].to_numpy(), expected):
        raise RuntimeError("Self-test year-pivot alignment failed.")

    print("Snow/VIIRS literature-method integration self-test: PASS", flush=True)


def fixed_nuts(config: dict) -> gpd.GeoDataFrame:
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release)
    id_field = config["nuts"]["id_field"]
    return nuts[[id_field, "geometry"]].copy()


def extract_snow_month(
    tif: Path,
    year: int,
    month: int,
    nuts: gpd.GeoDataFrame,
    config: dict,
    cache_dir: Path,
) -> pd.DataFrame:
    id_field = config["nuts"]["id_field"]
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"nuts2016_mod10cm_{year}_{month:02d}.parquet"
    sig = f"{SNOW_CLEANING_VERSION}:{tif.name}:{tif.stat().st_size}"

    if cache.exists():
        d = pd.read_parquet(cache)
        required = {
            id_field, "year", "month", "snow_mean_pct",
            "snow_valid_area_fraction", "source_signature",
        }
        if required.issubset(d.columns):
            if (
                d["source_signature"].dropna().astype(str).unique().tolist() == [sig]
                and set(d[id_field].astype(str)) == set(nuts[id_field].astype(str))
            ):
                return d

    with rasterio.open(tif) as ds:
        result = exact_extract(
            [
                SnowValueSource(ds, name="snow"),
                SnowValidSource(ds, name="snow_valid"),
            ],
            nuts,
            ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=[id_field],
            output="pandas",
            progress=True,
        )

    expected = {id_field, "snow_mean", "snow_valid_mean"}
    missing = expected - set(result.columns)
    if missing:
        raise RuntimeError(f"Unexpected MOD10CM exactextract output: missing {missing}")

    d = result.rename(columns={
        "snow_mean": "snow_mean_pct",
        "snow_valid_mean": "snow_valid_area_fraction",
    })
    d.insert(0, "month", int(month))
    d.insert(0, "year", int(year))
    d["source_signature"] = sig
    d.to_parquet(cache, index=False)
    return d


def discover_viirs_month(
    config: dict, year: int, month: int, tile: str
) -> tuple[Path, Path]:
    root = PROJECT_ROOT / f"data/raw/viirs/monthly_vcmslcfg_{year}/tif"
    matches = []
    for p in root.rglob("*.avg_rade9h.tif"):
        m = VIIRS_RE.match(p.name)
        if not m:
            continue
        if (
            int(m.group("start")[:4]) == year
            and int(m.group("start")[4:6]) == month
            and m.group("tile") == tile
        ):
            matches.append(p)
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one VIIRS VCMSLCFG raster for {year}-{month:02d}; "
            f"found {len(matches)}."
        )
    avg = matches[0]
    cf = Path(str(avg).replace(".avg_rade9h.tif", ".cf_cvg.tif"))
    if not cf.exists():
        raise FileNotFoundError(cf)
    return avg, cf


def extract_viirs_literature_metrics(
    year: int,
    month: int,
    avg: Path,
    cf: Path,
    support: Path,
    nuts: gpd.GeoDataFrame,
    config: dict,
    cache_dir: Path,
) -> pd.DataFrame:
    id_field = config["nuts"]["id_field"]
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"nuts2016_viirs_lit_metrics_{year}_{month:02d}.parquet"
    sig = "|".join(
        f"{p.name}:{p.stat().st_size}" for p in (avg, cf, support)
    )

    if cache.exists():
        d = pd.read_parquet(cache)
        req = {
            id_field, "year", "month", "cf_cvg_mean_common",
            "source_signature",
        } | {f"pct_lit_gt_{t}" for t in LIT_THRESHOLDS}
        if req.issubset(d.columns):
            if (
                d["source_signature"].dropna().astype(str).unique().tolist() == [sig]
                and set(d[id_field].astype(str)) == set(nuts[id_field].astype(str))
            ):
                return d

    avg_ds = rasterio.open(avg)
    cf_ds = rasterio.open(cf)
    support_ds = rasterio.open(support)
    try:
        sources = [
            CommonSupportThresholdSource(
                avg_ds, support_ds, t, name=f"lit_gt_{t}"
            )
            for t in LIT_THRESHOLDS
        ]
        sources.append(
            CommonSupportValueSource(
                cf_ds, support_ds, name="cf_cvg"
            )
        )
        result = exact_extract(
            sources,
            nuts,
            ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=[id_field],
            output="pandas",
            progress=True,
        )
    finally:
        avg_ds.close()
        cf_ds.close()
        support_ds.close()

    rename = {
        f"lit_gt_{t}_mean": f"pct_lit_gt_{t}"
        for t in LIT_THRESHOLDS
    }
    rename["cf_cvg_mean"] = "cf_cvg_mean_common"
    d = result.rename(columns=rename)

    req = {id_field, "cf_cvg_mean_common"} | {
        f"pct_lit_gt_{t}" for t in LIT_THRESHOLDS
    }
    missing = req - set(d.columns)
    if missing:
        raise RuntimeError(
            f"Unexpected VIIRS literature metrics schema: missing {missing}"
        )

    for t in LIT_THRESHOLDS:
        d[f"pct_lit_gt_{t}"] = 100.0 * d[f"pct_lit_gt_{t}"]
    d.insert(0, "month", int(month))
    d.insert(0, "year", int(year))
    d["source_signature"] = sig
    d.to_parquet(cache, index=False)
    return d


def load_viirs_monthly_regional(config: dict) -> pd.DataFrame:
    path = (
        configured_path(config, "processed")
        / "fixed_nuts2016_monthly_2021_localization_regional.csv"
    )
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run future-radiance-monthly-2021-localization first."
        )
    d = pd.read_csv(path)
    needed = {
        "NUTS_ID", "month", "common_support_area_fraction",
        "radiance_2020_common_nw", "radiance_2021_common_nw",
        "radiance_2022_common_nw", "anomaly_2021_vs_flanks_nw",
    }
    missing = needed - set(d.columns)
    if missing:
        raise RuntimeError(f"Monthly VIIRS regional file missing {missing}")
    return d


def viirs_long(viirs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    meta_cols = [
        "NUTS_ID", "month", "common_support_area_fraction",
        "CNTR_CODE", "density_band", "brightness_quintile",
        "baseline_raw_nw", "centroid_latitude", "centroid_longitude",
    ]
    meta_cols = [c for c in meta_cols if c in viirs.columns]
    for year in YEARS:
        x = viirs[meta_cols].copy()
        x["year"] = year
        x["radiance_common_nw"] = viirs[f"radiance_{year}_common_nw"].to_numpy()
        rows.append(x)
    return pd.concat(rows, ignore_index=True)


def in_continental_extent(d: pd.DataFrame, bounds: list[float]) -> pd.Series:
    west, south, east, north = map(float, bounds)
    return (
        d["centroid_longitude"].between(west, east)
        & d["centroid_latitude"].between(south, north)
    )


def spearman_clean(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    z = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(z) < 3 or z["x"].nunique() < 2 or z["y"].nunique() < 2:
        return np.nan, np.nan, len(z)
    r = spearmanr(z["x"], z["y"])
    return float(r.statistic), float(r.pvalue), len(z)


def build_primary_anomaly(
    merged: pd.DataFrame,
    viirs: pd.DataFrame,
    aux: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    min_snow = float(cfg["min_snow_valid_area_fraction"])
    min_common = float(cfg["min_viirs_common_area_fraction"])
    rows = []

    for month in PRIMARY_MONTHS:
        snow_m = merged[merged["month"] == month].copy()
        snow_p = snow_m.pivot(
            index="NUTS_ID", columns="year",
            values=["snow_mean_pct", "snow_valid_area_fraction"]
        )
        missing_years = set(YEARS) - set(snow_p.columns.get_level_values(1))
        if missing_years:
            raise RuntimeError(
                f"Month {month}: snow pivot is missing years {sorted(missing_years)}."
            )
        # Preserve NUTS_ID as a normal column and assign pivoted values by
        # position.  Assigning the Series directly would align its NUTS_ID
        # index against the DataFrame's RangeIndex and silently create NaNs.
        out = snow_p.index.to_frame(index=False, name="NUTS_ID")
        for year in YEARS:
            out[f"snow_{year}_pct"] = snow_p["snow_mean_pct"][year].to_numpy()
            out[f"snow_valid_{year}"] = (
                snow_p["snow_valid_area_fraction"][year].to_numpy()
            )

        out["snow_anomaly_2021_vs_flanks_pp"] = (
            out["snow_2021_pct"]
            - (out["snow_2020_pct"] + out["snow_2022_pct"]) / 2.0
        )
        out["min_snow_valid_fraction"] = out[
            [f"snow_valid_{y}" for y in YEARS]
        ].min(axis=1)

        v = viirs[viirs["month"] == month][[
            "NUTS_ID", "anomaly_2021_vs_flanks_nw",
            "common_support_area_fraction", "CNTR_CODE", "density_band",
            "brightness_quintile", "baseline_raw_nw",
            "centroid_latitude", "centroid_longitude",
        ]].copy()
        out = out.merge(v, on="NUTS_ID", how="inner", validate="one_to_one")

        a = aux[aux["month"] == month]
        cf = a.pivot(index="NUTS_ID", columns="year", values="cf_cvg_mean_common")
        missing_cf_years = set(YEARS) - set(cf.columns)
        if missing_cf_years:
            raise RuntimeError(
                f"Month {month}: cf_cvg pivot is missing years "
                f"{sorted(missing_cf_years)}."
            )
        cfout = cf.index.to_frame(index=False, name="NUTS_ID")
        for year in YEARS:
            cfout[f"cf_cvg_{year}"] = cf[year].to_numpy()
        cfout["cf_cvg_anomaly_2021_vs_flanks"] = (
            cfout["cf_cvg_2021"]
            - (cfout["cf_cvg_2020"] + cfout["cf_cvg_2022"]) / 2.0
        )
        out = out.merge(cfout, on="NUTS_ID", how="left", validate="one_to_one")

        out["month"] = month
        out["month_name"] = MONTH_NAMES[month - 1]
        out["eligible_primary"] = (
            (out["min_snow_valid_fraction"] >= min_snow)
            & (out["common_support_area_fraction"] >= min_common)
            & out["snow_anomaly_2021_vs_flanks_pp"].notna()
            & out["anomaly_2021_vs_flanks_nw"].notna()
            & out["cf_cvg_anomaly_2021_vs_flanks"].notna()
        )
        rows.append(out)

    return pd.concat(rows, ignore_index=True)


def primary_tests(anomaly: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    tests = []
    models = []

    for month in PRIMARY_MONTHS:
        g = anomaly[
            (anomaly["month"] == month) & anomaly["eligible_primary"]
        ].copy()
        if g.empty:
            diag = anomaly[anomaly["month"] == month]
            raise RuntimeError(
                f"Month {month}: zero eligible regions after primary filters. "
                f"Rows={len(diag)}, "
                f"snow_complete={int(diag['snow_anomaly_2021_vs_flanks_pp'].notna().sum())}, "
                f"viirs_complete={int(diag['anomaly_2021_vs_flanks_nw'].notna().sum())}, "
                f"cf_complete={int(diag['cf_cvg_anomaly_2021_vs_flanks'].notna().sum())}, "
                f"snow_valid={int((diag['min_snow_valid_fraction'] >= 0.5).sum())}, "
                f"viirs_support={int((diag['common_support_area_fraction'] >= 0.5).sum())}."
            )

        rho, p, n = spearman_clean(
            g["snow_anomaly_2021_vs_flanks_pp"],
            g["anomaly_2021_vs_flanks_nw"],
        )

        country = (
            g.groupby("CNTR_CODE", observed=True)
            .agg(
                snow_anomaly_pp=("snow_anomaly_2021_vs_flanks_pp", "median"),
                radiance_anomaly_nw=("anomaly_2021_vs_flanks_nw", "median"),
                latitude=("centroid_latitude", "median"),
                n_nuts3=("NUTS_ID", "nunique"),
            )
            .reset_index()
        )
        crho, cp, cn = spearman_clean(
            country["snow_anomaly_pp"], country["radiance_anomaly_nw"]
        )

        tests.append({
            "month": month,
            "month_name": MONTH_NAMES[month - 1],
            "role": (
                "primary" if month == 2
                else "supportive" if month == 1
                else "negative_control"
            ),
            "n_nuts3": n,
            "rho_region_snow_vs_radiance_anomaly": rho,
            "p_region_snow_vs_radiance_anomaly": p,
            "n_countries": cn,
            "rho_country_median_snow_vs_radiance_anomaly": crho,
            "p_country_median_snow_vs_radiance_anomaly": cp,
            "median_snow_anomaly_pp": float(
                g["snow_anomaly_2021_vs_flanks_pp"].median()
            ),
            "median_radiance_anomaly_nw": float(
                g["anomaly_2021_vs_flanks_nw"].median()
            ),
        })

        model_data = g[[
            "anomaly_2021_vs_flanks_nw",
            "snow_anomaly_2021_vs_flanks_pp",
            "baseline_raw_nw",
            "centroid_latitude",
            "CNTR_CODE",
            "cf_cvg_anomaly_2021_vs_flanks",
        ]].replace([np.inf, -np.inf], np.nan).dropna().copy()
        model_data["log10_baseline"] = np.log10(
            model_data["baseline_raw_nw"].clip(lower=1e-6)
        )
        model_data["log10_baseline_c"] = (
            model_data["log10_baseline"]
            - model_data["log10_baseline"].mean()
        )

        base_formula = (
            "anomaly_2021_vs_flanks_nw ~ log10_baseline_c + "
            "centroid_latitude + cf_cvg_anomaly_2021_vs_flanks + C(CNTR_CODE)"
        )
        full_formula = (
            "anomaly_2021_vs_flanks_nw ~ "
            "snow_anomaly_2021_vs_flanks_pp * log10_baseline_c + "
            "centroid_latitude + cf_cvg_anomaly_2021_vs_flanks + C(CNTR_CODE)"
        )
        base = smf.ols(base_formula, data=model_data).fit(cov_type="HC3")
        full = smf.ols(full_formula, data=model_data).fit(cov_type="HC3")

        snow_name = "snow_anomaly_2021_vs_flanks_pp"
        int_name = "snow_anomaly_2021_vs_flanks_pp:log10_baseline_c"

        try:
            joint = full.wald_test(
                f"{snow_name} = 0, {int_name} = 0", scalar=True
            )
            joint_stat = float(np.asarray(joint.statistic).squeeze())
            joint_p = float(np.asarray(joint.pvalue).squeeze())
        except Exception:
            joint_stat = np.nan
            joint_p = np.nan

        partial_r2 = (
            max(0.0, (base.ssr - full.ssr) / base.ssr)
            if base.ssr > 0 else np.nan
        )
        ci = full.conf_int()

        models.append({
            "month": month,
            "month_name": MONTH_NAMES[month - 1],
            "role": (
                "primary" if month == 2
                else "supportive" if month == 1
                else "negative_control"
            ),
            "n": int(full.nobs),
            "r2_base": float(base.rsquared),
            "r2_full": float(full.rsquared),
            "partial_r2_snow_block": partial_r2,
            "snow_coef_nw_per_pp": float(full.params[snow_name]),
            "snow_ci_low": float(ci.loc[snow_name, 0]),
            "snow_ci_high": float(ci.loc[snow_name, 1]),
            "snow_p": float(full.pvalues[snow_name]),
            "snow_x_brightness_coef": float(full.params[int_name]),
            "snow_x_brightness_p": float(full.pvalues[int_name]),
            "joint_snow_block_wald": joint_stat,
            "joint_snow_block_p": joint_p,
            "cf_cvg_anomaly_coef": float(
                full.params["cf_cvg_anomaly_2021_vs_flanks"]
            ),
            "cf_cvg_anomaly_p": float(
                full.pvalues["cf_cvg_anomaly_2021_vs_flanks"]
            ),
        })

    return pd.DataFrame(tests), pd.DataFrame(models)


def seasonal_region_correlations(
    merged: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    min_snow = float(cfg["min_snow_valid_area_fraction"])
    min_common = float(cfg["min_viirs_common_area_fraction"])
    min_months = int(cfg["min_seasonal_months"])

    rows = []
    for nuts_id, g in merged.groupby("NUTS_ID", observed=True):
        gv = g[
            (g["snow_valid_area_fraction"] >= min_snow)
            & (g["common_support_area_fraction"] >= min_common)
        ].copy()
        rho, p, n = spearman_clean(
            gv["snow_mean_pct"], gv["radiance_common_nw"]
        )
        if n < min_months:
            continue
        first = gv.iloc[0]
        rows.append({
            "NUTS_ID": nuts_id,
            "CNTR_CODE": first["CNTR_CODE"],
            "density_band": first["density_band"],
            "brightness_quintile": first["brightness_quintile"],
            "baseline_raw_nw": first["baseline_raw_nw"],
            "centroid_latitude": first["centroid_latitude"],
            "n_months": n,
            "rho_snow_vs_radiance": rho,
            "p_snow_vs_radiance": p,
        })
    return pd.DataFrame(rows)


def winter_summer_table(merged: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    min_snow = float(cfg["min_snow_valid_area_fraction"])
    min_common = float(cfg["min_viirs_common_area_fraction"])
    region_rows = []
    summary_rows = []

    for year in YEARS:
        g = merged[
            (merged["year"] == year)
            & merged["month"].isin(LEVIN_MONTHS)
        ].copy()
        p = g.pivot(index="NUTS_ID", columns="month", values=[
            "snow_mean_pct", "snow_valid_area_fraction",
            "radiance_common_nw", "common_support_area_fraction",
        ])
        if not {1, 7}.issubset(set(p.columns.get_level_values(1))):
            continue

        meta = g.drop_duplicates("NUTS_ID")[[
            "NUTS_ID", "CNTR_CODE", "centroid_latitude", "baseline_raw_nw"
        ]].set_index("NUTS_ID")

        d = pd.DataFrame(index=p.index)
        d["delta_snow_jan_minus_jul_pp"] = (
            p["snow_mean_pct"][1] - p["snow_mean_pct"][7]
        )
        d["delta_radiance_jan_minus_jul_nw"] = (
            p["radiance_common_nw"][1] - p["radiance_common_nw"][7]
        )
        d["min_snow_valid_fraction"] = pd.concat(
            [p["snow_valid_area_fraction"][1], p["snow_valid_area_fraction"][7]],
            axis=1,
        ).min(axis=1)
        d["min_common_support_fraction"] = pd.concat(
            [p["common_support_area_fraction"][1], p["common_support_area_fraction"][7]],
            axis=1,
        ).min(axis=1)
        d = d.join(meta)
        d = d[
            (d["min_snow_valid_fraction"] >= min_snow)
            & (d["min_common_support_fraction"] >= min_common)
        ].reset_index()
        d["year"] = year
        region_rows.append(d)

        rho, pval, n = spearman_clean(
            d["delta_snow_jan_minus_jul_pp"],
            d["delta_radiance_jan_minus_jul_nw"],
        )
        country = (
            d.groupby("CNTR_CODE", observed=True)
            .agg(
                dsnow=("delta_snow_jan_minus_jul_pp", "median"),
                drad=("delta_radiance_jan_minus_jul_nw", "median"),
            ).reset_index()
        )
        crho, cp, cn = spearman_clean(country["dsnow"], country["drad"])
        summary_rows.append({
            "year": year,
            "n_nuts3": n,
            "rho_region_jan_minus_jul_snow_vs_radiance": rho,
            "p_region": pval,
            "n_countries": cn,
            "rho_country_median_jan_minus_jul_snow_vs_radiance": crho,
            "p_country": cp,
        })

    return pd.concat(region_rows, ignore_index=True), pd.DataFrame(summary_rows)


def lit_area_winter_summer(
    aux: pd.DataFrame,
    snow: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    min_snow = float(cfg["min_snow_valid_area_fraction"])
    merged = aux.merge(
        snow[[
            "NUTS_ID", "year", "month", "snow_mean_pct",
            "snow_valid_area_fraction",
        ]],
        on=["NUTS_ID", "year", "month"],
        how="left",
        validate="many_to_one",
    )

    rows = []
    for year in YEARS:
        g = merged[
            (merged["year"] == year) & merged["month"].isin(LEVIN_MONTHS)
        ].copy()
        snow_p = g.pivot(index="NUTS_ID", columns="month", values=[
            "snow_mean_pct", "snow_valid_area_fraction"
        ])
        meta = g.drop_duplicates("NUTS_ID")[["NUTS_ID"]].set_index("NUTS_ID")

        base = meta.copy()
        base["delta_snow_jan_minus_jul_pp"] = (
            snow_p["snow_mean_pct"][1] - snow_p["snow_mean_pct"][7]
        )
        base["min_snow_valid_fraction"] = pd.concat(
            [
                snow_p["snow_valid_area_fraction"][1],
                snow_p["snow_valid_area_fraction"][7],
            ], axis=1
        ).min(axis=1)

        for threshold in LIT_THRESHOLDS:
            p = g.pivot(
                index="NUTS_ID", columns="month",
                values=f"pct_lit_gt_{threshold}",
            )
            d = base.copy()
            d["delta_lit_pct_jan_minus_jul"] = p[1] - p[7]
            d = d[d["min_snow_valid_fraction"] >= min_snow].dropna()
            rho, pval, n = spearman_clean(
                d["delta_snow_jan_minus_jul_pp"],
                d["delta_lit_pct_jan_minus_jul"],
            )
            rows.append({
                "year": year,
                "threshold_nw": threshold,
                "n_nuts3": n,
                "rho_snow_change_vs_percent_lit_change": rho,
                "p_value": pval,
            })
    return pd.DataFrame(rows)


def make_plots(
    anomaly: pd.DataFrame,
    tests: pd.DataFrame,
    models: pd.DataFrame,
    winter_summer: pd.DataFrame,
    seasonal: pd.DataFrame,
    lit_summary: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    feb = anomaly[
        (anomaly["month"] == 2) & anomaly["eligible_primary"]
    ]
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    ax.scatter(
        feb["snow_anomaly_2021_vs_flanks_pp"],
        feb["anomaly_2021_vs_flanks_nw"],
        s=12,
        alpha=0.45,
    )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel("February 2021 snow anomaly vs 2020/2022 (percentage points)")
    ax.set_ylabel("February 2021 VIIRS anomaly vs 2020/2022 (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Primary test: does anomalous snow track the February 2021 VIIRS pulse?")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "february_snow_vs_viirs_anomaly.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    ax.bar(tests["month_name"], tests["rho_region_snow_vs_radiance_anomaly"])
    ax.axhline(0, linewidth=0.8)
    ax.set_ylabel("Spearman rho")
    ax.set_xlabel("Pre-specified test month")
    ax.set_title("Snow anomaly versus VIIRS anomaly")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(output_dir / "snow_viirs_anomaly_rho_jan_feb_mar.pdf", dpi=dpi)
    plt.close(fig)

    ws21 = winter_summer[winter_summer["year"] == 2021]
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    ax.scatter(
        ws21["delta_snow_jan_minus_jul_pp"],
        ws21["delta_radiance_jan_minus_jul_nw"],
        s=12,
        alpha=0.45,
    )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel("Snow cover: Jan - Jul 2021 (percentage points)")
    ax.set_ylabel("VIIRS radiance: Jan - Jul 2021 (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Levin-style winter-summer contrast, 2021")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "levin_style_jan_july_snow_radiance_2021.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    vals = seasonal["rho_snow_vs_radiance"].dropna().to_numpy(float)
    ax.hist(vals, bins=30)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel("Within-NUTS3 Spearman rho: monthly snow vs VIIRS radiance")
    ax.set_ylabel("NUTS3 regions")
    ax.set_title("Levin-style seasonal snow-radiance relationship, 2020-2022")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(output_dir / "seasonal_snow_radiance_rho_distribution.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    for year, g in lit_summary.groupby("year", observed=True):
        ax.plot(
            g["threshold_nw"],
            g["rho_snow_change_vs_percent_lit_change"],
            marker="o",
            label=str(year),
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xticks(list(LIT_THRESHOLDS))
    ax.set_xticklabels([str(x) for x in LIT_THRESHOLDS])
    ax.set_xlabel("Published VIIRS lit-area threshold (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_ylabel("Spearman rho: Jan-Jul snow change vs percent-lit change")
    ax.set_title("Levin-style percent-lit-area sensitivity")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "levin_style_percent_lit_thresholds.pdf", dpi=dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    integration_self_test()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_snow_viirs_test"]
    raw_dir = PROJECT_ROOT / cfg["mod10cm_raw_dir"]
    if not raw_dir.exists():
        raise FileNotFoundError(raw_dir)

    hdfs = discover_mod10cm(raw_dir)
    cache_tif = configured_path(config, "processed") / "mod10cm_061_geotiff_cache"
    cache_snow = configured_path(config, "processed") / "mod10cm_061_nuts2016_cache"
    cache_viirs = (
        configured_path(config, "processed")
        / "viirs_monthly_literature_metrics_cache"
    )

    nuts = fixed_nuts(config)
    id_field = config["nuts"]["id_field"]
    nuts[id_field] = nuts[id_field].astype(str)

    snow_frames = []
    manifest_rows = []
    for (year, month), hdf in sorted(hdfs.items()):
        tif = ensure_snow_geotiff(hdf, cache_tif, config["nuts"]["crs"])
        print(f"Extracting MOD10CM {year}-{month:02d}", flush=True)
        snow_frames.append(
            extract_snow_month(
                tif, year, month, nuts, config, cache_snow
            )
        )
        manifest_rows.append({
            "year": year,
            "month": month,
            "month_name": MONTH_NAMES[month - 1],
            "hdf_file": hdf.name,
            "cached_geotiff": tif.name,
            "hdf_size_bytes": int(hdf.stat().st_size),
        })

    snow = pd.concat(snow_frames, ignore_index=True)
    # Product invariant: monthly mean snow cover is a percentage.  Any value
    # outside 0..100 means categorical MOD10CM codes have leaked into a mean.
    finite_snow = snow.loc[
        snow["snow_mean_pct"].notna(), "snow_mean_pct"
    ].astype(float)
    if ((finite_snow < 0) | (finite_snow > 100)).any():
        bad = snow.loc[
            snow["snow_mean_pct"].notna()
            & ((snow["snow_mean_pct"] < 0) | (snow["snow_mean_pct"] > 100)),
            ["NUTS_ID", "year", "month", "snow_mean_pct", "snow_valid_area_fraction"],
        ].head(10)
        raise RuntimeError(
            "MOD10CM invariant failed: snow_mean_pct outside 0..100. "
            f"First offending rows:\n{bad.to_string(index=False)}"
        )
    manifest = pd.DataFrame(manifest_rows)

    viirs = load_viirs_monthly_regional(config)
    viirs["NUTS_ID"] = viirs["NUTS_ID"].astype(str)
    vlong = viirs_long(viirs)

    snow["NUTS_ID"] = snow["NUTS_ID"].astype(str)
    merged = vlong.merge(
        snow[[
            "NUTS_ID", "year", "month",
            "snow_mean_pct", "snow_valid_area_fraction",
        ]],
        on=["NUTS_ID", "year", "month"],
        how="left",
        validate="one_to_one",
    )

    bounds = list(cfg["continental_europe_bounds_4326"])
    merged["in_continental_europe_extent"] = in_continental_extent(merged, bounds)
    merged = merged[merged["in_continental_europe_extent"]].copy()
    viirs_europe_ids = set(merged["NUTS_ID"].unique())
    viirs = viirs[viirs["NUTS_ID"].isin(viirs_europe_ids)].copy()

    support_dir = (
        configured_path(config, "processed")
        / "viirs_monthly_2021_localization_support"
    )
    tile = str(cfg["viirs_tile"])
    aux_frames = []
    for year in YEARS:
        for month in AUX_VIIRS_MONTHS:
            avg, cf = discover_viirs_month(config, year, month, tile)
            support = support_dir / f"common_cf_support_month_{month:02d}.tif"
            if not support.exists():
                raise FileNotFoundError(
                    f"Missing {support}. Run future-radiance-monthly-2021-localization."
                )
            print(
                f"Extracting Levin-style VIIRS metrics {year}-{month:02d}",
                flush=True,
            )
            aux_frames.append(
                extract_viirs_literature_metrics(
                    year, month, avg, cf, support, nuts, config, cache_viirs
                )
            )
    aux = pd.concat(aux_frames, ignore_index=True)
    aux["NUTS_ID"] = aux["NUTS_ID"].astype(str)
    aux = aux[aux["NUTS_ID"].isin(viirs_europe_ids)].copy()

    anomaly = build_primary_anomaly(merged, viirs, aux, cfg)
    tests, models = primary_tests(anomaly)
    seasonal = seasonal_region_correlations(merged, cfg)
    winter_summer, winter_summary = winter_summer_table(merged, cfg)
    lit_summary = lit_area_winter_summer(aux, snow, cfg)

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_snow_viirs_test"
    )
    processed.mkdir(parents=True, exist_ok=True)

    outputs = {
        "manifest": processed / "fixed_nuts2016_mod10cm_2020_2022_manifest.csv",
        "snow": processed / "fixed_nuts2016_mod10cm_snow_monthly_2020_2022.csv",
        "merged": processed / "fixed_nuts2016_snow_viirs_monthly_merged.csv",
        "anomaly": processed / "fixed_nuts2016_snow_2021_anomaly_regional.csv",
        "tests": processed / "fixed_nuts2016_snow_2021_anomaly_tests.csv",
        "models": processed / "fixed_nuts2016_snow_2021_anomaly_models.csv",
        "seasonal": processed / "fixed_nuts2016_snow_levin_seasonal_correlations.csv",
        "winter_region": processed / "fixed_nuts2016_snow_levin_jan_july_regional.csv",
        "winter_summary": processed / "fixed_nuts2016_snow_levin_jan_july_summary.csv",
        "lit": processed / "fixed_nuts2016_snow_levin_percent_lit_thresholds.csv",
        "summary": processed / "fixed_nuts2016_snow_viirs_test_summary.txt",
    }

    manifest.to_csv(outputs["manifest"], index=False)
    snow.to_csv(outputs["snow"], index=False)
    merged.to_csv(outputs["merged"], index=False)
    anomaly.to_csv(outputs["anomaly"], index=False)
    tests.to_csv(outputs["tests"], index=False)
    models.to_csv(outputs["models"], index=False)
    seasonal.to_csv(outputs["seasonal"], index=False)
    winter_summer.to_csv(outputs["winter_region"], index=False)
    winter_summary.to_csv(outputs["winter_summary"], index=False)
    lit_summary.to_csv(outputs["lit"], index=False)

    make_plots(
        anomaly, tests, models, winter_summer, seasonal,
        lit_summary, figures, int(cfg["figure_dpi"])
    )

    seasonal_rhos = seasonal["rho_snow_vs_radiance"].dropna()
    feb_test = tests.loc[tests["month"] == 2].iloc[0]
    feb_model = models.loc[models["month"] == 2].iloc[0]
    mar_test = tests.loc[tests["month"] == 3].iloc[0]

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 SNOW / VIIRS 2021 ANOMALY TEST\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            "Literature-led design: mean VIIRS radiance, MOD10CM mean snow cover, "
            "cloud-free coverage, published percent-lit thresholds, Spearman "
            "correlations, and January-vs-July contrasts following the analytical "
            "logic of Levin & Zhang (2017). NDVI/albedo are not included in this "
            "first snow-specific test because they have not yet been acquired.\n"
        )
        h.write(
            "Primary anomaly: 2021 - mean(2020, 2022). February is primary; "
            "January supportive; March negative control.\n"
        )
        h.write(
            f"MOD10CM valid values: 0-100%. Primary snow valid-area threshold: "
            f"{100*float(cfg['min_snow_valid_area_fraction']):.0f}%.\n"
        )
        h.write(
            f"Continental-Europe centroid extent: {bounds}. "
            f"Regions represented: {merged['NUTS_ID'].nunique()}.\n"
        )
        h.write(
            "October 2022 MOD10CM is structurally absent and is not imputed.\n\n"
        )

        h.write("PRIMARY / SUPPORTIVE / NEGATIVE-CONTROL SPEARMAN TESTS\n")
        h.write(tests.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nCOUNTRY/LATITUDE/BRIGHTNESS/CF_CVG-ADJUSTED MODELS\n")
        h.write(models.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nLEVIN-STYLE JANUARY-MINUS-JULY TESTS\n")
        h.write(winter_summary.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nLEVIN-STYLE PERCENT-LIT THRESHOLD TESTS\n")
        h.write(lit_summary.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nSEASONAL WITHIN-REGION SUMMARY\n")
        h.write(
            f"N regions with >= {int(cfg['min_seasonal_months'])} paired months: "
            f"{len(seasonal)}\n"
        )
        h.write(
            f"Median within-region snow-radiance Spearman rho: "
            f"{seasonal_rhos.median():+.4f}\n"
        )
        h.write(
            f"Fraction of regions with positive rho: "
            f"{100*(seasonal_rhos > 0).mean():.1f}%\n"
        )

        h.write("\nPRE-SPECIFIED INTERPRETATION CHECK\n")
        h.write(
            f"- February region-level rho: "
            f"{feb_test['rho_region_snow_vs_radiance_anomaly']:+.4f}; "
            f"country-median rho: "
            f"{feb_test['rho_country_median_snow_vs_radiance_anomaly']:+.4f}.\n"
        )
        h.write(
            f"- February adjusted snow coefficient: "
            f"{feb_model['snow_coef_nw_per_pp']:+.6f} nW per snow pp "
            f"(p={feb_model['snow_p']:.4g}); "
            f"snow-block partial R2={feb_model['partial_r2_snow_block']:.4f}; "
            f"joint p={feb_model['joint_snow_block_p']:.4g}.\n"
        )
        h.write(
            f"- February snow x baseline-brightness interaction p="
            f"{feb_model['snow_x_brightness_p']:.4g}.\n"
        )
        h.write(
            f"- March negative-control rho: "
            f"{mar_test['rho_region_snow_vs_radiance_anomaly']:+.4f}.\n"
        )
        h.write(
            "- Support for snow/albedo requires a positive February association "
            "that persists after country, latitude, baseline brightness and "
            "cloud-free coverage controls, with March materially weaker/absent. "
            "A positive snow x brightness interaction would strengthen the "
            "reflection mechanism.\n"
        )
        h.write(
            "- Concordant January-vs-July radiance and percent-lit responses "
            "provide literature-method replication; disagreement does not "
            "invalidate the anomaly test but would weaken generalization.\n"
        )
        h.write(
            "- This stage tests snow specifically. If supported, NDVI and surface "
            "albedo should be added next to distinguish snow presence from the "
            "broader surface-reflectance mechanism.\n"
        )

    print("\nFIXED NUTS2016 SNOW / VIIRS 2021 ANOMALY TEST")
    print("\nPRIMARY TESTS")
    print(tests.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nADJUSTED MODELS")
    print(models.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nLEVIN-STYLE JAN-JUL SUMMARY")
    print(winter_summary.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
