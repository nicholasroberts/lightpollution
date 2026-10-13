#!/usr/bin/env python3
"""Locate the 2021 VIIRS disturbance in monthly VCMSLCFG composites.

For each calendar month, the 2020, 2021 and 2022 VCMSLCFG composites are
compared on exactly the same pixel observation support: pixels with cf_cvg > 0
in all three years for that month. The primary regional anomaly is

    R_2021 - (R_2020 + R_2022) / 2

where each R is the physical-area-weighted mean radiance on that month-specific
common support.

Actual start/end dates are retained because some EOG monthly composites do not
span an identical number of days. This stage is diagnostic only and does not
modify annual production data or apply any correction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
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
from shapely.geometry import box

from diagnose_fixed_nuts2016_2021_reversal import add_latitude, load_locked_panel
from diagnose_fixed_nuts2016_mask_occupancy import attach_fixed_brightness
from project_config import PROJECT_ROOT, configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import load_fixed_nuts


YEARS = (2020, 2021, 2022)
MONTH_NAMES = ("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
AVG_RE = re.compile(
    r"^SVDNB_(?P<satellite>npp|j01)_(?P<start>\d{8})-(?P<end>\d{8})_"
    r"(?P<tile>[^_]+)_vcmslcfg_v10_.*\.avg_rade9h\.tif$",
    re.IGNORECASE,
)


def discover_files(cfg: dict) -> dict[int, dict[int, dict[str, Path]]]:
    tile = str(cfg["tile"])
    out: dict[int, dict[int, dict[str, Path]]] = {}
    for year in YEARS:
        root = PROJECT_ROOT / str(cfg["raw_dir_template"]).format(year=year)
        if not root.exists():
            raise FileNotFoundError(f"Missing monthly directory: {root}")
        months: dict[int, dict[str, Path]] = {}
        for avg in sorted(root.glob("*.avg_rade9h.tif")):
            m = AVG_RE.match(avg.name)
            if not m:
                continue
            if int(m.group("start")[:4]) != year or m.group("tile") != tile:
                continue
            month = int(m.group("start")[4:6])
            if month in months:
                raise RuntimeError(f"Duplicate {year}-{month:02d} radiance raster.")
            cf = Path(str(avg).replace(".avg_rade9h.tif", ".cf_cvg.tif"))
            cvg = Path(str(avg).replace(".avg_rade9h.tif", ".cvg.tif"))
            for p in (avg, cf, cvg):
                if not p.exists():
                    raise FileNotFoundError(f"Missing monthly companion raster: {p}")
            months[month] = {"avg": avg, "cf": cf, "cvg": cvg}
        missing = sorted(set(range(1, 13)) - set(months))
        if missing:
            raise RuntimeError(f"{year}: missing calendar months {missing}.")
        out[year] = months
    return out


def parse_dates(path: Path) -> tuple[str, str, int]:
    m = AVG_RE.match(path.name)
    if not m:
        raise RuntimeError(f"Cannot parse monthly filename: {path.name}")
    start = pd.Timestamp(m.group("start"))
    end = pd.Timestamp(m.group("end"))
    return start.date().isoformat(), end.date().isoformat(), int((end - start).days + 1)


def parse_satellite(path: Path) -> str:
    m = AVG_RE.match(path.name)
    if not m:
        raise RuntimeError(f"Cannot parse monthly filename: {path.name}")
    return m.group("satellite").lower()


def validate_alignment(files, expected_crs: str) -> pd.DataFrame:
    reference = None
    rows = []
    for year in YEARS:
        for month in range(1, 13):
            for layer, path in files[year][month].items():
                with rasterio.open(path) as ds:
                    if ds.count != 1:
                        raise RuntimeError(f"Expected one band: {path}")
                    if ds.crs is None or ds.crs.to_string() != expected_crs:
                        raise RuntimeError(
                            f"Unexpected CRS for {path.name}: {ds.crs}; expected {expected_crs}"
                        )
                    signature = (
                        ds.width, ds.height, tuple(ds.transform),
                        tuple(ds.bounds), ds.crs.to_string()
                    )
                    if reference is None:
                        reference = signature
                    elif signature != reference:
                        raise RuntimeError(f"Monthly rasters are not aligned: {path.name}")
                    rows.append({
                        "year": year,
                        "month": month,
                        "layer": layer,
                        "filename": path.name,
                        "size_bytes": int(path.stat().st_size),
                        "width": ds.width,
                        "height": ds.height,
                        "crs": ds.crs.to_string(),
                        "nodata": ds.nodata,
                    })
    return pd.DataFrame(rows)


def month_signature(paths: dict[int, dict[str, Path]]) -> str:
    text = "|".join(
        f"{year}:{layer}:{paths[year][layer].name}:{paths[year][layer].stat().st_size}"
        for year in YEARS for layer in ("avg", "cf", "cvg")
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def build_common_support_mask(
    month: int,
    paths: dict[int, dict[str, Path]],
    cache_dir: Path,
    expected_crs: str,
) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    mask_path = cache_dir / f"common_cf_support_month_{month:02d}.tif"
    meta_path = cache_dir / f"common_cf_support_month_{month:02d}.json"
    signature = month_signature(paths)

    if mask_path.exists() and meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            meta = {}
        if meta.get("source_signature") == signature:
            with rasterio.open(mask_path) as ds:
                if (
                    ds.count == 1
                    and ds.crs is not None
                    and ds.crs.to_string() == expected_crs
                    and ds.dtypes[0] == "uint8"
                ):
                    print(f"Using cached support mask for month {month:02d}", flush=True)
                    return mask_path

    sources = {year: rasterio.open(paths[year]["cf"]) for year in YEARS}
    try:
        ref = sources[2020]
        profile = ref.profile.copy()
        profile.update(
            dtype="uint8", count=1, nodata=None,
            compress="deflate", predictor=1, BIGTIFF="IF_SAFER"
        )
        tmp = mask_path.with_suffix(".tif.part")
        if tmp.exists():
            tmp.unlink()
        print(f"Building common 2020/2021/2022 support for month {month:02d}", flush=True)
        with rasterio.open(tmp, "w", **profile) as dst:
            for _, window in ref.block_windows(1):
                support = np.ones(
                    (int(window.height), int(window.width)), dtype=bool
                )
                for year in YEARS:
                    arr = sources[year].read(1, window=window, masked=True)
                    data = np.asarray(arr.data if np.ma.isMaskedArray(arr) else arr)
                    valid = (
                        ~np.ma.getmaskarray(arr)
                        if np.ma.isMaskedArray(arr)
                        else np.ones(arr.shape, dtype=bool)
                    )
                    support &= valid & np.isfinite(data) & (data > 0)
                dst.write(support.astype("uint8"), 1, window=window)
        tmp.replace(mask_path)
        meta_path.write_text(
            json.dumps({"month": month, "source_signature": signature}, indent=2),
            encoding="utf-8",
        )
    finally:
        for ds in sources.values():
            ds.close()
    return mask_path


class CommonSupportRadianceSource(RasterioRasterSource):
    """Expose radiance only where a precomputed common-support mask equals 1."""

    def __init__(self, avg_ds, support_ds, *, name: str):
        super().__init__(avg_ds, 1, name=name)
        self.avg = RasterioRasterSource(avg_ds, 1, name=f"{name}_avg")
        self.support = RasterioRasterSource(support_ds, 1, name=f"{name}_support")

    def read_window(self, x0, y0, nx, ny):
        avg = self.avg.read_window(x0, y0, nx, ny)
        sup = self.support.read_window(x0, y0, nx, ny)

        avg_data = np.asarray(avg.data if np.ma.isMaskedArray(avg) else avg, dtype=np.float32)
        sup_data = np.asarray(sup.data if np.ma.isMaskedArray(sup) else sup)

        avg_valid = (
            ~np.ma.getmaskarray(avg)
            if np.ma.isMaskedArray(avg)
            else np.ones(avg_data.shape, dtype=bool)
        )
        sup_valid = (
            ~np.ma.getmaskarray(sup)
            if np.ma.isMaskedArray(sup)
            else np.ones(sup_data.shape, dtype=bool)
        )
        valid = avg_valid & sup_valid & np.isfinite(avg_data) & (sup_data == 1)
        return np.ma.masked_array(avg_data, mask=~valid, copy=False)


def integration_self_test() -> None:
    support_mem = MemoryFile()
    avg_mems = []
    datasets = []
    try:
        support_ds = support_mem.open(
            driver="GTiff", height=2, width=2, count=1, dtype="uint8",
            crs="EPSG:4326", transform=from_origin(0, 2, 1, 1)
        )
        support_ds.write(np.array([[1, 0], [1, 0]], dtype="uint8"), 1)

        for value in (1.0, 1.5, 2.0):
            mem = MemoryFile()
            ds = mem.open(
                driver="GTiff", height=2, width=2, count=1, dtype="float32",
                crs="EPSG:4326", transform=from_origin(0, 2, 1, 1)
            )
            ds.write(np.full((2, 2), value, dtype="float32"), 1)
            avg_mems.append(mem)
            datasets.append(ds)

        gdf = gpd.GeoDataFrame(
            {"NUTS_ID": ["TEST"]},
            geometry=[box(0, 0, 2, 2)],
            crs="EPSG:4326",
        )
        sources = [RasterioRasterSource(support_ds, 1, name="common_support")]
        for year, ds in zip(YEARS, datasets):
            sources.append(
                CommonSupportRadianceSource(ds, support_ds, name=f"r{year}")
            )
        out = exact_extract(
            sources, gdf, ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=["NUTS_ID"], output="pandas", progress=False
        ).iloc[0]
        if not np.isclose(out["common_support_mean"], 0.5, atol=0.01):
            raise RuntimeError("Self-test common-support fraction failed.")
        if not np.isclose(out["r2020_mean"], 1.0, atol=1e-6):
            raise RuntimeError("Self-test 2020 common radiance failed.")
        if not np.isclose(out["r2021_mean"], 1.5, atol=1e-6):
            raise RuntimeError("Self-test 2021 common radiance failed.")
        if not np.isclose(out["r2022_mean"], 2.0, atol=1e-6):
            raise RuntimeError("Self-test 2022 common radiance failed.")
    finally:
        for ds in datasets:
            ds.close()
        for mem in avg_mems:
            mem.close()
        if "support_ds" in locals():
            support_ds.close()
        support_mem.close()
    print("Monthly 2021 localization integration self-test: PASS", flush=True)


def load_geometry_subset(config: dict, bounds, locked_ids: set[str]):
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release).copy()
    id_field = config["nuts"]["id_field"]
    nuts[id_field] = nuts[id_field].astype(str)
    nuts = nuts[nuts[id_field].isin(locked_ids)].copy()
    tile = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    fully_covered = nuts.geometry.apply(tile.covers)
    excluded = nuts.loc[~fully_covered, id_field].tolist()
    nuts = nuts.loc[fully_covered].copy()
    return nuts, excluded


def extract_month(month, paths, support_path, nuts, config, cache_dir):
    id_field = config["nuts"]["id_field"]
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"nuts2016_month_{month:02d}_common_support.parquet"
    signature = month_signature(paths)

    if cache.exists():
        d = pd.read_parquet(cache)
        needed = {
            id_field, "month", "source_signature",
            "common_support_area_fraction",
            "radiance_2020_common_nw", "radiance_2021_common_nw",
            "radiance_2022_common_nw",
        }
        if needed.issubset(d.columns):
            if (
                d["source_signature"].dropna().astype(str).unique().tolist() == [signature]
                and set(d[id_field].astype(str)) == set(nuts[id_field].astype(str))
            ):
                print(f"Using cached extraction for month {month:02d}", flush=True)
                return d

    avg_ds = {}
    support_ds = rasterio.open(support_path)
    try:
        for year in YEARS:
            avg_ds[year] = rasterio.open(paths[year]["avg"])
        sources = [RasterioRasterSource(support_ds, 1, name="common_support")]
        for year in YEARS:
            sources.append(
                CommonSupportRadianceSource(
                    avg_ds[year], support_ds, name=f"r{year}"
                )
            )
        print(f"Extracting month {month:02d} common-support radiance", flush=True)
        out = exact_extract(
            sources, nuts, ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=[id_field], output="pandas", progress=True
        )
    finally:
        for ds in avg_ds.values():
            ds.close()
        support_ds.close()

    out = out.rename(columns={
        "common_support_mean": "common_support_area_fraction",
        "r2020_mean": "radiance_2020_common_nw",
        "r2021_mean": "radiance_2021_common_nw",
        "r2022_mean": "radiance_2022_common_nw",
    })
    required = {
        id_field, "common_support_area_fraction",
        "radiance_2020_common_nw", "radiance_2021_common_nw",
        "radiance_2022_common_nw",
    }
    missing = required - set(out.columns)
    if missing:
        raise RuntimeError(
            f"Unexpected exactextract output for month {month:02d}; "
            f"missing {sorted(missing)}; found {list(out.columns)}"
        )
    d = out[[id_field] + sorted(required - {id_field})].copy()
    d.insert(0, "month", int(month))
    d["source_signature"] = signature
    d.to_parquet(cache, index=False)
    return d


def date_table(files) -> pd.DataFrame:
    rows = []
    for month in range(1, 13):
        row = {"month": month, "month_name": MONTH_NAMES[month - 1]}
        for year in YEARS:
            start, end, days = parse_dates(files[year][month]["avg"])
            row[f"start_date_{year}"] = start
            row[f"end_date_{year}"] = end
            row[f"n_calendar_days_{year}"] = days
            row[f"satellite_{year}"] = parse_satellite(files[year][month]["avg"])
        row["same_satellite_all_years"] = len(
            {row[f"satellite_{y}"] for y in YEARS}
        ) == 1
        row["day_count_range"] = (
            max(row[f"n_calendar_days_{y}"] for y in YEARS)
            - min(row[f"n_calendar_days_{y}"] for y in YEARS)
        )
        rows.append(row)
    return pd.DataFrame(rows)


def add_metrics(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    p = pd.to_numeric(d["common_support_area_fraction"], errors="coerce")
    if ((p < -1e-8) | (p > 1 + 1e-8)).any():
        raise RuntimeError("Invalid common-support area fraction.")
    d["common_support_area_fraction"] = p.clip(0, 1)
    d["common_support_area_pct"] = 100 * d["common_support_area_fraction"]
    d["delta_2020_2021_nw"] = d["radiance_2021_common_nw"] - d["radiance_2020_common_nw"]
    d["delta_2021_2022_nw"] = d["radiance_2022_common_nw"] - d["radiance_2021_common_nw"]
    d["anomaly_2021_vs_flanks_nw"] = (
        d["radiance_2021_common_nw"]
        - (d["radiance_2020_common_nw"] + d["radiance_2022_common_nw"]) / 2
    )
    return d


def spearman_clean(x: pd.Series, y: pd.Series):
    z = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(z) < 3 or z["x"].nunique() < 2 or z["y"].nunique() < 2:
        return np.nan, np.nan, len(z)
    r = spearmanr(z["x"], z["y"])
    return float(r.statistic), float(r.pvalue), len(z)


def country_summary(d, min_common, min_country_regions):
    rows = []
    for (month, country), g in d.groupby(["month", "CNTR_CODE"], observed=True):
        gv = g[
            (g["common_support_area_fraction"] >= min_common)
            & g["anomaly_2021_vs_flanks_nw"].notna()
        ]
        if len(gv) < min_country_regions:
            continue
        rows.append({
            "month": int(month),
            "month_name": MONTH_NAMES[int(month) - 1],
            "CNTR_CODE": str(country),
            "n_nuts3": int(len(gv)),
            "median_centroid_latitude": float(gv["centroid_latitude"].median()),
            "median_common_support_area_pct": float(gv["common_support_area_pct"].median()),
            "median_radiance_2020_common_nw": float(gv["radiance_2020_common_nw"].median()),
            "median_radiance_2021_common_nw": float(gv["radiance_2021_common_nw"].median()),
            "median_radiance_2022_common_nw": float(gv["radiance_2022_common_nw"].median()),
            "median_anomaly_2021_vs_flanks_nw": float(
                gv["anomaly_2021_vs_flanks_nw"].median()
            ),
        })
    return pd.DataFrame(rows)


def month_summary(d, countries, dates, min_common):
    rows = []
    for month, g in d.groupby("month", observed=True):
        gv = g[
            (g["common_support_area_fraction"] >= min_common)
            & g["anomaly_2021_vs_flanks_nw"].notna()
        ]
        c = countries[countries["month"] == int(month)]
        rho, p, n = spearman_clean(
            gv["centroid_latitude"], gv["anomaly_2021_vs_flanks_nw"]
        )
        rows.append({
            "month": int(month),
            "month_name": MONTH_NAMES[int(month) - 1],
            "n_valid_nuts3": int(len(gv)),
            "n_valid_countries": int(c["CNTR_CODE"].nunique()) if len(c) else 0,
            "median_common_support_area_pct": float(gv["common_support_area_pct"].median())
            if len(gv) else np.nan,
            "regional_median_anomaly_2021_vs_flanks_nw": float(
                gv["anomaly_2021_vs_flanks_nw"].median()
            ) if len(gv) else np.nan,
            "pct_regions_positive_anomaly": float(
                100 * (gv["anomaly_2021_vs_flanks_nw"] > 0).mean()
            ) if len(gv) else np.nan,
            "country_balanced_radiance_2020_common_nw": float(
                c["median_radiance_2020_common_nw"].mean()
            ) if len(c) else np.nan,
            "country_balanced_radiance_2021_common_nw": float(
                c["median_radiance_2021_common_nw"].mean()
            ) if len(c) else np.nan,
            "country_balanced_radiance_2022_common_nw": float(
                c["median_radiance_2022_common_nw"].mean()
            ) if len(c) else np.nan,
            "country_balanced_anomaly_2021_vs_flanks_nw": float(
                c["median_anomaly_2021_vs_flanks_nw"].mean()
            ) if len(c) else np.nan,
            "rho_latitude_vs_anomaly": rho,
            "p_latitude_vs_anomaly": p,
            "n_latitude_vs_anomaly": n,
        })
    out = pd.DataFrame(rows).merge(
        dates, on=["month", "month_name"], how="left", validate="one_to_one"
    )
    out["abs_country_balanced_anomaly_rank"] = np.nan
    primary = out["same_satellite_all_years"].fillna(False)
    out.loc[primary, "abs_country_balanced_anomaly_rank"] = (
        out.loc[primary, "country_balanced_anomaly_2021_vs_flanks_nw"]
        .abs().rank(method="min", ascending=False)
    )
    out["primary_same_satellite_comparison"] = primary
    return out.sort_values("month").reset_index(drop=True)


def density_summary(d, min_common):
    rows = []
    for (month, band), g in d.groupby(["month", "density_band"], observed=True):
        gv = g[
            (g["common_support_area_fraction"] >= min_common)
            & g["anomaly_2021_vs_flanks_nw"].notna()
        ]
        if gv.empty:
            continue
        rows.append({
            "month": int(month),
            "month_name": MONTH_NAMES[int(month) - 1],
            "density_band": str(band),
            "n_nuts3": int(len(gv)),
            "median_common_support_area_pct": float(gv["common_support_area_pct"].median()),
            "median_anomaly_2021_vs_flanks_nw": float(
                gv["anomaly_2021_vs_flanks_nw"].median()
            ),
        })
    return pd.DataFrame(rows)


def make_plots(monthly, countries, density, output_dir: Path, dpi: int):
    output_dir.mkdir(parents=True, exist_ok=True)
    x = monthly["month"]

    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    ax.plot(x, monthly["country_balanced_anomaly_2021_vs_flanks_nw"], marker="o")
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(range(1, 13), MONTH_NAMES)
    ax.set_xlabel("Calendar month")
    ax.set_ylabel("2021 anomaly vs mean of 2020/2022 (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Monthly localization of the 2021 VCMSLCFG disturbance")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "monthly_2021_anomaly_vs_flanks.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    for year in YEARS:
        ax.plot(
            x, monthly[f"country_balanced_radiance_{year}_common_nw"],
            marker="o", label=str(year)
        )
    ax.set_xticks(range(1, 13), MONTH_NAMES)
    ax.set_xlabel("Calendar month")
    ax.set_ylabel("Country-balanced radiance on common support (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Matched monthly VCMSLCFG radiance")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "monthly_matched_radiance_2020_2022.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    ax.plot(x, monthly["rho_latitude_vs_anomaly"], marker="o")
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(range(1, 13), MONTH_NAMES)
    ax.set_xlabel("Calendar month")
    ax.set_ylabel("Spearman rho: latitude vs 2021 anomaly")
    ax.set_title("Latitude structure of the monthly 2021 anomaly")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "monthly_anomaly_latitude_rho.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.4, 5.9))
    for band in ["<30", "30-100", "100-300", "300-1000", ">=1000"]:
        g = density[density["density_band"] == band].sort_values("month")
        if not g.empty:
            ax.plot(g["month"], g["median_anomaly_2021_vs_flanks_nw"], marker="o", label=band)
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(range(1, 13), MONTH_NAMES)
    ax.set_xlabel("Calendar month")
    ax.set_ylabel("Median 2021 anomaly vs flanks (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Monthly 2021 anomaly by fixed population-density band")
    ax.legend(title="people km$^{-2}$")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "monthly_2021_anomaly_by_density.pdf", dpi=dpi)
    plt.close(fig)

    if not countries.empty:
        order = (
            countries.groupby("CNTR_CODE", observed=True)["median_centroid_latitude"]
            .median().sort_values().index.tolist()
        )
        pivot = countries.pivot(
            index="CNTR_CODE", columns="month",
            values="median_anomaly_2021_vs_flanks_nw"
        ).reindex(order)
        arr = pivot.to_numpy(float)
        finite = np.abs(arr[np.isfinite(arr)])
        vmax = float(np.quantile(finite, 0.95)) if finite.size else 1.0
        vmax = max(vmax, 1e-6)
        fig, ax = plt.subplots(figsize=(10.0, 8.5))
        im = ax.imshow(
            arr, aspect="auto", interpolation="nearest",
            vmin=-vmax, vmax=vmax, cmap="coolwarm"
        )
        ax.set_xticks(np.arange(12), MONTH_NAMES)
        ax.set_yticks(np.arange(len(order)), order)
        ax.set_xlabel("Calendar month")
        ax.set_ylabel("Country, ordered south to north")
        ax.set_title("Country-median 2021 anomaly versus 2020/2022 flanks")
        fig.colorbar(im, ax=ax, label="nW cm$^{-2}$ sr$^{-1}$")
        fig.tight_layout()
        fig.savefig(output_dir / "monthly_country_anomaly_heatmap.pdf", dpi=dpi)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    integration_self_test()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_monthly_2021_localization"]
    min_common = float(cfg["min_common_area_fraction"])
    min_country_regions = int(cfg["min_country_regions"])

    files = discover_files(cfg)
    manifest = validate_alignment(files, config["nuts"]["crs"])
    dates = date_table(files)

    panel = add_latitude(attach_fixed_brightness(load_locked_panel(config)), config)
    id_field = config["nuts"]["id_field"]
    panel[id_field] = panel[id_field].astype(str)
    locked_ids = set(panel[id_field])

    with rasterio.open(files[2020][1]["avg"]) as ref:
        bounds = ref.bounds
    nuts, excluded = load_geometry_subset(config, bounds, locked_ids)

    print(
        f"Monthly tile fully contains {len(nuts)} of {len(locked_ids)} locked "
        f"NUTS3 regions; {len(excluded)} excluded at tile edge/outside.",
        flush=True,
    )

    processed = configured_path(config, "processed")
    support_cache = processed / "viirs_monthly_2021_localization_support"
    extract_cache = processed / "viirs_monthly_2021_localization_cache"

    frames = []
    for month in range(1, 13):
        paths = {year: files[year][month] for year in YEARS}
        support = build_common_support_mask(
            month, paths, support_cache, config["nuts"]["crs"]
        )
        frames.append(
            extract_month(month, paths, support, nuts, config, extract_cache)
        )

    extracted = pd.concat(frames, ignore_index=True)
    extracted[id_field] = extracted[id_field].astype(str)

    meta = panel[panel["year"] == 2020][[
        id_field, "CNTR_CODE", "density_band", "brightness_quintile",
        "baseline_raw_nw", "centroid_latitude", "centroid_longitude"
    ]].copy()
    regional = extracted.merge(meta, on=id_field, how="left", validate="many_to_one")
    if regional["CNTR_CODE"].isna().any():
        raise RuntimeError("Monthly extraction lost locked-panel metadata.")
    regional = add_metrics(regional)

    labels = list(config["analysis"]["density_temporal_change"]["density_band_labels"])
    regional["density_band"] = pd.Categorical(
        regional["density_band"], categories=labels, ordered=True
    )

    countries = country_summary(regional, min_common, min_country_regions)
    monthly = month_summary(regional, countries, dates, min_common)
    density = density_summary(regional, min_common)

    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_monthly_2021_localization"
    )
    processed.mkdir(parents=True, exist_ok=True)

    outputs = {
        "manifest": processed / "fixed_nuts2016_monthly_2021_localization_manifest.csv",
        "regional": processed / "fixed_nuts2016_monthly_2021_localization_regional.csv",
        "monthly": processed / "fixed_nuts2016_monthly_2021_localization_summary.csv",
        "country": processed / "fixed_nuts2016_monthly_2021_localization_by_country.csv",
        "density": processed / "fixed_nuts2016_monthly_2021_localization_by_density.csv",
        "summary": processed / "fixed_nuts2016_monthly_2021_localization_summary.txt",
    }

    manifest["tile_fully_covered_nuts3"] = len(nuts)
    manifest["tile_excluded_locked_nuts3"] = len(excluded)
    manifest.to_csv(outputs["manifest"], index=False)
    regional.to_csv(outputs["regional"], index=False)
    monthly.to_csv(outputs["monthly"], index=False)
    countries.to_csv(outputs["country"], index=False)
    density.to_csv(outputs["density"], index=False)

    make_plots(monthly, countries, density, figures, int(cfg["figure_dpi"]))

    ranked = monthly[
        monthly["primary_same_satellite_comparison"]
    ].sort_values("abs_country_balanced_anomaly_rank")
    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 MONTHLY LOCALIZATION OF THE 2021 VCMSLCFG ANOMALY\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            "Each calendar month uses pixels with cf_cvg > 0 in all of 2020, "
            "2021 and 2022 for that same month. The primary anomaly is "
            "2021 - mean(2020, 2022) on this fixed monthly support.\n"
        )
        h.write(
            f"Common-support summaries require >= {100*min_common:.1f}% of "
            f"NUTS3 area; country summaries require >= {min_country_regions} "
            "valid NUTS3 regions.\n"
        )
        h.write(
            f"The tile fully contains {len(nuts)} of {len(locked_ids)} locked "
            f"NUTS3 regions; {len(excluded)} are excluded at the tile edge/outside.\n"
        )
        h.write(
            "Actual monthly start/end dates, day counts and source satellite "
            "are retained because some EOG composites do not span the full "
            "calendar month and EOG may substitute NOAA-20/J01.\n"
        )
        h.write(
            "Primary monthly anomaly ranking is restricted to months using the "
            "same satellite in 2020, 2021 and 2022. Mixed-satellite months are "
            "retained only as explicitly flagged sensitivity results.\n\n"
        )
        h.write("MONTHLY SUMMARY\n")
        h.write(monthly.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nMONTHS RANKED BY ABSOLUTE COUNTRY-BALANCED 2021 ANOMALY\n")
        h.write(
            ranked[[
                "month", "month_name",
                "country_balanced_anomaly_2021_vs_flanks_nw",
                "abs_country_balanced_anomaly_rank",
                "n_valid_nuts3", "n_valid_countries",
                "rho_latitude_vs_anomaly", "day_count_range"
            ]].to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )
        h.write("\n\nDECISION RULE\n")
        h.write(
            "- One/few months with a large positive 2021 anomaly and strong "
            "latitude structure indicate a seasonal/product-specific disturbance.\n"
        )
        h.write(
            "- A step beginning in January/February and persisting across later "
            "months is more consistent with a broad calibration/LUT change.\n"
        )
        h.write(
            "- A summer-confined anomaly that strengthens markedly with latitude "
            "supports a stray-light/solar-illumination pathway.\n"
        )
        h.write(
            "- A diffuse year-round anomaly with weak seasonality/latitude points "
            "toward broader compositing or calibration behaviour.\n"
        )
        h.write(
            "- Different monthly date ranges remain an explicit sampling caveat. "
            "Mixed-satellite months (notably August 2022 J01/NOAA-20) are excluded "
            "from the primary anomaly ranking because a platform difference could "
            "mimic a temporal anomaly.\n"
        )
        h.write(
            "- This localization stage does not define an annual correction.\n"
        )

    print("\nFIXED NUTS2016 MONTHLY 2021 LOCALIZATION")
    print("\nMONTHLY SUMMARY")
    print(monthly.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nMONTHS RANKED BY ABSOLUTE COUNTRY-BALANCED ANOMALY")
    print(
        ranked[[
            "month", "month_name",
            "country_balanced_anomaly_2021_vs_flanks_nw",
            "abs_country_balanced_anomaly_rank",
            "n_valid_nuts3", "n_valid_countries",
            "rho_latitude_vs_anomaly", "day_count_range"
        ]].to_string(index=False, float_format=lambda x: f"{x:.6g}")
    )
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
