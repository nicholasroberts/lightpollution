#!/usr/bin/env python3
"""July VCMCFG versus VCMSLCFG pilot for the 2021 VIIRS anomaly.

Purpose
-------
Test the specific hypothesis that the strong northern 2021 disturbance seen in
Annual VNL is associated with the monthly stray-light-corrected pathway.

The pilot deliberately uses only July 2020, July 2021 and July 2022 and pairs
VCMCFG with VCMSLCFG *within the same year/date range*. This is important
because the downloaded 2022 product ends on 26 July whereas 2020 and 2021 run
through 31 July. Absolute cross-year July means are therefore descriptive only;
the primary diagnostics are within-year paired product contrasts.

For every NUTS3 and year the script measures:
- area with cloud-free observations in VCMCFG and VCMSLCFG (cf_cvg > 0);
- common observed area and VCMSLCFG-only added area;
- radiance on the exact common observed pixel support;
- VCMSLCFG-minus-VCMCFG radiance difference on that common support;
- the polygon-scale radiance contribution added by changed support versus
  changed radiance on common support; and
- area-weighted cf_cvg and cvg counts for sampling diagnostics.

Zeros in avg_rade9h are not treated as missing. Observation support is defined
explicitly from cf_cvg > 0 because these monthly rasters have nodata=None.

This is a diagnostic stage only. It does not modify annual production data or
apply a correction.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

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
import geopandas as gpd

from diagnose_fixed_nuts2016_2021_reversal import (
    add_latitude,
    load_locked_panel,
)
from diagnose_fixed_nuts2016_mask_occupancy import attach_fixed_brightness
from project_config import PROJECT_ROOT, configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import load_fixed_nuts


AVG_RE = re.compile(
    r"^SVDNB_npp_(?P<start>\d{8})-(?P<end>\d{8})_75N060W_"
    r"(?P<product>vcmcfg|vcmslcfg)_v10_.*\.avg_rade9h\.tif$"
)

DERIVED_MODES = [
    "vcm_observed_area",
    "vcmsl_observed_area",
    "common_observed_area",
    "vcm_observed_radiance",
    "vcmsl_observed_radiance",
    "vcm_common_radiance",
    "vcmsl_common_radiance",
]


def _as_data_valid(arr):
    if np.ma.isMaskedArray(arr):
        data = np.asarray(arr.data, dtype=np.float32)
        valid = ~np.ma.getmaskarray(arr)
    else:
        data = np.asarray(arr, dtype=np.float32)
        valid = np.ones(data.shape, dtype=bool)
    valid &= np.isfinite(data)
    return data, valid


class MonthlyPairSource(RasterioRasterSource):
    """Derived source from aligned VCMCFG/VCMSLCFG radiance and cf_cvg."""

    def __init__(self, datasets: dict[str, object], *, mode: str, name: str):
        if mode not in DERIVED_MODES:
            raise ValueError(f"Unknown mode: {mode}")
        super().__init__(datasets["vcmsl_avg"], 1, name=name)
        self.mode = mode
        self.sources = {
            key: RasterioRasterSource(ds, 1, name=key)
            for key, ds in datasets.items()
            if key in {"vcm_avg", "vcm_cf", "vcmsl_avg", "vcmsl_cf"}
        }

    def read_window(self, x0, y0, nx, ny):
        values = {}
        valid_all = None
        for key, src in self.sources.items():
            data, valid = _as_data_valid(src.read_window(x0, y0, nx, ny))
            values[key] = data
            valid_all = valid if valid_all is None else (valid_all & valid)

        vcm_obs = values["vcm_cf"] > 0
        sl_obs = values["vcmsl_cf"] > 0
        common = vcm_obs & sl_obs

        out = np.zeros(values["vcm_avg"].shape, dtype=np.float32)
        if self.mode == "vcm_observed_area":
            out[vcm_obs] = 1.0
        elif self.mode == "vcmsl_observed_area":
            out[sl_obs] = 1.0
        elif self.mode == "common_observed_area":
            out[common] = 1.0
        elif self.mode == "vcm_observed_radiance":
            out[vcm_obs] = values["vcm_avg"][vcm_obs]
        elif self.mode == "vcmsl_observed_radiance":
            out[sl_obs] = values["vcmsl_avg"][sl_obs]
        elif self.mode == "vcm_common_radiance":
            out[common] = values["vcm_avg"][common]
        elif self.mode == "vcmsl_common_radiance":
            out[common] = values["vcmsl_avg"][common]
        else:  # pragma: no cover
            raise RuntimeError(self.mode)

        if valid_all is not None and not valid_all.all():
            return np.ma.masked_array(out, mask=~valid_all, copy=False)
        return out


def discover_monthly_files(root: Path, years: list[int]) -> dict[int, dict[str, Path]]:
    found: dict[int, dict[str, Path]] = {}
    for path in sorted(root.glob("*.avg_rade9h.tif")):
        m = AVG_RE.match(path.name)
        if not m:
            continue
        year = int(m.group("start")[:4])
        if year not in years:
            continue
        product = m.group("product")
        if year in found and product in found[year]:
            raise RuntimeError(
                f"Duplicate {product} avg_rade9h raster for {year}: "
                f"{found[year][product]} and {path}"
            )
        found.setdefault(year, {})[product] = path

    missing = []
    for year in years:
        for product in ("vcmcfg", "vcmslcfg"):
            if product not in found.get(year, {}):
                missing.append(f"{year}:{product}")
    if missing:
        raise FileNotFoundError(
            "Missing monthly radiance rasters: " + ", ".join(missing)
        )

    result: dict[int, dict[str, Path]] = {}
    for year in years:
        result[year] = {}
        for product in ("vcmcfg", "vcmslcfg"):
            avg = found[year][product]
            cf = Path(str(avg).replace(".avg_rade9h.tif", ".cf_cvg.tif"))
            cvg = Path(str(avg).replace(".avg_rade9h.tif", ".cvg.tif"))
            for p in (avg, cf, cvg):
                if not p.exists():
                    raise FileNotFoundError(f"Missing monthly companion raster: {p}")
            result[year][f"{product}_avg"] = avg
            result[year][f"{product}_cf"] = cf
            result[year][f"{product}_cvg"] = cvg
    return result


def parse_date_range(path: Path) -> tuple[str, str, int]:
    m = AVG_RE.match(path.name)
    if not m:
        raise RuntimeError(f"Cannot parse monthly filename: {path.name}")
    start = pd.Timestamp(m.group("start"))
    end = pd.Timestamp(m.group("end"))
    return start.date().isoformat(), end.date().isoformat(), int((end - start).days + 1)


def validate_alignment(files: dict[int, dict[str, Path]], expected_crs: str):
    reference = None
    rows = []
    for year in sorted(files):
        for key, path in sorted(files[year].items()):
            with rasterio.open(path) as ds:
                if ds.count != 1:
                    raise RuntimeError(f"Expected one band: {path}")
                if ds.crs is None or ds.crs.to_string() != expected_crs:
                    raise RuntimeError(
                        f"Unexpected CRS for {path.name}: {ds.crs}; "
                        f"expected {expected_crs}"
                    )
                signature = (
                    ds.width,
                    ds.height,
                    tuple(ds.transform),
                    tuple(ds.bounds),
                    ds.crs.to_string(),
                )
                if reference is None:
                    reference = signature
                elif signature != reference:
                    raise RuntimeError(
                        "Monthly pilot rasters are not perfectly aligned. "
                        f"Mismatch at {path.name}."
                    )
                rows.append(
                    {
                        "year": year,
                        "layer_key": key,
                        "filename": path.name,
                        "size_bytes": int(path.stat().st_size),
                        "width": ds.width,
                        "height": ds.height,
                        "crs": ds.crs.to_string(),
                        "nodata": ds.nodata,
                    }
                )
    assert reference is not None
    return reference, pd.DataFrame(rows)


def integration_self_test():
    """Real exactextract test of support and common-pixel product contrast."""
    shape = (2, 2)
    arrays = {
        "vcm_avg": np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
        "vcm_cf": np.array([[1.0, 0.0], [1.0, 0.0]], dtype=np.float32),
        "vcm_cvg": np.ones(shape, dtype=np.float32) * 2,
        "vcmsl_avg": np.array([[1.5, 5.0], [2.5, 6.0]], dtype=np.float32),
        "vcmsl_cf": np.array([[1.0, 1.0], [1.0, 1.0]], dtype=np.float32),
        "vcmsl_cvg": np.ones(shape, dtype=np.float32) * 3,
    }
    memories = {}
    datasets = {}
    try:
        for key, arr in arrays.items():
            mem = MemoryFile()
            ds = mem.open(
                driver="GTiff",
                height=2,
                width=2,
                count=1,
                dtype="float32",
                crs="EPSG:4326",
                transform=from_origin(0, 2, 1, 1),
                nodata=None,
            )
            ds.write(arr, 1)
            memories[key] = mem
            datasets[key] = ds

        gdf = gpd.GeoDataFrame(
            {"NUTS_ID": ["TEST"]},
            geometry=[box(0, 0, 2, 2)],
            crs="EPSG:4326",
        )
        out = extract_pair_from_open_datasets(
            datasets, gdf, "NUTS_ID", progress=False
        ).iloc[0]
        if not np.isclose(out["vcm_observed_area_fraction"], 0.5, atol=0.01):
            raise RuntimeError("Self-test VCM observed-area fraction failed.")
        if not np.isclose(out["vcmsl_observed_area_fraction"], 1.0, atol=0.01):
            raise RuntimeError("Self-test VCMSL observed-area fraction failed.")
        if not np.isclose(out["common_observed_area_fraction"], 0.5, atol=0.01):
            raise RuntimeError("Self-test common-area fraction failed.")
        if not (out["vcmsl_added_area_fraction"] > 0.45):
            raise RuntimeError("Self-test VCMSL added support failed.")
        if not np.isclose(
            out["common_support_vcmsl_minus_vcm_nw"],
            0.0,
            atol=0.05,
        ):
            # Common pixels are 1.0->1.5 and 3.0->2.5, averaging to zero.
            raise RuntimeError("Self-test common-support radiance contrast failed.")
        if abs(out["product_difference_reconstruction_residual_nw"]) > 1e-6:
            raise RuntimeError("Self-test product decomposition failed.")
    finally:
        for ds in datasets.values():
            ds.close()
        for mem in memories.values():
            mem.close()

    print("July VCM/VCMSL integration self-test: PASS", flush=True)


def extract_pair_from_open_datasets(datasets, nuts, id_field, *, progress):
    sources = [
        MonthlyPairSource(datasets, mode=mode, name=mode)
        for mode in DERIVED_MODES
    ] + [
        RasterioRasterSource(datasets["vcm_cf"], 1, name="vcm_cf_cvg"),
        RasterioRasterSource(datasets["vcm_cvg"], 1, name="vcm_cvg"),
        RasterioRasterSource(datasets["vcmsl_cf"], 1, name="vcmsl_cf_cvg"),
        RasterioRasterSource(datasets["vcmsl_cvg"], 1, name="vcmsl_cvg"),
    ]

    result = exact_extract(
        sources,
        nuts,
        ["mean(coverage_weight=area_spherical_m2)"],
        include_cols=[id_field],
        output="pandas",
        progress=progress,
    )
    rename = {
        f"{mode}_mean": mode
        for mode in DERIVED_MODES
    }
    rename.update(
        {
            "vcm_cf_cvg_mean": "vcm_cf_cvg_mean",
            "vcm_cvg_mean": "vcm_cvg_mean",
            "vcmsl_cf_cvg_mean": "vcmsl_cf_cvg_mean",
            "vcmsl_cvg_mean": "vcmsl_cvg_mean",
        }
    )
    result = result.rename(columns=rename)

    required = {id_field} | set(DERIVED_MODES) | {
        "vcm_cf_cvg_mean",
        "vcm_cvg_mean",
        "vcmsl_cf_cvg_mean",
        "vcmsl_cvg_mean",
    }
    missing = required - set(result.columns)
    if missing:
        raise RuntimeError(
            f"Unexpected exactextract schema; missing {sorted(missing)}; "
            f"found {list(result.columns)}"
        )

    d = result[[id_field] + sorted(required - {id_field})].copy()
    d = derive_pair_metrics(d)
    return d


def _safe_divide(num: pd.Series, den: pd.Series) -> pd.Series:
    n = pd.to_numeric(num, errors="coerce").to_numpy(float)
    d = pd.to_numeric(den, errors="coerce").to_numpy(float)
    out = np.full(len(n), np.nan, dtype=float)
    ok = np.isfinite(n) & np.isfinite(d) & (d > 0)
    out[ok] = n[ok] / d[ok]
    return pd.Series(out, index=num.index)


def derive_pair_metrics(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    for col in DERIVED_MODES + [
        "vcm_cf_cvg_mean",
        "vcm_cvg_mean",
        "vcmsl_cf_cvg_mean",
        "vcmsl_cvg_mean",
    ]:
        d[col] = pd.to_numeric(d[col], errors="coerce")

    tol = 2e-6
    for col in [
        "vcm_observed_area",
        "vcmsl_observed_area",
        "common_observed_area",
    ]:
        bad = (d[col] < -tol) | (d[col] > 1 + tol)
        if bad.any():
            raise RuntimeError(f"Invalid area fraction in {col}.")
        d[col] = d[col].clip(0, 1)

    if (
        d["common_observed_area"]
        > d["vcm_observed_area"] + tol
    ).any() or (
        d["common_observed_area"]
        > d["vcmsl_observed_area"] + tol
    ).any():
        raise RuntimeError("Common monthly support exceeds product support.")

    d["vcm_observed_area_fraction"] = d["vcm_observed_area"]
    d["vcmsl_observed_area_fraction"] = d["vcmsl_observed_area"]
    d["common_observed_area_fraction"] = d["common_observed_area"]
    d["vcmsl_added_area_fraction"] = (
        d["vcmsl_observed_area"] - d["common_observed_area"]
    ).clip(lower=0)
    d["vcm_only_area_fraction"] = (
        d["vcm_observed_area"] - d["common_observed_area"]
    ).clip(lower=0)

    d["vcm_conditional_mean_nw"] = _safe_divide(
        d["vcm_observed_radiance"], d["vcm_observed_area"]
    )
    d["vcmsl_conditional_mean_nw"] = _safe_divide(
        d["vcmsl_observed_radiance"], d["vcmsl_observed_area"]
    )
    d["vcm_common_mean_nw"] = _safe_divide(
        d["vcm_common_radiance"], d["common_observed_area"]
    )
    d["vcmsl_common_mean_nw"] = _safe_divide(
        d["vcmsl_common_radiance"], d["common_observed_area"]
    )
    d["common_support_vcmsl_minus_vcm_nw"] = (
        d["vcmsl_common_mean_nw"] - d["vcm_common_mean_nw"]
    )

    d["vcmsl_added_radiance_contribution_nw"] = (
        d["vcmsl_observed_radiance"] - d["vcmsl_common_radiance"]
    )
    d["vcm_only_radiance_contribution_nw"] = (
        d["vcm_observed_radiance"] - d["vcm_common_radiance"]
    )
    d["vcmsl_added_mean_nw"] = _safe_divide(
        d["vcmsl_added_radiance_contribution_nw"],
        d["vcmsl_added_area_fraction"],
    )
    d["vcm_only_mean_nw"] = _safe_divide(
        d["vcm_only_radiance_contribution_nw"],
        d["vcm_only_area_fraction"],
    )

    d["common_radiance_difference_contribution_nw"] = (
        d["vcmsl_common_radiance"] - d["vcm_common_radiance"]
    )
    d["product_difference_contribution_nw"] = (
        d["vcmsl_observed_radiance"] - d["vcm_observed_radiance"]
    )
    d["product_difference_reconstructed_nw"] = (
        d["common_radiance_difference_contribution_nw"]
        + d["vcmsl_added_radiance_contribution_nw"]
        - d["vcm_only_radiance_contribution_nw"]
    )
    d["product_difference_reconstruction_residual_nw"] = (
        d["product_difference_contribution_nw"]
        - d["product_difference_reconstructed_nw"]
    )
    max_resid = np.nanmax(
        np.abs(d["product_difference_reconstruction_residual_nw"].to_numpy(float))
    )
    if max_resid > 1e-6:
        raise RuntimeError(
            f"Monthly product-difference decomposition failed; "
            f"max residual={max_resid:.3g}"
        )

    d["vcm_cloudfree_share"] = _safe_divide(
        d["vcm_cf_cvg_mean"], d["vcm_cvg_mean"]
    )
    d["vcmsl_cloudfree_share"] = _safe_divide(
        d["vcmsl_cf_cvg_mean"], d["vcmsl_cvg_mean"]
    )
    d["vcmsl_minus_vcm_observed_area_pp"] = 100 * (
        d["vcmsl_observed_area_fraction"] - d["vcm_observed_area_fraction"]
    )
    d["vcmsl_added_area_pp"] = 100 * d["vcmsl_added_area_fraction"]
    return d


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

    if nuts.empty:
        raise RuntimeError("No locked NUTS3 polygons are fully inside the monthly tile.")
    return nuts, excluded


def source_signature(paths: dict[str, Path]) -> str:
    return "|".join(
        f"{key}:{path.name}:{path.stat().st_size}"
        for key, path in sorted(paths.items())
    )


def extract_or_load_year(
    year: int,
    paths: dict[str, Path],
    nuts: gpd.GeoDataFrame,
    config: dict,
    cache_dir: Path,
):
    id_field = config["nuts"]["id_field"]
    cache = cache_dir / f"nuts2016_july_{year}_vcm_vcmsl.parquet"
    signature = source_signature(paths)

    if cache.exists():
        d = pd.read_parquet(cache)
        required = {
            id_field,
            "year",
            "source_signature",
            "common_support_vcmsl_minus_vcm_nw",
            "vcmsl_added_area_pp",
            "product_difference_contribution_nw",
        }
        if required.issubset(d.columns):
            sigs = d["source_signature"].dropna().astype(str).unique().tolist()
            if sigs == [signature] and set(d[id_field].astype(str)) == set(
                nuts[id_field].astype(str)
            ):
                print(f"Using cached July {year} extraction: {cache}", flush=True)
                return d

    datasets = {}
    try:
        key_map = {
            "vcm_avg": paths["vcmcfg_avg"],
            "vcm_cf": paths["vcmcfg_cf"],
            "vcm_cvg": paths["vcmcfg_cvg"],
            "vcmsl_avg": paths["vcmslcfg_avg"],
            "vcmsl_cf": paths["vcmslcfg_cf"],
            "vcmsl_cvg": paths["vcmslcfg_cvg"],
        }
        for key, path in key_map.items():
            datasets[key] = rasterio.open(path)

        print(f"Extracting paired July {year} VCMCFG/VCMSLCFG metrics", flush=True)
        d = extract_pair_from_open_datasets(
            datasets,
            nuts,
            id_field,
            progress=True,
        )
    finally:
        for ds in datasets.values():
            ds.close()

    d.insert(0, "year", int(year))
    start, end, days = parse_date_range(paths["vcmcfg_avg"])
    sl_start, sl_end, sl_days = parse_date_range(paths["vcmslcfg_avg"])
    if (start, end, days) != (sl_start, sl_end, sl_days):
        raise RuntimeError(f"VCM and VCMSL date ranges differ for {year}.")
    d["start_date"] = start
    d["end_date"] = end
    d["n_calendar_days"] = days
    d["vcm_file"] = paths["vcmcfg_avg"].name
    d["vcmsl_file"] = paths["vcmslcfg_avg"].name
    d["source_signature"] = signature

    cache_dir.mkdir(parents=True, exist_ok=True)
    d.to_parquet(cache, index=False)
    print(f"Cached July {year} extraction: {cache}", flush=True)
    return d


def _spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    z = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(z) < 3 or z["x"].nunique() < 2 or z["y"].nunique() < 2:
        return np.nan, np.nan, len(z)
    r = spearmanr(z["x"], z["y"])
    return float(r.statistic), float(r.pvalue), len(z)


def add_analysis_metadata(extracted: pd.DataFrame, panel: pd.DataFrame):
    meta = panel[panel["year"] == 2020][
        [
            "NUTS_ID",
            "CNTR_CODE",
            "density_band",
            "brightness_quintile",
            "baseline_raw_nw",
            "centroid_latitude",
            "centroid_longitude",
        ]
    ].copy()
    if meta["NUTS_ID"].duplicated().any():
        raise RuntimeError("Duplicate NUTS_ID in locked metadata.")

    d = extracted.merge(
        meta,
        on="NUTS_ID",
        how="left",
        validate="many_to_one",
    )
    if d["CNTR_CODE"].isna().any():
        raise RuntimeError("Monthly extraction failed to match locked panel metadata.")
    return d


def country_summary(d: pd.DataFrame, min_common_fraction: float):
    rows = []
    metrics = [
        "vcmsl_added_area_pp",
        "vcmsl_minus_vcm_observed_area_pp",
        "product_difference_contribution_nw",
        "common_radiance_difference_contribution_nw",
        "vcmsl_added_radiance_contribution_nw",
        "vcm_only_radiance_contribution_nw",
        "vcm_cf_cvg_mean",
        "vcmsl_cf_cvg_mean",
        "vcm_cloudfree_share",
        "vcmsl_cloudfree_share",
    ]
    for (year, country), g in d.groupby(["year", "CNTR_CODE"], observed=True):
        common_ok = g["common_observed_area_fraction"] >= min_common_fraction
        row = {
            "year": int(year),
            "CNTR_CODE": str(country),
            "n_nuts3": int(len(g)),
            "n_common_support_regions": int(common_ok.sum()),
            "median_centroid_latitude": float(g["centroid_latitude"].median()),
            "median_vcm_observed_area_pct": float(
                100 * g["vcm_observed_area_fraction"].median()
            ),
            "median_vcmsl_observed_area_pct": float(
                100 * g["vcmsl_observed_area_fraction"].median()
            ),
            "median_common_observed_area_pct": float(
                100 * g["common_observed_area_fraction"].median()
            ),
        }
        for col in metrics:
            row[f"median_{col}"] = float(g[col].median())
        row["median_common_support_vcmsl_minus_vcm_nw"] = (
            float(g.loc[common_ok, "common_support_vcmsl_minus_vcm_nw"].median())
            if common_ok.any()
            else np.nan
        )
        rows.append(row)
    return pd.DataFrame(rows)


def annual_summary(d: pd.DataFrame, countries: pd.DataFrame, min_common_fraction: float):
    rows = []
    for year, g in d.groupby("year", observed=True):
        c = countries[countries["year"] == int(year)].copy()
        common_ok = g["common_observed_area_fraction"] >= min_common_fraction

        rho_added, p_added, n_added = _spearman(
            g["centroid_latitude"], g["vcmsl_added_area_pp"]
        )
        rho_common, p_common, n_common = _spearman(
            g.loc[common_ok, "centroid_latitude"],
            g.loc[common_ok, "common_support_vcmsl_minus_vcm_nw"],
        )

        valid_country_common = c[
            c["median_common_support_vcmsl_minus_vcm_nw"].notna()
        ]
        rows.append(
            {
                "year": int(year),
                "start_date": g["start_date"].iloc[0],
                "end_date": g["end_date"].iloc[0],
                "n_calendar_days": int(g["n_calendar_days"].iloc[0]),
                "n_nuts3": int(len(g)),
                "n_common_support_regions": int(common_ok.sum()),
                "n_countries": int(g["CNTR_CODE"].nunique()),
                "n_countries_with_common_support": int(
                    valid_country_common["CNTR_CODE"].nunique()
                ),
                "country_balanced_median_vcmsl_added_area_pp": float(
                    c["median_vcmsl_added_area_pp"].mean()
                ),
                "country_balanced_median_product_difference_contribution_nw": float(
                    c["median_product_difference_contribution_nw"].mean()
                ),
                "country_balanced_median_common_difference_contribution_nw": float(
                    c["median_common_radiance_difference_contribution_nw"].mean()
                ),
                "country_balanced_median_vcmsl_added_contribution_nw": float(
                    c["median_vcmsl_added_radiance_contribution_nw"].mean()
                ),
                "country_balanced_common_support_vcmsl_minus_vcm_nw": float(
                    valid_country_common[
                        "median_common_support_vcmsl_minus_vcm_nw"
                    ].mean()
                )
                if len(valid_country_common)
                else np.nan,
                "rho_latitude_vs_vcmsl_added_area": rho_added,
                "p_latitude_vs_vcmsl_added_area": p_added,
                "n_latitude_vs_vcmsl_added_area": n_added,
                "rho_latitude_vs_common_support_difference": rho_common,
                "p_latitude_vs_common_support_difference": p_common,
                "n_latitude_vs_common_support_difference": n_common,
            }
        )
    return pd.DataFrame(rows).sort_values("year").reset_index(drop=True)


def density_summary(d: pd.DataFrame, min_common_fraction: float):
    rows = []
    for (year, band), g in d.groupby(["year", "density_band"], observed=True):
        common_ok = g["common_observed_area_fraction"] >= min_common_fraction
        rows.append(
            {
                "year": int(year),
                "density_band": str(band),
                "n_nuts3": int(len(g)),
                "n_common_support_regions": int(common_ok.sum()),
                "median_vcm_observed_area_pct": float(
                    100 * g["vcm_observed_area_fraction"].median()
                ),
                "median_vcmsl_observed_area_pct": float(
                    100 * g["vcmsl_observed_area_fraction"].median()
                ),
                "median_vcmsl_added_area_pp": float(
                    g["vcmsl_added_area_pp"].median()
                ),
                "median_product_difference_contribution_nw": float(
                    g["product_difference_contribution_nw"].median()
                ),
                "median_common_support_vcmsl_minus_vcm_nw": (
                    float(
                        g.loc[
                            common_ok,
                            "common_support_vcmsl_minus_vcm_nw",
                        ].median()
                    )
                    if common_ok.any()
                    else np.nan
                ),
                "median_vcm_cf_cvg_mean": float(g["vcm_cf_cvg_mean"].median()),
                "median_vcmsl_cf_cvg_mean": float(
                    g["vcmsl_cf_cvg_mean"].median()
                ),
            }
        )
    return pd.DataFrame(rows)


def make_plots(
    d: pd.DataFrame,
    annual: pd.DataFrame,
    countries: pd.DataFrame,
    output_dir: Path,
    dpi: int,
    min_common_fraction: float,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    ax.plot(
        annual["year"],
        annual["country_balanced_median_vcmsl_added_area_pp"],
        marker="o",
    )
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(annual["year"])
    ax.set_xlabel("July composite year")
    ax.set_ylabel("Country-balanced median VCMSLCFG-only observed area (pp)")
    ax.set_title("Extra July coverage admitted by VCMSLCFG")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "july_vcmsl_added_coverage.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    ax.plot(
        annual["year"],
        annual["country_balanced_common_support_vcmsl_minus_vcm_nw"],
        marker="o",
    )
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(annual["year"])
    ax.set_xlabel("July composite year")
    ax.set_ylabel(
        "Country-balanced VCMSLCFG - VCMCFG on common support "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title("Paired July radiance difference on identical observed pixels")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "july_common_support_radiance_difference.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    for year, g in d.groupby("year", observed=True):
        ax.scatter(
            g["centroid_latitude"],
            g["vcmsl_added_area_pp"],
            s=10,
            alpha=0.35,
            label=str(year),
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("NUTS3 centroid latitude (degrees N)")
    ax.set_ylabel("VCMSLCFG-only observed area (percentage points)")
    ax.set_title("Latitude dependence of July VCMSLCFG added support")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "july_added_coverage_vs_latitude.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    for year, g in d.groupby("year", observed=True):
        g = g[g["common_observed_area_fraction"] >= min_common_fraction]
        ax.scatter(
            g["centroid_latitude"],
            g["common_support_vcmsl_minus_vcm_nw"],
            s=10,
            alpha=0.35,
            label=str(year),
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("NUTS3 centroid latitude (degrees N)")
    ax.set_ylabel(
        "VCMSLCFG - VCMCFG on common support (nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title("Latitude dependence of paired July radiance difference")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "july_common_difference_vs_latitude.pdf", dpi=dpi)
    plt.close(fig)

    c21 = countries[countries["year"] == 2021].dropna(
        subset=["median_common_support_vcmsl_minus_vcm_nw"]
    ).sort_values("median_centroid_latitude")
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    ax.scatter(
        c21["median_centroid_latitude"],
        c21["median_common_support_vcmsl_minus_vcm_nw"],
        s=32,
    )
    for _, r in c21.iterrows():
        ax.annotate(
            r["CNTR_CODE"],
            (
                r["median_centroid_latitude"],
                r["median_common_support_vcmsl_minus_vcm_nw"],
            ),
            fontsize=7,
            xytext=(3, 3),
            textcoords="offset points",
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Country median NUTS3 centroid latitude (degrees N)")
    ax.set_ylabel(
        "2021 median VCMSLCFG - VCMCFG on common support "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title("July 2021 paired product contrast by country")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "july_2021_country_common_difference.pdf", dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    integration_self_test()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_july_vcm_vcmsl_pilot"]
    years = [int(y) for y in cfg["years"]]
    min_common_fraction = float(cfg["min_common_area_fraction"])
    root = PROJECT_ROOT / cfg["raw_dir"]

    if not root.exists():
        raise FileNotFoundError(
            f"Monthly pilot directory does not exist: {root}"
        )

    files = discover_monthly_files(root, years)
    signature, manifest = validate_alignment(files, config["nuts"]["crs"])

    panel = add_latitude(attach_fixed_brightness(load_locked_panel(config)), config)
    locked_ids = set(panel["NUTS_ID"].astype(str))

    ref_path = files[years[0]]["vcmcfg_avg"]
    with rasterio.open(ref_path) as ref:
        bounds = ref.bounds
    nuts, excluded_ids = load_geometry_subset(config, bounds, locked_ids)

    print(
        f"Monthly tile fully contains {len(nuts)} of {len(locked_ids)} locked "
        f"NUTS3 regions; {len(excluded_ids)} excluded at tile edges/outside.",
        flush=True,
    )

    processed = configured_path(config, "processed")
    cache_dir = processed / "viirs_monthly_july_pilot_cache"

    frames = []
    for year in years:
        frames.append(
            extract_or_load_year(
                year,
                files[year],
                nuts,
                config,
                cache_dir,
            )
        )
    extracted = pd.concat(frames, ignore_index=True)
    regional = add_analysis_metadata(extracted, panel)

    # Restore intended ordered density categories after CSV/parquet roundtrips.
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    regional["density_band"] = pd.Categorical(
        regional["density_band"],
        categories=labels,
        ordered=True,
    )

    countries = country_summary(regional, min_common_fraction)
    annual = annual_summary(regional, countries, min_common_fraction)
    density = density_summary(regional, min_common_fraction)

    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_july_vcm_vcmsl_pilot"
    )
    processed.mkdir(parents=True, exist_ok=True)

    outputs = {
        "manifest": processed / "fixed_nuts2016_july_vcm_vcmsl_manifest.csv",
        "regional": processed / "fixed_nuts2016_july_vcm_vcmsl_regional.csv",
        "annual": processed / "fixed_nuts2016_july_vcm_vcmsl_annual_summary.csv",
        "country": processed / "fixed_nuts2016_july_vcm_vcmsl_by_country.csv",
        "density": processed / "fixed_nuts2016_july_vcm_vcmsl_by_density.csv",
        "summary": processed / "fixed_nuts2016_july_vcm_vcmsl_summary.txt",
    }

    manifest["tile_fully_covered_nuts3"] = len(nuts)
    manifest["tile_excluded_locked_nuts3"] = len(excluded_ids)
    manifest.to_csv(outputs["manifest"], index=False)
    regional.to_csv(outputs["regional"], index=False)
    annual.to_csv(outputs["annual"], index=False)
    countries.to_csv(outputs["country"], index=False)
    density.to_csv(outputs["density"], index=False)

    make_plots(
        regional,
        annual,
        countries,
        figures,
        int(cfg["figure_dpi"]),
        min_common_fraction,
    )

    a = annual.set_index("year")
    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 JULY VCMCFG / VCMSLCFG PILOT\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            "Primary comparison: within-year VCMSLCFG versus VCMCFG. "
            "Observation support is cf_cvg > 0; avg_rade9h zeros are retained "
            "as data where cf_cvg > 0.\n"
        )
        h.write(
            "The 2022 archive ends 2022-07-26, so absolute cross-year July "
            "means/coverage counts are not assumed to have identical temporal "
            "sampling. Paired within-year product contrasts remain directly "
            "matched in date range.\n"
        )
        h.write(
            f"Monthly tile fully contains {len(nuts)} of {len(locked_ids)} "
            f"locked NUTS3 regions; {len(excluded_ids)} are excluded because "
            "their geometry is not wholly inside the downloaded 75N060W tile.\n"
        )
        h.write(
            f"Common-support conditional radiance summaries require at least "
            f"{100*min_common_fraction:.1f}% of NUTS3 area observed in both "
            "products.\n\n"
        )

        h.write("ANNUAL PAIRED PRODUCT SUMMARY\n")
        h.write(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))

        h.write("\n\n2021 COUNTRY MEDIANS (SORTED BY LATITUDE)\n")
        c21 = countries[countries["year"] == 2021].sort_values(
            "median_centroid_latitude"
        )
        keep = [
            "CNTR_CODE",
            "n_nuts3",
            "n_common_support_regions",
            "median_centroid_latitude",
            "median_vcm_observed_area_pct",
            "median_vcmsl_observed_area_pct",
            "median_vcmsl_added_area_pp",
            "median_common_support_vcmsl_minus_vcm_nw",
            "median_product_difference_contribution_nw",
        ]
        h.write(c21[keep].to_string(index=False, float_format=lambda x: f"{x:.6g}"))

        if all(y in a.index for y in [2020, 2021, 2022]):
            flank_common = np.nanmean(
                [
                    a.loc[2020, "country_balanced_common_support_vcmsl_minus_vcm_nw"],
                    a.loc[2022, "country_balanced_common_support_vcmsl_minus_vcm_nw"],
                ]
            )
            flank_added = np.nanmean(
                [
                    a.loc[2020, "country_balanced_median_vcmsl_added_area_pp"],
                    a.loc[2022, "country_balanced_median_vcmsl_added_area_pp"],
                ]
            )
            h.write("\n\nPRE-SPECIFIED 2021 CONTRAST\n")
            h.write(
                f"- 2021 country-balanced common-support VCMSLCFG-VCMCFG: "
                f"{a.loc[2021, 'country_balanced_common_support_vcmsl_minus_vcm_nw']:+.6f} nW.\n"
            )
            h.write(
                f"- Mean of 2020 and 2022 common-support contrasts: "
                f"{flank_common:+.6f} nW.\n"
            )
            h.write(
                f"- 2021 excess over flanking mean: "
                f"{a.loc[2021, 'country_balanced_common_support_vcmsl_minus_vcm_nw'] - flank_common:+.6f} nW.\n"
            )
            h.write(
                f"- 2021 country-balanced VCMSLCFG-only added area: "
                f"{a.loc[2021, 'country_balanced_median_vcmsl_added_area_pp']:+.3f} pp.\n"
            )
            h.write(
                f"- Mean of 2020 and 2022 added area: "
                f"{flank_added:+.3f} pp.\n"
            )
            h.write(
                f"- 2021 latitude correlation with added support: "
                f"rho={a.loc[2021, 'rho_latitude_vs_vcmsl_added_area']:+.4f}.\n"
            )
            h.write(
                f"- 2021 latitude correlation with common-support product "
                f"difference: "
                f"rho={a.loc[2021, 'rho_latitude_vs_common_support_difference']:+.4f}.\n"
            )

        h.write("\nDECISION RULE\n")
        h.write(
            "- A 2021-specific positive VCMSLCFG-VCMCFG radiance contrast on "
            "common observed pixels, especially if it strengthens with "
            "latitude, supports a stray-light-correction radiance pathway.\n"
        )
        h.write(
            "- A 2021-specific expansion of VCMSLCFG-only observed area and/or "
            "large radiance contribution from that added support supports a "
            "stray-light-related sampling/support pathway.\n"
        )
        h.write(
            "- If 2021 resembles both flanking years in both diagnostics, "
            "July stray-light processing is not a sufficient explanation and "
            "the next step should expand to other months or broader DNB "
            "calibration/compositing changes.\n"
        )
        h.write(
            "- This three-July pilot identifies whether the hypothesis merits "
            "a full monthly analysis; it does not itself define an annual "
            "correction.\n"
        )

    print("\nFIXED NUTS2016 JULY VCMCFG / VCMSLCFG PILOT")
    print(
        f"Tile-covered locked NUTS3: {len(nuts)}; "
        f"excluded outside/edge-crossing: {len(excluded_ids)}"
    )
    print("\nANNUAL PAIRED PRODUCT SUMMARY")
    print(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
