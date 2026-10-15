#!/usr/bin/env python3
"""Test whether surface albedo predicts the 2021 VIIRS anomaly and attenuates snow.

This stage is deliberately narrow and pre-specified.  It uses AppEEARS
VNP43MA3.002 daily visible broadband albedo for January-March 2020-2022 and
asks whether:

1. February 2021 surface-albedo anomaly predicts the matched-support VIIRS
   anomaly;
2. adding albedo materially improves the existing snow model; and
3. the snow block attenuates once albedo is included.

February is primary, January supportive and March a negative control.  White-
sky visible albedo is primary; black-sky visible albedo is a sensitivity.
Mandatory quality 0/1 is primary so the analysis does not condition winter
coverage on full BRDF inversion alone; quality 0 only is a stricter sensitivity.

The primary outcome is absolute VIIRS anomaly in nW cm-2 sr-1.  Relative
percentage anomaly is retained explicitly as context because large percentage
changes in dark regions can correspond to very small absolute radiance changes.

No annual VIIRS data are modified by this stage.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource
from rasterio.windows import Window
from scipy.stats import spearmanr
import statsmodels.formula.api as smf
from shapely.geometry import box

from project_config import PROJECT_ROOT, configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import load_fixed_nuts


YEARS = (2020, 2021, 2022)
MONTHS = (1, 2, 3)
MONTH_NAMES = {1: "Jan", 2: "Feb", 3: "Mar"}
ROLE = {1: "supportive", 2: "primary", 3: "negative_control"}

BAND_NAMES = {
    1: "wsa_qa01_mean",
    2: "bsa_qa01_mean",
    3: "wsa_qa0_mean",
    4: "bsa_qa0_mean",
    5: "wsa_qa01_day_fraction",
    6: "bsa_qa01_day_fraction",
    7: "wsa_qa0_day_fraction",
    8: "bsa_qa0_day_fraction",
    9: "wsa_qa01_valid_area",
    10: "bsa_qa01_valid_area",
    11: "wsa_qa0_valid_area",
    12: "bsa_qa0_valid_area",
}

DATE_PATTERNS = [
    re.compile(r"doy(?P<year>20\d{2})(?P<doy>\d{3})", re.IGNORECASE),
    re.compile(r"\.A(?P<year>20\d{2})(?P<doy>\d{3})\.", re.IGNORECASE),
    re.compile(r"(?P<year>20\d{2})[-_](?P<month>\d{2})[-_](?P<day>\d{2})"),
]


def integration_self_test() -> None:
    names = [
        "VNP43MA3.002_Albedo_WSA_vis_doy2020032_aid0001.tif",
        "VNP43MA3.A2021060.002.Albedo_BSA_vis.tif",
        "VNP43MA3_2022-03-15_BRDF_Albedo_Band_Mandatory_Quality_vis.tif",
    ]
    dates = [parse_date_from_name(x) for x in names]
    expected = [
        pd.Timestamp("2020-02-01"),
        pd.Timestamp("2021-03-01"),
        pd.Timestamp("2022-03-15"),
    ]
    if dates != expected:
        raise RuntimeError(f"VNP43MA3 date-parser self-test failed: {dates}")

    toy = pd.DataFrame({
        "NUTS_ID": ["A", "A", "A", "B", "B", "B"],
        "year": [2020, 2021, 2022, 2020, 2021, 2022],
        "albedo": [0.1, 0.3, 0.2, 0.4, 0.5, 0.6],
    })
    p = toy.pivot(index="NUTS_ID", columns="year", values="albedo")
    wide = p.index.to_frame(index=False, name="NUTS_ID")
    for year in YEARS:
        wide[f"a{year}"] = p[year].to_numpy()
    if not np.allclose(wide["a2020"].to_numpy(), [0.1, 0.4]):
        raise RuntimeError("VNP43MA3 pivot-alignment self-test failed.")

    before, after = 0.4, 0.1
    attenuation = 100.0 * (1.0 - after / before)
    if not np.isclose(attenuation, 75.0):
        raise RuntimeError("Attenuation self-test failed.")

    print("Albedo/VIIRS attenuation integration self-test: PASS", flush=True)


def parse_date_from_name(name: str) -> pd.Timestamp:
    for pattern in DATE_PATTERNS:
        m = pattern.search(name)
        if not m:
            continue
        gd = m.groupdict()
        year = int(gd["year"])
        if "doy" in gd and gd.get("doy"):
            return pd.Timestamp(f"{year}-01-01") + pd.Timedelta(days=int(gd["doy"]) - 1)
        return pd.Timestamp(year=year, month=int(gd["month"]), day=int(gd["day"]))
    raise ValueError(f"Could not parse date from AppEEARS filename: {name}")


def extract_zip_bundles(raw_dir: Path, cache_dir: Path) -> list[Path]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    roots = [raw_dir]
    for zpath in sorted(raw_dir.rglob("*.zip")):
        dest = cache_dir / zpath.stem
        sig_path = dest / ".source_signature.json"
        sig = {"name": zpath.name, "size": int(zpath.stat().st_size)}
        existing = None
        if sig_path.exists():
            try:
                existing = json.loads(sig_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                existing = None
        if existing != sig:
            if dest.exists():
                for p in sorted(dest.rglob("*"), reverse=True):
                    if p.is_file() or p.is_symlink():
                        p.unlink()
                    elif p.is_dir():
                        p.rmdir()
            dest.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(zpath) as z:
                z.extractall(dest)
            sig_path.write_text(json.dumps(sig, indent=2), encoding="utf-8")
        roots.append(dest)
    return roots


def classify_layer(name: str, cfg: dict) -> str | None:
    lower = name.lower()
    tokens = {
        "wsa": str(cfg["layer_wsa"]).lower(),
        "bsa": str(cfg["layer_bsa"]).lower(),
        "qa": str(cfg["layer_qa"]).lower(),
    }
    for key, token in tokens.items():
        if token in lower:
            return key
    return None


def discover_daily_files(raw_dir: Path, cache_extract: Path, cfg: dict):
    roots = extract_zip_bundles(raw_dir, cache_extract)
    found: dict[pd.Timestamp, dict[str, Path]] = defaultdict(dict)
    ignored = []
    for root in roots:
        for p in sorted(root.rglob("*.tif")):
            layer = classify_layer(p.name, cfg)
            if layer is None:
                ignored.append(p)
                continue
            try:
                date = parse_date_from_name(p.name)
            except ValueError:
                continue
            if date.year not in YEARS or date.month not in MONTHS:
                continue
            if layer in found[date]:
                raise RuntimeError(
                    f"Duplicate {layer} raster for {date.date()}: "
                    f"{found[date][layer]} and {p}"
                )
            found[date][layer] = p

    expected_dates = []
    for year in YEARS:
        for month in MONTHS:
            start = pd.Timestamp(year=year, month=month, day=1)
            end = start + pd.offsets.MonthEnd(0)
            expected_dates.extend(pd.date_range(start, end, freq="D").tolist())

    missing = []
    for date in expected_dates:
        layers = found.get(date, {})
        absent = [x for x in ("wsa", "bsa", "qa") if x not in layers]
        if absent:
            missing.append((date, absent))

    if missing:
        preview = "\n".join(
            f"  {d.date()}: missing {','.join(parts)}" for d, parts in missing[:20]
        )
        counts = {
            year: sum(1 for d in expected_dates if d.year == year and d in found)
            for year in YEARS
        }
        raise RuntimeError(
            "Incomplete VNP43MA3 AppEEARS inventory. The stage requires all "
            "January-March daily WSA/BSA/QA rasters for 2020-2022.\n"
            f"Dates with any discovered layer by year: {counts}\n"
            f"First missing dates/layers:\n{preview}\n"
            "Place all three completed AppEEARS JFM bundles under the configured "
            "raw directory; ZIP bundles may be left zipped."
        )

    return {d: found[d] for d in expected_dates}


def raster_signature(paths: dict[pd.Timestamp, dict[str, Path]]) -> dict:
    rows = []
    for date in sorted(paths):
        for layer in ("wsa", "bsa", "qa"):
            p = paths[date][layer]
            rows.append({
                "date": str(date.date()),
                "layer": layer,
                "name": p.name,
                "size": int(p.stat().st_size),
            })
    return {"files": rows}


def validate_alignment(datasets: list[rasterio.io.DatasetReader]) -> None:
    ref = datasets[0]
    for ds in datasets[1:]:
        if (
            ds.width != ref.width
            or ds.height != ref.height
            or ds.transform != ref.transform
            or ds.crs != ref.crs
        ):
            raise RuntimeError(
                "AppEEARS VNP43MA3 rasters are not on one aligned grid. "
                f"Reference={ref.name}; mismatch={ds.name}"
            )


def determine_albedo_scale(path: Path, fallback: float) -> tuple[float, float, str]:
    with rasterio.open(path) as ds:
        scale = float(ds.scales[0]) if ds.scales and np.isfinite(ds.scales[0]) else 1.0
        offset = float(ds.offsets[0]) if ds.offsets and np.isfinite(ds.offsets[0]) else 0.0
        if not np.isclose(scale, 1.0) or not np.isclose(offset, 0.0):
            return scale, offset, "raster_metadata"

        h = min(ds.height, 256)
        w = min(ds.width, 256)
        row0 = max(0, (ds.height - h) // 2)
        col0 = max(0, (ds.width - w) // 2)
        a = ds.read(1, window=Window(col0, row0, w, h), masked=True).astype(float)
        vals = a.compressed()
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            raise RuntimeError(f"Could not inspect albedo scaling in {path}")
        q99 = float(np.nanpercentile(vals, 99))
        if q99 <= 1.5:
            return 1.0, 0.0, "already_physical"
        scaled_q99 = q99 * fallback
        if scaled_q99 <= 1.5:
            return fallback, 0.0, "configured_fallback"
        raise RuntimeError(
            f"Ambiguous albedo scaling for {path.name}: sampled 99th percentile "
            f"is {q99:g}; raster scale metadata are unity and fallback "
            f"{fallback:g} would give {scaled_q99:g}. Inspect with gdalinfo."
        )


def build_monthly_composite(
    year: int,
    month: int,
    daily: dict[pd.Timestamp, dict[str, Path]],
    out: Path,
    meta_path: Path,
    cfg: dict,
) -> tuple[Path, dict]:
    month_dates = sorted(d for d in daily if d.year == year and d.month == month)
    expected_days = int((pd.Timestamp(year=year, month=month, day=1) + pd.offsets.MonthEnd(0)).day)
    if len(month_dates) != expected_days:
        raise RuntimeError(
            f"{year}-{month:02d}: expected {expected_days} daily rasters; "
            f"found {len(month_dates)}."
        )

    signature = {
        "year": year,
        "month": month,
        "min_pixel_valid_day_fraction": float(cfg["min_pixel_valid_day_fraction"]),
        "inputs": [
            {
                "date": str(d.date()),
                **{
                    k: {
                        "name": daily[d][k].name,
                        "size": int(daily[d][k].stat().st_size),
                    }
                    for k in ("wsa", "bsa", "qa")
                },
            }
            for d in month_dates
        ],
    }

    if out.exists() and meta_path.exists():
        try:
            old = json.loads(meta_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            old = None
        if old == signature:
            return out, signature

    fallback = float(cfg["fallback_albedo_scale_factor"])
    scale_info = {}
    for layer in ("wsa", "bsa"):
        p = daily[month_dates[0]][layer]
        scale_info[layer] = determine_albedo_scale(p, fallback)

    wsa_ds = [rasterio.open(daily[d]["wsa"]) for d in month_dates]
    bsa_ds = [rasterio.open(daily[d]["bsa"]) for d in month_dates]
    qa_ds = [rasterio.open(daily[d]["qa"]) for d in month_dates]
    all_ds = [*wsa_ds, *bsa_ds, *qa_ds]
    try:
        validate_alignment(all_ds)
        ref = wsa_ds[0]
        profile = ref.profile.copy()
        profile.update(
            driver="GTiff",
            dtype="float32",
            count=12,
            nodata=-9999.0,
            compress="deflate",
            predictor=3,
        )
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_suffix(".tmp.tif")
        min_day_frac = float(cfg["min_pixel_valid_day_fraction"])

        with rasterio.open(tmp, "w", **profile) as dst:
            for b, name in BAND_NAMES.items():
                dst.set_band_description(b, name)

            for _, window in ref.block_windows(1):
                shape = (int(window.height), int(window.width))
                sums = {
                    "wsa_qa01": np.zeros(shape, dtype=np.float64),
                    "bsa_qa01": np.zeros(shape, dtype=np.float64),
                    "wsa_qa0": np.zeros(shape, dtype=np.float64),
                    "bsa_qa0": np.zeros(shape, dtype=np.float64),
                }
                counts = {k: np.zeros(shape, dtype=np.uint16) for k in sums}

                for i in range(len(month_dates)):
                    qa = qa_ds[i].read(1, window=window, masked=True)
                    qdata = np.asarray(qa.data)
                    qmask = np.ma.getmaskarray(qa) | ~np.isfinite(qdata)

                    for layer, datasets in (("wsa", wsa_ds), ("bsa", bsa_ds)):
                        a = datasets[i].read(1, window=window, masked=True).astype(np.float64)
                        data = np.asarray(a.data)
                        mask = np.ma.getmaskarray(a) | ~np.isfinite(data)
                        scale, offset, _ = scale_info[layer]
                        physical = data * scale + offset
                        base_valid = (
                            ~mask
                            & ~qmask
                            & np.isfinite(physical)
                            & (physical >= 0.0)
                            & (physical <= 1.0)
                        )
                        valid01 = base_valid & (qdata >= 0) & (qdata <= 1)
                        valid0 = base_valid & (qdata == 0)
                        for suffix, valid in (("qa01", valid01), ("qa0", valid0)):
                            key = f"{layer}_{suffix}"
                            sums[key][valid] += physical[valid]
                            counts[key][valid] += 1

                means = {}
                fracs = {}
                valids = {}
                for key in sums:
                    frac = counts[key].astype(np.float64) / float(expected_days)
                    mean = np.full(shape, -9999.0, dtype=np.float32)
                    good = (counts[key] > 0) & (frac >= min_day_frac)
                    mean[good] = (sums[key][good] / counts[key][good]).astype(np.float32)
                    means[key] = mean
                    fracs[key] = frac.astype(np.float32)
                    valids[key] = good.astype(np.float32)

                dst.write(means["wsa_qa01"], 1, window=window)
                dst.write(means["bsa_qa01"], 2, window=window)
                dst.write(means["wsa_qa0"], 3, window=window)
                dst.write(means["bsa_qa0"], 4, window=window)
                dst.write(fracs["wsa_qa01"], 5, window=window)
                dst.write(fracs["bsa_qa01"], 6, window=window)
                dst.write(fracs["wsa_qa0"], 7, window=window)
                dst.write(fracs["bsa_qa0"], 8, window=window)
                dst.write(valids["wsa_qa01"], 9, window=window)
                dst.write(valids["bsa_qa01"], 10, window=window)
                dst.write(valids["wsa_qa0"], 11, window=window)
                dst.write(valids["bsa_qa0"], 12, window=window)

        tmp.replace(out)
        signature["scale_info"] = {
            k: {"scale": v[0], "offset": v[1], "mode": v[2]}
            for k, v in scale_info.items()
        }
        meta_path.write_text(json.dumps(signature, indent=2), encoding="utf-8")
        return out, signature
    finally:
        for ds in all_ds:
            ds.close()


def load_geometry(config: dict, allowed_ids: set[str], bounds: list[float]) -> gpd.GeoDataFrame:
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release).copy()
    id_field = config["nuts"]["id_field"]
    nuts[id_field] = nuts[id_field].astype(str)
    nuts = nuts[nuts[id_field].isin(allowed_ids)].copy()
    extent = box(*bounds)
    fully = nuts.geometry.apply(extent.covers)
    excluded = nuts.loc[~fully, id_field].tolist()
    if excluded:
        print(
            f"Excluding {len(excluded)} NUTS3 regions not fully covered by "
            "the AppEEARS rectangle.",
            flush=True,
        )
    return nuts.loc[fully, [id_field, "geometry"]].copy()


def extract_monthly_zonal(
    tif: Path,
    year: int,
    month: int,
    nuts: gpd.GeoDataFrame,
    config: dict,
    cache_dir: Path,
) -> pd.DataFrame:
    id_field = config["nuts"]["id_field"]
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"vnp43ma3_nuts2016_{year}_{month:02d}.parquet"
    signature = f"{tif.name}:{tif.stat().st_size}"

    if cache.exists():
        d = pd.read_parquet(cache)
        required = {id_field, "year", "month", "source_signature"}
        if required.issubset(d.columns):
            if (
                d["source_signature"].dropna().astype(str).unique().tolist() == [signature]
                and set(d[id_field].astype(str)) == set(nuts[id_field].astype(str))
            ):
                return d

    with rasterio.open(tif) as ds:
        sources = [
            RasterioRasterSource(ds, band, name=name)
            for band, name in BAND_NAMES.items()
        ]
        out = exact_extract(
            sources,
            nuts,
            ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=[id_field],
            output="pandas",
            progress=True,
        )

    rename = {f"{name}_mean": name for name in BAND_NAMES.values()}
    out = out.rename(columns=rename)
    out.insert(0, "month", int(month))
    out.insert(0, "year", int(year))
    out["source_signature"] = signature

    for layer in ("wsa", "bsa"):
        for mode in ("qa01", "qa0"):
            alb = f"{layer}_{mode}_mean"
            day = f"{layer}_{mode}_day_fraction"
            valid = f"{layer}_{mode}_valid_area"
            if alb in out:
                bad = out[alb].notna() & ((out[alb] < -1e-8) | (out[alb] > 1 + 1e-8))
                if bad.any():
                    raise RuntimeError(
                        f"Physical albedo invariant failed for {alb}: "
                        f"{out.loc[bad, [id_field, alb]].head().to_dict('records')}"
                    )
                out.loc[out[alb].notna(), alb] = out.loc[out[alb].notna(), alb].clip(0, 1)
            if day in out:
                out.loc[out[day].notna(), day] = out.loc[out[day].notna(), day].clip(0, 1)
            if valid in out:
                out.loc[out[valid].notna(), valid] = out.loc[out[valid].notna(), valid].clip(0, 1)

    out.to_parquet(cache, index=False)
    return out


def spearman_clean(x: pd.Series, y: pd.Series):
    d = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 3 or d["x"].nunique() < 2 or d["y"].nunique() < 2:
        return np.nan, np.nan, len(d)
    s = spearmanr(d["x"], d["y"])
    return float(s.statistic), float(s.pvalue), len(d)


def build_analysis(
    zonal: pd.DataFrame,
    snow_path: Path,
    monthly_viirs_path: Path,
    cfg: dict,
) -> pd.DataFrame:
    snow = pd.read_csv(snow_path)
    viirs = pd.read_csv(monthly_viirs_path)
    snow["NUTS_ID"] = snow["NUTS_ID"].astype(str)
    viirs["NUTS_ID"] = viirs["NUTS_ID"].astype(str)

    keep_viirs = [
        "NUTS_ID", "month",
        "radiance_2020_common_nw", "radiance_2021_common_nw",
        "radiance_2022_common_nw", "common_support_area_fraction",
    ]
    missing = set(keep_viirs) - set(viirs.columns)
    if missing:
        raise RuntimeError(f"Monthly VIIRS regional table missing {sorted(missing)}.")

    base = snow.merge(
        viirs[keep_viirs],
        on=["NUTS_ID", "month"],
        how="inner",
        validate="one_to_one",
    )
    flank = (base["radiance_2020_common_nw"] + base["radiance_2022_common_nw"]) / 2.0
    base["viirs_flank_mean_nw"] = flank
    base["viirs_anomaly_abs_nw"] = (
        base["radiance_2021_common_nw"] - base["viirs_flank_mean_nw"]
    )
    min_flank = float(cfg["relative_min_flank_radiance_nw"])
    base["viirs_anomaly_relative_pct"] = np.where(
        flank > min_flank,
        100.0 * base["viirs_anomaly_abs_nw"] / flank,
        np.nan,
    )

    rows = []
    min_area = float(cfg["min_region_valid_area_fraction"])
    for layer in ("wsa", "bsa"):
        for mode in ("qa01", "qa0"):
            value = f"{layer}_{mode}_mean"
            valid = f"{layer}_{mode}_valid_area"
            dayfrac = f"{layer}_{mode}_day_fraction"
            p = zonal.pivot(index="NUTS_ID", columns=["year", "month"])
            # Use explicit month subsets below; this pivot object is only used
            # to guarantee a complete year/month key set without label alignment.
            del p

            for month in MONTHS:
                m = zonal[zonal["month"] == month].copy()
                pv = m.pivot(index="NUTS_ID", columns="year", values=value)
                pa = m.pivot(index="NUTS_ID", columns="year", values=valid)
                pdays = m.pivot(index="NUTS_ID", columns="year", values=dayfrac)
                if set(YEARS) - set(pv.columns):
                    raise RuntimeError(
                        f"{layer}/{mode} month {month}: missing years in albedo zonal table."
                    )
                out = pv.index.to_frame(index=False, name="NUTS_ID")
                for year in YEARS:
                    out[f"albedo_{year}"] = pv[year].to_numpy()
                    out[f"albedo_valid_area_{year}"] = pa[year].to_numpy()
                    out[f"albedo_mean_day_fraction_{year}"] = pdays[year].to_numpy()
                out["albedo_anomaly_fraction"] = (
                    out["albedo_2021"]
                    - (out["albedo_2020"] + out["albedo_2022"]) / 2.0
                )
                out["albedo_anomaly_pp"] = 100.0 * out["albedo_anomaly_fraction"]
                out["min_albedo_valid_area_fraction"] = out[
                    [f"albedo_valid_area_{y}" for y in YEARS]
                ].min(axis=1)
                out["month"] = month
                out["albedo_type"] = layer.upper()
                out["quality_mode"] = mode
                out["eligible_albedo"] = (
                    (out["min_albedo_valid_area_fraction"] >= min_area)
                    & out["albedo_anomaly_pp"].notna()
                )
                merged = base[base["month"] == month].merge(
                    out, on=["NUTS_ID", "month"], how="inner", validate="one_to_one"
                )
                rows.append(merged)

    d = pd.concat(rows, ignore_index=True)
    d["month_name"] = d["month"].map(MONTH_NAMES)
    d["role"] = d["month"].map(ROLE)
    return d


def partial_r2(reduced, full) -> float:
    if reduced.ssr <= 0:
        return np.nan
    return max(0.0, float((reduced.ssr - full.ssr) / reduced.ssr))


def block_wald(model, terms: list[str]) -> tuple[float, float]:
    try:
        expr = ", ".join(f"{term} = 0" for term in terms)
        w = model.wald_test(expr, scalar=True)
        return float(np.asarray(w.statistic).squeeze()), float(np.asarray(w.pvalue).squeeze())
    except Exception:
        return np.nan, np.nan


def fit_models(d: pd.DataFrame, outcome: str) -> dict:
    cols = [
        outcome, "snow_anomaly_2021_vs_flanks_pp", "albedo_anomaly_pp",
        "baseline_raw_nw", "centroid_latitude", "CNTR_CODE",
        "cf_cvg_anomaly_2021_vs_flanks",
    ]
    x = d[cols].replace([np.inf, -np.inf], np.nan).dropna().copy()
    x["log10_baseline"] = np.log10(x["baseline_raw_nw"].clip(lower=1e-6))
    x["log10_baseline_c"] = x["log10_baseline"] - x["log10_baseline"].mean()

    base = (
        f"{outcome} ~ log10_baseline_c + centroid_latitude + "
        "cf_cvg_anomaly_2021_vs_flanks + C(CNTR_CODE)"
    )
    snow = (
        f"{outcome} ~ snow_anomaly_2021_vs_flanks_pp * log10_baseline_c + "
        "centroid_latitude + cf_cvg_anomaly_2021_vs_flanks + C(CNTR_CODE)"
    )
    alb = (
        f"{outcome} ~ albedo_anomaly_pp * log10_baseline_c + "
        "centroid_latitude + cf_cvg_anomaly_2021_vs_flanks + C(CNTR_CODE)"
    )
    both = (
        f"{outcome} ~ snow_anomaly_2021_vs_flanks_pp * log10_baseline_c + "
        "albedo_anomaly_pp * log10_baseline_c + centroid_latitude + "
        "cf_cvg_anomaly_2021_vs_flanks + C(CNTR_CODE)"
    )

    m0 = smf.ols(base, data=x).fit(cov_type="HC3")
    m1 = smf.ols(snow, data=x).fit(cov_type="HC3")
    m2 = smf.ols(alb, data=x).fit(cov_type="HC3")
    m3 = smf.ols(both, data=x).fit(cov_type="HC3")

    snow_terms = [
        "snow_anomaly_2021_vs_flanks_pp",
        "snow_anomaly_2021_vs_flanks_pp:log10_baseline_c",
    ]
    alb_terms = [
        "albedo_anomaly_pp",
        "albedo_anomaly_pp:log10_baseline_c",
    ]
    snow_wald = block_wald(m3, snow_terms)
    alb_wald = block_wald(m3, alb_terms)

    snow_before = partial_r2(m0, m1)
    snow_after = partial_r2(m2, m3)
    albedo_after_snow = partial_r2(m1, m3)
    attenuation = (
        100.0 * (1.0 - snow_after / snow_before)
        if np.isfinite(snow_before) and snow_before > 0 and np.isfinite(snow_after)
        else np.nan
    )
    snow_coef_before = float(m1.params["snow_anomaly_2021_vs_flanks_pp"])
    snow_coef_after = float(m3.params["snow_anomaly_2021_vs_flanks_pp"])
    coef_att = (
        100.0 * (1.0 - snow_coef_after / snow_coef_before)
        if snow_coef_before != 0 and np.sign(snow_coef_before) == np.sign(snow_coef_after)
        else np.nan
    )

    # Country-clustered sensitivity for the combined model.
    clustered = smf.ols(both, data=x).fit(
        cov_type="cluster", cov_kwds={"groups": x["CNTR_CODE"]}
    )
    c_alb_wald = block_wald(clustered, alb_terms)
    c_snow_wald = block_wald(clustered, snow_terms)

    return {
        "n": int(m3.nobs),
        "r2_m0_controls": float(m0.rsquared),
        "r2_m1_snow": float(m1.rsquared),
        "r2_m2_albedo": float(m2.rsquared),
        "r2_m3_snow_albedo": float(m3.rsquared),
        "partial_r2_snow_before_albedo": snow_before,
        "partial_r2_snow_after_albedo": snow_after,
        "snow_block_attenuation_pct": attenuation,
        "snow_main_coef_before": snow_coef_before,
        "snow_main_coef_after": snow_coef_after,
        "snow_main_coef_attenuation_pct": coef_att,
        "partial_r2_albedo_after_snow": albedo_after_snow,
        "albedo_main_coef_m3": float(m3.params["albedo_anomaly_pp"]),
        "albedo_main_p_m3_hc3": float(m3.pvalues["albedo_anomaly_pp"]),
        "albedo_x_brightness_coef_m3": float(
            m3.params["albedo_anomaly_pp:log10_baseline_c"]
        ),
        "albedo_x_brightness_p_m3_hc3": float(
            m3.pvalues["albedo_anomaly_pp:log10_baseline_c"]
        ),
        "snow_block_wald_m3_hc3": snow_wald[0],
        "snow_block_p_m3_hc3": snow_wald[1],
        "albedo_block_wald_m3_hc3": alb_wald[0],
        "albedo_block_p_m3_hc3": alb_wald[1],
        "snow_block_p_m3_country_cluster": c_snow_wald[1],
        "albedo_block_p_m3_country_cluster": c_alb_wald[1],
        "albedo_main_p_m3_country_cluster": float(
            clustered.pvalues["albedo_anomaly_pp"]
        ),
        "albedo_x_brightness_p_m3_country_cluster": float(
            clustered.pvalues["albedo_anomaly_pp:log10_baseline_c"]
        ),
    }


def tests_and_models(analysis: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    tests = []
    models = []
    for (month, atype, qmode), g in analysis.groupby(
        ["month", "albedo_type", "quality_mode"], observed=True
    ):
        gv = g[g["eligible_albedo"]].copy()
        if gv.empty:
            continue
        for outcome, scale in [
            ("viirs_anomaly_abs_nw", "absolute_nw"),
            ("viirs_anomaly_relative_pct", "relative_pct"),
        ]:
            rho, p, n = spearman_clean(gv["albedo_anomaly_pp"], gv[outcome])
            country = (
                gv.groupby("CNTR_CODE", observed=True)
                .agg(
                    albedo_anomaly_pp=("albedo_anomaly_pp", "median"),
                    outcome=(outcome, "median"),
                    n_nuts3=("NUTS_ID", "nunique"),
                )
                .reset_index()
            )
            crho, cp, cn = spearman_clean(
                country["albedo_anomaly_pp"], country["outcome"]
            )
            tests.append({
                "month": int(month),
                "month_name": MONTH_NAMES[int(month)],
                "role": ROLE[int(month)],
                "albedo_type": atype,
                "quality_mode": qmode,
                "outcome_scale": scale,
                "n_nuts3": n,
                "rho_region_albedo_vs_viirs": rho,
                "p_region_albedo_vs_viirs": p,
                "n_countries": cn,
                "rho_country_median_albedo_vs_viirs": crho,
                "p_country_median_albedo_vs_viirs": cp,
                "median_albedo_anomaly_pp": float(gv["albedo_anomaly_pp"].median()),
                "median_viirs_outcome": float(gv[outcome].median()),
            })
            fit = fit_models(gv, outcome)
            fit.update({
                "month": int(month),
                "month_name": MONTH_NAMES[int(month)],
                "role": ROLE[int(month)],
                "albedo_type": atype,
                "quality_mode": qmode,
                "outcome_scale": scale,
            })
            models.append(fit)
    return pd.DataFrame(tests), pd.DataFrame(models)


def density_context(analysis: pd.DataFrame) -> pd.DataFrame:
    g = analysis[
        (analysis["albedo_type"] == "WSA")
        & (analysis["quality_mode"] == "qa01")
        & analysis["eligible_albedo"]
    ].copy()
    rows = []
    order = ["<30", "30-100", "100-300", "300-1000", ">=1000"]
    for month in MONTHS:
        gm = g[g["month"] == month]
        for band in order:
            x = gm[gm["density_band"].astype(str) == band]
            if x.empty:
                continue
            rows.append({
                "month": month,
                "month_name": MONTH_NAMES[month],
                "role": ROLE[month],
                "density_band": band,
                "n_nuts3": int(len(x)),
                "median_baseline_raw_nw": float(x["baseline_raw_nw"].median()),
                "median_snow_anomaly_pp": float(
                    x["snow_anomaly_2021_vs_flanks_pp"].median()
                ),
                "median_albedo_anomaly_pp": float(x["albedo_anomaly_pp"].median()),
                "median_viirs_anomaly_abs_nw": float(x["viirs_anomaly_abs_nw"].median()),
                "median_viirs_anomaly_relative_pct": float(
                    x["viirs_anomaly_relative_pct"].median()
                ),
            })
    return pd.DataFrame(rows)


def make_plots(
    analysis: pd.DataFrame,
    tests: pd.DataFrame,
    models: pd.DataFrame,
    density: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    primary = analysis[
        (analysis["month"] == 2)
        & (analysis["albedo_type"] == "WSA")
        & (analysis["quality_mode"] == "qa01")
        & analysis["eligible_albedo"]
    ].copy()

    for outcome, filename, ylabel in [
        (
            "viirs_anomaly_abs_nw",
            "february_wsa_albedo_vs_viirs_anomaly_absolute.pdf",
            "February 2021 VIIRS anomaly (nW cm$^{-2}$ sr$^{-1}$)",
        ),
        (
            "viirs_anomaly_relative_pct",
            "february_wsa_albedo_vs_viirs_anomaly_relative.pdf",
            "February 2021 VIIRS anomaly (% of 2020/2022 flank mean)",
        ),
    ]:
        fig, ax = plt.subplots(figsize=(8.4, 5.8))
        ax.scatter(primary["albedo_anomaly_pp"], primary[outcome], s=12, alpha=0.45)
        ax.axhline(0, linewidth=0.8)
        ax.axvline(0, linewidth=0.8)
        ax.set_xlabel("February 2021 white-sky visible albedo anomaly (percentage points)")
        ax.set_ylabel(ylabel)
        ax.set_title("Primary albedo mechanism test")
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(output_dir / filename, dpi=dpi)
        plt.close(fig)

    m = models[
        (models["month"] == 2)
        & (models["albedo_type"] == "WSA")
        & (models["quality_mode"] == "qa01")
        & (models["outcome_scale"] == "absolute_nw")
    ]
    if len(m) == 1:
        r = m.iloc[0]
        vals = [
            r["partial_r2_snow_before_albedo"],
            r["partial_r2_snow_after_albedo"],
            r["partial_r2_albedo_after_snow"],
        ]
        labels = ["Snow\nbefore albedo", "Snow\nafter albedo", "Albedo\nafter snow"]
        fig, ax = plt.subplots(figsize=(7.5, 5.8))
        ax.bar(labels, vals)
        ax.set_ylabel("Partial $R^2$")
        ax.set_title("February block contributions: snow and albedo")
        ax.grid(alpha=0.2, axis="y")
        fig.tight_layout()
        fig.savefig(output_dir / "february_snow_albedo_partial_r2.pdf", dpi=dpi)
        plt.close(fig)

    t = tests[
        (tests["albedo_type"] == "WSA")
        & (tests["quality_mode"] == "qa01")
        & (tests["outcome_scale"] == "absolute_nw")
    ].sort_values("month")
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.bar(t["month_name"], t["rho_region_albedo_vs_viirs"])
    ax.axhline(0, linewidth=0.8)
    ax.set_ylabel("Spearman rho")
    ax.set_xlabel("Pre-specified test month")
    ax.set_title("Albedo anomaly versus absolute VIIRS anomaly")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(output_dir / "albedo_viirs_rho_jan_feb_mar.pdf", dpi=dpi)
    plt.close(fig)

    feb = density[density["month"] == 2].copy()
    order = ["<30", "30-100", "100-300", "300-1000", ">=1000"]
    feb["density_band"] = pd.Categorical(feb["density_band"], categories=order, ordered=True)
    feb = feb.sort_values("density_band")
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.4))
    axes[0].bar(feb["density_band"].astype(str), feb["median_viirs_anomaly_abs_nw"])
    axes[0].axhline(0, linewidth=0.8)
    axes[0].set_ylabel("Median VIIRS anomaly (nW cm$^{-2}$ sr$^{-1}$)")
    axes[0].set_title("Absolute anomaly")
    axes[0].tick_params(axis="x", rotation=25)
    axes[1].bar(
        feb["density_band"].astype(str),
        feb["median_viirs_anomaly_relative_pct"],
    )
    axes[1].axhline(0, linewidth=0.8)
    axes[1].set_ylabel("Median VIIRS anomaly (%)")
    axes[1].set_title("Relative anomaly")
    axes[1].tick_params(axis="x", rotation=25)
    fig.suptitle("February 2021 anomaly: absolute and proportional views")
    fig.tight_layout()
    fig.savefig(output_dir / "february_viirs_anomaly_absolute_relative_by_density.pdf", dpi=dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    integration_self_test()
    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_albedo_viirs_test"]

    raw_dir = PROJECT_ROOT / cfg["vnp43ma3_raw_dir"]
    if not raw_dir.exists():
        raise FileNotFoundError(
            f"{raw_dir}\nCreate this directory and place the three AppEEARS "
            "JFM 2020/2021/2022 downloads there. ZIP bundles may remain zipped."
        )

    processed = configured_path(config, "processed")
    cache_extract = processed / "vnp43ma3_002_appeears_extracted_cache"
    cache_monthly = processed / "vnp43ma3_002_monthly_composite_cache"
    cache_zonal = processed / "vnp43ma3_002_nuts2016_cache"

    daily = discover_daily_files(raw_dir, cache_extract, cfg)

    snow_path = processed / "fixed_nuts2016_snow_2021_anomaly_regional.csv"
    monthly_viirs_path = processed / "fixed_nuts2016_monthly_2021_localization_regional.csv"
    if not snow_path.exists():
        raise FileNotFoundError(
            f"{snow_path}. Run future-radiance-snow-2021-test first."
        )
    if not monthly_viirs_path.exists():
        raise FileNotFoundError(
            f"{monthly_viirs_path}. Run future-radiance-monthly-2021-localization first."
        )

    snow_ids = set(pd.read_csv(snow_path, usecols=["NUTS_ID"])["NUTS_ID"].astype(str))
    bounds = [float(x) for x in cfg["continental_europe_bounds_4326"]]
    nuts = load_geometry(config, snow_ids, bounds)

    manifest_rows = []
    zonal_frames = []
    for year in YEARS:
        for month in MONTHS:
            out = cache_monthly / f"vnp43ma3_002_{year}_{month:02d}_monthly.tif"
            meta = out.with_suffix(".json")
            print(f"Building/using VNP43MA3 monthly composite {year}-{month:02d}", flush=True)
            tif, sig = build_monthly_composite(
                year, month, daily, out, meta, cfg
            )
            print(f"Extracting VNP43MA3 NUTS3 means {year}-{month:02d}", flush=True)
            z = extract_monthly_zonal(
                tif, year, month, nuts, config, cache_zonal
            )
            zonal_frames.append(z)
            scale_info = sig.get("scale_info", {})
            manifest_rows.append({
                "year": year,
                "month": month,
                "month_name": MONTH_NAMES[month],
                "n_daily_files": int(
                    sum(1 for d in daily if d.year == year and d.month == month)
                ),
                "monthly_composite": tif.name,
                "wsa_scale": scale_info.get("wsa", {}).get("scale", np.nan),
                "wsa_scale_mode": scale_info.get("wsa", {}).get("mode", ""),
                "bsa_scale": scale_info.get("bsa", {}).get("scale", np.nan),
                "bsa_scale_mode": scale_info.get("bsa", {}).get("mode", ""),
            })

    zonal = pd.concat(zonal_frames, ignore_index=True)
    analysis = build_analysis(zonal, snow_path, monthly_viirs_path, cfg)
    tests, models = tests_and_models(analysis)
    density = density_context(analysis)
    manifest = pd.DataFrame(manifest_rows)

    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_albedo_viirs_test"
    )
    figures.mkdir(parents=True, exist_ok=True)

    outputs = {
        "manifest": processed / "fixed_nuts2016_vnp43ma3_manifest.csv",
        "monthly_zonal": processed / "fixed_nuts2016_vnp43ma3_monthly_2020_2022.csv",
        "analysis": processed / "fixed_nuts2016_albedo_snow_viirs_anomaly_regional.csv",
        "tests": processed / "fixed_nuts2016_albedo_2021_anomaly_tests.csv",
        "models": processed / "fixed_nuts2016_albedo_snow_attenuation_models.csv",
        "density": processed / "fixed_nuts2016_albedo_viirs_density_context.csv",
        "summary": processed / "fixed_nuts2016_albedo_viirs_test_summary.txt",
    }
    manifest.to_csv(outputs["manifest"], index=False)
    zonal.to_csv(outputs["monthly_zonal"], index=False)
    analysis.to_csv(outputs["analysis"], index=False)
    tests.to_csv(outputs["tests"], index=False)
    models.to_csv(outputs["models"], index=False)
    density.to_csv(outputs["density"], index=False)

    make_plots(
        analysis, tests, models, density, figures, int(cfg["figure_dpi"])
    )

    key_test = tests[
        (tests["month"] == 2)
        & (tests["albedo_type"] == "WSA")
        & (tests["quality_mode"] == "qa01")
        & (tests["outcome_scale"] == "absolute_nw")
    ].iloc[0]
    key_model = models[
        (models["month"] == 2)
        & (models["albedo_type"] == "WSA")
        & (models["quality_mode"] == "qa01")
        & (models["outcome_scale"] == "absolute_nw")
    ].iloc[0]

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 ALBEDO / SNOW / VIIRS 2021 ATTENUATION TEST\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            "Primary question: does visible surface-albedo anomaly predict the "
            "February 2021 VIIRS anomaly and attenuate the snow block?\n"
        )
        h.write(
            "Pre-specification: February primary; January supportive; March "
            "negative control. White-sky visible albedo primary; black-sky "
            "sensitivity. Mandatory QA 0/1 primary; QA 0-only stricter sensitivity.\n"
        )
        h.write(
            "Primary outcome is absolute radiance anomaly in nW cm-2 sr-1. "
            "Relative percentage anomaly is retained as context because large "
            "percentage changes in dark regions may represent small absolute "
            "radiance changes.\n\n"
        )
        h.write("APPeeARS / VNP43MA3 INVENTORY\n")
        h.write(manifest.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nALBEDO-VIIRS SPEARMAN TESTS\n")
        h.write(tests.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nSEQUENTIAL ATTENUATION MODELS\n")
        h.write(models.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nABSOLUTE / RELATIVE DENSITY CONTEXT\n")
        h.write(density.to_string(index=False, float_format=lambda x: f"{x:.6g}"))

        h.write("\n\nPRIMARY FEBRUARY WSA QA0/1 ABSOLUTE RESULT\n")
        h.write(
            f"- Region-level albedo vs VIIRS rho: "
            f"{key_test['rho_region_albedo_vs_viirs']:+.4f}; "
            f"country-median rho: "
            f"{key_test['rho_country_median_albedo_vs_viirs']:+.4f}.\n"
        )
        h.write(
            f"- M0 controls R2={key_model['r2_m0_controls']:.4f}; "
            f"M1 +snow R2={key_model['r2_m1_snow']:.4f}; "
            f"M2 +albedo R2={key_model['r2_m2_albedo']:.4f}; "
            f"M3 snow+albedo R2={key_model['r2_m3_snow_albedo']:.4f}.\n"
        )
        h.write(
            f"- Snow partial R2 before albedo="
            f"{key_model['partial_r2_snow_before_albedo']:.4f}; "
            f"after albedo={key_model['partial_r2_snow_after_albedo']:.4f}; "
            f"block attenuation="
            f"{key_model['snow_block_attenuation_pct']:.1f}%.\n"
        )
        h.write(
            f"- Albedo partial R2 after snow="
            f"{key_model['partial_r2_albedo_after_snow']:.4f}; "
            f"HC3 joint albedo-block p="
            f"{key_model['albedo_block_p_m3_hc3']:.4g}; "
            f"country-clustered p="
            f"{key_model['albedo_block_p_m3_country_cluster']:.4g}.\n"
        )
        h.write(
            f"- Albedo x baseline-brightness interaction HC3 p="
            f"{key_model['albedo_x_brightness_p_m3_hc3']:.4g}; "
            f"country-clustered p="
            f"{key_model['albedo_x_brightness_p_m3_country_cluster']:.4g}.\n"
        )
        h.write(
            "\nINTERPRETATION RULE\n"
            "- Strong support for the albedo pathway requires February albedo "
            "to predict absolute VIIRS anomaly, add explanatory power beyond "
            "snow, and materially reduce the snow block. January should support "
            "the same direction and March should remain materially weak/absent.\n"
            "- Relative-percent results are contextual/sensitivity results, not "
            "a substitute for the primary absolute-radiance analysis.\n"
            "- This stage does not yet correct the annual 2013-2024 series. If "
            "the mechanism is supported, the next task is to determine how to "
            "control the 2021 surface-reflectance disturbance and then refit the "
            "density-dependent historical response in both absolute and relative "
            "radiance space.\n"
        )

    print("\nFIXED NUTS2016 ALBEDO / SNOW / VIIRS ATTENUATION TEST")
    print("\nPRIMARY FEBRUARY RESULT")
    print(key_test.to_string())
    print("\nPRIMARY FEBRUARY MODEL")
    print(key_model.to_string())
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
