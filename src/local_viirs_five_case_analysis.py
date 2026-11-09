#!/usr/bin/env python3
"""Five predeclared local-lighting cases, monthly EOG VNL V1 2015–2024.

Read-only VIIRS source analysis: no annual corrections, refits, or PRC changes.
Requires the 75N060W VCMSLCFG avg_rade9h + cf_cvg GeoTIFF pairs for all
months 2015–2024. Strict monthly and grid validation; no silent data gaps.
Output: independent CSVs, event-aligned series, native-grid maps and GeoTIFFs.

    python src/local_viirs_five_case_analysis.py --self-test
    python src/local_viirs_five_case_analysis.py --inventory
    python src/local_viirs_five_case_analysis.py
    python src/local_viirs_five_case_analysis.py --cases echirolles portsmouth
"""
from __future__ import annotations

import argparse
import json
import math
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlencode

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask, geometry_window
from rasterio.windows import bounds as window_bounds
from shapely.geometry import mapping, shape

ROOT = Path(__file__).resolve().parents[1]
RAWDIR = ROOT / "data/raw/viirs"
BOUNDARIES = ROOT / "data/raw/local_lighting_boundaries"
NUTS2016 = ROOT / "data/raw/nuts/NUTS_RG_01M_2016_4326_LEVL_3.geojson"
OUT = ROOT / "data/processed/local_lighting_interventions/five_case_monthly"
FIG = ROOT / "figures/local_lighting_interventions/five_case_monthly"
PRODUCT = "vcmslcfg"
TILE = "75N060W"
YEARS = tuple(range(2015, 2025))
MONTHS = tuple(f"{y}-{m:02d}" for y in YEARS for m in range(1, 13))
UNIT = "nW cm$^{-2}$ sr$^{-1}$"
ONS2024 = (
    "https://services1.arcgis.com/ESMARspQHYMw9BZ9/arcgis/rest/services/"
    "Local_Authority_Districts_December_2024_Boundaries_UK_BGC/"
    "FeatureServer/0/query"
)
NAME = re.compile(
    r"^SVDNB_(?:npp|j01)_(?P<start>\d{8})-(?P<end>\d{8})_"
    r"75N060W_vcmslcfg_v10_.*\."
    r"(?P<layer>avg_rade9h|avg_rade9|cf_cvg)\.tif$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class Case:
    name: str
    authority_code: str
    nuts3: str
    transition_start: str
    transition_end: str
    mechanism: str
    contrast_label: str
    treatment_year: int
    reference_years: tuple[int, ...]
    compare_months: tuple[int, ...]
    shown_years: tuple[int, ...]
    comment: str


# These designs were set from independent published council timing, NOT from
# VIIRS outcomes. Transitions are shaded and excluded from primary contrasts.
CASES = {
    "echirolles": Case(
        "Échirolles", "38151", "FRK24", "2022-07", "2022-09",
        "Street-light switch-offs, some sites from July; roads Oct 2022",
        "Oct-Dec 2022 vs 2020/2021 same calendar months",
        2022, (2020, 2021), (10, 11, 12),
        (2019, 2020, 2021, 2022, 2023, 2024),
        "2023 LED conversion contaminates later follow-up, not Oct-Dec 2022"),
    "portsmouth": Case(
        "Portsmouth", "E06000044", "UKJ31", "2020-06", "2020-08",
        "LED conversion 2018-19; adaptive dimming summer 2020",
        "Sep-Dec 2020 vs 2018/2019 same calendar months",
        2020, (2018, 2019), (9, 10, 11, 12),
        (2017, 2018, 2019, 2020, 2021, 2022),
        "2018-19 reference already had gradual LED conversion; sensitivity essential"),
    "coventry": Case(
        "Coventry", "E08000026", "UKG33", "2024-05", "2024-07",
        "Part-night switching introduced May-Jul 2024",
        "Aug-Dec 2024 vs 2022/2023 same calendar months",
        2024, (2022, 2023), (8, 9, 10, 11, 12),
        (2020, 2021, 2022, 2023, 2024),
        "May-Jul implementation months excluded; no 2025 LED programme included"),
    "bristol": Case(
        "Bristol", "E06000023", "UKK11", "2023-03", "2024-12",
        "Rolling LED replacement and CMS dimming from Mar 2023",
        "Apr-Dec 2023 vs 2020-22 same calendar months",
        2023, (2020, 2021, 2022), tuple(range(4, 13)),
        (2019, 2020, 2021, 2022, 2023, 2024),
        "2023 is during rollout, NOT fully treated; LED and operation confounded"),
    "thanet": Case(
        "Thanet", "E07000114", "UKJ44", "2016-12", "2017-05",
        "First LED/all-night conversion Dec 2016–May 2017; 2018 second phase",
        "Jun-Nov 2017 vs 2015/2016 same calendar months",
        2017, (2015, 2016), (6, 7, 8, 9, 10, 11),
        (2015, 2016, 2017, 2018, 2019),
        "Later 2018 lighting changes cause mixed exposure; test sign not flux"),
}

# Deliberately descriptive second contrasts, not additional independent events.
SECONDARY = {
    "echirolles": (2023, (2020, 2021), (10, 11, 12),
                   "2023 follow-up, with LED conversions"),
    "portsmouth": (2021, (2018, 2019), (9, 10, 11, 12),
                   "2021 persistence; lighting/meteorology may have changed"),
    "bristol": (2024, (2020, 2021, 2022), tuple(range(4, 13)),
                "2024 continued rollout, NOT a clean isolated effect"),
    "thanet": (2018, (2015, 2016), (6, 7, 8, 9, 10, 11),
               "2018 second phase not individually dated"),
}
MONTHLY_INDEX = pd.period_range("2015-01", "2024-12", freq="M").astype(str)


def build_inventory():
    """Find original EOG observations, not _cYYYYMMDD creation-date tokens."""
    candidates: dict[tuple[str, str], list[tuple[Path, str]]] = {}
    for path in RAWDIR.rglob("*.tif"):
        match = NAME.match(path.name)
        if not match:
            continue
        date0, date1 = match["start"], match["end"]
        ym = date0[:4] + "-" + date0[4:6]
        if ym not in MONTHS:
            continue
        layer = "rad" if match["layer"].lower().startswith("avg_rade9") else "cf"
        candidates.setdefault((ym, layer), []).append((path, date1))
    result = {}
    missing = []
    for ym in MONTHS:
        chosen = {}
        for layer in ("rad", "cf"):
            paths = candidates.get((ym, layer), [])
            # Prefer the canonical per-year directory populated by the importer.
            canonical = f"monthly_vcmslcfg_{ym[:4]}/tif/"
            paths.sort(key=lambda x: (
                canonical in x[0].as_posix().lower(),
                "monthly_vcmslcfg_" in x[0].as_posix().lower()), reverse=True)
            if not paths:
                missing.append((ym, layer))
                continue
            if len(paths) > 1:
                a, b = paths[:2]
                s = lambda q: (
                    canonical in q[0].as_posix().lower(),
                    "monthly_vcmslcfg_" in q[0].as_posix().lower())
                if s(a) == s(b):
                    raise ValueError(f"Ambiguous {ym} {layer}: {a[0]} and {b[0]}")
            chosen[layer] = paths[0]
        if len(chosen) == 2:
            rad, end_rad = chosen["rad"]
            cov, end_cov = chosen["cf"]
            if end_rad != end_cov:
                raise ValueError(f"Nonidentical monthly observation intervals {ym}: {rad}, {cov}")
            result[ym] = (rad, cov)
    return result, missing


def fetch_json(url):
    req = urllib.request.Request(url, headers={
        "User-Agent": "Bristol-University-VIIRS-five-case-study/1.0",
        "Accept": "application/geo+json, application/json",
    })
    with urllib.request.urlopen(req, timeout=75) as response:
        return json.load(response)


def geojson_geometry(obj, code=None):
    if obj.get("type") == "FeatureCollection":
        features = obj.get("features", [])
        if len(features) != 1:
            raise ValueError(f"Expected one local-authority feature, received {len(features)}")
        obj = features[0]
    if obj.get("type") != "Feature":
        raise ValueError("Expected GeoJSON Feature, not an error page or ArcGIS JSON")
    if code and code.startswith("E"):
        actual = (obj.get("properties") or {}).get("LAD24CD")
        if actual is not None and actual != code:
            raise ValueError(f"Incorrect local authority: expected {code}, got {actual}")
    geom = shape(obj["geometry"])
    if geom.is_empty or not geom.is_valid or geom.geom_type not in ("Polygon", "MultiPolygon"):
        raise ValueError(f"Unusable {code} polygon from official source")
    return geom


def load_boundary(key, case, *, allow_nuts3=False):
    BOUNDARIES.mkdir(parents=True, exist_ok=True)
    if key == "echirolles":
        path = BOUNDARIES / "echirolles_insee_38151.geojson"
        url = ("https://geo.api.gouv.fr/communes/38151?"
               "fields=nom,code,contour&format=geojson&geometry=contour")
        label = "INSEE commune 38151"
    else:
        path = BOUNDARIES / f"ons2024_lad_{case.authority_code}.geojson"
        params = urlencode({
            "where": f"LAD24CD='{case.authority_code}'",
            "outFields": "*", "returnGeometry": "true",
            "outSR": 4326, "f": "geojson"})
        url = ONS2024 + "?" + params
        label = f"ONS 2024 LAD {case.authority_code} (20m generalized BGC)"
    try:
        if path.exists():
            obj = json.loads(path.read_text(encoding="utf-8"))
        else:
            obj = fetch_json(url)
        geom = geojson_geometry(obj, case.authority_code)
        if not path.exists():
            path.write_text(json.dumps(obj, ensure_ascii=False), encoding="utf-8")
            print(f"  Cached official boundary: {path}")
        return geom, label, str(path.relative_to(ROOT))
    except Exception as exc:
        if not allow_nuts3:
            raise RuntimeError(
                f"{key}: cannot obtain/validate official municipality boundary.\n"
                f"Source: {url}\nCache path: {path}\nReason: {exc}\n"
                "Do NOT substitute the full NUTS3 region without explicitly "
                "choosing --allow-nuts3-fallback.") from exc
        if not NUTS2016.exists():
            raise FileNotFoundError(NUTS2016) from exc
        nuts = gpd.read_file(NUTS2016)
        selected = nuts.loc[nuts.NUTS_ID == case.nuts3]
        if len(selected) != 1:
            raise ValueError(f"Missing NUTS3 {case.nuts3}") from exc
        geom = selected.to_crs("EPSG:4326").geometry.iloc[0]
        print(f"  WARNING: {key} uses NUTS3 {case.nuts3}, NOT the local authority")
        return geom, f"FALLBACK NUTS2016 {case.nuts3} (NOT municipality)", str(NUTS2016.relative_to(ROOT))


def grid_month(polygon, radfile, covfile):
    with rasterio.open(radfile) as rd, rasterio.open(covfile) as cv:
        if (rd.crs is None or rd.crs != cv.crs or rd.shape != cv.shape or
                rd.transform != cv.transform):
            raise ValueError(f"Inconsistent EOG grids: {radfile} and {covfile}")
        geom = gpd.GeoSeries([polygon], crs="EPSG:4326").to_crs(rd.crs).iloc[0]
        win = geometry_window(rd, [mapping(geom)])
        trans = rd.window_transform(win)
        inside = geometry_mask(
            [mapping(geom)], out_shape=(int(win.height), int(win.width)),
            transform=trans, invert=True, all_touched=False)
        if not inside.any():
            raise ValueError(f"Zero VIIRS centres inside polygon for {radfile}")
        rad = rd.read(1, window=win, masked=True)
        cf = cv.read(1, window=win, masked=True)
        a = np.asarray(rad.data, dtype=np.float32)
        b = np.asarray(cf.data, dtype=np.float32)
        valid = (inside & ~np.ma.getmaskarray(rad) &
                 ~np.ma.getmaskarray(cf) & np.isfinite(a) &
                 np.isfinite(b) & (b > 0))
        rv = np.full(a.shape, np.nan, dtype=np.float32)
        cvv = np.full(a.shape, np.nan, dtype=np.float32)
        rv[valid] = a[valid]
        cvv[inside & ~np.ma.getmaskarray(cf) & np.isfinite(b)] = b[
            inside & ~np.ma.getmaskarray(cf) & np.isfinite(b)]
        return {"rad": rv, "cf": cvv, "inside": inside,
                "transform": trans, "extent": window_bounds(win, rd.transform),
                "crs": rd.crs, "geometry": geom}


def grid_stats(key, case, ym, grid, label, paths):
    x = grid["rad"]
    cv = grid["cf"]
    good = np.isfinite(x)
    obs = cv[good]
    possible = int(grid["inside"].sum())
    values = x[good]
    result = {
        "case": key, "place": case.name, "support": label, "month": ym,
        "year": int(ym[:4]), "calendar_month": int(ym[5:]),
        "is_transition": case.transition_start <= ym <= case.transition_end,
        "possible_pixels": possible, "valid_pixels": int(good.sum()),
        "valid_percent": float(100 * good.sum() / possible),
        "mean_radiance": float(values.mean()) if len(values) else np.nan,
        "median_radiance": float(np.median(values)) if len(values) else np.nan,
        "q25_radiance": float(np.quantile(values, .25)) if len(values) else np.nan,
        "q75_radiance": float(np.quantile(values, .75)) if len(values) else np.nan,
        "median_cf_cvg": float(np.median(obs)) if len(obs) else np.nan,
        "q25_cf_cvg": float(np.quantile(obs, .25)) if len(obs) else np.nan,
        "q75_cf_cvg": float(np.quantile(obs, .75)) if len(obs) else np.nan,
        "n_cf3": int((obs >= 3).sum()),
        "n_cf5": int((obs >= 5).sum()),
        "mean_cf3": float(x[good & (cv >= 3)].mean())
                    if (good & (cv >= 3)).any() else np.nan,
        "mean_cf5": float(x[good & (cv >= 5)].mean())
                    if (good & (cv >= 5)).any() else np.nan,
        "radiance_tif": str(paths[0].relative_to(ROOT)),
        "cf_cvg_tif": str(paths[1].relative_to(ROOT)),
        "winter_2021_surface_caution": ym in ("2021-01", "2021-02", "2021-03"),
    }
    return result


def _outline(ax, geom):
    for polygon in geom.geoms if hasattr(geom, "geoms") else [geom]:
        if hasattr(polygon, "exterior"):
            ax.plot(*polygon.exterior.xy, color="black", linewidth=.85)


def draw_map(ax, values, base, title, *, norm=None, vmin=None, vmax=None, cmap="viridis"):
    left, bottom, right, top = base["extent"]
    im = ax.imshow(values, extent=(left, right, bottom, top),
                   origin="upper", interpolation="nearest",
                   norm=norm, vmin=vmin, vmax=vmax, cmap=cmap)
    _outline(ax, base["geometry"])
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    return im


def _save_fig(fig, stem):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".pdf"):
        fig.savefig(FIG / f"{stem}{ext}", dpi=180)
    plt.close(fig)


def _save_tif(path, arr, template):
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", driver="GTiff",
                       height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype="float32", crs=template["crs"],
                       transform=template["transform"], nodata=np.nan,
                       compress="deflate") as ds:
        ds.write(np.asarray(arr, dtype=np.float32), 1)


def monthly_contrast(panel, year, references, months):
    """Seasonally matched municipality-month comparisons; NO causal estimator."""
    result = []
    panel = panel.set_index(["year", "calendar_month"])
    for m in months:
        current = panel.loc[(year, m)]
        prior = panel.loc[[(y, m) for y in references]]
        avg = float(prior["mean_radiance"].mean())
        observed = float(current["mean_radiance"])
        delta = observed - avg
        result.append({
            "year": year, "month": m, "reference_years": ",".join(map(str, references)),
            "baseline_mean": avg, "outcome_mean": observed,
            "delta_abs": delta, "delta_percent": 100 * delta / avg if avg > .05 else np.nan,
            "valid_percent_outcome": current["valid_percent"],
            "valid_percent_reference": prior["valid_percent"].mean(),
            "median_cf_cvg_outcome": current["median_cf_cvg"],
            "median_cf_cvg_reference": prior["median_cf_cvg"].mean(),
            "delta_cf3": (float(current["mean_cf3"]) -
                          float(prior["mean_cf3"].mean())),
            "delta_cf5": (float(current["mean_cf5"]) -
                          float(prior["mean_cf5"].mean())),
        })
    return pd.DataFrame(result)


def paired_pixel_maps(grids, year, refs, months):
    """Per-calendar-month common-support baseline/outcome arrays and ΔR.

    A cell is matched only if that SAME cell is observed in the outcome and
    every reference year in the calendar month; a period map needs >=75% of
    requested months supported. No treatment/control matching is inferred.
    """
    months = tuple(months)
    before, after, diffs = [], [], []
    for m in months:
        cur = grids[f"{year}-{m:02d}"]["rad"]
        previous = [grids[f"{y}-{m:02d}"]["rad"] for y in refs]
        stack = np.stack(previous)
        common = np.isfinite(cur) & np.all(np.isfinite(stack), axis=0)
        b = np.full(cur.shape, np.nan, dtype=np.float32)
        a = np.full(cur.shape, np.nan, dtype=np.float32)
        baseline = np.mean(stack, axis=0)
        b[common], a[common] = baseline[common], cur[common]
        before.append(b)
        after.append(a)
        diffs.append(a - b)
    coverage = np.sum(np.isfinite(np.stack(diffs)), axis=0)
    threshold = max(1, math.ceil(len(months) * .75))
    def period_mean(arrays):
        stack = np.stack(arrays)
        denom = np.isfinite(stack).sum(axis=0)
        result = np.full(stack.shape[1:], np.nan, dtype=np.float32)
        np.divide(np.nansum(stack, axis=0), denom, out=result,
                  where=denom >= threshold)
        return result
    return period_mean(before), period_mean(after), period_mean(diffs), diffs, coverage


def plot_spatial(key, case, grids, year, refs, months, description):
    before, after, diff, monthdiff, coverage = paired_pixel_maps(
        grids, year, refs, months)
    base = grids[f"{year}-{months[0]:02d}"]
    finite = np.r_[before[np.isfinite(before)], after[np.isfinite(after)]]
    if not len(finite):
        raise ValueError(f"{key}: no common observed cells for spatial contrast")
    vmax = max(1e-3, float(np.quantile(finite, .995)))
    valid_diff = np.abs(diff[np.isfinite(diff)])
    dmax = max(1e-3, float(np.quantile(valid_diff, .98)))
    div = TwoSlopeNorm(vmin=-dmax, vcenter=0, vmax=dmax)
    stem = f"{key}_{year}_versus_{'-'.join(map(str, refs))}"

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    draw_map(axes[0], before, base, "Prior years, same calendar months",
             vmin=0, vmax=vmax)
    im = draw_map(axes[1], after, base, f"Outcome {year}: months {list(months)}",
                  vmin=0, vmax=vmax)
    imd = draw_map(axes[2], diff, base, f"{year} minus prior-year mean",
                   norm=div, cmap="coolwarm")
    fig.colorbar(im, ax=axes[:2], shrink=.75, label=UNIT)
    fig.colorbar(imd, ax=axes[2], shrink=.75, label=f"Δ radiance ({UNIT})")
    fig.suptitle(f"{case.name}: {description}; SAME monthly pixel support", fontsize=12)
    _save_fig(fig, stem + "_period_maps")

    n = len(months)
    cols = 3
    fig, axes = plt.subplots(math.ceil(n / cols), cols,
                             figsize=(13, 3.8 * math.ceil(n / cols)),
                             constrained_layout=True)
    flat = np.atleast_1d(axes).ravel()
    max_m = np.nanmax(np.abs(np.stack(monthdiff)))
    dnorm = TwoSlopeNorm(vmin=-max(max_m, 1e-3), vcenter=0,
                         vmax=max(max_m, 1e-3))
    for ax, m, z in zip(flat, months, monthdiff):
        im = draw_map(ax, z, base, f"{year}-{m:02d} minus reference",
                      norm=dnorm, cmap="coolwarm")
    for ax in flat[n:]:
        ax.set_visible(False)
    fig.colorbar(im, ax=[a for a in flat[:n]], shrink=.7, label=f"Δ ({UNIT})")
    fig.suptitle(f"{case.name}: calendar-month matched native-grid differences")
    _save_fig(fig, stem + "_monthly_maps")

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    im = draw_map(ax, np.where(base["inside"], coverage, np.nan),
                  base, "Number of paired months with valid observations",
                  vmin=0, vmax=n, cmap="cividis")
    fig.colorbar(im, ax=ax, label=f"Paired months, of {n}")
    fig.suptitle(f"{case.name}: spatial support for period contrast")
    _save_fig(fig, stem + "_matched_coverage_map")
    _save_tif(OUT / f"{stem}_delta.tif", diff, base)
    _save_tif(OUT / f"{stem}_matched_month_count.tif",
              np.where(base["inside"], coverage, np.nan), base)
    return int(np.isfinite(diff).sum()), float(np.nanmedian(diff))


def plot_monthly_series(key, case, panel):
    data = panel.sort_values("month").copy()
    data["date"] = pd.to_datetime(data["month"] + "-01")
    a, b = pd.Timestamp(case.transition_start + "-01"), (
        pd.Timestamp(case.transition_end + "-01") + pd.offsets.MonthBegin(1))
    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True,
                             constrained_layout=True)
    axes[0].plot(data.date, data.mean_radiance, "-",
                 linewidth=1.3, label="Mean radiance")
    axes[0].plot(data.date, data.median_radiance, "-", alpha=.7,
                 linewidth=1, label="Median")
    axes[0].legend()
    axes[0].set_ylabel(f"Radiance ({UNIT})")
    axes[1].plot(data.date, data.valid_percent)
    axes[1].set_ylabel("Valid area (%)")
    axes[1].set_ylim(-2, 102)
    axes[2].plot(data.date, data.median_cf_cvg)
    axes[2].set_ylabel("Median cf_cvg")
    for ax in axes:
        ax.axvspan(a, b, color="gray", alpha=.18)
        ax.grid(alpha=.25)
    axes[-1].xaxis.set_major_locator(mdates.YearLocator())
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axes[-1].tick_params(axis="x", rotation=35)
    axes[-1].set_xlabel("EOG month")
    fig.suptitle(f"{case.name}: VCMSLCFG monthly 2015–2024; shaded rollout")
    _save_fig(fig, key + "_2015_2024_series")

    # Keep the same-month seasonal trajectories separate from other months.
    fig, (ax, cov) = plt.subplots(2, 1, figsize=(11, 7), sharex=True,
                                constrained_layout=True)
    sub = data[data["year"].isin(case.shown_years)]
    for year, group in sub.groupby("year"):
        ax.plot(group.calendar_month, group.mean_radiance, "-o",
                markersize=3, label=str(year))
        cov.plot(group.calendar_month, group.median_cf_cvg, "-o",
                 markersize=3, label=str(year))
    ax.set_ylabel(f"Mean radiance ({UNIT})")
    ax.set_title("Same-calendar-month comparisons; 2021 snow caution")
    cov.set_ylabel("Median cf_cvg")
    cov.set_xlabel("Calendar month")
    cov.set_xticks(range(1, 13))
    for axis in (ax, cov):
        axis.grid(alpha=.25)
    ax.legend(ncol=3, loc="best")
    _save_fig(fig, key + "_seasonal_comparison")


def plot_contrast(key, case, contrast, year, refs, label, *, primary):
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True,
                             constrained_layout=True)
    x = contrast["month"]
    axes[0].plot(x, contrast["outcome_mean"], "o-", label=f"{year} outcome")
    axes[0].plot(x, contrast["baseline_mean"], "o-",
                 label="Mean of " + "/".join(map(str, refs)))
    axes[0].set_ylabel(f"Monthly radiance ({UNIT})")
    axes[0].legend()
    axes[1].bar(x, contrast["delta_abs"])
    axes[1].axhline(0, color="black", linewidth=.8)
    axes[1].set_ylabel(f"Absolute Δ ({UNIT})")
    axes[2].plot(x, contrast["median_cf_cvg_outcome"], "o-",
                 label=f"{year} cf_cvg")
    axes[2].plot(x, contrast["median_cf_cvg_reference"], "o-",
                 label="Reference cf_cvg")
    axes[2].set_ylabel("Median cf_cvg")
    axes[2].set_xlabel("Calendar month")
    axes[2].set_xticks(x)
    axes[2].legend()
    for ax in axes:
        ax.grid(alpha=.25)
    fig.suptitle(f"{case.name}: {label} — descriptive, NOT causal", fontsize=12)
    _save_fig(fig, f"{key}_{year}_{'primary' if primary else 'secondary'}_contrast")


def run_case(key, case, index, *, allow_nuts3):
    print(f"\n=== {case.name}: {case.mechanism} ===", flush=True)
    geom, label, boundary_file = load_boundary(
        key, case, allow_nuts3=allow_nuts3)
    print(f"Geographic support: {label}", flush=True)
    grids, records = {}, []
    signature = None
    for n, ym in enumerate(MONTHS):
        paths = index[ym]
        grid = grid_month(geom, *paths)
        this_signature = (grid["crs"], grid["transform"], grid["rad"].shape)
        if signature is None:
            signature = this_signature
        elif this_signature != signature:
            raise ValueError(f"{key}: monthly grid/clip changed in {ym} — abort")
        grids[ym] = grid
        records.append(grid_stats(key, case, ym, grid, label, paths))
        if (n + 1) % 24 == 0:
            print(f"  {n+1}/120 months extracted", flush=True)
    panel = pd.DataFrame(records)
    panel.to_csv(OUT / f"{key}_2015_2024_monthly.csv", index=False)
    plot_monthly_series(key, case, panel)
    specs = [("primary", case.treatment_year, case.reference_years,
              case.compare_months, case.contrast_label)]
    if key in SECONDARY:
        y, r, m, comment = SECONDARY[key]
        specs.append(("secondary", y, r, m, comment))
    contrasts = []
    summaries = []
    for kind, year, refs, months, label_event in specs:
        sub = monthly_contrast(panel, year, refs, months)
        sub.insert(0, "case", key)
        sub.insert(1, "kind", kind)
        sub.insert(2, "comment", label_event)
        contrasts.append(sub)
        plot_contrast(key, case, sub, year, refs, label_event,
                      primary=(kind == "primary"))
        observed = float(sub.outcome_mean.mean())
        expected = float(sub.baseline_mean.mean())
        summary = {
            "case": key, "kind": kind, "geographic_support": label,
            "boundary_file": boundary_file, "outcome_year": year,
            "reference_years": ",".join(map(str, refs)),
            "months": ",".join(f"{m:02d}" for m in months),
            "outcome_mean": observed, "reference_mean": expected,
            "period_delta_absolute": observed - expected,
            "period_delta_percent": (
                100 * (observed - expected) / expected if expected > .05 else np.nan),
            "mean_valid_percent_outcome": float(sub.valid_percent_outcome.mean()),
            "mean_cf_cvg_outcome": float(sub.median_cf_cvg_outcome.mean()),
            "mean_cf_cvg_reference": float(sub.median_cf_cvg_reference.mean()),
            "interpretation": label_event + ". " + case.comment,
        }
        if kind == "primary":
            n, median = plot_spatial(key, case, grids, year, refs, months, label_event)
            summary["matched_period_pixels"] = n
            summary["spatial_median_delta"] = median
        summaries.append(summary)
        print(f"  {kind} {year}: ΔR = {observed-expected:+.2f} "
              f"({summary['period_delta_percent']:+.1f}%), "
              f"refs={refs}, months={months}", flush=True)
    pd.concat(contrasts, ignore_index=True).to_csv(
        OUT / f"{key}_calendar_matched_contrasts.csv", index=False)
    return records, summaries


def self_test():
    # Unit-level tests, no real EOG data, no network or local boundary needed.
    source = Path("SVDNB_npp_20240701-20240731_75N060W_vcmslcfg_v10_c202408091300.avg_rade9h.tif")
    assert NAME.match(source.name)["start"] == "20240701"
    assert not NAME.match(source.name.replace("75N060W", "75N060E"))
    def make(x):
        return {"rad": np.array(x, dtype=np.float32)}
    z = {
        "2020-10": make([[2, 3], [2, np.nan]]),
        "2021-10": make([[4, 3], [2, 2]]),
        "2022-10": make([[1, 4], [np.nan, 3]]),
        "2020-11": make([[2, 3], [2, 2]]),
        "2021-11": make([[4, 3], [2, 2]]),
        "2022-11": make([[1, 4], [3, 3]]),
    }
    before, after, diff, monthly, count = paired_pixel_maps(
        z, 2022, (2020, 2021), (10, 11))
    assert np.isclose(diff[0, 0], -2)
    assert np.isclose(diff[0, 1], 1)
    assert np.isnan(diff[1, 0])  # one monthly observation missing => <75% support
    assert np.isnan(diff[1, 1])  # baseline missing in October
    assert count[1, 1] == 1
    assert len(monthly) == 2 and np.isclose(before[0, 0], 3)
    print("Self-test passed (dates, tile, matched pixels, coverage)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", choices=list(CASES),
                        default=list(CASES),
                        help="Subset of five cases; default all")
    parser.add_argument("--inventory", action="store_true",
                        help="Validate file presence without GIS or extraction")
    parser.add_argument("--self-test", action="store_true",
                        help="Run internal synthetic unit tests and quit")
    parser.add_argument("--allow-nuts3-fallback", action="store_true",
                        help="Explicitly permit approximate NUTS3 instead of official LAD")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    index, missing = build_inventory()
    if missing:
        print("MISSING EOG VCMSLCFG 75N060W monthly layers:")
        for ym, layer in missing:
            print(f"  {ym}: {layer}")
        raise RuntimeError(f"{len(missing)} required layers missing; run the "
                           "existing acquisition inventory/importer first.")
    print("Validated file presence: 120/120 monthly VCMSLCFG radiance + cf_cvg pairs")
    if args.inventory:
        return
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    all_records, all_summaries = [], []
    for key in args.cases:
        rows, summaries = run_case(key, CASES[key], index,
                                   allow_nuts3=args.allow_nuts3_fallback)
        all_records.extend(rows)
        all_summaries.extend(summaries)
    pd.DataFrame(all_records).to_csv(OUT / "five_case_monthly_2015_2024.csv",
                                     index=False)
    pd.DataFrame(all_summaries).to_csv(OUT / "five_case_contrast_summary.csv",
                                       index=False)
    print(f"\nCompleted: {len(args.cases)} case(s); "
          f"{len(all_records)} observations, {len(all_summaries)} contrasts")
    print(f"Tables and GeoTIFFs: {OUT}")
    print(f"PNG/PDF figures: {FIG}")
    print("DESCRIPTIVE ONLY: no causal effects, LED spectral correction or "
          "luminous-flux attribution. No locked annual/PRC analysis altered.")


if __name__ == "__main__":
    main()
