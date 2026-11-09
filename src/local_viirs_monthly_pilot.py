#!/usr/bin/env python3
"""Preliminary municipality-scale EOG VIIRS monthly radiance/coverage pilot.

Portsmouth: UKJ31 NUTS2016 polygon, 2020-2021.
Échirolles: INSEE 38151 official French commune polygon, 2020-2022.

Spatial maps: paired-month radiance differences, multi-month period contrasts,
and per-pixel cloud-free observation counts. Maps use the original EOG grid.

DOES NOT change the locked annual 2014-24 radiance or biological PRC data.
Run from the lightpollution repository root in .venv:
    python src/local_viirs_monthly_pilot.py
"""
from __future__ import annotations

import json
import re
import urllib.request
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask, geometry_window
from rasterio.windows import bounds as window_bounds
from matplotlib.colors import TwoSlopeNorm
from shapely.geometry import mapping, shape

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/viirs"
GEO = ROOT / "data/raw/local_lighting_boundaries"
OUT = ROOT / "data/processed/local_lighting_interventions/monthly_pilot"
FIG = ROOT / "figures/local_lighting_interventions/monthly_pilot"
PRODUCTS = ("vcmslcfg", "vcmcfg")
YM = re.compile(r"(?<!\d)(20\d{2})(0[1-9]|1[0-2])[0-3]\d(?!\d)")
CASES = {"portsmouth": ("2020-01", "2021-12"),
         "echirolles": ("2020-01", "2023-01")}


# Predeclared exploratory contrasts. Portsmouth's within-2020 comparison
# changes season and MUST NOT be interpreted as a policy effect.
SPATIAL_CASES = {
    "echirolles": {
        "periods": (("2021-10", "2021-12", "Oct-Dec 2021"),
                    ("2022-10", "2022-12", "Oct-Dec 2022")),
        "pairs": [(f"2022-{m}", f"2021-{m}") for m in ("09", "10", "11", "12")],
        "description": "Same calendar months; Oct 2022 road switch-off",
    },
    "portsmouth": {
        "periods": (("2020-01", "2020-05", "Jan-May 2020"),
                    ("2020-09", "2020-12", "Sep-Dec 2020")),
        "pairs": [(f"2020-{m}", f"2021-{m}") for m in ("09", "10", "11", "12")],
        "description": "Exploratory within-2020 contrast (seasonally confounded)",
    },
}


def find_months():
    """Pick one correctly labelled radiance and coverage raster per product-month."""
    groups = defaultdict(list)
    for p in INPUT.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in {".tif", ".tiff"}:
            continue
        m = YM.search(p.name)
        if not m:
            continue
        ym = m.group(1) + "-" + m.group(2)
        full = str(p).lower()
        product = ("vcmslcfg" if "vcmslcfg" in full else
                   "vcmcfg" if "vcmcfg" in full else None)
        name = p.name.lower()
        kind = ("rad" if "avg_rade9" in name else
                "cov" if "cf_cvg" in name else None)
        if product and kind:
            groups[(product, ym, kind)].append(p)

    def chosen(paths):
        if not paths:
            return None

        def rank(p):
            full = str(p).lower().replace("\\", "/")
            name = p.name.lower()
            return (4 if "global" in name else
                    3 if "75n060w" in name and "monthly_" in full else
                    2 if "75n060w" in name else 1,
                    0 if "monthly_july_pilot" in full else 1)

        ranks = [rank(p) for p in paths]
        best = max(ranks)
        pool = [p for p, score in zip(paths, ranks) if score == best]
        if len(pool) > 1:
            raise ValueError(f"Ambiguous monthly files: {[str(p) for p in pool]}")
        return pool[0]

    return {(product, ym): (chosen(groups[(product, ym, "rad")]),
                            chosen(groups[(product, ym, "cov")]))
            for product in PRODUCTS
            for ym in pd.period_range("2020-01", "2023-01", freq="M").astype(str)}


def polygons():
    nuts = ROOT / "data/raw/nuts/NUTS_RG_01M_2016_4326_LEVL_3.geojson"
    if not nuts.exists():
        raise FileNotFoundError(f"Missing fixed NUTS2016 geography: {nuts}")
    gdf = gpd.read_file(nuts)
    part = gdf.loc[gdf.NUTS_ID == "UKJ31"]
    if len(part) != 1:
        raise ValueError(f"Expected one Portsmouth UKJ31 polygon, found {len(part)}")
    portsmouth = part.to_crs("EPSG:4326").geometry.iloc[0]

    GEO.mkdir(parents=True, exist_ok=True)
    path = GEO / "echirolles_insee_38151.geojson"
    if not path.exists():
        url = ("https://geo.api.gouv.fr/communes/38151?"
               "fields=nom,code,contour&format=geojson&geometry=contour")
        req = urllib.request.Request(
            url, headers={"User-Agent": "Bristol-university-VIIRS-study/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=35) as response:
                obj = json.load(response)
        except Exception as exc:
            raise RuntimeError(
                "Could not obtain official Echirolles municipality boundary. "
                f"Save an authoritative INSEE 38151 GeoJSON Feature to {path}, "
                f"then rerun. Cause: {exc}") from exc
        if obj.get("type") not in {"Feature", "FeatureCollection"}:
            raise ValueError("French GIS API did not return GeoJSON")
        if obj.get("type") == "FeatureCollection":
            fs = obj["features"]
            if len(fs) != 1:
                raise ValueError("Expected exactly one Echirolles commune feature")
            obj = fs[0]
        path.write_text(json.dumps(obj, ensure_ascii=False), encoding="utf-8")
        print("Saved official Echirolles boundary:", path)
    obj = json.loads(path.read_text(encoding="utf-8"))
    if obj["type"] == "FeatureCollection":
        obj = obj["features"][0]
    echirolles = shape(obj["geometry"])
    if not echirolles.is_valid or echirolles.is_empty:
        raise ValueError("Echirolles GeoJSON has invalid/empty polygon")
    return {"portsmouth": portsmouth, "echirolles": echirolles}


def extract(city, polygon, product, ym, radfile, covfile):
    with rasterio.open(radfile) as rd, rasterio.open(covfile) as cv:
        if rd.crs != cv.crs or rd.transform != cv.transform or rd.shape != cv.shape:
            raise ValueError(f"Different grids for {radfile.name} and {covfile.name}")
        if rd.crs is None:
            raise ValueError(f"Missing CRS in {radfile}")
        geom = gpd.GeoSeries([polygon], crs="EPSG:4326").to_crs(rd.crs).iloc[0]
        try:
            win = geometry_window(rd, [mapping(geom)])
        except Exception as exc:
            raise ValueError(f"{city} outside raster {radfile}: {exc}") from exc
        inside = geometry_mask(
            [mapping(geom)],
            out_shape=(int(win.height), int(win.width)),
            transform=rd.window_transform(win),
            invert=True,
            all_touched=False)
        rad = rd.read(1, window=win, masked=True)
        cov = cv.read(1, window=win, masked=True)
        rr = np.asarray(rad.data, dtype=float)
        cc = np.asarray(cov.data, dtype=float)
        observed = (inside & ~np.ma.getmaskarray(rad) &
                    ~np.ma.getmaskarray(cov) & np.isfinite(rr) &
                    np.isfinite(cc) & (cc > 0))
        possible = int(inside.sum())
        n = int(observed.sum())
        if possible == 0:
            raise ValueError(f"No pixel centres inside {city} boundary for {ym}")
        valid_cf = inside & ~np.ma.getmaskarray(cov) & np.isfinite(cc)
        cc_all = np.where(valid_cf, cc, 0)
        z = rr[observed]
        return dict(
            city=city, product=product, month=ym,
            polygon_pixels=possible, valid_pixels=n,
            valid_percent=100*n/possible,
            mean_radiance=float(np.mean(z)) if n else np.nan,
            median_radiance=float(np.median(z)) if n else np.nan,
            mean_cf_cvg_all_pixels=float(cc_all[inside].mean()),
            median_cf_cvg_observed=float(np.median(cc[observed])) if n else np.nan,
            sum_cf_cvg=int(np.rint(cc_all[inside].sum())),
            rad_file=str(radfile.relative_to(ROOT)),
            cov_file=str(covfile.relative_to(ROOT)))


def plot_city(city, df):
    dat = df.copy().sort_values("date")
    fig, axes = plt.subplots(
        4, 1, figsize=(12, 10), sharex=True, constrained_layout=True)
    variants = [
        ("mean_radiance", "Mean radiance (nW cm$^{-2}$ sr$^{-1}$)"),
        ("valid_percent", "City polygon observed (%)"),
        ("median_cf_cvg_observed",
         "Median cloud-free observations / valid pixel"),
        ("sum_cf_cvg", "Total cloud-free observation counts across city")]
    for product, d in dat.groupby("product", sort=False):
        for ax, (var, _) in zip(axes, variants):
            ax.plot(d.date, d[var], "-o", markersize=3, label=product)
    for ax, (_, ylabel) in zip(axes, variants):
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
        if city == "echirolles":
            ax.axvline(pd.Timestamp("2022-10-01"),
                       color="gray", linestyle="--", linewidth=1)
        else:
            ax.axvspan(pd.Timestamp("2020-06-01"),
                       pd.Timestamp("2020-09-01"),
                       color="gray", alpha=0.14)
    axes[0].legend(loc="best")
    if city == "echirolles":
        fig.suptitle(
            "Échirolles (municipality): 2020–2022 VIIRS, October 2022 switch-off")
    else:
        fig.suptitle(
            "Portsmouth (NUTS2016 UKJ31): monthly VIIRS; "
            "summer 2020 dimming window approximate")
    axes[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    axes[-1].tick_params(axis="x", rotation=35)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{city}_monthly_radiance_and_cloudfree.png", dpi=180)
    fig.savefig(FIG / f"{city}_monthly_radiance_and_cloudfree.pdf")
    plt.close(fig)

    if city == "echirolles":
        fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
        plotdf = dat[dat["product"] == "vcmslcfg"]
        for year, d in plotdf.groupby(plotdf.date.dt.year):
            ax.plot(d.date.dt.month, d.mean_radiance, "-o", label=str(year))
        ax.axvspan(10, 12.5, alpha=0.11, color="gray")
        ax.set(xlim=(0.7, 12.5), xticks=list(range(1, 13)),
               xlabel="Calendar month",
               ylabel="Mean radiance (nW cm$^{-2}$ sr$^{-1}$)",
               title="Échirolles: 2020–2022 same-month reference "
                     "(October 2022 switch-off)")
        ax.grid(alpha=0.25)
        ax.legend()
        fig.savefig(FIG / "echirolles_same_calendar_month_comparison.png",
                    dpi=180)
        fig.savefig(FIG / "echirolles_same_calendar_month_comparison.pdf")
        plt.close(fig)



def _read_spatial(city, polygon, radfile, covfile):
    """Read native monthly grid pixels, clipping to municipality polygon centres."""
    with rasterio.open(radfile) as rd, rasterio.open(covfile) as cv:
        if rd.crs != cv.crs or rd.transform != cv.transform or rd.shape != cv.shape:
            raise ValueError(f"Raster grid mismatch: {radfile} vs {covfile}")
        geom = gpd.GeoSeries([polygon], crs="EPSG:4326").to_crs(rd.crs).iloc[0]
        window = geometry_window(rd, [mapping(geom)])
        wt = rd.window_transform(window)
        shape_xy = (int(window.height), int(window.width))
        inside = geometry_mask([mapping(geom)], out_shape=shape_xy,
                               transform=wt, invert=True, all_touched=False)
        rad = rd.read(1, window=window, masked=True)
        cf = cv.read(1, window=window, masked=True)
        valid = (inside & ~np.ma.getmaskarray(rad) &
                 ~np.ma.getmaskarray(cf) &
                 np.isfinite(np.asarray(rad.data)) &
                 np.isfinite(np.asarray(cf.data)) &
                 (np.asarray(cf.data) > 0))
        r = np.full(shape_xy, np.nan, dtype=np.float32)
        c = np.full(shape_xy, np.nan, dtype=np.float32)
        r[valid] = np.asarray(rad.data)[valid]
        c[inside & ~np.ma.getmaskarray(cf) &
          np.isfinite(np.asarray(cf.data))] = np.asarray(cf.data)[
              inside & ~np.ma.getmaskarray(cf) &
              np.isfinite(np.asarray(cf.data))]
        # A zero cloud-free count is a valid zero for coverage, not radiance.
        return {"rad": r, "cf": c, "inside": inside,
                "transform": wt, "crs": rd.crs, "boundary": geom,
                "extent": window_bounds(window, rd.transform)}


def _polygon_outline(ax, geom):
    parts = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
    for part in parts:
        if not hasattr(part, "exterior"):
            continue
        x, y = part.exterior.xy
        ax.plot(x, y, color="black", linewidth=0.9)


def _map(ax, a, reference, title, *, cmap="viridis",
         norm=None, vmin=None, vmax=None):
    l, b, r, t = reference["extent"]
    im = ax.imshow(a, extent=(l, r, b, t), origin="upper",
                   interpolation="nearest", cmap=cmap,
                   norm=norm, vmin=vmin, vmax=vmax)
    _polygon_outline(ax, reference["boundary"])
    ax.set_title(title, fontsize=10)
    ax.tick_params(labelsize=7)
    ax.set_xlabel("Raster x-coordinate", fontsize=8)
    ax.set_ylabel("Raster y-coordinate", fontsize=8)
    return im


def _mean_grids(grid_by_month, start, end, key):
    months = pd.period_range(start, end, freq="M").astype(str)
    arrays = [grid_by_month[m][key] for m in months
              if m in grid_by_month]
    if len(arrays) != len(months):
        missing = [m for m in months if m not in grid_by_month]
        raise ValueError(f"Cannot map {start} to {end}; missing: {missing}")
    stack = np.stack(arrays)
    good_count = np.sum(np.isfinite(stack), axis=0)
    summed = np.nansum(stack, axis=0)
    out = np.full(stack.shape[1:], np.nan, dtype=np.float32)
    np.divide(summed, good_count, out=out, where=good_count > 0)
    out[good_count == 0] = np.nan
    return out


def _save_native_geotiff(path, raster, reference):
    """Spatial change raster in original source CRS, with missing pixels as NaN."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", driver="GTiff", height=raster.shape[0],
                       width=raster.shape[1], count=1, dtype="float32",
                       crs=reference["crs"], transform=reference["transform"],
                       nodata=np.nan, compress="deflate") as dst:
        dst.write(np.asarray(raster, dtype="float32"), 1)


def plot_spatial_maps(city, polygon, index):
    """Plot city-month source-grid radiance, differences and observation coverage.

    Primary mapping source is VCMSLCFG because VCMCFG loses most summer months.
    Before/after means are compared on their **common valid pixel support**.
    Each monthly difference is also strictly matched pixel-to-pixel.
    """
    spec = SPATIAL_CASES[city]
    needed = set()
    for start, end, _ in spec["periods"]:
        needed.update(pd.period_range(start, end, freq="M").astype(str))
    for after, before in spec["pairs"]:
        needed.update((after, before))
    grids = {}
    for ym in sorted(needed):
        rad, cov = index.get(("vcmslcfg", ym), (None, None))
        if rad is None or cov is None:
            print(f"  MAP missing input for {city}: {ym}")
            continue
        grids[ym] = _read_spatial(city, polygon, rad, cov)
    if len(grids) != len(needed):
        print(f"  MAPS SKIPPED for {city}: missing some required months")
        return

    reference = next(iter(grids.values()))
    for ym, grid in grids.items():
        if (grid["rad"].shape != reference["rad"].shape or
            grid["crs"] != reference["crs"] or
            grid["transform"] != reference["transform"]):
            raise ValueError(f"{city}: nonmatching spatial grid in {ym}")

    (a0, a1, label0), (b0, b1, label1) = spec["periods"]
    before = _mean_grids(grids, a0, a1, "rad")
    after = _mean_grids(grids, b0, b1, "rad")
    shared = np.isfinite(before) & np.isfinite(after)
    before[~shared], after[~shared] = np.nan, np.nan
    delta = after - before
    extent = np.r_[before[np.isfinite(before)], after[np.isfinite(after)]]
    if not extent.size:
        raise ValueError(f"{city}: no common spatial radiance data")
    vmax = float(np.nanpercentile(extent, 99.5))
    scale = float(np.nanmax(np.abs(delta)))
    scale = max(scale, 1e-6)
    norm = TwoSlopeNorm(vmin=-scale, vcenter=0, vmax=scale)

    FIG.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5),
                             constrained_layout=True)
    _map(axes[0], before, reference, label0, vmin=0, vmax=vmax)
    m1 = _map(axes[1], after, reference, label1, vmin=0, vmax=vmax)
    m2 = _map(axes[2], delta, reference, label1 + " minus " + label0,
              cmap="coolwarm", norm=norm)
    fig.colorbar(m1, ax=axes[:2], shrink=0.75,
                 label="VIIRS nW cm$^{-2}$ sr$^{-1}$")
    fig.colorbar(m2, ax=axes[2], shrink=0.75, label="Δ radiance")
    fig.suptitle(f"{city.title()}: spatial radiance — {spec['description']}")
    name = f"{city}_vcmslcfg_period_mean_maps"
    for ext in (".png", ".pdf"):
        fig.savefig(FIG / (name + ext), dpi=200)
    plt.close(fig)
    _save_native_geotiff(OUT / f"{city}_vcmslcfg_period_difference.tif",
                         delta, reference)

    diffs = []
    for after_ym, before_ym in spec["pairs"]:
        diff = grids[after_ym]["rad"] - grids[before_ym]["rad"]
        diffs.append((after_ym, before_ym, diff))
    dscale = float(np.nanmax(np.abs(np.stack([a for _, _, a in diffs]))))
    dscale = max(dscale, 1e-6)
    dnorm = TwoSlopeNorm(vmin=-dscale, vcenter=0, vmax=dscale)
    fig, axes = plt.subplots(2, 2, figsize=(11, 10),
                             constrained_layout=True)
    for ax, (after_ym, before_ym, diff) in zip(axes.flat, diffs):
        im = _map(ax, diff, reference, after_ym + " minus " + before_ym,
                  cmap="coolwarm", norm=dnorm)
        ax.text(0.02, 0.02, f"Median Δ: {np.nanmedian(diff):.2f}",
                transform=ax.transAxes, fontsize=9, va="bottom",
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
    fig.colorbar(im, ax=axes, shrink=0.75,
                 label="VIIRS Δ nW cm$^{-2}$ sr$^{-1}$")
    fig.suptitle(f"{city.title()}: same-calendar-month radiance contrasts")
    name = f"{city}_vcmslcfg_matched_month_difference_maps"
    for ext in (".png", ".pdf"):
        fig.savefig(FIG / (name + ext), dpi=200)
    plt.close(fig)

    # Observation counts for each multi-month period, on the *same* grid.
    # The total is summed (not averaged): expected number of monthly looks.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5),
                             constrained_layout=True)
    cover = []
    for start, end, label in spec["periods"]:
        seq = [grids[m]["cf"] for m in pd.period_range(
            start, end, freq="M").astype(str)]
        cover.append((label, np.nansum(np.stack(seq), axis=0)))
    covmax = float(max(np.nanmax(a) for _, a in cover))
    for ax, (label, a) in zip(axes, cover):
        m = _map(ax, np.where(reference["inside"], a, np.nan),
                 reference, label, vmin=0, vmax=covmax, cmap="cividis")
    fig.colorbar(m, ax=axes, shrink=0.75,
                 label="Total cloud-free observations across period")
    fig.suptitle(f"{city.title()}: spatial cloud-free observation coverage")
    name = f"{city}_vcmslcfg_period_cloudfree_maps"
    for ext in (".png", ".pdf"):
        fig.savefig(FIG / (name + ext), dpi=200)
    plt.close(fig)
    print(f"  MAPS created for {city}: period radiance, matched monthly "
          "differences, cloud-free observation counts, native delta GeoTIFF")


def main():
    index = find_months()
    polygons_by_city = polygons()
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for city, (begin, end) in CASES.items():
        print(f"\n{city.upper()} ({begin} to {end})")
        for product in PRODUCTS:
            done, missing = [], []
            for ym in pd.period_range(begin, end, freq="M").astype(str):
                rad, cov = index.get((product, ym), (None, None))
                if rad is None or cov is None:
                    missing.append(ym)
                    continue
                item = extract(
                    city, polygons_by_city[city], product, ym, rad, cov)
                all_rows.append(item)
                done.append(ym)
            print(f"  {product}: {len(done)} usable month-pairs; "
                  f"missing {', '.join(missing) or 'none'}")
        city_rows = [r for r in all_rows if r["city"] == city]
        if not city_rows:
            print(f"  WARNING: no months available for {city}")
            continue
        df = pd.DataFrame(city_rows)
        df["date"] = pd.to_datetime(df.month + "-01")
        df.to_csv(
            OUT / f"{city}_monthly_radiance_cloudfree.csv", index=False)
        plot_city(city, df)
        plot_spatial_maps(city, polygons_by_city[city], index)
        for p, group in df.groupby("product"):
            print(
                f"  {p}: geometric cells "
                f"{sorted(group.polygon_pixels.unique().tolist())}; "
                f"observed area median {group.valid_percent.median():.1f}%, "
                f"min {group.valid_percent.min():.1f}%; "
                f"median cf_cvg {group.median_cf_cvg_observed.median():.1f}")
    pd.DataFrame(all_rows).to_csv(
        OUT / "combined_monthly_radiance_cloudfree.csv", index=False)
    print(f"\nCSV: {OUT}\nPlots: {FIG}")
    print("NOTE: descriptive plots only; no causal effect inferred "
          "and 2021 snow is uncorrected here.")


if __name__ == "__main__":
    main()
