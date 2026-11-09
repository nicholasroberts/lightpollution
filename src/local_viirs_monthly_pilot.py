#!/usr/bin/env python3
"""Preliminary municipality-scale EOG VIIRS monthly radiance/coverage pilot.

Portsmouth: UKJ31 NUTS2016 polygon, 2020-2021.
Échirolles: INSEE 38151 official French commune polygon, 2021-Jan 2023.

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
