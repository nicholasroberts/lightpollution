#!/usr/bin/env python3
"""LED-dominated counterexamples: matching locked annual + municipal monthly figures.

Designed BEFORE inspecting either case's VIIRS direction of change.

Rome:  ITI43 locked annual 2014-24; municipality EOG monthly 75N060W,
       2016-2019; event-year illustration 2018 vs same months 2017,
       when independent Rome records show near-constant *installed* lumens.
Richmond: UKI75 locked annual 2014-24; borough EOG monthly 75N060W,
       2017-2024; staged illustrative 2022 vs 2017/18 (no exact
       commissioning date, not a validated undimmed control).

Outputs are descriptive: annual NUTS3 VNL V2 and municipal monthly VNL V1
MUST NOT be numerically subtracted, nor interpreted as 'LED-only' effect sizes.
Rome LED rollout involved earlier reduction in installed lumens; Richmond
council explicitly reports nighttime dimming profiles.

Run:
  python src/plot_led_dominated_contrast_figures.py --self-test
  python src/plot_led_dominated_contrast_figures.py --inventory
  python src/plot_led_dominated_contrast_figures.py --annual-only
  python src/plot_led_dominated_contrast_figures.py --cases richmond
  python src/plot_led_dominated_contrast_figures.py --cases rome

Output figures: figures/local_lighting_interventions/led_dominated_contrasts/
Output CSVs: data/processed/local_lighting_interventions/led_dominated_contrasts/
"""
from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import urlencode
from types import SimpleNamespace

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import shape, mapping
from shapely.ops import unary_union

from export_coventry_echirolles_nuts3_bridge import (
    ROOT, ANNUAL, ENDPOINT, strict_snow_cache,
)
from local_viirs_five_case_analysis import (
    ONS2024, BOUNDARIES, RAWDIR, fetch_json, geojson_geometry, grid_month,
)

OUT = ROOT / "figures/local_lighting_interventions/led_dominated_contrasts"
TABLES = ROOT / "data/processed/local_lighting_interventions/led_dominated_contrasts"
UNIT = "nW cm$^{-2}$ sr$^{-1}$"
# Public water authority map service serving polygons of Roma municipi;
# union all valid subdivisions to create Rome *municipality*, not ITI43.
ROME_ADMIN = (
    "https://webgis.abdac.it/server/rest/services/"
    "T05_P5_Cf/MapServer/1/query"
)
SOURCE_NAME = re.compile(
    r"^SVDNB_(?:npp|j01)_(?P<start>20\d{6})-(?P<end>20\d{6})_"
    r"(?P<tile>75N060[WE])_vcmslcfg_v10_.*\."
    r"(?P<layer>avg_rade9h|cf_cvg)\.tif$", re.I
)


@dataclass(frozen=True)
class Candidate:
    key: str
    city: str
    nuts3: str
    region_name: str
    tile: str
    first_year: int
    last_year: int
    reference_years: tuple[int, ...]
    focus_year: int
    rollout_start: str
    rollout_end: str
    markers: tuple[tuple[str, str], ...]
    caution: str


CASES = {
    "richmond": Candidate(
        "richmond", "Richmond upon Thames", "UKI75",
        "Hounslow and Richmond upon Thames", "75N060W",
        2017, 2024, (2017, 2018), 2022,
        "2019-01-01", "2023-01-01",
        (("2019-12-01", "Dec 2019: financing (not start date)"),
         ("2022-12-01", "By end 2022: reported conversion progress")),
        ("LED retrofit with stated unchanged hours; borough also has "
         "dimming profiles and 2019+ rollout was phased. Do not call "
         "this an undimmed control or assert a known 2022 change point."),
    ),
    "rome": Candidate(
        "rome", "Roma Capitale", "ITI43",
        "Roma metropolitan area", "75N060W",
        2016, 2019, (2017,), 2018,
        "2015-01-01", "2020-01-01",
        (("2016-11-01", "Nov 2016: ~51k retrofit completed"),
         ("2018-12-01", "Dec 2018: ~170k retrofit completed"),
         ("2019-12-01", "End 2019: >180k installed")),
        ("Major HPS-to-LED rollout; installed luminous flux 2017→2018 "
         "~+1% while annual estimated energy fell ~27%. Total installed "
         "lumens fell greatly 2014→2018. At-overpass dimming is "
         "unverified; NOT an undimmed or constant-lumen programme."),
    ),
}


def dates(c: Candidate):
    return [f"{y}-{m:02d}" for y in range(c.first_year, c.last_year + 1)
            for m in range(1, 13)]


def raster_index():
    """Match EOG *observation* month, tile and product, never creation stamp."""
    entries = {}
    for path in RAWDIR.rglob("*.tif"):
        m = SOURCE_NAME.match(path.name)
        if not m:
            continue
        ym = m["start"][:4] + "-" + m["start"][4:6]
        key = (m["tile"].upper(), ym, "rad" if m["layer"] == "avg_rade9h"
               else "cf")
        entries.setdefault(key, []).append(path)
    # Original importer may put W and E files in the same per-year/tif dir.
    # Mirror the existing five-case canonical file-preference rules so
    # duplicate older pilot extracts do not invalidate all of Richmond.
    chosen = {}
    for (tile, ym, layer), paths in entries.items():
        canonical = f"/monthly_vcmslcfg_{ym[:4]}/tif/"
        def priority(p):
            normalized = "/" + p.as_posix().lower().lstrip("/")
            return (
                canonical in normalized,
                "/monthly_vcmslcfg_" in normalized,
            )
        paths = sorted(paths, key=priority, reverse=True)
        if len(paths) > 1 and priority(paths[0]) == priority(paths[1]):
            raise ValueError(
                f"Ambiguous equally preferred {tile} {ym} {layer}: "
                + ", ".join(str(p) for p in paths[:2])
            )
        chosen[(tile, ym, layer)] = paths[0]
    return chosen


def availability(c: Candidate, index):
    missing = [
        (ym, layer) for ym in dates(c) for layer in ("rad", "cf")
        if (c.tile, ym, layer) not in index
    ]
    return missing


def official_boundary(c: Candidate):
    """Cache only verified municipality/borough polygons; never use NUTS3."""
    BOUNDARIES.mkdir(parents=True, exist_ok=True)
    if c.key == "richmond":
        code = "E09000027"
        path = BOUNDARIES / f"ons2024_lad_{code}.geojson"
        url = ONS2024 + "?" + urlencode({
            "where": f"LAD24CD='{code}'", "outFields": "*",
            "returnGeometry": "true", "outSR": 4326, "f": "geojson",
        })
        if path.exists():
            obj = json.loads(path.read_text(encoding="utf-8"))
        else:
            obj = fetch_json(url)
        geom = geojson_geometry(obj, code)
        if not path.exists():
            path.write_text(json.dumps(obj), encoding="utf-8")
        return geom, f"ONS 2024 LAD {code}", str(path)

    path = BOUNDARIES / "rome_capitale_union_15_municipi.geojson"
    if path.exists():
        obj = json.loads(path.read_text(encoding="utf-8"))
        geom = geojson_geometry(obj)
    else:
        url = ROME_ADMIN + "?" + urlencode({
            "where": "1=1", "outFields": "*",
            "returnGeometry": "true", "outSR": 4326,
            "f": "geojson",
        })
        obj = fetch_json(url)
        if obj.get("type") != "FeatureCollection":
            raise ValueError("Rome boundary query did not return GeoJSON features")
        parts = obj.get("features", [])
        # Source layer 'Limiti del comune di Roma' depicts the municipal
        # subdistricts. Refuse incomplete/duplicated or NUTS3 proxies.
        if not 12 <= len(parts) <= 16:
            raise ValueError(
                f"Unexpected {len(parts)} Rome municipi; need full municipality")
        geoms = [shape(f["geometry"]) for f in parts]
        if not all(g.is_valid and not g.is_empty for g in geoms):
            raise ValueError("Invalid Rome municipal subdistricts")
        geom = unary_union(geoms)
        record = {"type": "Feature", "properties": {
            "name": "Roma Capitale, union of source municipal subdivisions",
            "source_url": url, "number_subdivisions": len(parts),
        }, "geometry": mapping(geom)}
    if not geom.is_valid or geom.is_empty:
        raise ValueError("Rome municipality boundary invalid")
    # Expected Rome municipality approx 1285 km2; reject mistakenly
    # returning metropolitan NUTS3/too-small one-district masks.
    square_km = float(gpd.GeoSeries(
        [geom], crs="EPSG:4326").to_crs("EPSG:3035").area.iloc[0] / 1e6)
    if not 900 < square_km < 1700:
        raise ValueError(
            f"Rome boundary area {square_km:.1f} km2 impossible for municipality; "
            f"check source {ROME_ADMIN}, do not substitute ITI43")
    if not path.exists():
        path.write_text(json.dumps(record), encoding="utf-8")
    return geom, "Rome Capitale municipality (15 subdistrict union)", str(path)


def frozen_annual():
    """Use ONLY the already locked European VNL V2 annual data and snow audit."""
    if not ANNUAL.exists() or not ENDPOINT.exists():
        raise FileNotFoundError(
            f"Frozen annual inputs absent: {ANNUAL} or {ENDPOINT}")
    snow = strict_snow_cache()
    required = ["NUTS_ID", "year", "radiance_mean_corrected"]
    ids = [c.nuts3 for c in CASES.values()]
    a = pd.read_parquet(ANNUAL, columns=required)
    a = a.loc[a.NUTS_ID.isin(ids)].copy()
    a["year"] = a.year.astype(int)
    if len(a) != 22 or a.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Need 11 unambiguous locked annual years in both NUTS3s")
    if not set(a.year) == set(range(2014, 2025)):
        raise ValueError("Annual 2014–24 years incomplete")
    if not np.isfinite(a.radiance_mean_corrected).all():
        raise ValueError("Nonfinite annual radiance")
    scols = ["NUTS_ID", "radiance_2021_original_calibrated",
             "radiance_2021_snow_corrected"]
    s = snow.loc[snow.NUTS_ID.isin(ids), scols]
    a = a.merge(s, on="NUTS_ID", validate="many_to_one")
    if a[scols[1]].isna().any():
        raise ValueError("Verified snow cache missing one of the NUTS3 regions")
    is2021 = a.year == 2021
    orig = a.loc[is2021, "radiance_mean_corrected"].to_numpy(float)
    if not np.allclose(
        orig,
        a.loc[is2021, "radiance_2021_original_calibrated"].to_numpy(float),
        atol=2e-5, rtol=2e-5,
    ):
        raise ValueError("Annual 2021 parquet differs from locked Figure 1 cache")
    a["radiance_locked"] = a.radiance_mean_corrected.astype(float)
    a.loc[is2021, "radiance_locked"] = a.loc[
        is2021, "radiance_2021_snow_corrected"].to_numpy(float)

    endpoints = pd.read_csv(ENDPOINT, dtype={"NUTS_ID": str})
    needed = ["NUTS_ID", "radiance_2014", "radiance_2024"]
    if any(n not in endpoints.columns for n in needed):
        raise ValueError("Missing locked endpoint radiance fields")
    e = endpoints.loc[endpoints.NUTS_ID.isin(ids), needed].copy()
    if len(e) != 2 or e.NUTS_ID.duplicated().any():
        raise ValueError("Locked endpoint missing/duplicate Rome or Richmond NUTS3")
    a = a.merge(e, on="NUTS_ID", validate="many_to_one")
    for year in (2014, 2024):
        subset = a[a.year.eq(year)]
        error = np.abs(
            subset.radiance_locked.to_numpy(float) -
            subset[f"radiance_{year}"].to_numpy(float)
        )
        if (error > .025).any():
            raise ValueError(f"Locked annual VNL V2 {year} vs PRC figure mismatch")
    return a.sort_values(["NUTS_ID", "year"])


def local_monthly(c: Candidate, index):
    """Recompute municipality means with original V1 raster and cf_cvg support."""
    geom, label, boundary_path = official_boundary(c)
    records = []
    signature = None
    for ym in dates(c):
        paths = [index[(c.tile, ym, layer)] for layer in ("rad", "cf")]
        grid = grid_month(geom, *paths)
        sig = (grid["crs"], grid["transform"], grid["rad"].shape)
        if signature is None:
            signature = sig
        elif sig != signature:
            raise ValueError(f"{c.key}: month-specific grid alignment changed")
        good = np.isfinite(grid["rad"])
        if not good.any():
            raise ValueError(f"{c.key} {ym}: no cloud-free pixels")
        v = grid["rad"][good]
        cf = grid["cf"]
        row = {
            "case": c.key, "nuts3": c.nuts3,
            "geography": label, "boundary_path": boundary_path,
            "tile": c.tile, "product": "VCMSLCFG VNL V1",
            "month": ym, "year": int(ym[:4]),
            "calendar_month": int(ym[5:]),
            "date": pd.Timestamp(ym+"-01"),
            "mean_radiance": float(v.mean()),
            "median_radiance": float(np.median(v)),
            "valid_percent": float(100*good.sum()/grid["inside"].sum()),
            "median_cf_cvg": float(np.median(cf[good])),
            "mean_cf3": (float(grid["rad"][good & (cf>=3)].mean())
                         if np.any(good & (cf>=3)) else np.nan),
            "mean_cf5": (float(grid["rad"][good & (cf>=5)].mean())
                         if np.any(good & (cf>=5)) else np.nan),
        }
        records.append(row)
    return pd.DataFrame(records)


def matching_year(c: Candidate, m):
    refs = (m[m.year.isin(c.reference_years)]
            .groupby("calendar_month").mean_radiance.mean())
    event = m[m.year.eq(c.focus_year)].copy()
    if len(event)!=12 or len(refs)!=12:
        raise ValueError(f"{c.key}: focus year/reference years incomplete")
    event["reference_mean"] = event.calendar_month.map(refs)
    event["anomaly"] = event.mean_radiance-event.reference_mean
    return event


def draw_annual(ax, a, c):
    part = a[a.NUTS_ID.eq(c.nuts3)].sort_values("year")
    yy = part.year.to_numpy(int)
    rr = part.radiance_locked.to_numpy(float)
    ax.plot(yy, rr, "-o", color="#24465f", lw=1.9, ms=4.5,
            label="Locked annual VNL V2 NUTS2016")
    start = int(c.rollout_start[:4])
    end = int(c.rollout_end[:4])
    ax.axvspan(start-.4, min(end-.4, 2024.5),
               color="#90b7d3", alpha=.24)
    ax.set(xlim=(2013.6, 2024.4), xticks=list(range(2014,2025)))
    ax.set_ylabel(f"Regional annual radiance\n({UNIT})")
    ax.grid(axis="y", alpha=.2)
    ax.set_title(f"A   {c.nuts3} | {c.region_name} — 2014–2024")
    ax.text(.015, .08, f"2014 {rr[0]:.2f} → 2024 {rr[-1]:.2f}\n"
            f"Δ {100*(rr[-1]/rr[0]-1):+.1f}%",
            transform=ax.transAxes, fontsize=9,
            va="bottom", color="#24465f")
    ax.legend(loc="upper right", fontsize=7.8)


def render(c: Candidate, a, m, folder: Path, *, formats=("png","pdf","svg")):
    event = matching_year(c,m)
    plt.rcParams.update({
        "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
        "font.size":9.5,"axes.spines.top":False,
        "axes.spines.right":False,
    })
    fig, axes = plt.subplots(3,1,figsize=(12.8,10.1),
                            constrained_layout=True,
                            gridspec_kw={"height_ratios":[1,1.18,.85]})
    aa,ab,ac=axes
    draw_annual(aa,a,c)
    ab.plot(m.date,m.mean_radiance,"-",lw=1.2,color="#969da2",
            label="Municipal monthly VNL V1")
    selected=m[m.year.eq(c.focus_year)]
    ab.plot(selected.date,selected.mean_radiance,"-o",lw=1.8,ms=4,
            color="#193d5a",label=f"{c.focus_year} municipality")
    ab.plot(event.date,event.reference_mean,"--",lw=1.5,
            color="#bf743c",label="Same-month "+"/".join(
                map(str,c.reference_years))+" mean")
    start,end=(pd.Timestamp(d) for d in (c.rollout_start,c.rollout_end))
    ab.axvspan(start,end,color="#90b7d3",alpha=.22)
    ac.axvspan(max(start,pd.Timestamp(f"{c.focus_year}-01-01")),
               min(end,pd.Timestamp(f"{c.focus_year+1}-01-01")),
               color="#90b7d3",alpha=.22)
    for date,label in c.markers:
        p=pd.Timestamp(date)
        ab.axvline(p,color="#456e87",linestyle=":",lw=1,alpha=.85)
        if c.key=="rome" or date=="2019-12-01":
            ab.annotate(label,
                        xy=(mdates.date2num(p),1),
                        xycoords=ab.get_xaxis_transform(),
                        xytext=(3,-10),textcoords="offset points",
                        rotation=90,ha="left",va="top",fontsize=7.2,
                        color="#435e70")
    ab.set_xlim(pd.Timestamp(f"{c.first_year}-01-01"),
                pd.Timestamp(f"{c.last_year+1}-01-01"))
    ab.set_ylabel(f"Municipal monthly radiance\n({UNIT})")
    ab.xaxis.set_major_locator(mdates.YearLocator())
    ab.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ab.set_title(
        f"B   {c.city} municipality | VNL V1 monthly "
        f"{c.first_year}–{c.last_year}\n"
        "       Pale blue: documented multi-year LED rollout, "
        "NOT verified undimmed operation"
    )
    ab.grid(axis="y",alpha=.2)
    ab.legend(fontsize=7.8,loc="upper left",ncol=2)
    ac.plot(event.date,event.anomaly,"-o",lw=1.8,ms=4.3,
            color="#2c6080")
    ac.axhline(0,color="#444",lw=.8)
    ac.xaxis.set_major_locator(mdates.MonthLocator())
    ac.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ac.set_xlim(pd.Timestamp(f"{c.focus_year}-01-01"),
                pd.Timestamp(f"{c.focus_year+1}-01-01"))
    ac.set_ylabel(f"Same-month anomaly\n({UNIT})")
    ac.set_title(
        f"C   {c.focus_year} vs {','.join(map(str,c.reference_years))}"
        "  |  descriptive, not an LED-only effect"
    )
    ac.grid(axis="y",alpha=.2)
    mean_a=float(event.mean_radiance.mean())
    mean_b=float(event.reference_mean.mean())
    ac.text(.98,.92,
            f"Mean ΔR: {mean_a-mean_b:+.2f} "
            f"({100*(mean_a/mean_b-1):+.1f}%)",
            transform=ac.transAxes,ha="right",va="top",
            bbox={"facecolor":"white","edgecolor":"#ddd","alpha":.85})
    fig.suptitle(
        f"{c.city}: LED-dominated conversion and VIIRS history\n"
        "Annual: locked NUTS3 VNL V2 | Monthly: municipality VNL V1 "
        "| NOT interchangeable",
        fontsize=13,weight="bold",y=1.04)
    fig.text(.01,-.035,c.caution+"  Not a control with proven constant "
             "output; annual NUTS3 != municipal monthly VIIRS.",
             fontsize=8,color="#555",ha="left",va="top")
    folder.mkdir(parents=True,exist_ok=True)
    stem=f"{c.key}_led_contrast_annual_and_municipal_monthly_2014_2024"
    for ext in formats:
        fig.savefig(folder/f"{stem}.{ext}",dpi=200,bbox_inches="tight")
    plt.close(fig)
    return event,stem


def render_annual_only(c,a,folder,*,formats=("png","pdf","svg")):
    fig,ax=plt.subplots(figsize=(11,4.1),constrained_layout=True)
    draw_annual(ax,a,c)
    fig.suptitle(
        f"{c.city}: provisional LED-dominated contrast, "
        "locked annual NUTS3 only\n"
        "NOT a municipal LED-only experiment",
        fontsize=11)
    folder.mkdir(parents=True,exist_ok=True)
    stem=f"{c.key}_led_contrast_nuts3_annual_only_2014_2024"
    for ext in formats:
        fig.savefig(folder/f"{stem}.{ext}",dpi=200,bbox_inches="tight")
    plt.close(fig)
    return stem


def self_test():
    # EOG Tile 2 (60W75N) contains both UK and Italy. Reject the prior
    # mistake that assigned Rome to Tile 3 (60E75N, East Asia/Pacific).
    assert CASES["rome"].tile == CASES["richmond"].tile == "75N060W"
    assert len(dates(CASES["rome"])) == 48
    assert len(dates(CASES["richmond"])) == 96
    rg=np.random.default_rng(2026)
    annual=pd.DataFrame(
        [{"NUTS_ID":c.nuts3,"year":y,
          "radiance_locked":(13. if c.key=="rome" else 5.)
          -.12*(y-2014)} for c in CASES.values()
         for y in range(2014,2025)])
    with TemporaryDirectory(prefix="viirs-led-contrasts-") as tmp:
        for c in CASES.values():
            months=pd.period_range(f"{c.first_year}-01",
                                   f"{c.last_year}-12",freq="M")
            m=pd.DataFrame({
                "year":months.year, "calendar_month":months.month,
                "date":months.to_timestamp(),
                "mean_radiance":13+rg.normal(0,1,len(months)),
                "median_cf_cvg":10,
            })
            result,stem=render(c,annual,m,Path(tmp),formats=("png",))
            assert len(result)==12
            assert (Path(tmp)/f"{stem}.png").is_file()
            annualstem=render_annual_only(c,annual,Path(tmp),formats=("png",))
            assert (Path(tmp)/f"{annualstem}.png").is_file()
            assert np.isfinite(result.anomaly).all()
    print("Self-test passed: Rome/Richmond three-panel and annual-only "
          "synthetic plot rendering and seasonal reference.")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--inventory",action="store_true")
    ap.add_argument("--annual-only",action="store_true")
    ap.add_argument("--strict",action="store_true",
                    help="Fail if any requested monthly comparison is incomplete")
    ap.add_argument("--cases",nargs="+",choices=list(CASES),
                    default=list(CASES))
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return
    index=raster_index() if not args.annual_only else {}
    if args.inventory:
        for key in args.cases:
            c=CASES[key]
            missing=availability(c,index)
            by_month=sorted(set(ym for ym,layer in missing))
            n=len(dates(c))
            print(f"{c.city}: {n-len(by_month)}/{n} monthly matched "
                  f"radiance/coverage pairs ready, tile {c.tile}.")
            if by_month:
                print("  Missing observation months:",", ".join(by_month))
        print("Both Rome and Richmond use EOG tile 75N060W (60W75N), "
              "which includes Europe and was already extracted for the "
              "original five-city study. No eastern-tile download required.")
        return
    annual=frozen_annual()
    TABLES.mkdir(parents=True,exist_ok=True)
    OUT.mkdir(parents=True,exist_ok=True)
    annual[["NUTS_ID","year","radiance_locked"]].to_csv(
        TABLES/"rome_richmond_locked_annual_NUTS3_2014_2024.csv",index=False)
    for key in args.cases:
        c=CASES[key]
        if args.annual_only:
            stem=render_annual_only(c,annual,OUT)
            print("Annual-only:",OUT/f"{stem}.png")
            continue
        missing=availability(c,index)
        if missing:
            by_month=sorted(set(ym for ym,_ in missing))
            print(
                f"SKIPPED {key}: {len(by_month)}/{len(dates(c))} "
                f"municipality monthly periods missing for {c.tile}."
                " Run --inventory for dates.")
            if args.strict:
                raise FileNotFoundError(
                    f"{key}: missing {by_month}; tile {c.tile}")
            continue
        monthly=local_monthly(c,index)
        monthly.to_csv(TABLES/f"{key}_municipality_monthly_vnlv1.csv",
                       index=False)
        event,stem=render(c,annual,monthly,OUT)
        event[["year","calendar_month","date","mean_radiance",
               "reference_mean","anomaly","median_cf_cvg",
               "valid_percent","mean_cf3","mean_cf5"]].to_csv(
            TABLES/f"{key}_monthly_same_month_anomaly.csv",index=False)
        print(f"Completed {key}: {OUT/(stem+'.png')}")
    print(f"Sources read-only; outputs {OUT} and {TABLES}.")


if __name__=="__main__":
    main()
