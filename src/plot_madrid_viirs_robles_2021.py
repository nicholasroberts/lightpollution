#!/usr/bin/env python3
"""Madrid: locked 2014–24 NUTS3 annual VIIRS, local monthly VNL V1, UCM sky colour.

Core ground photometry is directly transcribed from Robles et al. (2021),
Remote Sensing 13:1511, Table 6, AFTER MIDNIGHT, yearly medians 2011/15/19.
No intermediate UCM measurements are invented. The shaded 2015–16 retrofit
period is broadly documented but NOT a verified per-column switching date.

Sources: https://doi.org/10.3390/rs13081511
         https://doi.org/10.5281/zenodo.4633001 (full raw data)

No numerical mixing of upward VIIRS satellite radiance and zenith sky brightness.
    python src/plot_madrid_viirs_robles_2021.py --self-test
    python src/plot_madrid_viirs_robles_2021.py --inventory
    python src/plot_madrid_viirs_robles_2021.py --check-boundary
    python src/plot_madrid_viirs_robles_2021.py
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import urlencode

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Point, shape

from local_viirs_five_case_analysis import (
    ROOT, BOUNDARIES, build_inventory, grid_month, fetch_json,
)
from export_coventry_echirolles_nuts3_bridge import (
    ANNUAL, ENDPOINT, strict_snow_cache,
)

REGION = "ES300"
TILE = "75N060W"
# Madrid City Council's authoritative municipality polygon, Layer 4.
# This must not be replaced by the much larger ES300 Comunidad de Madrid.
BOUNDARY_SOURCE = (
    "https://services-eu1.arcgis.com/UrtTdLP2PYSxwtI6/ArcGIS/rest/services/"
    "L%C3%ADmites_administrativos_oficiales_del_Ayuntamiento_de_Madrid/"
    "FeatureServer/4/query"
)
BOUNDARY_CACHE = BOUNDARIES / "madrid_municipality_ayto_layer4.geojson"
PROCESSED = ROOT / "data/processed/local_lighting_interventions/madrid_ground_viirs"
FIGURES = ROOT / "figures/local_lighting_interventions/madrid_ground_viirs"
MONTHLY_CACHE = PROCESSED / "madrid_municipality_vcmslcfg_2015_2024.csv"
ANNUAL_CSV = PROCESSED / "madrid_ES300_locked_annual_vnlv2_2014_2024.csv"
GROUND_CSV = PROCESSED / "robles2021_table6_after_midnight_yearly_medians.csv"

# Published Table 6, after-midnight AstMon annual medians in mag/arcsec2 and
# their printed +/- values (NOT labelled as confidence intervals).
TABLE6 = {
    2011: {"B": (19.15, .17), "V": (17.82, .24), "R": (17.11, .25),
           "B-V": (1.32, .20), "V-R": (.78, .08)},
    2015: {"B": (19.27, .18), "V": (18.07, .22), "R": (17.42, .23),
           "B-V": (1.20, .11), "V-R": (.66, .04)},
    2019: {"B": (18.98, .14), "V": (17.95, .16), "R": (17.44, .18),
           "B-V": (1.02, .08), "V-R": (.51, .03)},
}


def ground_table():
    rows = []
    for year, bands in TABLE6.items():
        for band, (mag, spread) in bands.items():
            row = {
                "year": year, "band": band,
                "published_median_mag_per_arcsec2": mag,
                "published_plusminus_spread": spread,
                "photometer": "AstMon UCM zenith",
                "period": "after midnight; clear/moonless",
                "table": "Robles et al. 2021, Table 6",
                "doi": "10.3390/rs13081511",
            }
            if band in ("B", "V", "R"):
                baseline = TABLE6[2015][band][0]
                row["relative_linear_radiance_2015_100"] = (
                    100 * 10 ** (-.4 * (mag - baseline)))
            else:
                row["relative_linear_radiance_2015_100"] = np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def official_city_boundary(*, allow_network=True):
    """Require exactly one official Madrid municipality polygon and validate it."""
    BOUNDARIES.mkdir(parents=True, exist_ok=True)
    if BOUNDARY_CACHE.exists():
        obj = json.loads(BOUNDARY_CACHE.read_text(encoding="utf-8"))
        if (obj.get("type") != "Feature" or
                obj.get("properties", {}).get("source") != BOUNDARY_SOURCE):
            raise ValueError("Cached Madrid boundary has unverified provenance")
        geo = obj
    else:
        if not allow_network:
            raise FileNotFoundError(BOUNDARY_CACHE)
        url = BOUNDARY_SOURCE + "?" + urlencode({
            "where": "1=1", "outFields": "*", "returnGeometry": "true",
            "outSR": 4326, "f": "geojson",
        })
        result = fetch_json(url)
        if result.get("type") != "FeatureCollection" or len(
                result.get("features", [])) != 1:
            raise ValueError(
                "Madrid council municipality layer did not return exactly "
                "one GeoJSON feature. Do NOT substitute ES300.")
        geo = {
            "type": "Feature",
            "properties": {
                "name": "Madrid municipality",
                "source": BOUNDARY_SOURCE,
                "date_of_boundary": "current fixed support",
            },
            "geometry": result["features"][0]["geometry"],
        }
    polygon = shape(geo["geometry"])
    if (polygon.is_empty or not polygon.is_valid or
            polygon.geom_type not in ("Polygon", "MultiPolygon")):
        raise ValueError("Madrid municipality polygon invalid")
    west, south, east, north = polygon.bounds
    if not (-4.1 < west < -3.5 and -4.0 < east < -3.3 and
            40.2 < south < 40.55 and 40.4 < north < 40.7):
        raise ValueError(f"Madrid municipality bounding box implausible: {polygon.bounds}")
    if not polygon.covers(Point(-3.7038, 40.4168)):
        raise ValueError("Madrid municipality fails Puerta del Sol location test")
    area_km2 = float(gpd.GeoSeries([polygon], crs="EPSG:4326")
                     .to_crs("EPSG:3035").area.iloc[0] / 1e6)
    if not 550 < area_km2 < 665:
        raise ValueError(f"Madrid city boundary area {area_km2:.1f} km2 wrong")
    if not BOUNDARY_CACHE.exists():
        BOUNDARY_CACHE.write_text(json.dumps(geo), encoding="utf-8")
    return polygon, area_km2


def frozen_nuts3():
    """Use annual Figure-1-compatible 2021 snow patch, audit and endpoints."""
    if not ANNUAL.is_file() or not ENDPOINT.is_file():
        raise FileNotFoundError("Frozen NUTS3 annual parquet/endpoints missing")
    source = pd.read_parquet(
        ANNUAL, columns=["NUTS_ID", "year", "radiance_mean_corrected"])
    source = source[source.NUTS_ID.eq(REGION)].copy()
    source["year"] = source.year.astype(int)
    if (len(source) != 11 or source.year.nunique() != 11 or
            set(source.year) != set(range(2014, 2025))):
        raise ValueError("Expected exactly 11 locked Madrid ES300 annual rows")
    if not np.isfinite(source.radiance_mean_corrected).all():
        raise ValueError("Nonfinite Madrid annual values")
    snow = strict_snow_cache()
    patch = snow[snow.NUTS_ID.eq(REGION)]
    if len(patch) != 1:
        raise ValueError("Missing ES300 in audited 2021 snow patch")
    idx = source.year.eq(2021)
    original = float(source.loc[idx, "radiance_mean_corrected"].iloc[0])
    expected = float(patch.radiance_2021_original_calibrated.iloc[0])
    if not np.isclose(original, expected, atol=2e-5, rtol=2e-5):
        raise ValueError("2021 annual disagrees with approved Figure 1 snow cache")
    source["radiance_locked"] = source.radiance_mean_corrected.copy()
    source.loc[idx, "radiance_locked"] = float(
        patch.radiance_2021_snow_corrected.iloc[0])
    source["snow_corrected"] = idx
    endpoint = pd.read_csv(ENDPOINT, dtype={"NUTS_ID": str})
    endpoint = endpoint[endpoint.NUTS_ID.eq(REGION)]
    if len(endpoint) != 1:
        raise ValueError("Madrid ES300 absent or duplicated in locked PRC endpoints")
    for year in (2014, 2024):
        obs = float(source.loc[source.year.eq(year), "radiance_locked"].iloc[0])
        ref = float(endpoint[f"radiance_{year}"].iloc[0])
        if abs(obs-ref) > .025:
            raise ValueError(f"{year} Madrid annual differs from PRC endpoint")
    return source.sort_values("year")[
        ["NUTS_ID", "year", "radiance_locked", "snow_corrected"]]


def local_months(geom, index):
    """Original VNL V1 grid, municipality only; cf_cvg>0 source selection."""
    rows = []
    signature = None
    for year in range(2015, 2025):
        for mon in range(1, 13):
            ym = f"{year}-{mon:02d}"
            if ym not in index:
                raise FileNotFoundError(f"Missing VNL V1 European tile month {ym}")
            raster = grid_month(geom, *index[ym])
            sig = (raster["crs"], raster["transform"], raster["rad"].shape)
            if signature is None:
                signature = sig
            elif sig != signature:
                raise ValueError(f"Grid changed in {ym}")
            valid = np.isfinite(raster["rad"])
            if not valid.any():
                raise ValueError(f"{ym}: no valid cloud-free VIIRS pixels")
            obs = raster["rad"][valid]
            cf = raster["cf"][valid]
            rows.append({
                "year": year, "calendar_month": mon,
                "date": pd.Timestamp(f"{ym}-01"),
                "municipality_mean_vnlv1": float(np.mean(obs)),
                "municipality_median_vnlv1": float(np.median(obs)),
                "valid_pixels": int(valid.sum()),
                "inside_pixels": int(raster["inside"].sum()),
                "median_cf_cvg": float(np.median(cf)),
                "mean_cf3": float(np.mean(obs[cf >= 3]))
                            if (cf >= 3).any() else np.nan,
                "mean_cf5": float(np.mean(obs[cf >= 5]))
                            if (cf >= 5).any() else np.nan,
                "source_product": "EOG VNL V1 VCMSLCFG 75N060W",
            })
        print(f"Madrid municipality {year}: 12/12 months", flush=True)
    return pd.DataFrame(rows)


def draw_figure(annual, monthly, ground, out, *, monthly_end=2020):
    """Independent panels: frozen annual, city monthly, ground BVR, colour."""
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        4, 1, figsize=(12.2, 13.3), constrained_layout=True,
        gridspec_kw={"height_ratios": [1, 1.1, 1, .8]})
    a, b, c, d = axes
    a.plot(annual.year, annual.radiance_locked, "-o",
           color="#173d5b", lw=1.9, ms=5)
    a.axvspan(2014.8, 2016.9, alpha=.23, color="#dda963")
    a.set_xticks(range(2014, 2025))
    a.set_xlim(2013.7, 2024.4)
    a.set_ylabel("VIIRS annual radiance\nnW cm$^{-2}$ sr$^{-1}$")
    a.set_title("A  ES300 Madrid region | locked annual VNL V2, 2014–2024")
    a.text(.98, .95,
           "Major HPS-to-LED conversion 2015–2016\n"
           "Includes power and overnight dimming changes",
           transform=a.transAxes, ha="right", va="top", fontsize=8.5)
    part = monthly[monthly.year <= monthly_end]
    b.plot(part.date, part.municipality_mean_vnlv1, "-",
           lw=1.35, color="#486675")
    b.axvspan(pd.Timestamp("2015-01-01"), pd.Timestamp("2017-01-01"),
              color="#dda963", alpha=.23)
    b.set_xlim(pd.Timestamp("2015-01-01"),
               pd.Timestamp(f"{monthly_end+1}-01-01"))
    b.xaxis.set_major_locator(mdates.YearLocator())
    b.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    b.set_ylabel("VIIRS monthly radiance\nnW cm$^{-2}$ sr$^{-1}$")
    b.set_title("B  Madrid municipality | monthly VNL V1 VCMSLCFG")
    b.text(.99, .95, "2015–16: broad transition, not exact commissioning dates",
           transform=b.transAxes, ha="right", va="top", fontsize=8)
    colours = {"B": "#357cbb", "V": "#398668", "R": "#b4554b"}
    for band in ("B", "V", "R"):
        subset = ground[ground.band.eq(band)].sort_values("year")
        c.plot(subset.year, subset.relative_linear_radiance_2015_100,
               "o", label=f"Johnson {band}",
               color=colours[band], lw=1.6, ms=6)
    c.axhline(100, color="#666", alpha=.4, linestyle="--", lw=.8)
    c.axvspan(2014.8, 2016.9, alpha=.23, color="#dda963")
    c.set_xticks([2011, 2015, 2019]); c.set_xlim(2010.5, 2019.6)
    c.set_ylabel("Relative zenith radiance\n2015 = 100 within each band")
    c.set_title("C  UCM AstMon sky brightness AFTER midnight | Robles Table 6")
    c.legend(loc="upper left", ncol=3, fontsize=8)
    c.text(.99, .94,
           "2011 predates VIIRS annual record.\n"
           "Published 3-year anchors only; no inferred intermediate years.",
           transform=c.transAxes, ha="right", va="top", fontsize=8)
    for band, colour in (("B-V", "#5e6cab"), ("V-R", "#99724d")):
        sub = ground[ground.band.eq(band)].sort_values("year")
        d.errorbar(sub.year, sub.published_median_mag_per_arcsec2,
                   yerr=sub.published_plusminus_spread,
                   fmt="o", capsize=3, lw=1.6, color=colour, label=band)
    d.axvspan(2014.8, 2016.9, alpha=.23, color="#dda963")
    d.set_xticks([2011, 2015, 2019]); d.set_xlim(2010.5, 2019.6)
    d.set_ylabel("Colour index (mag)")
    d.set_title("D  UCM AstMon B−V / V−R | lower index = bluer sky")
    d.legend(fontsize=8)
    for ax in axes:
        ax.grid(axis="y", alpha=.2)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
    fig.suptitle(
        "Madrid | VIIRS upward radiance versus independent ground sky colour\n"
        "Different spectral sensitivities, viewing directions and geographies",
        fontsize=13.2, weight="bold", y=1.027)
    fig.text(
        .01, -.025,
        "Ground: Robles et al. (2021), Remote Sensing 13:1511, "
        "Table 6 (after midnight). Published ± values are not assumed CIs. "
        "VIIRS: locked NUTS3 annual / fixed municipality monthly. "
        "No causal separation.",
        fontsize=7.5, color="#555", ha="left")
    stem = f"madrid_ES300_annual_municipal_monthly_UCM_AstMon_table6_to_{monthly_end}"
    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(out / (stem+ext), dpi=190, bbox_inches="tight")
    plt.close(fig)
    return stem


def self_test():
    g = ground_table()
    assert len(g) == 15
    anchor = g[(g.band.eq("B")) & (g.year.eq(2015))].iloc[0]
    post = g[(g.band.eq("B")) & (g.year.eq(2019))].iloc[0]
    assert np.isclose(anchor.relative_linear_radiance_2015_100, 100)
    assert 129 < post.relative_linear_radiance_2015_100 < 132
    annual = pd.DataFrame({
        "year": range(2014, 2025),
        "radiance_locked": np.linspace(10, 7, 11),
    })
    months = pd.period_range("2015-01", "2024-12", freq="M")
    monthly = pd.DataFrame({
        "year": months.year, "date": months.to_timestamp(),
        "municipality_mean_vnlv1": 12 + np.sin(np.arange(120)/12),
    })
    with TemporaryDirectory(prefix="madrid-plot-") as folder:
        name = draw_figure(annual, monthly, g, Path(folder))
        for ext in (".png", ".pdf", ".svg"):
            assert (Path(folder) / (name+ext)).is_file()
    print("Self-test PASS: Table 6 magnitudes, relative-radiance conversion "
          "and PNG/SVG/PDF rendering. No real sources accessed.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--inventory", action="store_true")
    p.add_argument("--check-boundary", action="store_true")
    p.add_argument("--monthly-end", type=int, choices=range(2019, 2025),
                   default=2020, metavar="YEAR")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return
    if args.check_boundary:
        _, area = official_city_boundary()
        print(f"Validated official Madrid municipality: {area:.1f} km2, "
              f"cached {BOUNDARY_CACHE}")
        return
    index, missing = build_inventory()
    if args.inventory:
        print(f"Madrid tile {TILE}, {len(index)}/120 complete VCMSLCFG "
              f"month pairs. Missing individual layers: {len(missing)}.")
        print("Madrid city footprint != wider ES300 Comunidad de Madrid.")
        return
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} source layers missing in {TILE}; run inventory.")
    geom, area = official_city_boundary()
    annual = frozen_nuts3()
    monthly = local_months(geom, index)
    ground = ground_table()
    PROCESSED.mkdir(parents=True, exist_ok=True)
    annual.to_csv(ANNUAL_CSV, index=False)
    monthly.to_csv(MONTHLY_CACHE, index=False)
    ground.to_csv(GROUND_CSV, index=False)
    stem = draw_figure(annual, monthly, ground, FIGURES,
                       monthly_end=args.monthly_end)
    val = ground.loc[
        (ground.band.eq("B")) & (ground.year.eq(2019)),
        "relative_linear_radiance_2015_100"].iloc[0]
    print(f"Ground Table 6 Johnson B 2015→2019: {val-100:+.1f}%")
    print(f"Official Madrid municipal area: {area:.1f} km2")
    print(f"Figure: {FIGURES / (stem+'.png')}")
    print(f"CSV files: {PROCESSED}")
    print("Locked annual, snow cache, monthly rasters, PRC READ ONLY.")


if __name__ == "__main__":
    main()
