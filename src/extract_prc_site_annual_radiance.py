#!/usr/bin/env python3
"""Compare published Meah et al. Table S1 VIIRS 2019 values with local EOG 2018/2019.

Point extraction from existing Annual VNL V2.2 average_masked rasters only.
No NUTS3 aggregation, downloads, monthly products or NASA Black Marble.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np
import rasterio

from project_config import PROJECT_ROOT, configured_path, load_config


YEARS = (2018, 2019)


def locate_raster(root: Path, pattern: str, year: int) -> Path:
    candidates = sorted(p for p in root.rglob(pattern)
                        if p.is_file() and re.search(rf"(?<!\d){year}(?!\d)", p.name))
    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected one {year} average_masked GeoTIFF under {root}, "
            f"found {len(candidates)}:\n" +
            "\n".join(str(p) for p in candidates) +
            f"\nIf needed, specify --raster-{year} /path/to/file.tif"
        )
    return candidates[0]


def read_sites(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 14 or len({r["site_id"] for r in rows}) != 14:
        raise ValueError("Expected 14 unique sites from Meah et al. Table S1")
    for s in rows:
        s["latitude"] = float(s["latitude"])
        s["longitude"] = float(s["longitude"])
        s["published_viirs_2019_nw_cm2_sr"] = float(
            s["published_viirs_2019_nw_cm2_sr"]
        )
        if not (-90 <= s["latitude"] <= 90 and -180 <= s["longitude"] <= 180):
            raise ValueError(f"Invalid coordinates at site {s['site_id']}")
    return rows


def extract_points(path: Path, sites: list[dict], offset: float):
    points = {}
    with rasterio.open(path) as ds:
        if ds.count != 1 or ds.crs is None or ds.crs.to_epsg() != 4326:
            raise ValueError(f"Expected single-band EPSG:4326 EOG GeoTIFF: {path}")
        for s in sites:
            # VIIRS is indexed as (longitude, latitude), not (latitude, longitude).
            lon, lat = s["longitude"], s["latitude"]
            if not (ds.bounds.left <= lon < ds.bounds.right and
                    ds.bounds.bottom < lat <= ds.bounds.top):
                raise ValueError(f"{s['site_id']} is outside {path.name}")
            row, col = ds.index(lon, lat)
            sample = next(ds.sample([(lon, lat)], indexes=1, masked=True))
            invalid = bool(np.ma.getmaskarray(sample)[0])
            raw = float(sample[0]) if not invalid else math.nan
            if (not math.isfinite(raw) or
                    (ds.nodata is not None and raw == ds.nodata)):
                raw = math.nan
            # Same 2017+ positive-pixel calibration as apply_radiance_calibration.py.
            corrected = max(raw - offset, 0.0) if math.isfinite(raw) and raw > 0 else raw
            cx, cy = ds.xy(row, col)
            points[s["site_id"]] = (raw, corrected, row, col, cx, cy)
        metadata = (str(path), str(ds.crs), ds.width, ds.height,
                    ds.res[0], ds.res[1], ds.nodata)
    return points, metadata


def difference(a: float, b: float) -> float:
    return a - b if math.isfinite(a) and math.isfinite(b) else math.nan


def display(value: float) -> str:
    return f"{value:.3f}" if math.isfinite(value) else "NA"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", default="data/reference/prc_sites_meah_table_s1.csv")
    parser.add_argument("--raster-2018", type=Path)
    parser.add_argument("--raster-2019", type=Path)
    parser.add_argument("--output-dir",
                        default="data/processed/prc_site_radiance_2018_2019")
    args = parser.parse_args()
    config = load_config()
    offset = float(config["calibration"]["pixel_offset_nW_cm2_sr"])
    from_year = int(config["calibration"]["apply_from_year"])
    if from_year > min(YEARS):
        raise ValueError("Current production pixel-calibration start is later than 2018")
    base = configured_path(config, "viirs_raw")
    resolve = lambda x: Path(x) if Path(x).is_absolute() else PROJECT_ROOT / x
    sites = read_sites(resolve(args.sites))
    rasters = {}
    for year in YEARS:
        given = getattr(args, f"raster_{year}")
        path = resolve(given) if given else locate_raster(
            base, config["viirs"]["raster_pattern"], year
        )
        if "average_masked" not in path.name.lower():
            raise ValueError(f"Not an EOG average_masked raster: {path}")
        print(f"Sampling EOG {year}: {path}", flush=True)
        rasters[year] = extract_points(path, sites, offset)

    rows = []
    for s in sites:
        p18 = rasters[2018][0][s["site_id"]]
        p19 = rasters[2019][0][s["site_id"]]
        pub = s["published_viirs_2019_nw_cm2_sr"]
        year = int(s["field_date"][:4])
        field = p18 if year == 2018 else p19
        delta = difference(p19[0], p18[0])
        rows.append({
            **s,
            "eog_2018_raw": p18[0],
            "eog_2019_raw": p19[0],
            "eog_2018_pixel_corrected": p18[1],
            "eog_2019_pixel_corrected": p19[1],
            "eog_2019_raw_minus_website": difference(p19[0], pub),
            "eog_2019_corrected_minus_website": difference(p19[1], pub),
            "eog_raw_2019_minus_2018": delta,
            "eog_raw_2019_minus_2018_pct":
                100 * delta / p18[0] if math.isfinite(delta) and p18[0] > 0
                else math.nan,
            "eog_field_year_raw": field[0],
            "eog_field_year_corrected": field[1],
            "eog_2018_row": p18[2],
            "eog_2018_col": p18[3],
            "eog_2019_row": p19[2],
            "eog_2019_col": p19[3],
            "eog_2019_pixel_centre_lon": p19[4],
            "eog_2019_pixel_centre_lat": p19[5],
        })
    out = resolve(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv_out = out / "prc_site_annual_radiance_comparison.csv"
    with csv_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Meah et al. Table S1 website VIIRS 2019 values versus EOG annual pixels",
        "Units: nW cm^-2 sr^-1. Website provenance: viirs_2019 layer.",
        f"EOG calibration: positive pixels corrected by -{offset:g} after 2017.",
        "Raw EOG values are the first comparison; these may use different",
        "composite processing from the original lightpollutionmap.info website.",
        "No NUTS3 means; no Black Marble; no monthly or daily products.",
        "Printed Tesco date 21/04/2109 retained, interpreted as 2019-04-21.",
        "",
    ]
    for year in YEARS:
        m = rasters[year][1]
        lines.extend([f"{year}: {m[0]}", f"  CRS={m[1]}, size={m[2]}x{m[3]}, "
                      f"resolution=({m[4]}, {m[5]}), nodata={m[6]}"])
    lines.extend(["", "site | website 2019 | EOG 2018 raw | EOG 2019 raw | "
                  "EOG 2019 corrected | EOG 2019 raw - website | EOG 2019 raw - 2018 raw"])
    for r in rows:
        lines.append(" | ".join([
            r["site_id"] + " " + r["location"],
            f'{r["published_viirs_2019_nw_cm2_sr"]:.2f}',
            display(r["eog_2018_raw"]),
            display(r["eog_2019_raw"]),
            display(r["eog_2019_pixel_corrected"]),
            display(r["eog_2019_raw_minus_website"]),
            display(r["eog_raw_2019_minus_2018"]),
        ]))
    diffs = np.array([r["eog_2019_raw_minus_website"] for r in rows])
    valid = np.isfinite(diffs)
    lines.append(f"\nValid 2019 EOG pixel samples: {valid.sum()}/{len(rows)}")
    if valid.any():
        lines.append(f"Median absolute EOG 2019 versus website difference: "
                     f"{np.median(np.abs(diffs[valid])):.4f} nW cm^-2 sr^-1")
    summary = out / "prc_site_annual_radiance_comparison.txt"
    summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines), flush=True)
    print(f"\nCSV: {csv_out}\nSummary: {summary}")


if __name__ == "__main__":
    main()
