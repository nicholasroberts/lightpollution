#!/usr/bin/env python3
"""Revisit the original spider PRC fit with all 18 individually located field sites.

Use local calibrated 2018 and 2019 annual EOG VIIRS VNL average_masked pixels.
All 18 sites now have dates confirmed by the investigator; 3 in 2018
and 15 in 2019. Compare original multiyear VIIRS-website radiances against
our calibrated, field-year EOG pixels. Do not overwrite the frozen transfer.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from extract_prc_site_annual_radiance import extract_points, locate_raster
from fit_prc_observation_year import fit, predict
from project_config import PROJECT_ROOT, configured_path, load_config


def root_path(value):
    p = Path(value)
    return p if p.is_absolute() else PROJECT_ROOT / p


def read_records(path):
    d = pd.read_csv(path)
    required = {
        "sigg_index", "site_id", "source_filename", "location", "latitude",
        "longitude", "field_year", "prc_sigg",
        "original_multiyear_monthly_radiance_mean",
    }
    if missing := required - set(d.columns):
        raise ValueError(f"Missing reference columns: {missing}")
    if len(d) != 18 or not d.sigg_index.is_unique or not d.site_id.is_unique:
        raise ValueError("Expected 18 unique original PRC and site records")
    if not d.sigg_index.tolist() == list(range(1, 19)):
        raise ValueError("Reference is not in original R sigg order")
    years = pd.to_numeric(d.field_year, errors="coerce")
    if years.isna().any() or not years.isin([2018, 2019]).all():
        raise ValueError("Every site must have a confirmed field year of 2018 or 2019")
    if (int(years.eq(2018).sum()), int(years.eq(2019).sum())) != (3, 15):
        raise ValueError("Expected 3 field sites in 2018, 15 in 2019")
    date_years = pd.to_datetime(d.field_date, errors="raise").dt.year
    if not date_years.eq(years).all():
        raise ValueError("field_date and field_year disagree")
    d["field_year"] = years.astype(int)
    for field in ["latitude", "longitude", "prc_sigg",
                  "original_multiyear_monthly_radiance_mean"]:
        d[field] = pd.to_numeric(d[field], errors="raise")
    if not (d.latitude.between(-90, 90) &
            d.longitude.between(-180, 180)).all():
        raise ValueError("Invalid geographic coordinate")
    return d


def produce_comparisons(d, raster_paths, offset):
    site_dicts = d[["site_id", "latitude", "longitude"]].to_dict("records")
    samples = {}
    for year, path in sorted(raster_paths.items()):
        print(f"Reading EOG calibrated pixel radiance for {year}: {path}", flush=True)
        samples[year], _ = extract_points(path, site_dicts, offset)
        d[f"eog_{year}_corrected"] = [
            samples[year][site_id][1] for site_id in d.site_id
        ]
        d[f"eog_{year}_row"] = [
            samples[year][site_id][2] for site_id in d.site_id
        ]
        d[f"eog_{year}_col"] = [
            samples[year][site_id][3] for site_id in d.site_id
        ]
    d["eog_corrected_field_year"] = np.where(
        d.field_year.eq(2018), d.eog_2018_corrected, d.eog_2019_corrected)
    d["date_status"] = "confirmed"
    if not np.isfinite(d.eog_corrected_field_year.to_numpy(float)).all():
        raise ValueError("One or more site-year calibrated EOG values are unavailable")
    return d


def fit_record(name, frame, xname):
    x = frame[xname].to_numpy(float)
    y = frame.prc_sigg.to_numpy(float)
    if not np.all(np.isfinite(x)):
        raise RuntimeError(f"Invalid/absent EOG samples in {name}")
    info = fit(x, y)
    row = {k:v for k,v in info.items() if k != "z"}
    row.update({"model": name, "x_column": xname,
                "distinct_coordinates":len(frame)})
    return row, info


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sites",
        default="data/reference/prc_spider_18_sites_verified.csv")
    ap.add_argument("--raster-2018")
    ap.add_argument("--raster-2019")
    ap.add_argument("--output-dir",
        default="data/processed/prc_spider_18_site_refit_confirmed")
    a = ap.parse_args()
    config = load_config()
    offset = float(config["calibration"]["pixel_offset_nW_cm2_sr"])
    if int(config["calibration"]["apply_from_year"]) > 2018:
        raise ValueError("The 2017+ correction must apply in 2018 and 2019")
    src_root = configured_path(config, "viirs_raw")
    d = read_records(root_path(a.sites))
    rasters = {}
    for year in [2018, 2019]:
        given = getattr(a, f"raster_{year}")
        raster = (root_path(given) if given else
                  locate_raster(src_root, config["viirs"]["raster_pattern"],year))
        if "average_masked" not in raster.name.lower():
            raise ValueError(f"Expected EOG average_masked: {raster}")
        rasters[year] = raster
    d = produce_comparisons(d, rasters, offset)
    out = root_path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out/"spider18_site_eog_corrected_2018_2019.csv", index=False)

    original = "original_multiyear_monthly_radiance_mean"
    updated = "eog_corrected_field_year"
    # Primary comparison fits the same 18 spider measurements twice.
    # Published Table S1's 14 rows are supplementary sensitivity fits.
    core14 = d.loc[d.site_id.isin(
        [f"s{i:02d}" for i in range(1, 15)])].copy()
    if len(core14) != 14:
        raise ValueError("Could not identify all 14 Table S1 sites")
    specs = [
        ("Original multiyear website CSV, all 18", d, original),
        ("Corrected EOG observation year, all 18", d, updated),
        ("Original multiyear website CSV, Table S1 14", core14, original),
        ("Corrected EOG observation year, Table S1 14", core14, updated),
    ]
    model_rows, fits = [], {}
    for name, frame, col in specs:
        summary, info = fit_record(name, frame, col)
        model_rows.append(summary)
        fits[name] = info
    pd.DataFrame(model_rows).to_csv(out/"spider18_refit_comparison.csv",index=False)

    # Field-year assignment is now confirmed for all 18 sites. Do not
    # produce obsolete "unknown-year" sensitivity outputs.
    grid = np.linspace(0,65,500)
    fig, ax = plt.subplots(figsize=(10,6))
    for i, (name, frame, col) in enumerate(specs):
        fitobj = fits[name]
        ax.plot(grid,predict(grid,fitobj["z"]),label=
                f"{name} (R²={fitobj['R2']:.3f})",lw=2 if i in (0,2) else 1.5)
        if i in (0,2):
            ax.scatter(frame[col],frame.prc_sigg,
                       s=16 if i==0 else 32,alpha=.6,marker="x" if i==0 else "o")
    ax.set(xlim=(-1,65), ylim=(.09,.51),
           xlabel="Radiance (nW cm$^{-2}$ sr$^{-1}$)",
           ylabel="Spider PRC",title="Complete 18-site PRC calibration")
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)
    fig.savefig(out/"spider18_site_refit_comparison.pdf",bbox_inches="tight")
    plt.close(fig)

    # An unexpectedly shared physical pixel is flagged, not silently pooled.
    for year in [2018,2019]:
        dup = d[d.duplicated([f"eog_{year}_row",f"eog_{year}_col"],
                             keep=False)]
        if len(dup):
            print(f"{year} sites sharing VIIRS native pixels:",
                  dup[["site_id","location"]].to_string(index=False),flush=True)
    lines = [
        "Original 18 distinct R filename/PRC observations now have site coordinates.",
        "Original x = multi-year monthly CSV mean; published S1 was rounded.",
        "Original OpenSquare/BusinessPark monthly CSVs identical; audit source.",
        "All pixels are sampled independently from 2018 and 2019 EOG rasters.",
        "18 EOG observations individually geolocated, not 14 duplicated locations.",
        "All 18 measurement dates are confirmed: 3 sites in 2018, 15 in 2019.",
        "Four formerly unlisted Table S1 sites now have exact locations and dates.",
        "WClawn: (51.754256,-1.263855), distinct from WCpitch",
        "  (51.755998,-1.267212) retained from Table S1.",
        "Only calibrated EOG annual pixel data used for new x; no NUTS3 means.",
        "The original Open Square csv and Oxford Business Park csv are identical",
        "but their independently sampled EOG annual pixels may differ.",
        "The production PRC transfer function is NOT modified by this stage.",
        "",
    ]
    for info in model_rows:
        lines.append(f"{info['model']}: y0={info['y_0']:.6f}, "
                     f"yf={info['y_f']:.6f}, k={info['k']:.6f}, "
                     f"R2={info['R2']:.6f}, RMSE={info['RMSE']:.6f}")
    (out/"spider18_refit_summary.txt").write_text(
        "\n".join(lines)+"\n",encoding="utf-8")
    print("\n".join(lines),flush=True)
    print("Wrote:",out)


if __name__=="__main__":
    main()
