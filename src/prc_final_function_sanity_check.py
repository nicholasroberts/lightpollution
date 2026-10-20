#!/usr/bin/env python3
"""ONE final PRC calibration sanity check, then freeze the confirmed 18-site fit.

A. Inspect Slabhouse 2018/2019 raw, background-corrected annual EOG radiance,
   and cloud-free coverage when the companion cf_cvg raster is available.
B. Compare 2014->2024 area-weighted NUTS2016 PRC endpoints under the current
   provisional 16-observation vs confirmed 18-observation response functions.
All regions use identical VIIRS raster pixels, calibration offset and geography.
No 11-year atlas regeneration, no maps, no parameter optimization or overwrites.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

from extract_nuts3_radiance import discover_rasters
from historical_prc_pixel_nuts3 import extract_year, load_inputs, load_nuts
from project_config import PROJECT_ROOT, configured_path, load_config

OLD = (0.144845, 0.461175, 0.106483)  # yf, y0, k; current config checked below
FINAL = (0.1452426021, 0.4558190410, 0.1035007080)
SLABHOUSE_LAT_LON = (51.230935, -2.587688)
YEARS = (2014, 2024)


def companion_cf(path, root):
    name = path.name.replace(".average_masked.dat.tif", ".cf_cvg.dat.tif")
    if name == path.name:
        return None
    exact = list(root.rglob(name))
    if len(exact) == 1:
        return exact[0]
    possible = [p for p in root.rglob("*.tif") if
                f"_{2019 if '_2019_' in path.name else 2018}_global_" in p.name
                and ".cf_cvg." in p.name]
    return possible[0] if len(possible) == 1 else None


def slabhouse_check(rasters, config):
    raw_root = configured_path(config, "viirs_raw")
    offset = float(config["calibration"]["pixel_offset_nW_cm2_sr"])
    rows = []
    for year in (2018, 2019):
        yearpath = discover_rasters(config, [year])[year]
        lon, lat = SLABHOUSE_LAT_LON[1], SLABHOUSE_LAT_LON[0]
        with rasterio.open(yearpath) as ds:
            values = next(ds.sample([(lon, lat)], indexes=1, masked=True))
            raw = (float(values[0]) if not np.ma.getmaskarray(values)[0] else np.nan)
            corrected = (max(0.0, raw-offset) if np.isfinite(raw) and raw>0
                         else raw)
            pixelrow, pixelcol = ds.index(lon, lat)
        cfpath = companion_cf(yearpath, raw_root)
        coverage = np.nan
        cfpath_name = "(not present: coverage unverified)"
        if cfpath:
            with rasterio.open(cfpath) as cf:
                if cf.crs is None or cf.crs.to_epsg() != 4326:
                    raise ValueError(f"Cloud-free layer has unexpected CRS: {cfpath}")
                val = next(cf.sample([(lon, lat)], indexes=1, masked=True))
                if not np.ma.getmaskarray(val)[0]:
                    coverage = float(val[0])
                cfpath_name = cfpath.name
        rows.append({"year":year,"site":"Slabhouse, Wells",
                     "latitude":lat,"longitude":lon,
                     "eog_original_radiance":raw,
                     "eog_calibrated_radiance":corrected,
                     "annual_cloudfree_observation_count":coverage,
                     "viirs_row":pixelrow,"viirs_col":pixelcol,
                     "radiance_file":yearpath.name,"cf_cvg_file":cfpath_name})
    return pd.DataFrame(rows)


def compare_endpoints(config, rasters, min_coverage):
    nuts = load_nuts(config)
    _, current, offset, first = load_inputs(
        config, PROJECT_ROOT/"config/prc_transfer.yaml")
    if not np.allclose(current,OLD,atol=5e-7,rtol=0):
        raise RuntimeError("Existing provisional PRC parameters changed. "
                           "Do not compare against an unknown baseline.")
    pieces=[]
    for label,pars in [("provisional16",current),("confirmed18",FINAL)]:
        for year in YEARS:
            d=extract_year(rasters[year],year,nuts,config,pars,offset,first)
            d["function"]=label
            pieces.append(d)
    long=pd.concat(pieces,ignore_index=True)
    key=["NUTS_ID","year"]
    a=long[long.function.eq("provisional16")].drop(columns="function")
    b=long[long.function.eq("confirmed18")].drop(columns="function")
    pair=a.merge(b,on=key,validate="one_to_one",
                 suffixes=("_provisional16","_confirmed18"))
    da=np.abs(pair.valid_area_km2_provisional16-pair.valid_area_km2_confirmed18)
    if np.any(da > 1e-5):
        raise AssertionError("Both transfer functions must use IDENTICAL pixels")
    p14=pair[pair.year.eq(2014)].set_index("NUTS_ID")
    p24=pair[pair.year.eq(2024)].set_index("NUTS_ID")
    first_cols=["prc_mean_provisional16","prc_mean_confirmed18",
                "valid_area_km2_provisional16","valid_area_km2_confirmed18"]
    d=p14[first_cols].join(p24[first_cols],how="inner",
                            lsuffix="_2014",rsuffix="_2024")
    for suffix in ("provisional16","confirmed18"):
        d["delta_"+suffix]=(d["prc_mean_"+suffix+"_2024"]-
                             d["prc_mean_"+suffix+"_2014"])
        d["pct_"+suffix]=100*d["delta_"+suffix]/d["prc_mean_"+suffix+"_2014"]
    d["delta_difference"]=d.delta_confirmed18-d.delta_provisional16
    d["percent_difference"]=d.pct_confirmed18-d.pct_provisional16
    d["coverage_stable_endpoint"]=(
        (d.valid_area_km2_provisional16_2014>0)
        & (d.valid_area_km2_provisional16_2024>0)
        & (d.valid_area_km2_provisional16_2024 /
           d.valid_area_km2_provisional16_2014 >= min_coverage)
        & (d.valid_area_km2_provisional16_2024 /
           d.valid_area_km2_provisional16_2014 <= 1/min_coverage))
    all_year_panel=PROJECT_ROOT/(
        "data/processed/historical_prc_2014_2024/"
        "prc_nuts3_annual_cumulative.csv")
    if all_year_panel.exists():
        p=pd.read_csv(all_year_panel)
        if "coverage_stable" in p:
            v=p[p.year.eq(2024)].set_index("NUTS_ID").coverage_stable
            d["coverage_stable_full_11years"]=d.index.map(v).fillna(False).astype(bool)
            mask=d.coverage_stable_endpoint & d.coverage_stable_full_11years
        else:
            mask=d.coverage_stable_endpoint
    else:
        mask=d.coverage_stable_endpoint
    eligible=d[mask].dropna(subset=["delta_provisional16","delta_confirmed18"])
    if eligible.empty:
        raise RuntimeError("No comparable regions with stable pixel coverage")
    for name in ["provisional16","confirmed18"]:
        d["gain_rank_"+name]=d["delta_"+name].rank(ascending=False,method="min")
    # Europe-wide endpoint means use same regional valid areas as the original
    # historical PRC atlas. The two calibrations share identical raster masks.
    europe={}
    for year in YEARS:
        x=pair[pair.year.eq(year)]
        w=x.valid_area_km2_provisional16.to_numpy(float)
        valid=(np.isfinite(w)&(w>0)&np.isfinite(x.prc_mean_provisional16)
               &np.isfinite(x.prc_mean_confirmed18))
        europe[year]={name:float(np.average(x.loc[valid,"prc_mean_"+name],
                                             weights=w[valid]))
                      for name in ["provisional16","confirmed18"]}
    delta_e={name:europe[2024][name]-europe[2014][name]
             for name in ["provisional16","confirmed18"]}
    pct_e={name:100*delta_e[name]/europe[2014][name]
           for name in delta_e}
    a_top=set(eligible.nlargest(20,"delta_provisional16").index)
    b_top=set(eligible.nlargest(20,"delta_confirmed18").index)
    spear=eligible["delta_provisional16"].corr(eligible["delta_confirmed18"],
                                                method="spearman")
    info=[
        "FINAL PRC-SENSITIVITY SANITY CHECK — same native pixels for both functions",
        f"NUTS3 regions with stable endpoint/available full-series coverage: {len(eligible)}",
        f"2014 Europe mean PRC: provisional={europe[2014]['provisional16']:.7f} "
        f"confirmed={europe[2014]['confirmed18']:.7f}",
        f"2024 Europe mean PRC: provisional={europe[2024]['provisional16']:.7f} "
        f"confirmed={europe[2024]['confirmed18']:.7f}",
        f"2014->2024 Europe absolute PRC change: "
        f"provisional={delta_e['provisional16']:+.8f} "
        f"confirmed={delta_e['confirmed18']:+.8f}",
        f"2014->2024 Europe percent PRC change: "
        f"provisional={pct_e['provisional16']:+.5f}% "
        f"confirmed={pct_e['confirmed18']:+.5f}%",
        f"Region absolute-delta Spearman rank correlation={spear:.7f}",
        f"Top-20 region PRC gain overlap={len(a_top & b_top)}/20",
        f"Median absolute difference in region 2014->2024 delta="
        f"{np.median(np.abs(eligible.delta_difference)):.9f}",
        f"Maximum absolute difference in region 2014->2024 delta="
        f"{np.max(np.abs(eligible.delta_difference)):.9f}",
        "This is a two-endpoint sensitivity, not a new 11-year atlas.",
        "No 2021 corrections or calibration coefficient optimisation applied.",
        "Decision: use the confirmed full-18 field-year PRC function after review.",
    ]
    return d.reset_index(), info


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",default=
        "data/processed/prc_18site_final_sanity_check")
    p.add_argument("--min-coverage",type=float,default=.95)
    p.add_argument("--slabhouse-only",action="store_true",
                   help="Quick single-point EOG/cf_cvg check without NUTS extraction")
    args=p.parse_args()
    config=load_config()
    out=PROJECT_ROOT/args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    slab=slabhouse_check({},config)
    slab.to_csv(out/"slabhouse_annual_radiance_and_coverage.csv",index=False)
    txt=[
        "Slabhouse two-year point audit",
        slab[["year","eog_original_radiance","eog_calibrated_radiance",
              "annual_cloudfree_observation_count","cf_cvg_file"]].to_string(index=False),
        "If annual cf_cvg is not found, absence is reported rather than guessed.",
        "",
    ]
    if args.slabhouse_only:
        (out/"final_prc_sanity_summary.txt").write_text("\n".join(txt)+"\n")
        print("\n".join(txt))
        print("Outputs:",out)
        return
    rasters=discover_rasters(config,list(YEARS))
    region, lines=compare_endpoints(config,rasters,args.min_coverage)
    region.to_csv(out/"prc_2014_2024_function_sensitivity_nuts3.csv",index=False)
    txt+=lines
    (out/"final_prc_sanity_summary.txt").write_text("\n".join(txt)+"\n")
    print("\n".join(txt))
    print("Outputs:",out)


if __name__=="__main__":
    main()
