#!/usr/bin/env python3
"""Versioned UK-inclusive 2014–2020 Jan1 population / VIIRS / PRC cohort.

Extend—not overwrite—the existing frozen 1,273-region matched scientific
cohort with the 179 UK NUTS2016 regions, using:
  2014–2019: existing Eurostat Jan1 total counts, demo_r_pjanaggr3 JSON
  2020: official ONS 1-January 2020 Unified Demographic Collection workbook
  all years: frozen 2018–2020 NUTS2016 land area, calibrated VIIRS annual
             NUTS3 radiance, independently locked pixel-derived PRC means.

ONS Excel P02: C = both-sex total, D = male total, E = female total.
Require C == D+E for *every* NUTS3 region and exact 179 region IDs.
The site-specific ONS 2020 workbook is left untouched. The old 1,273 data
and all PRC/VIIRS source tables remain immutable. A fresh source-annotated
1,452-row cohort, 179-region UK annual panel and a new descriptive country
Spearman table are written under a distinct directory.

This is an explicitly versioned population-coverage extension, not a new
receptor calibration or a correction to any sensor radiance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from historical_prc_pixel_nuts3 import load_nuts
from inspect_ons_uk_jan1_2020_workbook import (
    FILE as ONS_XLSX, ONS_URL, sheet_rows, sheets, strings
)
from project_config import PROJECT_ROOT as ROOT, load_config
from reconstruct_population_density_jan1_nuts2016 import fixed_land_areas
from validate_population_fallback import prepare_population, prepare_area

YEARS = tuple(range(2014, 2021))
ENDPOINTS = (2014, 2019, 2020)
MATCHED_OLD = 1273
N_UK = 179
MATCHED_NEW = MATCHED_OLD + N_UK
OLD = (ROOT / "data/processed/jan1_population_radiance_prc_country_nl333"
       / "matched_region_change_2014_2020.csv")
RAW = ROOT / "data/raw/eurostat"
PRC = {
    t: ROOT / "data/processed" /
       f"historical_prc_2014_2024_{t}_snow_corrected" /
       "prc_nuts3_annual_cumulative.csv"
    for t in ("spider", "moth")
}
RADIANCE = ROOT / "data/processed/fixed_nuts2016_viirs_2014_2024.parquet"
OUT = ROOT / "data/processed/publication_figures/jan1_population_prc_with_uk_2014_2020"


def rootpath(arg):
    p = Path(arg)
    return p if p.is_absolute() else ROOT / p


def year_region_validation(frame, ids, label):
    if frame.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError(f"{label}: duplicate region/year")
    actual = set(zip(frame.NUTS_ID, frame.year))
    expected = {(i, y) for i in ids for y in YEARS}
    if actual != expected:
        missing, extra = expected - actual, actual - expected
        raise ValueError(
            f"{label}: expected {len(expected)} region-years; "
            f"missing={sorted(missing)[:7]}, extra={sorted(extra)[:7]}"
        )


def _as_int(value, label):
    try:
        n = int(str(value).strip())
    except (ValueError, TypeError) as exc:
        raise ValueError(f"ONS {label} not an integer: {value!r}") from exc
    if n <= 0:
        raise ValueError(f"ONS {label} is not a positive headcount: {n}")
    return n


def read_ons_2020_p02(path, ids):
    if not path.is_file():
        raise FileNotFoundError(f"Official ONS 1-Jan-2020 workbook absent: {path}")
    with zipfile.ZipFile(path) as archive:
        ss = strings(archive)
        table = dict(sheets(archive))
        if "P02" not in table:
            raise ValueError("ONS workbook lacks required NUTS3 P02 worksheet")
        items = list(sheet_rows(archive, table["P02"], ss))
    if not any("Population by age, sex and region (NUTS 3)" in
               str(value) for _, cells in items[:8] for value in cells.values()):
        raise ValueError("The ONS P02 workbook has an unexpected NUTS3 table heading")

    records = []
    for row_number, cells in items:
        code = str(cells.get(1, "")).strip().upper()
        if code not in ids:
            continue
        # Verify this is the all-ages BOTH-SEXES value by direct equality
        # with male and female totals in neighbouring P02 columns.
        total = _as_int(cells.get(3), f"{code} P02 C{row_number}")
        males = _as_int(cells.get(4), f"{code} P02 D{row_number}")
        females = _as_int(cells.get(5), f"{code} P02 E{row_number}")
        if total != males + females:
            raise ValueError(
                f"{code}: P02 all-sex total {total} != male+female "
                f"{males}+{females}={males+females}"
            )
        records.append({
            "NUTS_ID": code, "year": 2020,
            "population_jan1": total,
            "ons_p02_males": males, "ons_p02_females": females,
            "ons_sheet_row": row_number,
            "population_source": "ONS 1 Jan 2020 UDC P02 all ages both sexes",
        })
    ons = pd.DataFrame(records)
    if len(ons) != N_UK or not ons.NUTS_ID.is_unique or set(ons.NUTS_ID) != ids:
        raise ValueError(
            f"ONS P02 should contain exact {N_UK} UK NUTS2016 codes; "
            f"found {len(ons)} unique={ons.NUTS_ID.nunique()}"
        )
    # Confirm reproducibility with the investigator's first verified sheet
    # row; require explicitly reviewed workbook content, not any XLSX.
    if int(ons.loc[ons.NUTS_ID.eq("UKK11"), "population_jan1"].iloc[0]) != 467617:
        raise ValueError("ONS workbook is not the reviewed Jan-2020 UKK11=467617 version")
    national = int(ons.population_jan1.sum())
    if not 60_000_000 < national < 75_000_000:
        raise ValueError(f"Implausible ONS UK total, check age/sex categories: {national}")
    return ons.sort_values("NUTS_ID").reset_index(drop=True)


def uk_jan1_annual(source, ons, ids):
    with source.open(encoding="utf8") as stream:
        payload = json.load(stream)
    p = prepare_population(payload)
    p["NUTS_ID"] = p.NUTS_ID.astype(str)
    p = p[p.NUTS_ID.isin(ids) & p.year.between(2014, 2019)].copy()
    p["population_jan1"] = pd.to_numeric(p.population_jan1, errors="coerce")
    if not np.isfinite(p.population_jan1.to_numpy(float)).all() or (
        p.population_jan1 <= 0
    ).any():
        raise ValueError("Missing/non-positive Jan1 2014–2019 UK population")
    p = p[["NUTS_ID","year","population_jan1"]].copy()
    p["population_source"] = "Eurostat demo_r_pjanaggr3 Jan1 T TOTAL NR"
    append = ons[["NUTS_ID","year","population_jan1","population_source"]]
    a = pd.concat([p, append], ignore_index=True)
    year_region_validation(a, ids, "Jan1 UK annual")
    return a


def fixed_area(source, ids):
    with source.open(encoding="utf8") as stream:
        payload=json.load(stream)
    area=prepare_area(payload)
    land=fixed_land_areas(area)
    land=land[land.NUTS_ID.isin(ids)].copy()
    if len(land)!=N_UK or not land.NUTS_ID.is_unique:
        raise ValueError("Missing/duplicated UK fixed NUTS2016 areas")
    if (land.fixed_land_area_km2 <= 0).any() or (
        land.area_change_over_1pct.fillna(False)
    ).any():
        raise ValueError("UK land-area instability or invalid denominator")
    return land


def annual_prc(path, ids, taxon):
    cols=["NUTS_ID","year","prc_mean","snow_status",
          "cumulative_prc_change_absolute"]
    d=pd.read_csv(path,dtype={"NUTS_ID":str},usecols=cols)
    d=d[d.NUTS_ID.isin(ids)&d.year.isin(YEARS)].copy()
    year_region_validation(d,ids,f"{taxon} PRC")
    if not d.snow_status.eq("locked").all() or not (
        pd.to_numeric(d.prc_mean,errors="coerce").between(0,1).all()
    ):
        raise ValueError(f"{taxon} has invalid unlocked or nonfinite PRC values")
    return d[["NUTS_ID","year","prc_mean","cumulative_prc_change_absolute"]].rename(
        columns={
            "prc_mean": f"{taxon}_prc",
            "cumulative_prc_change_absolute": f"{taxon}_source_delta",
        }
    )


def annual_radiance(path, ids):
    d=pd.read_parquet(path,columns=[
        "NUTS_ID","year","radiance_mean_corrected"
    ])
    d["NUTS_ID"]=d.NUTS_ID.astype(str)
    d=d[d.NUTS_ID.isin(ids)&d.year.isin(YEARS)].copy()
    year_region_validation(d,ids,"VIIRS corrected radiance")
    x=pd.to_numeric(d.radiance_mean_corrected,errors="coerce")
    if not np.isfinite(x.to_numpy(float)).all() or (x<0).any():
        raise ValueError("Missing/invalid already calibrated UK VIIRS radiance")
    return d


def compile_uk_annual(pop,land,rad,prcs,ids):
    a=pop.merge(
        land[["NUTS_ID","fixed_land_area_km2","land_area_reference_year"]],
        on="NUTS_ID",validate="many_to_one"
    )
    a["density_jan1_fixed_land"] = (
        a.population_jan1 / a.fixed_land_area_km2
    )
    a=a.merge(rad,on=["NUTS_ID","year"],validate="one_to_one")
    for taxon in ("spider","moth"):
        a=a.merge(prcs[taxon],on=["NUTS_ID","year"],validate="one_to_one")
    year_region_validation(a,ids,"Combined UK yearly")
    if not (a.density_jan1_fixed_land > 0).all():
        raise ValueError("Bad UK derived density")
    baseline=a[a.year.eq(2014)].set_index("NUTS_ID")
    for taxon in ("spider","moth"):
        observed=(a[f"{taxon}_prc"]-
                  a.NUTS_ID.map(baseline[f"{taxon}_prc"]))
        if not np.allclose(
            a[f"{taxon}_source_delta"],observed,rtol=1e-5,atol=2e-6
        ):
            raise ValueError(f"{taxon} 2014 cumulative PRC deltas disagree")
    a=a.sort_values(["NUTS_ID","year"]).reset_index(drop=True)
    return a


def uk_endpoint_rows(a,old):
    b=a[a.year.isin(ENDPOINTS)].copy()
    wide=b.pivot(index="NUTS_ID",columns="year",
                 values=["population_jan1","density_jan1_fixed_land",
                         "radiance_mean_corrected","spider_prc","moth_prc"])
    if len(wide)!=N_UK:
        raise ValueError("UK endpoint pivot incomplete")
    x=pd.DataFrame(index=wide.index,columns=old.columns)
    x["NUTS_ID"]=x.index
    x["country"]="UK"
    ref=a.drop_duplicates("NUTS_ID").set_index("NUTS_ID")
    for col in ["fixed_land_area_km2","land_area_reference_year"]:
        if col not in x:
            raise ValueError(f"Original frozen matched data missing {col}")
        x[col]=ref[col]
    for y in ENDPOINTS:
        for src,stem in [
            ("population_jan1","population_jan1"),
            ("density_jan1_fixed_land","density_jan1_fixed_land"),
            ("radiance_mean_corrected","radiance")
        ]:
            col=f"{stem}_{y}"
            if col not in x:
                raise ValueError(f"Original frozen cohort lacks expected {col}")
            x[col]=wide[(src,y)]
        if y!=2014:
            y0=2014
            for stem in ["density_jan1_fixed_land","radiance"]:
                col=(f"change_{stem}_2014_{y}" if stem=="density_jan1_fixed_land"
                     else f"radiance_delta_2014_{y}")
                if col not in x:
                    raise ValueError(f"Original frozen cohort lacks expected {col}")
                x[col]=x[f"{stem}_{y}"]-x[f"{stem}_{y0}"]
            for taxon in ("spider","moth"):
                col=f"cumulative_delta_{taxon}_prc_{y}"
                if col not in x:
                    raise ValueError(f"Original frozen cohort lacks expected {col}")
                x[col]=wide[(f"{taxon}_prc",y)]-wide[(f"{taxon}_prc",2014)]
    # Keep legacy official-density fields NaN for UK 2019/2020 where the
    # original Eurostat DEMO_R_D3DENS series is demonstrably unavailable.
    # They are NOT variables used by the Jan1 figures.
    required=[
        "fixed_land_area_km2",
        "density_jan1_fixed_land_2014","density_jan1_fixed_land_2020",
        "change_density_jan1_fixed_land_2014_2020",
        "radiance_2014","radiance_2020","radiance_delta_2014_2020",
        "cumulative_delta_spider_prc_2020",
        "cumulative_delta_moth_prc_2020",
    ]
    if not np.isfinite(x[required].to_numpy(float)).all():
        raise ValueError("UK expanded endpoint row lacks a primary scientific field")
    return x.reset_index(drop=True)


def descriptive_correlations(d):
    rows=[]
    variable={
        "radiance":"radiance_delta_2014_2020",
        "spider":"cumulative_delta_spider_prc_2020",
        "moth":"cumulative_delta_moth_prc_2020",
    }
    predictor="change_density_jan1_fixed_land_2014_2020"
    for country,g in [("ALL",d),*list(d.groupby("country"))]:
        for outcome,col in variable.items():
            x=g[predictor].to_numpy(float)
            y=g[col].to_numpy(float)
            if len(x)<3 or np.unique(x).size<2 or np.unique(y).size<2:
                rho=float("nan")
            else:
                rho=float(spearmanr(x,y).statistic)
            rows.append({
                "country":country,"outcome":outcome,"n":len(g),
                "spearman_rho":rho
            })
    return pd.DataFrame(rows).sort_values(["outcome","country"])


def self_test():
    assert _as_int("467617","test")==467617
    try:_as_int("nan","test")
    except ValueError: pass
    else:raise AssertionError("Non-numeric cell was not rejected")
    ids={"UKK11","UKK12"}
    toy=pd.DataFrame({
        "NUTS_ID":["UKK11"]*7+["UKK12"]*7,
        "year":list(YEARS)*2
    })
    year_region_validation(toy,ids,"toy")
    assert MATCHED_OLD+N_UK==MATCHED_NEW==1452
    print("UK ONS C=D+E numeric validation and year-region join tests: PASS")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ons",type=Path,default=ONS_XLSX)
    p.add_argument("--existing-cohort",type=Path,default=OLD)
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    self_test()
    if args.self_test:
        return
    source={
        "ons":rootpath(args.ons),
        "old_frozen":rootpath(args.existing_cohort),
        "eurostat_2014_2019":RAW/"demo_r_pjanaggr3_nuts3.json",
        "area":RAW/"reg_area3_nuts3.json",
        "radiance":RADIANCE,
        **{f"prc_{k}":v for k,v in PRC.items()},
    }
    missing=[f"{name}: {path}" for name,path in source.items()
             if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required read-only source absent:\n"+"\n".join(missing))
    nuts=load_nuts(load_config())
    ids=set(nuts.loc[nuts.CNTR_CODE.eq("UK"),"NUTS_ID"].astype(str))
    if len(ids)!=N_UK:
        raise ValueError(f"Expected {N_UK} fixed UK NUTS2016 region IDs, found {len(ids)}")
    ons=read_ons_2020_p02(source["ons"],ids)
    jan1=uk_jan1_annual(source["eurostat_2014_2019"],ons,ids)
    area=fixed_area(source["area"],ids)
    rad=annual_radiance(source["radiance"],ids)
    prcs={t:annual_prc(source[f"prc_{t}"],ids,t) for t in ("spider","moth")}
    uk_annual=compile_uk_annual(jan1,area,rad,prcs,ids)
    old=pd.read_csv(source["old_frozen"],dtype={"NUTS_ID":str,"country":str})
    if len(old)!=MATCHED_OLD or not old.NUTS_ID.is_unique or (
        old.NUTS_ID.str.startswith("UK")
    ).any():
        raise ValueError("Frozen 1,273-region old source changed unexpectedly")
    uk=uk_endpoint_rows(uk_annual,old)
    if set(old.NUTS_ID)&set(uk.NUTS_ID):
        raise ValueError("Duplicated region IDs between original cohort and UK")
    expanded=pd.concat([old,uk],ignore_index=True)
    if len(expanded)!=MATCHED_NEW or not expanded.NUTS_ID.is_unique:
        raise ValueError("Expanded matching cohort must have exactly 1,452 unique regions")
    # Appending NA-only legacy fields may upcast an original integer column
    # without changing its values. Preserve *every* original cell exactly,
    # rather than rejecting an in-memory NumPy dtype promotion.
    pd.testing.assert_frame_equal(
        old, expanded.iloc[:len(old)][old.columns].reset_index(drop=True),
        check_dtype=False, check_exact=True
    )
    # Validate all 1,452 using the same published endpoint relationships.
    popvar="change_density_jan1_fixed_land_2014_2020"
    if not np.allclose(
        pd.to_numeric(expanded[popvar],errors="coerce"),
        pd.to_numeric(expanded.density_jan1_fixed_land_2020,errors="coerce")
        -pd.to_numeric(expanded.density_jan1_fixed_land_2014,errors="coerce"),
        atol=1e-8,rtol=1e-10
    ):
        raise ValueError("Expanded cohort population changes inconsistent")
    for col in [
        popvar,"radiance_delta_2014_2020",
        "cumulative_delta_spider_prc_2020",
        "cumulative_delta_moth_prc_2020",
    ]:
        if not np.isfinite(pd.to_numeric(expanded[col],errors="coerce")).all():
            raise ValueError(f"Expanded cohort has non-finite {col}")

    out=rootpath(args.out)
    out.mkdir(parents=True,exist_ok=True)
    ons.to_csv(out/"uk_jan1_2020_ons_p02_179_validated.csv",index=False)
    uk_annual.to_csv(out/"uk_nuts3_annual_2014_2020_validated.csv",index=False)
    expanded.to_csv(
        out/"matched_region_change_2014_2020_with_uk.csv",index=False
    )
    country=descriptive_correlations(expanded)
    country.to_csv(
        out/"within_country_correlations_2014_2020_with_uk.csv",index=False
    )
    br=uk_annual.loc[uk_annual.NUTS_ID.eq("UKK11") &
                     uk_annual.year.isin((2019,2020))]
    if list(br.population_jan1.astype(int))!=[467099,467617]:
        raise ValueError("UKK11 2019->2020 control value mismatch")
    manifest={
        "years":list(YEARS),"change_endpoints":[2014,2020],
        "original_cohort_regions":len(old),
        "uk_nuts2016_regions_added":len(uk),
        "expanded_regions":len(expanded),
        "original_cohort_rows_preserved_without_modification":True,
        "matched_geography":"fixed NUTS2016 NUTS3 polygons",
        "uk_population_source_2014_to_2019":
            "Eurostat demo_r_pjanaggr3, Jan1, TOTAL both sexes, original snapshot",
        "uk_population_source_2020":
            "ONS UDC Jan1 2020, P02 NUTS3, column C = D + E",
        "ons_official_url":ONS_URL,
        "ons_workbook_sha256":hashlib.sha256(
            source["ons"].read_bytes()
        ).hexdigest(),
        "ons_2020_total_uk":int(ons.population_jan1.sum()),
        "uk_landed_fixed_area_source":str(source["area"]),
        "uk_radiance_source":str(RADIANCE),
        "uk_taxon_prc_sources":{t:str(f) for t,f in PRC.items()},
        "density_definition":"population on Jan1 / same fixed land area",
        "prior_1_273_old_2020_source_does_not_get_changed":True,
        "legacy_official_density_absent_2019_2020_for_UK":True,
        "new_correlations_are_descriptive_and_VERSIONED":True,
        "source_paths":{k:str(v) for k,v in source.items()},
        "warning":[
          "ONS 2020 population is from official UDC workbook, not Eurostat 2026 snapshot",
          "Separate source vintages must be disclosed in publication methods",
          "No PRC transfer parameters, VIIRS radiance or original output overwritten",
          "Comparison may have a 2019/2020 source-revision discontinuity; investigate if needed",
          "No inference of spatial independence from NUTS3 rank correlations"
        ],
    }
    (out/"uk_inclusion_provenance.json").write_text(
        json.dumps(manifest,indent=2)+"\n",encoding="utf8"
    )
    print("\nUK-INCLUSIVE NUTS2016 COHORT VALIDATION",flush=True)
    print(f"  Original region rows unchanged: {len(old):,}",flush=True)
    print(f"  ONS UK 2020 matched region totals: {len(ons):,}",flush=True)
    print(f"  New UK annual region-years: {len(uk_annual):,}",flush=True)
    print(f"  New UK NUTS3 rows: {len(uk):,}",flush=True)
    print(f"  Expanded region cohort: {len(expanded):,}",flush=True)
    print(f"  UK Jan1 2020 total: {int(ons.population_jan1.sum()):,}",flush=True)
    print(f"  Bristol UKK11 2019->2020: 467099 -> 467617",flush=True)
    print("  Country correlation sample size:",
          country.loc[country.country.eq("UK"),["outcome","n","spearman_rho"]].
          to_string(index=False),flush=True)
    print("\nOUTPUTS",flush=True)
    for f in sorted(out.glob("*")):
        print(f"  {f} ({f.stat().st_size/1024:.1f} KiB)",flush=True)


if __name__=="__main__":
    main()
