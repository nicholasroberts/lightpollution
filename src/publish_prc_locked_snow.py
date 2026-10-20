#!/usr/bin/env python3
"""Publish native-pixel 2021 snow-adjusted PRC as standard spider/moth analyses."""
import argparse
import json
import numpy as np
import pandas as pd
import yaml
from project_config import PROJECT_ROOT as ROOT
from historical_prc_pixel_nuts3 import build_panel,europe_table,load_nuts,load_config,plot_region_atlases,plot_europe

SETS={
 "spider":("prc_2021_locked17_combined","historical_prc_2014_2024_locked17","config/prc_transfer.yaml"),
 "moth":("prc_2021_locked_moth16_combined","historical_prc_2014_2024_locked_moth16","config/prc_transfer_moth.yaml")}

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--self-test",action="store_true")
 p.add_argument("--no-maps",action="store_true")
 a=p.parse_args()
 model=yaml.safe_load((ROOT/"config/viirs_2021_snow_correction.yaml").read_text())
 if model["status"]!="locked" or model["applies_to_year"]!=2021:
  raise ValueError("Shared snow adjustment not locked")
 for _,_,config in SETS.values():
  if yaml.safe_load((ROOT/config).read_text())["status"]!="locked":
   raise ValueError("Taxon PRC not locked")
 if a.self_test:
  print("Both PRCs and shared 2021 VIIRS snow correction locked: PASS")
  return
 for taxon,(candidate,archive,_) in SETS.items():
  old=pd.read_csv(ROOT/"data/processed"/archive/"prc_nuts3_annual_cumulative.csv",dtype={"NUTS_ID":str})
  new=pd.read_csv(ROOT/"data/processed"/candidate/"prc_nuts3_annual_cumulative_candidate.csv",dtype={"NUTS_ID":str})
  keys=["NUTS_ID","year"]
  j=old[keys+["prc_mean","valid_area_km2"]].merge(
    new[keys+["prc_mean","valid_area_km2"]],on=keys,validate="one_to_one",suffixes=("_old","_new"))
  if len(j)!=16742 or len(old)!=len(new):
   raise ValueError("Unmatched NUTS3 years")
  unchanged=j.year.ne(2021)
  if not np.array_equal(j.loc[unchanged,"prc_mean_old"],j.loc[unchanged,"prc_mean_new"]):
   raise ValueError("Non-2021 PRC modified")
  if not np.allclose(j.valid_area_km2_old,j.valid_area_km2_new,atol=.05):
   raise ValueError("Pixel coverage changed")
  gain=j.loc[~unchanged,"prc_mean_new"]-j.loc[~unchanged,"prc_mean_old"]
  if not (gain.ge(-5e-6).all() and gain.gt(1e-7).any()):
   raise ValueError("Snow adjustment absent or reversed")
  out=ROOT/"data/processed"/("historical_prc_2014_2024_"+taxon+"_snow_corrected")
  figs=ROOT/"figures"/("historical_prc_2014_2024_"+taxon+"_snow_corrected")
  out.mkdir(parents=True,exist_ok=True)
  figs.mkdir(parents=True,exist_ok=True)
  panel=build_panel(new,list(range(2014,2025)),.95)
  panel["snow_status"]="locked"
  panel.to_csv(out/"prc_nuts3_annual_cumulative.csv",index=False)
  panel.to_parquet(out/"prc_nuts3_annual_cumulative.parquet",index=False)
  europe=europe_table(panel)
  europe.to_csv(out/"prc_europe_annual_cumulative.csv",index=False)
  last=panel.loc[panel.year.eq(2024)&panel.coverage_stable].copy()
  last["rank_largest_gain"]=last.cumulative_prc_change_absolute.rank(ascending=False)
  last["rank_largest_loss"]=last.cumulative_prc_change_absolute.rank(ascending=True)
  last.to_csv(out/"prc_nuts3_rankings_2014_2024.csv",index=False)
  if not a.no_maps:
   nuts=load_nuts(load_config())
   plot_region_atlases(panel,nuts,figs,list(range(2014,2025)))
   plot_europe(europe,figs)
  rawmean=europe_table(old)
  before=float(rawmean.loc[rawmean.year.eq(2021),"area_weighted_prc"].iloc[0])
  after=float(europe.loc[europe.year.eq(2021),"area_weighted_prc"].iloc[0])
  info={"status":"locked","taxon":taxon,"2021_prc_original":before,
        "2021_prc_corrected":after,"2021_increment":after-before,
        "unmodified_years":True,"nuts3_regions":1522}
  (out/"snow_correction_summary.json").write_text(json.dumps(info,indent=2)+"\n")
  print(taxon,info,flush=True)

if __name__=="__main__":
 main()
