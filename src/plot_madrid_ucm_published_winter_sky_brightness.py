#!/usr/bin/env python3
"""Madrid winter 01:30 AstMon published sky brightness vs local VIIRS.

Standalone, read-only comparison; requires existing --viirs output from
madrid_ucm_local_ground_bridge.py. Does NOT use the Zenodo unfiltered
AstMon_evol.csv (which has NO cloud/moon quality fields).

AstMon data: Robles et al. (2021) Tables 3-5, winter Nov-Jan, *published*
binned 01:30 measurements, 2015 and 2019 only:
 https://doi.org/10.3390/rs13081511

VIIRS data: 2015 vs 2019 *calendar-year* Jan/Nov/Dec mean, 3km/10km
centred on UCM, VCMSLCFG cf_cvg>=3, same 3 months in each year.
CAUTION: publication's "winter 2015" grouping may cross calendar years;
the VIIRS month selection is a similar-season sensitivity test, NOT an
exact time-paired or equivalent night set. Winter may be stratified
differently than the author's source data.

Run:
  python src/plot_madrid_ucm_published_winter_sky_brightness.py --self-test
  python src/plot_madrid_ucm_published_winter_sky_brightness.py

Outputs: figures/local_lighting_interventions/madrid_ground_viirs/local_ucm/
         data/processed/local_lighting_interventions/madrid_ground_viirs/local_ucm/
"""
from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data/processed/local_lighting_interventions/madrid_ground_viirs/local_ucm"
VIIRS = DATA / "madrid_ucm_3km_10km_2015_2019_monthly.csv"
FIGURES = ROOT / "figures/local_lighting_interventions/madrid_ground_viirs/local_ucm"

# Robles et al. 2021 Tables 3–5, "Binned Winter Evolution", 01:30.
# Values are mag arcsec-2 with reported ± spreads (not confidence intervals).
# Paper defines "winter" as Nov–Jan. Don't invent intermediate ground years.
ASTMON = [
    {"band":"B","mag_2015":19.30,"spread_2015":0.07,
     "mag_2019":18.97,"spread_2019":0.06},
    {"band":"V","mag_2015":18.07,"spread_2015":0.11,
     "mag_2019":17.96,"spread_2019":0.02},
    {"band":"R","mag_2015":17.42,"spread_2015":0.07,
     "mag_2019":17.45,"spread_2019":0.03},
]


def build_changes(local):
    required = {"support","year","month","cf_threshold","mean_radiance",
                "n_valid_pixels","n_possible_pixels"}
    if not required.issubset(local):
        raise ValueError(f"VIIRS local source lacks fields {required-set(local)}")
    sub = local[local.cf_threshold.eq(3) &
                local.year.isin([2015,2019])].copy()
    # Reject duplicate calendar months or absent months.
    for support in ("UCM_3km","UCM_10km"):
        for year in (2015,2019):
            part = sub[sub.support.eq(support) & sub.year.eq(year)]
            if len(part)!=12 or set(part.month)!=set(range(1,13)):
                raise ValueError(f"Incomplete/duplicated {support} {year}")
            if (part.n_valid_pixels < part.n_possible_pixels).any():
                raise ValueError("VIIRS support changes; use strict paired pixels")
            if not (part.mean_radiance > 0).all():
                raise ValueError("Non-positive radiance")
    rows = []
    # Same year Jan, Nov, Dec is a DESCRIPTIVE proxy, not guaranteed
    # the exact cross-year Nov-Jan winter used by the authors.
    for support in ("UCM_3km","UCM_10km"):
        d=sub[sub.support.eq(support)&sub.month.isin([1,11,12])]
        means=d.groupby("year").mean_radiance.mean()
        rows.append({
            "measurement":support,"source":"VIIRS VNL V1",
            "temporal_support":"calendar-year Jan, Nov, Dec",
            "baseline_2015":float(means.loc[2015]),
            "outcome_2019":float(means.loc[2019]),
            "percent_change":float(100*(means.loc[2019]/means.loc[2015]-1)),
            "units":"nW cm^-2 sr^-1",
            "quality":"cf_cvg>=3; all pixels valid; not overpass-night matched"})
    for record in ASTMON:
        a,b=record["mag_2015"],record["mag_2019"]
        rows.append({
            "measurement":"AstMon Johnson "+record["band"],
            "source":"ground zenith, Robles Tables 3-5",
            "temporal_support":"published Nov-Jan winter, 01:30 bin",
            "baseline_2015":a,"outcome_2019":b,
            "percent_change":float(100*(10**(-.4*(b-a))-1)),
            "units":"mag arcsec^-2",
            "quality":"original paper clear/moonless filtered; published ± not CI"})
    return pd.DataFrame(rows)


def figure(frame, out, formats=("png","pdf","svg")):
    plt.rcParams.update({
        "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none",
        "font.size":10, "axes.spines.top":False,
        "axes.spines.right":False})
    fig,ax=plt.subplots(figsize=(11.8,5.6),constrained_layout=True)
    x=np.array([0,1,3,4,5])
    col=["#33536d","#6793af","#377cb3","#439166","#b05b4b"]
    values=frame.percent_change.to_numpy(float)
    ax.bar(x,values,color=col,width=.64,zorder=3)
    for xx,yy in zip(x,values):
        ax.annotate(f"{yy:+.1f}%",(xx,yy),
                    xytext=(0,6 if yy>0 else -18),textcoords="offset points",
                    ha="center",fontweight="bold")
    ax.axhline(0,color="#444",linewidth=1)
    ax.axvline(2,color="#aaa",linestyle=":",linewidth=1)
    ax.set_xticks(x,[
        "VIIRS UCM\n3 km","VIIRS UCM\n10 km",
        "AstMon\nJohnson B","AstMon\nJohnson V","AstMon\nJohnson R"])
    ax.set_ylabel("2015→2019 % change (positive = brighter)")
    ax.set_title(
        "Madrid winter sky brightness at 01:30 versus UCM-local VIIRS\n"
        "Ground data are published filtered observations, not raw Zenodo CSV",
        fontsize=13,weight="bold")
    ax.grid(axis="y",alpha=.18,zorder=0)
    ax.text(.19,.96,
            "Orbital upward radiance\nJan/Nov/Dec calendar-month mean",
            transform=ax.transAxes,ha="center",va="top",fontsize=8)
    ax.text(.77,.96,
            "UCM zenith sky brightness\npublished Nov–Jan winter 01:30 bin",
            transform=ax.transAxes,ha="center",va="top",fontsize=8)
    fig.text(.01,-.03,
            "Robles et al. 2021, Tables 3–5; ground winter Jan/Nov/Dec "
            "may cross calendar-year boundaries. VIIRS same-year "
            "Jan/Nov/Dec is only a seasonal sensitivity proxy. Different "
            "spectra, views, year groupings and spatial source footprints. "
            "Not a simultaneous observation or causal LED test.",
            fontsize=8,ha="left",color="#555")
    out.mkdir(parents=True,exist_ok=True)
    stem="madrid_ucm_winter_0130_published_sky_brightness_vs_viirs"
    for ext in formats:
        fig.savefig(out/f"{stem}.{ext}",dpi=180,bbox_inches="tight")
    plt.close(fig)
    return stem


def self_test():
    period = pd.period_range("2015-01","2019-12",freq="M")
    fake=pd.DataFrame([
        {"support":s,"year":p.year,"month":p.month,"cf_threshold":3,
         "mean_radiance":float(40 if p.year==2015 else 38),
         "n_valid_pixels":n,"n_possible_pixels":n}
        for s,n in (("UCM_3km",172),("UCM_10km",1913)) for p in period])
    out=build_changes(fake)
    assert len(out)==5
    assert np.isclose(out.percent_change.iloc[0],-5)
    assert 35 < out.percent_change.iloc[2] < 36
    with TemporaryDirectory() as tmp:
        stem=figure(out,Path(tmp),formats=("png",))
        assert (Path(tmp)/f"{stem}.png").exists()
    print("Self-test PASS: approved published 01:30 Table 3–5 values, "
          "exact month validity, mag conversion and plot rendering.")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:
        self_test()
        return
    if not VIIRS.is_file():
        raise FileNotFoundError(
            f"{VIIRS}; first run "
            "python src/madrid_ucm_local_ground_bridge.py --viirs")
    changes=build_changes(pd.read_csv(VIIRS))
    DATA.mkdir(parents=True,exist_ok=True)
    changes.to_csv(DATA/"madrid_ucm_winter_0130_published_vs_viirs.csv",
                   index=False)
    stem=figure(changes,FIGURES)
    print(changes.to_string(index=False))
    print("Figure:",FIGURES/(stem+".png"))


if __name__=="__main__":
    main()
