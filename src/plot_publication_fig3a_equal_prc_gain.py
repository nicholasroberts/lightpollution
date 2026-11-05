#!/usr/bin/env python3
"""Illustrate R reductions necessary for +0.05 PRC at R=80 vs R=5.

Presentation-only companion to Figure 3A, using independently LOCKED spider
and moth functions from config/prc_transfer{,_moth}.yaml. The inversion of
the exponential function is analytic; no data reconstruction or refitting.
Creates a 16:9 editable SVG, vector PDF, preview PNG, value CSV and JSON.
Neither Figure 3A/3B nor the original scientific inputs are overwritten.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from figure_style import apply_figure_style, SERIES_COLOURS, save_figure
from plot_publication_figure3_prc import transfer_functions, response, RAD_UNIT
from project_config import PROJECT_ROOT as ROOT, load_config

DEFAULT_OUT=ROOT/"figures/publication/figure3_prc"
NAME="fig3A_example_equal_0p05_prc_gain_80_vs_5"
DEFAULT_GAIN=.05
STARTS=(80.,5.)
TAXON_LABEL={"spider":"Spider","moth":"Moth"}
ORANGE="#E58A2B"
TEAL="#148774"
NAVY="#21334D"
BG="#FBFCFE"


def inverse_prc(target, parameters):
    """Analytic inverse of yf+(y0-yf)*exp(-k*R)."""
    yf,y0,k=parameters
    if not (yf<target<=y0):
        raise ValueError("Requested target PRC outside locked function range")
    result=-math.log((target-yf)/(y0-yf))/k
    if not math.isfinite(result) or result<0:
        raise ValueError("Invalid inverted radiance")
    return result


def required_reduction(taxon, start, gain, parameters):
    original=float(response(np.asarray([start]),parameters)[0])
    target=original+gain
    final=inverse_prc(target,parameters)
    if not (0<=final<start):
        raise ValueError("No valid radiance reduction achieves this PRC increase")
    if not np.isclose(
        float(response(np.asarray([final]),parameters)[0]),
        target,rtol=1e-12,atol=1e-12
    ):
        raise ValueError("Inverse function verification failed")
    decrease=start-final
    return {
        "taxon":taxon,"starting_radiance_nw":start,
        "starting_prc":original,"target_prc":target,
        "absolute_prc_increase":gain,
        "final_radiance_nw":final,"radiance_reduction_nw":decrease,
        "radiance_reduction_pct":100*decrease/start,
    }


def plot_panel(ax,taxon,start,params,result,zoom,ylim,gain):
    xmax=7. if zoom else 100.
    r=np.linspace(0,xmax,1401)
    ax.plot(r,response(r,params),lw=2.6,
            color=SERIES_COLOURS[taxon],zorder=3)
    final=result["final_radiance_nw"]
    orig_prc=result["starting_prc"]
    new_prc=result["target_prc"]
    segment=np.linspace(final,start,221)
    ax.plot(segment,response(segment,params),color=TEAL,
            lw=4.6,solid_capstyle="round",zorder=5)
    ax.plot(start,orig_prc,marker="o",markersize=9.5,
            markerfacecolor=ORANGE,markeredgecolor="white",
            markeredgewidth=1.5,zorder=7)
    ax.plot(final,new_prc,marker="D",markersize=8.5,
            markerfacecolor=TEAL,markeredgecolor="white",
            markeredgewidth=1.3,zorder=7)
    r1=final+.56*(start-final)
    r2=final+.36*(start-final)
    ax.annotate("",
        xy=(r2,float(response(np.asarray([r2]),params)[0])),
        xytext=(r1,float(response(np.asarray([r1]),params)[0])),
        arrowprops=dict(arrowstyle="-|>",mutation_scale=14,lw=1.6,
                        color=TEAL),zorder=8)
    bracket=6.2 if zoom else 90.
    ax.annotate("",xy=(bracket,new_prc),xytext=(bracket,orig_prc),
        arrowprops=dict(arrowstyle="<->",lw=1.35,color=NAVY,
                        shrinkA=0,shrinkB=0,mutation_scale=9),zorder=5)
    ax.text(bracket+(.4 if zoom else 1.5),.5*(orig_prc+new_prc),
            f"+{gain:g} PRC",ha="left",va="center",fontsize=9.5,
            color=NAVY,bbox=dict(boxstyle="round,pad=.12",
            facecolor="white",edgecolor="none",alpha=.85),zorder=9)
    for mark,color in ((final,TEAL),(start,ORANGE)):
        ax.axvline(mark,color=color,lw=.9,linestyle=":",alpha=.7,zorder=2)
    ax.set_xscale("linear")
    ax.set_xlim(0.,xmax)
    ax.set_ylim(*ylim)
    ax.set_xticks(list(range(8)) if zoom else [0,20,40,60,80,100])
    ax.tick_params(labelsize=10.5,colors="#435469",length=3.5)
    ax.grid(alpha=.16,linewidth=.65)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_facecolor("white")
    ax.set_title(f"{TAXON_LABEL[taxon]}  •  starting at {start:g}",
                 fontsize=13.2,fontweight="bold",loc="left",
                 color=NAVY,pad=12)
    ax.set_xlabel(f"Calibrated radiance ({RAD_UNIT})",
                  fontsize=10.8,labelpad=6)
    ax.set_ylabel("PRC",fontsize=11)
    assert ax.get_xscale()=="linear"
    assert np.allclose(ax.get_xlim(),(0.,xmax))


def make_figure(parameters,rows,out,gain):
    apply_figure_style()
    plt.rcParams.update({"font.size":10.4,"svg.fonttype":"none",
                         "pdf.fonttype":42,"savefig.facecolor":"white"})
    fig=plt.figure(figsize=(16,9),facecolor=BG)
    fig.text(.055,.965,
        f"How much must radiance fall to recover +{gain:g} PRC?",
        fontsize=20.5,fontweight="bold",color=NAVY,ha="left",va="top")
    fig.text(.055,.922,
        "Locked Figure 3A spider and moth response functions  •  "
        "Same PRC gain, very different starting radiance",
        fontsize=11.3,color="#54657A",ha="left",va="top")
    panels=[
        ("spider",80.,[.070,.552,.425,.285],False,(.127,.47)),
        ("spider",5., [.562,.552,.373,.285],True,(.295,.475)),
        ("moth",80., [.070,.145,.425,.285],False,(.082,.72)),
        ("moth",5.,  [.562,.145,.373,.285],True,(.40,.70)),
    ]
    for taxon,start,pos,zoom,ylim in panels:
        result=rows[(taxon,start)]
        ax=fig.add_axes(pos)
        plot_panel(ax,taxon,start,parameters[taxon],result,zoom,ylim,gain)
        fig.text(pos[0],pos[1]-.068,
            f"{start:g} → {result['final_radiance_nw']:.2f}  |  "
            f"−{result['radiance_reduction_nw']:.2f} nW  "
            f"({result['radiance_reduction_pct']:.0f}% reduction)",
            fontsize=12.3,fontweight="bold",color=TEAL)
    handles=[
        Line2D([],[],marker="o",linestyle="None",color=ORANGE,
               markeredgecolor="white",markeredgewidth=1.1,
               markersize=10,label="Starting radiance"),
        Line2D([],[],marker="D",linestyle="None",color=TEAL,
               markeredgecolor="white",markeredgewidth=1.1,
               markersize=9,label=f"Radiance after +{gain:g} PRC"),
    ]
    fig.legend(handles=handles,loc="lower left",bbox_to_anchor=(.059,.031),
               ncol=2,frameon=False,fontsize=10.3,handlelength=1.5,
               columnspacing=2.4)
    fig.text(.935,.033,
        "Illustrative inverse calculation, not measured lighting "
        "restoration; PRC is dimensionless.",
        ha="right",va="center",fontsize=9.4,color="#5F6D7D")
    save_figure(fig,out/NAME,pad_inches=.12)
    fig.savefig(out/f"{NAME}.png",dpi=220,format="png",
                bbox_inches="tight",pad_inches=.12,facecolor=BG)
    plt.close(fig)


def self_test():
    p=(.1453159020633674,.45612835620720266,.1037115607084461)
    r=required_reduction("spider",80.,.05,p)
    assert abs(r["final_radiance_nw"]-17.60284263)<1e-6
    r=required_reduction("spider",5.,.05,p)
    assert abs(r["final_radiance_nw"]-2.69388333)<1e-6
    p=(.10128334523974113,.6870948746846356,.0782312036944599)
    r=required_reduction("moth",80.,.05,p)
    assert abs(r["final_radiance_nw"]-31.17421787)<1e-6
    r=required_reduction("moth",5.,.05,p)
    assert abs(r["final_radiance_nw"]-3.48069856)<1e-6
    print("Figure 3A inverse PRC/radiance calibration test: PASS",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",default=DEFAULT_OUT,type=Path)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    self_test()
    if args.self_test:
        return
    parameters,model_details=transfer_functions(load_config())
    results={
        (taxon,start):required_reduction(taxon,start,DEFAULT_GAIN,p)
        for taxon,p in parameters.items() for start in STARTS
    }
    output=args.out if args.out.is_absolute() else ROOT/args.out
    output.mkdir(parents=True,exist_ok=True)
    rows=[results[(taxon,start)] for taxon in ("spider","moth")
          for start in STARTS]
    pd.DataFrame(rows).to_csv(
        output/"fig3A_example_equal_0p05_prc_gain_calculations.csv",index=False)
    (output/"fig3A_example_equal_0p05_prc_gain_provenance.json").write_text(
        json.dumps({
            "source":"locked Figure 3A functions, no refitting",
            "transfer_functions":model_details,
            "absolute_prc_increase":DEFAULT_GAIN,
            "initial_radiance_nw":list(STARTS),
            "calculations":rows,
            "formula":"Rfinal=-ln((PRC_start+gain-yf)/(y0-yf))/k",
            "radiance_axis":"linear 0–100 plus linear 0–7 close-up",
            "original_figures_3A_and_3B_unchanged":True,
            "scientific_data_modified":False,
        },indent=2)+"\n",encoding="utf8")
    make_figure(parameters,results,output,DEFAULT_GAIN)
    print("Rendered Figure 3A annotated inverse-gain example:",output,flush=True)
    for row in rows:
        print(row["taxon"],row["starting_radiance_nw"],"->",
              round(row["final_radiance_nw"],5),"reduction",
              round(row["radiance_reduction_nw"],5),"nW",flush=True)


if __name__=="__main__":
    main()
