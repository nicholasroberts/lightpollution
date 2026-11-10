#!/usr/bin/env python3
"""Madrid: absolute rather than percentage-only comparison of ground vs VIIRS.

Uses measured local monthly VNL V1 csv from madrid_ucm_local_ground_bridge.py
and the independently published quality-filtered Robles et al (2021)
winter 01:30 Johnson B/V/R figures in Tables 3–5.
https://doi.org/10.3390/rs13081511

WARNING: AstMon measures ground zenith sky radiance, VIIRS measures upward
sensor-weighted band radiance. Same SI dimensions DOES NOT imply equivalent
fluxes or a photon budget. VIIRS spectrum is unknown: 700 nm photon
equivalent is illustrative only. Vega filter conversions approximate.
Figure panels use DIFFERENT y-axis scales. Do not sum the B/V/R bands.
No frozen NUTS/VIIRS/PRC sources are altered.

    python src/plot_madrid_absolute_photon_radiance_2015_2019.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BASE = "local_lighting_interventions/madrid_ground_viirs/local_ucm"
RAW = ROOT / "data/processed" / BASE
FIG = ROOT / "figures" / BASE
VIIRS = RAW / "madrid_ucm_3km_10km_2015_2019_monthly.csv"
H = 6.62607015e-34
C = 299792458.
ONE_ARCSEC2_SR = (np.pi/(180.*3600.))**2
# Published 2015,2019 winter 01:30 mag/arcsec²; approximate Vega zero W m^-2 nm^-1;
# nominal band wavelength nm and FWHM nm. R is original Johnson R, NOT Cousins R.
BANDS = {
    "Johnson B": (19.30, 18.97, 6.31e-11, 436., 89.),
    "Johnson V": (18.07, 17.96, 3.607e-11, 545., 84.),
    "Johnson R": (17.42, 17.45, 1.76e-11, 700., 220.),
}

def run():
    d = pd.read_csv(VIIRS)
    required = {"support", "year", "month", "cf_threshold",
                "mean_radiance", "n_valid_pixels", "n_possible_pixels"}
    if not required.issubset(d.columns):
        raise ValueError("Missing original monthly Madrid UCM source columns")
    out = []
    ground = {}
    for name, (m15, m19, zero, wl, width) in BANDS.items():
        phot = {}
        for year, mag in ((2015, m15), (2019, m19)):
            L = zero * 10**(-.4 * mag) * width / ONE_ARCSEC2_SR
            phot[year] = (L, L*wl*1e-9/(H*C))
        ground[name] = phot
        a, b = phot[2015][0], phot[2019][0]
        out.append({
            "measure": name, "instrument": "ground UCM zenith",
            "time_basis": "published 01:30 winter Nov-Jan",
            "radiance_2015_uW_m2_sr": a*1e6,
            "radiance_2019_uW_m2_sr": b*1e6,
            "delta_uW_m2_sr": (b-a)*1e6,
            "delta_photons_1e12_m2_s_sr": (phot[2019][1]-phot[2015][1])/1e12,
            "percent_change": 100*(b/a-1),
            "wavelength_nm_assumed": wl,
            "note": "approx Johnson Vega zero point & FWHM"})
    viirs = {}
    for area in ("UCM_3km", "UCM_10km"):
        sub = d.loc[d.support.eq(area) & d.year.isin((2015,2019)) &
                    d.month.isin((1,11,12)) & d.cf_threshold.eq(3)]
        if (len(sub) != 6 or sub.duplicated(["year","month"]).any() or
                not (sub.n_valid_pixels == sub.n_possible_pixels).all()):
            raise ValueError(f"Missing/invalid 2015 or 2019 winter VIIRS {area}")
        annual = sub.groupby("year").mean_radiance.mean()
        a, b = (float(annual.loc[y])*1e-5 for y in (2015,2019))
        viirs[area] = {2015:a, 2019:b}
        out.append({
            "measure": area, "instrument": "VIIRS VNL V1 upward",
            "time_basis": "calendar-year Jan/Nov/Dec; not same winters",
            "radiance_2015_uW_m2_sr": a*1e6,
            "radiance_2019_uW_m2_sr": b*1e6,
            "delta_uW_m2_sr": (b-a)*1e6,
            "delta_photons_1e12_m2_s_sr": (b-a)*700e-9/(H*C)/1e12,
            "percent_change": 100*(b/a-1),
            "wavelength_nm_assumed": 700.,
            "note": "DNB 700nm photon equivalent, NOT measured spectrum"})
    table = pd.DataFrame(out)
    RAW.mkdir(parents=True, exist_ok=True)
    table.to_csv(RAW/"madrid_2015_2019_absolute_radiance_photon_comparison.csv",
                 index=False)
    plt.rcParams.update({"font.size":10, "axes.spines.top":False,
                         "axes.spines.right":False, "pdf.fonttype":42,
                         "ps.fonttype":42, "svg.fonttype":"none"})
    fig, axes = plt.subplots(2,2,figsize=(13.6,10.2),
                            gridspec_kw={"height_ratios":[1.08,1]},
                            constrained_layout=True)
    ag, av, apg, apv = axes.ravel()
    gnames = list(BANDS)
    vkeys = ["UCM_3km","UCM_10km"]
    def compare(ax,names,left,right,title):
        x=np.arange(len(names))
        b1=ax.bar(x-.17,left,.34,label="2015")
        b2=ax.bar(x+.17,right,.34,label="2019")
        ax.bar_label(b1,fmt="%.1f",padding=3,fontsize=8.5)
        ax.bar_label(b2,fmt="%.1f",padding=3,fontsize=8.5)
        ax.set_xticks(x,names)
        ax.set_ylim(0,max(max(left),max(right))*1.3)
        ax.set_ylabel("Radiance (µW m⁻² sr⁻¹)")
        ax.set_title(title,loc="left",fontweight="bold")
        ax.legend(frameon=False,ncol=2)
        ax.grid(axis="y",alpha=.16)
    compare(ag,gnames,[ground[n][2015][0]*1e6 for n in gnames],
            [ground[n][2019][0]*1e6 for n in gnames],
            "A  Ground zenith sky radiance")
    compare(av,["UCM 3 km","UCM 10 km"],
            [viirs[n][2015]*1e6 for n in vkeys],
            [viirs[n][2019]*1e6 for n in vkeys],
            "B  Orbital upward VIIRS radiance")
    def photon(ax,names,keys,title):
        x=np.arange(len(names))
        values=np.array([table.set_index("measure").loc[k,"delta_photons_1e12_m2_s_sr"]
                         for k in keys],dtype=float)
        perc=np.array([table.set_index("measure").loc[k,"percent_change"]
                       for k in keys],dtype=float)
        ax.bar(x,values)
        ax.axhline(0,lw=1,color="0.3")
        ax.set_xticks(x,names)
        ax.set_ylabel("Δ photon equivalent (10¹² m⁻² s⁻¹ sr⁻¹)")
        ax.set_title(title,loc="left",fontweight="bold")
        for xx,v,p in zip(x,values,perc):
            ax.annotate(f"{v:+.1f}\n({p:+.1f}%)",(xx,v),
                        xytext=(0,7 if v>=0 else -28),
                        textcoords="offset points",ha="center",
                        fontsize=9,fontweight="bold")
        ax.grid(axis="y",alpha=.16)
        return values
    photon(apg,gnames,gnames,"C  Ground: estimated photon-radiance change")
    photon(apv,["UCM 3 km","UCM 10 km"],vkeys,
           "D  VIIRS: illustrative 700 nm photon equivalent")
    # Represent unknown DNB wavelength as hypothetical 500–900 nm range,
    # NOT measurement error bars and NOT an invertible spectral estimate.
    for i,k in enumerate(vkeys):
        delta=viirs[k][2019]-viirs[k][2015]
        vals=np.array([delta*wl*1e-9/(H*C)/1e12
                       for wl in (500.,900.)])
        apv.vlines(i,min(vals),max(vals),color="0.25",lw=2,alpha=.5)
        apv.hlines(vals,i-.09,i+.09,color="0.25",lw=2,alpha=.5)
    apv.set_ylim(-180,22)
    apv.text(.02,.04,
             "Grey whiskers: 500–900 nm assumed wavelength,\n"
             "NOT measured spectral uncertainty",
             transform=apv.transAxes,fontsize=8)
    fig.suptitle(
        "Madrid 2015–2019  |  Absolute radiance and photon-equivalent changes\n"
        "Same dimensions; different wavelengths, viewpoints and source areas",
        fontsize=15,weight="bold")
    fig.text(.03,-.018,
        "Ground: Robles et al. (2021), Tables 3–5, filtered winter 01:30. "
        "Johnson zero points and nominal FWHM are approximate (especially R).\n"
        "VIIRS: measured 3/10 km VNL V1, equal means of Jan, Nov, Dec "
        "within each calendar year. DNB is spectrally response-weighted; "
        "photons at 700 nm are an EQUIVALENT, not a measured quantity.\n"
        "Do not subtract ground zenith radiance from upward orbital radiance "
        "to infer total emitted photons. B/V/R filter responses overlap; "
        "do not sum the bands.",
        fontsize=8.2,color="0.3",va="top")
    FIG.mkdir(parents=True,exist_ok=True)
    stem="madrid_absolute_radiance_photons_2015_2019"
    for suffix in ("png","pdf","svg"):
        fig.savefig(FIG/f"{stem}.{suffix}",dpi=185,bbox_inches="tight")
    plt.close(fig)
    print(table.to_string(index=False))
    print("Figure:", FIG/f"{stem}.png")

if __name__ == "__main__":
    run()
