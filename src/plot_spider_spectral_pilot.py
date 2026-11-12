"""Summarize visible spectral scenarios and the separate unvalidated UV proxy."""
from pathlib import Path
import argparse,csv
from datetime import datetime
import numpy as np
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spider_spectral_lunar import pigment_proxy


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['spectra','hps','led','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    a.output.mkdir(parents=True)
    fig,axes=plt.subplots(2,2,figsize=(12,9))
    ax=axes[0,0]
    for label,name in [('Natural','HPS_visible_input.csv'),('HPS','HPS_visible_input.csv'),('LED 2958–3969K group','LED_mid_visible_input.csv')]:
        d=np.genfromtxt(a.spectra/name,delimiter=',',names=True);wave=d['wavelength_nm']
        values=d['natural'] if label=='Natural' else d['alan']
        ax.plot(wave,values/trapezoid(values,wave),label=label)
    ax.set(title='Source + power-law transfer scenarios; equal band energy',xlabel='Wavelength (nm)',ylabel='Normalized energy spectrum (per nm)');ax.legend(fontsize=8)
    ax=axes[0,1];wave=np.arange(300.,751.)
    for amplitude in [0,.29,1]:ax.plot(wave,pigment_proxy(wave,uv_amplitude=amplitude),label=f'UV amplitude {amplitude:g}')
    ax.axvspan(450,750,color='gray',alpha=.1,label='Current PRC calculation band')
    ax.set(title='510nm pigment proxy; UV strength is a scenario',xlabel='Wavelength (nm)',ylabel='Relative photon sensitivity');ax.legend(fontsize=8)
    for ax,site in zip(axes[1],['Cottage','oxfordbussinesspark']):
        for label,folder in [('HPS',a.hps),('LED 2958–3969K group',a.led)]:
            rows=[r for r in csv.DictReader((folder/'prc_ranges.csv').open()) if r['site']==site and float(r['atmospheric_factor'])==.8 and float(r['alan_band_ratio'])==1 and float(r['body_tilt_deg'])==0]
            # Convert UTC to each site's local hour, retaining hours >24.
            from zoneinfo import ZoneInfo
            zone=ZoneInfo('Europe/Stockholm' if site=='Cottage' else 'Europe/London')
            times=[datetime.fromisoformat(r['utc']).astimezone(zone) for r in rows];day=times[0].date()
            hours=[24*(t.date()-day).days+t.hour+t.minute/60 for t in times]
            ax.plot(hours,[float(r['median_abs_prc']) for r in rows],label=label)
        ax.set(title=site,xlabel='Local hour (24 = midnight)',ylabel='Median absolute PRC over body headings',ylim=(0,.7));ax.legend(fontsize=8)
    fig.suptitle('VISIBLE-ONLY pilot: horizontal body, rotating 360deg; KS91 changing moonlight\nBottom panels: atmospheric factor0.8; fixed ALAN = capture-time zenith natural band intensity\nGeneric source spectra + assumed transfer; not measured site mixtures or UV-inclusive predictions',fontsize=10)
    fig.tight_layout(rect=[0,0,1,.9]);fig.savefig(a.output/'visible_spectral_nightly_comparison.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()
