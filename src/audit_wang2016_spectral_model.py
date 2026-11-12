"""Audit corrected supplement coefficients, spectral DoLP and AoP covariance.

No clamping or empirical fitting. Full equation raw outputs and failures are
reported explicitly; this is a validation gate before lunar interpretation.
"""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wang2016_sky import (complex_pattern,source_meridian_pattern,wang_spectral_inputs,
 horizon_factor,turbidity_multiplier,zenith_chromaticity)
from reconstruct_spider_eye_contrast import digest


def covariance_audit():
    theta=np.linspace(.05,1.4,50);phi=np.linspace(0,2*np.pi,50)
    rows=[]
    for name,fn in [('printed_eq4_5',complex_pattern),('source_meridian_spin2',source_meridian_pattern)]:
        a=fn(theta,phi,.7,.2,24)
        b=fn(theta,phi+.7,.7,.9,24)
        rows.append(dict(convention=name,mean_abs_dop_covariance_error=float(np.mean(abs(abs(a)-abs(b)))),
            mean_complex_spin2_error=float(np.mean(abs(b-a*np.exp(2j*.7))))))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source-elevation-deg',type=float,default=45)
    p.add_argument('--source-azimuth-deg',type=float,default=30)
    p.add_argument('--neutral-separation-deg',type=float,required=True)
    p.add_argument('--turbidity',type=float,default=3)
    a=p.parse_args()
    if a.output.exists():p.error('Use a new output folder')
    y,x=np.indices((181,181));dx=(x-90)/89;dy=(y-90)/89;r=np.hypot(dx,dy);sky=r<.995
    theta=2*np.arcsin(np.clip(r/np.sqrt(2),0,1));phi=np.arctan2(dy,dx)
    th=theta[sky];ph=phi[sky];zs=np.deg2rad(90-a.source_elevation_deg);az=np.deg2rad(a.source_azimuth_deg)
    fig,axes=plt.subplots(2,4,figsize=(14,8));results=[]
    art=plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap='turbo')
    for row,convention in enumerate(['supplement_as_printed','corrected_preetham']):
        for col,wl in enumerate([380,450,550,780]):
            try:
                sd,s90,ssource,s=wang_spectral_inputs(th,ph,zs,az,a.turbidity,wl,convention)
                w=complex_pattern(th,ph,zs,az,a.neutral_separation_deg)
                dop=abs(w)*(th*horizon_factor(th)+(np.pi/2-th)*s)*turbidity_multiplier(a.turbidity)
                valid=np.isfinite(dop)&(dop>=0)&(dop<=1)&np.isfinite(sd)&(sd>0)
                item=dict(chromaticity_convention=convention,wavelength_nm=wl,
                    zenith_xy=zenith_chromaticity(zs,a.turbidity,convention),
                    min_spectral_s=float(np.nanmin(s)),max_spectral_s=float(np.nanmax(s)),
                    min_raw_dop=float(np.nanmin(dop)),max_raw_dop=float(np.nanmax(dop)),
                    nonphysical_pixels=int((~valid).sum()),sky_pixels=int(sky.sum()),
                    reference_s90=float(s90),reference_ssource=float(ssource))
                arr=np.full(theta.shape,np.nan);arr[sky]=np.where(valid,dop,np.nan)
                art=axes[row,col].imshow(arr,cmap='turbo',vmin=0,vmax=1)
                axes[row,col].set_title(f'{convention}\n{wl}nm; invalid {item["nonphysical_pixels"]}')
            except ValueError as error:
                item=dict(chromaticity_convention=convention,wavelength_nm=wl,error=str(error))
                axes[row,col].text(.5,.5,str(error),ha='center',va='center',wrap=True,transform=axes[row,col].transAxes)
            axes[row,col].axis('off');results.append(item)
    fig.suptitle('Full printed eq17 spectral audit: nonphysical pixels masked and counted, not clipped\nNo lunar fit; reference ray and chromaticity convention explicit')
    fig.subplots_adjust(left=.02,right=.91,bottom=.02,top=.84,hspace=.28,wspace=.08)
    cax=fig.add_axes([.93,.15,.015,.6]);fig.colorbar(art,cax=cax,label='DoLP (0..1)')
    a.output.mkdir(parents=True)
    fig.savefig(a.output/'spectral_coefficient_audit.png',dpi=170);plt.close(fig)
    # AoP comparison diagnostic with consistent [-90,90] scale.
    fig,axes=plt.subplots(1,2,figsize=(10,5))
    for ax,(name,fn) in zip(axes,[('Printed eq4/5',complex_pattern),('Source-meridian spin2 comparison',source_meridian_pattern)]):
        angle=np.full(theta.shape,np.nan);angle[sky]=np.rad2deg(.5*np.angle(fn(th,ph,zs,az,a.neutral_separation_deg)))
        art=ax.imshow(angle,cmap='twilight',vmin=-90,vmax=90);ax.set_title(name);ax.axis('off');fig.colorbar(art,ax=ax,shrink=.7,label='AoP (degrees)')
    fig.suptitle('Coordinate-convention audit, not validation against measured skies');fig.tight_layout()
    fig.savefig(a.output/'aop_coordinate_audit.png',dpi=170);plt.close(fig)
    summary=dict(status='MODEL_AUDIT_NOT_VALIDATED_LUNAR_SIMULATION',parameters=vars(a)|{'output':str(a.output)},
        covariance_checks=covariance_audit(),spectral_checks=results,script_sha256=digest(__file__),
        coefficient_note='Corrected separate x/y matrices from Preetham author manuscript Appendix A.2; printed duplicate preserved',
        limitations=['Eq17 may yield nonphysical values for explicit reference choices; failures are not hidden',
        'AoP basis must be resolved before fixed-eye dynamics interpreted',
        'CIE daylight basis is not a validated lunar spectrum',
        '380..780nm only; no ultraviolet or near-infrared extrapolation'])
    (a.output/'wang_spectral_audit.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
