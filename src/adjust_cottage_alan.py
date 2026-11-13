"""User-authorized Cottage scenario adjustment; never read experimental DoLP.

Set mean local, proxy-weighted sky DoLP in zenith30deg cap at capture time for
HPS at a specified scenario level. Apply same physical source scaling to LED.
"""
from pathlib import Path
import argparse,json
import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import brentq
from model_spider_nightly_prc import body_field
from wang2016_sky import berry_normalized_pattern,horizon_factor
from spider_spectral_lunar import ks91_brightness


def mean_dop(pnatural,intensity,alpha,ni,ai):
    return float(np.mean(pnatural*intensity*ni/(intensity*ni+alpha*ai)))


def target_addition(pnatural,intensity,target,ni=1.,ai=1.):
    if not np.isfinite(target) or not 0<target<1:raise ValueError('Target DoLP must be in(0,1)')
    baseline=mean_dop(pnatural,intensity,0,ni,ai)
    if target>baseline:raise ValueError(f'Target {target:g} exceeds no-ALAN mean {baseline:g}; ALAN cannot increase DoLP')
    if abs(target-baseline)<1e-12:return 0.
    upper=1.
    while mean_dop(pnatural,intensity,upper,ni,ai)>target:
        upper*=2
        if upper>1e9:raise ValueError('Failed to bracket ALAN addition')
    return float(brentq(lambda v:mean_dop(pnatural,intensity,v,ni,ai)-target,0,upper,xtol=1e-12))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--configs',type=Path,required=True);p.add_argument('--spectra',type=Path,required=True)
    p.add_argument('--moon-anchors',type=Path,required=True);p.add_argument('--target-dop',type=float,required=True)
    p.add_argument('--target-level',type=float,default=5.);p.add_argument('--atmospheric-factor',type=float,default=.8)
    p.add_argument('--extinction-v',type=float,default=.2);p.add_argument('--neutral-separation-deg',type=float,default=24.)
    a=p.parse_args()
    if not np.isfinite(a.target_level) or a.target_level<=0 or not np.isfinite(a.atmospheric_factor) or not 0<a.atmospheric_factor<=1:p.error('Positive level and atmospheric factor in(0,1] required')
    paths=[a.configs/(s+'_site_alan.json') for s in ['HPS','LED_mid']];configs=[json.loads(f.read_text()) for f in paths]
    if any(c['normalization_site']!='Cottage' for c in configs):raise ValueError('Cottage must be common natural reference')
    anchor=json.loads(a.moon_anchors.read_text())['sites']['Cottage'];el=anchor['moon_elevation_deg'];az=anchor['moon_compass_azimuth_deg']
    rays,_=body_field(0,0,8192);th=np.arccos(rays[:,2]);ph=np.arctan2(rays[:,1],rays[:,0]);zs=np.deg2rad(90-el);ps=np.deg2rad(90-az)
    rho=np.arccos(np.clip(np.cos(th)*np.cos(zs)+np.sin(th)*np.sin(zs)*np.cos(ph-ps),-1,1))
    intensity=ks91_brightness(th,el,rho,0.,a.extinction_v)/ks91_brightness(0.,el,zs,0.,a.extinction_v)
    spectrum=np.genfromtxt(a.spectra/'HPS_visible_input.csv',delimiter=',',names=True);w=spectrum['wavelength_nm']
    n=spectrum['natural']/trapezoid(spectrum['natural'],w);art=spectrum['alan']/trapezoid(spectrum['alan'],w)
    scale=spectrum['dop_scale'] if 'dop_scale' in spectrum.dtype.names else np.ones(len(w))
    ni=trapezoid(n*spectrum['sensitivity'],w);ai=trapezoid(art*spectrum['sensitivity'],w);nq=trapezoid(n*spectrum['sensitivity']*scale,w)
    pn=a.atmospheric_factor*abs(berry_normalized_pattern(th,ph,zs,ps,a.neutral_separation_deg))*horizon_factor(th)*nq/ni
    alpha=target_addition(pn,intensity,a.target_dop,ni,ai)
    old=configs[0]['alan_scale_by_site']['Cottage'];new=alpha/(a.target_level*configs[0]['spectral_band_factor'])
    if new<old-1e-12:raise ValueError('Requested target would reduce Cottage ALAN; this adjustment is increase-only')
    adjustment=dict(status='USER_REQUESTED_SCENARIO_NOT_FIT_TO_EXPERIMENT',site='Cottage',target_mean_local_dop=a.target_dop,reference_spectrum='HPS',reference_level=a.target_level,field='Zenith-centred60deg full field; mean of local DoLP, not DoLP of spatially integrated Stokes',capture_time=anchor['local_time'],atmospheric_factor=a.atmospheric_factor,extinction_v=a.extinction_v,neutral_separation_deg=a.neutral_separation_deg,old_cottage_scale=old,new_cottage_scale=new,alan_increase_factor=None if old==0 else new/old,hps_local_capture_alan_ratio_at_target=alpha,achieved_mean_dop=mean_dop(pn,intensity,alpha,ni,ai),no_alan_mean_dop=mean_dop(pn,intensity,0,ni,ai),effective_oxford_over_cottage_source_ratio=configs[0]['alan_scale_by_site']['oxfordbussinesspark']/new,samples=len(th),note='Cottage site-contrast assumption explicitly superseded; Oxford scale unchanged; LED retains DNB-overlap factor; no experimental map read')
    for path,config in zip(paths,configs):
        config['alan_scale_by_site']['Cottage']=new;config['cottage_adjustment']=adjustment;config['status']='USER_ADJUSTED_COTTAGE_SCENARIO'
        config['site_contrast']['proxy_superseded_for_cottage']=True
        path.write_text(json.dumps(config,indent=2)+'\n')
    (a.configs/'cottage_adjustment.json').write_text(json.dumps(adjustment,indent=2)+'\n')
    print(f'Cottage HPS level{a.target_level:g}: mean DoLP{adjustment["achieved_mean_dop"]:.4f}; ALAN increase{adjustment["alan_increase_factor"]:.3g}x; Oxford unchanged',flush=True)

if __name__=='__main__':main()
