"""Run a declared reduced Wang angular/horizon pilot with additive ALAN.

--demo is an elevation sweep, NOT a real lunar ephemeris. For actual lunar
positions supply --start-utc/--end-utc and verified location. No photo fitting,
HPS/LED spectral conversion, or full eq17 claim is made by this pilot.
"""
from pathlib import Path
from datetime import datetime,timedelta,timezone
import argparse,csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wang2016_sky import natural_stokes,add_unpolarized,lunar_track
from reconstruct_spider_eye_contrast import eye_result,digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--demo',action='store_true')
    p.add_argument('--start-utc');p.add_argument('--end-utc')
    p.add_argument('--latitude',type=float);p.add_argument('--longitude',type=float)
    p.add_argument('--step-minutes',type=int,default=15)
    p.add_argument('--neutral-separation-deg',type=float,required=True)
    p.add_argument('--turbidity',type=float,default=3.)
    p.add_argument('--replace-turbidity-with-alan',action='store_true',help='Set natural amplitude multiplier to 1 and model dilution only via explicit ALAN')
    a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    if a.step_minutes<=0:p.error('Step must be positive')
    if a.demo:
        elevation=np.arange(5,86,5);azimuth=np.full(elevation.shape,30.)
        labels=[f'elevation_{e:g}' for e in elevation];kind='synthetic_source_elevation_sweep_not_lunar_track'
    else:
        if None in (a.start_utc,a.end_utc,a.latitude,a.longitude):p.error('Provide verified explicit UTC range and location, or --demo')
        start=datetime.fromisoformat(a.start_utc.replace('Z','+00:00'));end=datetime.fromisoformat(a.end_utc.replace('Z','+00:00'))
        if start.tzinfo is None or end.tzinfo is None:p.error('Times need explicit UTC timezone (+00:00 or Z)')
        if start.utcoffset()!=timedelta(0) or end.utcoffset()!=timedelta(0):p.error('Use UTC timestamps, not local clock times')
        if not timedelta(0)<end-start<=timedelta(days=2):p.error('UTC interval must be positive and at most two days')
        times=[start+timedelta(minutes=m) for m in range(0,int((end-start).total_seconds()/60)+1,a.step_minutes)]
        labels=[t.astimezone(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S') for t in times]
        elevation,azimuth=lunar_track(labels,a.latitude,a.longitude);kind='lunar_geometry_solar_model_extension'
    # Equal-area/equisolid pixel grid; sky azimuth zero on +x and +90 on +y.
    y,x=np.indices((201,201),dtype=float);dx=(x-100)/99;dy=(y-100)/99
    radial=np.hypot(dx,dy);sky=radial<.999
    theta=2*np.arcsin(np.clip(radial/np.sqrt(2),0,1));phi=np.arctan2(dy,dx)
    # Evaluate only upper hemisphere; arrays outside sky masked for eye_result.
    theta=np.where(sky,theta,0);phi=np.where(sky,phi,0)
    mask=sky&(theta<=np.deg2rad(30))
    rows=[];snapshots=[]
    for label,el,az in zip(labels,elevation,azimuth):
        print(f'{label}: source elevation {el:.2f}deg, azimuth {az:.2f}deg',flush=True)
        if el<=0:continue
        i,q,u=natural_stokes(theta,phi,np.deg2rad(90-el),np.deg2rad(az),
            a.neutral_separation_deg,np.ones_like(theta),turbidity=a.turbidity,mode='eq6_horizon',
            depolarization_multiplier=1. if a.replace_turbidity_with_alan else None)
        curves={}
        for addition in [0.,.5,1.,3.]:
            ii,qq,uu=add_unpolarized(i,q,u,addition)
            result=eye_result(ii,qq,uu,mask,invalid_policy='error')
            contrast=[]
            for heading in range(0,180,2):
                absorption_p=10**(-.4);absorption_t=10**(-.04)
                b=1-(absorption_p+absorption_t)/2;gain=(absorption_t-absorption_p)/2
                angle=np.deg2rad(heading)
                mean=b*ii[mask].sum();mod=-gain*(qq[mask].sum()*np.cos(2*angle)+uu[mask].sum()*np.sin(2*angle))
                c=float(np.log((mean+mod)/(mean-mod)));contrast.append(c)
                rows.append(dict(frame=label,source_elevation_deg=float(el),source_azimuth_deg=float(az),
                    alan_relative_intensity=addition,heading_deg=heading,signed_log_contrast=c,
                    maximum_abs_contrast=result['prc']))
            curves[addition]=contrast
        snapshots.append((label,el,az,q,u,curves))
    if not rows:raise ValueError('Moon/source below horizon throughout interval')
    a.output.mkdir(parents=True)
    with (a.output/'heading_responses.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for addition in [0.,.5,1.,3.]:
        values=[next(r['maximum_abs_contrast'] for r in rows if r['frame']==s[0] and r['alan_relative_intensity']==addition) for s in snapshots]
        axes[0].plot([s[1] for s in snapshots],values,'o-',label=f'ALAN / natural I = {addition:g}')
    axes[0].set(xlabel='Source elevation (deg)',ylabel='Maximum absolute eye log contrast',title='Fixed zenith field: 60deg full width');axes[0].legend(fontsize=8)
    snap=snapshots[len(snapshots)//2]
    for addition,curve in snap[5].items():axes[1].plot(range(0,180,2),curve,label=f'ALAN {addition:g}')
    axes[1].set(xlabel='Body heading relative to model axes (deg)',ylabel='Signed eye log contrast',title=f'Example: source elevation {snap[1]:.1f}deg');axes[1].legend(fontsize=8)
    fig.suptitle('Wang eq6 + eq8 reduced pilot; uniform relative intensity; not full spectral model')
    fig.tight_layout();fig.savefig(a.output/'wang_spider_pilot.png',dpi=170);plt.close(fig)
    (a.output/'provenance.json').write_text(json.dumps(dict(status='REDUCED_MODEL_PILOT_NOT_EMPIRICAL_VALIDATION',
        track_kind=kind,parameters=vars(a)|{'output':str(a.output)},model='Wang eqs1/3/4/5/6/7/8; no eq17 spectral factor',
        paper_doi='10.1088/2040-8978/18/6/065601',script_sha256=digest(__file__),
        limitations=['Uniform natural I=1 retained at every time: lunar brightness/phase not simulated',
        'ALAN is relative unpolarized intensity, not VIIRS radiance or HPS/LED conversion',
        'Solar model extended to lunar positions provisionally',
        'Neutral separation explicitly supplied, not recovered from paper',
        'No original photograph fit; no timestamp validation; no terrain mask',
        'Paper stereographic AoP convention retained; image/lab basis correspondence unverified']),indent=2)+'\n')

if __name__=='__main__':main()
