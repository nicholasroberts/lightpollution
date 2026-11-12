"""Test a 60deg full eye field at zenith and at the DoLP-band elevation.

Angular spherical caps use nominal equisolid lens mapping. Pixels contribute
with equal weight, as instructed for the historical eye model. Band position
is chosen by mean local DoLP, never by PRC or agreement with reference. Azimuth
is image-relative: zero right, +90 down in original TIFF coordinates. The
azimuth of a field centre is distinct from optimized polarization-axis angle.
"""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reconstruct_spider_eye_contrast import eye_result,digest


def sky_vectors(shape,center,r90):
    y,x=np.indices(shape,dtype=float)
    dx=(x-center[0])/r90;dy=(y-center[1])/r90
    r=np.hypot(dx,dy)
    theta=2*np.arcsin(np.clip(r/np.sqrt(2),0,1))
    az=np.arctan2(dy,dx)
    return np.array([np.sin(theta)*np.cos(az),np.sin(theta)*np.sin(az),np.cos(theta)]),r<=1


def cap(vectors,sky,zenith_deg,azimuth_deg):
    z,a=np.deg2rad([zenith_deg,azimuth_deg])
    direction=np.array([np.sin(z)*np.cos(a),np.sin(z)*np.sin(a),np.cos(z)])
    dot=np.einsum('i,ijk->jk',direction,vectors)
    return sky & (dot>=np.cos(np.deg2rad(30)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run=a.run.resolve();out=a.output.resolve()
    if out.exists():p.error('Use a fresh output folder')
    if out==run or run in out.parents or out in run.parents:p.error('Output must be separate from run')
    report=json.loads((run/'rotation_analysis.json').read_text())
    center=report['reference_center_xy_px'];r90=report['radius_60_px']*np.sqrt(2)
    arrays={}
    for method in ['linear_stokes','legacy_dop_aop']:
        with np.load(run/f'{method}.npz',allow_pickle=False) as z:
            arrays[method]=[z[k] for k in ['I','Q','U']]
    shape=arrays['linear_stokes'][0].shape
    # Downsample only the position-selection scan; final azimuth results use all pixels.
    stride=max(1,int(np.ceil(max(shape)/500)))
    sampled=[x[::stride,::stride] for x in arrays['linear_stokes']]
    v,sky=sky_vectors(sampled[0].shape,np.array(center)/stride,r90/stride)
    i,q,u=sampled
    with np.errstate(divide='ignore',invalid='ignore'):d=np.hypot(q,u)/i
    physical=np.isfinite(d)&(i>0)&(d<=1+1e-10)
    scan=[]
    for zenith in range(0,61,5):
        for azimuth in ([0] if zenith==0 else range(0,360,10)):
            mask=cap(v,sky,zenith,azimuth)
            valid=mask&physical
            fraction=float(valid.sum()/mask.sum()) if mask.any() else 0
            mean=float(np.mean(d[valid])) if fraction>=.95 else None
            scan.append(dict(zenith_deg=zenith,elevation_deg=90-zenith,azimuth_deg=azimuth,
                mean_dop=mean,physical_fraction=fraction))
    eligible=[s for s in scan if s['mean_dop'] is not None]
    if not eligible:raise ValueError('No field with >=95% physical pixels')
    band=max(eligible,key=lambda s:s['mean_dop'])
    print(f'Band selected by mean DoLP: elevation {band["elevation_deg"]}deg, image azimuth {band["azimuth_deg"]}deg',flush=True)
    v,sky=sky_vectors(shape,center,r90)
    rows=[]
    for zenith,azimuth,label in [(0,0,'zenith')]+[(band['zenith_deg'],az,'band_elevation_scan') for az in range(0,360,5)]:
        mask=cap(v,sky,zenith,azimuth)
        for method,(i,q,u) in arrays.items():
            with np.errstate(divide='ignore',invalid='ignore'):d=np.hypot(q,u)/i
            physical=mask&np.isfinite(d)&(i>0)&(d<=1+1e-10)
            fraction=float(physical.sum()/mask.sum()) if mask.any() else 0
            for mode in ['stokes_intensity','normalized_direction']:
                result=eye_result(i,q,u,mask,response_mode=mode,invalid_policy='exclude')
                rows.append(dict(field=label,method=method,response_mode=mode,
                    elevation_deg=90-zenith,azimuth_deg=azimuth,mean_dop=float(np.mean(d[physical])),
                    physical_fraction=fraction,prc=result['prc'],max_axis_deg=result['max_axis_deg_mod_180']))
    out.mkdir(parents=True)
    with (out/'azimuth_scan.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    selected=[r for r in rows if r['field']=='zenith' or r['azimuth_deg']==band['azimuth_deg']]
    payload=dict(status='directional_field_test_not_historical_validation',site_id=report['site_id'],
        selection='Maximum mean physical DoLP on 5deg zenith / 10deg azimuth grid; direct Stokes; >=95% physical pixels',
        selected_band=band,zenith_and_band_results=selected,position_scan=scan,
        angular_radius_deg=30,full_field_deg=60,image_relative_azimuth='0 right, 90 down in original TIFF; not compass bearing',
        source_report_sha256=digest(run/'rotation_analysis.json'),script_sha256=digest(__file__),
        source_array_sha256={m:digest(run/f'{m}.npz') for m in arrays},raw_images_modified=False,
        limitations=['Equisolid projection and transferred optical centre remain approximate',
        'Global image-plane polarization axes retained; full curved-eye angular geometry not introduced',
        'Residual alignment errors reused',
        'Band location is a computational mean-DoLP definition, not proof of historical viewing direction',
        'No terrain mask: ground within a field is not automatically identified',
        'Equal pixel contributions retained; no solid-angle weighting'])
    (out/'view_direction_summary.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    i,q,u=arrays['linear_stokes']
    with np.errstate(divide='ignore',invalid='ignore'):d=np.hypot(q,u)/i
    im=axes[0].imshow(d,cmap='turbo',vmin=0,vmax=1)
    for z,az,color,label in [(0,0,'white','Zenith'),(band['zenith_deg'],band['azimuth_deg'],'black','Band')]:
        mask=cap(v,sky,z,az)
        axes[0].contour(mask,levels=[.5],colors=color,linewidths=1.2)
        axes[0].plot([],[],color=color,label=label)
    axes[0].legend();axes[0].set_title('DoLP and 60deg full angular fields');axes[0].axis('off');fig.colorbar(im,ax=axes[0],shrink=.7)
    for method,color in [('linear_stokes','tab:blue'),('legacy_dop_aop','tab:orange')]:
        chosen=[r for r in rows if r['field']=='band_elevation_scan' and r['method']==method and r['response_mode']=='stokes_intensity']
        axes[1].plot([r['azimuth_deg'] for r in chosen],[r['prc'] for r in chosen],label=method,color=color)
        baseline=next(r['prc'] for r in rows if r['field']=='zenith' and r['method']==method and r['response_mode']=='stokes_intensity')
        axes[1].axhline(baseline,color=color,ls='--',alpha=.7)
    axes[1].axvline(band['azimuth_deg'],color='black',ls=':',label='DoLP-selected band direction')
    axes[1].set(xlabel='Image-relative field-centre azimuth (deg)',ylabel='Eye-to-eye log contrast',
        title=f'Field-centre elevation {band["elevation_deg"]}deg; dashed = zenith')
    axes[1].legend(fontsize=8);fig.tight_layout();fig.savefig(out/'band_vs_azimuth.png',dpi=170);plt.close(fig)
    print(f'Wrote view_direction_summary.json, azimuth_scan.csv and band_vs_azimuth.png in {out}',flush=True)

if __name__=='__main__':main()
