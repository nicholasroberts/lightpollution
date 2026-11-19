"""Unfitted common-pair AoP check of saved fixed-geometry mixture fits.

Reports both stored-U and conjugated-U conventions without choosing either.
No analyzer zero, rotation, source weight or transmission is fitted here.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit_phase(z):
    norm=np.abs(z)
    return np.divide(z,norm,out=np.full(z.shape,np.nan+1j*np.nan),where=norm>1e-12)


def common_mask(withheld,observed,phases,dolps,minimum=.05):
    mask=np.asarray(withheld,dtype=bool)&np.isfinite(observed)&(observed>minimum)
    for name,phase in phases.items():
        mask&=np.isfinite(phase)&np.isfinite(dolps[name])&(dolps[name]>minimum)
    return mask


def change_error(measured,predicted,mask,axis,lag):
    if lag<1 or lag>=mask.shape[axis]:raise ValueError('Invalid pair spacing')
    lo=[slice(None)]*2;hi=lo.copy();lo[axis]=slice(None,-lag);hi[axis]=slice(lag,None)
    lo=tuple(lo);hi=tuple(hi)
    m=measured[hi]*np.conj(measured[lo]);p=predicted[hi]*np.conj(predicted[lo])
    valid=mask[hi]&mask[lo]&np.isfinite(m)&np.isfinite(p)
    return np.where(valid,np.rad2deg(.5*np.angle(m*np.conj(p))),np.nan)


def metrics(measured,predicted,mask,lag):
    values=np.concatenate([change_error(measured,predicted,mask,axis,lag).ravel() for axis in [0,1]])
    values=values[np.isfinite(values)]
    if not len(values):return dict(pair_count=0,median_absolute_deg=None,p95_absolute_deg=None,rms_deg=None)
    return dict(pair_count=int(len(values)),median_absolute_deg=float(np.median(abs(values))),
        p95_absolute_deg=float(np.percentile(abs(values),95)),rms_deg=float(np.sqrt(np.mean(values**2))))


def run(a):
    from audit_spider_lunar_aop import block_mean,rays,fixed_moon_in_camera
    from fit_spider_twilight_moon import component,mixture,moon_plane_reflection
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report_path=a.fit/'twilight_moon_fit.json';arrays_path=a.fit/'twilight_moon_fit_arrays.npz'
    fit=json.loads(report_path.read_text());site_name=fit['site']
    if digest(a.reference_lock)!=fit['source_lock_sha256']:raise ValueError('Fit and registration lock differ')
    lock=json.loads(a.reference_lock.read_text())['sites'][site_name]
    stokes_path=a.reference_lock.parent/site_name/'registration/linear_stokes.npz'
    checks=json.loads((a.reference_lock.parent/'reference_sha256.json').read_text())
    for path in [a.reference_lock,stokes_path]:
        if digest(path)!=checks[str(path.relative_to(a.reference_lock.parent))]:raise ValueError('Locked checksum mismatch')
    if digest(stokes_path)!=fit['stokes_sha256']:raise ValueError('Fit and camera Stokes differ')
    with np.load(stokes_path) as z:I,Q,U=[block_mean(z[k],a.bin_size) for k in ['I','Q','U']]
    with np.load(arrays_path) as z:arrays={k:np.array(z[k]) for k in z.files}
    with np.errstate(divide='ignore',invalid='ignore'):observed=np.hypot(Q,U)/I
    if observed.shape!=arrays['measured_dolp'].shape or not np.allclose(observed,arrays['measured_dolp'],equal_nan=True):raise ValueError('Bin size or stored measured arrays differ from fit')
    y,x=np.indices(I.shape);xy=np.stack([x*a.bin_size+(a.bin_size-1)/2,y*a.bin_size+(a.bin_size-1)/2],axis=-1)
    center=np.asarray(lock['reference_center_xy_px']);radius=lock['rim_radius_px']
    view=rays(xy,center,radius);rotation=lock['camera_to_north_xy']
    moon=fixed_moon_in_camera(lock['astronomical_anchor'],rotation)
    solar=fit['solar'];sun=fixed_moon_in_camera(dict(moon_elevation_deg=solar['sun_elevation_deg'],moon_compass_azimuth_deg=solar['sun_compass_azimuth_deg']),rotation)
    if fit['handedness_hypothesis']=='reversed':
        reflection=moon_plane_reflection(center,lock['moon_pixel_xy'])
        moon=np.r_[reflection@moon[:2],moon[2]];sun=np.r_[reflection@sun[:2],sun[2]]
    elif fit['handedness_hypothesis']!='original':raise ValueError('Unknown parity hypothesis')
    m=component(view,moon,fit['fixed_pmax']);s=component(view,sun,fit['fixed_pmax'])
    r,b=fit['models']['mixed']['parameters'];mixed=mixture(m,s,r,b)
    models={'moon':(m[0]+fit['models']['moon']['parameters'][0],m[1]),
            'sun':(s[0]+fit['models']['sun']['parameters'][0],s[1]),'mixed':mixed}
    phases={};dolps={}
    for name,(intensity,z) in models.items():
        dolps[name]=abs(z)/intensity;phases[name]=unit_phase(z)
        if not np.allclose(dolps[name],arrays[name+'_dolp'],equal_nan=True):raise ValueError('Recreated model differs from saved fit: '+name)
    if not all(np.allclose(v,arrays[k],equal_nan=True) for k,v in [('mixed_I',mixed[0]),('mixed_Q',mixed[1].real),('mixed_U',mixed[1].imag)]):raise ValueError('Saved mixed Stokes differ')
    measured=unit_phase(Q+1j*U)
    mask=common_mask(arrays['withheld'],observed,phases,dolps)&np.isfinite(measured)
    if not mask.any():raise ValueError('No common AoP samples above threshold')
    predictions=dict(phases,constant=np.ones(I.shape,dtype=complex));result={}
    for label,phase in [('stored_U',measured),('conjugated_U_diagnostic',np.conj(measured))]:
        result[label]={str(lag*a.bin_size):{name:metrics(phase,predicted,mask,lag) for name,predicted in predictions.items()} for lag in a.lag_bins}
        for spacing,values in result[label].items():
            if len({v['pair_count'] for v in values.values()})!=1:raise ValueError('AoP comparison pairs differ')
    a.output.mkdir(parents=True,exist_ok=False)
    fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
    phase_panels=[('Measured AoP; unknown constant zero',measured),('Fixed-Moon AoP',phases['moon']),('Fixed mixed AoP',phases['mixed'])]
    for ax,(title,phase) in zip(axes[0],phase_panels):
        image=ax.imshow(np.where(mask,np.rad2deg(.5*np.angle(phase)),np.nan),cmap='twilight',vmin=-90,vmax=90)
        ax.set_title(title);fig.colorbar(image,ax=ax,shrink=.7,label='Degrees modulo 180');ax.axis('off')
    lag=a.lag_bins[0]
    for ax,name in zip(axes[1],['moon','sun','mixed']):
        v=change_error(measured,phases[name],mask,1,lag)
        image=ax.imshow(v,cmap='RdBu_r',vmin=-15,vmax=15);ax.set_title(name+': horizontal AoP-change error')
        ax.axis('off');fig.colorbar(image,ax=ax,shrink=.7,label='Degrees; display ±15')
    fig.suptitle(site_name+' — camera coordinates; identical withheld pairs for every model\nNo refit, zero-angle adjustment or convention selection. Horizontal panels: '+str(lag*a.bin_size)+' px spacing.')
    fig.savefig(a.output/'common_pair_aop.png',dpi=160);plt.close(fig)
    summary=dict(site=site_name,camera_channel=fit['camera_channel'],handedness_hypothesis=fit['handedness_hypothesis'],
        common_withheld_pixel_count=int(mask.sum()),minimum_measured_and_each_model_dolp=.05,
        statistics=result,bin_size=a.bin_size,source_fit_sha256=digest(report_path),source_fit_arrays_sha256=digest(arrays_path),
        source_lock_sha256=digest(a.reference_lock),source_stokes_sha256=digest(stokes_path),script_sha256=digest(__file__),
        fit_or_correction_applied=False,limitations=['Stored-U and conjugated-U are reported, not selected by sky agreement.',
        'Spatial AoP changes cancel a constant analyzer zero; absolute celestial AoP remains uncalibrated.',
        'The shared mask depends on all three saved models and restricts comparison to common support.',
        'Pair differences combine horizontal and vertical directions; pairs are correlated, not independent observations.',
        'Geometric parity, projection, analyzer tilt transport and instrument polarization remain separate uncertainties.'])
    (a.output/'common_pair_aop.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary['statistics'],indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--fit',type=Path,required=True)
    p.add_argument('--bin-size',type=int,default=8)
    p.add_argument('--lag-bins',type=int,nargs='+',default=[4,8])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.bin_size<1:p.error('Bin size must be positive')
    try:run(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))


if __name__=='__main__':main()
