"""Camera-band integrated PRC comparison of checksum-locked Sun/Moon fits.

No refit or analyzer-zero calibration. Uses the same included pixels for all
models; reports orientation-invariant amplitudes separately from native axes.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def integrated_prc(I,Q,U,mask,parallel=.4,perpendicular=.04):
    """Established absorbance response; sum Stokes before receptor contrast."""
    if not np.isfinite([parallel,perpendicular]).all() or not 0<=perpendicular<parallel:
        raise ValueError('Require finite 0 <= perpendicular < parallel absorbance')
    mask=np.asarray(mask,dtype=bool)
    if I.shape!=Q.shape or I.shape!=U.shape or I.shape!=mask.shape or not mask.any():
        raise ValueError('Stokes/mask shapes differ or support is empty')
    values=[np.asarray(v,dtype=float)[mask] for v in [I,Q,U]]
    if not all(np.isfinite(v).all() for v in values):raise ValueError('Nonfinite included Stokes')
    i,q,u=values
    if np.any(i<=0) or np.any(np.hypot(q,u)>i*(1+1e-10)):
        raise ValueError('Nonphysical included Stokes')
    ti,tq,tu=map(float,[i.sum(),q.sum(),u.sum()])
    b=1-(10**(-parallel)+10**(-perpendicular))/2
    a=(10**(-perpendicular)-10**(-parallel))/2
    k=a*np.hypot(tq,tu)/(b*ti)
    if not 0<=k<1:raise ValueError('Nonpositive receptor response')
    return dict(sum_I=ti,sum_Q=tq,sum_U=tu,pixel_count=int(mask.sum()),
        integrated_dolp=float(np.hypot(tq,tu)/ti),modulation_fraction=float(k),
        maximum_absolute_prc=float(2*np.arctanh(k)),
        continuous_heading_median_absolute_prc=float(2*np.arctanh(k/np.sqrt(2))),
        maximum_axis_native_deg=float(np.rad2deg(.5*np.arctan2(tu,tq))%180))


def curve(result,angle_deg,relative_to_peak=False):
    phase=0 if relative_to_peak else result['maximum_axis_native_deg']
    return 2*np.arctanh(result['modulation_fraction']*np.cos(2*np.deg2rad(np.asarray(angle_deg)-phase)))


def run(a):
    from audit_spider_lunar_aop import block_mean,rays,fixed_moon_in_camera
    from fit_spider_twilight_moon import component,mixture,moon_plane_reflection
    fitpath=a.fit/'twilight_moon_fit.json';arraypath=a.fit/'twilight_moon_fit_arrays.npz'
    fit=json.loads(fitpath.read_text());name=fit['site'];base=a.reference_lock.parent
    if digest(a.reference_lock)!=fit['source_lock_sha256']:raise ValueError('Fit and lock differ')
    lock=json.loads(a.reference_lock.read_text())['sites'][name]
    if not lock.get('analyzer_stokes_convention'):raise ValueError('Explicit signed analyzer convention required')
    path=base/name/'registration/linear_stokes.npz'
    checks=json.loads((base/'reference_sha256.json').read_text())
    for f in [a.reference_lock,path]:
        if checks.get(str(f.relative_to(base)))!=digest(f):raise ValueError('Locked checksum mismatch')
    if digest(path)!=fit['stokes_sha256']:raise ValueError('Fit and measured Stokes differ')
    with np.load(path,allow_pickle=False) as z:I,Q,U=[block_mean(z[k],a.bin_size) for k in ['I','Q','U']]
    with np.load(arraypath,allow_pickle=False) as z:arrays={k:np.array(z[k]) for k in z.files}
    with np.errstate(divide='ignore',invalid='ignore'):observed=np.hypot(Q,U)/I
    if observed.shape!=arrays['measured_dolp'].shape or not np.allclose(observed,arrays['measured_dolp'],equal_nan=True):
        raise ValueError('Bin size or measured fit arrays differ')
    yy,xx=np.indices(I.shape);xy=np.stack([xx*a.bin_size+(a.bin_size-1)/2,yy*a.bin_size+(a.bin_size-1)/2],axis=-1)
    center=np.asarray(lock['reference_center_xy_px']);radius=lock['rim_radius_px']
    view=rays(xy,center,radius)
    moon=fixed_moon_in_camera(lock['astronomical_anchor'],lock['camera_to_north_xy'])
    solar=fit['solar'];sun=fixed_moon_in_camera(dict(moon_elevation_deg=solar['sun_elevation_deg'],
        moon_compass_azimuth_deg=solar['sun_compass_azimuth_deg']),lock['camera_to_north_xy'])
    if fit['handedness_hypothesis']=='reversed':
        reflection=moon_plane_reflection(center,lock['moon_pixel_xy'])
        moon=np.r_[reflection@moon[:2],moon[2]];sun=np.r_[reflection@sun[:2],sun[2]]
    elif fit['handedness_hypothesis']!='original':raise ValueError('Unknown parity')
    m=component(view,moon,fit['fixed_pmax']);s=component(view,sun,fit['fixed_pmax'])
    mixed=mixture(m,s,*fit['models']['mixed']['parameters'])
    models={'moon':(m[0]+fit['models']['moon']['parameters'][0],m[1]),
            'sun':(s[0]+fit['models']['sun']['parameters'][0],s[1]),'mixed':mixed}
    stokes={'measured':(I,Q,U)}
    for key,(intensity,z) in models.items():
        if not np.allclose(abs(z)/intensity,arrays[key+'_dolp'],equal_nan=True):raise ValueError('Saved model differs: '+key)
        stokes[key]=(intensity,z.real,z.imag)
    for key,value in zip(['mixed_I','mixed_Q','mixed_U'],stokes['mixed']):
        if not np.allclose(value,arrays[key],equal_nan=True):raise ValueError('Saved mixed Stokes differ')
    circle=np.linalg.norm(xy-center,axis=-1)<lock['eye_radius_px']
    fit_support=arrays['train']|arrays['withheld']
    support=circle&fit_support
    for i,q,u in stokes.values():
        support&=np.isfinite(i)&np.isfinite(q)&np.isfinite(u)&(i>0)&(np.hypot(q,u)<=i*(1+1e-10))
    if not support.any():raise ValueError('No common physical central-field support')
    regions={'common_central_field':support,'withheld_blocks_only':support&arrays['withheld']}
    results={region:{key:integrated_prc(*values,mask,a.parallel,a.perpendicular)
        for key,values in stokes.items()} for region,mask in regions.items()}
    for values in results.values():
        reference=values['measured']
        for key,v in values.items():
            v['maximum_prc_minus_measured']=v['maximum_absolute_prc']-reference['maximum_absolute_prc']
            v['maximum_prc_ratio_to_measured']=(v['maximum_absolute_prc']/reference['maximum_absolute_prc']
                if reference['maximum_absolute_prc']>0 else None)
    report=dict(site=name,camera_channel=fit['camera_channel'],status='camera_band_PRC_comparison_not_production_spectral_model',
        results=results,bin_size=a.bin_size,full_central_circle_pixels=int(circle.sum()),
        common_central_pixels=int(support.sum()),common_circle_coverage_fraction=float(support.sum()/circle.sum()),
        fitted_support_pixels_outside_central_circle=int((fit_support&~circle).sum()),
        absorbance_parallel=a.parallel,absorbance_perpendicular=a.perpendicular,
        prc_definition='ln(R1/R2); orthogonal receptors; R1=b*sum(I)+a*(sum(Q)*cos(2theta)+sum(U)*sin(2theta)); R2 reverses modulation',
        spatial_weighting='equal included binned image pixels; provisional equisolid projection',
        analyzer_convention=lock['analyzer_stokes_convention'],handedness_hypothesis=fit['handedness_hypothesis'],
        source_lock_sha256=digest(a.reference_lock),source_stokes_sha256=digest(path),
        source_fit_sha256=digest(fitpath),source_fit_arrays_sha256=digest(arraypath),script_sha256=digest(__file__),
        fit_or_zero_calibration_applied=False,
        limitations=['Primary comparison integrates the common included portion of the central 60-degree field; inspect coverage.',
        'Withheld-block integration samples a fragmented field, not a second full-eye measurement.',
        'Primary PRC uses regions used to fit DoLP; it is a derived comparison, not independent validation.',
        'Maximum and continuous heading-median absolute PRC are invariant to one constant analyzer-zero rotation.',
        'Native-axis curves retain unknown measured analyzer zero; they are not compass-heading comparisons.',
        'Peak-relative curves compare amplitudes only; each has its own zero, with no fitted angular offset.',
        'Measured and modeled intensity weighting are separate: model weights come from saved DoLP-only fits.',
        'Camera-plane receptor proxy uses existing absorbances 0.4/0.04 and ideal lens transport; not calibrated spider retinal transport.',
        'Camera-band result; no pigment spectral integration, calibrated ALAN spectrum or absolute source radiometry.',
        'Registration geometry, camera parity and off-axis analyzer calibration retain their stated uncertainties.'])
    a.output.mkdir(parents=True,exist_ok=False)
    (a.output/'prc_comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    angles=np.arange(360)
    with (a.output/'prc_orientation_curves.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['region','model','angle_deg','signed_prc_native_axis_unknown_zero','absolute_prc_relative_to_own_peak'])
        for region,values in results.items():
            for key,v in values.items():
                for theta,native,peak in zip(angles,curve(v,angles),abs(curve(v,angles,True))):
                    writer.writerow([region,key,int(theta),float(native),float(peak)])
    np.savez_compressed(a.output/'prc_support.npz',common_central_field=support,withheld_blocks_only=regions['withheld_blocks_only'])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    values=results['common_central_field'];labels=list(values);positions=np.arange(len(labels));colors=['black','C0','C1','C2']
    axes[0].bar(positions-.18,[values[k]['maximum_absolute_prc'] for k in labels],width=.36,label='Maximum |PRC|')
    axes[0].bar(positions+.18,[values[k]['continuous_heading_median_absolute_prc'] for k in labels],width=.36,label='Heading-median |PRC|')
    axes[0].set_xticks(positions);axes[0].set_xticklabels(labels);axes[0].set_ylabel('Integrated log receptor contrast');axes[0].legend()
    for key,color in zip(labels,colors):axes[1].plot(angles[:181],abs(curve(values[key],angles[:181],True)),label=key,color=color)
    axes[1].set_xlabel('Receptor axis relative to each pattern\'s own peak (deg)')
    axes[1].set_ylabel('|PRC|; amplitude comparison only');axes[1].legend()
    fig.suptitle(name+' — '+fit['camera_channel']+' camera band; matched central-field support\n'
        f"Coverage {report['common_circle_coverage_fraction']:.1%}; Stokes integrated first; no analyzer-zero or sky fit")
    fig.savefig(a.output/'prc_comparison.png',dpi=160);plt.close(fig)
    print(json.dumps(dict(site=name,coverage=report['common_circle_coverage_fraction'],results=values),indent=2))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--fit',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bin-size',type=int,default=8)
    p.add_argument('--parallel',type=float,default=.4)
    p.add_argument('--perpendicular',type=float,default=.04)
    a=p.parse_args()
    if a.bin_size<1:p.error('Bin size must be positive')
    try:run(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))


if __name__=='__main__':main()
