"""Subfigure 1C: measured museum and saved lunar/solar/combined patterns.

Reuse the mixture fit; estimate one global analyzer-basis offset on training pixels.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cmcrameri import cm
from reconstruct_spider_eye_contrast import digest
from model_paper1_lunar_baseline import lunar_component
from spider_fisheye_registration import rays
from fit_paper1_solar_lunar_angle import source_vector


def analyzer_basis_alignment(arr):
    """One spin-2 offset from combined model to measured analyzer basis."""
    measured=arr['measured_Q']+1j*arr['measured_U']
    model=arr['combined_Q']+1j*arr['combined_U']
    mask=arr['fit_mask'].astype(bool)
    for name in ['measured','lunar','solar','combined']:
        z=arr[name+'_Q']+1j*arr[name+'_U']
        mask&=np.isfinite(z)&(abs(z)>1e-10)
    held=arr['withheld_pixel_partition'].astype(bool)
    train=mask&~held;test=mask&held
    if min(train.sum(),test.sum())<20:
        raise ValueError('Insufficient common training/withheld analyzer samples')
    delta=(measured[train]/abs(measured[train]))*np.conj(model[train]/abs(model[train]))
    mean=delta.mean();coherence=float(abs(mean))
    if coherence<1e-6:
        raise ValueError('Global analyzer offset is undefined: incoherent training angles')
    offset=.5*np.angle(mean)
    spin=np.exp(2j*offset)
    diagnostics={}
    for name in ['lunar','solar','combined']:
        z=arr[name+'_Q']+1j*arr[name+'_U']
        diagnostics[name]={}
        for label,selection in [('training',train),('withheld',test)]:
            residual=np.rad2deg(.5*np.angle(measured[selection]*np.conj(z[selection]*spin)))
            raw=np.rad2deg(.5*np.angle(measured[selection]*np.conj(z[selection])))
            diagnostics[name][label]=dict(pixel_count=int(selection.sum()),
                median_absolute_deg=float(np.median(abs(residual))),
                rms_deg=float(np.sqrt(np.mean(residual**2))),
                p95_absolute_deg=float(np.percentile(abs(residual),95)),
                unaligned_rms_deg=float(np.sqrt(np.mean(raw**2))))
    return spin,dict(model_to_measured_offset_deg=float(np.rad2deg(offset)),
        training_circular_coherence=coherence,absolute_AoP_residuals=diagnostics,
        method='Half-angle of mean unit spin-2 measured/combined phase difference on training pixels only; equal weights. Same offset applied to every theoretical component.',
        calibration_status='Pattern-estimated global analyzer-basis registration, not independently measured polarizer calibration',
        sky_rotation_fitted=False,mixture_refitted=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fit',type=Path,required=True)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--north',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    rp=a.fit/'solar_lunar_fit.json';ap=a.fit/'solar_lunar_fit_arrays.npz'
    fit=json.loads(rp.read_text());site=fit['site']
    for path,key in [(a.reference_lock,'source_lock_sha256'),(a.north/'north_stokes.npz','source_north_arrays_sha256'),(a.north/'north_rotation.json','source_north_report_sha256')]:
        if digest(path)!=fit[key]:p.error('Source changed: '+str(path))
    for name,key in [('wang2016_sky.py','kernel_sha256'),('model_paper1_lunar_baseline.py','component_module_sha256')]:
        if digest(Path(__file__).with_name(name))!=fit[key]:p.error('Model code changed since fit: '+name)
    lock=json.loads(a.reference_lock.read_text())['sites'][site]
    c=np.asarray(lock['reference_center_xy_px']);radius=lock['rim_radius_px']
    with np.load(ap,allow_pickle=False) as z:arr={k:np.array(z[k]) for k in z.files}
    with np.load(a.north/'north_stokes.npz',allow_pickle=False) as z:
        observed_moon=np.array(z['moon_xy'])
    x,y=arr['x_px'],arr['y_px']
    view=rays(np.stack([x,2*c[1]-y],-1),c,radius)
    anchor=fit['astronomical_anchor'];solar=fit['solar']
    moon=source_vector(anchor['moon_elevation_deg'],anchor['moon_compass_azimuth_deg'])
    sun=source_vector(solar['sun_elevation_deg'],solar['sun_compass_azimuth_deg'])
    for name,source in [('lunar',moon),('solar',sun)]:
        i,z=lunar_component(view,source,fit['neutral_separation_deg'])
        for key,v in [('I',i),('Q',z.real),('U',z.imag)]:
            if not np.allclose(arr[name+'_'+key],v,equal_nan=True):p.error('Saved component differs: '+name)
    f=fit['solar_fraction']
    for k in ['I','Q','U']:
        expected=(1-f)*arr['lunar_'+k]+f*arr['solar_'+k]
        if not np.allclose(expected,arr['combined_'+k],equal_nan=True):p.error('Saved mixture differs')
    spin,alignment=analyzer_basis_alignment(arr)
    names=['measured','lunar','solar','combined']
    labels=['Measured green; north-up','Lunar only','Solar only',f'Best combined\nSolar {100*f:.1f}% / lunar {100*(1-f):.1f}%']
    disc=np.hypot(x-c[0],y-c[1])<.985*radius
    field=np.hypot(x-c[0],y-c[1])<lock['eye_radius_px']
    fig,axes=plt.subplots(2,4,figsize=(14,7.8),layout='constrained')
    stats={}
    for col,(name,label) in enumerate(zip(names,labels)):
        i,q,u=[arr[name+'_'+k] for k in ['I','Q','U']]
        if name!='measured':
            rotated=(q+1j*u)*spin
            q,u=rotated.real,rotated.imag
        with np.errstate(divide='ignore',invalid='ignore'):d=np.hypot(q,u)/i
        valid=disc&np.isfinite(d)&(i>0)&(d>=0)&(d<=1)
        angle=np.rad2deg(.5*np.arctan2(u,q))
        angle_valid=valid&(np.hypot(q,u)>1e-10)
        for row,(values,mask,cmap,lo,hi) in enumerate([(d,valid,cm.batlow,0,1),(angle,angle_valid,cm.romaO,-90,90)]):
            ax=axes[row,col]
            im=ax.pcolormesh(x,y,np.where(mask,values,np.nan),shading='nearest',cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
            ax.add_patch(plt.Circle(c,lock['eye_radius_px'],fill=False,color='white',lw=.8))
            ax.plot(*observed_moon,'kx',ms=4)
            mr=radius*np.sqrt(1-moon[2])
            predicted=c+np.array([moon[0],-moon[1]])*mr/max(np.linalg.norm(moon[:2]),1e-12)
            ax.plot(*predicted,'r+',ms=4)
            for tx,ty,t in [(.5,.99,'N'),(.01,.5,'E'),(.5,.01,'S'),(.99,.5,'W')]:
                ax.text(tx,ty,t,transform=ax.transAxes,ha='center',va='center',fontsize=8)
            ax.set(xlim=(c[0]-1.07*radius,c[0]+1.07*radius),ylim=(c[1]+1.07*radius,c[1]-1.07*radius),aspect='equal')
            ax.axis('off')
            if row==0:ax.set_title(label,fontsize=10)
        v=d[field&valid]
        stats[name]=dict(mean_dolp=float(v.mean()),spatial_sd=float(v.std(ddof=1)),included_pixels=int(v.size),
            support='Central60degree binned grid, physical pixels for this field. Measured statistics are north-resampled diagnostics, not replacement camera-frame publication values.')
    fig.colorbar(axes[0,0].collections[0],ax=axes[0,:].tolist(),shrink=.65,label='DoLP')
    fig.colorbar(axes[1,0].collections[0],ax=axes[1,:].tolist(),shrink=.65,label='AoP (degrees modulo 180)')
    fig.suptitle(site+' — fixed-geometry solar/lunar comparison\nNo ALAN; one global analyzer-basis offset fitted on training pixels',fontsize=13)
    fig.text(.012,.985,'C',ha='left',va='top',fontsize=18,fontweight='bold')
    a.output.mkdir(parents=True)
    for ext in ['png','pdf','svg']:fig.savefig(a.output/('fig01c_'+site+'.'+ext),dpi=300)
    plt.close(fig)
    report=dict(site=site,solar_fraction=f,lunar_fraction=1-f,solar=solar,astronomical_anchor=anchor,
        source_fit_report_sha256=digest(rp),source_fit_arrays_sha256=digest(ap),script_sha256=digest(__file__),
        central_statistics=stats,colormaps=dict(DoLP='batlow',AoP='romaO'),
        fit_changed=False,background_added=0,absolute_analyzer_zero_fitted=True,analyzer_basis_alignment=alignment,
        limitations=['Observed north-up data, not ALAN-corrected measurements.',
        'Full-disc model is context outside the fit mask; no claim that foreground is predicted.',
        'One global analyzer offset estimated from training pattern, not independent instrument calibration; assess absolute residuals on withheld pixels.',
        'Mixture inferred on central spatial pairs, not independent-night validation.',
        'Below-horizon solar component is an exploratory geometric extension.'])
    (a.output/'analyzer_basis_alignment.json').write_text(json.dumps(dict(site=site,source_fit_report_sha256=digest(rp),**alignment),indent=2)+'\n')
    (a.output/'figure_values.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(alignment,indent=2))
    print('Saved',a.output)


if __name__=='__main__':main()
