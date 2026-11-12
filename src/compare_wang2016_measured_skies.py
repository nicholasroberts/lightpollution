"""Compare Cottage and Oxford Business Park Stokes patterns to Wang model.

Visible-source geometry is provisional until marked Moon candidate reviewed.
No exposure UTC, compass bearing, sky-only terrain mask or spectral ALAN ratio
is inferred. One constant AoP offset is estimated; DoLP amplitude is not fitted.
"""
from pathlib import Path
import argparse,json
import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from wang2016_sky import (berry_normalized_pattern,horizon_factor,alan_relative_from_depth)
from reconstruct_spider_eye_contrast import digest


def axial_offset_and_residual(observed,model,mask):
    """Angles are half-complex-phase; estimates only one global basis offset."""
    selected=mask&np.isfinite(observed)&np.isfinite(model)&(abs(observed)>.02)&(abs(model)>.02)
    if selected.sum()<10:raise ValueError('Too few polarized pixels for AoP comparison')
    offset=.5*np.angle(np.sum(observed[selected]*np.conj(model[selected])))
    residual=.5*np.angle(observed[selected]*np.conj(model[selected]*np.exp(2j*offset)))
    return float(offset),dict(pixels=int(selected.sum()),
        median_abs_deg=float(np.rad2deg(np.median(abs(residual)))),
        p95_abs_deg=float(np.rad2deg(np.quantile(abs(residual),.95))))


def source_candidate(intensity,sky):
    """Compact local brightness peak; explicitly not established Moon identity."""
    im=np.where(sky,np.nan_to_num(intensity),0)
    prominence=ndimage.gaussian_filter(im,2)-ndimage.gaussian_filter(im,18)
    permitted=ndimage.binary_erosion(sky,iterations=6)
    score=np.where(permitted,prominence,-np.inf)
    y,x=np.unravel_index(np.argmax(score),score.shape)
    return [float(x),float(y)],float(score[y,x])


def ephemeris_source_pixel(anchor,center,r90):
    if anchor.get('candidate_pixel_xy') is None:
        raise ValueError('No image bearing available for Moon anchor')
    elevation=float(anchor['moon_elevation_deg'])
    if not 0<elevation<90:raise ValueError('Moon elevation must be within (0,90)')
    vector=np.array(anchor['candidate_pixel_xy'],dtype=float)-np.array(center)
    length=np.linalg.norm(vector)
    if length<=0:raise ValueError('Moon image bearing cannot be inferred from image centre')
    radius=r90*np.sin(np.deg2rad((90-elevation)/2))/np.sin(np.pi/4)
    return (np.array(center)+radius*vector/length).tolist()


def compare_site(directory,source_xy,separation,turbidity=None,wavelength=None,anchor=None,alan_relative=0.):
    if turbidity is not None or wavelength is not None:
        raise ValueError("This comparison uses a turbidity-free angular/horizon model; spectral/turbidity arguments are unsupported")
    if not np.isfinite(alan_relative) or alan_relative<0:
        raise ValueError("ALAN relative intensity must be finite and nonnegative")
    report_path=directory/'rotation_analysis.json';report=json.loads(report_path.read_text())
    path=directory/'linear_stokes.npz'
    with np.load(path,allow_pickle=False) as z:i,q,u=[z[k] for k in ['I','Q','U']]
    center=np.array(report['reference_center_xy_px']);r90=report['radius_60_px']*np.sqrt(2)
    y,x=np.indices(i.shape);dx=(x-center[0])/r90;dy=(y-center[1])/r90;r=np.hypot(dx,dy)
    sky=r<.985;theta=2*np.arcsin(np.clip(r/np.sqrt(2),0,1));phi=np.arctan2(dy,dx)
    if anchor is not None:
        if anchor['source_report_sha256']!=digest(report_path):raise ValueError('Moon anchor belongs to another alignment run')
        if source_xy is not None:raise ValueError('Choose either explicit Moon pixels or ephemeris anchors')
        source_xy=ephemeris_source_pixel(anchor,center,r90);prominence=None
        selection='ephemeris_elevation_with_reviewed_image_bearing'
    elif source_xy is None:source_xy,prominence=source_candidate(i,sky);selection='automatic_bright_source_candidate_UNVERIFIED'
    else:prominence=None;selection='user_supplied_source_pixel'
    delta=(np.array(source_xy)-center)/r90;sr=np.linalg.norm(delta)
    if sr>=1 or sr<0:raise ValueError('Source outside nominal sky circle')
    zs=2*np.arcsin(sr/np.sqrt(2));az=np.arctan2(delta[1],delta[0])
    if not 0<zs<np.pi/2:raise ValueError('Source must be strictly above horizon and below zenith')
    th=theta[sky];ph=phi[sky]
    w=berry_normalized_pattern(th,ph,zs,az,separation)
    pangular=abs(w)*horizon_factor(th)
    # I_total=I_natural*(1+a), Q/U unchanged: no atmospheric attenuation.
    pfull=pangular/(1+alan_relative)
    model=np.full(i.shape,np.nan+1j*np.nan,dtype=complex)
    physicalmodel=np.isfinite(pfull)&(pfull>=0)&(pfull<=1)
    model[sky]=np.where(physicalmodel,pfull*np.exp(1j*np.angle(w)),np.nan)
    with np.errstate(divide='ignore',invalid='ignore'):observed=(q+1j*u)/i
    field=sky&(theta<=np.deg2rad(30))
    physical=field&np.isfinite(observed)&(i>0)&(abs(observed)<=1)&np.isfinite(model)
    offset,residual=axial_offset_and_residual(observed,model,physical)
    model*=np.exp(2j*offset)
    dob=abs(observed[physical]);dmod=abs(model[physical])
    result=dict(site_id=report['site_id'],source_selection=selection,source_pixel_xy=source_xy,
        source_candidate_prominence=prominence,source_elevation_deg=float(90-np.rad2deg(zs)),
        source_azimuth_image_deg=float(np.rad2deg(az)%360),global_aop_offset_deg=float(np.rad2deg(offset)),
        aop_residual_in_zenith_field=residual,field_pixels=int(field.sum()),compared_pixels=int(physical.sum()),
        observed_mean_dop=float(np.mean(dob)),model_mean_dop=float(np.mean(dmod)),
        dop_mean_absolute_error=float(np.mean(abs(dob-dmod))),
        model_nonphysical_sky_pixels=int((~physicalmodel).sum()),
        angular_reduced_max_dop=float(np.max(pangular)),full_model_max_raw_dop=float(np.nanmax(pfull)),
        model_structure="Berry eq4.2 (-4 prefactor) + Wang horizon eq8; no atmospheric turbidity or spectral extension",
        alan_relative_intensity=float(alan_relative),alan_depth=float(np.log1p(alan_relative)),
        rotation_report_sha256=digest(report_path),stokes_sha256=digest(path),
        flagged_exposures=sum(bool(r['flags']) for r in report['rotations']))
    if anchor is not None:
        result['ephemeris_anchor']=anchor
        result['geometry_note']='Ephemeris elevation fixes radial coordinate; image candidate fixes bearing; level-camera assumption retained'
    return result,(i,observed,model,sky,field,center,r90)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--moon-anchors',type=Path,help='Reviewed moon_ephemeris_anchors.json')
    p.add_argument('--cottage-moon-xy',nargs=2,type=float)
    p.add_argument('--oxford-moon-xy',nargs=2,type=float)
    p.add_argument('--neutral-separation-deg',type=float,required=True)
    group=p.add_mutually_exclusive_group()
    group.add_argument('--alan-relative',type=float,default=None,help='I_ALAN/I_natural; default 0')
    group.add_argument('--alan-depth',type=float,help='ln(1+I_ALAN/I_natural); not atmospheric turbidity')
    a=p.parse_args();batch=a.batch.resolve();out=a.output.resolve()
    addition=float(alan_relative_from_depth(a.alan_depth)) if a.alan_depth is not None else (a.alan_relative or 0.)
    if not np.isfinite(addition) or addition<0:p.error('ALAN addition must be finite and nonnegative')
    if out.exists():p.error('Use fresh output directory')
    if out==batch or batch in out.parents or out in batch.parents:p.error('Output must be separate from existing batch')
    anchors=json.loads(a.moon_anchors.read_text())['sites'] if a.moon_anchors else {}
    results={};fig,axes=plt.subplots(2,4,figsize=(15,8))
    for row,(site,source) in enumerate([('Cottage',a.cottage_moon_xy),('oxfordbussinesspark',a.oxford_moon_xy)]):
        print(f'Comparing {site}',flush=True)
        stats,(i,observed,model,sky,field,center,r90)=compare_site(batch/site/'registration',source,
            a.neutral_separation_deg,anchor=anchors.get(site),alan_relative=addition)
        results[site]=stats
        for ax,arr,title,lo,hi,cmap in zip(axes[row],
            [abs(observed),abs(model),np.rad2deg(.5*np.angle(observed)),np.rad2deg(.5*np.angle(model))],
            ['Measured DoLP','Wang DoLP (no amplitude fit)','Measured AoP','Wang AoP (global offset fitted)'],
            [0,0,-90,-90],[1,1,90,90],['turbo','turbo','twilight','twilight']):
            cm=plt.get_cmap(cmap).copy();cm.set_bad('lightgray')
            arr=np.where(sky,arr,np.nan)
            art=ax.imshow(arr,cmap=cm,vmin=lo,vmax=hi)
            ax.add_patch(Circle(center,r90*np.sin(np.deg2rad(15))/np.sin(np.deg2rad(45)),fill=False,color='white',lw=.8))
            ax.plot(*stats['source_pixel_xy'],'x',color='black',ms=7)
            ax.set_title(f'{site}\n{title}',fontsize=10);ax.axis('off');fig.colorbar(art,ax=ax,shrink=.65)
    fig.suptitle(f'Measured versus Wang source-meridian model: Berry-normalized angular/horizon model; no turbidity; ALAN relative I = {addition:g}\nBlack cross = selected Moon geometry; white circle = 60deg full animal field')
    fig.tight_layout(rect=[0,0,1,.93]);out.mkdir(parents=True)
    fig.savefig(out/'cottage_oxford_model_comparison.png',dpi=170);plt.close(fig)
    payload=dict(status='PROVISIONAL_SOURCE_GEOMETRY_AND_PATTERN_COMPARISON',sites=results,
        parameters=vars(a)|{'batch':str(batch),'output':str(out),'moon_anchors':str(a.moon_anchors) if a.moon_anchors else None},script_sha256=digest(__file__),
        limitations=['Automatic bright source may be artificial light: verify each black cross before interpretation',
        'No UTC, Moon ephemeris or compass orientation inferred',
        'Reduced angular/horizon model: spectral eq17 omitted to avoid hidden atmospheric turbidity',
        'ALAN addition is proportional to natural intensity across the sky; not uniform absolute skyglow',
        'Berry unity-peak polarization intensity used as a DoLP envelope assumption, not radiometric derivation',
        'Berry minus sign rotates raw AoP by 90deg; one global fitted basis offset compensates in comparison',
        'No amplitude calibration or ALAN amount inferred from between-site differences',
        'No sky/terrain segmentation; primary comparison restricted to zenith 30deg radius',
        'Residual image registration and uncertain geometry reused',
        'Source-meridian convention is explicit; constant AoP offset is registration, not validation',
        'Cottage low ALAN does not mean ALAN-free; neutral separation remains an assumption'])
    (out/'cottage_oxford_model_summary.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print(f'Wrote comparison PNG and summary JSON in {out}',flush=True)

if __name__=='__main__':main()
