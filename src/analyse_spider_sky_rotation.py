"""Register camera-roll sky exposures; write diagnostic, unreviewed PRC previews.

No TIFFs are altered. Uses the accepted 5% Cottage circle fits, rotates each
exposure about its own centre into image 1 coordinates, without resizing.
Scene clockwise roll is the source-sampling angle, not an eye-axis rotation.
Analyzer convention remains nominal 0,45,90,135 in common coordinates; reversing
all analyzer angles changes U sign but not maximized integrated contrast.
This is the historical global image-plane model, not a calibrated full fisheye
polarimetric ray model. Registration does not normalize image brightness.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np
from scipy import ndimage, optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from reconstruct_spider_eye_contrast import (tiff_channel, digest,
    combine_replicates, run, GROUPS, circle_mask, eye_result)


def polar_features(image, center, radius, angles=720):
    # Local spatial contrast reduces broad analyzer-dependent brightness changes.
    image=np.asarray(image,dtype=float)
    detail=image-ndimage.gaussian_filter(image,6)
    radial=np.linspace(.15*radius,.985*radius,200)
    theta=np.arange(angles)*2*np.pi/angles
    x=center[0]+radial[:,None]*np.cos(theta)
    y=center[1]+radial[:,None]*np.sin(theta)
    polar=ndimage.map_coordinates(detail,[y,x],order=1,mode='constant',cval=0)
    polar-=polar.mean(axis=1,keepdims=True)
    norms=np.sqrt(np.mean(polar**2,axis=1,keepdims=True))
    # Limit leverage of near-empty annuli while retaining stars and skyline.
    polar/=np.maximum(norms,np.percentile(norms,30)+1e-12)
    return polar


def estimate_roll(reference, moving):
    """Positive angle: scene in moving image is clockwise from reference."""
    n=reference.shape[1]
    if moving.shape!=reference.shape: raise ValueError('Polar dimensions differ')
    denom=np.linalg.norm(reference)*np.linalg.norm(moving)
    if denom<=1e-12: raise ValueError('No structural signal for rotation')
    corr=np.fft.irfft(np.conj(np.fft.rfft(reference,axis=1))*
        np.fft.rfft(moving,axis=1),n=n,axis=1).sum(axis=0)/denom
    best=int(np.argmax(corr))
    # Subsample phase interpolation, bounded around the global peak.
    def score(shift):
        shifted=np.fft.irfft(np.fft.rfft(moving,axis=1)*
            np.exp(2j*np.pi*np.fft.rfftfreq(n)*shift),n=n,axis=1)
        return float(np.sum(reference*shifted)/denom)
    fit=optimize.minimize_scalar(lambda s:-score(s),bounds=(best-1,best+1),method='bounded')
    angle=(fit.x*360/n+180)%360-180
    distance=abs((np.arange(n)-fit.x+n/2)%n-n/2)*360/n
    alternate=float(corr[distance>15].max())
    return dict(scene_roll_clockwise_deg=float(angle),correlation=score(fit.x),
        alternate_peak_outside_15deg=alternate,peak_margin=score(fit.x)-alternate)


def warp_to_reference(image, source_center, target_center, roll_deg):
    y,x=np.indices(image.shape,dtype=float)
    dx=x-target_center[0]; dy=y-target_center[1]
    a=np.deg2rad(roll_deg)
    sx=source_center[0]+np.cos(a)*dx-np.sin(a)*dy
    sy=source_center[1]+np.sin(a)*dx+np.cos(a)*dy
    return ndimage.map_coordinates(np.asarray(image,dtype=float),[sy,sx],
        order=1,mode='constant',cval=np.nan,prefilter=False)


def refine_registration(reference, moving, reference_center, source_center, radius, initial_roll):
    """Local rotation search with FFT translation, using inner-sky structure.

    Returns a source-coordinate centre correction. Intensities for scientific
    output remain untouched. Translation is estimated independently of rim fit.
    """
    step=max(1,int(np.ceil(max(reference.shape)/600)))
    ref=reference[::step,::step]
    mov=moving[::step,::step]
    rc=np.array(reference_center)/step; sc=np.array(source_center)/step
    y,x=np.indices(ref.shape)
    r=np.hypot(x-rc[0],y-rc[1])
    # Exclude horizon; the biological integration region drives registration.
    window=np.clip((radius/step-r)/8,0,1)
    def features(im):
        im=np.nan_to_num(im)
        f=im-ndimage.gaussian_filter(im,2)
        # Bound bright Moon/hot-pixel leverage; retain signed local structure.
        limit=np.percentile(abs(f)[window>.5],99.5)+1e-12
        return np.clip(f,-limit,limit)*window
    rf=features(ref)
    maxshift=max(3,int(radius/step*.15))
    def evaluate(angle):
        warped=warp_to_reference(mov,sc,rc,angle)
        mf=features(warped)
        corr=np.fft.ifft2(np.conj(np.fft.fft2(rf))*np.fft.fft2(mf)).real
        yy=np.arange(-maxshift,maxshift+1); xx=yy.copy()
        local=corr[np.ix_(yy%ref.shape[0],xx%ref.shape[1])]
        iy,ix=np.unravel_index(np.argmax(local),local.shape)
        delta=np.array([xx[ix],yy[iy]],dtype=float)
        shifted=ndimage.shift(np.nan_to_num(warped),[-delta[1],-delta[0]],order=1)
        ff=features(shifted)
        score=np.sum(rf*ff)/(np.linalg.norm(rf)*np.linalg.norm(ff)+1e-12)
        return float(score),delta
    candidates=np.arange(initial_roll-12,initial_roll+12.001,1.)
    coarse=[evaluate(a) for a in candidates]
    best_angle=float(candidates[np.argmax([v[0] for v in coarse])])
    candidates=np.arange(best_angle-1,best_angle+1.001,.1)
    fine=[evaluate(a) for a in candidates]
    k=int(np.argmax([v[0] for v in fine])); angle=float(candidates[k])
    score,delta=fine[k]
    # Continuous refinement starts from the global local-search maximum.
    def cost(p):
        a,dx,dy=p
        c=np.cos(np.deg2rad(a)); t=np.sin(np.deg2rad(a))
        center=sc+np.array([c*dx-t*dy,t*dx+c*dy])
        mf=features(warp_to_reference(mov,center,rc,a))
        return -np.sum(rf*mf)/(np.linalg.norm(rf)*np.linalg.norm(mf)+1e-12)
    fit=optimize.minimize(cost,[angle,*delta],method='Powell',
        bounds=[(angle-1,angle+1),(delta[0]-2,delta[0]+2),(delta[1]-2,delta[1]+2)],
        options=dict(maxiter=30,xtol=.01,ftol=1e-5))
    angle,dx,dy=fit.x
    c=np.cos(np.deg2rad(angle)); t=np.sin(np.deg2rad(angle))
    correction=step*np.array([c*dx-t*dy,t*dx+c*dy])
    return dict(scene_roll_clockwise_deg=float(angle),inner_sky_correlation=float(-fit.fun),
        source_center_correction_xy_px=correction.tolist(),
        refined_source_center_xy_px=(np.array(source_center)+correction).tolist(),
        translation_limit_reached=bool(abs(dx)>=maxshift-1 or abs(dy)>=maxshift-1),
        rotation_search_limit_reached=bool(abs(angle-initial_roll)>11)), correction



def refine_multistart(reference, moving, reference_center, source_center, radius,
                      initial_roll, nominal):
    """Compare independent starts by structural image correlation, never PRC."""
    seeds=[float(initial_roll),float(nominal),float(-nominal)]
    unique=[]
    for seed in seeds:
        if all(abs((seed-old+180)%360-180)>2 for old in unique):
            unique.append(seed)
    trials=[]
    for seed in unique:
        fit,correction=refine_registration(reference,moving,reference_center,
                                           source_center,radius,seed)
        trials.append((fit,correction,seed))
    best=max(trials,key=lambda t:t[0]['inner_sky_correlation'])
    fit,correction,seed=best
    fit['selected_seed_clockwise_deg']=seed
    fit['multistart_trials']=[dict(seed_clockwise_deg=t[2],
        scene_roll_clockwise_deg=t[0]['scene_roll_clockwise_deg'],
        inner_sky_correlation=t[0]['inner_sky_correlation'],
        rotation_search_limit_reached=t[0]['rotation_search_limit_reached'])
        for t in trials]
    return fit,correction


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--geometry',type=Path,required=True)
    p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--site-id',default='s14')
    p.add_argument('--multistart',action='store_true',help='Retry polar and both signed nominal rotation starts; select by structural correlation')
    a=p.parse_args()
    source=a.folder.resolve(); target=a.output.resolve()
    if source==target or source in target.parents or target in source.parents:
        p.error('Output must be separate from TIFF input folder')
    if target.exists(): p.error('Output already exists; use a fresh folder')
    geometry=json.loads(a.geometry.read_text())
    records={Path(r['file']).name:r for r in geometry['records']}
    groups=GROUPS if len(records)==12 else [[1,2,3,4]] if len(records)==4 else None
    if groups is None: p.error('Expected twelve or four geometry records')
    channel=geometry['channel']
    images={}; centres={}; radii={}; inputs=[]
    for num in sorted(n for g in groups for n in g):
        name=f'{num}.tif'; rec=records[name]
        path=source/name
        image,meta=tiff_channel(path,channel)
        if meta['sha256']!=rec['sha256']: raise ValueError(f'Geometry hash mismatch: {name}')
        f=next(f for f in rec['fits'] if f['threshold_fraction']==.05 and 'error' not in f)
        images[num]=image; centres[num]=[f['center_x_px'],f['center_y_px']]
        radii[num]=f['rim_radius_px']; inputs.append(meta)
    shape=images[1].shape
    if any(im.shape!=shape for im in images.values()): raise ValueError('Image dimensions differ')
    target.mkdir(parents=True)
    ref=polar_features(images[1],centres[1],radii[1])
    aligned={}; rotations=[]
    radius60=radii[1]/np.sqrt(2) # historical report compatibility only
    eye_radius=radii[1]*np.sin(np.deg2rad(15))/np.sin(np.deg2rad(45))
    fig,axes=plt.subplots(int(np.ceil(len(images)/4)),4,figsize=(16,4*np.ceil(len(images)/4)),squeeze=False)
    ref_edge=ndimage.gaussian_gradient_magnitude(images[1],1.5)
    for (num,image),ax in zip(images.items(),axes.flat):
        estimate=estimate_roll(ref,polar_features(image,centres[num],radii[num]))
        if num==1: estimate['scene_roll_clockwise_deg']=0.
        initial_roll=estimate['scene_roll_clockwise_deg']
        group_index=next(k for g in groups for k,n in enumerate(g) if n==num)
        nominal=[0,45,90,135][group_index]
        if num==1:
            refinement=dict(scene_roll_clockwise_deg=0.,inner_sky_correlation=1.,
                source_center_correction_xy_px=[0.,0.],refined_source_center_xy_px=centres[1],
                translation_limit_reached=False,rotation_search_limit_reached=False)
        else:
            if a.multistart:
                refinement,_=refine_multistart(images[1],image,centres[1],centres[num],radius60,initial_roll,nominal)
            else:
                refinement,_=refine_registration(images[1],image,centres[1],centres[num],radius60,initial_roll)
        estimate['initial_polar_roll_clockwise_deg']=initial_roll
        estimate['refinement']=refinement
        roll=refinement['scene_roll_clockwise_deg']
        estimate['scene_roll_clockwise_deg']=roll
        aligned[num]=warp_to_reference(image,refinement['refined_source_center_xy_px'],centres[1],roll)
        group_index=next(k for g in groups for k,n in enumerate(g) if n==num)
        nominal=[0,45,90,135][group_index]
        error=min(abs((roll-s*nominal+180)%360-180) for s in [-1,1])
        estimate.update(filename=f'{num}.tif',nominal_analyzer_deg=nominal,
            distance_to_either_signed_nominal_deg=error,
            flags=(['rotation differs from either nominal sign by >5deg'] if error>5 else [])+
                  (['weak or ambiguous structural match'] if estimate['correlation']<.2 or estimate['peak_margin']<.03 else []))
        if refinement['inner_sky_correlation']<.2:
            estimate['flags'].append('refined inner-sky match remains weak')
        if refinement['translation_limit_reached'] or refinement['rotation_search_limit_reached']:
            estimate['flags'].append('registration search limit reached')
        rotations.append(estimate)
        im=aligned[num]
        # Display normalization is confined to this visual diagnostic.
        def display(x):
            lo,hi=np.nanpercentile(x,[10,99]); return np.nan_to_num(np.clip((x-lo)/(hi-lo+1e-12),0,1))
        edge=ndimage.gaussian_gradient_magnitude(np.nan_to_num(im),1.5)
        rgb=np.stack([display(ref_edge),display(edge),display(edge)],axis=-1)
        ax.imshow(rgb); ax.add_patch(Circle(centres[1],eye_radius,fill=False,color='yellow',lw=.7))
        ax.set_title(f'{num}.tif: roll {roll:.2f}deg; inner corr {refinement["inner_sky_correlation"]:.2f}')
        ax.axis('off')
    for ax in list(axes.flat)[len(images):]: ax.axis('off')
    fig.suptitle('ALIGNMENT REVIEW: red = reference edges; cyan = aligned exposure edges\nWhite overlap = agreement; yellow = 60deg full field (30deg radius); display scales only')
    fig.tight_layout(rect=[0,0,1,.94]); fig.savefig(target/'rotation_overlay.png',dpi=150); plt.close(fig)
    report=dict(status='unreviewed_registration_and_prc_preview',site_id=a.site_id,
        geometry_sha256=digest(a.geometry),script_sha256=digest(__file__),inputs=inputs,
        multistart_enabled=a.multistart,geometry_threshold=.05,reference_center_xy_px=centres[1],radius_60_px=radius60,
        eye_angular_radius_deg=30,eye_full_field_deg=60,eye_radius_px=eye_radius,
        rotations=rotations,raw_images_modified=False,
        limitations=['Structural rotation estimation can be biased by polarization brightness gradients',
        'Circle centres initialize translations; inner-sky structure refines them; residual parallax, tilt or distortion not corrected',
        'Nominal analyzer angles assumed; registration angles are not automatically substituted for physical analyzer calibration',
        'Common global image-plane polarization axes follow historical model; no full ray-dependent Mueller calibration',
        'No exposure or gamma correction; TIFF radiometric assumptions remain those of historical pipeline'])
    (target/'rotation_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    # Two explicitly separate reconstructions; no historical value is overwritten.
    report['prc_diagnostics']=[]
    for method in ['linear_stokes','legacy_dop_aop']:
        i,q,u,qa=combine_replicates(([aligned[n] for n in g] for g in groups),method)
        for mode in ['stokes_intensity','normalized_direction']:
            for factor in [.95,1.,1.05]:
                mask=circle_mask(shape,centres[1],eye_radius*factor)
                result=eye_result(i,q,u,mask,response_mode=mode,invalid_policy='exclude')
                report['prc_diagnostics'].append(dict(method=method,response_mode=mode,
                    radius_factor=factor,prc=result['prc'],excluded_fraction=result['excluded_fraction']))
        mask=circle_mask(shape,centres[1],eye_radius)
        report.setdefault('analyzer_pair_closure',{})[method]=[
            dict(median_abs=float(np.nanmedian(abs(rec['closure'][mask]))),
                 p95_abs=float(np.nanpercentile(abs(rec['closure'][mask]),95))) for rec in qa]
        np.savez_compressed(target/f'{method}.npz',I=i,Q=q,U=u)
        config=dict(schema_version=1,site_id=a.site_id,geometry_reviewed=True,
            registration_reviewed=False,input=dict(kind='stokes_npz',path=f'{method}.npz'),
            fov=dict(angular_radius_deg=30,center_xy_px=centres[1],radius_px=eye_radius),
            absorbance_parallel=.4,absorbance_perpendicular=.04,
            response_mode='stokes_intensity',invalid_policy='exclude',
            processing_note=method,replicate_count=len(groups))
        manifest=target/f'{method}_manifest.json';manifest.write_text(json.dumps(config,indent=2)+'\n')
        try: run(manifest,target/f'{method}_preview',preview_only=True)
        except ValueError as e:
            (target/f'{method}_preview_error.txt').write_text(str(e)+'\n')
    (target/'rotation_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'Review {target}/rotation_overlay.png and rotation_analysis.json before accepting PRC previews')

if __name__=='__main__': main()
