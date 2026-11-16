"""Register original exposures to 1.tif and freeze reproducible site coordinates.

Only static image landmarks determine exposure transforms. Astronomical Moon
coordinates remain fixed; an identified image Moon sets a display rotation.
Numerical acceptance does not replace review of the saved edge overlays.
"""
import argparse,json,math,shutil
from pathlib import Path
import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from reconstruct_spider_eye_contrast import digest,tiff_channel,combine_replicates,GROUPS

ROOT=Path(__file__).resolve().parents[1]
# Approximate image-source seeds from the 16 Nov camera-coordinate audit.
# Offsets are relative to reference centre, in units of the saved rim radius.
# They do not use astronomical elevation, DoLP or AoP to choose a source.
MOON_SEEDS={'minifactory':[.70,-.05],'museum':[-.57,-.06],
    'opensquare':[-.61,-.03],'southhorrington':[-.91,0.],
    'stannes':[.61,-.04],'WCpitch':[-.58,0.],'zoology':[.65,-.02]}
INHERITED={'Cottage','oxfordbussinesspark','Skillinge','slabhouse'}


def rigid_fit(source,target):
    """Fit target = R source + t, with unit scale and no reflection."""
    source=np.asarray(source,dtype=float);target=np.asarray(target,dtype=float)
    if source.shape!=target.shape or source.ndim!=2 or source.shape[1]!=2 or len(source)<2:
        raise ValueError('Need corresponding 2D points')
    a=source-source.mean(0);b=target-target.mean(0)
    u,_,vt=np.linalg.svd(a.T@b)
    d=np.eye(2);d[-1,-1]=np.linalg.det(vt.T@u.T)
    rotation=vt.T@d@u.T
    return np.column_stack([rotation,target.mean(0)-rotation@source.mean(0)])


def transform_points(points,matrix):
    return np.asarray(points)@matrix[:,:2].T+matrix[:,2]


def robust_rigid(source,target,tolerance=2.5,trials=1500):
    """Deterministic rigid RANSAC; no affine distortion or scale fitting."""
    source=np.asarray(source);target=np.asarray(target)
    if len(source)<8:raise ValueError('Fewer than eight training landmark matches')
    rng=np.random.default_rng(20261116);best=None
    for _ in range(trials):
        ids=rng.choice(len(source),2,replace=False)
        if np.linalg.norm(source[ids[0]]-source[ids[1]])<30:continue
        matrix=rigid_fit(source[ids],target[ids])
        error=np.linalg.norm(transform_points(source,matrix)-target,axis=1)
        mask=error<tolerance
        score=(int(mask.sum()),-float(np.median(error[mask])) if mask.any() else -np.inf)
        if best is None or score>best[0]:best=(score,mask)
    if best is None or best[1].sum()<8:raise ValueError('Insufficient consistent rigid landmarks')
    mask=best[1]
    for _ in range(3):
        matrix=rigid_fit(source[mask],target[mask])
        mask=np.linalg.norm(transform_points(source,matrix)-target,axis=1)<tolerance
        if mask.sum()<8:raise ValueError('Unstable rigid fit')
    return rigid_fit(source[mask],target[mask]),mask


def warp(image,source_to_reference,shape=None):
    """One inverse bilinear resampling; unsupported borders are NaN."""
    if isinstance(source_to_reference,dict):
        from spider_fisheye_registration import warp as fisheye_warp
        return fisheye_warp(image,source_to_reference,shape)
    shape=image.shape if shape is None else shape
    rotation=source_to_reference[:,:2];translation=source_to_reference[:,2]
    inverse=rotation.T
    y,x=np.indices(shape,dtype=float)
    sx=inverse[0,0]*(x-translation[0])+inverse[0,1]*(y-translation[1])
    sy=inverse[1,0]*(x-translation[0])+inverse[1,1]*(y-translation[1])
    return ndimage.map_coordinates(np.asarray(image,dtype=float),[sy,sx],order=1,
        mode='constant',cval=np.nan,prefilter=False)


def landmark_features(image,center,radius):
    import cv2
    y,x=np.indices(image.shape);r=np.hypot(x-center[0],y-center[1])
    valid=np.isfinite(image)&(image>0)
    # Include static skyline/buildings, exclude most inner sky and black rim.
    mask=valid&(r>.55*radius)&(r<1.03*radius)
    values=image[mask]
    if values.size<1000:raise ValueError('Insufficient outer-field pixels')
    lo,hi=np.percentile(values,[1,99.5])
    scaled=np.clip((np.nan_to_num(image)-lo)/max(hi-lo,1e-12),0,1)
    display=(255*np.log1p(30*scaled)/np.log(31)).astype(np.uint8)
    sift=cv2.SIFT_create(nfeatures=6000,contrastThreshold=.015,edgeThreshold=15)
    keys,desc=sift.detectAndCompute(display,(mask*255).astype(np.uint8))
    if desc is None or len(keys)<20:raise ValueError('Too few static SIFT landmarks')
    return np.array([k.pt for k in keys]),desc


def match_landmarks(reference,moving):
    import cv2
    rp,rd=reference;mp,md=moving;matcher=cv2.BFMatcher(cv2.NORM_L2)
    def good(a,b):
        return {pair[0].queryIdx:pair[0].trainIdx for pair in matcher.knnMatch(a,b,k=2)
                if len(pair)==2 and pair[0].distance<.72*pair[1].distance}
    forward=good(md,rd);reverse=good(rd,md)
    pairs=[(m,r) for m,r in forward.items() if reverse.get(r)==m]
    if len(pairs)<16:raise ValueError('Fewer than sixteen reciprocal landmark matches')
    return np.array([mp[m] for m,r in pairs]),np.array([rp[r] for m,r in pairs])


def register_landmarks(reference_features,moving_features,center,radius):
    source,target=match_landmarks(reference_features,moving_features)
    # Withhold every fourth correspondence before fitting or model selection.
    held=np.arange(len(source))%4==0
    matrix,inside=robust_rigid(source[~held],target[~held])
    error=np.linalg.norm(transform_points(source,matrix)-target,axis=1)
    validation=error[held];confirmed=validation<3.
    good_target=target[~held][inside]
    sectors=((np.arctan2(good_target[:,1]-center[1],good_target[:,0]-center[0])+np.pi)*8/(2*np.pi)).astype(int)%8
    span=float(np.linalg.norm(np.ptp(good_target,axis=0)))
    quality=dict(training_matches=int((~held).sum()),training_inliers=int(inside.sum()),
        training_median_px=float(np.median(error[~held][inside])),
        heldout_matches=int(held.sum()),heldout_inliers=int(confirmed.sum()),
        heldout_inlier_fraction=float(confirmed.mean()),
        heldout_median_all_matches_px=float(np.median(validation)),
        heldout_inlier_p95_px=float(np.percentile(validation[confirmed],95)) if confirmed.any() else None,
        occupied_reference_sectors=int(len(np.unique(sectors))),landmark_span_px=span)
    passed=(inside.sum()>=8 and confirmed.sum()>=4 and confirmed.mean()>=.6 and
        np.median(validation)<2. and len(np.unique(sectors))>=3 and span>.5*radius)
    quality['numerical_pass']=bool(passed)
    return matrix,quality,source,target,error


def moon_coordinate(image,center,radius,site,anchor,overrides):
    """Localize only around a visually nominated pixel, never around ephemeris."""
    supplied=overrides.get(site,{})
    if supplied.get('moon_pixel_xy') is not None:
        xy=np.asarray(supplied['moon_pixel_xy'],dtype=float)
        status='investigator_supplied_image_coordinate'
        if not supplied.get('note'):raise ValueError('Provide a note for supplied Moon identity')
    else:
        if site in INHERITED and anchor and anchor.get('candidate_pixel_xy') is not None:
            seed=np.asarray(anchor['candidate_pixel_xy'],dtype=float)
        elif site in MOON_SEEDS:seed=np.asarray(center)+radius*np.asarray(MOON_SEEDS[site])
        else:return None,dict(status='unresolved_identity_no_north_rotation')
        y,x=np.indices(image.shape);valid=np.isfinite(image)&(image>0)
        patch=valid&((x-seed[0])**2+(y-seed[1])**2<(.08*radius)**2)
        detail=ndimage.gaussian_filter(np.nan_to_num(image),2)-ndimage.gaussian_filter(np.nan_to_num(image),15)
        if not patch.any():return None,dict(status='empty_moon_seed_patch')
        iy,ix=np.unravel_index(np.argmax(np.where(patch,detail,-np.inf)),image.shape)
        if detail[iy,ix]<=0 or np.linalg.norm(np.array([ix,iy])-seed)>.07*radius:
            return None,dict(status='moon_peak_missing_or_near_patch_boundary')
        xy=np.array([ix,iy],dtype=float);status='audit_seed_raw_image_localization_provisional'
    if xy.shape!=(2,) or not np.isfinite(xy).all() or not(0<=xy[0]<image.shape[1] and 0<=xy[1]<image.shape[0]):
        raise ValueError('Invalid Moon pixel')
    rr=float(np.linalg.norm(xy-center))
    implied=90-float(np.rad2deg(2*np.arcsin(np.clip(rr/(radius*np.sqrt(2)),0,1))))
    meta=dict(status=status,pixel_xy=xy.tolist(),image_implied_elevation_deg=implied,
        note=supplied.get('note','Pixel selected from raw 1.tif near visually nominated source; identity remains provisional'))
    if anchor:
        meta['elevation_residual_deg']=implied-anchor['moon_elevation_deg']
        meta['elevation_consistent']=bool(abs(meta['elevation_residual_deg'])<=6)
    return xy,meta


def north_rotation(center,moon_pixel,azimuth):
    """Pure image rotation: observed Moon bearing -> fixed compass azimuth."""
    direction=np.asarray(moon_pixel)-center
    if np.linalg.norm(direction)<1:raise ValueError('Moon too close to zenith to determine bearing')
    observed=math.atan2(direction[1],direction[0])
    alpha=math.radians(azimuth-90)-observed
    c,s=math.cos(alpha),math.sin(alpha);rotation=np.array([[c,-s],[s,c]])
    return np.column_stack([rotation,np.asarray(center)-rotation@center]),float(np.rad2deg(alpha))


def edge_display(image):
    edge=ndimage.gaussian_gradient_magnitude(np.nan_to_num(image),1.2)
    hi=np.percentile(edge,99.5)
    return np.clip(edge/max(hi,1e-12),0,1)


def process_site(site,batch,out,anchor,overrides,saved_site=None,registration_model='rigid',linear_raw_manifest=None):
    oldpath=batch/site/'registration/rotation_analysis.json'
    old=json.loads(oldpath.read_text());center=np.asarray(old['reference_center_xy_px'],dtype=float)
    if saved_site is not None and saved_site.get('source_report_sha256')!=digest(oldpath):
        raise ValueError('Saved landmark lock uses a different original registration report')
    radius=float(old['radius_60_px'])*np.sqrt(2)
    oldinputs={Path(r['path']).name:r for r in old['inputs']}
    numbers=sorted(int(Path(n).stem) for n in oldinputs)
    groups=GROUPS if numbers==list(range(1,13)) else [[1,2,3,4]] if numbers==[1,2,3,4] else None
    if groups is None:raise ValueError('Unsupported explicit exposure grouping')
    channel=oldinputs['1.tif'].get('channel')
    raw_set=None;raw_frames={};geometry_transfer=None
    if linear_raw_manifest is not None:
        if saved_site is not None:raise ValueError('New RAW images require fresh landmarks, not --from-lock')
        raw_set=json.loads(Path(linear_raw_manifest).read_text())
        raw_channel=raw_set.get('camera_channel','blue')
        if raw_channel not in ('blue','green'):raise ValueError('Unsupported RAW camera channel')
        if raw_set.get('status')!='linear_raw_'+raw_channel+'_prepared_not_registered_complete':raise ValueError('Linear RAW conversion is incomplete')
        settings=raw_set.get('processing',{})
        if settings.get('method')!='visible_Bayer_'+raw_channel+'_samples_minus_camera_reported_black':raise ValueError('Unknown linear RAW preparation method')
        for key in ['white_balance_applied','gamma_applied','color_matrix_applied','display_profile_or_look_applied','auto_brightness_applied']:
            if settings.get(key) is not False:raise ValueError('RAW preparation must disable '+key)
        raw_frames={int(f['exposure']):f for f in raw_set['frames']}
        if len(raw_frames)!=len(raw_set['frames']) or 1 not in raw_frames:raise ValueError('Invalid RAW exposure list or missing reference')
        if not set(raw_frames)<=set(numbers):raise ValueError('RAW exposure absent from original batch')
        numbers=sorted(raw_frames);groups=raw_set.get('analyzer_groups')
        expected=[g for g in GROUPS if set(g)<=set(numbers)]
        if groups!=expected or not groups:raise ValueError('RAW analyzer groups must be exactly the available complete quartets')
        if set(raw_set.get('excluded_exposures',[]))!=set(range(1,13))-set(numbers):raise ValueError('RAW exclusions do not match available frames')
    images={};inputs=[]
    for n in numbers:
        record=oldinputs[f'{n}.tif'];path=Path(record['path'])
        if digest(path)!=record['sha256']:raise ValueError('Raw exposure hash changed: '+str(path))
        if raw_set is not None:
            from PIL import Image
            with Image.open(path) as original:
                timestamp=str(original.getexif().get_ifd(34665).get(36867,''))
            frame=raw_frames[n]
            if frame.get('datetime_original')!=timestamp:raise ValueError('RAW/TIFF capture timestamps differ for exposure '+str(n))
            path=Path(linear_raw_manifest).parent/f'{n}.tif'
            if digest(path)!=frame['output_sha256']:raise ValueError('Linear RAW hash mismatch: '+str(path))
            image,meta=tiff_channel(path,None)
            if not np.isfinite(image).all():raise ValueError('Nonfinite linear RAW pixels')
            meta.update(radiometry='black_subtracted_linear_camera_'+raw_channel,source_raw_sha256=frame['source_sha256'],
                datetime_original=frame['datetime_original'])
        else:image,meta=tiff_channel(path,channel)
        images[n]=image;inputs.append(meta)
    if any(im.shape!=images[1].shape for im in images.values()):raise ValueError('Exposure dimensions differ')
    if raw_set is not None:
        old_shape=oldinputs['1.tif'].get('shape')
        if old_shape is None:old_shape=tiff_channel(Path(oldinputs['1.tif']['path']),channel)[0].shape
        scale=np.array([images[1].shape[1]/old_shape[1],images[1].shape[0]/old_shape[0]])
        if abs(scale[0]/scale[1]-1)>.01:raise ValueError('RAW canvas aspect differs substantially; independent geometry required')
        center=center*scale;radius=radius*float(np.sqrt(np.prod(scale)))
        geometry_transfer=dict(status='provisional_canvas_scaled_initializer_fresh_landmarks_required',
            old_shape=list(old_shape[:2]),new_shape=list(images[1].shape),scale_xy=scale.tolist(),
            independent_angular_calibration=False)
    destination=out/site;destination.mkdir();reg=destination/'registration';reg.mkdir()
    reference=None if saved_site is not None else landmark_features(images[1],center,radius)
    shared_fits={};shared_geometry=None;site_pairs={}
    if registration_model=='fisheye-shared':
        from spider_shared_fisheye_registration import fit_site
        for n in numbers:
            if n==1:continue
            if saved_site is not None:
                saved=next(t for t in saved_site['transforms'] if t['filename']==f'{n}.tif')
                pairs=saved['landmark_correspondences']
                site_pairs[n]=(np.asarray(pairs['source_xy_px']),np.asarray(pairs['reference_xy_px']))
            else:
                site_pairs[n]=match_landmarks(reference,landmark_features(images[n],center,radius))
        shared_fits,shared_geometry=fit_site(site_pairs,center,radius)
    aligned={1:images[1]};transforms=[]
    fig,axes=plt.subplots(int(np.ceil(len(numbers)/4)),4,figsize=(16,4*int(np.ceil(len(numbers)/4))),squeeze=False,layout='constrained')
    reference_edge=edge_display(images[1]);all_pass=True
    for n,ax in zip(numbers,axes.flat):
        try:
            if n==1:
                matrix=np.array([[1.,0.,0.],[0.,1.,0.]])
                quality=dict(numerical_pass=True,reference_identity=True)
            else:
                if saved_site is not None:
                    saved=next(t for t in saved_site['transforms'] if t['filename']==f'{n}.tif')
                    pairs=saved['landmark_correspondences'];sp=np.asarray(pairs['source_xy_px']);tp=np.asarray(pairs['reference_xy_px'])
                elif registration_model=='fisheye':
                    sp,tp=match_landmarks(reference,landmark_features(images[n],center,radius))
                if registration_model=='fisheye-shared':
                    sp,tp=site_pairs[n];matrix,quality,error=shared_fits[n]
                elif registration_model=='fisheye':
                    if saved_site is not None and saved['quality'].get('numerical_pass') and saved.get('source_to_reference_xy') is not None:
                        matrix=np.asarray(saved['source_to_reference_xy']);quality=dict(saved['quality'],preserved_passing_rigid_transform=True)
                        error=np.linalg.norm(transform_points(sp,matrix)-tp,axis=1)
                    else:
                        from spider_fisheye_registration import fit_correspondences
                        matrix,quality,error=fit_correspondences(sp,tp,center,radius)
                else:
                    matrix,quality,sp,tp,error=register_landmarks(reference,landmark_features(images[n],center,radius),center,radius)
            if not quality['numerical_pass']:all_pass=False
            if n!=1:aligned[n]=warp(images[n],matrix)
            if isinstance(matrix,dict):
                rot=np.asarray(matrix['ray_rotation_source_to_reference']);angle=float(np.rad2deg(np.arctan2(rot[1,0],rot[0,0])))
                source_center=np.array(matrix['source_center_xy_px'])
            else:
                angle=float(np.rad2deg(np.arctan2(matrix[1,0],matrix[0,0])))
                source_center=transform_points([center],np.column_stack([matrix[:,:2].T,-matrix[:,:2].T@matrix[:,2]]))[0]
            entry=dict(filename=f'{n}.tif',source_to_reference_xy=None if isinstance(matrix,dict) else matrix.tolist(),
                image_rotation_clockwise_deg=angle,source_reference_center_xy_px=source_center.tolist(),quality=quality)
            if isinstance(matrix,dict):entry['fisheye_mapping']=matrix
            if n!=1:
                if registration_model=='fisheye-shared':
                    from spider_shared_fisheye_registration import heldout_cells
                    held=heldout_cells(tp)
                else:held=np.arange(len(sp))%4==0
                entry['landmark_correspondences']=dict(source_xy_px=sp.tolist(),reference_xy_px=tp.tolist(),
                    residual_px=error.tolist(),heldout=held.tolist())
            transforms.append(entry)
            edge=edge_display(aligned[n]);ax.imshow(np.stack([reference_edge,edge,edge],axis=-1))
            ax.set_title(f'{n}.tif: '+('PASS numerical' if quality['numerical_pass'] else 'PENDING')+f'; rotation {angle:.2f}°',fontsize=9)
        except (ValueError,KeyError,StopIteration) as exc:
            all_pass=False;transforms.append(dict(filename=f'{n}.tif',quality=dict(numerical_pass=False,error=str(exc))))
            ax.text(.5,.5,str(exc),ha='center',wrap=True,transform=ax.transAxes);ax.set_title(f'{n}.tif: PENDING')
        ax.add_patch(Circle(center,radius,fill=False,color='yellow',lw=.5));ax.axis('off')
    for ax in list(axes.flat)[len(numbers):]:ax.axis('off')
    fig.suptitle(site+': red reference / cyan registered exposure; white means overlapping edges\n'+registration_model+' geometry; no reflection; verify buildings, skyline and stars before accepting',fontsize=13)
    fig.savefig(destination/'registration_edges.png',dpi=150);plt.close(fig)
    lock=dict(site=site,site_id=old['site_id'],status='numerically_locked_awaiting_visual_review' if all_pass else 'pending_registration_not_for_comparison',
        reference_exposure='1.tif',reference_center_xy_px=center.tolist(),rim_radius_px=radius,
        eye_radius_px=radius*np.sin(np.deg2rad(15))/np.sin(np.pi/4),eye_full_field_deg=60,
        geometry_status='Saved reference-image geometry frozen; inherited calibration not independently validated here',
        raw_inputs=inputs,transforms=transforms,source_report_sha256=digest(oldpath),
        astronomical_anchor=anchor,moon_pixel_xy=None,north_orientation_status='unresolved',north_up_dop_path=None,
        coordinate_convention='zero-based original pixels, x right, y down; matrices map source xy to reference xy',
        polarimetry_note='Nominal analyzer grouping and camera-basis Stokes equations unchanged; exposure registration is not a full fisheye Mueller calibration',
        registration_method=registration_model+' from static outer-field SIFT correspondences and withheld landmark checks; no sky polarization objective')
    if shared_geometry is not None:lock['shared_image_geometry']=shared_geometry
    if raw_set is not None:
        lock.update(linear_raw_manifest_sha256=digest(linear_raw_manifest),analyzer_groups=groups,
            excluded_exposures=raw_set['excluded_exposures'],radiometry_status='black_subtracted_linear_camera_'+raw_channel,camera_channel=raw_channel,
            geometry_transfer=geometry_transfer,
            geometry_status='Provisional canvas-scaled reference geometry; newly fitted RAW landmarks; not independent angular calibration')
    if all_pass:
        i,q,u,qa=combine_replicates(([aligned[n] for n in g] for g in groups),'linear_stokes')
        np.savez_compressed(reg/'linear_stokes.npz',I=i,Q=q,U=u)
        report=dict(old,status=lock['status'],reference_center_xy_px=center.tolist(),inputs=inputs,
                    locked_transforms=transforms,registration_lock_schema=1,
                    limitations=old.get('limitations',[])+['Rebuilt from one-pass '+registration_model+' landmark resampling; visual review pending'])
        (reg/'rotation_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
        with np.errstate(divide='ignore',invalid='ignore'):dop=np.hypot(q,u)/i
        y,x=np.indices(i.shape);valid=(np.hypot(x-center[0],y-center[1])<.985*radius)&(i>0)&(dop<=1)&np.isfinite(dop)
        dop=np.where(valid,dop,np.nan)
        moon,meta=moon_coordinate(images[1],center,radius,site,anchor,overrides)
        lock['moon_localization']=meta;lock['moon_pixel_xy']=None if moon is None else moon.tolist()
        north_dop=None
        if moon is not None and anchor and meta.get('elevation_consistent',False):
            matrix,angle=north_rotation(center,moon,anchor['moon_compass_azimuth_deg'])
            north_dop=warp(dop,matrix);north_intensity=warp(i,matrix)
            np.savez_compressed(reg/'north_up_dop.npz',DoLP=north_dop,I=north_intensity)
            lock.update(north_orientation_status='pure_rotation_from_observed_moon_and_fixed_ephemeris_provisional',
                camera_to_north_xy=matrix.tolist(),north_rotation_clockwise_deg=angle,
                north_up_dop_path=f'{site}/registration/north_up_dop.npz',
                handedness='orientation preserving; no mirror fitted; independent handedness verification pending',
                north_frame_note='Scalar DoLP/intensity only; camera-basis Q/U remain in linear_stokes.npz, not relabelled celestial Stokes')
        closure=[dict(median_abs=float(np.nanmedian(abs(r['closure'][valid]))),
                      p95_abs=float(np.nanpercentile(abs(r['closure'][valid]),95))) for r in qa]
        lock['analyzer_pair_closure']=closure
        # Capture-time audit uses original 1.tif, not a ghosted composite.
        fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        lo,hi=np.percentile(images[1][np.isfinite(images[1])],[5,99.8])
        axes[0].imshow(images[1],cmap='gray',vmin=lo,vmax=hi)
        if moon is not None:axes[0].plot(*moon,'rx')
        axes[0].set_title('Raw reference 1.tif; image Moon candidate')
        axes[1].imshow(dop,cmap='turbo',vmin=0,vmax=1);axes[1].set_title('Rebuilt camera-frame DoLP')
        if north_dop is not None:
            axes[2].imshow(north_dop,cmap='turbo',vmin=0,vmax=1);axes[2].set_title('Experimental DoLP rotated to north\nProvisional Moon identity; no mirroring')
            axes[2].text(.03,.97,'N ↑',transform=axes[2].transAxes,va='top')
        else:axes[2].text(.5,.5,'North unresolved',ha='center',transform=axes[2].transAxes)
        for ax in axes:ax.axis('off')
        fig.suptitle(site+' — fixed astronomical Moon coordinates; no DoLP/AoP fit')
        fig.savefig(destination/'registered_dop_and_north.png',dpi=150);plt.close(fig)
    (destination/'reference_coords.json').write_text(json.dumps(lock,indent=2,allow_nan=False)+'\n')
    return lock


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--model-run',type=Path,required=True,help='all_sites directory with fixed astronomical anchors')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--sites',nargs='+',help='Default all sites present in saved batch')
    p.add_argument('--moon-coordinates',type=Path,help='Optional investigator coordinates: {site:{moon_pixel_xy:[x,y],note:...}}')
    p.add_argument('--registration-model',choices=['rigid','fisheye','fisheye-shared'],default='rigid')
    p.add_argument('--linear-raw-manifest',type=Path,help='Prepared linear camera-blue RAW manifest for one explicitly selected site')
    p.add_argument('--from-lock',type=Path,help='Reuse saved image landmark coordinates, not previous transforms')
    a=p.parse_args()
    if a.from_lock and a.registration_model not in ('fisheye','fisheye-shared'):p.error('--from-lock requires a fisheye registration model')
    if a.linear_raw_manifest and (a.from_lock or not a.sites or len(a.sites)!=1):p.error('--linear-raw-manifest requires one --sites entry and fresh landmarks (no --from-lock)')
    if a.output.exists():p.error('Use a fresh versioned output directory')
    try:import cv2
    except ImportError:p.error('Install requirements-registration.txt in the active venv')
    prov=json.loads((a.model_run/'hps/provenance.json').read_text())
    anchors={r['site']:r['anchor'] for r in prov['sites']}
    sites=a.sites or sorted(f.name for f in a.batch.iterdir() if (f/'registration/rotation_analysis.json').exists())
    if not sites:p.error('No saved batch sites found')
    overrides=json.loads(a.moon_coordinates.read_text()) if a.moon_coordinates else {}
    saved=json.loads(a.from_lock.read_text())['sites'] if a.from_lock else {}
    a.output.mkdir(parents=True);locks={}
    for index,site in enumerate(sites):
        print(f'[{index+1}/{len(sites)}] Registering {site} to raw 1.tif',flush=True)
        try:
            if a.from_lock and site not in saved:raise ValueError('Site absent from source landmark lock')
            locks[site]=process_site(site,a.batch,a.output,anchors.get(site),overrides,saved.get(site),a.registration_model,a.linear_raw_manifest)
        except (OSError,ValueError,KeyError) as exc:
            locks[site]=dict(site=site,status='pending_inputs_or_registration',reason=str(exc),moon_pixel_xy=None)
            print('  PENDING: '+str(exc),flush=True)
        manifest=dict(schema_version=1,status='FROZEN_NUMERICAL_CANDIDATES_VISUAL_REVIEW_PENDING',
            sites=locks,raw_images_modified=False,model_parameters_changed=False,
            script_sha256=digest(__file__),opencv_version=cv2.__version__,
            registration_module_sha256={name:digest(Path(__file__).with_name(name)) for name in
                ['spider_fisheye_registration.py','spider_shared_fisheye_registration.py'] if Path(__file__).with_name(name).exists()})
        (a.output/'locked_site_references.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
        summary={name:dict(status=r['status'],north_orientation_status=r.get('north_orientation_status'),
            reason=r.get('reason'),exposures=[dict(filename=t['filename'],quality=t['quality']) for t in r.get('transforms',[])]) for name,r in locks.items()}
        (a.output/'registration_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    checks={str(f.relative_to(a.output)):digest(f) for f in a.output.rglob('*') if f.is_file()}
    (a.output/'reference_sha256.json').write_text(json.dumps(checks,indent=2)+'\n')
    print('Frozen references written. Inspect registration_edges.png before using experimental comparisons.',flush=True)


if __name__=='__main__':main()


