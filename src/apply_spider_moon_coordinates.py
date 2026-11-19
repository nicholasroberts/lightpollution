"""Apply investigator Moon pixels to saved RAW registrations without refitting.

Writes a fresh derived lock, unchanged camera Stokes, and Moon geometry panels.
The inherited zenith/projection and camera handedness remain provisional.
"""
import argparse
import copy
import hashlib
import json
import math
import shutil
from pathlib import Path
import numpy as np
from scipy import ndimage


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require_hash(path, expected):
    if not expected or digest(path)!=expected:
        raise ValueError('Source hash mismatch: '+str(path))


def geometry(center,radius,moon,anchor):
    center=np.asarray(center,dtype=float);moon=np.asarray(moon,dtype=float)
    delta=moon-center;distance=float(np.linalg.norm(delta))
    if radius<=0 or not np.isfinite(radius) or distance<1:
        raise ValueError('Invalid radius or Moon too near zenith')
    implied=90-np.rad2deg(2*np.arcsin(np.clip(distance/(radius*np.sqrt(2)),0,1)))
    residual=float(implied-anchor['moon_elevation_deg'])
    alpha=math.radians(anchor['moon_compass_azimuth_deg']-90)-math.atan2(delta[1],delta[0])
    c,s=math.cos(alpha),math.sin(alpha);rotation=np.array([[c,-s],[s,c]])
    matrix=np.column_stack([rotation,center-rotation@center])
    expected_radius=radius*np.sqrt(2)*np.sin(np.deg2rad((90-anchor['moon_elevation_deg'])/2))
    predicted=center+delta/distance*expected_radius
    return matrix,dict(image_implied_elevation_deg=float(implied),elevation_residual_deg=residual,
        elevation_consistent=abs(residual)<=6,moon_pixel_separation=float(abs(expected_radius-distance)),
        predicted_moon_camera_xy=predicted.tolist(),north_rotation_clockwise_deg=float(np.rad2deg(alpha)))


def warp(image,matrix):
    y,x=np.indices(image.shape,dtype=float);r=matrix[:,:2];t=matrix[:,2]
    sx=r[0,0]*(x-t[0])+r[1,0]*(y-t[1]);sy=r[0,1]*(x-t[0])+r[1,1]*(y-t[1])
    return ndimage.map_coordinates(image,[sy,sx],order=1,mode='constant',cval=np.nan,prefilter=False)


def apply(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import tifffile
    summary_path=a.registration/'registration_batch_summary.json'
    rows={r['site']:r for r in json.loads(summary_path.read_text())['sites']}
    review=json.loads(a.visual_review.read_text())
    require_hash(summary_path,review['registration_summary_sha256'])
    reviews={r['site']:r for r in review['sites']}
    overrides=json.loads(a.moon_coordinates.read_text());prepared=[]
    for site in a.sites:
        if Path(site).name!=site or site in ('.','..'):raise ValueError('Invalid site name')
        r=rows[site];o=overrides[site]
        if reviews[site]['alignment']!='acceptable':raise ValueError('Visual alignment not accepted: '+site)
        if r['status']!='numerically_locked_awaiting_visual_review':raise ValueError('Numerical registration pending: '+site)
        source=a.registration/r['lock_path'];require_hash(source,r['lock_sha256'])
        require_hash(source,o['registration_lock_sha256'])
        lock=copy.deepcopy(json.loads(source.read_text())['sites'][site])
        checks=json.loads((source.parent/'reference_sha256.json').read_text())
        # Verify all saved products before copying them unchanged.
        for rel,sha in checks.items():
            f=source.parent/rel
            if not f.resolve().is_relative_to(source.parent.resolve()):raise ValueError('Unsafe checksum path')
            require_hash(f,sha)
        raw=a.raw_root/site/'1.tif';manifest=raw.parent/'linear_raw_manifest.json'
        require_hash(raw,o['reference_tiff_sha256']);require_hash(manifest,o['linear_raw_manifest_sha256'])
        if lock['linear_raw_manifest_sha256']!=o['linear_raw_manifest_sha256']:raise ValueError('Lock/RAW manifest mismatch')
        image=tifffile.imread(raw)
        if list(image.shape)!=o['image_shape_yx'] or image.ndim!=2:raise ValueError('Reference image dimensions changed')
        moon=np.asarray(o['moon_pixel_xy'],dtype=float)
        if moon.shape!=(2,) or not np.isfinite(moon).all() or not (0<=moon[0]<image.shape[1] and 0<=moon[1]<image.shape[0]):raise ValueError('Moon pixel out of bounds')
        if not o.get('note'):raise ValueError('Moon identity note missing')
        stokes=source.parent/site/'registration/linear_stokes.npz'
        with np.load(stokes) as z:
            i,q,u=(np.array(z[k],dtype=float) for k in ['I','Q','U'])
        if any(v.shape!=image.shape for v in (i,q,u)):raise ValueError('Stokes/reference shape mismatch')
        matrix,meta=geometry(lock['reference_center_xy_px'],lock['rim_radius_px'],moon,lock['astronomical_anchor'])
        prepared.append((site,source,lock,image,moon,matrix,meta,i,q,u,stokes))
    a.output.mkdir(parents=True,exist_ok=False);locks={};results=[]
    for site,source,lock,image,moon,matrix,meta,i,q,u,stokes in prepared:
        shutil.copytree(source.parent/site,a.output/site)
        destination=a.output/site;reg=destination/'registration'
        for name in ['north_up_dop.npz']:
            (reg/name).unlink(missing_ok=True)
        o=overrides[site]
        lock.update(moon_pixel_xy=moon.tolist(),moon_localization=dict(meta,status='investigator_supplied_image_coordinate',
            pixel_xy=moon.tolist(),note=o['note']),visual_alignment_status='investigator_accepted',
            source_registration_lock_sha256=digest(source),moon_overrides_sha256=digest(a.moon_coordinates),
            visual_review_sha256=digest(a.visual_review),north_up_dop_path=None,north_orientation_status='pending_projection_elevation_check')
        for key in ['camera_to_north_xy','north_rotation_clockwise_deg']:lock.pop(key,None)
        with np.errstate(divide='ignore',invalid='ignore'):dop=np.hypot(q,u)/i
        y,x=np.indices(image.shape);center=np.array(lock['reference_center_xy_px'])
        valid=np.isfinite(dop)&(i>0)&(dop<=1)&(np.hypot(x-center[0],y-center[1])<.985*lock['rim_radius_px'])
        dop=np.where(valid,dop,np.nan)
        north=None
        if meta['elevation_consistent']:
            north=warp(dop,matrix)
            np.savez_compressed(reg/'north_up_dop.npz',DoLP=north,I=warp(i,matrix))
            lock.update(camera_to_north_xy=matrix.tolist(),north_rotation_clockwise_deg=meta['north_rotation_clockwise_deg'],
                north_orientation_status='pure_rotation_from_observed_moon_and_fixed_ephemeris_provisional',
                north_up_dop_path=f'{site}/registration/north_up_dop.npz')
        lock['north_frame_note']='Moon identity supplied by investigator; inherited zenith/projection and handedness provisional. Proper rotation with east right; camera Q/U unchanged. No polarization fit.'
        # Preserve alignment and Stokes exactly, not just within tolerance.
        require_hash(reg/'linear_stokes.npz',digest(stokes))
        if lock['transforms']!=json.loads(source.read_text())['sites'][site]['transforms']:raise ValueError('Unexpected transform mutation')
        fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        positive=image[image>0];lo,hi=np.percentile(positive,[1,99.9])
        display=np.log1p(100*np.clip((image-lo)/max(hi-lo,1e-12),0,1))/np.log(101)
        axes[0].imshow(display,cmap='gray');axes[0].plot(*moon,'rx');axes[0].set_title('Reference 1.tif; investigator Moon')
        axes[1].imshow(dop,cmap='viridis',vmin=0,vmax=1);axes[1].plot(*moon,'kx');axes[1].plot(*meta['predicted_moon_camera_xy'],'r+');axes[1].set_title('Camera DoLP; fixed-elevation check')
        if north is not None:
            axes[2].imshow(north,cmap='viridis',vmin=0,vmax=1)
            for point,style in [(moon,'kx'),(meta['predicted_moon_camera_xy'],'r+')]:
                xy=matrix[:,:2]@point+matrix[:,2];axes[2].plot(*xy,style)
            axes[2].set_title('Moon-anchored north; parity provisional')
            axes[2].text(.02,.96,'N ↑',transform=axes[2].transAxes,va='top');axes[2].text(.9,.5,'E →',transform=axes[2].transAxes)
        else:axes[2].text(.5,.5,'Elevation/projection check pending',ha='center',transform=axes[2].transAxes)
        for ax in axes:ax.axis('off')
        fig.suptitle(site+f" — unchanged alignment and camera Stokes; elevation residual {meta['elevation_residual_deg']:.2f}°\nBlack ×: observed Moon; red +: fixed astronomical elevation under inherited projection. No position or pattern fit.")
        fig.savefig(destination/'registered_dop_and_north.png',dpi=160);plt.close(fig)
        (destination/'reference_coords.json').write_text(json.dumps(lock,indent=2,allow_nan=False)+'\n')
        locks[site]=lock;results.append(dict(site=site,moon_pixel_xy=moon.tolist(),**meta,
            north_orientation_status=lock['north_orientation_status'],registration_refitted=False,camera_stokes_unchanged=True))
    (a.output/'locked_site_references.json').write_text(json.dumps(dict(schema_version=1,status='DERIVED_MOON_COORDINATES_PROJECTION_AND_PARITY_PROVISIONAL',sites=locks,raw_images_modified=False,model_parameters_changed=False),indent=2,allow_nan=False)+'\n')
    (a.output/'moon_coordinate_summary.json').write_text(json.dumps(dict(sites=results,unapplied_override_sites=sorted(set(overrides)-set(a.sites))),indent=2)+'\n')
    checks={str(f.relative_to(a.output)):digest(f) for f in a.output.rglob('*') if f.is_file()}
    (a.output/'reference_sha256.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(results,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registration',type=Path,required=True)
    p.add_argument('--raw-root',type=Path,required=True)
    p.add_argument('--moon-coordinates',type=Path,required=True)
    p.add_argument('--visual-review',type=Path,required=True)
    p.add_argument('--sites',nargs='+',default=['WClawn','oxfordbussinesspark'])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    try:apply(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))


if __name__=='__main__':main()
