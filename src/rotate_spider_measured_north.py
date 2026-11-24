"""Rotate audited measured Stokes into a Moon-anchored north-up sky view."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def coordinate_matrices(center, moon, azimuth, handedness):
    v = np.asarray(moon, float)-center
    if np.linalg.norm(v) < 1:
        raise ValueError('Moon too near centre to anchor orientation')
    alpha = np.deg2rad(azimuth-90)-np.arctan2(v[1], v[0])
    co, si = np.cos(alpha), np.sin(alpha)
    rotation = np.array([[co, -si], [si, co]])
    h = np.eye(2)
    if handedness == 'reversed':
        direction = v/np.linalg.norm(v)
        h = 2*np.outer(direction, direction)-np.eye(2)
    # Cartesian sky coordinates: +x west/right, +y north/up.
    cartesian = -rotation@h
    # Raster coordinates: right and down.
    raster = np.diag([1., -1.])@cartesian
    return cartesian, raster


def tensor_transform(q, u, a):
    return ((a[0,0]**2-a[0,1]**2)*q+2*a[0,0]*a[0,1]*u,
            (a[0,0]*a[1,0]-a[0,1]*a[1,1])*q+
            (a[0,0]*a[1,1]+a[0,1]*a[1,0])*u)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock', type=Path, required=True)
    p.add_argument('--audit', type=Path, required=True)
    p.add_argument('--site', default='museum')
    p.add_argument('--handedness', choices=['original','reversed'], required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        p.error('Use a fresh output directory')
    report = json.loads((args.audit/'stokes_audit.json').read_text())
    if report['site'] != args.site or report['source_lock_sha256'] != digest(args.reference_lock):
        p.error('Audit site or source lock differs')
    lock = json.loads(args.reference_lock.read_text())['sites'][args.site]
    center = np.asarray(lock['reference_center_xy_px'], float)
    if lock.get('moon_pixel_xy') is None or not lock.get('astronomical_anchor'):
        p.error('Observed Moon and fixed astronomical anchor required')
    az = lock['astronomical_anchor']['moon_compass_azimuth_deg']
    a, b = coordinate_matrices(center, lock['moon_pixel_xy'], az, args.handedness)
    source = args.audit/'measured_stokes_and_polarization.npz'
    with np.load(source, allow_pickle=False) as z:
        intensity, q, u = [np.array(z[k]) for k in ['I_pair','Q','U']]
    q, u = tensor_transform(q, u, a)
    y, x = np.indices(intensity.shape, dtype=float)
    # Inverse of orthogonal source->destination raster transform.
    dx, dy = x-center[0], y-center[1]
    sx = b[0,0]*dx+b[1,0]*dy+center[0]
    sy = b[0,1]*dx+b[1,1]*dy+center[1]
    intensity, q, u = [map_coordinates(v, [sy,sx], order=1,
        mode='constant', cval=np.nan, prefilter=False) for v in [intensity,q,u]]
    with np.errstate(divide='ignore', invalid='ignore'):
        dolp = np.hypot(q,u)/intensity
    aop = np.rad2deg(.5*np.arctan2(u,q))
    radius = lock['rim_radius_px']
    valid = (np.hypot(dx,dy)<.985*radius)&(intensity>0)&np.isfinite(dolp)&(dolp<=1)
    angle_valid = valid&(np.hypot(q,u)>0)
    dolp = np.where(valid,dolp,np.nan)
    aop = np.where(angle_valid,aop,np.nan)
    predicted_moon=lock.get('moon_position_status','').startswith('predicted_')
    moon = center+b@(np.asarray(lock['moon_pixel_xy'])-center)
    args.output.mkdir(parents=True)
    for name, values in [('DoLP_north',dolp),('AoP_deg_north',aop)]:
        tifffile.imwrite(args.output/(name+'.tif'),values.astype(np.float32),metadata=None)
    np.savez_compressed(args.output/'north_stokes.npz',I=intensity,Q=q,U=u,
        DoLP=dolp,AoP_deg=aop,valid=valid,camera_to_sky_cartesian=a,
        camera_to_raster=b,center_xy=center,moon_xy=moon)
    fig, axes = plt.subplots(1,2,figsize=(10,5),layout='constrained')
    for ax,v,cmap,lo,hi,title in zip(axes,[dolp,aop],['viridis','twilight'],
            [0,-90],[1,90],['Measured green DoLP','AoP; unknown analyzer zero']):
        im=ax.imshow(v,cmap=cmap,vmin=lo,vmax=hi)
        ax.plot(*moon,'r+' if predicted_moon else 'kx',ms=7)
        for tx,ty,t in [(.5,.98,'N'),(.02,.5,'E'),(.5,.02,'S'),(.98,.5,'W')]:
            ax.text(tx,ty,t,transform=ax.transAxes,ha='center',va='center')
        ax.set_xlim(center[0]-1.08*radius,center[0]+1.08*radius)
        ax.set_ylim(center[1]+1.08*radius,center[1]-1.08*radius)
        ax.set_title(title);ax.axis('off');fig.colorbar(im,ax=ax,shrink=.8)
    fig.suptitle(args.site+(': star-anchored; predicted Moon' if predicted_moon else ': observed-Moon anchored')+'; north up, east left')
    fig.savefig(args.output/'DoLP_AoP_north.png',dpi=200)
    plt.close(fig)
    provenance=dict(site=args.site,source_lock_sha256=digest(args.reference_lock),
        source_audit_sha256=digest(args.audit/'stokes_audit.json'),
        source_arrays_sha256=digest(source),script_sha256=digest(__file__),
        moon_position_status=lock.get('moon_position_status','unspecified'),celestial_orientation=lock.get('celestial_orientation'),
        handedness_hypothesis=args.handedness,moon_compass_azimuth_deg=az,
        camera_to_sky_cartesian=a.tolist(),camera_to_raster=b.tolist(),
        note='Pair estimator retained. Transform I,Q,U then calculate DoLP/AoP. One bilinear spatial resampling; no interpolation of angles. No sky-pattern fit. AoP basis +x west/right, +y north/up; unknown constant analyzer zero and off-axis calibration remain. Reversed parity is explicit and separate from analyzer sequence sign. Moon radius not adjusted to match astronomical elevation.')
    (args.output/'north_rotation.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Saved',args.output)


if __name__=='__main__':
    main()
