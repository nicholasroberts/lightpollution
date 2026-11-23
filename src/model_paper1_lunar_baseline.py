"""Figure 1B lunar-only geometric baseline, without ALAN or fitted dilution."""
import argparse
import json
from pathlib import Path
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wang2016_sky import natural_stokes
from spider_fisheye_registration import rays
from rotate_spider_measured_north import digest


def lunar_component(view, source, separation_deg):
    """Reuse the earlier Berry-normalized Wang angular/horizon kernel."""
    theta = np.arccos(np.clip(view[...,2], -1, 1))
    phi = np.arctan2(view[...,1], view[...,0])
    source_zenith = np.arccos(np.clip(source[2], -1, 1))
    source_azimuth = np.arctan2(source[1], source[0])
    i,q,u = natural_stokes(theta,phi,source_zenith,source_azimuth,
        separation_deg,np.ones_like(theta),mode='eq6_horizon',
        depolarization_multiplier=1.)
    return i,q+1j*u


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock', type=Path, required=True)
    p.add_argument('--north', type=Path, required=True)
    p.add_argument('--site', default='museum')
    p.add_argument('--neutral-separation-deg', type=float, required=True, help='Use the explicitly recorded prior Berry/Wang neutral-point separation')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists() or not 0<a.neutral_separation_deg<90:
        p.error('Fresh output and 0 < neutral separation < 90 degrees required')
    provenance = json.loads((a.north/'north_rotation.json').read_text())
    if provenance['site'] != a.site or provenance['source_lock_sha256'] != digest(a.reference_lock):
        p.error('North product and reference lock differ')
    lock = json.loads(a.reference_lock.read_text())['sites'][a.site]
    with np.load(a.north/'north_stokes.npz', allow_pickle=False) as z:
        measured = np.array(z['DoLP'])
        measured_angle = np.array(z['AoP_deg'])
        valid = np.array(z['valid'])
        center = np.array(z['center_xy'])
        observed_moon = np.array(z['moon_xy'])
    y,x = np.indices(measured.shape, dtype=float)
    # Same equisolid projection, expressed in Cartesian sky basis x west,y north.
    points = np.stack([x, 2*center[1]-y], -1)
    view = rays(points, center, lock['rim_radius_px'])
    el,az = np.deg2rad([lock['astronomical_anchor']['moon_elevation_deg'],
                       lock['astronomical_anchor']['moon_compass_azimuth_deg']])
    moon = np.array([-np.cos(el)*np.sin(az),np.cos(el)*np.cos(az),np.sin(el)])
    intensity,pol = lunar_component(view,moon,a.neutral_separation_deg)
    predicted = abs(pol)/intensity
    angle = np.rad2deg(.5*np.angle(pol))
    disc = np.hypot(x-center[0],y-center[1])<.985*lock['rim_radius_px']
    predicted = np.where(disc,predicted,np.nan)
    angle = np.where(disc&(abs(pol)>1e-10),angle,np.nan)
    # Predicted Moon position uses ephemeris elevation, not observed radial distance.
    mr = lock['rim_radius_px']*np.sqrt(1-moon[2])
    norm = np.linalg.norm(moon[:2])
    predicted_moon = center+np.array([moon[0],-moon[1]])*mr/max(norm,1e-12)
    field = np.hypot(x-center[0],y-center[1])<lock['eye_radius_px']
    common = field&valid&np.isfinite(predicted)
    if common.sum()<2:
        p.error('Insufficient common central-field pixels')
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'lunar_only_stokes.npz',I=intensity,Q=pol.real,U=pol.imag,
        DoLP=predicted,AoP_deg=angle,common_central_field=common)
    for name,v in [('DoLP_lunar',predicted),('AoP_deg_lunar',angle)]:
        tifffile.imwrite(a.output/(name+'.tif'),v.astype(np.float32),metadata=None)
    fig,axes=plt.subplots(2,2,figsize=(10,9),layout='constrained')
    for row,(d,psi,label) in enumerate([(measured,measured_angle,'Measured green'),
                                       (predicted,angle,'Berry/Wang lunar only')]):
        for col,(v,cmap,lo,hi,title) in enumerate([(d,'viridis',0,1,'DoLP'),(psi,'twilight',-90,90,'AoP')]):
            ax=axes[row,col]
            im=ax.imshow(v,cmap=cmap,vmin=lo,vmax=hi)
            ax.add_patch(plt.Circle(center,lock['eye_radius_px'],fill=False,color='white',lw=.8))
            ax.plot(*observed_moon,'kx',ms=6)
            ax.plot(*predicted_moon,'r+',ms=6)
            ax.set_title(label+' '+title)
            ax.text(.5,.98,'N',transform=ax.transAxes,ha='center')
            ax.text(.02,.5,'E',transform=ax.transAxes)
            r=lock['rim_radius_px']
            ax.set_xlim(center[0]-1.06*r,center[0]+1.06*r)
            ax.set_ylim(center[1]+1.06*r,center[1]-1.06*r)
            ax.axis('off');fig.colorbar(im,ax=ax,shrink=.8)
    fig.suptitle(a.site+': lunar-only baseline; no dilution or sky fit\nMeasured AoP has unknown constant zero; absolute colors are not calibrated agreement')
    for ext in ['png','pdf','svg']:
        fig.savefig(a.output/('fig01b_lunar_only.'+ext),dpi=220)
    plt.close(fig)
    stats = {name:dict(mean=float(v[common].mean()),spatial_sd=float(v[common].std(ddof=1)))
             for name,v in [('measured_north_resampled',measured),('lunar',predicted)]}
    report=dict(site=a.site,source_lock_sha256=digest(a.reference_lock),
        source_north_arrays_sha256=digest(a.north/'north_stokes.npz'),
        source_north_report_sha256=digest(a.north/'north_rotation.json'),
        script_sha256=digest(__file__),kernel_sha256=digest(Path(__file__).with_name('wang2016_sky.py')),
        projection_module_sha256=digest(Path(__file__).with_name('spider_fisheye_registration.py')),
        neutral_separation_deg=a.neutral_separation_deg,solar_ratio=0,unpolarized_background=0,horizon_depolarization='Wang eq8: cos(theta)^(1/10)',
        atmosphere='Berry-normalized geometric envelope times Wang eq8; atmospheric multiplier unity. No full radiative transfer or Wang eq17 spectral term.',
        spectral_assumption='Achromatic polarization and separable common lunar spectral template; common response integral cancels in DoLP/AoP. Camera-green comparator, not claimed animal-weighted prediction.',
        projection='Zenith-axis equisolid; rim radius at 90 degree zenith angle',
        astronomical_anchor=lock['astronomical_anchor'],common_central_pixel_count=int(common.sum()),
        central_statistics=stats,observed_moon_xy=observed_moon.tolist(),predicted_moon_xy=predicted_moon.tolist(),
        limitations=['Full-disc model is geometric context; no automatic foreground segmentation.',
        'Central support is physical measured pixels, not an independently reviewed clear-sky mask.',
        'Absolute measured analyzer zero unknown; compare spatial AoP structure separately.',
        'Earlier Berry source-meridian phase convention retained; full camera-basis correspondence is not independently calibrated.',
        'Resampled statistics are diagnostics; published camera-frame mean remains separately recorded.',
        'No imposed pmax or clipping. Berry geometric intensity used as DoLP envelope is the earlier declared assumption. Uniform natural I=1. No fit performed.'])
    (a.output/'lunar_only.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(stats,indent=2))
    print('Saved',a.output)


if __name__=='__main__':
    main()
