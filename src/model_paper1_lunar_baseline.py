"""Figure 1B lunar-only geometric baseline, without ALAN or fitted dilution."""
import argparse
import json
from pathlib import Path
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_spider_lunar_aop import lunar_phase
from spider_fisheye_registration import rays
from rotate_spider_measured_north import digest


def lunar_component(view, source, pmax):
    mu = np.clip(np.sum(view*source, axis=-1), -1, 1)
    intensity = (1+mu**2)/(1+source[2]**2)
    amplitude = pmax*(1-mu**2)/(1+source[2]**2)
    phase = lunar_phase(view, source)
    pol = np.where(amplitude < 1e-12, 0j, amplitude*phase)
    return intensity, pol


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock', type=Path, required=True)
    p.add_argument('--north', type=Path, required=True)
    p.add_argument('--site', default='museum')
    p.add_argument('--pmax', type=float, required=True, help='0.8 pilot baseline; 1 ideal Rayleigh sensitivity')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists() or not 0<a.pmax<=1:
        p.error('Fresh output and 0 < pmax <= 1 required')
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
    intensity,pol = lunar_component(view,moon,a.pmax)
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
                                       (predicted,angle,'Lunar only; pmax='+str(a.pmax))]):
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
        script_sha256=digest(__file__),phase_module_sha256=digest(Path(__file__).with_name('audit_spider_lunar_aop.py')),
        projection_module_sha256=digest(Path(__file__).with_name('spider_fisheye_registration.py')),
        pmax=a.pmax,solar_ratio=0,unpolarized_background=0,horizon_depolarization='none',
        atmosphere='No extinction or multiple-scattering radiative transfer; fixed pmax is a declared phenomenological amplitude.',
        spectral_assumption='Achromatic polarization and separable common lunar spectral template; common response integral cancels in DoLP/AoP. Camera-green comparator, not claimed animal-weighted prediction.',
        projection='Zenith-axis equisolid; rim radius at 90 degree zenith angle',
        astronomical_anchor=lock['astronomical_anchor'],common_central_pixel_count=int(common.sum()),
        central_statistics=stats,observed_moon_xy=observed_moon.tolist(),predicted_moon_xy=predicted_moon.tolist(),
        limitations=['Full-disc model is geometric context; no automatic foreground segmentation.',
        'Central support is physical measured pixels, not an independently reviewed clear-sky mask.',
        'Absolute measured analyzer zero unknown; compare spatial AoP structure separately.',
        'Ideal meridional polarization transport, not measured lens Mueller calibration.',
        'Resampled statistics are diagnostics; published camera-frame mean remains separately recorded.',
        'pmax=0.8 is not an independently calibrated atmospheric estimate. No fit performed.'])
    (a.output/'lunar_only.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(stats,indent=2))
    print('Saved',a.output)


if __name__=='__main__':
    main()
