"""Offset-free AoP geometry checks of a checksum-locked linear RAW reconstruction.

No sky-pattern rotation, analyzer-zero fit, mirroring or ALAN fit is performed.
The primary prediction uses fixed astronomical Moon coordinates and the saved
camera-to-north rotation, with a provisional zenith-axis/equisolid assumption.
An observed-Moon-direction sensitivity check is reported separately.
"""
import argparse,json,hashlib,csv
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spider_fisheye_registration import rays,pixels


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def block_mean(image,size):
    h,w=image.shape;h=h//size*size;w=w//size*size
    blocks=image[:h,:w].reshape(h//size,size,w//size,size)
    count=np.isfinite(blocks).sum(axis=(1,3));total=np.nansum(blocks,axis=(1,3))
    return np.divide(total,count,out=np.full(total.shape,np.nan),where=count>0)


def image_tangent_basis(vectors):
    """Ideal meridional ray transport to the reference rear-analyzer plane."""
    x,y,z=np.moveaxis(vectors,-1,0);rho=np.hypot(x,y)
    cp=np.divide(x,rho,out=np.ones_like(x),where=rho>1e-12)
    sp=np.divide(y,rho,out=np.zeros_like(y),where=rho>1e-12)
    meridian=np.stack([z*cp,z*sp,-rho],axis=-1)
    azimuthal=np.stack([-sp,cp,np.zeros_like(x)],axis=-1)
    return cp[...,None]*meridian-sp[...,None]*azimuthal,sp[...,None]*meridian+cp[...,None]*azimuthal


def lunar_phase(vectors,moon):
    """exp(2i AoP) for ideal single scattering; no intensity/DoLP fit."""
    field=np.cross(moon,vectors);ex,ey=image_tangent_basis(vectors)
    field_xy=np.sum(field*ex,axis=-1)+1j*np.sum(field*ey,axis=-1)
    norm=np.abs(field_xy)
    unit=np.divide(field_xy,norm,out=np.full(field_xy.shape,np.nan+1j*np.nan),where=norm>1e-9)
    return unit**2


def fixed_moon_in_camera(anchor,camera_to_north):
    el,az=np.deg2rad([anchor['moon_elevation_deg'],anchor['moon_compass_azimuth_deg']])
    # North-up image coordinates: x east, y south/down, z zenith.
    north=np.array([np.cos(el)*np.sin(az),-np.cos(el)*np.cos(az),np.sin(el)])
    rotation=np.asarray(camera_to_north,dtype=float)[:,:2]
    if not np.allclose(rotation.T@rotation,np.eye(2),atol=1e-6) or np.linalg.det(rotation)<.999999:
        raise ValueError('Saved camera-to-north mapping must be a proper rotation')
    return np.r_[rotation.T@north[:2],north[2]]


def solar_geometry(anchor,site_name,reference_path,camera_to_north):
    from astropy.utils import iers
    from astropy.coordinates import EarthLocation,AltAz,get_sun
    from astropy.time import Time
    import astropy.units as units
    with Path(reference_path).open() as stream:rows=list(csv.DictReader(stream))
    matches=[r for r in rows if Path(r['source_filename']).stem==site_name]
    if len(matches)!=1:raise ValueError('Need one verified site reference row for the solar control')
    row=matches[0]
    if anchor.get('site_id') and row['site_id']!=anchor['site_id']:raise ValueError('Solar-control site identity differs from Moon anchor')
    utc=datetime.fromisoformat(anchor['utc'].replace('Z','+00:00'))
    if utc.tzinfo is None:utc=utc.replace(tzinfo=timezone.utc)
    iers.conf.auto_download=False
    location=EarthLocation(lat=float(row['latitude'])*units.deg,lon=float(row['longitude'])*units.deg)
    moment=Time(utc.astimezone(timezone.utc))
    sun=get_sun(moment).transform_to(AltAz(obstime=moment,location=location,pressure=0*units.hPa))
    solar=dict(moon_elevation_deg=float(sun.alt.deg),moon_compass_azimuth_deg=float(sun.az.deg))
    vector=fixed_moon_in_camera(solar,camera_to_north)
    return vector,dict(sun_elevation_deg=float(sun.alt.deg),sun_compass_azimuth_deg=float(sun.az.deg),
        utc=utc.isoformat(),latitude_deg=float(row['latitude']),longitude_deg=float(row['longitude']),
        site_reference_sha256=digest(reference_path),pressure_hpa=0,source='astropy.get_sun; offline IERS',
        interpretation='Ideal solar scattering-plane control only; no twilight brightness or Moon/Sun mixing fit')


def phase_changes(phase,axis,lag):
    if lag<1 or lag>=phase.shape[axis]:raise ValueError('Lag must fit inside the diagnostic grid')
    lo=[slice(None)]*2;hi=lo.copy();lo[axis]=slice(None,-lag);hi[axis]=slice(lag,None)
    return phase[tuple(hi)]*np.conj(phase[tuple(lo)])


def change_residual(measured,predicted,mask,axis,lag):
    m=phase_changes(measured,axis,lag);p=phase_changes(predicted,axis,lag)
    lo=[slice(None)]*2;hi=lo.copy();lo[axis]=slice(None,-lag);hi[axis]=slice(lag,None)
    good=mask[tuple(lo)]&mask[tuple(hi)]&np.isfinite(m)&np.isfinite(p)
    result=.5*np.angle(m*np.conj(p))
    return np.where(good,result,np.nan)


def residual_summary(measured,predicted,mask,lag):
    values=np.concatenate([change_residual(measured,predicted,mask,axis,lag).ravel() for axis in [0,1]])
    values=values[np.isfinite(values)]
    if not len(values):return dict(pair_count=0,median_abs_change_difference_deg=None,p95_abs_change_difference_deg=None)
    return dict(pair_count=len(values),median_abs_change_difference_deg=float(np.rad2deg(np.median(abs(values)))),
                p95_abs_change_difference_deg=float(np.rad2deg(np.percentile(abs(values),95))))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--site',default='museum')
    p.add_argument('--reference',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/prc_spider_18_sites_verified.csv',
        help='Verified site coordinates for the independent solar geometry control')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bin-size',type=int,default=8,help='Average I/Q/U blocks before calculating angles')
    p.add_argument('--lag-bins',type=int,default=4,help='Pair spacing on the binned grid')
    p.add_argument('--min-dolp',type=float,default=.05)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    if not a.reference.is_file():p.error('Verified site coordinates missing; supply --reference')
    if a.bin_size<1 or a.lag_bins<1:p.error('Bin size and lag must be positive')
    if not 0<a.min_dolp<.95:p.error('Choose min DoLP between 0 and 0.95')
    base=a.reference_lock.parent;checks=json.loads((base/'reference_sha256.json').read_text())
    for path in [a.reference_lock,base/a.site/'registration/linear_stokes.npz']:
        key=str(path.relative_to(base))
        if key not in checks or digest(path)!=checks[key]:p.error('Checksum mismatch or unlisted input: '+key)
    site=json.loads(a.reference_lock.read_text())['sites'][a.site]
    if site.get('radiometry_status')!='black_subtracted_linear_camera_blue':p.error('Use the new linear RAW registration lock')
    if site['status']!='numerically_locked_awaiting_visual_review':p.error('Registration has not passed its numerical checks')
    if not site.get('astronomical_anchor',{}).get('utc'):p.error('The locked capture UTC is required')
    if not site.get('camera_to_north_xy') or site.get('moon_pixel_xy') is None:p.error('A locked image Moon and north rotation are required')
    c=np.asarray(site['reference_center_xy_px']);radius=site['rim_radius_px'];moon_xy=np.asarray(site['moon_pixel_xy'])
    with np.load(base/a.site/'registration/linear_stokes.npz') as z:
        I,Q,U=[block_mean(z[k],a.bin_size) for k in ['I','Q','U']]
    if a.lag_bins>=min(I.shape):p.error('Lag exceeds image dimensions')
    yy,xx=np.indices(I.shape,dtype=float);xx=xx*a.bin_size+(a.bin_size-1)/2;yy=yy*a.bin_size+(a.bin_size-1)/2
    points=np.stack([xx,yy],axis=-1);vectors=rays(points,c,radius)
    primary_moon=fixed_moon_in_camera(site['astronomical_anchor'],site['camera_to_north_xy'])
    observed_moon=rays(moon_xy[None,:],c,radius)[0]
    predictions={'fixed_ephemeris_zenith_axis':lunar_phase(vectors,primary_moon),
                 'locked_image_direction_sensitivity':lunar_phase(vectors,observed_moon),
                 'constant_aop_control':np.ones(I.shape,dtype=complex)}
    solar_metadata=None
    if a.reference.is_file() and site['astronomical_anchor'].get('utc'):
        sun_vector,solar_metadata=solar_geometry(site['astronomical_anchor'],a.site,a.reference,site['camera_to_north_xy'])
        predictions['fixed_sun_geometry_control']=lunar_phase(vectors,sun_vector)
    with np.errstate(invalid='ignore',divide='ignore'):
        dolp=np.hypot(Q,U)/I;measured=(Q+1j*U)/np.hypot(Q,U)
    rho=np.hypot(xx-c[0],yy-c[1]);eye_radius=site.get('eye_radius_px',radius*np.sin(np.pi/12)/np.sin(np.pi/4))
    usable=np.isfinite(I)&np.isfinite(measured)&(I>0)&(dolp>=a.min_dolp)&(dolp<=.95)
    # Mask both source singularities; bright Moon/nearby structures are not angle references.
    for moon in [primary_moon,observed_moon]:
        usable&=abs(np.sum(vectors*moon,axis=-1))<np.cos(np.deg2rad(10))
    mask=usable&(rho<eye_radius)
    if mask.sum()<20:p.error('Too few usable central-field angle samples')
    metrics={}
    for sign,phase in [('stored_U',measured),('reversed_U_diagnostic_only',np.conj(measured))]:
        metrics[sign]={}
        for name,pred in predictions.items():
            metrics[sign][name]={str(lag*a.bin_size):residual_summary(phase,pred,mask,lag)
                for lag in sorted(set([2,a.lag_bins,8])) if lag<min(I.shape)}
    primary=predictions['fixed_ephemeris_zenith_axis'];residual=change_residual(measured,primary,mask,1,a.lag_bins)
    rx=xx[:,:-a.lag_bins]+a.bin_size*a.lag_bins/2;ry=yy[:,:-a.lag_bins]
    observed_angle=.5*np.angle(measured);expected_angle=.5*np.angle(primary)
    measured_change=.5*np.angle(phase_changes(measured,1,a.lag_bins))
    expected_change=.5*np.angle(phase_changes(primary,1,a.lag_bins))
    fig,axes=plt.subplots(2,3,figsize=(15,10),layout='constrained')
    display=usable&(rho<.985*radius)
    for ax,data,title,cmap,vmin,vmax in [
        (axes[0,0],dolp,'Measured linear-RAW DoLP','turbo',0,1),
        (axes[0,1],np.rad2deg(observed_angle),'Measured AoP: zero axis not calibrated','twilight',-90,90),
        (axes[0,2],np.rad2deg(expected_angle),'Fixed-Moon ideal single-scattering AoP','twilight',-90,90)]:
        art=ax.pcolormesh(xx,yy,np.where(display,data,np.nan),shading='nearest',cmap=cmap,vmin=vmin,vmax=vmax)
        fig.colorbar(art,ax=ax,shrink=.65,label='DoLP' if ax is axes[0,0] else 'degrees, modulo 180')
    valid=np.isfinite(residual)
    for ax,data,title in [(axes[1,0],np.rad2deg(measured_change),'Measured horizontal AoP change'),
                          (axes[1,1],np.rad2deg(expected_change),'Predicted horizontal AoP change'),
                          (axes[1,2],np.rad2deg(residual),'Change difference; global zero cancels')]:
        art=ax.pcolormesh(rx,ry,np.where(valid,data,np.nan),shading='nearest',cmap='RdBu_r',vmin=-15,vmax=15)
        fig.colorbar(art,ax=ax,shrink=.65,label='degrees; display clipped at ±15')
    predicted_xy=pixels(primary_moon[None,:],c,radius)[0]
    for ax,title in zip(axes.flat,['Measured linear-RAW DoLP','Measured AoP: zero axis not calibrated',
        'Fixed-Moon ideal single-scattering AoP','Measured horizontal AoP change',
        'Predicted horizontal AoP change','Change difference; global zero cancels']):
        ax.set_title(title,fontsize=10);ax.plot(*moon_xy,'kx',ms=5);ax.plot(*predicted_xy,'r+',ms=6)
        ax.add_patch(plt.Circle(c,eye_radius,fill=False,color='white',lw=.8))
        ax.set(xlim=(c[0]-radius*1.06,c[0]+radius*1.06),ylim=(c[1]+radius*1.06,c[1]-radius*1.06),aspect='equal');ax.axis('off')
    fig.suptitle(a.site+' — AoP geometry audit; no rotation/zero-angle fit to the pattern\n'
        f'Camera frame. Black ×: observed Moon; red +: fixed-ephemeris prediction with provisional zenith axis. '
        f'White circle: central field; pair spacing {a.lag_bins*a.bin_size} px.',fontsize=12)
    a.output.mkdir(parents=True);fig.savefig(a.output/'aop_geometry_audit.png',dpi=160);plt.close(fig)
    if solar_metadata is not None:
        fig,axs=plt.subplots(2,3,figsize=(15,10),layout='constrained')
        solar=predictions['fixed_sun_geometry_control']
        for ax,phase,title in zip(axs[0],[measured,primary,solar],
            ['Measured AoP (uncalibrated constant zero)','Fixed-Moon scattering-plane AoP','Fixed-Sun scattering-plane AoP']):
            art=ax.pcolormesh(xx,yy,np.where(display,np.rad2deg(.5*np.angle(phase)),np.nan),
                shading='nearest',cmap='twilight',vmin=-90,vmax=90)
            ax.set_title(title,fontsize=10);fig.colorbar(art,ax=ax,shrink=.65,label='degrees, modulo 180')
        for ax,pred,title in zip(axs[1],[primary,solar,np.ones(I.shape,dtype=complex)],
            ['Measured-minus-Moon AoP change','Measured-minus-Sun AoP change','Measured-minus-constant AoP change']):
            error=change_residual(measured,pred,mask,1,a.lag_bins)
            art=ax.pcolormesh(rx,ry,np.rad2deg(error),shading='nearest',cmap='RdBu_r',vmin=-15,vmax=15)
            ax.set_title(title,fontsize=10);fig.colorbar(art,ax=ax,shrink=.65,label='degrees; display clipped at ±15')
        for ax in axs.flat:
            ax.add_patch(plt.Circle(c,eye_radius,fill=False,color='white',lw=.8));ax.plot(*moon_xy,'kx',ms=5)
            ax.set(xlim=(c[0]-radius*1.06,c[0]+radius*1.06),ylim=(c[1]+radius*1.06,c[1]-radius*1.06),aspect='equal');ax.axis('off')
        fig.suptitle(a.site+f" — Moon/Sun geometry control; Sun altitude {solar_metadata['sun_elevation_deg']:.2f}°, azimuth {solar_metadata['sun_compass_azimuth_deg']:.2f}°"+
            '\nNo fitted analyzer zero, sky rotation, source intensity or Moon/Sun mixing. Absolute AoP colors have an unknown constant zero.',fontsize=12)
        fig.savefig(a.output/'aop_moon_sun_control.png',dpi=160);plt.close(fig)
    summary=dict(site=a.site,status='diagnostic_only_no_pattern_fit',reference_lock_sha256=digest(a.reference_lock),
        script_sha256=digest(__file__),astronomical_anchor=site['astronomical_anchor'],
        observed_moon_pixel_xy=moon_xy.tolist(),fixed_ephemeris_predicted_pixel_xy=predicted_xy.tolist(),
        moon_pixel_separation=float(np.linalg.norm(predicted_xy-moon_xy)),
        angular_separation_of_two_moon_directions_deg=float(np.rad2deg(np.arccos(np.clip(primary_moon@observed_moon,-1,1)))),
        solar_control=solar_metadata,bin_size=a.bin_size,usable_central_bins=int(mask.sum()),central_bin_count=int((rho<eye_radius).sum()),
        min_dolp=a.min_dolp,excluded_source_angle_deg=10,metrics=metrics,
        fitted_pattern_rotation_deg=None,fitted_analyzer_zero_deg=None,reflection_applied=False,
        limitations=['Unknown fixed analyzer zero prevents absolute AoP color comparison; pair changes cancel that zero.',
            'Primary prediction assumes optical axis at zenith and provisional equisolid projection.',
            'Observed-Moon predictor is a separately labelled geometry sensitivity check, not an alignment fit.',
            'Rear-analyzer transport is idealized; no measured lens Mueller calibration.',
            'Differences can reflect reconstruction, instrument effects or real polarized sky light; not a unique diagnosis.',
            'Masks and nominal analyzer reconstruction remain provisional; no model or ALAN values changed.'])
    (a.output/'aop_geometry_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print('Fixed-Moon / stored-U check:',metrics['stored_U']['fixed_ephemeris_zenith_axis'][str(a.lag_bins*a.bin_size)],flush=True)
    if solar_metadata is not None:print('Fixed-Sun / stored-U check:',metrics['stored_U']['fixed_sun_geometry_control'][str(a.lag_bins*a.bin_size)],'Sun altitude:',solar_metadata['sun_elevation_deg'],flush=True)
    print('Saved '+str(a.output/'aop_geometry_audit.png')+' and aop_geometry_summary.json')


if __name__=='__main__':main()
