"""Lunar-only Berry/Wang maps from an original TIFF capture timestamp.

Theory only: does not decode measured polarization or require a NEF lock.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from anchor_moon_from_local_time import capture_time, local_to_utc
from wang2016_sky import lunar_track
from model_paper1_lunar_baseline import lunar_component
from reconstruct_spider_eye_contrast import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,help='Existing SITE/registration/rotation_analysis.json collection')
    p.add_argument('--reference-tiff',type=Path,help='Explicit original reference TIFF instead of saved batch record')
    p.add_argument('--site',default='Cottage')
    p.add_argument('--timezone',default='Europe/Stockholm')
    p.add_argument('--reference',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/prc_spider_18_sites_verified.csv')
    p.add_argument('--neutral-separation-deg',type=float,required=True)
    p.add_argument('--grid-size',type=int,default=801)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or not 0<a.neutral_separation_deg<90 or a.grid_size<101:
        p.error('Require fresh output, 0 < neutral separation < 90 and grid size >=101')
    report_path=None
    if a.reference_tiff:
        source=a.reference_tiff
    elif a.batch:
        report_path=a.batch/a.site/'registration/rotation_analysis.json'
        records=json.loads(report_path.read_text())['inputs']
        candidates=[r for r in records if Path(r['path']).name.lower()=='1.tif']
        if len(candidates)!=1:
            p.error('Need one reference 1.tif record; supply --reference-tiff')
        source=Path(candidates[0]['path'])
        if digest(source)!=candidates[0]['sha256']:
            p.error('Saved reference TIFF checksum differs')
    else:
        p.error('Supply --batch or --reference-tiff')
    stamp,origin=capture_time(source)
    # Export/edit timestamps are not silently substituted for capture times.
    if origin not in ('EXIF DateTimeOriginal','EXIF 36867'):
        p.error('DateTimeOriginal required; found '+origin)
    utc,offset=local_to_utc(stamp,a.timezone)
    rows=[r for r in csv.DictReader(a.reference.open()) if Path(r['source_filename']).stem==a.site]
    if len(rows)!=1:
        p.error('Need exactly one site reference row')
    ref=rows[0]
    if stamp[:10].replace(':','-')!=ref['field_date']:
        p.error('TIFF date differs from verified site date; reconcile before calculation')
    lat,lon=float(ref['latitude']),float(ref['longitude'])
    elevations,azimuths=lunar_track([utc],lat,lon)
    elevation,azimuth=float(elevations[0]),float(azimuths[0])
    if not 0<elevation<90:
        p.error('Reduced lunar kernel requires Moon above horizon and below zenith')
    el,az=np.deg2rad([elevation,azimuth])
    moon=np.array([-np.cos(el)*np.sin(az),np.cos(el)*np.cos(az),np.sin(el)])
    # Equisolid hemisphere; Cartesian sky axes x west/right, y north/up.
    y,x=np.indices((a.grid_size,a.grid_size),dtype=float)
    centre=(a.grid_size-1)/2;radius=centre-2
    xx,yy=(x-centre)/radius,-(y-centre)/radius
    rho=np.hypot(xx,yy);sky=rho<.999
    theta=2*np.arcsin(np.minimum(rho/np.sqrt(2),1))
    phi=np.arctan2(yy,xx)
    theta=np.where(sky,theta,0);phi=np.where(sky,phi,0)
    view=np.stack([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)],-1)
    i,z=lunar_component(view,moon,a.neutral_separation_deg)
    dop=np.where(sky,abs(z)/i,np.nan)
    angle=np.where(sky&(abs(z)>1e-10),np.rad2deg(.5*np.angle(z)),np.nan)
    field=sky&(theta<np.pi/6)
    mr=np.sqrt(1-moon[2]);mn=np.linalg.norm(moon[:2])
    mp=np.array([moon[0],-moon[1]])*mr/max(mn,1e-12)*radius+centre
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'lunar_only_stokes.npz',I=np.where(sky,i,np.nan),
        Q=np.where(sky,z.real,np.nan),U=np.where(sky,z.imag,np.nan),
        DoLP=dop,AoP_deg=angle,central_60deg_field=field)
    for name,v in [('DoLP_lunar',dop),('AoP_deg_lunar',angle)]:
        tifffile.imwrite(a.output/(name+'.tif'),v.astype(np.float32),metadata=None)
    fig,axes=plt.subplots(1,2,figsize=(10,5),layout='constrained')
    for ax,v,cmap,lo,hi,title in zip(axes,[dop,angle],['viridis','twilight'],[0,-90],[1,90],['Lunar-only DoLP','Lunar-only AoP']):
        im=ax.imshow(v,cmap=cmap,vmin=lo,vmax=hi)
        ax.add_patch(plt.Circle((centre,centre),radius*np.sqrt(2)*np.sin(np.pi/12),fill=False,color='white',lw=.8))
        ax.plot(*mp,'r+',ms=7)
        for tx,ty,t in [(.5,.99,'N'),(.01,.5,'E'),(.5,.01,'S'),(.99,.5,'W')]:
            ax.text(tx,ty,t,transform=ax.transAxes,ha='center',va='center')
        ax.set_title(title);ax.axis('off');fig.colorbar(im,ax=ax,shrink=.8)
    fig.suptitle(a.site+'; Berry/Wang theory only; no dilution or sky fit\n'+stamp+' '+a.timezone+'; Moon '+f'{elevation:.2f} deg elevation, {azimuth:.2f} deg azimuth')
    for ext in ['png','pdf','svg']:
        fig.savefig(a.output/('fig01b_lunar_only.'+ext),dpi=220)
    plt.close(fig)
    result=dict(site=a.site,status='theory_from_TIFF_capture_time_not_measured_reconstruction',
        local_time=stamp,timezone=a.timezone,utc_offset_hours=offset,utc=utc+'Z',
        timestamp_source=origin,reference_field_date=ref['field_date'],
        latitude_deg=lat,longitude_deg=lon,moon_elevation_deg=elevation,moon_compass_azimuth_deg=azimuth,
        source_tiff=str(source),source_tiff_sha256=digest(source),
        source_batch_report_sha256=digest(report_path) if report_path else None,
        site_reference_sha256=digest(a.reference),script_sha256=digest(__file__),
        kernel_sha256=digest(Path(__file__).with_name('wang2016_sky.py')),
        component_module_sha256=digest(Path(__file__).with_name('model_paper1_lunar_baseline.py')),
        neutral_separation_deg=a.neutral_separation_deg,
        theoretical_central_mean_dolp=float(dop[field].mean()),
        theoretical_central_spatial_sd=float(dop[field].std(ddof=1)),
        central_pixel_count=int(field.sum()),grid_size=a.grid_size,
        assumptions=['Same reduced Berry-normalized Wang angular/horizon component as museum.',
        'Uniform natural I=1; unity atmospheric multiplier; no imposed pmax, no ALAN, no solar contribution.',
        'Zenith-axis equisolid; N up E left; ephemeris Moon marker, not an observed pixel.',
        'No absolute experimental camera orientation or polarization calibration established by this calculation.',
        'Achromatic polarization; not an animal-specific spectral prediction.',
        'Equal-grid-pixel theoretical mean/SD; not measured Cottage statistics.'])
    (a.output/'lunar_only.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['local_time','utc','moon_elevation_deg','moon_compass_azimuth_deg','theoretical_central_mean_dolp','theoretical_central_spatial_sd']},indent=2))
    print('Saved',a.output)


if __name__=='__main__':
    main()
